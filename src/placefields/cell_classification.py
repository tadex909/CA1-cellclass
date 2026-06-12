from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_CLASSIFICATION_COLUMNS = ("session_id", "cell_id", "pred_type")
CONFIDENT_PRED_TYPE_COLUMN = "Sure (P(pred_type) > 0.6)"
COMPACT_CLASSIFICATION_COLUMNS = (
    "session_id",
    "cell_id",
    "pred_type",
    "p_pred_type",
    CONFIDENT_PRED_TYPE_COLUMN,
    "age_group",
)


def _normalize_pred_type(value: object) -> str:
    return str(value).strip().lower()


def _add_pred_type_probability(out: pd.DataFrame) -> pd.DataFrame:
    if "p_pred_type" in out.columns:
        out["p_pred_type"] = pd.to_numeric(out["p_pred_type"], errors="coerce")
        out[CONFIDENT_PRED_TYPE_COLUMN] = out["p_pred_type"] > 0.6
        return out

    p = pd.Series(np.nan, index=out.index, dtype=np.float64)
    if {"gmm_p_interneuron", "gmm_p_pyramidal"}.issubset(out.columns):
        p_interneuron = pd.to_numeric(out["gmm_p_interneuron"], errors="coerce")
        p_pyramidal = pd.to_numeric(out["gmm_p_pyramidal"], errors="coerce")
        is_interneuron = out["pred_type"].astype(str) == "interneuron"
        is_pyramidal = out["pred_type"].astype(str) == "pyramidal"
        p.loc[is_interneuron] = p_interneuron.loc[is_interneuron]
        p.loc[is_pyramidal] = p_pyramidal.loc[is_pyramidal]

    if "gmm_pmax" in out.columns:
        pmax = pd.to_numeric(out["gmm_pmax"], errors="coerce")
        p = p.fillna(pmax)

    out["p_pred_type"] = p
    out[CONFIDENT_PRED_TYPE_COLUMN] = out["p_pred_type"] > 0.6
    return out


def _prepare_classification_frame(d: pd.DataFrame, *, source_label: str) -> pd.DataFrame:
    missing = [c for c in REQUIRED_CLASSIFICATION_COLUMNS if c not in d.columns]
    if missing:
        raise KeyError(f"{source_label} missing required column(s): {missing}")

    out = d.copy()
    out["session_id"] = out["session_id"].astype(str).str.strip()
    out["cell_id"] = pd.to_numeric(out["cell_id"], errors="coerce").astype("Int64")
    out["pred_type"] = out["pred_type"].map(_normalize_pred_type)
    out = out.dropna(subset=["session_id", "cell_id", "pred_type"])
    out = out.loc[out["session_id"] != ""].copy()
    out = out.loc[out["pred_type"] != ""].copy()
    out["cell_id"] = out["cell_id"].astype(np.int64)
    out = _add_pred_type_probability(out)
    return out


def load_cell_classification_table(source: str | Path) -> pd.DataFrame:
    src = Path(source)
    if src.is_file():
        frames = [_prepare_classification_frame(pd.read_csv(src), source_label=str(src))]
    elif src.is_dir():
        csv_paths = sorted(src.rglob("*_classification_info.csv"))
        if not csv_paths:
            raise FileNotFoundError(f"No *_classification_info.csv files found under {src}")
        frames = [
            _prepare_classification_frame(pd.read_csv(csv_path), source_label=str(csv_path))
            for csv_path in csv_paths
        ]
    else:
        raise FileNotFoundError(f"Classification source does not exist: {src}")

    merged = pd.concat(frames, ignore_index=True)
    if merged.empty:
        return pd.DataFrame(columns=[c for c in COMPACT_CLASSIFICATION_COLUMNS if c != "age_group"])

    keep_cols = [c for c in COMPACT_CLASSIFICATION_COLUMNS if c in merged.columns]
    merged = merged[keep_cols].copy()

    conflicts = (
        merged.groupby(["session_id", "cell_id"], sort=False)["pred_type"]
        .nunique(dropna=True)
        .reset_index(name="n_pred_type")
    )
    bad = conflicts.loc[conflicts["n_pred_type"] > 1, ["session_id", "cell_id"]]
    if not bad.empty:
        sample = bad.head(10).to_dict(orient="records")
        raise ValueError(
            "Conflicting pred_type labels found for the same session_id/cell_id. "
            f"Examples: {sample}"
        )

    merged = merged.drop_duplicates(subset=["session_id", "cell_id"], keep="first")
    merged = merged.sort_values(["session_id", "cell_id"], kind="stable").reset_index(drop=True)
    return merged


def filter_cell_ids_by_pred_type(
    *,
    session_id: str,
    cell_ids: np.ndarray,
    classification_table: pd.DataFrame,
    pred_types: tuple[str, ...],
) -> np.ndarray:
    pred_keep = {_normalize_pred_type(v) for v in pred_types if str(v).strip()}
    if not pred_keep:
        return np.asarray(cell_ids, dtype=np.int64).copy()

    session_rows = classification_table.loc[
        classification_table["session_id"].astype(str) == str(session_id)
    ]
    if session_rows.empty:
        return np.empty(0, dtype=np.int64)

    pred_by_cell = {
        int(cell_id): _normalize_pred_type(pred_type)
        for cell_id, pred_type in zip(
            session_rows["cell_id"].astype(np.int64),
            session_rows["pred_type"].astype(str),
        )
    }
    keep_mask = np.asarray(
        [pred_by_cell.get(int(cell_id), "") in pred_keep for cell_id in np.asarray(cell_ids, dtype=np.int64)],
        dtype=bool,
    )
    return np.asarray(cell_ids, dtype=np.int64)[keep_mask]
