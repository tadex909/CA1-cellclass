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
    "u_type",
    "p_pred_type",
    CONFIDENT_PRED_TYPE_COLUMN,
    "age_group",
)


def read_classification_frame(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix == ".csv":
        return pd.read_csv(path)
    raise ValueError(f"Unsupported classification file type: {path}")


def _normalize_pred_type(value: object) -> str:
    return str(value).strip().lower()


def _normalize_u_type_series(values: pd.Series) -> pd.Series:
    """
    Normalize legacy allcel type_u labels to readable cell-type names.

    The original allcel convention is:
      0 = interneuron
      1 = pyramidal
    """
    if pd.api.types.is_numeric_dtype(values):
        numeric = pd.to_numeric(values, errors="coerce")
        out = numeric.map({0: "interneuron", 1: "pyramidal"})
        return out.astype("string")

    text = values.astype(str).str.strip().str.lower()
    out = text.map(
        {
            "0": "interneuron",
            "interneuron": "interneuron",
            "int": "interneuron",
            "1": "pyramidal",
            "pyramidal": "pyramidal",
            "pyr": "pyramidal",
        }
    )
    return out.astype("string")


def _add_legacy_u_type(out: pd.DataFrame) -> pd.DataFrame:
    u_type = pd.Series(pd.NA, index=out.index, dtype="string")

    if "u_type" in out.columns:
        u_type = u_type.fillna(_normalize_u_type_series(out["u_type"]))

    for source_col in ("type_u_type", "type_u_binary", "allcel__type_u"):
        if source_col in out.columns:
            u_type = u_type.fillna(_normalize_u_type_series(out[source_col]))

    if u_type.notna().any():
        out["u_type"] = u_type
        return out

    return out


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
    out = _add_legacy_u_type(out)
    return out


def _merge_classification_frames(frames: list[pd.DataFrame]) -> pd.DataFrame:
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


def compact_cell_classification_table(
    df: pd.DataFrame,
    *,
    source_label: str = "<dataframe>",
) -> pd.DataFrame:
    return _merge_classification_frames(
        [_prepare_classification_frame(df, source_label=source_label)]
    )


def _comparison_frame_from_parquet(path: Path) -> pd.DataFrame:
    d = pd.read_parquet(path)
    if "age_group" not in d.columns and path.name.endswith("_gmm2_vs_type_u.parquet"):
        d["age_group"] = path.parent.name
    return d


def _classification_frames_from_directory(src: Path) -> list[pd.DataFrame]:
    all_comparison_path = src / "all_age_groups_gmm2_vs_type_u.parquet"
    if all_comparison_path.exists():
        return [
            _prepare_classification_frame(
                _comparison_frame_from_parquet(all_comparison_path),
                source_label=str(all_comparison_path),
            )
        ]

    comparison_paths = [
        path
        for path in sorted(src.rglob("*_gmm2_vs_type_u.parquet"))
        if path.name != "all_age_groups_gmm2_vs_type_u.parquet"
    ]
    if comparison_paths:
        return [
            _prepare_classification_frame(
                _comparison_frame_from_parquet(path),
                source_label=str(path),
            )
            for path in comparison_paths
        ]

    csv_paths = sorted(src.rglob("*_classification_info.csv"))
    if csv_paths:
        return [
            _prepare_classification_frame(pd.read_csv(csv_path), source_label=str(csv_path))
            for csv_path in csv_paths
        ]

    compact_csv = src / "cell_classification_table.csv"
    if compact_csv.exists():
        return [
            _prepare_classification_frame(pd.read_csv(compact_csv), source_label=str(compact_csv))
        ]

    raise FileNotFoundError(
        "No classification files found under "
        f"{src}. Expected all_age_groups_gmm2_vs_type_u.parquet, "
        "*_gmm2_vs_type_u.parquet, *_classification_info.csv, or "
        "cell_classification_table.csv."
    )


def load_cell_classification_table(source: str | Path) -> pd.DataFrame:
    src = Path(source)
    if src.is_file():
        frames = [
            _prepare_classification_frame(read_classification_frame(src), source_label=str(src))
        ]
    elif src.is_dir():
        frames = _classification_frames_from_directory(src)
    else:
        raise FileNotFoundError(f"Classification source does not exist: {src}")

    return _merge_classification_frames(frames)


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
