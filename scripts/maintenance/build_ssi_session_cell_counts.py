from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from cellclass.config import DEFAULT_CELL_CLASSIFICATION_TABLE


def parse_bool_series(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values.fillna(False)
    return values.map(lambda value: str(value).strip().lower() in {"true", "1", "yes"})


def normalize_cell_ids(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["session_id"] = out["session_id"].astype(str)
    out["cell_id"] = pd.to_numeric(out["cell_id"], errors="coerce")
    out = out.dropna(subset=["cell_id"]).copy()
    out["cell_id"] = out["cell_id"].astype(np.int64)
    return out


def normalize_session_key(value: object) -> str:
    return " ".join(str(value).strip().split())


def age_group_from_age(age: int) -> str:
    if age in (16, 17, 18):
        return "P16-18"
    if age in (19, 20, 21):
        return "P19-21"
    if age in (22, 23, 24):
        return "P22-24"
    return ""


def load_age_group_map(schedule_xlsx: Path, sheet_name: str) -> dict[str, str]:
    if not schedule_xlsx.exists():
        return {}
    d = pd.read_excel(schedule_xlsx, sheet_name=sheet_name)
    required = {"SessionName", "Age"}
    missing = required.difference(d.columns)
    if missing:
        raise KeyError(f"{schedule_xlsx} [{sheet_name}] missing required columns: {sorted(missing)}")
    d = d[["SessionName", "Age"]].copy()
    d["session_key"] = d["SessionName"].map(normalize_session_key)
    d["Age"] = pd.to_numeric(d["Age"], errors="coerce")
    d["age_group"] = d["Age"].map(lambda age: age_group_from_age(int(age)) if pd.notna(age) else "")
    return {
        str(row.session_key): str(row.age_group)
        for row in d.itertuples(index=False)
        if str(row.session_key) and str(row.age_group)
    }


def find_processed_ssi_csvs(processed_root: Path) -> list[Path]:
    return sorted(processed_root.glob("*/ssi/*/ssi_classification.csv"))


def decode_npz_json_scalar(value: np.ndarray) -> object:
    raw = np.asarray(value)
    if raw.shape == ():
        raw = raw.item()
    if isinstance(raw, bytes):
        raw = raw.decode("utf-8")
    return json.loads(str(raw))


def read_allcel_cell_ids(z: np.lib.npyio.NpzFile, allcel_path: Path) -> np.ndarray:
    if "allcel__id_cel" in z:
        raw = z["allcel__id_cel"]
    elif "allcel__id_cel__json" in z:
        raw = decode_npz_json_scalar(z["allcel__id_cel__json"])
    else:
        raise KeyError(f"{allcel_path} is missing allcel__id_cel / allcel__id_cel__json")
    return np.asarray(raw).ravel().astype(np.int64)


def pred_type_from_allcel(allcel_path: Path) -> pd.DataFrame:
    session_id = allcel_path.stem.removesuffix("_allcel")
    with np.load(allcel_path, allow_pickle=False) as z:
        cell_ids = read_allcel_cell_ids(z, allcel_path)

        pred = np.full(cell_ids.shape, "", dtype=object)
        if "allcel__pyr_uc" in z:
            pyr_mask = np.asarray(z["allcel__pyr_uc"])
            if pyr_mask.ndim > 1:
                pyr_mask = np.any(pyr_mask.astype(bool), axis=1)
            pyr_mask = np.asarray(pyr_mask).ravel().astype(bool)
            if pyr_mask.size == cell_ids.size:
                pred[pyr_mask] = "pyramidal"
        if "allcel__itn_uc" in z:
            itn_mask = np.asarray(z["allcel__itn_uc"])
            if itn_mask.ndim > 1:
                itn_mask = np.any(itn_mask.astype(bool), axis=1)
            itn_mask = np.asarray(itn_mask).ravel().astype(bool)
            if itn_mask.size == cell_ids.size:
                pred[itn_mask] = "interneuron"

        if "allcel__type_u" in z:
            type_u = np.asarray(z["allcel__type_u"]).ravel()
            if type_u.size == cell_ids.size:
                unresolved = pred == ""
                pred[unresolved & (type_u == 1)] = "pyramidal"
                pred[unresolved & (type_u == 0)] = "interneuron"

    out = pd.DataFrame(
        {
            "session_id": session_id,
            "cell_id": cell_ids,
            "pred_type": pred,
            "p_pred_type": np.nan,
            "Sure (P(pred_type) > 0.6)": False,
            "age_group": "",
        }
    )
    out = out.loc[out["pred_type"].astype(str) != ""].copy()
    return out


def fallback_cell_table_from_interim(session_id: str, interim_root: Path) -> pd.DataFrame:
    matches = sorted(interim_root.rglob(f"{session_id}_allcel.npz"))
    if not matches:
        return pd.DataFrame()
    return normalize_cell_ids(pred_type_from_allcel(matches[0]))


def age_group_for_session(cell_df: pd.DataFrame) -> str:
    if "age_group" not in cell_df.columns:
        return ""
    ages = sorted(str(age) for age in cell_df["age_group"].dropna().unique())
    return ages[0] if len(ages) == 1 else ";".join(ages)


def build_counts_for_session(
    ssi_path: Path,
    cell_table: pd.DataFrame,
    *,
    interim_root: Path,
    age_group_by_session: dict[str, str],
) -> dict[str, object]:
    ssi_df = pd.read_csv(ssi_path)
    if "session_id" in ssi_df.columns and not ssi_df.empty:
        session_id = str(ssi_df["session_id"].dropna().astype(str).iloc[0])
    else:
        session_id = ssi_path.parent.name

    session_cells = cell_table[cell_table["session_id"] == session_id].copy()
    if session_cells.empty:
        session_cells = fallback_cell_table_from_interim(session_id, interim_root)
    if session_cells.empty and {"session_id", "cell_id", "cell_type"}.issubset(ssi_df.columns):
        session_cells = ssi_df[["session_id", "cell_id", "cell_type"]].drop_duplicates().rename(
            columns={"cell_type": "pred_type"}
        )
        if "classification_certainty" in ssi_df.columns:
            certainty = (
                ssi_df[["session_id", "cell_id", "classification_certainty"]]
                .drop_duplicates(["session_id", "cell_id"])
            )
            session_cells = session_cells.merge(certainty, on=["session_id", "cell_id"], how="left")
            session_cells["Sure (P(pred_type) > 0.6)"] = (
                session_cells["classification_certainty"].astype(str).str.lower() == "sure"
            )
        else:
            session_cells["Sure (P(pred_type) > 0.6)"] = False
        session_cells = normalize_cell_ids(session_cells)

    pred_type = session_cells["pred_type"].astype(str).str.strip().str.lower()
    is_sure = parse_bool_series(session_cells["Sure (P(pred_type) > 0.6)"])

    sure_pyr_ids = session_cells.loc[(pred_type == "pyramidal") & is_sure, "cell_id"]
    uns_pyr_ids = session_cells.loc[(pred_type == "pyramidal") & ~is_sure, "cell_id"]

    interneuron_ids = set(session_cells.loc[pred_type == "interneuron", "cell_id"].astype(np.int64))

    if "SM" in ssi_df.columns:
        ssi_norm = normalize_cell_ids(ssi_df[["session_id", "cell_id", "SM"]])
        sm_ids = set(ssi_norm.loc[parse_bool_series(ssi_norm["SM"]), "cell_id"].astype(np.int64))
    else:
        sm_ids = set()
    sm_interneuron_ids = interneuron_ids.intersection(sm_ids)

    n_sure_pyr = int(sure_pyr_ids.nunique())
    n_uns_pyr = int(uns_pyr_ids.nunique())
    n_interneuron_with_significant_ssi = int(len(sm_interneuron_ids))
    age_group = age_group_for_session(session_cells)
    if not age_group:
        age_group = age_group_by_session.get(normalize_session_key(session_id), "")

    return {
        "session_id": session_id,
        "age_group": age_group,
        "n_sure_pyr": n_sure_pyr,
        "n_uns_pyr": n_uns_pyr,
        "n_interneuron_with_significant_SSI": n_interneuron_with_significant_ssi,
        "sum_of_the_previous_3_columns": (
            n_sure_pyr + n_uns_pyr + n_interneuron_with_significant_ssi
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Summarize processed SSI outputs by session, counting sure pyramidal cells, "
            "unsure pyramidal cells, and interneurons that are spatially modulated in at "
            "least one condition-direction."
        )
    )
    ap.add_argument(
        "--processed_root",
        type=str,
        default="data/processed",
        help="Root containing <mouse>/ssi/<session_id>/ssi_classification.csv outputs.",
    )
    ap.add_argument(
        "--interim_root",
        type=str,
        default="data/interim",
        help="Fallback source for allcel cell-type labels when a session is absent from the classification table.",
    )
    ap.add_argument(
        "--cell_classification_table",
        type=str,
        default=DEFAULT_CELL_CLASSIFICATION_TABLE,
        help="Canonical session_id/cell_id cell classification table.",
    )
    ap.add_argument(
        "--schedule_xlsx",
        type=str,
        default="data/schedule.xlsx",
        help="Optional Excel file containing SessionName/Age metadata for age_group fallback.",
    )
    ap.add_argument(
        "--schedule_sheet",
        type=str,
        default="VINCA",
        help="Sheet name in --schedule_xlsx.",
    )
    ap.add_argument(
        "--out_csv",
        type=str,
        default="results/tables/ssi_session_cell_counts.csv",
        help="Output summary CSV path.",
    )
    return ap


def main() -> None:
    args = build_parser().parse_args()
    processed_root = Path(args.processed_root)
    interim_root = Path(args.interim_root)
    cell_table_path = Path(args.cell_classification_table)
    schedule_xlsx = Path(args.schedule_xlsx)
    out_csv = Path(args.out_csv)

    ssi_paths = find_processed_ssi_csvs(processed_root)
    if not ssi_paths:
        raise FileNotFoundError(f"No processed SSI classification CSVs found under {processed_root}")

    cell_table = pd.read_csv(cell_table_path)
    cell_table = normalize_cell_ids(cell_table)
    age_group_by_session = load_age_group_map(schedule_xlsx, str(args.schedule_sheet))

    rows = [
        build_counts_for_session(
            path,
            cell_table,
            interim_root=interim_root,
            age_group_by_session=age_group_by_session,
        )
        for path in ssi_paths
    ]
    out = pd.DataFrame(rows).sort_values("session_id", kind="stable").reset_index(drop=True)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_csv, index=False)
    print(f"Wrote {len(out)} sessions to {out_csv}")
    print(out.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
