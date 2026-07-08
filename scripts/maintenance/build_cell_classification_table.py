from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from cellclass.config import DEFAULT_CELL_CLASSIFICATION_TABLE, DEFAULT_TYPE_U_COMPARISON_ROOT
from placefields import load_cell_classification_table


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Aggregate per-age-group *_classification_info.csv files into a canonical "
            "session_id/cell_id/pred_type/u_type/p_pred_type table with a low-confidence flag."
        )
    )
    ap.add_argument(
        "--classification_source",
        type=str,
        default=DEFAULT_TYPE_U_COMPARISON_ROOT,
        help=(
            "CSV, parquet, or directory root containing canonical GMM comparison "
            "outputs or legacy *_classification_info.csv files."
        ),
    )
    ap.add_argument(
        "--out_csv",
        type=str,
        default=DEFAULT_CELL_CLASSIFICATION_TABLE,
        help="Output CSV path.",
    )
    return ap


def main() -> None:
    args = build_parser().parse_args()
    table = load_cell_classification_table(args.classification_source)
    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_csv, index=False)
    print(f"Wrote {len(table)} rows to {out_csv}")
    if "pred_type" in table.columns:
        counts = table["pred_type"].value_counts(dropna=False).to_dict()
        print(f"pred_type counts: {counts}")


if __name__ == "__main__":
    main()
