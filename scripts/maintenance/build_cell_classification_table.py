from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


THIS_DIR = Path(__file__).resolve().parent
root = THIS_DIR
while root != root.parent and not (root / "src" / "placefields").is_dir():
    root = root.parent
src_dir = root / "src"
if not (src_dir / "placefields").is_dir():
    raise RuntimeError(f"Could not find src/placefields starting from {THIS_DIR}")
sys.path.insert(0, str(src_dir))

from placefields import load_cell_classification_table  # noqa: E402


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Aggregate per-age-group *_classification_info.csv files into a canonical "
            "session_id/cell_id/pred_type/p_pred_type table with a low-confidence flag."
        )
    )
    ap.add_argument(
        "--classification_source",
        type=str,
        default="results/type_u_comparison_valero_feats_3",
        help="CSV or directory root containing *_classification_info.csv files.",
    )
    ap.add_argument(
        "--out_csv",
        type=str,
        default="results/type_u_comparison_valero_feats_3/cell_classification_table.csv",
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
