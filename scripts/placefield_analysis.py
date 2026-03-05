from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd

# Allow running from repository root without editable install.
THIS_DIR = Path(__file__).resolve().parent
root = THIS_DIR
while root != root.parent and not (root / "src" / "placefields").is_dir():
    root = root.parent
src_dir = root / "src"
if not (src_dir / "placefields").is_dir():
    raise RuntimeError(f"Could not find src/placefields starting from {THIS_DIR}")
sys.path.insert(0, str(src_dir))

from placefields import PlaceFieldConfig, analyze_session_placefields, load_session_placefield_arrays


def parse_csv_list(raw: str | None) -> list[str]:
    if not raw:
        return []
    return [x.strip() for x in raw.split(",") if x.strip()]


def find_npz_inputs(input_root: Path, sessions: list[str]) -> list[Path]:
    if input_root.is_file():
        return [input_root]

    files = sorted(input_root.rglob("*.npz"))
    if not sessions:
        return files

    out: list[Path] = []
    session_set = set(sessions)
    for p in files:
        stem = p.stem
        if stem in session_set or any(stem.startswith(f"{s}_") for s in session_set):
            out.append(p)
    return out


def analyze_one(npz_path: Path, out_root: Path, cfg: PlaceFieldConfig) -> dict[str, str]:
    arrays = load_session_placefield_arrays(npz_path)
    df, summary = analyze_session_placefields(arrays, cfg)

    session_id = npz_path.stem
    session_dir = out_root / session_id
    session_dir.mkdir(parents=True, exist_ok=True)

    per_cell_path = session_dir / f"{session_id}_placefield_summary.parquet"
    summary_path = session_dir / f"{session_id}_run_summary.json"
    df.to_parquet(per_cell_path, index=False)
    summary_path.write_text(json.dumps(summary, indent=2))

    return {
        "session_id": session_id,
        "source_npz": str(npz_path),
        "per_cell_summary": str(per_cell_path),
        "run_summary": str(summary_path),
    }


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Scaffold place-field analysis runner. "
            "Input NPZ files must contain position_cm, dt_s, spike_positions_cm, "
            "spike_cell_ids, and cell_ids."
        )
    )
    ap.add_argument(
        "--input_root",
        type=str,
        default="data/placefields/interim",
        help="Folder with session NPZ files or one NPZ file path.",
    )
    ap.add_argument(
        "--out_root",
        type=str,
        default="results/placefields",
        help="Output root for place-field summaries.",
    )
    ap.add_argument(
        "--sessions",
        type=str,
        default="",
        help="Optional comma-separated session IDs/stems to process.",
    )
    ap.add_argument("--bin_size_cm", type=float, default=2.0)
    ap.add_argument("--min_occupancy_s", type=float, default=0.1)
    ap.add_argument("--smooth_sigma_bins", type=float, default=1.5)
    ap.add_argument("--field_threshold_ratio", type=float, default=0.2)
    ap.add_argument("--min_field_bins", type=int, default=3)
    return ap


def main() -> None:
    args = build_parser().parse_args()

    cfg = PlaceFieldConfig(
        bin_size_cm=args.bin_size_cm,
        min_occupancy_s=args.min_occupancy_s,
        smooth_sigma_bins=args.smooth_sigma_bins,
        field_threshold_ratio=args.field_threshold_ratio,
        min_field_bins=args.min_field_bins,
    )
    cfg.validate()

    input_root = Path(args.input_root)
    if not input_root.exists():
        raise FileNotFoundError(f"Input path not found: {input_root}")

    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    sessions = parse_csv_list(args.sessions)
    npz_files = find_npz_inputs(input_root, sessions)
    if not npz_files:
        print("No matching NPZ files found.")
        return

    rows: list[dict[str, str]] = []
    for npz_path in npz_files:
        try:
            row = analyze_one(npz_path, out_root, cfg)
            rows.append(row)
            print(f"OK: {npz_path}")
        except Exception as exc:
            print(f"FAIL: {npz_path} -> {exc}")

    if rows:
        index_path = out_root / "run_index.csv"
        pd.DataFrame(rows).to_csv(index_path, index=False)
        print(f"Wrote index: {index_path}")
    print("Done.")


if __name__ == "__main__":
    main()
