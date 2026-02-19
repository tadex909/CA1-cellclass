from __future__ import annotations

import argparse
from pathlib import Path
import sys
import traceback

# --- sys.path hack (so scripts/ works without pip install -e .)
THIS_DIR = Path(__file__).resolve().parent
root = THIS_DIR
while root != root.parent and not (root / "src" / "cellclass").is_dir():
    root = root.parent
src_dir = root / "src"
if not (src_dir / "cellclass").is_dir():
    raise RuntimeError(f"Could not find src/cellclass starting from {THIS_DIR}")
sys.path.insert(0, str(src_dir))

# Import your single-file processor
# IMPORTANT: adjust the import to match your filename/module:
# If your file is scripts/one_file_processing.py, import from there.
from one_file_processing import extract_one, parse_session_id  # adjust if needed


def iter_npz_files(interim_root: Path):
    # Your structure: data/interim/VS57/session_date/*.npz
    yield from sorted(interim_root.rglob("*.npz"))


def outputs_exist(npz_path: Path, processed_root: Path) -> bool:
    meta = parse_session_id(npz_path)
    session_id = meta["session_id"]
    mouse = meta["mouse"]
    feat = processed_root / mouse / "features" / f"{session_id}_features.parquet"
    acg = processed_root / mouse / "acg" / f"{session_id}_acg_counts_bin1_win50.npz"
    qc  = processed_root / mouse / "qc" / f"{session_id}_manifest.json"
    return feat.exists() and acg.exists() and qc.exists()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--interim_root", type=str, default="data/interim", help="Root folder containing interim npz files")
    ap.add_argument("--processed_root", type=str, default="data/processed", help="Where to write processed outputs")
    ap.add_argument("--skip_existing", action="store_true", help="Skip sessions already processed")
    ap.add_argument("--dry_run", action="store_true", help="List files without processing")
    # pass-through params for extract_one
    ap.add_argument("--acg_bin_ms", type=float, default=1.0)
    ap.add_argument("--acg_window_ms", type=float, default=50.0)
    ap.add_argument("--acg_smooth_ms", type=float, default=2.0)
    ap.add_argument("--waveform_fs_hz", type=float, default=25000.0)
    ap.add_argument("--waveform_trim", type=int, default=5)
    ap.add_argument("--qc_min_spikes", type=int, default=200)

    args = ap.parse_args()

    interim_root = Path(args.interim_root)
    processed_root = Path(args.processed_root)

    files = list(iter_npz_files(interim_root))
    print(f"Found {len(files)} npz files under {interim_root}")

    n_ok = 0
    n_skip = 0
    n_fail = 0

    for i, f in enumerate(files, 1):
        meta = parse_session_id(f)
        tag = f"[{i}/{len(files)}] {meta['session_id']}"
        if args.skip_existing and outputs_exist(f, processed_root):
            print(tag, "SKIP (already processed)")
            n_skip += 1
            continue

        if args.dry_run:
            print(tag, "DRY RUN:", f)
            continue

        try:
            extract_one(
                npz_path=f,
                processed_root=processed_root,
                acg_bin_ms=args.acg_bin_ms,
                acg_window_ms=args.acg_window_ms,
                acg_smooth_ms=args.acg_smooth_ms,
                waveform_fs_hz=args.waveform_fs_hz,
                waveform_trim=args.waveform_trim,
                qc_min_spikes=args.qc_min_spikes,
            )
            print(tag, "OK")
            n_ok += 1
        except Exception as e:
            print(tag, "FAIL:", repr(e))
            traceback.print_exc()
            n_fail += 1

    print("\nDone.")
    print("  OK  :", n_ok)
    print("  SKIP:", n_skip)
    print("  FAIL:", n_fail)


if __name__ == "__main__":
    main()
