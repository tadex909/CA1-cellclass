from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from cellclass.config import (
    DEFAULT_TYPE_U_COMPARISON_FEATURES,
    DEFAULT_TYPE_U_COMPARISON_N_INIT,
    DEFAULT_TYPE_U_COMPARISON_ROOT,
    csv_join,
)


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def quote_cmd(cmd: list[str]) -> str:
    return " ".join(f'"{part}"' if any(ch.isspace() for ch in part) else part for part in cmd)


def run_step(name: str, cmd: list[str], *, cwd: Path, dry_run: bool) -> None:
    print(f"\n== {name} ==")
    print(quote_cmd(cmd))
    if dry_run:
        return
    subprocess.run(cmd, cwd=cwd, check=True)


def add_if(cmd: list[str], condition: bool, *args: str) -> None:
    if condition:
        cmd.extend(args)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Run the canonical cell-classification workflow from raw ratemap MAT files "
            "to GMM-vs-type_u comparison outputs."
        )
    )
    ap.add_argument("--raw_root", type=str, default="data/raw")
    ap.add_argument("--interim_root", type=str, default="data/interim")
    ap.add_argument("--processed_root", type=str, default="data/processed")
    ap.add_argument("--schedule", type=str, default="data/schedule.xlsx")
    ap.add_argument("--results_root", type=str, default="results")
    ap.add_argument("--comparison_out", type=str, default=DEFAULT_TYPE_U_COMPARISON_ROOT)

    ap.add_argument(
        "--no_recursive_raw",
        action="store_true",
        help="Do not recursively search --raw_root for Ratemap MAT files.",
    )
    ap.add_argument("--overwrite_interim", action="store_true")
    ap.add_argument("--include_spikes", action="store_true")
    ap.add_argument(
        "--processed_pattern",
        type=str,
        default="*_allcel.npz",
        help="Interim NPZ glob consumed by interim_to_processed.",
    )
    ap.add_argument(
        "--reprocess_existing",
        action="store_true",
        help="Recompute processed feature outputs even when they already exist.",
    )

    ap.add_argument("--no_qc", action="store_true", help="Write clean units without QC filtering.")
    ap.add_argument("--qc_strict", action="store_true")
    ap.add_argument("--age_groups", type=str, default="")
    ap.add_argument("--features", type=str, default=csv_join(DEFAULT_TYPE_U_COMPARISON_FEATURES))
    ap.add_argument("--no_log_fr", action="store_true")
    ap.add_argument("--no_standardize", action="store_true")

    ap.add_argument("--random_state", type=int, default=0)
    ap.add_argument("--n_init", type=int, default=DEFAULT_TYPE_U_COMPARISON_N_INIT)
    ap.add_argument(
        "--covariance_type",
        type=str,
        default="full",
        choices=["full", "tied", "diag", "spherical"],
    )
    ap.add_argument("--min_units", type=int, default=10)
    ap.add_argument("--dry-run", "--dry_run", dest="dry_run", action="store_true")
    return ap


def main() -> None:
    args = build_parser().parse_args()
    root = repo_root()
    python = sys.executable

    mat_to_npz = [
        python,
        "-m",
        "cellclass.mat_to_npz",
        "--mode",
        "ratemap",
        "--input",
        args.raw_root,
        "--output",
        args.interim_root,
    ]
    add_if(mat_to_npz, not args.no_recursive_raw, "--recursive")
    add_if(mat_to_npz, args.overwrite_interim, "--overwrite")
    add_if(mat_to_npz, args.include_spikes, "--include-spikes")

    interim_to_processed = [
        python,
        "scripts/pipelines/interim_to_processed.py",
        "--interim_root",
        args.interim_root,
        "--processed_root",
        args.processed_root,
        "--pattern",
        args.processed_pattern,
    ]
    add_if(interim_to_processed, not args.reprocess_existing, "--skip_existing")

    aggregate_by_age = [
        python,
        "scripts/pipelines/aggregate_by_age.py",
        "--processed_root",
        args.processed_root,
        "--excel",
        args.schedule,
        "--outdir",
        args.results_root,
    ]
    add_if(aggregate_by_age, not args.no_qc, "--qc")
    add_if(aggregate_by_age, args.qc_strict, "--qc_strict")
    add_if(aggregate_by_age, bool(args.features), "--features", args.features)
    add_if(aggregate_by_age, args.no_log_fr, "--no_log_fr")
    add_if(aggregate_by_age, args.no_standardize, "--no_standardize")

    compare_type_u = [
        python,
        "-m",
        "models.compare_type_u",
        "--results_root",
        args.results_root,
        "--out_root",
        args.comparison_out,
        "--random_state",
        str(args.random_state),
        "--n_init",
        str(args.n_init),
        "--covariance_type",
        args.covariance_type,
        "--min_units",
        str(args.min_units),
    ]
    add_if(compare_type_u, bool(args.age_groups), "--age_groups", args.age_groups)
    add_if(compare_type_u, bool(args.features), "--features", args.features)
    add_if(compare_type_u, args.no_log_fr, "--no_log_fr")
    add_if(compare_type_u, args.no_standardize, "--no_standardize")

    run_step("raw MAT -> interim allcel NPZ", mat_to_npz, cwd=root, dry_run=args.dry_run)
    run_step(
        "interim allcel NPZ -> processed features",
        interim_to_processed,
        cwd=root,
        dry_run=args.dry_run,
    )
    run_step("processed features -> age-group tables", aggregate_by_age, cwd=root, dry_run=args.dry_run)
    run_step(
        "age-group tables -> type_u comparison",
        compare_type_u,
        cwd=root,
        dry_run=args.dry_run,
    )

    if args.dry_run:
        print("\nDry run complete.")
    else:
        print(f"\nDone. Classification comparison outputs are under {args.comparison_out}.")


if __name__ == "__main__":
    main()
