from __future__ import annotations

import argparse
import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path
import sys
from typing import Any

import numpy as np
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

from placefields import (
    RatemapConfig,
    build_default_xbin,
    build_ratemap_from_trials,
    build_trial_info_from_traj,
    ifreq_swap,
    matlab_1b_to_python_0b,
    normalize_x_to_100,
)
from placefields.interim_io import (
    find_pairs,
    load_allcel_spikes,
    load_traj_fields,
    save_ratemap_pack,
    stitch_trial_series,
)


def parse_csv_list(raw: str | None) -> list[str]:
    if not raw:
        return []
    return [x.strip() for x in raw.split(",") if x.strip()]


def format_param_token(name: str, value: float | int | str) -> str:
    txt = str(value).strip().replace("-", "m").replace(".", "p").replace(" ", "_")
    return f"{name}_{txt}"


def default_run_id(args: argparse.Namespace) -> str:
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    smooth_tag = format_param_token("smooth", f"{float(args.smooth_sigma_bins):g}")
    xrem_tag = format_param_token("xrem", int(args.xbin_rem))
    speed_val = "none" if not np.isfinite(args.min_speed) else f"{float(args.min_speed):g}"
    speed_tag = format_param_token("minspeed", speed_val)
    return f"ratemap__{smooth_tag}__{xrem_tag}__{speed_tag}__{ts}"


def get_git_head(repo_root: Path) -> str | None:
    try:
        proc = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=str(repo_root),
            check=True,
            capture_output=True,
            text=True,
        )
        return proc.stdout.strip() or None
    except Exception:
        return None


def write_run_config(path: Path, cfg: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(cfg, f, indent=2)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Pair *_allcel.npz and *_trajdata.npz from data/interim and build ratemap "
            "tensors using placefields.build_ratemap_from_trials."
        )
    )
    ap.add_argument("--interim_root", type=str, default="data/interim")
    ap.add_argument("--out_root", type=str, default="results/ratemap")
    ap.add_argument("--sessions", type=str, default="", help="Optional comma-separated session ids.")
    ap.add_argument("--spike_freq_hz", type=float, default=25000.0)
    ap.add_argument("--behavior_freq_hz", type=float, default=1000.0)
    ap.add_argument("--smooth_sigma_bins", type=float, default=2.8)
    ap.add_argument("--xbin_rem", type=int, default=10)
    ap.add_argument("--min_speed", type=float, default=2.0)
    ap.add_argument("--no_normalize_x", action="store_true")
    ap.add_argument(
        "--run_id",
        type=str,
        default="",
        help=(
            "Optional run identifier. If omitted and run subdirectories are enabled, "
            "a timestamped id is generated."
        ),
    )
    ap.add_argument(
        "--no_run_subdir",
        action="store_true",
        help="Write outputs directly into out_root (legacy behavior).",
    )
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--dry_run", action="store_true")
    return ap


def main() -> None:
    args = build_parser().parse_args()

    interim_root = Path(args.interim_root)
    out_root_base = Path(args.out_root)
    out_root_base.mkdir(parents=True, exist_ok=True)
    run_id = args.run_id.strip() or default_run_id(args)
    if args.no_run_subdir:
        out_root = out_root_base
        run_id = ""
    else:
        out_root = out_root_base / run_id
        out_root.mkdir(parents=True, exist_ok=True)

    min_speed_cfg = None if not np.isfinite(args.min_speed) else float(args.min_speed)
    cfg = {
        "script": "scripts/pipelines/build_ratemap_from_interim.py",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": sys.version,
        "platform": platform.platform(),
        "git_head": get_git_head(root),
        "run_id": run_id,
        "run_subdir_enabled": bool(not args.no_run_subdir),
        "out_root_base": str(out_root_base),
        "out_root": str(out_root),
        "args": {
            "interim_root": args.interim_root,
            "out_root": args.out_root,
            "sessions": args.sessions,
            "spike_freq_hz": float(args.spike_freq_hz),
            "behavior_freq_hz": float(args.behavior_freq_hz),
            "smooth_sigma_bins": float(args.smooth_sigma_bins),
            "xbin_rem": int(args.xbin_rem),
            "min_speed": min_speed_cfg,
            "no_normalize_x": bool(args.no_normalize_x),
            "overwrite": bool(args.overwrite),
            "dry_run": bool(args.dry_run),
            "run_id": args.run_id,
            "no_run_subdir": bool(args.no_run_subdir),
        },
    }
    write_run_config(out_root / "run_config.json", cfg)
    if not args.no_run_subdir:
        (out_root_base / "LATEST_RUN.txt").write_text(f"{run_id}\n", encoding="utf-8")
    print(f"Output root: {out_root}")

    sessions = set(parse_csv_list(args.sessions))
    pairs = find_pairs(interim_root, sessions)
    if not pairs:
        print("No allcel+trajdata pairs found.")
        return

    rows: list[dict[str, Any]] = []
    for session, allcel_path, traj_path in pairs:
        rel_parent = allcel_path.parent.relative_to(interim_root)
        out_path = out_root / rel_parent / f"{session}_rmap.npz"

        if out_path.exists() and not args.overwrite:
            print(f"SKIP: {session} (exists)")
            rows.append(
                {
                    "session_id": session,
                    "status": "skip_exists",
                    "allcel_npz": str(allcel_path),
                    "traj_npz": str(traj_path),
                    "rmap_npz": str(out_path),
                }
            )
            continue

        if args.dry_run:
            print(f"[DRY] {session}: {allcel_path} + {traj_path} -> {out_path}")
            continue

        try:
            itime_25k, id_spk, id_cel = load_allcel_spikes(allcel_path)
            cond, wb, start, stop, vr_list, speed_list = load_traj_fields(traj_path)

            x = stitch_trial_series(start, stop, vr_list)
            if not args.no_normalize_x:
                x = normalize_x_to_100(x)
            speed = (
                stitch_trial_series(start, stop, speed_list)
                if speed_list is not None
                else None
            )

            trials = build_trial_info_from_traj(
                cond=cond,
                wb=wb,
                start_1b=start,
                stop_1b=stop,
                n_samples=x.size,
            )
            xbin = build_default_xbin(x)
            spike_idx_1b = ifreq_swap(itime_25k, args.spike_freq_hz, args.behavior_freq_hz)
            spike_idx_0b = matlab_1b_to_python_0b(spike_idx_1b)

            min_speed = None if not np.isfinite(args.min_speed) else float(args.min_speed)
            cfg = RatemapConfig(
                freq_hz=float(args.behavior_freq_hz),
                smooth_sigma_bins=float(args.smooth_sigma_bins),
                xbin_rem=int(args.xbin_rem),
                nb_cond=None,
                min_speed=min_speed,
            )
            pack = build_ratemap_from_trials(
                position_x=x,
                spike_indices_0b=spike_idx_0b,
                spike_cell_ids=id_spk,
                cell_ids=id_cel,
                trials=trials,
                xbin_edges=xbin,
                cfg=cfg,
                speed=speed,
            )

            meta = {
                "session_id": session,
                "run_id": run_id,
                "run_out_root": str(out_root),
                "source_allcel_npz": str(allcel_path),
                "source_traj_npz": str(traj_path),
                "freq_spike_hz": float(args.spike_freq_hz),
                "freq_behavior_hz": float(args.behavior_freq_hz),
                "normalize_x_to_100": bool(not args.no_normalize_x),
                "xbin_rem": int(args.xbin_rem),
                "smooth_sigma_bins": float(args.smooth_sigma_bins),
                "min_speed": (None if min_speed is None else float(min_speed)),
                "n_cells": int(pack.cell_ids.size),
                "n_trials": int(len(pack.trial_info)),
                "n_bins": int(pack.xbin_centers.size),
                "nb_cond": int(pack.nb_cond),
            }
            save_ratemap_pack(pack, out_path, meta)
            print(f"OK: {session} -> {out_path}")
            rows.append(
                {
                    "session_id": session,
                    "status": "ok",
                    "allcel_npz": str(allcel_path),
                    "traj_npz": str(traj_path),
                    "rmap_npz": str(out_path),
                    "n_cells": int(pack.cell_ids.size),
                    "n_trials": int(len(pack.trial_info)),
                    "n_bins": int(pack.xbin_centers.size),
                    "nb_cond": int(pack.nb_cond),
                }
            )
        except Exception as exc:
            print(f"FAIL: {session} -> {exc}")
            rows.append(
                {
                    "session_id": session,
                    "status": "error",
                    "allcel_npz": str(allcel_path),
                    "traj_npz": str(traj_path),
                    "rmap_npz": str(out_path),
                    "error": str(exc),
                }
            )

    if rows:
        index_path = out_root / "run_index.csv"
        pd.DataFrame(rows).to_csv(index_path, index=False)
        print(f"Wrote index: {index_path}")
    print("Done.")


if __name__ == "__main__":
    main()
