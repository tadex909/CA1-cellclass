from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from placefields import (
    PopulationGeometryConfig,
    build_population_geometry_from_saved_ratemap,
    canonical_condition_name,
    filter_cell_ids_by_pred_type,
    load_cell_classification_table,
    load_saved_ratemap_pack,
    load_traj_condition_names,
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
    nbin_tag = format_param_token("gbin", int(args.n_geom_bins))
    occ_tag = format_param_token("occ", f"{float(args.min_occupancy_s):g}")
    norm_tag = format_param_token("norm", "-".join(parse_csv_list(args.normalizations)))
    pred_types = parse_csv_list(getattr(args, "cell_pred_types", ""))
    if pred_types:
        pred_tag = format_param_token("cells", "-".join(pred_types))
        return f"population_geometry__{nbin_tag}__{occ_tag}__{norm_tag}__{pred_tag}__{ts}"
    return f"population_geometry__{nbin_tag}__{occ_tag}__{norm_tag}__{ts}"


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


def find_repo_root(start: Path) -> Path:
    root = start.resolve()
    if root.is_file():
        root = root.parent
    while root != root.parent:
        if (root / "pyproject.toml").exists() or (root / ".git").exists():
            return root
        root = root.parent
    return start.resolve().parent


def write_run_config(path: Path, cfg: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(cfg, f, indent=2)


def find_ratemap_files(ratemap_root: Path, sessions: set[str]) -> list[tuple[str, Path]]:
    out: list[tuple[str, Path]] = []
    for rmap_path in sorted(ratemap_root.rglob("*_rmap.npz")):
        stem = rmap_path.stem
        if not stem.endswith("_rmap"):
            continue
        session_id = stem[: -len("_rmap")]
        if sessions and session_id not in sessions:
            continue
        out.append((session_id, rmap_path))
    return out


def build_condition_name_map(traj_path: Path | None) -> dict[int, str]:
    if traj_path is None or not traj_path.exists():
        return {}

    cond_raw, names_raw = load_traj_condition_names(traj_path)
    out: dict[int, str] = {}
    for cond in np.unique(cond_raw.astype(np.int64, copy=False)):
        names = [
            canonical_condition_name(name)
            for name in names_raw[cond_raw == cond].tolist()
            if canonical_condition_name(name)
        ]
        if not names:
            continue
        counts = Counter(names)
        best = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
        out[int(cond)] = str(best)
    return out


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Build session-level population geometry matrices from saved *_rmap.npz files."
        )
    )
    ap.add_argument("--ratemap_root", type=str, default="results/ratemap")
    ap.add_argument("--out_root", type=str, default="results/population_geometry")
    ap.add_argument("--sessions", type=str, default="", help="Optional comma-separated session ids.")
    ap.add_argument("--n_geom_bins", type=int, default=10)
    ap.add_argument("--min_occupancy_s", type=float, default=0.05)
    ap.add_argument(
        "--normalizations",
        type=str,
        default="raw,mean_rate",
        help="Comma-separated list drawn from {raw,mean_rate}.",
    )
    ap.add_argument("--n_splits", type=int, default=100)
    ap.add_argument("--max_exact_splits", type=int, default=128)
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument(
        "--cell_classification_source",
        type=str,
        default="",
        help=(
            "Optional CSV or directory containing *_classification_info.csv files "
            "with session_id, cell_id, pred_type."
        ),
    )
    ap.add_argument(
        "--cell_pred_types",
        type=str,
        default="",
        help="Optional comma-separated pred_type filter, e.g. pyramidal.",
    )
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
    root = find_repo_root(Path(__file__))

    ratemap_root = Path(args.ratemap_root)
    if not ratemap_root.exists():
        raise FileNotFoundError(f"ratemap_root does not exist: {ratemap_root}")

    out_root_base = Path(args.out_root)
    out_root_base.mkdir(parents=True, exist_ok=True)
    run_id = args.run_id.strip() or default_run_id(args)
    if args.no_run_subdir:
        out_root = out_root_base
        run_id = ""
    else:
        out_root = out_root_base / run_id
        out_root.mkdir(parents=True, exist_ok=True)

    normalizations = tuple(parse_csv_list(args.normalizations))
    cell_pred_types = tuple(parse_csv_list(args.cell_pred_types))
    geom_cfg = PopulationGeometryConfig(
        n_geom_bins=int(args.n_geom_bins),
        min_occupancy_s=float(args.min_occupancy_s),
        normalizations=normalizations,
        n_splits=int(args.n_splits),
        max_exact_splits=int(args.max_exact_splits),
        seed=(None if args.seed is None else int(args.seed)),
    )
    geom_cfg.validate()

    cfg = {
        "script": "scripts/pipelines/build_population_geometry_from_ratemap.py",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": sys.version,
        "platform": platform.platform(),
        "git_head": get_git_head(root),
        "run_id": run_id,
        "run_subdir_enabled": bool(not args.no_run_subdir),
        "out_root_base": str(out_root_base),
        "out_root": str(out_root),
        "args": {
            "ratemap_root": args.ratemap_root,
            "out_root": args.out_root,
            "sessions": args.sessions,
            "n_geom_bins": int(args.n_geom_bins),
            "min_occupancy_s": float(args.min_occupancy_s),
            "normalizations": list(geom_cfg.normalizations),
            "n_splits": int(args.n_splits),
            "max_exact_splits": int(args.max_exact_splits),
            "seed": (None if args.seed is None else int(args.seed)),
            "cell_classification_source": args.cell_classification_source,
            "cell_pred_types": list(cell_pred_types),
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

    if cell_pred_types and not str(args.cell_classification_source).strip():
        raise ValueError("--cell_pred_types requires --cell_classification_source")

    classification_table: pd.DataFrame | None = None
    if str(args.cell_classification_source).strip():
        classification_table = load_cell_classification_table(args.cell_classification_source)

    sessions = set(parse_csv_list(args.sessions))
    ratemaps = find_ratemap_files(ratemap_root, sessions)
    if not ratemaps:
        print("No ratemap files found.")
        return

    index_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    for session_id, rmap_path in ratemaps:
        rel_parent = rmap_path.parent.relative_to(ratemap_root)
        out_path = out_root / rel_parent / f"{session_id}_geom.npz"

        if out_path.exists() and not args.overwrite:
            print(f"SKIP: {session_id} (exists)")
            index_rows.append(
                {
                    "session_id": session_id,
                    "status": "skip_exists",
                    "rmap_npz": str(rmap_path),
                    "geom_npz": str(out_path),
                }
            )
            continue

        if args.dry_run:
            print(f"[DRY] {session_id}: {rmap_path} -> {out_path}")
            continue

        try:
            saved = load_saved_ratemap_pack(rmap_path)
            traj_path: Path | None = None
            source_traj = saved.meta.get("source_traj_npz")
            if source_traj:
                traj_path = Path(str(source_traj))
                if not traj_path.is_absolute():
                    traj_path = root / traj_path

            condition_name_map = build_condition_name_map(traj_path)
            selected_cell_ids: np.ndarray | None = None
            if classification_table is not None and cell_pred_types:
                selected_cell_ids = filter_cell_ids_by_pred_type(
                    session_id=session_id,
                    cell_ids=saved.cell_ids,
                    classification_table=classification_table,
                    pred_types=cell_pred_types,
                )
                if selected_cell_ids.size == 0:
                    print(f"SKIP: {session_id} (no cells match pred_type filter {cell_pred_types})")
                    index_rows.append(
                        {
                            "session_id": session_id,
                            "status": "skip_no_matching_cells",
                            "rmap_npz": str(rmap_path),
                            "geom_npz": str(out_path),
                            "n_cells_total": int(saved.cell_ids.size),
                            "n_cells_used": 0,
                            "cell_pred_types": ",".join(cell_pred_types),
                        }
                    )
                    continue

            session_result = build_population_geometry_from_saved_ratemap(
                saved=saved,
                cfg=geom_cfg,
                condition_names_by_base=condition_name_map,
                selected_cell_ids=selected_cell_ids,
            )

            meta = dict(saved.meta)
            meta.update(
                {
                    "session_id": session_id,
                    "run_id": run_id,
                    "run_out_root": str(out_root),
                    "source_rmap_npz": str(rmap_path),
                    "n_geom_bins": int(session_result.xbin_centers.size),
                    "min_occupancy_s": float(geom_cfg.min_occupancy_s),
                    "normalizations": list(geom_cfg.normalizations),
                    "n_splits": int(geom_cfg.n_splits),
                    "max_exact_splits": int(geom_cfg.max_exact_splits),
                    "seed": (None if geom_cfg.seed is None else int(geom_cfg.seed)),
                    "condition_name_map": {str(k): str(v) for k, v in condition_name_map.items()},
                    "cell_classification_source": (
                        str(Path(args.cell_classification_source))
                        if str(args.cell_classification_source).strip()
                        else ""
                    ),
                    "cell_pred_types": list(cell_pred_types),
                    "n_cells_total": int(saved.cell_ids.size),
                    "n_cells_used": int(session_result.cell_ids_u.size),
                }
            )

            payload = session_result.to_payload()
            payload["meta_json"] = np.array(json.dumps(meta), dtype=np.string_)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(out_path, **payload)

            summary_rows.extend(session_result.summary_rows(session_id=session_id))
            print(f"OK: {session_id} -> {out_path}")
            index_rows.append(
                {
                    "session_id": session_id,
                    "status": "ok",
                    "rmap_npz": str(rmap_path),
                    "geom_npz": str(out_path),
                    "n_groups": int(session_result.condway_1b_g.size),
                    "n_normalizations": int(len(session_result.normalizations)),
                    "n_geom_bins": int(session_result.xbin_centers.size),
                    "n_cells_total": int(saved.cell_ids.size),
                    "n_cells_used": int(session_result.cell_ids_u.size),
                    "cell_pred_types": ",".join(cell_pred_types),
                }
            )
        except Exception as exc:
            print(f"FAIL: {session_id} -> {exc}")
            index_rows.append(
                {
                    "session_id": session_id,
                    "status": "error",
                    "rmap_npz": str(rmap_path),
                    "geom_npz": str(out_path),
                    "error": str(exc),
                }
            )

    if index_rows:
        index_path = out_root / "run_index.csv"
        pd.DataFrame(index_rows).to_csv(index_path, index=False)
        print(f"Wrote index: {index_path}")
    if summary_rows:
        summary_path = out_root / "geometry_summary.csv"
        pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
        print(f"Wrote geometry summary: {summary_path}")
    print("Done.")


if __name__ == "__main__":
    main()
