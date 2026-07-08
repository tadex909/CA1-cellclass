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
    BayesianDecoderConfig,
    build_regular_xbin,
    build_trial_info_from_traj,
    canonical_condition_name,
    decode_bayesian_position_from_trials,
    decode_condway,
    ifreq_swap,
    load_traj_condition_names,
    matlab_1b_to_python_0b,
    normalize_x_to_100,
)
from placefields.interim_io import (
    find_pairs,
    load_allcel_spikes,
    load_traj_fields,
    stitch_trial_series,
)


SUMMARY_COLUMNS = [
    "session_id",
    "condway",
    "condition_label",
    "decode_groupby",
    "group_condition_families",
    "train_all_laps",
    "train_group",
    "train_group_label",
    "n_windows",
    "n_cells_used",
    "n_laps_decoded",
    "median_error_cm",
    "mean_error_cm",
    "mean_prob_actual",
]


def parse_csv_list(raw: str | None) -> list[str]:
    if not raw:
        return []
    return [x.strip() for x in raw.split(",") if x.strip()]


def parse_int_csv(raw: str | None) -> set[int]:
    out: set[int] = set()
    for tok in parse_csv_list(raw):
        out.add(int(tok))
    return out


def format_param_token(name: str, value: float | int | str) -> str:
    txt = str(value).strip().replace("-", "m").replace(".", "p").replace(" ", "_")
    return f"{name}_{txt}"


def default_run_id(args: argparse.Namespace) -> str:
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    tau_tag = format_param_token("tau", f"{float(args.tau_s):g}")
    group_tag = format_param_token("group", str(args.decode_groupby))
    family_tag = "family" if bool(args.group_condition_families) else "exact"
    train_tag = "train_all_laps" if bool(args.train_all_laps) else "leave_one_lap"
    cell_tag = format_param_token("cells", str(args.cell_selection))
    bin_tag = format_param_token("bin", f"{float(args.bin_size_cm):g}")
    smooth_tag = format_param_token("smooth", f"{float(args.smooth_sigma_bins):g}")
    xrem_tag = format_param_token("xrem", int(args.xbin_rem))
    speed_val = "none" if not np.isfinite(args.min_speed) else f"{float(args.min_speed):g}"
    speed_tag = format_param_token("minspeed", speed_val)
    return (
        f"bayes_decode__{group_tag}__{family_tag}__{train_tag}__{cell_tag}__{tau_tag}__{bin_tag}__"
        f"{smooth_tag}__{xrem_tag}__{speed_tag}__{ts}"
    )


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


def np_bytes_json(obj: Any) -> np.ndarray:
    return np.asarray(json.dumps(obj), dtype="S")


def decode_meta_json(raw: np.ndarray) -> dict[str, Any]:
    value = np.asarray(raw).item()
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    decoded = json.loads(str(value))
    return decoded if isinstance(decoded, dict) else {}


def build_condition_name_map(traj_path: Path | None) -> dict[int, str]:
    if traj_path is None or not traj_path.exists():
        return {}
    try:
        cond_raw, names_raw = load_traj_condition_names(traj_path)
    except (KeyError, ValueError):
        return {}

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


def condition_label_for_condway(condway: int, condition_name_map: dict[int, str]) -> str:
    base_condition_1b, direction = decode_condway(int(condway))
    name = condition_name_map.get(int(base_condition_1b), f"cond{base_condition_1b}")
    return f"{name} {direction}"


def parse_bool_series(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values.fillna(False)
    return values.map(lambda value: str(value).strip().lower() in {"true", "1", "yes"})


def processed_ssi_classification_path(session_id: str, processed_root: Path) -> Path:
    mouse = str(session_id).split("_", 1)[0]
    return processed_root / mouse / "ssi" / str(session_id) / "ssi_classification.csv"


def pyr_plus_sm_interneuron_cell_ids(
    *,
    session_id: str,
    all_cell_ids: np.ndarray,
    processed_root: Path,
) -> tuple[np.ndarray, dict[str, Any]]:
    ssi_path = processed_ssi_classification_path(session_id, processed_root)
    if not ssi_path.exists():
        raise FileNotFoundError(f"Missing processed SSI classification for {session_id}: {ssi_path}")

    d = pd.read_csv(ssi_path)
    required = {"cell_id", "SM"}
    missing = required.difference(d.columns)
    if missing:
        raise KeyError(f"{ssi_path} missing required columns: {sorted(missing)}")

    type_col = "cell_type" if "cell_type" in d.columns else "pred_type" if "pred_type" in d.columns else None
    if type_col is None:
        raise KeyError(f"{ssi_path} must contain cell_type or pred_type for PYR+SMINTER selection")

    d = d.copy()
    d["cell_id"] = pd.to_numeric(d["cell_id"], errors="coerce")
    d = d.dropna(subset=["cell_id"]).copy()
    d["cell_id"] = d["cell_id"].astype(np.int64)
    d = d[np.isin(d["cell_id"], np.asarray(all_cell_ids, dtype=np.int64))]
    d["cell_type_norm"] = d[type_col].astype(str).str.strip().str.lower()
    d["SM_bool"] = parse_bool_series(d["SM"])

    pyr_ids = set(d.loc[d["cell_type_norm"] == "pyramidal", "cell_id"].astype(np.int64))
    sm_inter_ids = set(
        d.loc[
            (d["cell_type_norm"] == "interneuron") & d["SM_bool"],
            "cell_id",
        ].astype(np.int64)
    )
    selected = np.asarray(sorted(pyr_ids.union(sm_inter_ids)), dtype=np.int64)
    meta = {
        "ssi_classification_csv": str(ssi_path),
        "n_pyramidal_selected": int(len(pyr_ids)),
        "n_sm_interneuron_selected": int(len(sm_inter_ids)),
    }
    return selected, meta


def selected_cell_ids_for_run(
    *,
    session_id: str,
    all_cell_ids: np.ndarray,
    requested_cell_ids: set[int],
    max_cells: int,
    cell_selection: str,
    processed_root: Path,
) -> tuple[np.ndarray, dict[str, Any]]:
    ids = np.asarray(all_cell_ids, dtype=np.int64).ravel()
    mode = str(cell_selection).strip().lower()
    selection_meta: dict[str, Any] = {"cell_selection": mode or "all"}
    if mode in {"pyr+sminter", "pyr+sm_inter", "pyr_sm_inter", "pyr+sm-inter"}:
        ids, selection_meta = pyr_plus_sm_interneuron_cell_ids(
            session_id=session_id,
            all_cell_ids=ids,
            processed_root=processed_root,
        )
        selection_meta["cell_selection"] = "PYR+SMINTER"
    elif mode not in {"all", ""}:
        raise ValueError("cell_selection must be 'all' or 'pyr+sm_inter'")
    if requested_cell_ids:
        ids = ids[np.isin(ids, list(requested_cell_ids))]
    if max_cells > 0:
        ids = ids[: int(max_cells)]
    if ids.size == 0:
        raise ValueError("No cells selected for decoding")
    selection_meta["n_cells_selected_after_limits"] = int(ids.size)
    return ids, selection_meta


def load_summary_rows_from_decode_npz(npz_path: Path, session_id: str) -> list[dict[str, Any]]:
    with np.load(npz_path, allow_pickle=False) as z:
        required = [
            "cell_ids",
            "decode__condway_w",
            "decode__train_group_w",
            "decode__trial_index_w",
            "decode__error_cm_w",
            "decode__prob_actual_w",
        ]
        missing = [k for k in required if k not in z.files]
        if missing:
            raise KeyError(f"Missing keys in {npz_path}: {missing}")

        cell_ids = np.asarray(z["cell_ids"], dtype=np.int64)
        condway_w = np.asarray(z["decode__condway_w"], dtype=np.int64)
        train_group_w = np.asarray(z["decode__train_group_w"], dtype=np.int64)
        if "decode__train_group_label_w" in z.files:
            train_group_label_w = np.asarray(z["decode__train_group_label_w"]).astype(str)
        else:
            train_group_label_w = np.asarray([str(v) for v in train_group_w], dtype=np.str_)
        trial_index_w = np.asarray(z["decode__trial_index_w"], dtype=np.int64)
        error_w = np.asarray(z["decode__error_cm_w"], dtype=np.float64)
        prob_w = np.asarray(z["decode__prob_actual_w"], dtype=np.float64)
        condition_name_map: dict[int, str] = {}
        if "meta_json" in z.files:
            try:
                meta = decode_meta_json(z["meta_json"])
                raw_map = meta.get("condition_name_map", {})
                if isinstance(raw_map, dict):
                    condition_name_map = {int(k): str(v) for k, v in raw_map.items()}
                decode_groupby = str(meta.get("decode_groupby", "unknown"))
                group_condition_families = bool(meta.get("group_condition_families", False))
                train_all_laps = bool(meta.get("train_all_laps", False))
            except Exception:
                condition_name_map = {}
                decode_groupby = "unknown"
                group_condition_families = False
                train_all_laps = False
        else:
            decode_groupby = "unknown"
            group_condition_families = False
            train_all_laps = False
        if "decode__condition_label_w" in z.files:
            condition_label_w = np.asarray(z["decode__condition_label_w"]).astype(str)
        else:
            condition_label_w = np.asarray(
                [condition_label_for_condway(c, condition_name_map) for c in condway_w],
                dtype=np.str_,
            )

    rows: list[dict[str, Any]] = []
    for condway in sorted(set(condway_w.tolist())):
        idx = condway_w == int(condway)
        labels = [v for v in sorted(set(condition_label_w[idx].tolist())) if str(v)]
        rows.append(
            {
                "session_id": str(session_id),
                "condway": int(condway),
                "condition_label": ",".join(labels)
                or condition_label_for_condway(int(condway), condition_name_map),
                "decode_groupby": decode_groupby,
                "group_condition_families": bool(group_condition_families),
                "train_all_laps": bool(train_all_laps),
                "train_group": ",".join(str(v) for v in sorted(set(train_group_w[idx].tolist()))),
                "train_group_label": ",".join(
                    v for v in sorted(set(train_group_label_w[idx].tolist())) if v
                ),
                "n_windows": int(np.sum(idx)),
                "n_cells_used": int(cell_ids.size),
                "n_laps_decoded": int(np.unique(trial_index_w[idx]).size),
                "median_error_cm": float(np.nanmedian(error_w[idx])),
                "mean_error_cm": float(np.nanmean(error_w[idx])),
                "mean_prob_actual": float(np.nanmean(prob_w[idx])),
            }
        )
    return rows


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Build leave-one-lap-out memoryless Bayesian position decoder outputs "
            "from interim allcel+trajdata pairs."
        )
    )
    ap.add_argument("--interim_root", type=str, default="data/interim")
    ap.add_argument("--out_root", type=str, default="results/position_decoding")
    ap.add_argument("--sessions", type=str, default="", help="Optional comma-separated session ids.")

    ap.add_argument("--spike_freq_hz", type=float, default=25000.0)
    ap.add_argument("--behavior_freq_hz", type=float, default=1000.0)
    ap.add_argument("--tau_s", type=float, default=0.150)
    ap.add_argument("--bin_size_cm", type=float, default=2.0)
    ap.add_argument("--smooth_sigma_bins", type=float, default=2.8)
    ap.add_argument("--xbin_rem", type=int, default=0)
    ap.add_argument("--min_speed", type=float, default=2.0)
    ap.add_argument(
        "--decode_groupby",
        type=str,
        default="condway",
        choices=["condway", "condition", "global"],
        help=(
            "Which held-out lap group defines the training map. condway keeps "
            "condition+direction separate; condition pools W/B within condition; "
            "global pools all laps except the held-out lap."
        ),
    )
    ap.add_argument(
        "--group_condition_families",
        action="store_true",
        help=(
            "Group condition variants by family before decoding. For example, "
            "PO/PO2/PO3/PONM share one group and POM/POMB share one group. "
            "With --decode_groupby condway, direction is still kept separate."
        ),
    )
    ap.add_argument(
        "--train_all_laps",
        action="store_true",
        help=(
            "Train each decoded lap on all laps in the same decode group, including "
            "the lap being decoded. Default is leave-one-lap-out."
        ),
    )
    ap.add_argument("--min_valid_window_fraction", type=float, default=0.5)
    ap.add_argument("--rate_floor_hz", type=float, default=1e-12)
    ap.add_argument("--no_normalize_x", action="store_true")

    ap.add_argument(
        "--cell_selection",
        type=str,
        default="all",
        choices=["all", "pyr+sm_inter"],
        help=(
            "Cell selection mode. 'all' uses all recorded cells. 'pyr+sm_inter' uses "
            "all pyramidal cells plus interneurons with SM=True in at least one "
            "condition-direction from data/processed/<mouse>/ssi/<session_id>/ssi_classification.csv."
        ),
    )
    ap.add_argument(
        "--processed_root",
        type=str,
        default="data/processed",
        help="Root containing processed SSI outputs used by --cell_selection pyr+sm_inter.",
    )
    ap.add_argument("--cell_ids", type=str, default="", help="Optional comma-separated cell ids.")
    ap.add_argument(
        "--max_cells",
        type=int,
        default=0,
        help="If >0, process at most this many selected cells (first by file order).",
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
    requested_cell_ids = parse_int_csv(args.cell_ids)
    processed_root = Path(args.processed_root)
    decoder_cfg = BayesianDecoderConfig(
        freq_hz=float(args.behavior_freq_hz),
        tau_s=float(args.tau_s),
        bin_size_cm=float(args.bin_size_cm),
        min_speed=min_speed_cfg,
        smooth_sigma_bins=float(args.smooth_sigma_bins),
        xbin_rem=int(args.xbin_rem),
        min_valid_window_fraction=float(args.min_valid_window_fraction),
        rate_floor_hz=float(args.rate_floor_hz),
        decode_groupby=str(args.decode_groupby),
        group_condition_families=bool(args.group_condition_families),
        train_all_laps=bool(args.train_all_laps),
    )
    decoder_cfg.validate()

    cfg = {
        "script": "scripts/pipelines/build_bayesian_decoder_from_interim.py",
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
            "tau_s": float(args.tau_s),
            "bin_size_cm": float(args.bin_size_cm),
            "smooth_sigma_bins": float(args.smooth_sigma_bins),
            "xbin_rem": int(args.xbin_rem),
            "min_speed": min_speed_cfg,
            "decode_groupby": str(args.decode_groupby),
            "group_condition_families": bool(args.group_condition_families),
            "train_all_laps": bool(args.train_all_laps),
            "min_valid_window_fraction": float(args.min_valid_window_fraction),
            "rate_floor_hz": float(args.rate_floor_hz),
            "no_normalize_x": bool(args.no_normalize_x),
            "cell_selection": str(args.cell_selection),
            "processed_root": str(args.processed_root),
            "cell_ids": args.cell_ids,
            "max_cells": int(args.max_cells),
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

    index_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []

    for session, allcel_path, traj_path in pairs:
        rel_parent = allcel_path.parent.relative_to(interim_root)
        out_path = out_root / rel_parent / f"{session}_bayes_decode.npz"

        if out_path.exists() and not args.overwrite:
            print(f"SKIP: {session} (exists)")
            try:
                summary_rows.extend(load_summary_rows_from_decode_npz(out_path, session))
            except Exception as exc:
                print(f"WARN: could not collect decoding summary from existing file ({out_path}): {exc}")
            index_rows.append(
                {
                    "session_id": session,
                    "status": "skip_exists",
                    "allcel_npz": str(allcel_path),
                    "traj_npz": str(traj_path),
                    "decode_npz": str(out_path),
                }
            )
            continue

        if args.dry_run:
            print(f"[DRY] {session}: {allcel_path} + {traj_path} -> {out_path}")
            continue

        try:
            itime_25k, id_spk, id_cel = load_allcel_spikes(allcel_path)
            cond, wb, start, stop, vr_list, speed_list = load_traj_fields(traj_path)
            condition_name_map = build_condition_name_map(traj_path)

            x = stitch_trial_series(start, stop, vr_list)
            if not args.no_normalize_x:
                x = normalize_x_to_100(x)

            if speed_list is None and min_speed_cfg is not None:
                raise ValueError(
                    "Trajectory file has no XSpeed/Speed field; use --min_speed nan to decode without speed filtering."
                )
            speed = stitch_trial_series(start, stop, speed_list) if speed_list is not None else None

            trials = build_trial_info_from_traj(
                cond=cond,
                wb=wb,
                start_1b=start,
                stop_1b=stop,
                n_samples=int(x.size),
            )
            xbin_edges = build_regular_xbin(x, float(args.bin_size_cm))
            spike_idx_1b = ifreq_swap(itime_25k, args.spike_freq_hz, args.behavior_freq_hz)
            spike_idx_0b = matlab_1b_to_python_0b(spike_idx_1b)
            cell_ids_decode, cell_selection_meta = selected_cell_ids_for_run(
                session_id=session,
                all_cell_ids=id_cel,
                requested_cell_ids=requested_cell_ids,
                max_cells=int(args.max_cells),
                cell_selection=str(args.cell_selection),
                processed_root=processed_root,
            )

            result = decode_bayesian_position_from_trials(
                position_x=x,
                spike_indices_0b=spike_idx_0b,
                spike_cell_ids=id_spk,
                cell_ids=cell_ids_decode,
                trials=trials,
                xbin_edges=xbin_edges,
                cfg=decoder_cfg,
                speed=speed,
                condition_names_by_base=condition_name_map,
            )

            meta = {
                "session_id": session,
                "run_id": run_id,
                "run_out_root": str(out_root),
                "source_allcel_npz": str(allcel_path),
                "source_traj_npz": str(traj_path),
                "freq_spike_hz": float(args.spike_freq_hz),
                "freq_behavior_hz": float(args.behavior_freq_hz),
                "tau_s": float(args.tau_s),
                "bin_size_cm": float(args.bin_size_cm),
                "smooth_sigma_bins": float(args.smooth_sigma_bins),
                "xbin_rem": int(args.xbin_rem),
                "min_speed": min_speed_cfg,
                "decode_groupby": str(args.decode_groupby),
                "group_condition_families": bool(args.group_condition_families),
                "train_all_laps": bool(args.train_all_laps),
                "min_valid_window_fraction": float(args.min_valid_window_fraction),
                "rate_floor_hz": float(args.rate_floor_hz),
                "normalize_x_to_100": bool(not args.no_normalize_x),
                "cell_selection": str(args.cell_selection),
                "cell_selection_meta": cell_selection_meta,
                "requested_cell_ids": sorted(requested_cell_ids),
                "max_cells": int(args.max_cells),
                "n_cells_requested": int(cell_ids_decode.size),
                "n_cells_used": int(result.cell_ids.size),
                "n_trials": int(len(trials)),
                "n_windows": int(result.posterior_wx.shape[0]),
                "n_bins": int(result.xbin_centers.size),
                "cv": (
                    f"train_all_laps_within_{args.decode_groupby}"
                    if bool(args.train_all_laps)
                    else f"leave_one_lap_within_{args.decode_groupby}"
                ),
                "condition_name_map": {str(k): v for k, v in sorted(condition_name_map.items())},
            }

            payload = result.to_payload()
            payload["decode__condition_label_w"] = np.asarray(
                [
                    condition_label_for_condway(int(condway), condition_name_map)
                    for condway in result.condway_w
                ],
                dtype=np.str_,
            )
            payload["meta_json"] = np_bytes_json(meta)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(out_path, **payload)

            for row in result.summary_rows(session_id=session):
                row["decode_groupby"] = str(args.decode_groupby)
                row["group_condition_families"] = bool(args.group_condition_families)
                row["train_all_laps"] = bool(args.train_all_laps)
                row["condition_label"] = condition_label_for_condway(
                    int(row["condway"]),
                    condition_name_map,
                )
                summary_rows.append(row)
            print(f"OK: {session} -> {out_path}")
            index_rows.append(
                {
                    "session_id": session,
                    "status": "ok",
                    "allcel_npz": str(allcel_path),
                    "traj_npz": str(traj_path),
                    "decode_npz": str(out_path),
                    "n_cells_requested": int(cell_ids_decode.size),
                    "n_cells_used": int(result.cell_ids.size),
                    "n_trials": int(len(trials)),
                    "n_windows": int(result.posterior_wx.shape[0]),
                    "n_bins": int(result.xbin_centers.size),
                    "decode_groupby": str(args.decode_groupby),
                    "group_condition_families": bool(args.group_condition_families),
                    "train_all_laps": bool(args.train_all_laps),
                    "cell_selection": str(args.cell_selection),
                    **cell_selection_meta,
                }
            )
        except Exception as exc:
            print(f"FAIL: {session} -> {exc}")
            index_rows.append(
                {
                    "session_id": session,
                    "status": "error",
                    "allcel_npz": str(allcel_path),
                    "traj_npz": str(traj_path),
                    "decode_npz": str(out_path),
                    "error": str(exc),
                }
            )

    if index_rows:
        index_path = out_root / "run_index.csv"
        pd.DataFrame(index_rows).to_csv(index_path, index=False)
        print(f"Wrote index: {index_path}")

    summary_path = out_root / "decoding_summary.csv"
    pd.DataFrame(summary_rows, columns=SUMMARY_COLUMNS).to_csv(summary_path, index=False)
    print(f"Wrote decoding summary: {summary_path}")
    print("Done.")


if __name__ == "__main__":
    main()
