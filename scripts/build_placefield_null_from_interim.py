from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
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
    NullBootstrapConfig,
    RatemapConfig,
    build_default_xbin,
    build_ratemap_from_trials,
    build_trial_info_from_traj,
    compute_condition_occupancy,
    empirical_pvalue_from_null,
    empirical_pval_cx,
    empirical_pval_tx,
    ifreq_swap,
    matlab_1b_to_python_0b,
    normalize_x_to_100,
    simulate_null_fr_s_txrep,
    ssi_null_by_condition,
    ssi_observed_by_condition,
)
from placefields.interim_io import (
    find_pairs,
    load_allcel_spikes,
    load_traj_fields,
    stitch_trial_series,
)


def parse_csv_list(raw: str | None) -> list[str]:
    if not raw:
        return []
    return [x.strip() for x in raw.split(",") if x.strip()]


def parse_int_csv(raw: str | None) -> set[int]:
    out: set[int] = set()
    for tok in parse_csv_list(raw):
        out.add(int(tok))
    return out


def build_ssi_classification_rows(
    *,
    session_id: str,
    cell_ids: np.ndarray,
    ssi_obs_cu: np.ndarray,
    ssi_pval_cu: np.ndarray,
    sm_alpha: float,
    source: str,
) -> list[dict[str, Any]]:
    if ssi_obs_cu.shape != ssi_pval_cu.shape:
        raise ValueError(
            f"Shape mismatch: ssi_obs {ssi_obs_cu.shape} vs ssi_pval {ssi_pval_cu.shape}"
        )
    if ssi_obs_cu.shape[0] != int(cell_ids.size):
        raise ValueError(
            f"Cell axis mismatch: {ssi_obs_cu.shape[0]} rows vs {int(cell_ids.size)} cell ids"
        )

    rows: list[dict[str, Any]] = []
    n_cond = int(ssi_obs_cu.shape[1])
    for u, cid in enumerate(cell_ids.astype(np.int64, copy=False)):
        for c in range(n_cond):
            ssi_obs = float(ssi_obs_cu[u, c])
            pval = float(ssi_pval_cu[u, c])
            is_sm = bool(np.isfinite(pval) and (pval < float(sm_alpha)))
            rows.append(
                {
                    "session_id": str(session_id),
                    "cell_id": int(cid),
                    "condition_1b": int(c + 1),
                    "ssi_obs": ssi_obs,
                    "p_value": pval,
                    "SM": is_sm,
                    "sm_alpha": float(sm_alpha),
                    "source": str(source),
                }
            )
    return rows


def load_ssi_arrays_from_pfnull(npz_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with np.load(npz_path, allow_pickle=False) as z:
        required = ["cell_ids", "pfnull__ssi_obs_cu", "pfnull__ssi_pval_cu"]
        missing = [k for k in required if k not in z.files]
        if missing:
            raise KeyError(f"Missing keys in {npz_path}: {missing}")
        cell_ids = z["cell_ids"].astype(np.int64, copy=False)
        ssi_obs_cu = z["pfnull__ssi_obs_cu"].astype(np.float64, copy=False)
        ssi_pval_cu = z["pfnull__ssi_pval_cu"].astype(np.float64, copy=False)
    return cell_ids, ssi_obs_cu, ssi_pval_cu


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Build MATLAB-like null placefield maps and empirical p-values from interim allcel+traj pairs."
        )
    )
    ap.add_argument("--interim_root", type=str, default="data/interim")
    ap.add_argument("--out_root", type=str, default="results/placefield_null")
    ap.add_argument("--sessions", type=str, default="", help="Optional comma-separated session ids.")

    ap.add_argument("--spike_freq_hz", type=float, default=25000.0)
    ap.add_argument("--behavior_freq_hz", type=float, default=1000.0)
    ap.add_argument("--smooth_sigma_bins", type=float, default=2.8)
    ap.add_argument("--xbin_rem", type=int, default=10)
    ap.add_argument("--min_speed", type=float, default=2.0)
    ap.add_argument("--no_normalize_x", action="store_true")

    ap.add_argument(
        "--null_method",
        type=str,
        default="random",
        choices=["random", "poisson", "circular_shift"],
    )
    ap.add_argument("--nb_rep", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--min_shift_frac", type=float, default=0.10)
    ap.add_argument("--min_shift_samples", type=int, default=1)

    ap.add_argument("--cell_ids", type=str, default="", help="Optional comma-separated cell ids.")
    ap.add_argument(
        "--max_cells",
        type=int,
        default=0,
        help="If >0 and --cell_ids is empty, process at most this many cells (first by file order).",
    )
    ap.add_argument("--save_null_maps", action="store_true", help="Save fr_s_txrep arrays (large files).")
    ap.add_argument(
        "--save_ssi_null",
        action="store_true",
        help="Save SSI null distributions per cell/condition (n_cells, n_cond, n_rep).",
    )
    ap.add_argument(
        "--sm_alpha",
        type=float,
        default=0.05,
        help="Threshold used to classify SM=True when SSI p-value < sm_alpha.",
    )
    ap.add_argument(
        "--classification_csv",
        type=str,
        default="ssi_classification.csv",
        help=(
            "Output CSV (relative to out_root unless absolute path) with one row "
            "per session/cell/condition and columns including p_value and SM."
        ),
    )
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--dry_run", action="store_true")
    return ap


def main() -> None:
    args = build_parser().parse_args()
    if not (0.0 < float(args.sm_alpha) < 1.0):
        raise ValueError("--sm_alpha must be strictly between 0 and 1.")

    interim_root = Path(args.interim_root)
    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    sessions = set(parse_csv_list(args.sessions))
    requested_cell_ids = parse_int_csv(args.cell_ids)

    pairs = find_pairs(interim_root, sessions)
    if not pairs:
        print("No allcel+trajdata pairs found.")
        return

    rows: list[dict[str, Any]] = []
    classification_rows: list[dict[str, Any]] = []
    for session, allcel_path, traj_path in pairs:
        rel_parent = allcel_path.parent.relative_to(interim_root)
        out_path = out_root / rel_parent / f"{session}_pfnull.npz"

        if out_path.exists() and not args.overwrite:
            print(f"SKIP: {session} (exists)")
            try:
                cell_ids_skip, ssi_obs_skip, ssi_pval_skip = load_ssi_arrays_from_pfnull(out_path)
                classification_rows.extend(
                    build_ssi_classification_rows(
                        session_id=session,
                        cell_ids=cell_ids_skip,
                        ssi_obs_cu=ssi_obs_skip,
                        ssi_pval_cu=ssi_pval_skip,
                        sm_alpha=float(args.sm_alpha),
                        source="skip_exists",
                    )
                )
            except Exception as exc:
                print(f"WARN: could not collect SSI classification from existing file ({out_path}): {exc}")
            rows.append(
                {
                    "session_id": session,
                    "status": "skip_exists",
                    "allcel_npz": str(allcel_path),
                    "traj_npz": str(traj_path),
                    "pfnull_npz": str(out_path),
                }
            )
            continue

        if args.dry_run:
            print(f"[DRY] {session}: {allcel_path} + {traj_path} -> {out_path}")
            continue

        try:
            # Load spikes and trajectory.
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
            n_samples = int(x.size)

            trials = build_trial_info_from_traj(
                cond=cond,
                wb=wb,
                start_1b=start,
                stop_1b=stop,
                n_samples=n_samples,
            )
            xbin = build_default_xbin(x)

            # Build ratemap tensors once (for dwell + observed fr_s used in p-values).
            spike_idx_1b = ifreq_swap(itime_25k, args.spike_freq_hz, args.behavior_freq_hz)
            spike_idx_0b = matlab_1b_to_python_0b(spike_idx_1b)

            min_speed = None if not np.isfinite(args.min_speed) else float(args.min_speed)
            rmap_cfg = RatemapConfig(
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
                cfg=rmap_cfg,
                speed=speed,
            )

            # Rebuild sample validity mask to match ratemap filtering.
            sample_keep = np.isfinite(x)
            if (speed is not None) and (min_speed is not None):
                sample_keep &= np.isfinite(speed) & (speed >= float(min_speed))

            # Select cells.
            cell_ids_all = pack.cell_ids.astype(np.int64)
            if requested_cell_ids:
                sel_mask = np.isin(cell_ids_all, list(requested_cell_ids))
                sel = np.where(sel_mask)[0]
            else:
                sel = np.arange(cell_ids_all.size, dtype=np.int64)
                if args.max_cells > 0:
                    sel = sel[: int(args.max_cells)]

            if sel.size == 0:
                raise RuntimeError("No cells selected for null simulation.")

            nb_sel = int(sel.size)
            n_trials = int(len(pack.trial_info))
            n_bins = int(pack.xbin_centers.size)
            nb_cond = int(pack.nb_cond)
            occupancy_cx = compute_condition_occupancy(
                pack.dwell_tx_x,
                pack.idcond_t,
                pack.nb_cond,
            )

            pval_tx_ux = np.full((nb_sel, n_trials, n_bins), np.nan, dtype=np.float64)
            pval_cx_ux = np.full((nb_sel, nb_cond, n_bins), np.nan, dtype=np.float64)
            ssi_obs_cu = np.full((nb_sel, nb_cond), np.nan, dtype=np.float64)
            ssi_pval_cu = np.full((nb_sel, nb_cond), np.nan, dtype=np.float64)
            ssi_null_cur = [] if args.save_ssi_null else None
            fr_s_txrep_u = [] if args.save_null_maps else None

            base_seed = args.seed
            for j, u in enumerate(sel):
                cid = int(cell_ids_all[u])
                spike_idx_cell = spike_idx_0b[id_spk.astype(np.int64) == cid]

                seed_j = None if base_seed is None else int(base_seed) + int(j)
                null_cfg = NullBootstrapConfig(
                    method=args.null_method,
                    nb_rep=int(args.nb_rep),
                    min_shift_frac=float(args.min_shift_frac),
                    min_shift_samples=int(args.min_shift_samples),
                    seed=seed_j,
                )

                fr_s_txrep = simulate_null_fr_s_txrep(
                    position_x=x,
                    spike_indices_cell_0b=spike_idx_cell,
                    trials=pack.trial_info,
                    xbin_edges=pack.xbin_edges,
                    dwell_s_tx_x=pack.dwell_s_tx_x,
                    smooth_sigma_bins=float(args.smooth_sigma_bins),
                    cfg=null_cfg,
                    valid_sample_mask=sample_keep,
                )

                pval_tx_ux[j] = empirical_pval_tx(fr_s_txrep, pack.fr_s_tx_ux[u])
                pval_cx_ux[j] = empirical_pval_cx(
                    fr_s_txrep,
                    pack.fr_s_cx_ux[u],
                    pack.idcond_t,
                    pack.nb_cond,
                )
                ssi_obs = ssi_observed_by_condition(
                    pack.fr_s_cx_ux[u],
                    occupancy_cx,
                )
                ssi_null = ssi_null_by_condition(
                    fr_s_txrep,
                    pack.idcond_t,
                    pack.nb_cond,
                    occupancy_cx,
                )
                ssi_obs_cu[j] = ssi_obs
                ssi_pval_cu[j] = empirical_pvalue_from_null(
                    ssi_obs,
                    ssi_null,
                    alternative="greater",
                    add_one=True,
                )

                if fr_s_txrep_u is not None:
                    fr_s_txrep_u.append(fr_s_txrep.astype(np.float32, copy=False))
                if ssi_null_cur is not None:
                    ssi_null_cur.append(ssi_null.astype(np.float32, copy=False))

            meta = {
                "session_id": session,
                "source_allcel_npz": str(allcel_path),
                "source_traj_npz": str(traj_path),
                "null_method": args.null_method,
                "nb_rep": int(args.nb_rep),
                "seed": (None if args.seed is None else int(args.seed)),
                "min_shift_frac": float(args.min_shift_frac),
                "min_shift_samples": int(args.min_shift_samples),
                "freq_spike_hz": float(args.spike_freq_hz),
                "freq_behavior_hz": float(args.behavior_freq_hz),
                "normalize_x_to_100": bool(not args.no_normalize_x),
                "min_speed": (None if min_speed is None else float(min_speed)),
                "smooth_sigma_bins": float(args.smooth_sigma_bins),
                "xbin_rem": int(args.xbin_rem),
                "n_cells_selected": int(nb_sel),
                "n_trials": int(n_trials),
                "n_bins": int(n_bins),
                "nb_cond": int(nb_cond),
                "ssi_alternative": "greater",
                "sm_alpha": float(args.sm_alpha),
            }

            payload: dict[str, Any] = {
                "meta_json": np.array(json.dumps(meta), dtype=np.string_),
                "cell_ids": cell_ids_all[sel],
                "idcond_t": pack.idcond_t.astype(np.int64, copy=False),
                "xbin_edges": pack.xbin_edges.astype(np.float64, copy=False),
                "xbin_centers": pack.xbin_centers.astype(np.float64, copy=False),
                # Observed dwell maps (seconds): trial-level and condition-level.
                "dwell_tx_x": pack.dwell_tx_x.astype(np.float32, copy=False),
                "dwell_cx_x": pack.dwell_cx_x.astype(np.float32, copy=False),
                "occupancy_cx": occupancy_cx.astype(np.float32, copy=False),
                # Observed maps for selected cells (same xbin_rem/smoothing as this bootstrap run).
                "pfnull__fr_tx_ux": pack.fr_tx_ux[sel].astype(np.float32, copy=False),
                "pfnull__fr_s_tx_ux": pack.fr_s_tx_ux[sel].astype(np.float32, copy=False),
                "pfnull__fr_cx_ux": pack.fr_cx_ux[sel].astype(np.float32, copy=False),
                "pfnull__fr_s_cx_ux": pack.fr_s_cx_ux[sel].astype(np.float32, copy=False),
                "pfnull__pval_tx_ux": pval_tx_ux.astype(np.float32, copy=False),
                "pfnull__pval_cx_ux": pval_cx_ux.astype(np.float32, copy=False),
                "pfnull__ssi_obs_cu": ssi_obs_cu.astype(np.float32, copy=False),
                "pfnull__ssi_pval_cu": ssi_pval_cu.astype(np.float32, copy=False),
            }
            if fr_s_txrep_u is not None:
                payload["pfnull__fr_s_txrep_uxtr"] = np.asarray(fr_s_txrep_u, dtype=np.float32)
            if ssi_null_cur is not None:
                payload["pfnull__ssi_null_cur"] = np.asarray(ssi_null_cur, dtype=np.float32)

            out_path.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(out_path, **payload)

            classification_rows.extend(
                build_ssi_classification_rows(
                    session_id=session,
                    cell_ids=cell_ids_all[sel],
                    ssi_obs_cu=ssi_obs_cu,
                    ssi_pval_cu=ssi_pval_cu,
                    sm_alpha=float(args.sm_alpha),
                    source="computed",
                )
            )

            print(f"OK: {session} -> {out_path}")
            rows.append(
                {
                    "session_id": session,
                    "status": "ok",
                    "allcel_npz": str(allcel_path),
                    "traj_npz": str(traj_path),
                    "pfnull_npz": str(out_path),
                    "n_cells_selected": int(nb_sel),
                    "n_trials": int(n_trials),
                    "n_bins": int(n_bins),
                    "nb_cond": int(nb_cond),
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
                    "pfnull_npz": str(out_path),
                    "error": str(exc),
                }
            )

    if rows:
        index_path = out_root / "run_index.csv"
        pd.DataFrame(rows).to_csv(index_path, index=False)
        print(f"Wrote index: {index_path}")
    if classification_rows:
        class_path = Path(args.classification_csv)
        if not class_path.is_absolute():
            class_path = out_root / class_path
        class_path.parent.mkdir(parents=True, exist_ok=True)
        dcls = pd.DataFrame(classification_rows).sort_values(
            ["session_id", "cell_id", "condition_1b"], kind="stable"
        )
        dcls.to_csv(class_path, index=False)
        print(f"Wrote SSI classification: {class_path}")
    print("Done.")


if __name__ == "__main__":
    main()
