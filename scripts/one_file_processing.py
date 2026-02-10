#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Tuple

import numpy as np
import pandas as pd

# Adjust these imports to your package structure
from cellclass.processing import compute_acg, smooth_acg_boxcar
from cellclass.features import (
    compute_cv2,
    refractory_ms_derivative_sd,
    burst_index_from_acg,
    acg_peak_latency_ms,
    waveformfeature_vinca,  # <-- or whatever you name it
)

# -------------------------
# Helpers
# -------------------------

def parse_session_id(npz_path: Path) -> Dict[str, str]:
    """
    Example filename:
      VS57_2022-12-18_18-51-04_Ratemap_final_thèse_allcel.npz
    We'll take the first 3 underscore-separated fields as:
      mouse, date, time
    and build session_id = mouse_date_time
    """
    stem = npz_path.stem
    parts = stem.split("_")
    out = {
        "mouse": parts[0] if len(parts) > 0 else "unknown",
        "date": parts[1] if len(parts) > 1 else "unknown",
        "time": parts[2] if len(parts) > 2 else "unknown",
    }
    out["session_id"] = f"{out['mouse']}_{out['date']}_{out['time']}"
    return out


def ensure_dirs(processed_root: Path) -> Dict[str, Path]:
    d = {
        "features": processed_root / "features",
        "acg": processed_root / "acg",
        "qc": processed_root / "qc",
    }
    for p in d.values():
        p.mkdir(parents=True, exist_ok=True)
    return d


def safe_get(z: np.lib.npyio.NpzFile, key: str) -> np.ndarray:
    if key not in z:
        raise KeyError(f"Missing key '{key}' in npz. Available keys: {list(z.keys())[:20]} ...")
    return z[key]


# -------------------------
# Main extraction
# -------------------------

def extract_one(
    npz_path: Path,
    processed_root: Path,
    *,
    # ACG params
    acg_bin_ms: float = 1.0,
    acg_window_ms: float = 50.0,
    acg_smooth_ms: float = 5.0,
    # Waveform params
    waveform_fs_hz: float = 25000.0,
    waveform_trim: int = 5,
    # QC params
    qc_min_spikes: int = 200,
) -> Tuple[Path, Path, Path]:
    """
    Returns paths to (features_parquet, acg_npz, manifest_json).
    """
    session_meta = parse_session_id(npz_path)
    dirs = ensure_dirs(processed_root)

    # Load interim npz
    z = np.load(npz_path, allow_pickle=False)

    # Required fields (names as you described)
    spike_times_s = safe_get(z, "allcel__time_spk").astype(np.float64)      # seconds
    spike_cluster_ids = safe_get(z, "allcel__id_spk").astype(np.int64)      # cluster id per spike
    cell_ids = safe_get(z, "allcel__id_cel").astype(np.int64)              # list of cell cluster ids

    bestswaveforms = safe_get(z, "allcel__bestswaveforms").astype(np.float64)  # (51, 100, n_cells)
    # optional
    # bestwaveform = z.get("allcel__bestwaveform", None)

    # -------------------------
    # Core per-cell counts
    # -------------------------
    # n_spikes per cell
    n_spikes = np.array([(spike_cluster_ids == cid).sum() for cid in cell_ids], dtype=np.int64)

    # firing rate requires session duration; if you don't have it, estimate from spike times
    # (This is a reasonable first pass; later you may prefer e.g. ephys recording duration.)
    if spike_times_s.size > 0:
        t0 = float(np.nanmin(spike_times_s))
        t1 = float(np.nanmax(spike_times_s))
        duration_s = max(t1 - t0, 1e-9)
    else:
        duration_s = np.nan

    fr_hz = n_spikes / duration_s

    # -------------------------
    # ACG (counts) + smoothing
    # -------------------------
    acg_res = compute_acg(
        spike_times_s=spike_times_s,
        spike_cluster_ids=spike_cluster_ids,
        cell_ids=cell_ids,
        bin_ms=acg_bin_ms,
        window_ms=acg_window_ms,
        normalize="count",
    )

    acg_counts = acg_res.acg               # (n_lags, n_cells)
    lags_ms = acg_res.bin_centers_ms       # (n_lags,)

    acg_counts_smooth = smooth_acg_boxcar(
        acg_counts, window_ms=acg_smooth_ms, bin_ms=acg_bin_ms, axis=0
    )

    # Refractory (Royer derivative method) — you can choose raw or smoothed
    refr = refractory_ms_derivative_sd(acg_counts, lags_ms)
    # Burst index: use smoothed counts by default (more stable)
    burst = burst_index_from_acg(acg_counts_smooth, lags_ms)

    # Peak latency as “rise time” proxy
    peak_lat_ms, peak_val = acg_peak_latency_ms(acg_counts_smooth, lags_ms, search_ms=(1.0, acg_window_ms))

    # -------------------------
    # ISI / CV2
    # -------------------------
    _, cv2 = compute_cv2(spike_times_s, spike_cluster_ids, cell_ids)

    # -------------------------
    # Waveform features (Vinca)
    # -------------------------
    wf_feat = waveformfeature_vinca(
        bestswaveforms=bestswaveforms,
        fs_hz=waveform_fs_hz,
        trim=waveform_trim,
        baseline_subtract=False,  # match lab default
    )

    # -------------------------
    # QC flags
    # -------------------------
    qc_min_spikes_flag = n_spikes >= qc_min_spikes
    qc_refractory = refr.valid
    qc_waveform = wf_feat.valid
    # You can add more later:
    # qc_acg_peak = peak_val >= some_threshold

    # -------------------------
    # Build feature table
    # -------------------------
    df = pd.DataFrame({
        "session_id": session_meta["session_id"],
        "mouse": session_meta["mouse"],
        "date": session_meta["date"],
        "time": session_meta["time"],
        "cell_id": cell_ids,

        "n_spikes": n_spikes,
        "fr_hz": fr_hz,
        "cv2": cv2,

        "refractory_ms_center": refr.refractory_ms_center,
        "refractory_ms_edge": refr.refractory_ms_edge,
        "refractory_peak_ms": refr.peak_ms_center,

        "burst_index": burst,

        "acg_peak_latency_ms": peak_lat_ms,
        "acg_peak_value": peak_val,

        "spk_duration_ms": wf_feat.spk_duration_ms,
        "spk_peaktrough_ms": wf_feat.spk_peaktrough_ms,
        "spk_asymmetry": wf_feat.spk_asymmetry,

        "qc_min_spikes": qc_min_spikes_flag,
        "qc_refractory": qc_refractory,
        "qc_waveform": qc_waveform,
    })

    # -------------------------
    # Save outputs
    # -------------------------
    session_id = session_meta["session_id"]

    features_path = dirs["features"] / f"{session_id}_features.parquet"
    acg_path = dirs["acg"] / f"{session_id}_acg_counts_bin{int(acg_bin_ms)}_win{int(acg_window_ms)}.npz"
    manifest_path = dirs["qc"] / f"{session_id}_manifest.json"

    df.to_parquet(features_path, index=False)

    np.savez_compressed(
        acg_path,
        lags_ms=lags_ms,
        cell_ids=cell_ids,
        acg_counts=acg_counts,
        acg_counts_smooth_boxcar_ms=acg_smooth_ms,
        acg_counts_smooth=acg_counts_smooth,
        acg_bin_ms=acg_bin_ms,
        acg_window_ms=acg_window_ms,
    )

    manifest: Dict[str, Any] = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "source_npz": str(npz_path),
        "session": session_meta,
        "params": {
            "acg_bin_ms": acg_bin_ms,
            "acg_window_ms": acg_window_ms,
            "acg_smooth_ms": acg_smooth_ms,
            "waveform_fs_hz": waveform_fs_hz,
            "waveform_trim": waveform_trim,
            "qc_min_spikes": qc_min_spikes,
        },
        "outputs": {
            "features_parquet": str(features_path),
            "acg_npz": str(acg_path),
        },
        "counts": {
            "n_cells": int(len(cell_ids)),
            "n_spikes_total": int(len(spike_times_s)),
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2))

    return features_path, acg_path, manifest_path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", type=str, required=True, help="Path to one interim .npz file")
    ap.add_argument("--processed_root", type=str, default="data/processed", help="Output processed root folder")

    ap.add_argument("--acg_bin_ms", type=float, default=1.0)
    ap.add_argument("--acg_window_ms", type=float, default=50.0)
    ap.add_argument("--acg_smooth_ms", type=float, default=5.0)

    ap.add_argument("--waveform_fs_hz", type=float, default=25000.0)
    ap.add_argument("--waveform_trim", type=int, default=5)

    ap.add_argument("--qc_min_spikes", type=int, default=200)

    args = ap.parse_args()

    features_path, acg_path, manifest_path = extract_one(
        npz_path=Path(args.npz),
        processed_root=Path(args.processed_root),
        acg_bin_ms=args.acg_bin_ms,
        acg_window_ms=args.acg_window_ms,
        acg_smooth_ms=args.acg_smooth_ms,
        waveform_fs_hz=args.waveform_fs_hz,
        waveform_trim=args.waveform_trim,
        qc_min_spikes=args.qc_min_spikes,
    )

    print("Wrote:")
    print("  ", features_path)
    print("  ", acg_path)
    print("  ", manifest_path)


if __name__ == "__main__":
    main()
