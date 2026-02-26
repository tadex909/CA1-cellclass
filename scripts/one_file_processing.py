from __future__ import annotations
import argparse
import json
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
import sys
from typing import Dict, Any, Tuple
import numpy as np
import pandas as pd
THIS_DIR = Path(__file__).resolve().parent
root = THIS_DIR
while root != root.parent and not (root / "src").is_dir():
    root = root.parent

src_dir = root / "src"
if not src_dir.is_dir():
    raise RuntimeError(f"Could not find 'src' dir starting from {THIS_DIR}")

sys.path.insert(0, str(src_dir))

from cellclass.processing import compute_acg, smooth_acg_boxcar, smooth_acg_gaussian
from cellclass.features import (
    compute_cv2,
    refractory_ms_derivative_sd,
    burst_index_from_acg,
    acg_peak_latency_ms,
    waveform_features_from_bestswaveforms,
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


def ensure_dirs(processed_root: Path, mouse: str) -> Dict[str, Path]:
    base = processed_root / mouse
    d = {
        "features": base / "features",
        "acg": base / "acg",
        "qc": base / "qc",
    }
    for p in d.values():
        p.mkdir(parents=True, exist_ok=True)
    return d


def safe_get(z: np.lib.npyio.NpzFile, key: str) -> np.ndarray:
    if key not in z:
        raise KeyError(f"Missing key '{key}' in npz. Available keys: {list(z.keys())[:20]} ...")
    return z[key]


def safe_get_any(z: np.lib.npyio.NpzFile, keys: list[str]) -> np.ndarray:
    for key in keys:
        if key in z:
            return z[key]
    raise KeyError(f"Missing keys {keys!r} in npz. Available keys: {list(z.keys())[:20]} ...")


def normalize_per_cell_field(arr: np.ndarray, n_cells: int, *, name: str) -> np.ndarray:
    """
    Ensure a per-cell field is 1D and aligned to cell_ids length.
    """
    out = np.asarray(arr)
    out = np.squeeze(out)
    if out.ndim != 1:
        raise ValueError(f"{name} must be 1D after squeeze; got shape {arr.shape} -> {out.shape}")
    if out.size != n_cells:
        raise ValueError(f"{name} length mismatch: expected {n_cells}, got {out.size}")
    if out.dtype.kind == "S":
        out = np.char.decode(out, "utf-8", errors="ignore")
    return out


def normalize_per_cell_matrix_field(arr: np.ndarray, n_cells: int, *, name: str) -> np.ndarray:
    """
    Ensure a per-cell matrix-like field is shaped as (n_cells, n_tasks).
    Accepts (n_cells, n_tasks), (n_tasks, n_cells), or object-per-cell vectors.
    """
    out = np.asarray(arr)
    out = np.squeeze(out)

    if out.ndim == 1 and out.size == n_cells and out.dtype == object:
        rows = [np.asarray(v).ravel() for v in out]
        if not rows:
            raise ValueError(f"{name} is empty")
        max_tasks = max(r.size for r in rows)
        if max_tasks == 0:
            out = np.zeros((n_cells, 1), dtype=np.float64)
        else:
            out_pad = np.full((n_cells, max_tasks), np.nan, dtype=np.float64)
            for i, r in enumerate(rows):
                if r.size > 0:
                    out_pad[i, : r.size] = r
            out = out_pad
    elif out.ndim == 1:
        if out.size == n_cells:
            out = out.reshape(n_cells, 1)
        elif out.size % n_cells == 0:
            out = out.reshape(n_cells, out.size // n_cells)
        else:
            raise ValueError(
                f"{name} length mismatch: expected {n_cells} or a multiple, got {out.size}"
            )
    elif out.ndim == 2:
        if out.shape[0] == n_cells:
            pass
        elif out.shape[1] == n_cells:
            out = out.T
        elif out.size % n_cells == 0:
            out = out.reshape(n_cells, out.size // n_cells)
        else:
            raise ValueError(
                f"{name} shape mismatch: expected one axis={n_cells}, got {out.shape}"
            )
    else:
        raise ValueError(f"{name} must be <=2D after squeeze; got shape {arr.shape} -> {out.shape}")

    return np.asarray(out, dtype=np.float64)


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
    acg_smooth_ms: float = 2.0,
    theta_acg_bin_ms: float = 10.0,
    theta_acg_window_ms: float = 700.0,
    theta_acg_smooth_ms: float | None = None,
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
    dirs = ensure_dirs(processed_root, session_meta["mouse"])

    # Load interim npz
    z = np.load(npz_path, allow_pickle=False)

    # Required fields (names as you described)
    spike_times_s = safe_get(z, "allcel__time_spk").astype(np.float64)      # seconds
    spike_cluster_ids = safe_get(z, "allcel__id_spk").astype(np.int64)      # cluster id per spike
    cell_ids = safe_get(z, "allcel__id_cel").astype(np.int64)              # list of cell cluster ids
    type_u_raw = safe_get(z, "allcel__type_u")                               # previous classification per cell
    burst_u_raw = safe_get_any(z, ["allcel__burst_u", "burst_u"])
    fr_u_raw = safe_get_any(z, ["allcel__fr_u", "fr_u"])
    duration_u_raw = safe_get_any(z, ["allcel__duration_u", "duration_u"])
    asymmetry_u_raw = safe_get_any(z, ["allcel__asymmetry_u", "asymmetry_u"])
    pyr_sm_uc_raw = safe_get_any(z, ["allcel__pyr_sm_uc", "pyr_sm_uc"])
    itn_sm_uc_raw = safe_get_any(z, ["allcel__itn_sm_uc", "itn_sm_uc"])

    bestswaveforms = safe_get(z, "allcel__bestswaveforms").astype(np.float64)  # (51, 100, n_cells)
    # optional
    # bestwaveform = z.get("allcel__bestwaveform", None)
    type_u = normalize_per_cell_field(type_u_raw, cell_ids.size, name="allcel__type_u")
    burst_u = normalize_per_cell_field(burst_u_raw, cell_ids.size, name="allcel__burst_u").astype(np.float64)
    fr_u = normalize_per_cell_field(fr_u_raw, cell_ids.size, name="allcel__fr_u").astype(np.float64)
    duration_u = normalize_per_cell_field(duration_u_raw, cell_ids.size, name="allcel__duration_u").astype(np.float64)
    asymmetry_u = normalize_per_cell_field(
        asymmetry_u_raw,
        cell_ids.size,
        name="allcel__asymmetry_u",
    ).astype(np.float64)
    pyr_sm_uc = normalize_per_cell_matrix_field(pyr_sm_uc_raw, cell_ids.size, name="allcel__pyr_sm_uc")
    itn_sm_uc = normalize_per_cell_matrix_field(itn_sm_uc_raw, cell_ids.size, name="allcel__itn_sm_uc")

    n_tasks = int(max(pyr_sm_uc.shape[1], itn_sm_uc.shape[1]))

    def _pad_tasks(arr: np.ndarray, n_cols: int) -> np.ndarray:
        if arr.shape[1] == n_cols:
            return arr
        out = np.full((arr.shape[0], n_cols), np.nan, dtype=np.float64)
        out[:, : arr.shape[1]] = arr
        return out

    pyr_sm_uc = _pad_tasks(pyr_sm_uc, n_tasks)
    itn_sm_uc = _pad_tasks(itn_sm_uc, n_tasks)

    type_u_num = pd.to_numeric(pd.Series(type_u), errors="coerce").to_numpy(dtype=np.float64)
    sm_u = np.where(type_u_num.reshape(-1, 1) == 1.0, pyr_sm_uc, itn_sm_uc)
    sm_u_any = np.nansum(sm_u > 0, axis=1) > 0

    # -------------------------
    # Core per-cell counts
    # -------------------------
    # n_spikes per cell
    n_spikes = np.array([(spike_cluster_ids == cid).sum() for cid in cell_ids], dtype=np.int64)

    #Average firing rate of the cell in the session is defined as the total number of spikes divided by the duration of the recording session. The duration is calculated as the difference between the maximum and minimum spike times. 
    #If there are no spikes, we set the duration to NaN to avoid division by zero.
    
    if spike_times_s.size > 0:
        t0 = float(np.nanmin(spike_times_s))
        t1 = float(np.nanmax(spike_times_s))
        duration_s = max(t1 - t0, 1e-9)
    else:
        duration_s = np.nan

    fr_hz_session = n_spikes / duration_s
    # Use the dataset-provided firing rate for compatibility with legacy features.
    fr_hz = fr_u

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
    acg_counts_smooth = smooth_acg_gaussian(acg_counts, sigma_ms= acg_smooth_ms, bin_ms=acg_bin_ms, axis=0) 

    theta_sigma_ms = acg_smooth_ms if theta_acg_smooth_ms is None else theta_acg_smooth_ms
    theta_acg_res = compute_acg(
        spike_times_s=spike_times_s,
        spike_cluster_ids=spike_cluster_ids,
        cell_ids=cell_ids,
        bin_ms=theta_acg_bin_ms,
        window_ms=theta_acg_window_ms,
        normalize="count",
    )
    theta_acg_counts = theta_acg_res.acg
    theta_lags_ms = theta_acg_res.bin_centers_ms
    theta_acg_counts_smooth = smooth_acg_gaussian(
        theta_acg_counts,
        sigma_ms=theta_sigma_ms,
        bin_ms=theta_acg_bin_ms,
        axis=0,
    )

    # Refractory (Royer derivative method) — you can choose raw or smoothed
    refr = refractory_ms_derivative_sd(acg_counts, lags_ms)
    # Burst index: use smoothed counts by default (more stable)
    burst = burst_index_from_acg(acg_counts_smooth, lags_ms)

    # Peak latency as “rise time” proxy
    #REMEMBER THAT IT'S GOING TO DEPEND ON THE SMOOTHING, SO YOU MIGHT WANT TO TUNE acg_smooth_ms
    peak_lat_ms, peak_val = acg_peak_latency_ms(acg_counts_smooth, lags_ms, search_ms=(1.0, acg_window_ms))

    # -------------------------
    # ISI / CV2
    # -------------------------
    _, cv2 = compute_cv2(spike_times_s, spike_cluster_ids, cell_ids)

    # -------------------------
    # Waveform features (Vinca)
    # -------------------------
    wf_feat = waveform_features_from_bestswaveforms(
    bestswaveforms=bestswaveforms,
    fs_hz=waveform_fs_hz,
    trim=waveform_trim,
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
        "allcel__type_u": type_u,
        "allcel__burst_u": burst_u,
        "allcel__fr_u": fr_u,
        "allcel__duration_u": duration_u,
        "allcel__asymmetry_u": asymmetry_u,
        "allcel__sm_u_any": sm_u_any,

        "n_spikes": n_spikes,
        "fr_hz": fr_hz,
        "fr_hz_session": fr_hz_session,
        "cv2": cv2,

        "refractory_ms_center": refr.refractory_ms_center,
        "refractory_ms_edge": refr.refractory_ms_edge,
        "refractory_peak_ms_center": refr.peak_ms_center,
        "refractory_peak_ms_edge": refr.peak_ms_edge,
        "refractory_bin_ms": refr.bin_ms,
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

    for k in range(sm_u.shape[1]):
        df[f"allcel__sm_u_task{k+1}"] = sm_u[:, k]

    # -------------------------
    # Save outputs
    # -------------------------
    session_id = session_meta["session_id"]

    features_path = dirs["features"] / f"{session_id}_features.parquet"
    acg_path = dirs["acg"] / f"{session_id}_acg_counts_bin{int(acg_bin_ms)}_win{int(acg_window_ms)}.npz"
    theta_acg_path = dirs["acg"] / f"{session_id}_acg_counts_bin{int(theta_acg_bin_ms)}_win{int(theta_acg_window_ms)}.npz"
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
    np.savez_compressed(
        theta_acg_path,
        lags_ms=theta_lags_ms,
        cell_ids=cell_ids,
        acg_counts=theta_acg_counts,
        acg_counts_smooth_boxcar_ms=theta_sigma_ms,
        acg_counts_smooth=theta_acg_counts_smooth,
        acg_bin_ms=theta_acg_bin_ms,
        acg_window_ms=theta_acg_window_ms,
    )

    manifest: Dict[str, Any] = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "source_npz": str(npz_path),
        "session": session_meta,
        "params": {
            "acg_bin_ms": acg_bin_ms,
            "acg_window_ms": acg_window_ms,
            "acg_smooth_ms": acg_smooth_ms,
            "theta_acg_bin_ms": theta_acg_bin_ms,
            "theta_acg_window_ms": theta_acg_window_ms,
            "theta_acg_smooth_ms": theta_sigma_ms,
            "waveform_fs_hz": waveform_fs_hz,
            "waveform_trim": waveform_trim,
            "qc_min_spikes": qc_min_spikes,
        },
        "outputs": {
            "features_parquet": str(features_path),
            "acg_npz": str(acg_path),
            "theta_acg_npz": str(theta_acg_path),
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
    ap.add_argument("--acg_smooth_ms", type=float, default=2.0)
    ap.add_argument("--theta_acg_bin_ms", type=float, default=10.0)
    ap.add_argument("--theta_acg_window_ms", type=float, default=700.0)
    ap.add_argument("--theta_acg_smooth_ms", type=float, default=None)

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
        theta_acg_bin_ms=args.theta_acg_bin_ms,
        theta_acg_window_ms=args.theta_acg_window_ms,
        theta_acg_smooth_ms=args.theta_acg_smooth_ms,
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
