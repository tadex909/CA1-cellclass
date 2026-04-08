from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np


def decode_npz_json_scalar(arr: np.ndarray) -> Any:
    return json.loads(arr.tobytes().decode("utf-8", errors="ignore"))


def stitch_trial_series(
    start_1b: np.ndarray,
    stop_1b: np.ndarray,
    series_by_trial: list[Any],
) -> np.ndarray:
    """
    Build a session-level vector from per-trial vectors.

    Placement rule:
    - if len(vals) == span: place on [start-1, stop)
    - if len(vals) == span-2: place on interior [start, stop-1)
    - otherwise: place as much as possible from trial start
    """
    s = np.asarray(start_1b, dtype=np.int64).ravel()
    e = np.asarray(stop_1b, dtype=np.int64).ravel()
    if s.size != e.size:
        raise ValueError(f"start/stop size mismatch: {s.size} vs {e.size}")
    if len(series_by_trial) != s.size:
        raise ValueError(
            f"series_by_trial length mismatch: expected {s.size}, got {len(series_by_trial)}"
        )

    n_samples = int(np.max(e)) if e.size else 0
    out = np.full(n_samples, np.nan, dtype=np.float64)

    for i in range(s.size):
        s0 = int(s[i] - 1)
        e0 = int(e[i])  # exclusive in Python after MATLAB conversion
        span = max(0, e0 - s0)
        vals = np.asarray(series_by_trial[i], dtype=np.float64).ravel()
        if span <= 0 or vals.size == 0:
            continue

        if vals.size == span:
            a = s0
            b = e0
        elif vals.size == (span - 2) and span >= 2:
            a = s0 + 1
            b = e0 - 1
        else:
            a = s0
            b = min(e0, s0 + vals.size)

        if b <= a:
            continue
        n = b - a
        out[a:b] = vals[:n]

    return out


def find_pairs(interim_root: Path, sessions: set[str]) -> list[tuple[str, Path, Path]]:
    out: list[tuple[str, Path, Path]] = []
    for allcel_path in sorted(interim_root.rglob("*_allcel.npz")):
        stem = allcel_path.stem
        if not stem.endswith("_allcel"):
            continue
        session = stem[: -len("_allcel")]
        if sessions and session not in sessions:
            continue
        traj_path = allcel_path.with_name(f"{session}_trajdata.npz")
        if traj_path.exists():
            out.append((session, allcel_path, traj_path))
    return out


def load_allcel_spikes(allcel_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with np.load(allcel_path, allow_pickle=False) as za:
        itime_25k = np.asarray(za["allcel__itime_spk"]).ravel()
        id_spk = np.asarray(za["allcel__id_spk"]).ravel().astype(np.int64)
        id_cel = np.asarray(za["allcel__id_cel"]).ravel().astype(np.int64)
    return itime_25k, id_spk, id_cel


def load_traj_fields(
    traj_path: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[Any], list[Any] | None]:
    with np.load(traj_path, allow_pickle=False) as zt:
        cond = np.asarray(zt["traj__Cond"]).ravel()
        wb = np.asarray(zt["traj__WB"]).ravel()
        start = np.asarray(zt["traj__start"]).ravel()
        stop = np.asarray(zt["traj__stop"]).ravel()
        vr_list = decode_npz_json_scalar(zt["traj__VRtraj__json"])
        # Prefer XSpeed for velocity masking; fallback to Speed for older files.
        if "traj__XSpeed__json" in zt:
            speed_list = decode_npz_json_scalar(zt["traj__XSpeed__json"])
        elif "traj__Speed__json" in zt:
            speed_list = decode_npz_json_scalar(zt["traj__Speed__json"])
        else:
            speed_list = None
    return cond, wb, start, stop, vr_list, speed_list


def save_ratemap_pack(pack: Any, out_path: Path, meta: dict[str, Any]) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {
        "meta_json": np.array(json.dumps(meta), dtype=np.string_),
        "cell_ids": pack.cell_ids,
        "idcond_t": pack.idcond_t,
        "xbin_edges": pack.xbin_edges,
        "xbin_centers": pack.xbin_centers,
        "rmap__nbspk_tx_ux": pack.nbspk_tx_ux,
        "rmap__dwell_tx_x": pack.dwell_tx_x,
        "rmap__fr_tx_ux": pack.fr_tx_ux,
        "rmap__nbspk_s_tx_ux": pack.nbspk_s_tx_ux,
        "rmap__dwell_s_tx_x": pack.dwell_s_tx_x,
        "rmap__fr_s_tx_ux": pack.fr_s_tx_ux,
        "rmap__nbspk_cx_ux": pack.nbspk_cx_ux,
        "rmap__dwell_cx_x": pack.dwell_cx_x,
        "rmap__fr_cx_ux": pack.fr_cx_ux,
        "rmap__nbspk_s_cx_ux": pack.nbspk_s_cx_ux,
        "rmap__dwell_s_cx_x": pack.dwell_s_cx_x,
        "rmap__fr_s_cx_ux": pack.fr_s_cx_ux,
    }
    np.savez_compressed(out_path, **payload)
