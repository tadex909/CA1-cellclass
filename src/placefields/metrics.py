from __future__ import annotations

from typing import Optional

import numpy as np


def make_bin_edges(track_min_cm: float, track_max_cm: float, bin_size_cm: float) -> np.ndarray:
    if track_max_cm <= track_min_cm:
        raise ValueError("track_max_cm must be > track_min_cm")
    if bin_size_cm <= 0:
        raise ValueError("bin_size_cm must be > 0")

    span = track_max_cm - track_min_cm
    n_bins = max(1, int(np.ceil(span / bin_size_cm)))
    return np.linspace(track_min_cm, track_max_cm, n_bins + 1, dtype=np.float64)


def compute_occupancy_s(position_cm: np.ndarray, dt_s: np.ndarray, bin_edges_cm: np.ndarray) -> np.ndarray:
    pos = np.asarray(position_cm, dtype=np.float64).ravel()
    dt = np.asarray(dt_s, dtype=np.float64).ravel()
    if pos.size != dt.size:
        raise ValueError(f"position_cm and dt_s size mismatch: {pos.size} vs {dt.size}")

    valid = np.isfinite(pos) & np.isfinite(dt) & (dt > 0)
    if not np.any(valid):
        return np.zeros(bin_edges_cm.size - 1, dtype=np.float64)
    occ, _ = np.histogram(pos[valid], bins=bin_edges_cm, weights=dt[valid])
    return occ.astype(np.float64)


def compute_spike_counts(spike_position_cm: np.ndarray, bin_edges_cm: np.ndarray) -> np.ndarray:
    spk_pos = np.asarray(spike_position_cm, dtype=np.float64).ravel()
    valid = np.isfinite(spk_pos)
    if not np.any(valid):
        return np.zeros(bin_edges_cm.size - 1, dtype=np.float64)
    counts, _ = np.histogram(spk_pos[valid], bins=bin_edges_cm)
    return counts.astype(np.float64)


def compute_rate_map_hz(
    spike_counts: np.ndarray,
    occupancy_s: np.ndarray,
    *,
    min_occupancy_s: float,
) -> np.ndarray:
    c = np.asarray(spike_counts, dtype=np.float64).ravel()
    occ = np.asarray(occupancy_s, dtype=np.float64).ravel()
    if c.size != occ.size:
        raise ValueError(f"spike_counts and occupancy_s size mismatch: {c.size} vs {occ.size}")

    rate = np.full_like(occ, np.nan, dtype=np.float64)
    valid = occ >= float(min_occupancy_s)
    rate[valid] = c[valid] / occ[valid]
    return rate


def _gaussian_kernel_1d(sigma_bins: float) -> np.ndarray:
    if sigma_bins <= 0:
        return np.array([1.0], dtype=np.float64)
    radius = max(1, int(np.ceil(4.0 * sigma_bins)))
    x = np.arange(-radius, radius + 1, dtype=np.float64)
    kernel = np.exp(-(x**2) / (2.0 * sigma_bins**2))
    kernel /= kernel.sum()
    return kernel


def smooth_rate_map_hz(rate_map_hz: np.ndarray, sigma_bins: float) -> np.ndarray:
    r = np.asarray(rate_map_hz, dtype=np.float64).ravel()
    if sigma_bins <= 0 or r.size == 0:
        return r.copy()

    finite = np.isfinite(r)
    if not np.any(finite):
        return r.copy()

    filled = r.copy()
    filled[~finite] = 0.0
    weights = finite.astype(np.float64)
    k = _gaussian_kernel_1d(float(sigma_bins))
    num = np.convolve(filled, k, mode="same")
    den = np.convolve(weights, k, mode="same")

    out = np.full_like(r, np.nan, dtype=np.float64)
    valid = den > 1e-12
    out[valid] = num[valid] / den[valid]
    return out


def compute_spatial_information_bits_per_spike(
    rate_map_hz: np.ndarray,
    occupancy_s: np.ndarray,
) -> Optional[float]:
    rate = np.asarray(rate_map_hz, dtype=np.float64).ravel()
    occ = np.asarray(occupancy_s, dtype=np.float64).ravel()
    if rate.size != occ.size:
        raise ValueError(f"rate_map_hz and occupancy_s size mismatch: {rate.size} vs {occ.size}")

    valid = np.isfinite(rate) & np.isfinite(occ) & (occ > 0) & (rate >= 0)
    if not np.any(valid):
        return None

    p_i = occ[valid] / occ[valid].sum()
    r_i = rate[valid]
    r_bar = float(np.sum(p_i * r_i))
    if r_bar <= 0:
        return None

    ratio = r_i / r_bar
    keep = ratio > 0
    if not np.any(keep):
        return 0.0
    info = float(np.sum(p_i[keep] * ratio[keep] * np.log2(ratio[keep])))
    return info


def detect_place_fields(
    rate_map_hz: np.ndarray,
    *,
    threshold_ratio: float,
    min_bins: int,
) -> list[tuple[int, int]]:
    if not (0 < threshold_ratio <= 1):
        raise ValueError("threshold_ratio must be in (0, 1]")
    if min_bins < 1:
        raise ValueError("min_bins must be >= 1")

    r = np.asarray(rate_map_hz, dtype=np.float64).ravel()
    if r.size == 0:
        return []

    finite = np.isfinite(r)
    if not np.any(finite):
        return []

    peak = float(np.nanmax(r))
    if not np.isfinite(peak) or peak <= 0:
        return []

    thr = peak * float(threshold_ratio)
    mask = np.isfinite(r) & (r >= thr)
    fields: list[tuple[int, int]] = []

    start: Optional[int] = None
    for i, on in enumerate(mask):
        if on and start is None:
            start = i
        if (not on) and (start is not None):
            end = i - 1
            if (end - start + 1) >= min_bins:
                fields.append((start, end))
            start = None
    if start is not None:
        end = mask.size - 1
        if (end - start + 1) >= min_bins:
            fields.append((start, end))

    return fields
