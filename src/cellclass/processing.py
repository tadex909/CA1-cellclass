from __future__ import annotations
#The following script takes script from data/interim and processes it to create the final NPZ files in data/processed. It also adds metadata about the source .mat file and the processing time.
from dataclasses import dataclass
from typing import Literal, Optional, Tuple

import numpy as np


NormalizeMode = Literal["count", "per_spike", "rate_hz"]


@dataclass
class ACGResult:
    bin_centers_ms: np.ndarray          # shape (2*n_bins + 1,)
    acg: np.ndarray                     # shape (2*n_bins + 1, n_cells)
    cell_ids: np.ndarray                # shape (n_cells,)
    normalize: str
    bin_ms: float
    window_ms: float


@dataclass
class CCGResult:
    bin_centers_ms: np.ndarray          # shape (2*n_bins + 1,)
    ccg: np.ndarray                     # shape (2*n_bins + 1, n_cells, n_cells)
    cell_ids: np.ndarray                # shape (n_cells,)
    normalize: str
    bin_ms: float
    window_ms: float


def _acg_one_unit(times_s: np.ndarray, bin_s: float, window_s: float) -> np.ndarray:
    """
    Positive-lag histogram for one unit.
    Returns hist_pos with length n_bins, where bin k counts dt in [k*bin_s, (k+1)*bin_s).
    """
    times_s = np.asarray(times_s, dtype=np.float64)
    if times_s.size < 2:
        n_bins = int(np.ceil(window_s / bin_s))
        return np.zeros(n_bins, dtype=np.int64)

    times_s = np.sort(times_s)
    n_bins = int(np.ceil(window_s / bin_s))
    hist_pos = np.zeros(n_bins, dtype=np.int64)

    # For each spike i, count spikes in (t_i, t_i + window]
    for i in range(times_s.size - 1):
        t0 = times_s[i]
        stop = np.searchsorted(times_s, t0 + window_s, side="right")
        if stop <= i + 1:
            continue
        dt = times_s[i + 1 : stop] - t0  # strictly positive
        # Bin indices
        b = np.floor(dt / bin_s).astype(np.int64)
        # Safety: dt==window can land in bin n_bins (rare float edge); clip
        b = b[(b >= 0) & (b < n_bins)]
        if b.size:
            hist_pos += np.bincount(b, minlength=n_bins).astype(np.int64)

    return hist_pos


def _ccg_one_pair(
    times_a_s: np.ndarray,
    times_b_s: np.ndarray,
    bin_s: float,
    window_s: float,
) -> np.ndarray:
    """
    Positive-lag histogram from A to B.
    Returns hist_pos with length n_bins, where bin k counts dt in [k*bin_s, (k+1)*bin_s).
    """
    times_a_s = np.asarray(times_a_s, dtype=np.float64)
    times_b_s = np.asarray(times_b_s, dtype=np.float64)
    n_bins = int(np.ceil(window_s / bin_s))
    if times_a_s.size == 0 or times_b_s.size == 0:
        return np.zeros(n_bins, dtype=np.int64)

    times_a_s = np.sort(times_a_s)
    times_b_s = np.sort(times_b_s)
    hist_pos = np.zeros(n_bins, dtype=np.int64)

    # For each spike in A, count spikes in B in (t_a, t_a + window]
    for i in range(times_a_s.size):
        t0 = times_a_s[i]
        start = np.searchsorted(times_b_s, t0, side="right")
        stop = np.searchsorted(times_b_s, t0 + window_s, side="right")
        if stop <= start:
            continue
        dt = times_b_s[start:stop] - t0
        b = np.floor(dt / bin_s).astype(np.int64)
        b = b[(b >= 0) & (b < n_bins)]
        if b.size:
            hist_pos += np.bincount(b, minlength=n_bins).astype(np.int64)

    return hist_pos


def compute_acg(
    spike_times_s: np.ndarray,
    spike_cluster_ids: np.ndarray,
    cell_ids: Optional[np.ndarray] = None,
    *,
    bin_ms: float = 1.0,
    window_ms: float = 50.0,
    normalize: NormalizeMode = "rate_hz",
) -> ACGResult:
    """
    Compute autocorrelograms for each cell id.

    Parameters
    ----------
    spike_times_s : (N,) float
        Spike times in seconds.
    spike_cluster_ids : (N,) int
        Cluster id per spike.
    cell_ids : (C,) int, optional
        Which cluster ids to compute ACG for. If None, uses sorted unique cluster ids.
    bin_ms : float
        Bin width in milliseconds.
    window_ms : float
        Half-window in milliseconds. Output covers [-window_ms, +window_ms].
    normalize : {"count","per_spike","rate_hz"}
        Normalization for the ACG values.

    Returns
    -------
    ACGResult
    """
    spike_times_s = np.asarray(spike_times_s, dtype=np.float64).reshape(-1)
    spike_cluster_ids = np.asarray(spike_cluster_ids).reshape(-1)
    assert spike_times_s.shape[0] == spike_cluster_ids.shape[0], "times and ids must match length"

    if cell_ids is None:
        cell_ids = np.unique(spike_cluster_ids)
    cell_ids = np.asarray(cell_ids).reshape(-1)

    bin_s = bin_ms / 1000.0
    window_s = window_ms / 1000.0
    n_bins = int(np.ceil(window_s / bin_s))

    # centers: negative, zero, positive
    pos_centers = (np.arange(n_bins) + 0.5) * bin_ms
    bin_centers_ms = np.concatenate([-pos_centers[::-1], [0.0], pos_centers])

    acg = np.zeros((2 * n_bins + 1, cell_ids.size), dtype=np.float64)

    for j, cid in enumerate(cell_ids):
        mask = spike_cluster_ids == cid
        t = spike_times_s[mask]
        hist_pos = _acg_one_unit(t, bin_s=bin_s, window_s=window_s)  # int64

        # Mirror to negative lags; set zero lag to 0
        hist_full = np.concatenate([hist_pos[::-1], np.array([0], dtype=np.int64), hist_pos])

        if normalize == "count":
            acg[:, j] = hist_full
        elif normalize == "per_spike":
            denom = max(t.size, 1)
            acg[:, j] = hist_full / denom
        elif normalize == "rate_hz":
            denom = max(t.size, 1) * bin_s
            acg[:, j] = hist_full / denom
        else:
            raise ValueError(f"Unknown normalize={normalize}")

    return ACGResult(
        bin_centers_ms=bin_centers_ms,
        acg=acg,
        cell_ids=cell_ids,
        normalize=normalize,
        bin_ms=bin_ms,
        window_ms=window_ms,
    )


def compute_ccg(
    spike_times_s: np.ndarray,
    spike_cluster_ids: np.ndarray,
    cell_ids: Optional[np.ndarray] = None,
    *,
    bin_ms: float = 1.0,
    window_ms: float = 50.0,
    normalize: NormalizeMode = "rate_hz",
    include_auto: bool = False,
) -> CCGResult:
    """
    Compute cross-correlograms for all cell pairs.

    Parameters
    ----------
    spike_times_s : (N,) float
        Spike times in seconds.
    spike_cluster_ids : (N,) int
        Cluster id per spike.
    cell_ids : (C,) int, optional
        Which cluster ids to compute CCG for. If None, uses sorted unique cluster ids.
    bin_ms : float
        Bin width in milliseconds.
    window_ms : float
        Half-window in milliseconds. Output covers [-window_ms, +window_ms].
    normalize : {"count","per_spike","rate_hz"}
        Normalization for the CCG values (per spike of the reference cell).
    include_auto : bool
        If True, fills diagonal with ACG counts. If False, diagonal is zeros.

    Returns
    -------
    CCGResult
    """
    spike_times_s = np.asarray(spike_times_s, dtype=np.float64).reshape(-1)
    spike_cluster_ids = np.asarray(spike_cluster_ids).reshape(-1)
    assert spike_times_s.shape[0] == spike_cluster_ids.shape[0], "times and ids must match length"

    if cell_ids is None:
        cell_ids = np.unique(spike_cluster_ids)
    cell_ids = np.asarray(cell_ids).reshape(-1)

    bin_s = bin_ms / 1000.0
    window_s = window_ms / 1000.0
    n_bins = int(np.ceil(window_s / bin_s))

    pos_centers = (np.arange(n_bins) + 0.5) * bin_ms
    bin_centers_ms = np.concatenate([-pos_centers[::-1], [0.0], pos_centers])

    n_cells = cell_ids.size
    ccg = np.zeros((2 * n_bins + 1, n_cells, n_cells), dtype=np.float64)

    # Pre-split spike times by cell for reuse
    times_by_cell = [spike_times_s[spike_cluster_ids == cid] for cid in cell_ids]

    for i in range(n_cells):
        t_i = times_by_cell[i]
        for j in range(n_cells):
            if i == j and not include_auto:
                continue
            t_j = times_by_cell[j]
            if i == j and include_auto:
                hist_pos = _acg_one_unit(t_i, bin_s=bin_s, window_s=window_s)
                hist_full = np.concatenate([hist_pos[::-1], np.array([0], dtype=np.int64), hist_pos])
            else:
                hist_pos_ij = _ccg_one_pair(t_i, t_j, bin_s=bin_s, window_s=window_s)
                hist_pos_ji = _ccg_one_pair(t_j, t_i, bin_s=bin_s, window_s=window_s)
                hist_full = np.concatenate([hist_pos_ji[::-1], np.array([0], dtype=np.int64), hist_pos_ij])

            if normalize == "count":
                ccg[:, i, j] = hist_full
            elif normalize == "per_spike":
                denom = max(t_i.size, 1)
                ccg[:, i, j] = hist_full / denom
            elif normalize == "rate_hz":
                denom = max(t_i.size, 1) * bin_s
                ccg[:, i, j] = hist_full / denom
            else:
                raise ValueError(f"Unknown normalize={normalize}")

    return CCGResult(
        bin_centers_ms=bin_centers_ms,
        ccg=ccg,
        cell_ids=cell_ids,
        normalize=normalize,
        bin_ms=bin_ms,
        window_ms=window_ms,
    )


def compute_ccg_from_npz(
    npz_path: str,
    *,
    bin_ms: float = 1.0,
    window_ms: float = 50.0,
    normalize: NormalizeMode = "rate_hz",
    include_auto: bool = False,
    time_key: str = "allcel__time_spk",
    id_key: str = "allcel__id_spk",
    cell_ids_key: str = "allcel__id_cel",
) -> CCGResult:
    """
    Load an interim NPZ and compute CCGs between neurons.
    """
    z = np.load(npz_path, allow_pickle=False)
    spike_times_s = z[time_key].astype(np.float64)
    spike_cluster_ids = z[id_key].astype(np.int64)
    if cell_ids_key in z:
        cell_ids = z[cell_ids_key].astype(np.int64)
    else:
        cell_ids = None

    return compute_ccg(
        spike_times_s=spike_times_s,
        spike_cluster_ids=spike_cluster_ids,
        cell_ids=cell_ids,
        bin_ms=bin_ms,
        window_ms=window_ms,
        normalize=normalize,
        include_auto=include_auto,
    )

#-----smoothing of the ACG-----
def smooth_acg_boxcar(
    acg: np.ndarray,
    *,
    window_ms: float = 5.0,
    bin_ms: float = 1.0,
    axis: int = 0,
    mode: str = "reflect",
) -> np.ndarray:
    """
    Smooth an ACG with a boxcar (moving average) window.

    Parameters
    ----------
    acg : ndarray
        ACG array, e.g. (n_lags, n_cells) or (n_lags,).
    window_ms : float
        Smoothing window in ms.
    bin_ms : float
        Bin size in ms.
    axis : int
        Axis corresponding to lag bins.
    mode : str
        Padding mode for np.pad: "reflect", "edge", "constant", ...

    Returns
    -------
    smoothed : ndarray
        Same shape as acg.
    """
    x = np.asarray(acg, dtype=np.float64)
    w = int(round(window_ms / bin_ms))
    w = max(1, w)
    if w == 1:
        return x.copy()

    # make w odd so smoothing is centered
    if w % 2 == 0:
        w += 1
    half = w // 2

    # move axis to front
    x0 = np.moveaxis(x, axis, 0)
    pad_width = [(half, half)] + [(0, 0)] * (x0.ndim - 1)
    xp = np.pad(x0, pad_width=pad_width, mode=mode)

    kernel = np.ones(w, dtype=np.float64) / w
    # convolve along first axis
    out = np.empty_like(x0)
    for idx in np.ndindex(x0.shape[1:]):
        out[(slice(None),) + idx] = np.convolve(xp[(slice(None),) + idx], kernel, mode="valid")

    return np.moveaxis(out, 0, axis)

def smooth_acg_gaussian(
    acg: np.ndarray,
    *,
    sigma_ms: float = 2.0,
    bin_ms: float = 1.0,
    axis: int = 0,
    truncate: float = 3.0,
    mode: str = "reflect",
) -> np.ndarray:
    """
    Smooth an ACG with a Gaussian kernel.

    sigma_ms ~ 2ms gives ~5ms-ish smoothing (since kernel spans ~6*sigma).
    """
    x = np.asarray(acg, dtype=np.float64)
    sigma_bins = sigma_ms / bin_ms
    if sigma_bins <= 0:
        return x.copy()

    radius = int(np.ceil(truncate * sigma_bins))
    kx = np.arange(-radius, radius + 1)
    kernel = np.exp(-0.5 * (kx / sigma_bins) ** 2)
    kernel /= kernel.sum()

    x0 = np.moveaxis(x, axis, 0)
    pad_width = [(radius, radius)] + [(0, 0)] * (x0.ndim - 1)
    xp = np.pad(x0, pad_width=pad_width, mode=mode)

    out = np.empty_like(x0)
    for idx in np.ndindex(x0.shape[1:]):
        out[(slice(None),) + idx] = np.convolve(xp[(slice(None),) + idx], kernel, mode="valid")

    return np.moveaxis(out, 0, axis)
