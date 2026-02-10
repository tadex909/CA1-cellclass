#The following script takes script from data/interim and processes it to create the final NPZ files in data/processed. It also adds metadata about the source .mat file and the processing time.
from __future__ import annotations

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
