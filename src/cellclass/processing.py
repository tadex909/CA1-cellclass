#The following script takes script from data/interim and processes it to create the final NPZ files in data/processed. It also adds metadata about the source .mat file and the processing time.
from __future__ import annotations

from dataclasses import dataclass
import warnings
import importlib
import importlib.util
from typing import Callable, Literal, Optional, Tuple

import numpy as np


NormalizeMode = Literal["count", "per_spike", "rate_hz"]
BackendMode = Literal["auto", "numpy", "numba"]
CCGBackendMode = Literal["numpy", "numba"]


_NUMBA_ACG_KERNEL: Optional[Callable] = None
_NUMBA_CCG_KERNEL: Optional[Callable] = None


def _get_numba_acg_kernel() -> Optional[Callable]:
    """Return a compiled numba kernel for ACG, or None if numba is unavailable."""
    global _NUMBA_ACG_KERNEL
    if _NUMBA_ACG_KERNEL is not None:
        return _NUMBA_ACG_KERNEL

    if importlib.util.find_spec("numba") is None:
        return None

    numba = importlib.import_module("numba")

    @numba.njit(cache=True)
    def _acg_one_unit_numba(
        times_s: np.ndarray,
        bin_s: float,
        window_s: float,
        n_bins: int,
    ) -> np.ndarray:
        n = times_s.shape[0]
        hist_pos = np.zeros(n_bins, dtype=np.int64)
        if n < 2:
            return hist_pos

        for i in range(n - 1):
            t0 = times_s[i]
            j = i + 1
            while j < n:
                dt = times_s[j] - t0
                if dt > window_s:
                    break
                b = int(dt / bin_s)
                if 0 <= b < n_bins:
                    hist_pos[b] += 1
                j += 1
        return hist_pos

    _NUMBA_ACG_KERNEL = _acg_one_unit_numba
    return _NUMBA_ACG_KERNEL


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
    backend: BackendMode = "auto",
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
    backend : {"auto", "numpy", "numba"}
        Backend used for histogram accumulation. "auto" uses numba when available,
        otherwise falls back to numpy.

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

    numba_kernel: Optional[Callable] = None
    if backend not in {"auto", "numpy", "numba"}:
        raise ValueError(f"Unknown backend={backend}")
    if backend in {"auto", "numba"}:
        numba_kernel = _get_numba_acg_kernel()
        if backend == "numba" and numba_kernel is None:
            raise ImportError(
                "backend='numba' requested but numba is not installed. "
                "Install numba or use backend='numpy' (or 'auto')."
            )

    # centers: negative, zero, positive
    pos_centers = (np.arange(n_bins) + 0.5) * bin_ms
    bin_centers_ms = np.concatenate([-pos_centers[::-1], [0.0], pos_centers])

    acg = np.zeros((2 * n_bins + 1, cell_ids.size), dtype=np.float64)

    for j, cid in enumerate(cell_ids):
        mask = spike_cluster_ids == cid
        t = spike_times_s[mask]
        if numba_kernel is not None:
            t_sorted = np.sort(t.astype(np.float64, copy=False))
            hist_pos = numba_kernel(t_sorted, bin_s, window_s, n_bins)  # int64
        else:
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


@dataclass
class CCGResult:
    bin_centers_ms: np.ndarray          # shape (2*n_bins + 1,)
    ccg: np.ndarray                     # shape (2*n_bins + 1, n_ref_cells, n_target_cells)
    ref_cell_ids: np.ndarray            # shape (n_ref_cells,)
    target_cell_ids: np.ndarray         # shape (n_target_cells,)
    normalize: str
    bin_ms: float
    window_ms: float
    backend: str


def _get_numba_ccg_kernel() -> Optional[Callable]:
    """Return a compiled numba kernel for CCG pair histograms, or None if unavailable."""
    global _NUMBA_CCG_KERNEL
    if _NUMBA_CCG_KERNEL is not None:
        return _NUMBA_CCG_KERNEL

    if importlib.util.find_spec("numba") is None:
        return None

    numba = importlib.import_module("numba")

    @numba.njit(cache=True)
    def _ccg_one_pair_numba(
        ref_times: np.ndarray,
        target_times: np.ndarray,
        bin_s: float,
        window_s: float,
        n_full_bins: int,
        remove_zero: bool,
    ) -> np.ndarray:
        hist = np.zeros(n_full_bins, dtype=np.int64)
        if ref_times.size == 0 or target_times.size == 0:
            return hist

        for i in range(ref_times.size):
            t0 = ref_times[i]
            left = np.searchsorted(target_times, t0 - window_s, side="left")
            right = np.searchsorted(target_times, t0 + window_s, side="right")
            for j in range(left, right):
                dt = target_times[j] - t0
                if remove_zero and dt == 0.0:
                    continue
                b = int((dt + window_s) / bin_s)
                if 0 <= b < n_full_bins:
                    hist[b] += 1
        return hist

    _NUMBA_CCG_KERNEL = _ccg_one_pair_numba
    return _NUMBA_CCG_KERNEL


def _ccg_one_pair_numpy(
    ref_times: np.ndarray,
    target_times: np.ndarray,
    *,
    bin_s: float,
    window_s: float,
    n_full_bins: int,
    remove_zero: bool,
) -> np.ndarray:
    """Event-based CCG histogram for one ref/target pair using searchsorted windows."""
    ref_times = np.asarray(ref_times, dtype=np.float64)
    target_times = np.asarray(target_times, dtype=np.float64)
    hist = np.zeros(n_full_bins, dtype=np.int64)

    if ref_times.size == 0 or target_times.size == 0:
        return hist

    ref_times = np.sort(ref_times)
    target_times = np.sort(target_times)

    for t0 in ref_times:
        left = np.searchsorted(target_times, t0 - window_s, side="left")
        right = np.searchsorted(target_times, t0 + window_s, side="right")
        if right <= left:
            continue

        dt = target_times[left:right] - t0
        if remove_zero:
            dt = dt[dt != 0.0]
        if dt.size == 0:
            continue

        bins = np.floor((dt + window_s) / bin_s).astype(np.int64)
        bins = bins[(bins >= 0) & (bins < n_full_bins)]
        if bins.size:
            hist += np.bincount(bins, minlength=n_full_bins).astype(np.int64)

    return hist


def compute_ccg(
    spike_times_s: np.ndarray,
    spike_cluster_ids: np.ndarray,
    ref_cell_ids: Optional[np.ndarray] = None,
    target_cell_ids: Optional[np.ndarray] = None,
    *,
    bin_ms: float = 1.0,
    window_ms: float = 50.0,
    normalize: NormalizeMode = "rate_hz",
    backend: CCGBackendMode = "numpy",
) -> CCGResult:
    """
    Compute cross-correlograms for ref/target cell pairs using event windows and searchsorted.

    Parameters
    ----------
    spike_times_s : (N,) float
        Spike times in seconds.
    spike_cluster_ids : (N,) int
        Cluster id per spike.
    ref_cell_ids : (R,) int, optional
        Reference cluster ids. If None, uses sorted unique cluster ids.
    target_cell_ids : (T,) int, optional
        Target cluster ids. If None, uses sorted unique cluster ids.
    bin_ms : float
        Bin width in milliseconds.
    window_ms : float
        Half-window in milliseconds. Output covers [-window_ms, +window_ms].
    normalize : {"count","per_spike","rate_hz"}
        Normalization for the CCG values based on reference spike count.
    backend : {"numpy", "numba"}
        Histogram backend. If numba is unavailable, automatically falls back to numpy.
    """
    spike_times_s = np.asarray(spike_times_s, dtype=np.float64).reshape(-1)
    spike_cluster_ids = np.asarray(spike_cluster_ids).reshape(-1)
    assert spike_times_s.shape[0] == spike_cluster_ids.shape[0], "times and ids must match length"

    if ref_cell_ids is None:
        ref_cell_ids = np.unique(spike_cluster_ids)
    if target_cell_ids is None:
        target_cell_ids = np.unique(spike_cluster_ids)

    ref_cell_ids = np.asarray(ref_cell_ids).reshape(-1)
    target_cell_ids = np.asarray(target_cell_ids).reshape(-1)

    if backend not in {"numpy", "numba"}:
        raise ValueError(f"Unknown backend={backend}")

    numba_kernel: Optional[Callable] = None
    backend_used = backend
    if backend == "numba":
        numba_kernel = _get_numba_ccg_kernel()
        if numba_kernel is None:
            warnings.warn(
                "backend='numba' requested but numba is not installed; falling back to numpy.",
                RuntimeWarning,
                stacklevel=2,
            )
            backend_used = "numpy"

    bin_s = bin_ms / 1000.0
    window_s = window_ms / 1000.0
    n_side_bins = int(np.ceil(window_s / bin_s))
    n_full_bins = 2 * n_side_bins + 1

    bin_centers_ms = np.arange(-n_side_bins, n_side_bins + 1, dtype=np.float64) * bin_ms
    ccg = np.zeros((n_full_bins, ref_cell_ids.size, target_cell_ids.size), dtype=np.float64)

    spike_times_by_id = {
        cid: np.sort(spike_times_s[spike_cluster_ids == cid].astype(np.float64, copy=False))
        for cid in np.union1d(ref_cell_ids, target_cell_ids)
    }

    for i, ref_cid in enumerate(ref_cell_ids):
        ref_t = spike_times_by_id[ref_cid]
        for j, tgt_cid in enumerate(target_cell_ids):
            tgt_t = spike_times_by_id[tgt_cid]
            remove_zero = ref_cid == tgt_cid

            if numba_kernel is not None:
                hist = numba_kernel(ref_t, tgt_t, bin_s, window_s, n_full_bins, remove_zero)
            else:
                hist = _ccg_one_pair_numpy(
                    ref_t,
                    tgt_t,
                    bin_s=bin_s,
                    window_s=window_s,
                    n_full_bins=n_full_bins,
                    remove_zero=remove_zero,
                )

            if normalize == "count":
                ccg[:, i, j] = hist
            elif normalize == "per_spike":
                ccg[:, i, j] = hist / max(ref_t.size, 1)
            elif normalize == "rate_hz":
                ccg[:, i, j] = hist / (max(ref_t.size, 1) * bin_s)
            else:
                raise ValueError(f"Unknown normalize={normalize}")

    return CCGResult(
        bin_centers_ms=bin_centers_ms,
        ccg=ccg,
        ref_cell_ids=ref_cell_ids,
        target_cell_ids=target_cell_ids,
        normalize=normalize,
        bin_ms=bin_ms,
        window_ms=window_ms,
        backend=backend_used,
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
