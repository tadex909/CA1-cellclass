from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, Tuple
import numpy as np
from typing import Tuple

@dataclass
class RefractoryResult:
    refractory_ms: np.ndarray     # (n_cells,)
    baseline: np.ndarray          # (n_cells,)
    threshold: np.ndarray         # (n_cells,)
    valid: np.ndarray             # (n_cells,) bool


def refractory_ms_from_acg(
    acg: np.ndarray,
    bin_centers_ms: np.ndarray,
    *,
    baseline_ms: Tuple[float, float] = (20.0, 50.0),
    frac_of_baseline: float = 0.2,
    min_baseline: float = 1e-9,
) -> RefractoryResult:
    """
    Compute refractory recovery time from an autocorrelogram.

    Parameters
    ----------
    acg : (n_lags, n_cells) array
        Autocorrelogram values (counts, per_spike, or rate_hz).
        Must include both negative and positive lags and a 0-lag bin.
    bin_centers_ms : (n_lags,) array
        Lag bin centers in ms, matching acg rows.
    baseline_ms : (lo, hi)
        Range of positive lags (ms) used to estimate baseline.
    frac_of_baseline : float
        Threshold as fraction of baseline.
    min_baseline : float
        If baseline < min_baseline, mark invalid (insufficient activity).

    Returns
    -------
    RefractoryResult
    """
    acg = np.asarray(acg, dtype=np.float64)
    lags = np.asarray(bin_centers_ms, dtype=np.float64)
    assert acg.shape[0] == lags.shape[0]

    # Select positive lags only (exclude 0)
    pos_mask = lags > 0
    lags_pos = lags[pos_mask]
    acg_pos = acg[pos_mask, :]  # (n_pos, n_cells)

    # Baseline mask within positive lags
    b0, b1 = baseline_ms
    base_mask = (lags_pos >= b0) & (lags_pos <= b1)
    if not np.any(base_mask):
        raise ValueError("baseline_ms range does not overlap positive lags in ACG.")

    baseline = np.nanmean(acg_pos[base_mask, :], axis=0)  # (n_cells,)
    valid = np.isfinite(baseline) & (baseline >= min_baseline)

    threshold = baseline * frac_of_baseline

    refractory = np.full(acg.shape[1], np.nan, dtype=np.float64)

    # For valid cells: first lag where ACG crosses threshold
    for j in np.where(valid)[0]:
        above = acg_pos[:, j] >= threshold[j]
        if np.any(above):
            idx = np.argmax(above)  # first True
            refractory[j] = lags_pos[idx]
        else:
            refractory[j] = np.nan  # never recovers within window

    return RefractoryResult(
        refractory_ms=refractory,
        baseline=baseline,
        threshold=threshold,
        valid=valid,
    )

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np


@dataclass
class CV2Result:
    cv2_mean: np.ndarray      # (n_cells,)
    cv2_values: list          # list of arrays, each (n_isi-1,) per cell
    cell_ids: np.ndarray      # (n_cells,)


def compute_cv2(
    spike_times_s: np.ndarray,
    spike_cluster_ids: np.ndarray,
    cell_ids: Optional[np.ndarray] = None,
    *,
    min_spikes: int = 3,
    eps: float = 1e-12,
    return_values: bool = False,
) -> Tuple[np.ndarray, np.ndarray] | CV2Result:
    """
    Compute CV2 per cell.

    Parameters
    ----------
    spike_times_s : (N,) float
        Spike times in seconds.
    spike_cluster_ids : (N,) int
        Cluster id per spike.
    cell_ids : (C,) int, optional
        Which cluster ids to compute CV2 for. If None, uses sorted unique ids.
    min_spikes : int
        Need at least 3 spikes to have 2 ISIs and 1 CV2 value.
    eps : float
        Numerical stabilizer for division.
    return_values : bool
        If True, also return per-interval CV2 arrays per cell.

    Returns
    -------
    If return_values=False:
        (cell_ids, cv2_mean) where cv2_mean is (C,) with NaN if insufficient spikes.
    If return_values=True:
        CV2Result with cv2_values list (one array per cell).
    """
    spike_times_s = np.asarray(spike_times_s, dtype=np.float64).reshape(-1)
    spike_cluster_ids = np.asarray(spike_cluster_ids).reshape(-1)
    assert spike_times_s.shape[0] == spike_cluster_ids.shape[0]

    if cell_ids is None:
        cell_ids = np.unique(spike_cluster_ids)
    cell_ids = np.asarray(cell_ids).reshape(-1)

    cv2_mean = np.full(cell_ids.size, np.nan, dtype=np.float64)
    cv2_values = [] if return_values else None

    for j, cid in enumerate(cell_ids):
        t = spike_times_s[spike_cluster_ids == cid]
        if t.size < min_spikes:
            if return_values:
                cv2_values.append(np.array([], dtype=np.float64))
            continue

        t = np.sort(t)
        isi = np.diff(t)
        # Need consecutive isi pairs: (isi[i], isi[i+1])
        if isi.size < 2:
            if return_values:
                cv2_values.append(np.array([], dtype=np.float64))
            continue

        a = isi[:-1]
        b = isi[1:]
        cv2_i = 2.0 * np.abs(b - a) / (a + b + eps)

        cv2_mean[j] = np.nanmean(cv2_i)
        if return_values:
            cv2_values.append(cv2_i)

    if return_values:
        return CV2Result(cv2_mean=cv2_mean, cv2_values=cv2_values, cell_ids=cell_ids)

    return cell_ids, cv2_mean

def burst_index_from_acg(
    acg: np.ndarray,
    bin_centers_ms: np.ndarray,
    *,
    early_ms: Tuple[float, float] = (0.0, 10.0),
    late_ms: Tuple[float, float] = (40.0, 50.0),
    use_positive_only: bool = True,
    eps: float = 1e-12,
) -> np.ndarray:
    """
    Compute burst index from ACG using the user's definition.

    a = sum(ACG in early window)
    b = sum(ACG in late window)
    burst = (a-b)/a if a>b else (a-b)/b

    Returns (n_cells,) with NaN if windows are empty or both a and b ~ 0.
    """
    acg = np.asarray(acg, dtype=np.float64)
    lags = np.asarray(bin_centers_ms, dtype=np.float64)
    assert acg.shape[0] == lags.shape[0]

    if use_positive_only:
        mask = lags >= 0  # include 0; your "0:10" includes near-zero
        lags = lags[mask]
        acg = acg[mask, :]

    e0, e1 = early_ms
    l0, l1 = late_ms

    early_mask = (lags >= e0) & (lags < e1)
    late_mask = (lags >= l0) & (lags < l1)

    if not np.any(early_mask):
        raise ValueError("early_ms window does not overlap available lags.")
    if not np.any(late_mask):
        raise ValueError("late_ms window does not overlap available lags.")

    a = np.nansum(acg[early_mask, :], axis=0)
    b = np.nansum(acg[late_mask, :], axis=0)

    out = np.full(acg.shape[1], np.nan, dtype=np.float64)

    # Avoid division by ~0
    denom_a = np.abs(a) > eps
    denom_b = np.abs(b) > eps

    # where both are zero-ish, keep NaN
    usable = denom_a | denom_b
    idx = np.where(usable)[0]

    for j in idx:
        if a[j] > b[j] and denom_a[j]:
            out[j] = (a[j] - b[j]) / a[j]
        elif denom_b[j]:
            out[j] = (a[j] - b[j]) / b[j]
        else:
            out[j] = np.nan

    return out

@dataclass
class RefractoryDerivResult:
    refractory_ms: np.ndarray   # (n_cells,)
    peak_ms: np.ndarray         # (n_cells,)
    deriv_sd: np.ndarray        # (n_cells,)
    valid: np.ndarray           # (n_cells,) bool


def refractory_ms_derivative_sd(
    acg_counts: np.ndarray,
    bin_centers_ms: np.ndarray,
    *,
    search_peak_ms: Tuple[float, float] = (1.0, 50.0),
    require_peak_ge: int = 1,
) -> RefractoryDerivResult:
    """
    Royer/Buzsáki-style refractory period from ACG counts with 1 ms bins.

    Steps (per cell):
      1) take positive-lag ACG including 0
      2) find the peak time within search_peak_ms
      3) compute derivative d[k] = acg[k] - acg[k-1] from 0 to peak
      4) compute sd = std(d)
      5) refractory = first bin where d[k] > sd

    Notes:
      - Works best with bin_ms = 1.0
      - If peak is too small or undefined -> returns NaN, valid=False
    """
    acg_counts = np.asarray(acg_counts)
    lags = np.asarray(bin_centers_ms, dtype=np.float64)
    assert acg_counts.shape[0] == lags.shape[0]

    if acg_counts.dtype.kind not in ("i", "u", "f"):
        raise ValueError("acg_counts must be numeric.")

    # Use positive lags including 0
    pos_mask = lags >= 0
    lags_pos = lags[pos_mask]
    acg_pos = np.asarray(acg_counts[pos_mask, :], dtype=np.float64)  # (n_pos, n_cells)

    # sanity: ensure 0 ms exists as first bin
    # If your bin centers are 0.5, 1.5, ... you need to adapt. With our compute_acg it’s 0 at center.
    if lags_pos.size == 0 or lags_pos[0] != 0.0:
        # We can still proceed, but interpretation changes.
        # Better to enforce our compute_acg convention.
        raise ValueError("Expected bin_centers_ms to include 0.0 at the first non-negative bin.")

    lo, hi = search_peak_ms
    peak_search = (lags_pos >= lo) & (lags_pos <= hi)
    if not np.any(peak_search):
        raise ValueError("search_peak_ms does not overlap available positive lags.")

    n_cells = acg_pos.shape[1]
    refractory = np.full(n_cells, np.nan, dtype=np.float64)
    peak_ms = np.full(n_cells, np.nan, dtype=np.float64)
    deriv_sd = np.full(n_cells, np.nan, dtype=np.float64)
    valid = np.zeros(n_cells, dtype=bool)

    # indices in the positive-lag array
    idx_candidates = np.where(peak_search)[0]

    for j in range(n_cells):
        y = acg_pos[:, j]

        # find peak within the search window
        ywin = y[idx_candidates]
        if not np.any(np.isfinite(ywin)):
            continue

        peak_rel = int(np.nanargmax(ywin))
        peak_idx = idx_candidates[peak_rel]
        peak_val = y[peak_idx]

        if peak_val < require_peak_ge:
            # too few counts / too flat -> unreliable
            continue

        # derivative from 0 up to peak: d[1..peak_idx] (since d[k]=y[k]-y[k-1])
        if peak_idx < 1:
            continue
        d = np.diff(y[: peak_idx + 1])

        sd = float(np.nanstd(d))
        deriv_sd[j] = sd
        peak_ms[j] = lags_pos[peak_idx]

        if not np.isfinite(sd) or sd <= 0:
            continue

        # first k where derivative exceeds 1 sd
        above = d > sd
        if np.any(above):
            k = int(np.argmax(above)) + 1  # +1 because d is between bins (k-1 -> k)
            refractory[j] = lags_pos[k]
            valid[j] = True

    return RefractoryDerivResult(
        refractory_ms=refractory,
        peak_ms=peak_ms,
        deriv_sd=deriv_sd,
        valid=valid,
    )

@dataclass
class WaveformFeatResult:
    width_ms: np.ndarray         # (n_cells,)
    asym: np.ndarray             # (n_cells,)
    pt_ratio: np.ndarray         # (n_cells,)
    valid: np.ndarray            # (n_cells,) bool


def _trimmed_mean(x: np.ndarray, trim: int) -> float:
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan
    x = np.sort(x)
    if 2 * trim >= x.size:
        return np.nan
    return float(np.mean(x[trim : x.size - trim]))


def waveform_features_from_bests(
    bestswaveforms: np.ndarray,
    fs_hz: float,
    *,
    peak_search_ms: float = 1.6,
    trim: int = 5,
    eps: float = 1e-12,
) -> WaveformFeatResult:
    """
    Compute waveform width (trough-to-peak latency) and asymmetry per cell,
    using a trimmed mean across waveforms.

    Parameters
    ----------
    bestswaveforms : array (t_step, n_waveforms, n_cells)
        Waveforms on the best channel, multiple sampled spikes per cell.
    fs_hz : float
        Sampling frequency in Hz (e.g., 25000).
    peak_search_ms : float
        Search window after trough to find the peak.
    trim : int
        Number of extreme waveforms to discard at each tail when averaging.
    """
    W = np.asarray(bestswaveforms, dtype=np.float64)
    if W.ndim != 3:
        raise ValueError(f"Expected bestswaveforms ndim=3, got shape {W.shape}")

    t_step, n_wf, n_cells = W.shape
    peak_search_samp = int(round((peak_search_ms / 1000.0) * fs_hz))
    peak_search_samp = max(1, peak_search_samp)

    width_ms = np.full(n_cells, np.nan, dtype=np.float64)
    asym = np.full(n_cells, np.nan, dtype=np.float64)
    pt_ratio = np.full(n_cells, np.nan, dtype=np.float64)
    valid = np.zeros(n_cells, dtype=bool)

    for c in range(n_cells):
        widths = np.full(n_wf, np.nan, dtype=np.float64)
        asyms = np.full(n_wf, np.nan, dtype=np.float64)
        ptrs = np.full(n_wf, np.nan, dtype=np.float64)

        for k in range(n_wf):
            w = W[:, k, c]
            w = w - np.mean(w[:5])

            if not np.any(np.isfinite(w)):
                continue

            tmin = int(np.nanargmin(w))
            # search peak AFTER trough within window
            t1 = min(t_step, tmin + peak_search_samp + 1)
            if t1 <= tmin + 1:
                continue

            post = w[tmin:t1]
            tmax_rel = int(np.nanargmax(post))
            tmax = tmin + tmax_rel

            A_trough = w[tmin]
            A_peak = w[tmax]

            # width in ms
            widths[k] = (tmax - tmin) * 1000.0 / fs_hz

            # peak-to-trough ratio
            ptrs[k] = (A_peak + eps) / (abs(A_trough) + eps)

            # asymmetry in [-1,1]
            asyms[k] = (A_peak - abs(A_trough)) / (A_peak + abs(A_trough) + eps)

        w_mean = _trimmed_mean(widths, trim=trim)
        a_mean = _trimmed_mean(asyms, trim=trim)
        p_mean = _trimmed_mean(ptrs, trim=trim)

        width_ms[c] = w_mean
        asym[c] = a_mean
        pt_ratio[c] = p_mean
        valid[c] = np.isfinite(w_mean)

    return WaveformFeatResult(width_ms=width_ms, asym=asym, pt_ratio=pt_ratio, valid=valid)


@dataclass(frozen=True)
class SpikeShape:
    baseline_min: float
    baseline_idx: int

    pre_peak_val: float
    pre_peak_idx: int

    trough_val: float
    trough_idx: int

    post_peak_val: float
    post_peak_idx: int

    fs_hz: float

    @property
    def pre_peak_t(self) -> float:
        return (self.pre_peak_idx + 1) / self.fs_hz

    @property
    def trough_t(self) -> float:
        return (self.trough_idx + 1) / self.fs_hz

    @property
    def post_peak_t(self) -> float:
        return (self.post_peak_idx + 1) / self.fs_hz

    @property
    def peak_to_trough_ms(self) -> float:
        return (self.trough_t - self.pre_peak_t) * 1000.0

    @property
    def trough_to_rebound_ms(self) -> float:
        return (self.post_peak_t - self.trough_t) * 1000.0

    @property
    def asymmetry(self) -> float:
        a1 = self.pre_peak_val - self.baseline_min
        a2 = self.post_peak_val - self.baseline_min
        return (a1 - a2) / (a1 + a2 + 1e-12)
    
def extract_spike_shape(waveform: np.ndarray, fs_hz: float) -> SpikeShape:
    w = np.asarray(waveform, dtype=np.float64)

    trough_idx = int(np.argmin(w))
    trough_val = w[trough_idx]

    # pre-trough peak
    pre_seg = w[: trough_idx + 1]
    pre_peak_idx = int(np.argmax(pre_seg))
    pre_peak_val = w[pre_peak_idx]

    # post-trough rebound: first crossing ≥ pre_peak_val
    post_seg = w[trough_idx:]
    above = np.flatnonzero(post_seg >= pre_peak_val)
    if above.size > 0:
        post_peak_idx = trough_idx + int(above[0])
    else:
        post_peak_idx = trough_idx + int(np.argmax(post_seg))
    post_peak_val = w[post_peak_idx]

    # baseline minimum before pre_peak
    base_seg = w[: pre_peak_idx + 1]
    baseline_idx = int(np.argmin(base_seg))
    baseline_min = w[baseline_idx]

    return SpikeShape(
        baseline_min=baseline_min,
        baseline_idx=baseline_idx,
        pre_peak_val=pre_peak_val,
        pre_peak_idx=pre_peak_idx,
        trough_val=trough_val,
        trough_idx=trough_idx,
        post_peak_val=post_peak_val,
        post_peak_idx=post_peak_idx,
        fs_hz=fs_hz,
    )

