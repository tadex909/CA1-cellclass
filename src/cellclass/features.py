from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, Tuple
import numpy as np

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
    refractory_ms_center: np.ndarray  # (n_cells,)
    refractory_ms_edge: np.ndarray    # (n_cells,)
    peak_ms_center: np.ndarray        # (n_cells,)
    peak_ms_edge: np.ndarray          # (n_cells,)
    deriv_sd: np.ndarray              # (n_cells,)
    valid: np.ndarray                 # (n_cells,) bool
    bin_ms: float


def refractory_ms_derivative_sd(
    acg_counts: np.ndarray,
    bin_centers_ms: np.ndarray,
    *,
    search_peak_ms: Tuple[float, float] = (1.0, 50.0),
    require_peak_ge: int = 1,
) -> RefractoryDerivResult:
    """
    Royer/Buzsáki-style refractory period from ACG counts.

    Works with bin centers at 0,1,2,... OR 0.5,1.5,2.5,... etc.
    We return both:
      - center convention: the bin center time
      - edge convention: center - bin_ms/2 (clipped at 0)

    Steps (per cell):
      1) use non-negative lags (>=0)
      2) find ACG peak within search_peak_ms
      3) compute derivative d[k] = y[k] - y[k-1] from start to peak
      4) sd = std(d)
      5) refractory = first bin where d[k] > sd
    """
    acg_counts = np.asarray(acg_counts)
    lags = np.asarray(bin_centers_ms, dtype=np.float64)
    if acg_counts.shape[0] != lags.shape[0]:
        raise ValueError("acg_counts and bin_centers_ms must have same first dimension")

    # Estimate bin size robustly
    if lags.size < 2:
        raise ValueError("Need at least 2 lag bins")
    bin_ms = float(np.median(np.diff(lags)))
    if not np.isfinite(bin_ms) or bin_ms <= 0:
        raise ValueError("Could not infer a valid bin_ms from bin_centers_ms")

    # Use non-negative lags
    pos_mask = lags >= 0
    lags_pos = lags[pos_mask]
    acg_pos = np.asarray(acg_counts[pos_mask, :], dtype=np.float64)  # (n_pos, n_cells)
    if lags_pos.size < 2:
        raise ValueError("Not enough non-negative lag bins to compute derivative")

    lo, hi = search_peak_ms
    peak_search = (lags_pos >= lo) & (lags_pos <= hi)
    idx_candidates = np.where(peak_search)[0]
    if idx_candidates.size == 0:
        raise ValueError("search_peak_ms does not overlap available non-negative lags")

    n_cells = acg_pos.shape[1]
    refractory_center = np.full(n_cells, np.nan, dtype=np.float64)
    peak_center = np.full(n_cells, np.nan, dtype=np.float64)
    deriv_sd = np.full(n_cells, np.nan, dtype=np.float64)
    valid = np.zeros(n_cells, dtype=bool)

    for j in range(n_cells):
        y = acg_pos[:, j]

        ywin = y[idx_candidates]
        if not np.any(np.isfinite(ywin)):
            continue

        peak_rel = int(np.nanargmax(ywin))
        peak_idx = int(idx_candidates[peak_rel])
        peak_val = y[peak_idx]
        if not np.isfinite(peak_val) or peak_val < require_peak_ge:
            continue

        if peak_idx < 1:
            continue

        # derivative from first bin up to peak bin
        d = np.diff(y[: peak_idx + 1])
        sd = float(np.nanstd(d))
        deriv_sd[j] = sd
        peak_center[j] = lags_pos[peak_idx]

        if not np.isfinite(sd) or sd <= 0:
            continue

        above = d > sd
        if np.any(above):
            k = int(np.argmax(above)) + 1  # +1 because d indexes transitions to bin k
            refractory_center[j] = lags_pos[k]
            valid[j] = True

    # Convert center -> edge convention
    refractory_edge = np.maximum(refractory_center - bin_ms / 2.0, 0.0)
    peak_edge = np.maximum(peak_center - bin_ms / 2.0, 0.0)

    return RefractoryDerivResult(
        refractory_ms_center=refractory_center,
        refractory_ms_edge=refractory_edge,
        peak_ms_center=peak_center,
        peak_ms_edge=peak_edge,
        deriv_sd=deriv_sd,
        valid=valid,
        bin_ms=bin_ms,
    )



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

@dataclass
class WaveformFeatResult:
    spk_duration_ms: np.ndarray        # (n_cells,)
    spk_peaktrough_ms: np.ndarray      # (n_cells,)
    spk_asymmetry: np.ndarray          # (n_cells,)
    valid: np.ndarray                  # (n_cells,) bool

@dataclass
class WaveformFeatResult:
    spk_duration_ms: np.ndarray        # (n_cells,)
    spk_peaktrough_ms: np.ndarray      # (n_cells,)
    spk_asymmetry: np.ndarray          # (n_cells,)
    valid: np.ndarray                  # (n_cells,) bool


def _trimmed_mean(x: np.ndarray, trim: int) -> float:
    """Drop `trim` values at each tail after sorting; return mean of remaining."""
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan
    x = np.sort(x)
    if 2 * trim >= x.size:
        return np.nan
    return float(np.mean(x[trim : x.size - trim]))


def waveform_features_from_bestswaveforms(
    bestswaveforms: np.ndarray,
    fs_hz: float,
    *,
    trim: int = 5,
) -> WaveformFeatResult:
    """
    Compute waveform features per cell from bestswaveforms using extract_spike_shape.

    Parameters
    ----------
    bestswaveforms : np.ndarray
        Shape (t_step, n_waveforms, n_cells).
    fs_hz : float
        Waveform sampling rate (e.g., 25000).
    trim : int
        Number of waveforms to drop at each tail for robust averaging
        (Vinca uses 5: mean(sorted_vals(5:end-5))).

    Returns
    -------
    WaveformFeatResult
    """
    W = np.asarray(bestswaveforms, dtype=np.float64)
    if W.ndim != 3:
        raise ValueError(f"Expected bestswaveforms with ndim=3, got shape {W.shape}")

    t_step, n_wf, n_cells = W.shape
    if fs_hz <= 0:
        raise ValueError("fs_hz must be > 0")

    spk_duration_ms = np.full(n_cells, np.nan, dtype=np.float64)
    spk_peaktrough_ms = np.full(n_cells, np.nan, dtype=np.float64)
    spk_asymmetry = np.full(n_cells, np.nan, dtype=np.float64)
    valid = np.zeros(n_cells, dtype=bool)

    for c in range(n_cells):
        dur = np.full(n_wf, np.nan, dtype=np.float64)
        pt = np.full(n_wf, np.nan, dtype=np.float64)
        asym = np.full(n_wf, np.nan, dtype=np.float64)

        for j in range(n_wf):
            w = W[:, j, c]
            if not np.any(np.isfinite(w)):
                continue

            shape = extract_spike_shape(w, fs_hz)

            # These names come from the SpikeShape class we defined earlier
            dur[j] = shape.trough_to_rebound_ms
            pt[j] = shape.peak_to_trough_ms
            asym[j] = shape.asymmetry

        spk_duration_ms[c] = _trimmed_mean(dur, trim=trim)
        spk_peaktrough_ms[c] = _trimmed_mean(pt, trim=trim)
        spk_asymmetry[c] = _trimmed_mean(asym, trim=trim)

        valid[c] = np.isfinite(spk_duration_ms[c]) and np.isfinite(spk_peaktrough_ms[c]) and np.isfinite(spk_asymmetry[c])

    return WaveformFeatResult(
        spk_duration_ms=spk_duration_ms,
        spk_peaktrough_ms=spk_peaktrough_ms,
        spk_asymmetry=spk_asymmetry,
        valid=valid,
    )

def acg_peak_latency_ms(
    acg_counts: np.ndarray,
    bin_centers_ms: np.ndarray,
    *,
    search_ms: Tuple[float, float] = (1.0, 50.0),
    use_positive_only: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Returns:
      peak_latency_ms: (n_cells,)
      peak_value:      (n_cells,)
    """
    acg = np.asarray(acg_counts, dtype=np.float64)
    lags = np.asarray(bin_centers_ms, dtype=np.float64)
    assert acg.shape[0] == lags.shape[0]

    if use_positive_only:
        m = lags > 0
        lags = lags[m]
        acg = acg[m, :]

    lo, hi = search_ms
    m = (lags >= lo) & (lags <= hi)
    if not np.any(m):
        raise ValueError("search_ms does not overlap available lags.")

    l = lags[m]
    a = acg[m, :]  # (n_bins, n_cells)

    peak_idx = np.nanargmax(a, axis=0)
    peak_latency = l[peak_idx]
    peak_value = a[peak_idx, np.arange(a.shape[1])]

    return peak_latency, peak_value