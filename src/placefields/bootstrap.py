from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Sequence

import numpy as np

from .metrics import smooth_last_axis
from .trials import TrialInfo

NullMethod = Literal["random", "poisson", "circular_shift"]


@dataclass(frozen=True)
class NullBootstrapConfig:
    """
    Parameters for MATLAB-like null map simulation (`fct_placefield_simtrain`-style).
    """

    method: NullMethod = "random"
    nb_rep: int = 1000
    min_shift_frac: float = 0.10
    min_shift_samples: int = 1
    seed: int | None = None

    def validate(self) -> None:
        if self.method not in {"random", "poisson", "circular_shift"}:
            raise ValueError("method must be one of {'random','poisson','circular_shift'}")
        if self.nb_rep < 1:
            raise ValueError("nb_rep must be >= 1")
        if not np.isfinite(self.min_shift_frac) or self.min_shift_frac < 0:
            raise ValueError("min_shift_frac must be finite and >= 0")
        if self.min_shift_samples < 0:
            raise ValueError("min_shift_samples must be >= 0")


def _bin_indices(x: np.ndarray, edges: np.ndarray) -> np.ndarray:
    idx = np.searchsorted(edges, x, side="right") - 1
    return idx.astype(np.int64, copy=False)


def _draw_nontrivial_shift(n: int, *, cfg: NullBootstrapConfig, rng: np.random.Generator) -> int:
    if n <= 1:
        return 0

    min_shift = max(int(cfg.min_shift_samples), int(np.floor(cfg.min_shift_frac * n)))
    if min_shift <= 0:
        min_shift = 1

    shifts = np.arange(n, dtype=np.int64)
    mask = (shifts >= min_shift) & (shifts <= (n - min_shift))
    valid = shifts[mask]

    if valid.size == 0:
        valid = shifts[1:]
    return int(valid[rng.integers(0, valid.size)])


def _counts_random(
    *,
    nb_rep: int,
    nb_spk: int,
    candidate_bin_idx: np.ndarray,
    n_bins: int,
    rng: np.random.Generator,
) -> np.ndarray:
    out = np.zeros((nb_rep, n_bins), dtype=np.float64)
    n_candidates = int(candidate_bin_idx.size)
    if nb_spk <= 0 or n_candidates <= 0:
        return out

    draw = rng.integers(0, n_candidates, size=(nb_rep, nb_spk))
    cols = candidate_bin_idx[draw.ravel()]
    rows = np.repeat(np.arange(nb_rep, dtype=np.int64), nb_spk)
    np.add.at(out, (rows, cols), 1.0)
    return out


def _counts_poisson(
    *,
    nb_rep: int,
    nb_spk: int,
    candidate_bin_idx: np.ndarray,
    n_bins: int,
    rng: np.random.Generator,
) -> np.ndarray:
    out = np.zeros((nb_rep, n_bins), dtype=np.float64)
    n_candidates = int(candidate_bin_idx.size)
    if nb_spk <= 0 or n_candidates <= 0:
        return out

    lam = float(nb_spk) / float(n_candidates)
    sampled = rng.poisson(lam, size=(nb_rep, n_candidates)).astype(np.float64, copy=False)
    cols = np.tile(candidate_bin_idx, nb_rep)
    rows = np.repeat(np.arange(nb_rep, dtype=np.int64), n_candidates)
    np.add.at(out, (rows, cols), sampled.ravel())
    return out


def _counts_circular_shift(
    *,
    nb_rep: int,
    spike_local_idx: np.ndarray,
    candidate_bin_idx: np.ndarray,
    n_bins: int,
    cfg: NullBootstrapConfig,
    rng: np.random.Generator,
) -> np.ndarray:
    out = np.zeros((nb_rep, n_bins), dtype=np.float64)
    n_candidates = int(candidate_bin_idx.size)
    if spike_local_idx.size == 0 or n_candidates <= 0:
        return out

    for r in range(nb_rep):
        shift = _draw_nontrivial_shift(n_candidates, cfg=cfg, rng=rng)
        shifted = (spike_local_idx + shift) % n_candidates
        bins = candidate_bin_idx[shifted]
        out[r] = np.bincount(bins, minlength=n_bins).astype(np.float64, copy=False)
    return out


def simulate_null_fr_s_txrep(
    *,
    position_x: np.ndarray,
    spike_indices_cell_0b: np.ndarray,
    trials: Sequence[TrialInfo],
    xbin_edges: np.ndarray,
    dwell_s_tx_x: np.ndarray,
    smooth_sigma_bins: float,
    cfg: NullBootstrapConfig,
    valid_sample_mask: np.ndarray | None = None,
) -> np.ndarray:
    """
    Build null smoothed rate maps for one cell, analogous to MATLAB `fct_placefield_simtrain`.

    Parameters
    ----------
    position_x
        Session position vector (1D, behavior frequency).
    spike_indices_cell_0b
        Spike sample indices for one cell (Python 0-based, same frequency as `position_x`).
    trials
        Trial boundaries/labels.
    xbin_edges
        Spatial bin edges used for ratemap construction.
    dwell_s_tx_x
        Smoothed dwell per trial/bin, shape (n_trials, n_bins).
    smooth_sigma_bins
        Gaussian sigma in bins used for both count and FR smoothing.
    cfg
        Null bootstrap configuration (`random`, `poisson`, `circular_shift`).
    valid_sample_mask
        Optional boolean mask over session samples (e.g. speed thresholding). If None, all finite
        samples are considered valid candidates.

    Returns
    -------
    np.ndarray
        `fr_s_txrep` with shape (n_trials, n_bins, nb_rep).
    """

    cfg.validate()
    x = np.asarray(position_x, dtype=np.float64).ravel()
    n_samples = int(x.size)
    if n_samples == 0:
        raise ValueError("position_x is empty")

    spk = np.asarray(spike_indices_cell_0b, dtype=np.int64).ravel()
    if spk.size:
        spk = spk[(spk >= 0) & (spk < n_samples)]

    edges = np.asarray(xbin_edges, dtype=np.float64).ravel()
    if edges.size < 2 or not np.all(np.diff(edges) > 0):
        raise ValueError("xbin_edges must be strictly increasing with at least 2 elements")
    n_bins = int(edges.size - 1)

    n_trials = len(trials)
    dwell_s = np.asarray(dwell_s_tx_x, dtype=np.float64)
    if dwell_s.shape != (n_trials, n_bins):
        raise ValueError(
            f"dwell_s_tx_x shape mismatch: expected {(n_trials, n_bins)}, got {dwell_s.shape}"
        )

    if valid_sample_mask is None:
        sample_valid = np.ones(n_samples, dtype=bool)
    else:
        sample_valid = np.asarray(valid_sample_mask, dtype=bool).ravel()
        if sample_valid.size != n_samples:
            raise ValueError(
                f"valid_sample_mask size mismatch: expected {n_samples}, got {sample_valid.size}"
            )

    sample_valid &= np.isfinite(x)
    sample_valid &= (x >= edges[0]) & (x <= edges[-1])

    out = np.zeros((n_trials, n_bins, cfg.nb_rep), dtype=np.float64)
    rng = np.random.default_rng(cfg.seed)

    for t, tr in enumerate(trials):
        s = max(0, int(tr.start_idx_0b))
        e = min(n_samples, int(tr.stop_idx_0b_exclusive))
        if e <= s:
            continue

        trial_idx = np.arange(s, e, dtype=np.int64)
        idxnonan = trial_idx[sample_valid[s:e]]

        # Match MATLAB behavior: if no candidate sample, run with N=1 and no spikes.
        if idxnonan.size == 0:
            candidate_bin_idx = np.array([], dtype=np.int64)
            nb_spk = 0
        else:
            candidate_bin_idx = _bin_indices(x[idxnonan], edges)
            keep_cand = (candidate_bin_idx >= 0) & (candidate_bin_idx < n_bins)
            idxnonan = idxnonan[keep_cand]
            candidate_bin_idx = candidate_bin_idx[keep_cand]

            in_trial = spk[(spk >= s) & (spk < e)]
            if in_trial.size:
                keep_spk = sample_valid[in_trial]
                spk_valid = in_trial[keep_spk]
            else:
                spk_valid = np.array([], dtype=np.int64)
            nb_spk = int(spk_valid.size)

        if cfg.method == "random":
            nbspk_repx = _counts_random(
                nb_rep=cfg.nb_rep,
                nb_spk=nb_spk,
                candidate_bin_idx=candidate_bin_idx,
                n_bins=n_bins,
                rng=rng,
            )
        elif cfg.method == "poisson":
            nbspk_repx = _counts_poisson(
                nb_rep=cfg.nb_rep,
                nb_spk=nb_spk,
                candidate_bin_idx=candidate_bin_idx,
                n_bins=n_bins,
                rng=rng,
            )
        else:
            in_trial = spk[(spk >= s) & (spk < e)]
            if in_trial.size and idxnonan.size:
                keep_spk = sample_valid[in_trial]
                spk_valid = in_trial[keep_spk]
                local_idx = np.searchsorted(idxnonan, spk_valid)
                good = (local_idx < idxnonan.size) & (idxnonan[local_idx] == spk_valid)
                local_idx = local_idx[good].astype(np.int64, copy=False)
            else:
                local_idx = np.array([], dtype=np.int64)

            nbspk_repx = _counts_circular_shift(
                nb_rep=cfg.nb_rep,
                spike_local_idx=local_idx,
                candidate_bin_idx=candidate_bin_idx,
                n_bins=n_bins,
                cfg=cfg,
                rng=rng,
            )

        nbspk_s_repx = smooth_last_axis(nbspk_repx, smooth_sigma_bins)
        fr_repx = np.divide(
            nbspk_s_repx,
            dwell_s[t][None, :],
            out=np.zeros_like(nbspk_s_repx),
            where=dwell_s[t][None, :] > 0,
        )
        fr_s_repx = smooth_last_axis(fr_repx, smooth_sigma_bins)
        out[t, :, :] = fr_s_repx.T

    return out


def empirical_pval_tx(fr_s_txrep: np.ndarray, fr_s_tx_x: np.ndarray) -> np.ndarray:
    """
    Trial-level empirical p-values:
        p(t,x) = P(null >= observed) estimated over repetitions.
    """
    null = np.asarray(fr_s_txrep, dtype=np.float64)
    obs = np.asarray(fr_s_tx_x, dtype=np.float64)
    if null.ndim != 3:
        raise ValueError(f"fr_s_txrep must be 3D (t,x,rep), got ndim={null.ndim}")
    if obs.shape != null.shape[:2]:
        raise ValueError(f"fr_s_tx_x shape {obs.shape} must match fr_s_txrep[:2] {null.shape[:2]}")
    return np.mean(null >= obs[:, :, None], axis=2)


def empirical_pval_cx(
    fr_s_txrep: np.ndarray,
    fr_s_cx_x: np.ndarray,
    idcond_t: np.ndarray,
    nb_cond: int,
) -> np.ndarray:
    """
    Condition-level empirical p-values, matching MATLAB logic:
    - average null maps over laps in each condition
    - compare to observed condition map
    """
    null = np.asarray(fr_s_txrep, dtype=np.float64)
    obs_c = np.asarray(fr_s_cx_x, dtype=np.float64)
    idc = np.asarray(idcond_t, dtype=np.int64).ravel()

    if null.ndim != 3:
        raise ValueError(f"fr_s_txrep must be 3D (t,x,rep), got ndim={null.ndim}")
    n_trials, n_bins, _ = null.shape
    if idc.size != n_trials:
        raise ValueError(f"idcond_t length {idc.size} must match n_trials {n_trials}")
    if obs_c.shape != (nb_cond, n_bins):
        raise ValueError(f"fr_s_cx_x shape {obs_c.shape} must be {(nb_cond, n_bins)}")

    out = np.full((nb_cond, n_bins), np.nan, dtype=np.float64)
    for c in range(1, nb_cond + 1):
        idx = idc == c
        if not np.any(idx):
            continue
        mean_null = np.nanmean(null[idx, :, :], axis=0)  # (x,rep)
        out[c - 1, :] = np.mean(mean_null >= obs_c[c - 1, :, None], axis=1)
    return out

