from __future__ import annotations

import numpy as np


def spatial_selectivity_index(rate_x: np.ndarray, occupancy_x: np.ndarray) -> float:
    """
    Compute Spatial Selectivity Index (SSI) for one map:

        SSI = sum_i p_i * (lambda_i / lambda_bar) * log2(lambda_i / lambda_bar)

    where:
    - p_i is normalized occupancy in bin i
    - lambda_i is response in bin i
    - lambda_bar = sum_i p_i * lambda_i

    Returns NaN if SSI is undefined (e.g. no valid occupancy or non-positive mean response).
    """
    r = np.asarray(rate_x, dtype=np.float64).ravel()
    p = np.asarray(occupancy_x, dtype=np.float64).ravel()
    if r.size != p.size:
        raise ValueError(f"rate/occupancy size mismatch: {r.size} vs {p.size}")

    valid = np.isfinite(r) & np.isfinite(p) & (p > 0)
    if not np.any(valid):
        return float("nan")

    rv = r[valid]
    pv = p[valid]
    ps = float(np.sum(pv))
    if ps <= 0:
        return float("nan")
    pv = pv / ps

    lam_bar = float(np.sum(pv * rv))
    if not np.isfinite(lam_bar) or lam_bar <= 0:
        return float("nan")

    ratio = rv / lam_bar
    good = np.isfinite(ratio) & (ratio > 0)
    if not np.any(good):
        return 0.0
    ssi = np.sum(pv[good] * ratio[good] * np.log2(ratio[good]))
    return float(ssi)


def spatial_selectivity_index_reps(rate_x_rep: np.ndarray, occupancy_x: np.ndarray) -> np.ndarray:
    """
    Vectorized SSI for many repetitions.

    Parameters
    ----------
    rate_x_rep
        Array with shape (n_bins, n_rep).
    occupancy_x
        Occupancy vector with shape (n_bins,).
    """
    A = np.asarray(rate_x_rep, dtype=np.float64)
    if A.ndim != 2:
        raise ValueError(f"rate_x_rep must be 2D (n_bins, n_rep), got ndim={A.ndim}")
    n_bins, n_rep = A.shape
    p = np.asarray(occupancy_x, dtype=np.float64).ravel()
    if p.size != n_bins:
        raise ValueError(f"occupancy size mismatch: expected {n_bins}, got {p.size}")

    out = np.full(n_rep, np.nan, dtype=np.float64)
    for r in range(n_rep):
        out[r] = spatial_selectivity_index(A[:, r], p)
    return out


def compute_condition_occupancy(
    dwell_tx_x: np.ndarray,
    idcond_t: np.ndarray,
    nb_cond: int,
) -> np.ndarray:
    """
    Build normalized occupancy per condition from trial dwell.

    Returns
    -------
    np.ndarray
        Shape (nb_cond, n_bins), each row normalized to sum=1 when possible.
    """
    dwell = np.asarray(dwell_tx_x, dtype=np.float64)
    if dwell.ndim != 2:
        raise ValueError(f"dwell_tx_x must be 2D (n_trials, n_bins), got ndim={dwell.ndim}")
    n_trials, n_bins = dwell.shape
    idc = np.asarray(idcond_t, dtype=np.int64).ravel()
    if idc.size != n_trials:
        raise ValueError(f"idcond_t length {idc.size} must match n_trials {n_trials}")

    occ = np.full((nb_cond, n_bins), np.nan, dtype=np.float64)
    for c in range(1, nb_cond + 1):
        idx = idc == c
        if not np.any(idx):
            continue
        occ_c = np.nansum(dwell[idx, :], axis=0)
        occ_c = np.where(np.isfinite(occ_c) & (occ_c > 0), occ_c, 0.0)
        s = float(np.sum(occ_c))
        if s <= 0:
            continue
        occ[c - 1, :] = occ_c / s
    return occ


def ssi_observed_by_condition(
    fr_s_cx_x: np.ndarray,
    occupancy_cx: np.ndarray,
) -> np.ndarray:
    """
    Compute observed SSI per condition.
    """
    fr = np.asarray(fr_s_cx_x, dtype=np.float64)
    occ = np.asarray(occupancy_cx, dtype=np.float64)
    if fr.ndim != 2 or occ.ndim != 2:
        raise ValueError("fr_s_cx_x and occupancy_cx must both be 2D (n_cond, n_bins)")
    if fr.shape != occ.shape:
        raise ValueError(f"shape mismatch: fr={fr.shape} vs occ={occ.shape}")

    n_cond = fr.shape[0]
    out = np.full(n_cond, np.nan, dtype=np.float64)
    for c in range(n_cond):
        out[c] = spatial_selectivity_index(fr[c], occ[c])
    return out


def ssi_null_by_condition(
    fr_s_txrep: np.ndarray,
    idcond_t: np.ndarray,
    nb_cond: int,
    occupancy_cx: np.ndarray,
) -> np.ndarray:
    """
    Compute null SSI distribution per condition.

    Parameters
    ----------
    fr_s_txrep
        Null maps with shape (n_trials, n_bins, n_rep).
    """
    null = np.asarray(fr_s_txrep, dtype=np.float64)
    if null.ndim != 3:
        raise ValueError(f"fr_s_txrep must be 3D (n_trials, n_bins, n_rep), got ndim={null.ndim}")
    n_trials, n_bins, n_rep = null.shape
    idc = np.asarray(idcond_t, dtype=np.int64).ravel()
    if idc.size != n_trials:
        raise ValueError(f"idcond_t length {idc.size} must match n_trials {n_trials}")
    occ = np.asarray(occupancy_cx, dtype=np.float64)
    if occ.shape != (nb_cond, n_bins):
        raise ValueError(f"occupancy_cx shape {occ.shape} must be {(nb_cond, n_bins)}")

    out = np.full((nb_cond, n_rep), np.nan, dtype=np.float64)
    for c in range(1, nb_cond + 1):
        idx = idc == c
        if not np.any(idx):
            continue
        mean_null = np.nanmean(null[idx, :, :], axis=0)  # (n_bins, n_rep)
        out[c - 1, :] = spatial_selectivity_index_reps(mean_null, occ[c - 1, :])
    return out


def empirical_pvalue_from_null(
    observed: np.ndarray,
    null_dist: np.ndarray,
    *,
    alternative: str = "greater",
    add_one: bool = True,
) -> np.ndarray:
    """
    Empirical p-values from null distribution.

    Parameters
    ----------
    observed
        Shape (n_cond,).
    null_dist
        Shape (n_cond, n_rep).
    alternative
        "greater" (default), "less", or "two-sided".
    add_one
        If True, use (k+1)/(n+1) correction.
    """
    obs = np.asarray(observed, dtype=np.float64).ravel()
    null = np.asarray(null_dist, dtype=np.float64)
    if null.ndim != 2:
        raise ValueError(f"null_dist must be 2D (n_cond, n_rep), got ndim={null.ndim}")
    n_cond, n_rep = null.shape
    if obs.size != n_cond:
        raise ValueError(f"observed length {obs.size} must match null_dist n_cond {n_cond}")

    out = np.full(n_cond, np.nan, dtype=np.float64)
    for c in range(n_cond):
        o = obs[c]
        nd = null[c]
        nd = nd[np.isfinite(nd)]
        if not np.isfinite(o) or nd.size == 0:
            continue

        if alternative == "greater":
            k = int(np.sum(nd >= o))
        elif alternative == "less":
            k = int(np.sum(nd <= o))
        elif alternative == "two-sided":
            c0 = float(np.nanmedian(nd))
            k = int(np.sum(np.abs(nd - c0) >= abs(o - c0)))
        else:
            raise ValueError("alternative must be one of {'greater','less','two-sided'}")

        n = int(nd.size)
        if add_one:
            out[c] = (k + 1.0) / (n + 1.0)
        else:
            out[c] = k / max(1.0, float(n))
    return out

