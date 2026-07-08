from __future__ import annotations

from typing import Any

import numpy as np

from .decoding import BayesianDecodingResult


def _window_indices(
    result: BayesianDecodingResult,
    window_idx: np.ndarray | list[int] | list[bool] | None,
) -> np.ndarray:
    n_windows = int(np.asarray(result.posterior_wx).shape[0])
    if window_idx is None:
        return np.arange(n_windows, dtype=np.int64)

    idx = np.asarray(window_idx)
    if idx.dtype == bool:
        if idx.ndim != 1 or idx.size != n_windows:
            raise ValueError(
                f"Boolean window_idx must have length {n_windows}, got shape {idx.shape}"
            )
        return np.flatnonzero(idx).astype(np.int64, copy=False)

    if not np.issubdtype(idx.dtype, np.integer):
        raise TypeError("window_idx must be None, a boolean mask, or integer indices")

    idx = idx.astype(np.int64, copy=False).ravel()
    if idx.size and (np.any(idx < 0) or np.any(idx >= n_windows)):
        raise IndexError(f"window_idx contains values outside [0, {n_windows})")
    return idx


def _as_int_radius(name: str, value: int) -> int:
    out = int(value)
    if out != value or out < 0:
        raise ValueError(f"{name} must be a non-negative integer")
    return out


def _as_nonnegative_float(name: str, value: float) -> float:
    out = float(value)
    if not np.isfinite(out) or out < 0:
        raise ValueError(f"{name} must be finite and >= 0")
    return out


def _as_positive_int(name: str, value: int) -> int:
    out = int(value)
    if out != value or out < 1:
        raise ValueError(f"{name} must be a positive integer")
    return out


def _nanmean(arr: np.ndarray) -> float:
    values = np.asarray(arr, dtype=np.float64)
    if values.size == 0 or not np.any(np.isfinite(values)):
        return float("nan")
    return float(np.nanmean(values))


def _nanmedian(arr: np.ndarray) -> float:
    values = np.asarray(arr, dtype=np.float64)
    if values.size == 0 or not np.any(np.isfinite(values)):
        return float("nan")
    return float(np.nanmedian(values))


def _posterior_rows(result: BayesianDecodingResult, idx: np.ndarray) -> np.ndarray:
    posterior = np.asarray(result.posterior_wx, dtype=np.float64)
    if posterior.ndim != 2:
        raise ValueError(f"posterior_wx must be 2D, got ndim={posterior.ndim}")

    rows = posterior[idx]
    normalized = np.full(rows.shape, np.nan, dtype=np.float64)
    if rows.size == 0:
        return normalized

    finite_rows = np.all(np.isfinite(rows), axis=1)
    row_sums = np.sum(rows, axis=1)
    valid = finite_rows & np.isfinite(row_sums) & (row_sums > 0)
    normalized[valid] = rows[valid] / row_sums[valid, None]
    return normalized


def _check_xbin_centers(result: BayesianDecodingResult) -> np.ndarray:
    x = np.asarray(result.xbin_centers, dtype=np.float64).ravel()
    n_bins = int(np.asarray(result.posterior_wx).shape[1])
    if x.size != n_bins:
        raise ValueError(f"xbin_centers length {x.size} must match posterior bins {n_bins}")
    return x


def _included_position_bins(
    result: BayesianDecodingResult,
    *,
    include_bins: np.ndarray | list[int] | list[bool] | None,
    exclude_edge_bins: int,
    x_min_cm: float | None,
    x_max_cm: float | None,
) -> np.ndarray:
    x = _check_xbin_centers(result)
    n_bins = int(x.size)
    include = np.ones(n_bins, dtype=bool)

    if include_bins is not None:
        raw = np.asarray(include_bins)
        if raw.dtype == bool:
            if raw.ndim != 1 or raw.size != n_bins:
                raise ValueError(
                    f"Boolean include_bins must have length {n_bins}, got shape {raw.shape}"
                )
            include &= raw
        elif np.issubdtype(raw.dtype, np.integer):
            idx = raw.astype(np.int64, copy=False).ravel()
            if idx.size and (np.any(idx < 0) or np.any(idx >= n_bins)):
                raise IndexError(f"include_bins contains values outside [0, {n_bins})")
            mask = np.zeros(n_bins, dtype=bool)
            mask[idx] = True
            include &= mask
        else:
            raise TypeError("include_bins must be None, a boolean mask, or integer bin indices")

    edge = _as_int_radius("exclude_edge_bins", exclude_edge_bins)
    if edge:
        include[:edge] = False
        include[max(0, n_bins - edge) :] = False

    if x_min_cm is not None:
        xmin = float(x_min_cm)
        if not np.isfinite(xmin):
            raise ValueError("x_min_cm must be finite when provided")
        include &= x >= xmin
    if x_max_cm is not None:
        xmax = float(x_max_cm)
        if not np.isfinite(xmax):
            raise ValueError("x_max_cm must be finite when provided")
        include &= x <= xmax
    if x_min_cm is not None and x_max_cm is not None and float(x_min_cm) > float(x_max_cm):
        raise ValueError("x_min_cm must be <= x_max_cm")

    return include


def decoding_error_cm(
    result: BayesianDecodingResult,
    window_idx: np.ndarray | list[int] | list[bool] | None = None,
) -> np.ndarray:
    """
    Return per-window argmax decoding error in centimeters.
    """

    idx = _window_indices(result, window_idx)
    return np.asarray(result.error_cm_w, dtype=np.float64)[idx]


def hard_decoding_correct(
    result: BayesianDecodingResult,
    *,
    tolerance_bins: int = 0,
    tolerance_cm: float | None = None,
    window_idx: np.ndarray | list[int] | list[bool] | None = None,
) -> np.ndarray:
    """
    Return whether each decoded argmax is correct within a bin or cm tolerance.
    """

    tol_bins = _as_int_radius("tolerance_bins", tolerance_bins)
    if tolerance_cm is not None and tol_bins != 0:
        raise ValueError("Pass either tolerance_bins or tolerance_cm, not both")

    idx = _window_indices(result, window_idx)
    if tolerance_cm is not None:
        tol_cm = _as_nonnegative_float("tolerance_cm", tolerance_cm)
        actual_x = np.asarray(result.actual_x_w, dtype=np.float64)[idx]
        decoded_x = np.asarray(result.decoded_x_w, dtype=np.float64)[idx]
        valid = np.isfinite(actual_x) & np.isfinite(decoded_x)
        return valid & (np.abs(decoded_x - actual_x) <= tol_cm)

    actual_bin = np.asarray(result.actual_bin_w, dtype=np.int64)[idx]
    decoded_bin = np.asarray(result.decoded_bin_w, dtype=np.int64)[idx]
    valid = (actual_bin >= 0) & (decoded_bin >= 0)
    return valid & (np.abs(decoded_bin - actual_bin) <= tol_bins)


def hard_decoding_accuracy(
    result: BayesianDecodingResult,
    *,
    tolerance_bins: int = 0,
    tolerance_cm: float | None = None,
    window_idx: np.ndarray | list[int] | list[bool] | None = None,
) -> float:
    """
    Return mean argmax decoding correctness.
    """

    correct = hard_decoding_correct(
        result,
        tolerance_bins=tolerance_bins,
        tolerance_cm=tolerance_cm,
        window_idx=window_idx,
    )
    if correct.size == 0:
        return float("nan")
    return float(np.mean(correct.astype(np.float64)))


def local_decoding_probability(
    result: BayesianDecodingResult,
    *,
    radius_bins: int = 0,
    radius_cm: float | None = None,
    window_idx: np.ndarray | list[int] | list[bool] | None = None,
) -> np.ndarray:
    """
    Return posterior mass assigned to the actual location or local neighborhood.
    """

    rad_bins = _as_int_radius("radius_bins", radius_bins)
    if radius_cm is not None and rad_bins != 0:
        raise ValueError("Pass either radius_bins or radius_cm, not both")

    idx = _window_indices(result, window_idx)
    posterior = _posterior_rows(result, idx)
    out = np.full(idx.size, np.nan, dtype=np.float64)
    if idx.size == 0:
        return out

    n_bins = int(posterior.shape[1])
    valid_rows = np.any(np.isfinite(posterior), axis=1)

    if radius_cm is not None:
        rad_cm = _as_nonnegative_float("radius_cm", radius_cm)
        x = _check_xbin_centers(result)
        actual_x = np.asarray(result.actual_x_w, dtype=np.float64)[idx]
        for i, actual in enumerate(actual_x):
            if not valid_rows[i] or not np.isfinite(actual):
                continue
            out[i] = float(np.nansum(posterior[i, np.abs(x - actual) <= rad_cm]))
        return out

    actual_bin = np.asarray(result.actual_bin_w, dtype=np.int64)[idx]
    for i, actual in enumerate(actual_bin):
        if not valid_rows[i] or actual < 0 or actual >= n_bins:
            continue
        start = max(0, int(actual) - rad_bins)
        stop = min(n_bins, int(actual) + rad_bins + 1)
        out[i] = float(np.nansum(posterior[i, start:stop]))
    return out


def bayesian_decoding_accuracy(
    result: BayesianDecodingResult,
    *,
    radius_bins: int = 0,
    radius_cm: float | None = None,
    window_idx: np.ndarray | list[int] | list[bool] | None = None,
) -> float:
    """
    Return mean posterior mass assigned to the actual location.
    """

    return _nanmean(
        local_decoding_probability(
            result,
            radius_bins=radius_bins,
            radius_cm=radius_cm,
            window_idx=window_idx,
        )
    )


def local_decoding_probability_by_actual_bin(
    result: BayesianDecodingResult,
    *,
    radius_bins: int = 0,
    radius_cm: float | None = None,
    window_idx: np.ndarray | list[int] | list[bool] | None = None,
    include_bins: np.ndarray | list[int] | list[bool] | None = None,
    exclude_edge_bins: int = 0,
    x_min_cm: float | None = None,
    x_max_cm: float | None = None,
    min_windows_per_bin: int = 1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Return mean local decoding probability for each actual position bin.

    Excluded or insufficiently sampled bins have NaN mean probability.
    """

    min_windows = _as_positive_int("min_windows_per_bin", min_windows_per_bin)
    idx = _window_indices(result, window_idx)
    x = _check_xbin_centers(result)
    include = _included_position_bins(
        result,
        include_bins=include_bins,
        exclude_edge_bins=exclude_edge_bins,
        x_min_cm=x_min_cm,
        x_max_cm=x_max_cm,
    )

    mean_probability_x = np.full(x.size, np.nan, dtype=np.float64)
    n_windows_x = np.zeros(x.size, dtype=np.int64)
    if idx.size == 0:
        return x.copy(), mean_probability_x, n_windows_x

    probability_w = local_decoding_probability(
        result,
        radius_bins=radius_bins,
        radius_cm=radius_cm,
        window_idx=idx,
    )
    actual_bin = np.asarray(result.actual_bin_w, dtype=np.int64)[idx]
    valid = np.isfinite(probability_w) & (actual_bin >= 0) & (actual_bin < x.size)

    for actual in np.flatnonzero(include):
        in_bin = valid & (actual_bin == int(actual))
        n_windows = int(np.sum(in_bin))
        n_windows_x[actual] = n_windows
        if n_windows >= min_windows:
            mean_probability_x[actual] = float(np.nanmean(probability_w[in_bin]))

    return x.copy(), mean_probability_x, n_windows_x


def position_balanced_bayesian_decoding_accuracy(
    result: BayesianDecodingResult,
    *,
    radius_bins: int = 0,
    radius_cm: float | None = None,
    window_idx: np.ndarray | list[int] | list[bool] | None = None,
    include_bins: np.ndarray | list[int] | list[bool] | None = None,
    exclude_edge_bins: int = 0,
    x_min_cm: float | None = None,
    x_max_cm: float | None = None,
    min_windows_per_bin: int = 1,
) -> float:
    """
    Return local posterior accuracy averaged equally over actual position bins.
    """

    _, probability_x, _ = local_decoding_probability_by_actual_bin(
        result,
        radius_bins=radius_bins,
        radius_cm=radius_cm,
        window_idx=window_idx,
        include_bins=include_bins,
        exclude_edge_bins=exclude_edge_bins,
        x_min_cm=x_min_cm,
        x_max_cm=x_max_cm,
        min_windows_per_bin=min_windows_per_bin,
    )
    return _nanmean(probability_x)


def posterior_peak_probability(
    result: BayesianDecodingResult,
    window_idx: np.ndarray | list[int] | list[bool] | None = None,
) -> np.ndarray:
    """
    Return the maximum posterior probability per decoded window.
    """

    idx = _window_indices(result, window_idx)
    posterior = _posterior_rows(result, idx)
    out = np.full(idx.size, np.nan, dtype=np.float64)
    valid = np.any(np.isfinite(posterior), axis=1)
    if np.any(valid):
        out[valid] = np.nanmax(posterior[valid], axis=1)
    return out


def posterior_entropy_bits(
    result: BayesianDecodingResult,
    *,
    normalized: bool = False,
    window_idx: np.ndarray | list[int] | list[bool] | None = None,
) -> np.ndarray:
    """
    Return Shannon entropy of the posterior per decoded window.
    """

    idx = _window_indices(result, window_idx)
    posterior = _posterior_rows(result, idx)
    out = np.full(idx.size, np.nan, dtype=np.float64)
    valid = np.any(np.isfinite(posterior), axis=1)
    if not np.any(valid):
        return out

    p = np.where(np.isfinite(posterior[valid]) & (posterior[valid] > 0), posterior[valid], 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        entropy = -np.sum(np.where(p > 0, p * np.log2(p), 0.0), axis=1)

    if normalized:
        n_bins = int(posterior.shape[1])
        entropy = np.zeros_like(entropy) if n_bins <= 1 else entropy / np.log2(float(n_bins))

    out[valid] = entropy
    return out


def posterior_mean_x(
    result: BayesianDecodingResult,
    window_idx: np.ndarray | list[int] | list[bool] | None = None,
) -> np.ndarray:
    """
    Return posterior expected position per decoded window.
    """

    idx = _window_indices(result, window_idx)
    posterior = _posterior_rows(result, idx)
    x = _check_xbin_centers(result)
    out = np.full(idx.size, np.nan, dtype=np.float64)
    valid = np.any(np.isfinite(posterior), axis=1)
    if np.any(valid):
        out[valid] = np.sum(posterior[valid] * x[None, :], axis=1)
    return out


def posterior_std_cm(
    result: BayesianDecodingResult,
    window_idx: np.ndarray | list[int] | list[bool] | None = None,
) -> np.ndarray:
    """
    Return posterior spatial standard deviation per decoded window.
    """

    idx = _window_indices(result, window_idx)
    posterior = _posterior_rows(result, idx)
    x = _check_xbin_centers(result)
    mean_x = posterior_mean_x(result, window_idx=idx)
    out = np.full(idx.size, np.nan, dtype=np.float64)
    valid = np.any(np.isfinite(posterior), axis=1) & np.isfinite(mean_x)
    if np.any(valid):
        out[valid] = np.sqrt(
            np.sum(posterior[valid] * (x[None, :] - mean_x[valid, None]) ** 2, axis=1)
        )
    return out


def posterior_expected_abs_error_cm(
    result: BayesianDecodingResult,
    window_idx: np.ndarray | list[int] | list[bool] | None = None,
) -> np.ndarray:
    """
    Return posterior-expected absolute distance from the actual position.
    """

    idx = _window_indices(result, window_idx)
    posterior = _posterior_rows(result, idx)
    x = _check_xbin_centers(result)
    actual_x = np.asarray(result.actual_x_w, dtype=np.float64)[idx]
    out = np.full(idx.size, np.nan, dtype=np.float64)
    valid = np.any(np.isfinite(posterior), axis=1) & np.isfinite(actual_x)
    if np.any(valid):
        out[valid] = np.sum(posterior[valid] * np.abs(x[None, :] - actual_x[valid, None]), axis=1)
    return out


def summarize_decoding_measures(
    result: BayesianDecodingResult,
    *,
    window_idx: np.ndarray | list[int] | list[bool] | None = None,
    hard_tolerance_bins: int = 0,
    hard_tolerance_cm: float | None = None,
    local_radius_bins: int = 0,
    local_radius_cm: float | None = None,
) -> dict[str, Any]:
    """
    Return compact scalar measures for a Bayesian decoding result.
    """

    idx = _window_indices(result, window_idx)
    errors = decoding_error_cm(result, idx)
    return {
        "n_windows": int(idx.size),
        "median_error_cm": _nanmedian(errors),
        "mean_error_cm": _nanmean(errors),
        "hard_accuracy": hard_decoding_accuracy(
            result,
            tolerance_bins=hard_tolerance_bins,
            tolerance_cm=hard_tolerance_cm,
            window_idx=idx,
        ),
        "bayesian_accuracy": bayesian_decoding_accuracy(
            result,
            radius_bins=local_radius_bins,
            radius_cm=local_radius_cm,
            window_idx=idx,
        ),
        "position_balanced_bayesian_accuracy": position_balanced_bayesian_decoding_accuracy(
            result,
            radius_bins=local_radius_bins,
            radius_cm=local_radius_cm,
            window_idx=idx,
        ),
        "mean_peak_probability": _nanmean(posterior_peak_probability(result, idx)),
        "mean_posterior_entropy_bits": _nanmean(posterior_entropy_bits(result, window_idx=idx)),
        "mean_posterior_std_cm": _nanmean(posterior_std_cm(result, idx)),
        "mean_posterior_expected_abs_error_cm": _nanmean(
            posterior_expected_abs_error_cm(result, idx)
        ),
    }
