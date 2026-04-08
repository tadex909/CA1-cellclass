from __future__ import annotations

import numpy as np


def gaussian_kernel_1d(sigma_bins: float) -> np.ndarray:
    if sigma_bins <= 0:
        return np.array([1.0], dtype=np.float64)
    radius = max(1, int(np.ceil(4.0 * sigma_bins)))
    x = np.arange(-radius, radius + 1, dtype=np.float64)
    k = np.exp(-(x**2) / (2.0 * sigma_bins**2))
    k /= np.sum(k)
    return k


def _convolve_same_as_input(x: np.ndarray, k: np.ndarray) -> np.ndarray:
    """
    MATLAB-like conv(...,'same') size behavior:
    output length equals input signal length, even when kernel is longer.
    """
    v = np.asarray(x, dtype=np.float64).ravel()
    kk = np.asarray(k, dtype=np.float64).ravel()
    full = np.convolve(v, kk, mode="full")
    start = (kk.size - 1) // 2
    end = start + v.size
    return full[start:end]


def smooth_last_axis(arr: np.ndarray, sigma_bins: float) -> np.ndarray:
    if sigma_bins <= 0:
        return np.asarray(arr, dtype=np.float64).copy()
    k = gaussian_kernel_1d(float(sigma_bins))
    x = np.asarray(arr, dtype=np.float64)
    if x.ndim == 1:
        return _convolve_same_as_input(x, k)
    if x.ndim == 2:
        out = np.empty_like(x, dtype=np.float64)
        for i in range(x.shape[0]):
            out[i] = _convolve_same_as_input(x[i], k)
        return out
    if x.ndim == 3:
        out = np.empty_like(x, dtype=np.float64)
        for i in range(x.shape[0]):
            for j in range(x.shape[1]):
                out[i, j] = _convolve_same_as_input(x[i, j], k)
        return out
    raise ValueError(f"Unsupported ndim for smoothing: {x.ndim}")


def divide_with_nan(num: np.ndarray, den: np.ndarray) -> np.ndarray:
    out = np.full(np.broadcast_shapes(num.shape, den.shape), np.nan, dtype=np.float64)
    np.divide(num, den, out=out, where=(den > 0))
    return out
