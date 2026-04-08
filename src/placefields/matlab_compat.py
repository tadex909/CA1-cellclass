from __future__ import annotations

import numpy as np


def ifreq_swap(indices: np.ndarray, f1_hz: float, f2_hz: float) -> np.ndarray:
    """
    MATLAB-compatible frequency index conversion.

    Equivalent to `fct_ifreq_swap.m`:
        if2 = floor((f2 / f1) * (if1 - 1)) + 1

    Parameters
    ----------
    indices
        1-based sample indices at frequency `f1_hz`.
    f1_hz
        Source frequency in Hz.
    f2_hz
        Target frequency in Hz.

    Returns
    -------
    np.ndarray
        Converted 1-based indices at frequency `f2_hz` (int64), same shape as input.
    """
    if f1_hz <= 0 or f2_hz <= 0:
        raise ValueError("f1_hz and f2_hz must be > 0")

    idx = np.asarray(indices, dtype=np.float64)
    out = np.floor((float(f2_hz) / float(f1_hz)) * (idx - 1.0)) + 1.0
    out = out.astype(np.int64)
    out[out < 1] = 1
    return out.reshape(np.asarray(indices).shape)


def matlab_1b_to_python_0b(indices_1b: np.ndarray) -> np.ndarray:
    """
    Convert MATLAB 1-based indices to Python 0-based indices.
    """
    idx = np.asarray(indices_1b, dtype=np.int64)
    out = idx - 1
    out[out < 0] = 0
    return out
