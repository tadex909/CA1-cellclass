from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .matlab_compat import matlab_1b_to_python_0b


@dataclass(frozen=True)
class TrialInfo:
    """
    One trial descriptor, using Python indexing conventions.

    Notes
    -----
    - `start_idx_0b` is inclusive.
    - `stop_idx_0b_exclusive` is exclusive (ready for Python slicing).
    """
    trial_index: int
    cond: int
    wb: str
    condway: int
    start_idx_0b: int
    stop_idx_0b_exclusive: int


def build_condway(cond: np.ndarray, wb: np.ndarray) -> np.ndarray:
    """
    Build condition-direction labels equivalent to the MATLAB logic in Rmap_G.m.

    For standard W/B coding:
    - condition i + W -> 2*i - 1
    - condition i + B -> 2*i
    """
    c = np.asarray(cond, dtype=np.int64).ravel()
    wbs = np.asarray(wb).astype(str).ravel()
    if c.size != wbs.size:
        raise ValueError(f"cond and wb size mismatch: {c.size} vs {wbs.size}")

    out = np.full(c.size, -1, dtype=np.int64)
    is_w = np.char.upper(wbs) == "W"
    is_b = np.char.upper(wbs) == "B"
    bad = ~(is_w | is_b)
    if np.any(bad):
        bad_vals = sorted(set(wbs[bad].tolist()))
        raise ValueError(f"Unsupported WB values: {bad_vals}. Expected only 'W'/'B'.")

    out[is_w] = 2 * c[is_w] - 1
    out[is_b] = 2 * c[is_b]
    return out


def build_trial_info_from_traj(
    cond: np.ndarray,
    wb: np.ndarray,
    start_1b: np.ndarray,
    stop_1b: np.ndarray,
    *,
    n_samples: int,
) -> list[TrialInfo]:
    """
    Construct per-trial metadata with Python-native indexing.

    `start_1b`/`stop_1b` are expected to come from MATLAB `Traj.start/Traj.stop`,
    where `stop` is inclusive in 1-based indexing.
    """
    c = np.asarray(cond, dtype=np.int64).ravel()
    w = np.asarray(wb).astype(str).ravel()
    s1 = np.asarray(start_1b, dtype=np.int64).ravel()
    e1 = np.asarray(stop_1b, dtype=np.int64).ravel()

    n = c.size
    if not (w.size == n and s1.size == n and e1.size == n):
        raise ValueError(
            f"size mismatch: cond={c.size}, wb={w.size}, start={s1.size}, stop={e1.size}"
        )
    if n_samples <= 0:
        raise ValueError("n_samples must be > 0")

    condway = build_condway(c, w)
    s0 = matlab_1b_to_python_0b(s1)
    # MATLAB inclusive 1b stop -> Python exclusive 0b stop keeps same integer value.
    e0_excl = np.asarray(e1, dtype=np.int64).copy()
    e0_excl[e0_excl < 0] = 0
    e0_excl[e0_excl > n_samples] = n_samples
    s0[s0 > n_samples] = n_samples

    trials: list[TrialInfo] = []
    for i in range(n):
        trials.append(
            TrialInfo(
                trial_index=int(i),
                cond=int(c[i]),
                wb=str(w[i]),
                condway=int(condway[i]),
                start_idx_0b=int(s0[i]),
                stop_idx_0b_exclusive=int(e0_excl[i]),
            )
        )
    return trials
