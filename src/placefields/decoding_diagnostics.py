from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np
import pandas as pd

from .cue_zones import (
    cue_zone_layout_for_condition_on_track,
    label_positions_by_zone,
    label_positions_by_zone_component,
)
from .decoding import BayesianDecodingResult
from .population_geometry import decode_condway
from .trials import canonical_condition_name


def condition_label_for_condway(
    condway: int,
    condition_names_by_base: Mapping[int, str] | None = None,
) -> str:
    base_condition, direction = decode_condway(int(condway))
    condition_names_by_base = condition_names_by_base or {}
    raw_name = condition_names_by_base.get(int(base_condition), f"cond{int(base_condition)}")
    condition_name = canonical_condition_name(raw_name) or str(raw_name)
    return f"{condition_name} {direction}"


def bayesian_decoding_result_to_frame(
    result: BayesianDecodingResult,
    *,
    freq_hz: float,
    condition_names_by_base: Mapping[int, str] | None = None,
    include_cue_zones: bool = True,
    cue_xbin_edges: np.ndarray | None = None,
) -> pd.DataFrame:
    """
    Convert a BayesianDecodingResult into one tidy row per decoded time window.
    """

    if not np.isfinite(freq_hz) or float(freq_hz) <= 0:
        raise ValueError("freq_hz must be finite and > 0")

    condways = np.asarray(result.condway_w, dtype=np.int64)
    decoded = [decode_condway(int(condway)) for condway in condways]
    base_conditions = [int(base) for base, _ in decoded]
    directions = [str(direction) for _, direction in decoded]
    condition_names_by_base = condition_names_by_base or {}

    condition_names: list[str] = []
    condition_labels: list[str] = []
    for base_condition, direction in zip(base_conditions, directions):
        raw_name = condition_names_by_base.get(int(base_condition), f"cond{int(base_condition)}")
        condition_name = canonical_condition_name(raw_name) or str(raw_name)
        condition_names.append(condition_name)
        condition_labels.append(f"{condition_name} {direction}")

    frame = pd.DataFrame(
        {
            "window_index": np.arange(condways.size, dtype=np.int64),
            "trial_index": result.trial_index_w,
            "condway": result.condway_w,
            "base_condition": base_conditions,
            "condition_name": condition_names,
            "condition_label": condition_labels,
            "direction": directions,
            "window_start": result.window_start_w,
            "window_stop": result.window_stop_w,
            "window_mid_s": (
                0.5
                * (np.asarray(result.window_start_w, dtype=np.float64) + np.asarray(result.window_stop_w, dtype=np.float64))
                / float(freq_hz)
            ),
            "actual_x": result.actual_x_w,
            "decoded_x": result.decoded_x_w,
            "actual_bin": result.actual_bin_w,
            "decoded_bin": result.decoded_bin_w,
            "error_cm": result.error_cm_w,
            "prob_actual": result.prob_actual_w,
            "n_spikes": result.n_spikes_w,
            "n_train_laps": result.n_train_laps_w,
            "train_group": result.train_group_w,
            "train_group_label": result.train_group_label_w.astype(str),
        }
    )
    if include_cue_zones:
        edges = result.xbin_edges if cue_xbin_edges is None else cue_xbin_edges
        actual_x = np.asarray(result.actual_x_w, dtype=np.float64)
        base_arr = np.asarray(base_conditions, dtype=np.int64)
        cue_zone = np.full(condways.size, "", dtype="<U5")
        cue_zone_component = np.full(condways.size, "", dtype="<U16")

        for base_condition in sorted(set(base_conditions)):
            row_mask = base_arr == int(base_condition)
            if not np.any(row_mask):
                continue
            raw_name = condition_names_by_base.get(
                int(base_condition),
                f"cond{int(base_condition)}",
            )
            condition_name = canonical_condition_name(raw_name) or str(raw_name)
            layout = cue_zone_layout_for_condition_on_track(
                condition_name,
                xbin_edges=edges,
            )
            cue_zone[row_mask] = label_positions_by_zone(
                actual_x[row_mask],
                layout=layout,
            )
            cue_zone_component[row_mask] = label_positions_by_zone_component(
                actual_x[row_mask],
                layout=layout,
            )

        frame["cue_zone"] = cue_zone
        frame["cue_zone_component"] = cue_zone_component
    return frame


def _unique_join(values: pd.Series) -> str:
    unique = sorted({str(value) for value in values.dropna().astype(str) if str(value)})
    return ",".join(unique)


def _summarize_group(group: pd.DataFrame, group_keys: dict[str, Any]) -> dict[str, Any]:
    row: dict[str, Any] = dict(group_keys)
    row.update(
        {
            "n_windows": int(len(group)),
            "n_laps_decoded": int(group["trial_index"].nunique()) if "trial_index" in group else 0,
            "median_error_cm": float(group["error_cm"].median()) if "error_cm" in group else np.nan,
            "mean_error_cm": float(group["error_cm"].mean()) if "error_cm" in group else np.nan,
            "mean_prob_actual": float(group["prob_actual"].mean()) if "prob_actual" in group else np.nan,
            "mean_spikes_per_window": float(group["n_spikes"].mean()) if "n_spikes" in group else np.nan,
            "mean_train_laps": float(group["n_train_laps"].mean()) if "n_train_laps" in group else np.nan,
        }
    )
    if "train_group" in group:
        row["train_group"] = _unique_join(group["train_group"])
    if "train_group_label" in group:
        row["train_group_label"] = _unique_join(group["train_group_label"])
    return row


def summarize_decoding_by_condition(
    decode_df: pd.DataFrame,
    *,
    group_cols: list[str] | None = None,
) -> pd.DataFrame:
    """
    Summarize decoded windows by condition/condition-direction.
    """

    if group_cols is None:
        group_cols = [col for col in ["condway", "condition_label", "direction"] if col in decode_df.columns]
    if not group_cols:
        raise ValueError("No grouping columns available. Pass group_cols explicitly.")

    rows: list[dict[str, Any]] = []
    for keys, group in decode_df.groupby(group_cols, sort=True, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        rows.append(_summarize_group(group, dict(zip(group_cols, keys))))
    return pd.DataFrame(rows).reset_index(drop=True)


def summarize_decoding_by_trial(
    decode_df: pd.DataFrame,
    *,
    group_cols: list[str] | None = None,
) -> pd.DataFrame:
    """
    Summarize decoded windows by trial/lap.
    """

    if group_cols is None:
        group_cols = [col for col in ["trial_index", "condway", "condition_label"] if col in decode_df.columns]
    if not group_cols:
        raise ValueError("No grouping columns available. Pass group_cols explicitly.")

    rows: list[dict[str, Any]] = []
    for keys, group in decode_df.groupby(group_cols, sort=True, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = _summarize_group(group, dict(zip(group_cols, keys)))
        if "error_cm" in group:
            row["total_abs_error_cm"] = float(group["error_cm"].sum())
            row["mean_abs_error_cm"] = float(group["error_cm"].mean())
            row["median_abs_error_cm"] = float(group["error_cm"].median())
        rows.append(row)
    return pd.DataFrame(rows).reset_index(drop=True)


def soft_confusion_matrix_for_windows(
    result: BayesianDecodingResult,
    window_idx: np.ndarray | list[int],
    *,
    direction: str | None = None,
    flip_backward: bool = False,
) -> tuple[np.ndarray, int]:
    """
    Build a row-normalized soft confusion matrix from posterior probabilities.

    Rows are actual position bins. Columns are decoded position bins. Each
    decoded window contributes its full posterior vector to the row
    corresponding to the actual bin.
    """

    n_bins = int(result.xbin_centers.size)
    idx = np.asarray(window_idx)
    if idx.dtype == bool:
        idx = np.flatnonzero(idx)
    idx = idx.astype(np.int64, copy=False).ravel()

    posterior = np.asarray(result.posterior_wx[idx], dtype=np.float64)
    actual_bin = np.asarray(result.actual_bin_w[idx], dtype=np.int64)

    valid = (actual_bin >= 0) & (actual_bin < n_bins)
    valid &= np.all(np.isfinite(posterior), axis=1)
    row_sums = np.sum(posterior, axis=1)
    valid &= np.isfinite(row_sums) & (row_sums > 0)

    posterior = posterior[valid]
    actual_bin = actual_bin[valid]
    row_sums = row_sums[valid]
    if posterior.shape[0] == 0:
        return np.full((n_bins, n_bins), np.nan, dtype=np.float64), 0

    posterior = posterior / row_sums[:, None]

    if flip_backward and direction is not None and str(direction).upper() == "B":
        posterior = posterior[:, ::-1]
        actual_bin = (n_bins - 1) - actual_bin

    accum = np.zeros((n_bins, n_bins), dtype=np.float64)
    np.add.at(accum, actual_bin, posterior)

    actual_row_sums = accum.sum(axis=1, keepdims=True)
    matrix = np.full_like(accum, np.nan, dtype=np.float64)
    valid_rows = actual_row_sums[:, 0] > 0
    matrix[valid_rows] = accum[valid_rows] / actual_row_sums[valid_rows]
    return matrix, int(posterior.shape[0])


def average_soft_confusions_by_condition(
    result: BayesianDecodingResult,
    decode_df: pd.DataFrame,
    *,
    label_col: str = "condition_label",
    trial_col: str = "trial_index",
    direction_col: str = "direction",
    flip_backward: bool = False,
) -> tuple[dict[str, np.ndarray], pd.DataFrame]:
    """
    Average lap-level soft confusion matrices within each condition label.
    """

    if trial_col not in decode_df.columns:
        raise KeyError(f"decode_df is missing required column: {trial_col}")
    if label_col not in decode_df.columns:
        label_col = "condway"
    if label_col not in decode_df.columns:
        raise KeyError("decode_df must contain condition_label or condway")

    work = decode_df.copy()
    if "window_index" not in work.columns:
        work = work.reset_index().rename(columns={"index": "window_index"})

    group_cols = [trial_col, label_col]
    if direction_col in work.columns:
        group_cols.append(direction_col)

    matrices_by_label: dict[str, list[np.ndarray]] = {}
    trial_rows: list[dict[str, Any]] = []
    for keys, group in work.groupby(group_cols, sort=False, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        key_map = dict(zip(group_cols, keys))
        label = str(key_map[label_col])
        direction = str(key_map[direction_col]) if direction_col in key_map else None
        matrix, n_windows = soft_confusion_matrix_for_windows(
            result,
            group["window_index"].to_numpy(dtype=np.int64),
            direction=direction,
            flip_backward=flip_backward,
        )
        if n_windows == 0:
            continue
        matrices_by_label.setdefault(label, []).append(matrix)
        trial_rows.append(
            {
                "trial_index": int(key_map[trial_col]),
                "condition_label": label,
                "direction": "" if direction is None else direction,
                "n_windows": int(n_windows),
                "n_actual_bins": int(np.sum(np.any(np.isfinite(matrix), axis=1))),
            }
        )

    averages: dict[str, np.ndarray] = {}
    for label, matrices in matrices_by_label.items():
        stack = np.stack(matrices, axis=0)
        counts = np.sum(np.isfinite(stack), axis=0)
        summed = np.nansum(stack, axis=0)
        avg = np.full_like(summed, np.nan, dtype=np.float64)
        avg[counts > 0] = summed[counts > 0] / counts[counts > 0]
        averages[label] = avg

    return averages, pd.DataFrame(trial_rows)


def summarize_soft_confusions(
    average_confusion_by_condition: Mapping[str, np.ndarray],
    trial_confusion_df: pd.DataFrame,
    xbin_centers: np.ndarray,
    *,
    near_diagonal_cm: float = 10.0,
) -> pd.DataFrame:
    """
    Summarize average soft confusion matrices per condition.
    """

    x = np.asarray(xbin_centers, dtype=np.float64).ravel()
    distance_cm = np.abs(x[:, None] - x[None, :])
    near_mask = distance_cm <= float(near_diagonal_cm)
    rows: list[dict[str, Any]] = []

    for label, matrix_raw in average_confusion_by_condition.items():
        matrix = np.asarray(matrix_raw, dtype=np.float64)
        if matrix.shape != distance_cm.shape:
            raise ValueError(
                f"Confusion matrix for {label!r} has shape {matrix.shape}; "
                f"expected {distance_cm.shape}"
            )
        valid_rows = np.any(np.isfinite(matrix), axis=1)
        matrix_filled = np.nan_to_num(matrix, nan=0.0)
        diag = np.diag(matrix)
        near_mass_by_row = np.sum(matrix_filled * near_mask, axis=1)
        expected_abs_error_by_row = np.sum(matrix_filled * distance_cm, axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            entropy_by_row = -np.sum(
                np.where(matrix_filled > 0, matrix_filled * np.log2(matrix_filled), 0.0),
                axis=1,
            )

        label_trials = (
            trial_confusion_df[trial_confusion_df["condition_label"].astype(str) == str(label)]
            if "condition_label" in trial_confusion_df.columns
            else pd.DataFrame()
        )
        rows.append(
            {
                "condition_label": str(label),
                "n_trials": int(label_trials["trial_index"].nunique()) if "trial_index" in label_trials else 0,
                "n_windows": int(label_trials["n_windows"].sum()) if "n_windows" in label_trials else 0,
                "mean_diagonal_mass": (
                    float(np.nanmean(diag)) if np.any(np.isfinite(diag)) else np.nan
                ),
                f"mean_near_diagonal_mass_{near_diagonal_cm:g}cm": (
                    float(np.mean(near_mass_by_row[valid_rows])) if np.any(valid_rows) else np.nan
                ),
                "expected_abs_error_cm_from_matrix": (
                    float(np.mean(expected_abs_error_by_row[valid_rows])) if np.any(valid_rows) else np.nan
                ),
                "mean_row_entropy_bits": (
                    float(np.mean(entropy_by_row[valid_rows])) if np.any(valid_rows) else np.nan
                ),
            }
        )

    return pd.DataFrame(rows).sort_values("condition_label").reset_index(drop=True)
