from __future__ import annotations

"""Boundary validators for cellclass data products."""

from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from cellclass.config import DEFAULT_AGE_GROUPS


class ValidationError(ValueError):
    """Raised when a pipeline data object violates its expected schema."""


ALLCEL_REQUIRED_KEYS: tuple[str, ...] = (
    "allcel__time_spk",
    "allcel__id_spk",
    "allcel__id_cel",
    "allcel__type_u",
    "allcel__bestswaveforms",
    "allpf__ispf_cxu",
)

ALLCEL_PER_CELL_ALIASES: dict[str, tuple[str, ...]] = {
    "burst": ("allcel__burst_u", "burst_u"),
    "firing_rate": ("allcel__fr_u", "fr_u"),
    "duration": ("allcel__duration_u", "duration_u"),
    "asymmetry": ("allcel__asymmetry_u", "asymmetry_u"),
}

FEATURE_TABLE_REQUIRED_COLUMNS: tuple[str, ...] = (
    "session_id",
    "mouse",
    "date",
    "time",
    "cell_id",
    "allcel__type_u",
    "n_spikes",
    "fr_hz",
    "fr_hz_session",
    "cv2",
    "refractory_ms_center",
    "refractory_ms_edge",
    "burst_index",
    "acg_peak_latency_ms",
    "spk_duration_ms",
    "spk_peaktrough_ms",
    "spk_asymmetry",
    "qc_min_spikes",
    "qc_refractory",
    "qc_waveform",
)

AGE_GROUP_TABLE_REQUIRED_COLUMNS: tuple[str, ...] = (
    "session_id",
    "cell_id",
    "unit_uid",
    "Age",
    "age_group",
)

IDENTIFIER_COLUMNS: tuple[str, ...] = ("session_id", "cell_id")

NUMERIC_FEATURE_COLUMNS: tuple[str, ...] = (
    "n_spikes",
    "fr_hz",
    "fr_hz_session",
    "cv2",
    "refractory_ms_center",
    "refractory_ms_edge",
    "burst_index",
    "acg_peak_latency_ms",
    "spk_duration_ms",
    "spk_peaktrough_ms",
    "spk_asymmetry",
)


def _source_name(source: str | Path | None, default: str) -> str:
    return str(source) if source is not None else default


def _raise_if_errors(errors: list[str], *, source: str | Path | None, default: str) -> None:
    if not errors:
        return
    label = _source_name(source, default)
    detail = "\n  - ".join(errors)
    raise ValidationError(f"{label} failed validation:\n  - {detail}")


def _keys(data: Mapping[str, Any] | np.lib.npyio.NpzFile) -> set[str]:
    return set(data.keys())


def _as_squeezed_array(data: Mapping[str, Any] | np.lib.npyio.NpzFile, key: str) -> np.ndarray:
    return np.squeeze(np.asarray(data[key]))


def _require_keys(
    keys: set[str],
    required: Iterable[str],
    errors: list[str],
) -> None:
    missing = [key for key in required if key not in keys]
    if missing:
        errors.append(f"missing required key(s): {', '.join(missing)}")


def _resolve_alias(
    keys: set[str],
    aliases: tuple[str, ...],
    *,
    label: str,
    errors: list[str],
) -> str | None:
    for key in aliases:
        if key in keys:
            return key
    errors.append(f"missing required {label} key; expected one of {aliases}")
    return None


def _check_1d_length(
    arr: np.ndarray,
    *,
    name: str,
    length: int | None,
    errors: list[str],
) -> None:
    if arr.ndim != 1:
        errors.append(f"{name} must be 1D after squeeze; got shape {arr.shape}")
        return
    if length is not None and arr.size != length:
        errors.append(f"{name} length must be {length}; got {arr.size}")


def _check_numeric_array(
    arr: np.ndarray,
    *,
    name: str,
    errors: list[str],
    integer: bool = False,
    finite: bool = True,
) -> None:
    if integer and arr.dtype.kind not in "iu":
        errors.append(f"{name} must have integer dtype; got {arr.dtype}")
        return
    if not integer and arr.dtype.kind not in "biuf":
        errors.append(f"{name} must have numeric dtype; got {arr.dtype}")
        return
    if finite and not np.all(np.isfinite(arr.astype(np.float64, copy=False))):
        errors.append(f"{name} contains non-finite values")


def validate_allcel_npz(
    data: Mapping[str, Any] | np.lib.npyio.NpzFile,
    *,
    source: str | Path | None = None,
    n_condition_pairs: int = 5,
) -> None:
    """
    Validate the interim ratemap/allcel NPZ consumed by cellclass.pipeline.

    The contract is intentionally about the fields used by feature extraction:
    spike times and IDs align, per-cell arrays align with allcel__id_cel,
    waveforms have cell axis last, and allpf__ispf_cxu has condition/cell axes.
    """
    errors: list[str] = []
    keys = _keys(data)
    _require_keys(keys, ALLCEL_REQUIRED_KEYS, errors)

    alias_keys = {
        label: _resolve_alias(keys, aliases, label=label, errors=errors)
        for label, aliases in ALLCEL_PER_CELL_ALIASES.items()
    }

    if errors:
        _raise_if_errors(errors, source=source, default="allcel NPZ")

    spike_times = _as_squeezed_array(data, "allcel__time_spk")
    spike_ids = _as_squeezed_array(data, "allcel__id_spk")
    cell_ids = _as_squeezed_array(data, "allcel__id_cel")
    type_u = _as_squeezed_array(data, "allcel__type_u")

    _check_1d_length(spike_times, name="allcel__time_spk", length=None, errors=errors)
    _check_1d_length(spike_ids, name="allcel__id_spk", length=spike_times.size, errors=errors)
    _check_1d_length(cell_ids, name="allcel__id_cel", length=None, errors=errors)
    _check_1d_length(type_u, name="allcel__type_u", length=cell_ids.size, errors=errors)

    if cell_ids.ndim == 1 and cell_ids.size == 0:
        errors.append("allcel__id_cel must contain at least one cell")
    if cell_ids.ndim == 1 and np.unique(cell_ids).size != cell_ids.size:
        errors.append("allcel__id_cel contains duplicate cell IDs")

    if spike_times.ndim == 1:
        _check_numeric_array(spike_times, name="allcel__time_spk", errors=errors)
    if spike_ids.ndim == 1:
        _check_numeric_array(spike_ids, name="allcel__id_spk", errors=errors, integer=True)
    if cell_ids.ndim == 1:
        _check_numeric_array(cell_ids, name="allcel__id_cel", errors=errors, integer=True)

    for label, key in alias_keys.items():
        if key is None:
            continue
        arr = _as_squeezed_array(data, key)
        _check_1d_length(arr, name=key, length=cell_ids.size, errors=errors)
        if arr.ndim == 1:
            _check_numeric_array(arr, name=key, errors=errors)

    waveforms = np.asarray(data["allcel__bestswaveforms"])
    if waveforms.ndim != 3:
        errors.append(f"allcel__bestswaveforms must have shape (time, waveform, cell); got {waveforms.shape}")
    elif waveforms.shape[2] != cell_ids.size:
        errors.append(
            "allcel__bestswaveforms cell axis must match allcel__id_cel length "
            f"({cell_ids.size}); got shape {waveforms.shape}"
        )
    elif waveforms.dtype.kind not in "biuf":
        errors.append(f"allcel__bestswaveforms must be numeric; got {waveforms.dtype}")

    ispf = _as_squeezed_array(data, "allpf__ispf_cxu")
    expected_conditions = 2 * n_condition_pairs
    if ispf.ndim == 2:
        if cell_ids.size != 1:
            errors.append(
                f"allpf__ispf_cxu is 2D but allcel__id_cel has {cell_ids.size} cells; expected 3D"
            )
        if expected_conditions not in ispf.shape:
            errors.append(
                f"allpf__ispf_cxu must include a condition axis of length {expected_conditions}; got {ispf.shape}"
            )
    elif ispf.ndim == 3:
        if cell_ids.size not in ispf.shape:
            errors.append(
                f"allpf__ispf_cxu must include a cell axis of length {cell_ids.size}; got {ispf.shape}"
            )
        if expected_conditions not in ispf.shape:
            errors.append(
                f"allpf__ispf_cxu must include a condition axis of length {expected_conditions}; got {ispf.shape}"
            )
    else:
        errors.append(f"allpf__ispf_cxu must be 2D/3D after squeeze; got {ispf.shape}")
    if ispf.dtype.kind not in "biuf":
        errors.append(f"allpf__ispf_cxu must be numeric; got {ispf.dtype}")

    _raise_if_errors(errors, source=source, default="allcel NPZ")


def validate_allcel_npz_file(path: str | Path, *, n_condition_pairs: int = 5) -> None:
    with np.load(path, allow_pickle=False) as data:
        validate_allcel_npz(data, source=path, n_condition_pairs=n_condition_pairs)


def _check_required_columns(
    df: pd.DataFrame,
    required: Iterable[str],
    errors: list[str],
) -> None:
    missing = [col for col in required if col not in df.columns]
    if missing:
        errors.append(f"missing required column(s): {', '.join(missing)}")


def _source_feature_columns(features: Iterable[str]) -> list[str]:
    out: list[str] = []
    for feature in features:
        if feature in {"log_fr_hz", "log10_fr_hz"}:
            source = "fr_hz"
        elif feature in {"log_fr_hz_session", "log10_fr_hz_session"}:
            source = "fr_hz_session"
        else:
            source = feature
        if source not in out:
            out.append(source)
    return out


def _check_identifier_columns(df: pd.DataFrame, errors: list[str]) -> None:
    if not set(IDENTIFIER_COLUMNS).issubset(df.columns):
        return
    for col in IDENTIFIER_COLUMNS:
        if df[col].isna().any():
            errors.append(f"{col} contains missing values")
    dup = df.duplicated(list(IDENTIFIER_COLUMNS))
    if dup.any():
        errors.append(
            "duplicate (session_id, cell_id) rows found: "
            f"{int(dup.sum())} duplicate row(s)"
        )


def _check_numeric_columns(
    df: pd.DataFrame,
    columns: Iterable[str],
    errors: list[str],
) -> None:
    for col in columns:
        if col not in df.columns:
            continue
        converted = pd.to_numeric(df[col], errors="coerce")
        bad = converted.isna() & df[col].notna()
        if bad.any():
            errors.append(f"{col} contains non-numeric value(s)")


def validate_feature_table(
    df: pd.DataFrame,
    *,
    source: str | Path | None = None,
    required_features: Iterable[str] | None = None,
) -> None:
    """Validate one processed per-session feature table."""
    errors: list[str] = []
    if df.empty:
        errors.append("feature table is empty")
    _check_required_columns(df, FEATURE_TABLE_REQUIRED_COLUMNS, errors)
    if required_features is not None:
        _check_required_columns(df, _source_feature_columns(required_features), errors)
    _check_identifier_columns(df, errors)
    _check_numeric_columns(df, NUMERIC_FEATURE_COLUMNS, errors)
    _raise_if_errors(errors, source=source, default="processed feature table")


def validate_age_group_table(
    df: pd.DataFrame,
    *,
    source: str | Path | None = None,
    age_group: str | None = None,
    required_features: Iterable[str] | None = None,
    require_type_u: bool = False,
) -> None:
    """Validate an age-group table consumed by model scripts."""
    errors: list[str] = []
    if df.empty:
        errors.append("age-group table is empty")
    _check_required_columns(df, AGE_GROUP_TABLE_REQUIRED_COLUMNS, errors)
    if require_type_u:
        _check_required_columns(df, ("allcel__type_u",), errors)
    if required_features is not None:
        _check_required_columns(df, _source_feature_columns(required_features), errors)

    _check_identifier_columns(df, errors)

    if "age_group" in df.columns:
        valid = set(DEFAULT_AGE_GROUPS)
        observed = set(df["age_group"].dropna().astype(str))
        invalid = sorted(observed - valid)
        if invalid:
            errors.append(
                "age_group contains unsupported value(s): "
                f"{', '.join(invalid)}; expected one of {', '.join(DEFAULT_AGE_GROUPS)}"
            )
        if age_group is not None:
            wrong = df["age_group"].notna() & (df["age_group"].astype(str) != age_group)
            if wrong.any():
                errors.append(
                    f"age_group column must be {age_group!r}; found {int(wrong.sum())} mismatched row(s)"
                )

    numeric_cols = list(NUMERIC_FEATURE_COLUMNS)
    if "Age" in df.columns:
        numeric_cols.append("Age")
    if required_features is not None:
        numeric_cols.extend(_source_feature_columns(required_features))
    _check_numeric_columns(df, numeric_cols, errors)
    _raise_if_errors(errors, source=source, default="age-group table")
