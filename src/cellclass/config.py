from __future__ import annotations

"""Shared configuration and label conventions for the cellclass workflow."""

from collections.abc import Iterable

import numpy as np
import pandas as pd


DEFAULT_AGE_GROUPS: tuple[str, ...] = ("P16-18", "P19-21", "P22-24")

AGE_GROUP_BINS: tuple[tuple[str, tuple[int, ...]], ...] = (
    ("P16-18", (16, 17, 18)),
    ("P19-21", (19, 20, 21)),
    ("P22-24", (22, 23, 24)),
)

DEFAULT_MODEL_FEATURES: tuple[str, ...] = (
    "fr_hz",
    "burst_index",
    "cv2",
    "spk_duration_ms",
    "spk_peaktrough_ms",
    "spk_asymmetry",
    "refractory_ms_edge",
    "acg_peak_latency_ms",
)

DEFAULT_AGE_AGGREGATION_FEATURES: tuple[str, ...] = (
    *DEFAULT_MODEL_FEATURES,
    "allcel__sm_u_any",
    "allcel__sm_u_task1",
    "allcel__sm_u_task2",
    "allcel__sm_u_task3",
    "allcel__sm_u_task4",
    "allcel__sm_u_task5",
)

DEFAULT_MOUSE_AGGREGATION_FEATURES: tuple[str, ...] = (
    "fr_hz",
    "burst_index",
    "cv2",
    "spk_duration_ms",
    "spk_peaktrough_ms",
    "spk_asymmetry",
    "refractory_ms_center",
    "acg_peak_latency_ms",
    "type_u",
)

FEATURE_SET_EXPERIMENTS: dict[str, tuple[str, ...]] = {
    "all_features": DEFAULT_MODEL_FEATURES,
    "some_features": (
        "fr_hz",
        "cv2",
        "acg_peak_latency_ms",
        "spk_duration_ms",
        "spk_asymmetry",
    ),
    "few_features": (
        "fr_hz",
        "spk_duration_ms",
        "refractory_ms_edge",
    ),
    "valero_features": (
        "cv2",
        "acg_peak_latency_ms",
        "spk_duration_ms",
        "spk_asymmetry",
        "fr_hz",
    ),
}

TYPE_U_INTERNEURON = 0
TYPE_U_PYRAMIDAL = 1
CELL_TYPE_INTERNEURON = "interneuron"
CELL_TYPE_PYRAMIDAL = "pyramidal"

TYPE_U_TO_CELL_TYPE: dict[int, str] = {
    TYPE_U_INTERNEURON: CELL_TYPE_INTERNEURON,
    TYPE_U_PYRAMIDAL: CELL_TYPE_PYRAMIDAL,
}

TYPE_U_TEXT_TO_BINARY: dict[str, int] = {
    "0": TYPE_U_INTERNEURON,
    CELL_TYPE_INTERNEURON: TYPE_U_INTERNEURON,
    "int": TYPE_U_INTERNEURON,
    "1": TYPE_U_PYRAMIDAL,
    CELL_TYPE_PYRAMIDAL: TYPE_U_PYRAMIDAL,
    "pyr": TYPE_U_PYRAMIDAL,
}


def parse_csv_list(raw: str | None) -> list[str]:
    if not raw:
        return []
    return [x.strip() for x in raw.split(",") if x.strip()]


def csv_join(values: Iterable[str]) -> str:
    return ",".join(values)


def age_group_from_age(age: int) -> str | None:
    for group, ages in AGE_GROUP_BINS:
        if age in ages:
            return group
    return None


def normalize_type_u(series: pd.Series) -> pd.Series:
    """
    Standardize external type_u labels to nullable integers:
      0 = interneuron
      1 = pyramidal
    """
    s = series.copy()
    if pd.api.types.is_numeric_dtype(s):
        out = pd.to_numeric(s, errors="coerce")
        out = out.where(out.isin([TYPE_U_INTERNEURON, TYPE_U_PYRAMIDAL]), np.nan)
        return out.astype("Int64")

    out = s.astype(str).str.strip().str.lower().map(TYPE_U_TEXT_TO_BINARY)
    return out.astype("Int64")


def type_u_binary_to_name(series: pd.Series) -> pd.Series:
    return series.map(TYPE_U_TO_CELL_TYPE).astype("string")
