from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .config import PlaceFieldConfig
from .metrics import (
    compute_occupancy_s,
    compute_rate_map_hz,
    compute_spatial_information_bits_per_spike,
    compute_spike_counts,
    detect_place_fields,
    make_bin_edges,
    smooth_rate_map_hz,
)


def load_session_placefield_arrays(npz_path: Path) -> dict[str, np.ndarray]:
    """
    Load session-level arrays needed for place-field analysis.

    Expected keys:
    - `position_cm`: shape (n_samples,)
    - `dt_s`: shape (n_samples,)
    - `spike_positions_cm`: shape (n_spikes,)
    - `spike_cell_ids`: shape (n_spikes,)
    - `cell_ids`: shape (n_cells,)
    """
    with np.load(npz_path, allow_pickle=False) as z:
        needed = [
            "position_cm",
            "dt_s",
            "spike_positions_cm",
            "spike_cell_ids",
            "cell_ids",
        ]
        missing = [k for k in needed if k not in z]
        if missing:
            raise KeyError(f"Missing keys in {npz_path}: {missing}")

        return {
            "position_cm": np.asarray(z["position_cm"], dtype=np.float64).ravel(),
            "dt_s": np.asarray(z["dt_s"], dtype=np.float64).ravel(),
            "spike_positions_cm": np.asarray(z["spike_positions_cm"], dtype=np.float64).ravel(),
            "spike_cell_ids": np.asarray(z["spike_cell_ids"], dtype=np.int64).ravel(),
            "cell_ids": np.asarray(z["cell_ids"], dtype=np.int64).ravel(),
        }


def analyze_session_placefields(
    arrays: dict[str, np.ndarray],
    cfg: PlaceFieldConfig,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    cfg.validate()

    position_cm = arrays["position_cm"]
    dt_s = arrays["dt_s"]
    spike_positions_cm = arrays["spike_positions_cm"]
    spike_cell_ids = arrays["spike_cell_ids"]
    cell_ids = arrays["cell_ids"]

    valid_pos = np.isfinite(position_cm)
    if not np.any(valid_pos):
        raise ValueError("No finite position samples in session")

    track_min = float(np.nanmin(position_cm[valid_pos]))
    track_max = float(np.nanmax(position_cm[valid_pos]))
    if not np.isfinite(track_min) or not np.isfinite(track_max) or track_max <= track_min:
        raise ValueError("Invalid track bounds from position_cm")

    bin_edges = make_bin_edges(track_min, track_max, cfg.bin_size_cm)
    occupancy_s = compute_occupancy_s(position_cm, dt_s, bin_edges)

    rows: list[dict[str, Any]] = []
    for cell_id in cell_ids:
        mask = spike_cell_ids == int(cell_id)
        counts = compute_spike_counts(spike_positions_cm[mask], bin_edges)
        rate_raw = compute_rate_map_hz(
            counts,
            occupancy_s,
            min_occupancy_s=cfg.min_occupancy_s,
        )
        rate_smooth = smooth_rate_map_hz(rate_raw, cfg.smooth_sigma_bins)
        fields = detect_place_fields(
            rate_smooth,
            threshold_ratio=cfg.field_threshold_ratio,
            min_bins=cfg.min_field_bins,
        )
        info_bits = compute_spatial_information_bits_per_spike(rate_smooth, occupancy_s)

        rows.append(
            {
                "cell_id": int(cell_id),
                "pf_n_fields": int(len(fields)),
                "pf_has_field": bool(len(fields) > 0),
                "pf_peak_rate_hz": (
                    float(np.nanmax(rate_smooth))
                    if np.isfinite(rate_smooth).any()
                    else np.nan
                ),
                "pf_spatial_info_bits_per_spike": (
                    float(info_bits) if info_bits is not None else np.nan
                ),
                "pf_field_bins": repr(fields),
            }
        )

    summary = {
        "n_cells": int(cell_ids.size),
        "n_position_samples": int(position_cm.size),
        "n_spikes": int(spike_positions_cm.size),
        "track_min_cm": track_min,
        "track_max_cm": track_max,
        "n_bins": int(bin_edges.size - 1),
        "bin_size_cm": float(cfg.bin_size_cm),
    }
    return pd.DataFrame(rows), summary
