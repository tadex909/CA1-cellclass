from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PlaceFieldConfig:
    bin_size_cm: float = 2.0
    min_occupancy_s: float = 0.1
    smooth_sigma_bins: float = 1.5
    field_threshold_ratio: float = 0.2
    min_field_bins: int = 3

    def validate(self) -> None:
        if self.bin_size_cm <= 0:
            raise ValueError("bin_size_cm must be > 0")
        if self.min_occupancy_s < 0:
            raise ValueError("min_occupancy_s must be >= 0")
        if self.smooth_sigma_bins < 0:
            raise ValueError("smooth_sigma_bins must be >= 0")
        if not (0 < self.field_threshold_ratio <= 1):
            raise ValueError("field_threshold_ratio must be in (0, 1]")
        if self.min_field_bins < 1:
            raise ValueError("min_field_bins must be >= 1")
