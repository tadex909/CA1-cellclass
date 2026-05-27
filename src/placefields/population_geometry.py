from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
import math
from typing import Any, Mapping

import numpy as np

from .cue_zones import CueZoneLayout, cue_zone_layout_for_condition, label_xbin_centers_by_zone
from .interim_io import SavedRatemapPack, subset_saved_ratemap_pack_cells
from .trials import condition_family_name


VALID_NORMALIZATIONS = ("raw", "mean_rate")


@dataclass(frozen=True)
class DisplacementProfile:
    delta_bin: np.ndarray
    delta_x: np.ndarray
    mean: np.ndarray
    n_pairs: np.ndarray


@dataclass(frozen=True)
class ZoneDisplacementProfiles:
    zone_labels_k: np.ndarray
    all_profile: DisplacementProfile
    anchor_profiles: dict[str, DisplacementProfile]
    within_profiles: dict[str, DisplacementProfile]
    layout: CueZoneLayout | None = None


@dataclass(frozen=True)
class PopulationGeometryConfig:
    n_geom_bins: int = 10
    min_occupancy_s: float = 0.05
    normalizations: tuple[str, ...] = ("raw", "mean_rate")
    n_splits: int = 100
    max_exact_splits: int = 128
    seed: int | None = None

    def validate(self) -> None:
        if self.n_geom_bins < 1:
            raise ValueError("n_geom_bins must be >= 1")
        if not np.isfinite(self.min_occupancy_s) or self.min_occupancy_s < 0:
            raise ValueError("min_occupancy_s must be finite and >= 0")
        if self.n_splits < 1:
            raise ValueError("n_splits must be >= 1")
        if self.max_exact_splits < 1:
            raise ValueError("max_exact_splits must be >= 1")

        seen: set[str] = set()
        norms: list[str] = []
        for raw in self.normalizations:
            norm = str(raw).strip()
            if not norm:
                continue
            if norm not in VALID_NORMALIZATIONS:
                raise ValueError(
                    f"Unsupported normalization {norm!r}. Expected one of {VALID_NORMALIZATIONS}."
                )
            if norm in seen:
                continue
            seen.add(norm)
            norms.append(norm)
        if not norms:
            raise ValueError("normalizations must contain at least one supported option")
        object.__setattr__(self, "normalizations", tuple(norms))


@dataclass(frozen=True)
class PopulationGeometryGroupResult:
    condway_1b: int
    base_condition_1b: int
    direction: str
    condition_name: str
    condition_family: str
    n_laps: int
    n_splits: int
    cv_status: str
    G_nkk: np.ndarray
    S_cv_nkk: np.ndarray
    valid_mean_nk: np.ndarray
    valid_cv_nk: np.ndarray


@dataclass(frozen=True)
class PopulationGeometrySessionResult:
    xbin_edges: np.ndarray
    xbin_centers: np.ndarray
    cell_ids_u: np.ndarray
    normalizations: tuple[str, ...]
    condway_1b_g: np.ndarray
    base_condition_1b_g: np.ndarray
    direction_g: np.ndarray
    condition_name_g: np.ndarray
    condition_family_g: np.ndarray
    G_gnkk: np.ndarray
    S_cv_gnkk: np.ndarray
    valid_mean_gnk: np.ndarray
    valid_cv_gnk: np.ndarray
    n_laps_g: np.ndarray
    n_splits_g: np.ndarray
    cv_status_g: np.ndarray

    def to_payload(self) -> dict[str, Any]:
        return {
            "xbin_edges": self.xbin_edges.astype(np.float64, copy=False),
            "xbin_centers": self.xbin_centers.astype(np.float64, copy=False),
            "geom__cell_ids_u": self.cell_ids_u.astype(np.int64, copy=False),
            "geom__normalizations_n": np.asarray(self.normalizations, dtype=np.str_),
            "geom__condway_1b_g": self.condway_1b_g.astype(np.int64, copy=False),
            "geom__base_condition_1b_g": self.base_condition_1b_g.astype(np.int64, copy=False),
            "geom__direction_g": self.direction_g.astype(np.str_),
            "geom__condition_name_g": self.condition_name_g.astype(np.str_),
            "geom__condition_family_g": self.condition_family_g.astype(np.str_),
            "geom__G_gnkk": self.G_gnkk.astype(np.float32, copy=False),
            "geom__S_cv_gnkk": self.S_cv_gnkk.astype(np.float32, copy=False),
            "geom__valid_mean_gnk": self.valid_mean_gnk.astype(bool, copy=False),
            "geom__valid_cv_gnk": self.valid_cv_gnk.astype(bool, copy=False),
            "geom__n_laps_g": self.n_laps_g.astype(np.int64, copy=False),
            "geom__n_splits_g": self.n_splits_g.astype(np.int64, copy=False),
            "geom__cv_status_g": self.cv_status_g.astype(np.str_),
        }

    def summary_rows(self, *, session_id: str) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        for g in range(int(self.condway_1b_g.size)):
            for n, norm in enumerate(self.normalizations):
                rows.append(
                    {
                        "session_id": str(session_id),
                        "condway_1b": int(self.condway_1b_g[g]),
                        "base_condition_1b": int(self.base_condition_1b_g[g]),
                        "direction": str(self.direction_g[g]),
                        "condition_name": str(self.condition_name_g[g]),
                        "condition_family": str(self.condition_family_g[g]),
                        "normalization": str(norm),
                        "n_cells_used": int(self.cell_ids_u.size),
                        "n_laps": int(self.n_laps_g[g]),
                        "n_splits": int(self.n_splits_g[g]),
                        "cv_status": str(self.cv_status_g[g]),
                        "n_geom_bins": int(self.xbin_centers.size),
                        "n_valid_mean_bins": int(np.sum(self.valid_mean_gnk[g, n])),
                        "n_valid_cv_bins": int(np.sum(self.valid_cv_gnk[g, n])),
                    }
                )
        return rows


def decode_condway(condway: int) -> tuple[int, str]:
    c = int(condway)
    if c < 1:
        raise ValueError(f"condway must be >= 1, got {c}")
    direction = "W" if (c % 2) == 1 else "B"
    base_condition_1b = (c + 1) // 2
    return base_condition_1b, direction


def rebin_trial_maps(
    *,
    counts_lkn: np.ndarray,
    dwell_lk: np.ndarray,
    xbin_edges: np.ndarray,
    n_geom_bins: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    counts = np.asarray(counts_lkn, dtype=np.float64)
    dwell = np.asarray(dwell_lk, dtype=np.float64)
    edges = np.asarray(xbin_edges, dtype=np.float64).ravel()

    if counts.ndim != 3:
        raise ValueError(f"counts_lkn must be 3D (lap, bin, cell), got ndim={counts.ndim}")
    if dwell.ndim != 2:
        raise ValueError(f"dwell_lk must be 2D (lap, bin), got ndim={dwell.ndim}")
    if counts.shape[:2] != dwell.shape:
        raise ValueError(f"counts/dwell shape mismatch: {counts.shape[:2]} vs {dwell.shape}")

    n_laps, n_bins_fine, n_cells = counts.shape
    if edges.size != (n_bins_fine + 1):
        raise ValueError(
            f"xbin_edges length {edges.size} must equal n_bins + 1 = {n_bins_fine + 1}"
        )
    if n_geom_bins < 1:
        raise ValueError("n_geom_bins must be >= 1")

    n_bins_geom = min(int(n_geom_bins), int(n_bins_fine))
    groups = [np.asarray(idx, dtype=np.int64) for idx in np.array_split(np.arange(n_bins_fine), n_bins_geom)]
    if any(g.size == 0 for g in groups):
        raise RuntimeError("Internal error: empty geometry bin after rebinning")

    counts_geom = np.zeros((n_laps, n_bins_geom, n_cells), dtype=np.float64)
    dwell_geom = np.zeros((n_laps, n_bins_geom), dtype=np.float64)
    edges_geom = np.empty(n_bins_geom + 1, dtype=np.float64)

    for k, group in enumerate(groups):
        counts_geom[:, k, :] = np.sum(counts[:, group, :], axis=1)
        dwell_geom[:, k] = np.sum(dwell[:, group], axis=1)
        edges_geom[k] = float(edges[int(group[0])])
    edges_geom[-1] = float(edges[int(groups[-1][-1]) + 1])

    centers_geom = 0.5 * (edges_geom[:-1] + edges_geom[1:])
    return counts_geom, dwell_geom, edges_geom, centers_geom


def compute_normalization_scales(
    *,
    counts_lkn: np.ndarray,
    dwell_lk: np.ndarray,
    min_occupancy_s: float,
    normalizations: tuple[str, ...],
) -> dict[str, np.ndarray]:
    counts = np.asarray(counts_lkn, dtype=np.float64)
    dwell = np.asarray(dwell_lk, dtype=np.float64)
    if counts.ndim != 3 or dwell.ndim != 2 or counts.shape[:2] != dwell.shape:
        raise ValueError("counts_lkn and dwell_lk must have aligned (lap, bin) axes")

    _, _, n_cells = counts.shape
    scales: dict[str, np.ndarray] = {
        "raw": np.ones(n_cells, dtype=np.float64),
    }

    if "mean_rate" in normalizations:
        rates = np.divide(
            counts,
            dwell[:, :, None],
            out=np.full_like(counts, np.nan, dtype=np.float64),
            where=dwell[:, :, None] > 0,
        )
        valid_lk = np.isfinite(dwell) & (dwell >= float(min_occupancy_s))
        valid_count = int(np.sum(valid_lk))
        if valid_count > 0:
            scale = (
                np.sum(np.where(valid_lk[:, :, None], rates, 0.0), axis=(0, 1))
                / float(valid_count)
            )
        else:
            scale = np.full(n_cells, np.nan, dtype=np.float64)
        scale = np.where(np.isfinite(scale) & (scale > 0), scale, 1.0)
        scales["mean_rate"] = scale.astype(np.float64, copy=False)
    return scales


def _apply_neuron_scale(rates_kn: np.ndarray, scale_n: np.ndarray) -> np.ndarray:
    scale = np.asarray(scale_n, dtype=np.float64).ravel()
    if rates_kn.shape[1] != scale.size:
        raise ValueError(f"scale length {scale.size} must match n_cells {rates_kn.shape[1]}")
    safe_scale = np.where(np.isfinite(scale) & (scale > 0), scale, 1.0)
    return rates_kn / safe_scale[None, :]


def _pooled_rates_from_laps(
    counts_lkn: np.ndarray,
    dwell_lk: np.ndarray,
    lap_idx: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    pooled_counts = np.sum(counts_lkn[lap_idx, :, :], axis=0)
    pooled_dwell = np.sum(dwell_lk[lap_idx, :], axis=0)
    pooled_rates = np.divide(
        pooled_counts,
        pooled_dwell[:, None],
        out=np.full_like(pooled_counts, np.nan, dtype=np.float64),
        where=pooled_dwell[:, None] > 0,
    )
    return pooled_rates, pooled_dwell


def _valid_vector_mask(
    rates_kn: np.ndarray,
    pooled_dwell_k: np.ndarray,
    *,
    min_occupancy_s: float,
) -> np.ndarray:
    row_finite = np.all(np.isfinite(rates_kn), axis=1)
    norms = np.linalg.norm(np.where(np.isfinite(rates_kn), rates_kn, 0.0), axis=1)
    return (
        np.isfinite(pooled_dwell_k)
        & (pooled_dwell_k >= float(min_occupancy_s))
        & row_finite
        & np.isfinite(norms)
        & (norms > 0)
    )


def _cosine_matrix(
    a_kn: np.ndarray,
    b_kn: np.ndarray,
    valid_a_k: np.ndarray,
    valid_b_k: np.ndarray,
) -> np.ndarray:
    a = np.asarray(a_kn, dtype=np.float64)
    b = np.asarray(b_kn, dtype=np.float64)
    valid_a = np.asarray(valid_a_k, dtype=bool).ravel()
    valid_b = np.asarray(valid_b_k, dtype=bool).ravel()
    if a.ndim != 2 or b.ndim != 2:
        raise ValueError("a_kn and b_kn must be 2D (bin, cell)")
    if a.shape[1] != b.shape[1]:
        raise ValueError(f"cell-axis mismatch: {a.shape[1]} vs {b.shape[1]}")
    if valid_a.size != a.shape[0] or valid_b.size != b.shape[0]:
        raise ValueError("valid mask length must match number of bins")

    out = np.full((a.shape[0], b.shape[0]), np.nan, dtype=np.float64)
    idx_a = np.flatnonzero(valid_a)
    idx_b = np.flatnonzero(valid_b)
    if idx_a.size == 0 or idx_b.size == 0:
        return out

    a_sel = a[idx_a, :]
    b_sel = b[idx_b, :]
    a_norm = np.linalg.norm(a_sel, axis=1)
    b_norm = np.linalg.norm(b_sel, axis=1)
    a_unit = a_sel / a_norm[:, None]
    b_unit = b_sel / b_norm[:, None]
    block = np.clip(a_unit @ b_unit.T, -1.0, 1.0)
    out[np.ix_(idx_a, idx_b)] = block
    return out


def _canonical_subset_key(subset: tuple[int, ...], *, n_laps: int) -> tuple[int, ...]:
    if (n_laps % 2) == 1 or 0 in subset:
        return subset
    subset_set = set(subset)
    complement = tuple(i for i in range(n_laps) if i not in subset_set)
    return complement


def _enumerate_exact_splits(n_laps: int) -> list[tuple[np.ndarray, np.ndarray]]:
    size_a = n_laps // 2
    splits: list[tuple[np.ndarray, np.ndarray]] = []
    for subset in combinations(range(n_laps), size_a):
        key = _canonical_subset_key(subset, n_laps=n_laps)
        if key != subset:
            continue
        subset_set = set(subset)
        other = tuple(i for i in range(n_laps) if i not in subset_set)
        splits.append(
            (
                np.asarray(subset, dtype=np.int64),
                np.asarray(other, dtype=np.int64),
            )
        )
    return splits


def _sample_splits(
    *,
    n_laps: int,
    n_splits: int,
    seed: int | None,
) -> list[tuple[np.ndarray, np.ndarray]]:
    size_a = n_laps // 2
    target = int(n_splits)
    rng = np.random.default_rng(seed)
    seen: set[tuple[int, ...]] = set()
    out: list[tuple[np.ndarray, np.ndarray]] = []
    max_attempts = max(200, 50 * target)

    attempts = 0
    while len(out) < target and attempts < max_attempts:
        subset = tuple(sorted(rng.choice(n_laps, size=size_a, replace=False).tolist()))
        key = _canonical_subset_key(subset, n_laps=n_laps)
        if key in seen:
            attempts += 1
            continue
        seen.add(key)
        key_set = set(key)
        other = tuple(i for i in range(n_laps) if i not in key_set)
        out.append(
            (
                np.asarray(key, dtype=np.int64),
                np.asarray(other, dtype=np.int64),
            )
        )
        attempts += 1
    return out


def _build_group_splits(cfg: PopulationGeometryConfig, n_laps: int) -> tuple[list[tuple[np.ndarray, np.ndarray]], str]:
    if n_laps < 2:
        return [], "too_few_laps_for_cv"

    total_splits = math.comb(n_laps, n_laps // 2)
    if (n_laps % 2) == 0:
        total_splits //= 2

    if total_splits <= int(cfg.max_exact_splits):
        return _enumerate_exact_splits(n_laps), "ok_exact"

    target = min(int(cfg.n_splits), int(total_splits))
    splits = _sample_splits(n_laps=n_laps, n_splits=target, seed=cfg.seed)
    return splits, "ok_sampled"


def _coerce_bin_mask(mask_k: np.ndarray | None, *, n_bins: int, name: str) -> np.ndarray:
    if mask_k is None:
        return np.ones(n_bins, dtype=bool)
    mask = np.asarray(mask_k, dtype=bool).ravel()
    if mask.size != n_bins:
        raise ValueError(f"{name} length {mask.size} must match n_bins {n_bins}")
    return mask


def compute_displacement_profile(
    similarity_kk: np.ndarray,
    *,
    xbin_centers: np.ndarray | None = None,
    valid_mask_k: np.ndarray | None = None,
    anchor_mask_k: np.ndarray | None = None,
    target_mask_k: np.ndarray | None = None,
    max_lag_bins: int | None = None,
) -> DisplacementProfile:
    sim = np.asarray(similarity_kk, dtype=np.float64)
    if sim.ndim != 2 or sim.shape[0] != sim.shape[1]:
        raise ValueError(
            f"similarity_kk must be a square 2D matrix, got shape {sim.shape}"
        )

    n_bins = int(sim.shape[0])
    if n_bins < 2:
        return DisplacementProfile(
            delta_bin=np.zeros(0, dtype=np.int64),
            delta_x=np.zeros(0, dtype=np.float64),
            mean=np.zeros(0, dtype=np.float64),
            n_pairs=np.zeros(0, dtype=np.int64),
        )

    if xbin_centers is None:
        centers = None
    else:
        centers = np.asarray(xbin_centers, dtype=np.float64).ravel()
        if centers.size != n_bins:
            raise ValueError(
                f"xbin_centers length {centers.size} must match n_bins {n_bins}"
            )

    valid_mask = _coerce_bin_mask(valid_mask_k, n_bins=n_bins, name="valid_mask_k")
    anchor_mask = _coerce_bin_mask(anchor_mask_k, n_bins=n_bins, name="anchor_mask_k")
    target_mask = _coerce_bin_mask(target_mask_k, n_bins=n_bins, name="target_mask_k")

    max_lag = n_bins - 1 if max_lag_bins is None else min(int(max_lag_bins), n_bins - 1)
    if max_lag < 1:
        raise ValueError("max_lag_bins must be >= 1 when provided")

    delta_bin = np.arange(1, max_lag + 1, dtype=np.int64)
    delta_x = np.full(max_lag, np.nan, dtype=np.float64)
    mean = np.full(max_lag, np.nan, dtype=np.float64)
    n_pairs = np.zeros(max_lag, dtype=np.int64)

    for lag in delta_bin.tolist():
        i = np.arange(0, n_bins - lag, dtype=np.int64)
        j = i + lag
        if centers is None:
            delta_x[lag - 1] = float(lag)
        else:
            delta_x[lag - 1] = float(np.mean(centers[j] - centers[i]))

        pair_mask = (
            valid_mask[i]
            & valid_mask[j]
            & anchor_mask[i]
            & target_mask[j]
            & np.isfinite(sim[i, j])
        )
        if not np.any(pair_mask):
            continue

        values = sim[i[pair_mask], j[pair_mask]]
        mean[lag - 1] = float(np.mean(values))
        n_pairs[lag - 1] = int(values.size)

    return DisplacementProfile(
        delta_bin=delta_bin,
        delta_x=delta_x,
        mean=mean,
        n_pairs=n_pairs,
    )


def compute_zone_displacement_profiles(
    similarity_kk: np.ndarray,
    *,
    zone_labels_k: np.ndarray,
    xbin_centers: np.ndarray | None = None,
    valid_mask_k: np.ndarray | None = None,
    zone_names: tuple[str, ...] = ("rich", "poor"),
    max_lag_bins: int | None = None,
    layout: CueZoneLayout | None = None,
) -> ZoneDisplacementProfiles:
    sim = np.asarray(similarity_kk, dtype=np.float64)
    if sim.ndim != 2 or sim.shape[0] != sim.shape[1]:
        raise ValueError(
            f"similarity_kk must be a square 2D matrix, got shape {sim.shape}"
        )

    labels = np.asarray(zone_labels_k, dtype=np.str_).ravel()
    n_bins = int(sim.shape[0])
    if labels.size != n_bins:
        raise ValueError(f"zone_labels_k length {labels.size} must match n_bins {n_bins}")

    all_profile = compute_displacement_profile(
        sim,
        xbin_centers=xbin_centers,
        valid_mask_k=valid_mask_k,
        max_lag_bins=max_lag_bins,
    )

    anchor_profiles: dict[str, DisplacementProfile] = {}
    within_profiles: dict[str, DisplacementProfile] = {}
    for zone_name in zone_names:
        zone_mask = labels == str(zone_name)
        anchor_profiles[str(zone_name)] = compute_displacement_profile(
            sim,
            xbin_centers=xbin_centers,
            valid_mask_k=valid_mask_k,
            anchor_mask_k=zone_mask,
            max_lag_bins=max_lag_bins,
        )
        within_profiles[str(zone_name)] = compute_displacement_profile(
            sim,
            xbin_centers=xbin_centers,
            valid_mask_k=valid_mask_k,
            anchor_mask_k=zone_mask,
            target_mask_k=zone_mask,
            max_lag_bins=max_lag_bins,
        )

    return ZoneDisplacementProfiles(
        zone_labels_k=labels,
        all_profile=all_profile,
        anchor_profiles=anchor_profiles,
        within_profiles=within_profiles,
        layout=layout,
    )


def compute_condition_zone_displacement_profiles(
    similarity_kk: np.ndarray,
    *,
    xbin_centers: np.ndarray,
    valid_mask_k: np.ndarray | None = None,
    condition_name: str | None = "",
    condition_family: str | None = "",
    max_lag_bins: int | None = None,
) -> ZoneDisplacementProfiles:
    layout = cue_zone_layout_for_condition(
        condition_name,
        condition_family=condition_family,
    )
    zone_labels = label_xbin_centers_by_zone(
        xbin_centers,
        layout=layout,
    )
    return compute_zone_displacement_profiles(
        similarity_kk,
        zone_labels_k=zone_labels,
        xbin_centers=xbin_centers,
        valid_mask_k=valid_mask_k,
        max_lag_bins=max_lag_bins,
        layout=layout,
    )


def build_population_geometry_for_group(
    *,
    counts_lkn: np.ndarray,
    dwell_lk: np.ndarray,
    cfg: PopulationGeometryConfig,
    condway_1b: int,
    condition_name: str = "",
    normalization_scales: Mapping[str, np.ndarray] | None = None,
) -> PopulationGeometryGroupResult:
    cfg.validate()

    counts = np.asarray(counts_lkn, dtype=np.float64)
    dwell = np.asarray(dwell_lk, dtype=np.float64)
    if counts.ndim != 3:
        raise ValueError(f"counts_lkn must be 3D (lap, bin, cell), got ndim={counts.ndim}")
    if dwell.ndim != 2:
        raise ValueError(f"dwell_lk must be 2D (lap, bin), got ndim={dwell.ndim}")
    if counts.shape[:2] != dwell.shape:
        raise ValueError(f"counts/dwell shape mismatch: {counts.shape[:2]} vs {dwell.shape}")

    n_laps, n_bins, n_cells = counts.shape
    if normalization_scales is None:
        scales_by_norm = compute_normalization_scales(
            counts_lkn=counts,
            dwell_lk=dwell,
            min_occupancy_s=float(cfg.min_occupancy_s),
            normalizations=cfg.normalizations,
        )
    else:
        scales_by_norm = {
            str(k): np.asarray(v, dtype=np.float64).ravel()
            for k, v in normalization_scales.items()
        }

    splits, cv_status = _build_group_splits(cfg, n_laps)
    n_norm = len(cfg.normalizations)

    G_nkk = np.full((n_norm, n_bins, n_bins), np.nan, dtype=np.float64)
    S_cv_nkk = np.full((n_norm, n_bins, n_bins), np.nan, dtype=np.float64)
    valid_mean_nk = np.zeros((n_norm, n_bins), dtype=bool)
    valid_cv_nk = np.zeros((n_norm, n_bins), dtype=bool)

    pooled_rates_kn, pooled_dwell_k = _pooled_rates_from_laps(
        counts,
        dwell,
        np.arange(n_laps, dtype=np.int64),
    )

    for n_idx, norm in enumerate(cfg.normalizations):
        scale = scales_by_norm.get(norm)
        if scale is None:
            if norm == "raw":
                scale = np.ones(n_cells, dtype=np.float64)
            else:
                raise KeyError(f"Missing normalization scale for {norm!r}")

        pooled_norm_kn = _apply_neuron_scale(pooled_rates_kn, scale)
        valid_mean = _valid_vector_mask(
            pooled_norm_kn,
            pooled_dwell_k,
            min_occupancy_s=float(cfg.min_occupancy_s),
        )
        G = _cosine_matrix(pooled_norm_kn, pooled_norm_kn, valid_mean, valid_mean)
        diag_idx = np.flatnonzero(valid_mean)
        if diag_idx.size:
            G[diag_idx, diag_idx] = 1.0
        G_nkk[n_idx] = G
        valid_mean_nk[n_idx] = valid_mean

        if not splits:
            continue

        split_stack = np.full((len(splits), n_bins, n_bins), np.nan, dtype=np.float64)
        for s_idx, (lap_a, lap_b) in enumerate(splits):
            rates_a_kn, dwell_a_k = _pooled_rates_from_laps(counts, dwell, lap_a)
            rates_b_kn, dwell_b_k = _pooled_rates_from_laps(counts, dwell, lap_b)
            rates_a_norm_kn = _apply_neuron_scale(rates_a_kn, scale)
            rates_b_norm_kn = _apply_neuron_scale(rates_b_kn, scale)
            valid_a = _valid_vector_mask(
                rates_a_norm_kn,
                dwell_a_k,
                min_occupancy_s=float(cfg.min_occupancy_s),
            )
            valid_b = _valid_vector_mask(
                rates_b_norm_kn,
                dwell_b_k,
                min_occupancy_s=float(cfg.min_occupancy_s),
            )
            mat_ab = _cosine_matrix(rates_a_norm_kn, rates_b_norm_kn, valid_a, valid_b)
            mat_ba = _cosine_matrix(rates_b_norm_kn, rates_a_norm_kn, valid_b, valid_a)
            split_stack[s_idx] = 0.5 * (mat_ab + mat_ba)

        finite = np.isfinite(split_stack)
        count = np.sum(finite, axis=0)
        summed = np.nansum(np.where(finite, split_stack, 0.0), axis=0)
        mean_mat = np.full((n_bins, n_bins), np.nan, dtype=np.float64)
        valid_entry = count > 0
        mean_mat[valid_entry] = summed[valid_entry] / count[valid_entry]
        S_cv_nkk[n_idx] = mean_mat
        valid_cv_nk[n_idx] = np.diag(count) > 0

    base_condition_1b, direction = decode_condway(condway_1b)
    return PopulationGeometryGroupResult(
        condway_1b=int(condway_1b),
        base_condition_1b=int(base_condition_1b),
        direction=str(direction),
        condition_name=str(condition_name),
        condition_family=condition_family_name(condition_name) if condition_name else "",
        n_laps=int(n_laps),
        n_splits=int(len(splits)),
        cv_status=str(cv_status),
        G_nkk=G_nkk,
        S_cv_nkk=S_cv_nkk,
        valid_mean_nk=valid_mean_nk,
        valid_cv_nk=valid_cv_nk,
    )


def build_population_geometry_from_saved_ratemap(
    *,
    saved: SavedRatemapPack,
    cfg: PopulationGeometryConfig,
    condition_names_by_base: Mapping[int, str] | None = None,
    selected_cell_ids: np.ndarray | None = None,
) -> PopulationGeometrySessionResult:
    cfg.validate()

    saved_use = saved
    if selected_cell_ids is not None:
        keep_ids = np.asarray(selected_cell_ids, dtype=np.int64).ravel()
        if keep_ids.size == 0:
            raise ValueError("selected_cell_ids is empty")
        saved_use = subset_saved_ratemap_pack_cells(saved, keep_ids)

    counts_lkn_fine = np.transpose(np.asarray(saved_use.nbspk_tx_ux, dtype=np.float64), (1, 2, 0))
    dwell_lk_fine = np.asarray(saved_use.dwell_tx_x, dtype=np.float64)
    counts_lkn, dwell_lk, edges_geom, centers_geom = rebin_trial_maps(
        counts_lkn=counts_lkn_fine,
        dwell_lk=dwell_lk_fine,
        xbin_edges=saved_use.xbin_edges,
        n_geom_bins=int(cfg.n_geom_bins),
    )
    scales_by_norm = compute_normalization_scales(
        counts_lkn=counts_lkn,
        dwell_lk=dwell_lk,
        min_occupancy_s=float(cfg.min_occupancy_s),
        normalizations=cfg.normalizations,
    )

    cond_names = condition_names_by_base or {}
    condways = np.unique(np.asarray(saved_use.idcond_t, dtype=np.int64))
    condways = condways[np.isfinite(condways)]
    condways = np.sort(condways.astype(np.int64, copy=False))
    groups: list[PopulationGeometryGroupResult] = []
    for condway in condways.tolist():
        idx = np.asarray(saved_use.idcond_t, dtype=np.int64) == int(condway)
        if not np.any(idx):
            continue
        base_condition_1b, _ = decode_condway(int(condway))
        groups.append(
            build_population_geometry_for_group(
                counts_lkn=counts_lkn[idx, :, :],
                dwell_lk=dwell_lk[idx, :],
                cfg=cfg,
                condway_1b=int(condway),
                condition_name=str(cond_names.get(int(base_condition_1b), "")),
                normalization_scales=scales_by_norm,
            )
        )

    n_groups = len(groups)
    n_norm = len(cfg.normalizations)
    n_bins = int(centers_geom.size)

    condway_1b_g = np.asarray([g.condway_1b for g in groups], dtype=np.int64)
    base_condition_1b_g = np.asarray([g.base_condition_1b for g in groups], dtype=np.int64)
    direction_g = np.asarray([g.direction for g in groups], dtype=np.str_)
    condition_name_g = np.asarray([g.condition_name for g in groups], dtype=np.str_)
    condition_family_g = np.asarray([g.condition_family for g in groups], dtype=np.str_)
    n_laps_g = np.asarray([g.n_laps for g in groups], dtype=np.int64)
    n_splits_g = np.asarray([g.n_splits for g in groups], dtype=np.int64)
    cv_status_g = np.asarray([g.cv_status for g in groups], dtype=np.str_)

    G_gnkk = np.full((n_groups, n_norm, n_bins, n_bins), np.nan, dtype=np.float64)
    S_cv_gnkk = np.full((n_groups, n_norm, n_bins, n_bins), np.nan, dtype=np.float64)
    valid_mean_gnk = np.zeros((n_groups, n_norm, n_bins), dtype=bool)
    valid_cv_gnk = np.zeros((n_groups, n_norm, n_bins), dtype=bool)

    for g_idx, group in enumerate(groups):
        G_gnkk[g_idx] = group.G_nkk
        S_cv_gnkk[g_idx] = group.S_cv_nkk
        valid_mean_gnk[g_idx] = group.valid_mean_nk
        valid_cv_gnk[g_idx] = group.valid_cv_nk

    return PopulationGeometrySessionResult(
        xbin_edges=edges_geom,
        xbin_centers=centers_geom,
        cell_ids_u=np.asarray(saved_use.cell_ids, dtype=np.int64),
        normalizations=cfg.normalizations,
        condway_1b_g=condway_1b_g,
        base_condition_1b_g=base_condition_1b_g,
        direction_g=direction_g,
        condition_name_g=condition_name_g,
        condition_family_g=condition_family_g,
        G_gnkk=G_gnkk,
        S_cv_gnkk=S_cv_gnkk,
        valid_mean_gnk=valid_mean_gnk,
        valid_cv_gnk=valid_cv_gnk,
        n_laps_g=n_laps_g,
        n_splits_g=n_splits_g,
        cv_status_g=cv_status_g,
    )
