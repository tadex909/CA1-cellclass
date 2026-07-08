from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .trials import canonical_condition_name, condition_family_name


DEFAULT_TRACK_START = 0.0
DEFAULT_TRACK_END = 100.0


@dataclass(frozen=True)
class CueZoneLayout:
    """
    Canonical cue-zone layout on the normalized 0-100 track.

    Notes
    -----
    - `rich_spans` are stored as half-open intervals `[start, end)`.
    - `excluded_spans` are omitted from both rich and poor zone labeling.
    - `poor_spans` are computed as the complement of `rich_spans ∪ excluded_spans`
      on the track.
    - `PNO` has no cue-rich spans by definition.
    """

    condition_name: str
    condition_family: str
    track_start: float
    track_end: float
    rich_spans: tuple[tuple[float, float], ...]
    excluded_spans: tuple[tuple[float, float], ...]
    object_centers: tuple[float, ...]

    @property
    def poor_spans(self) -> tuple[tuple[float, float], ...]:
        return complement_spans(
            tuple(self.rich_spans) + tuple(self.excluded_spans),
            track_start=float(self.track_start),
            track_end=float(self.track_end),
        )


def complement_spans(
    spans: tuple[tuple[float, float], ...] | list[tuple[float, float]],
    *,
    track_start: float = DEFAULT_TRACK_START,
    track_end: float = DEFAULT_TRACK_END,
) -> tuple[tuple[float, float], ...]:
    lo = float(track_start)
    hi = float(track_end)
    if hi <= lo:
        return tuple()

    clipped: list[tuple[float, float]] = []
    for start, end in sorted((float(a), float(b)) for a, b in spans):
        start_c = max(lo, min(hi, start))
        end_c = max(lo, min(hi, end))
        if end_c <= start_c:
            continue
        clipped.append((start_c, end_c))

    if not clipped:
        return ((lo, hi),)

    out: list[tuple[float, float]] = []
    cursor = lo
    for start_c, end_c in clipped:
        if start_c > cursor:
            out.append((cursor, start_c))
        cursor = max(cursor, end_c)
    if cursor < hi:
        out.append((cursor, hi))
    return tuple(out)


def _span_mask(values: np.ndarray, start: float, end: float, track_end: float) -> np.ndarray:
    vals = np.asarray(values, dtype=np.float64)
    stop = float(end)
    if np.isclose(stop, float(track_end)):
        return np.isfinite(vals) & (vals >= float(start)) & (vals <= stop)
    return np.isfinite(vals) & (vals >= float(start)) & (vals < stop)


def _track_bounds_from_xbin_edges(xbin_edges: np.ndarray) -> tuple[float, float]:
    edges = np.asarray(xbin_edges, dtype=np.float64).ravel()
    if edges.size < 2:
        raise ValueError("xbin_edges must contain at least 2 edges")
    if not np.all(np.isfinite(edges)):
        raise ValueError("xbin_edges must be finite")
    if not np.all(np.diff(edges) > 0):
        raise ValueError("xbin_edges must be strictly increasing")
    return float(edges[0]), float(edges[-1])


def _scale_value(value: float, *, source_start: float, source_end: float, track_start: float, track_end: float) -> float:
    scale = (float(track_end) - float(track_start)) / (float(source_end) - float(source_start))
    return float(track_start) + (float(value) - float(source_start)) * scale


def _scale_spans(
    spans: tuple[tuple[float, float], ...],
    *,
    source_start: float,
    source_end: float,
    track_start: float,
    track_end: float,
) -> tuple[tuple[float, float], ...]:
    return tuple(
        (
            _scale_value(start, source_start=source_start, source_end=source_end, track_start=track_start, track_end=track_end),
            _scale_value(end, source_start=source_start, source_end=source_end, track_start=track_start, track_end=track_end),
        )
        for start, end in spans
    )


def scale_cue_zone_layout(
    layout: CueZoneLayout,
    *,
    track_start: float,
    track_end: float,
) -> CueZoneLayout:
    """
    Scale a canonical cue-zone layout onto a target track coordinate system.
    """

    source_start = float(layout.track_start)
    source_end = float(layout.track_end)
    target_start = float(track_start)
    target_end = float(track_end)
    if not np.isfinite(target_start) or not np.isfinite(target_end) or target_end <= target_start:
        raise ValueError("track_start/track_end must be finite with track_end > track_start")
    if not np.isfinite(source_start) or not np.isfinite(source_end) or source_end <= source_start:
        raise ValueError("layout track_start/track_end must be finite with track_end > track_start")

    return CueZoneLayout(
        condition_name=layout.condition_name,
        condition_family=layout.condition_family,
        track_start=target_start,
        track_end=target_end,
        rich_spans=_scale_spans(
            layout.rich_spans,
            source_start=source_start,
            source_end=source_end,
            track_start=target_start,
            track_end=target_end,
        ),
        excluded_spans=_scale_spans(
            layout.excluded_spans,
            source_start=source_start,
            source_end=source_end,
            track_start=target_start,
            track_end=target_end,
        ),
        object_centers=tuple(
            _scale_value(
                center,
                source_start=source_start,
                source_end=source_end,
                track_start=target_start,
                track_end=target_end,
            )
            for center in layout.object_centers
        ),
    )


_CUE_LAYOUTS_BY_FAMILY: dict[str, dict[str, tuple[tuple[float, float], ...] | tuple[float, ...]]] = {
    "PO": {
        "rich_spans": ((0.0, 43.0), (81.0, 100.0)),
        "excluded_spans": tuple(),
        "object_centers": (20.0, 36.0, 88.0),
    },
    "POM": {
        "rich_spans": ((0.0, 28.0), (57.0, 100.0)),
        "excluded_spans": tuple(),
        "object_centers": (20.0, 64.0, 88.0),
    },
    "PNO": {
        "rich_spans": tuple(),
        "excluded_spans": tuple(),
        "object_centers": tuple(),
    },
}


def cue_zone_layout_for_condition(
    condition_name: str | None = "",
    *,
    condition_family: str | None = "",
    track_start: float = DEFAULT_TRACK_START,
    track_end: float = DEFAULT_TRACK_END,
) -> CueZoneLayout:
    cname = canonical_condition_name(condition_name)
    family = canonical_condition_name(condition_family) or condition_family_name(cname)
    spec = _CUE_LAYOUTS_BY_FAMILY.get(
        family,
        {
            "rich_spans": tuple(),
            "excluded_spans": tuple(),
            "object_centers": tuple(),
        },
    )
    return CueZoneLayout(
        condition_name=cname,
        condition_family=family,
        track_start=float(track_start),
        track_end=float(track_end),
        rich_spans=tuple((float(a), float(b)) for a, b in spec["rich_spans"]),
        excluded_spans=tuple((float(a), float(b)) for a, b in spec["excluded_spans"]),
        object_centers=tuple(float(x) for x in spec["object_centers"]),
    )


def cue_zone_layout_for_condition_on_track(
    condition_name: str | None = "",
    *,
    condition_family: str | None = "",
    xbin_edges: np.ndarray | None = None,
    track_start: float | None = None,
    track_end: float | None = None,
) -> CueZoneLayout:
    """
    Return a cue-zone layout scaled from canonical 0-100 units to a target track.
    """

    if xbin_edges is not None:
        if track_start is not None or track_end is not None:
            raise ValueError("Pass either xbin_edges or track_start/track_end, not both")
        track_start, track_end = _track_bounds_from_xbin_edges(xbin_edges)
    elif track_start is None and track_end is None:
        track_start, track_end = DEFAULT_TRACK_START, DEFAULT_TRACK_END
    elif track_start is None or track_end is None:
        raise ValueError("track_start and track_end must be provided together")

    canonical = cue_zone_layout_for_condition(
        condition_name,
        condition_family=condition_family,
    )
    return scale_cue_zone_layout(
        canonical,
        track_start=float(track_start),
        track_end=float(track_end),
    )


def _base_condition_for_condway(condway: int) -> int:
    c = int(condway)
    if c < 1:
        raise ValueError(f"condway must be >= 1, got {c}")
    return (c + 1) // 2


def cue_zone_layouts_by_condway(
    condways: np.ndarray | list[int],
    condition_names_by_base: dict[int, str] | None = None,
    *,
    xbin_edges: np.ndarray | None = None,
    track_start: float | None = None,
    track_end: float | None = None,
) -> dict[int, CueZoneLayout]:
    """
    Build scaled cue-zone layouts keyed by condition-direction id.
    """

    condition_names_by_base = condition_names_by_base or {}
    out: dict[int, CueZoneLayout] = {}
    for condway in sorted(set(np.asarray(condways, dtype=np.int64).ravel().tolist())):
        base_condition = _base_condition_for_condway(int(condway))
        condition_name = condition_names_by_base.get(
            int(base_condition),
            f"cond{int(base_condition)}",
        )
        out[int(condway)] = cue_zone_layout_for_condition_on_track(
            condition_name,
            xbin_edges=xbin_edges,
            track_start=track_start,
            track_end=track_end,
        )
    return out


def label_positions_by_zone(
    positions: np.ndarray,
    *,
    layout: CueZoneLayout,
    rich_label: str = "rich",
    poor_label: str = "poor",
) -> np.ndarray:
    values = np.asarray(positions, dtype=np.float64).ravel()

    max_label_len = max(1, len(str(rich_label)), len(str(poor_label)))
    labels = np.full(values.size, "", dtype=f"<U{max_label_len}")
    in_track = (
        np.isfinite(values)
        & (values >= float(layout.track_start))
        & (values <= float(layout.track_end))
    )
    labels[in_track] = poor_label

    for start, end in layout.rich_spans:
        labels[_span_mask(values, start, end, float(layout.track_end))] = rich_label

    for start, end in layout.excluded_spans:
        labels[_span_mask(values, start, end, float(layout.track_end))] = ""
    return labels


def label_xbin_centers_by_zone(
    xbin_centers: np.ndarray,
    *,
    layout: CueZoneLayout | None = None,
    condition_name: str | None = "",
    condition_family: str | None = "",
    rich_label: str = "rich",
    poor_label: str = "poor",
) -> np.ndarray:
    if layout is None:
        layout = cue_zone_layout_for_condition(
            condition_name,
            condition_family=condition_family,
        )
    return label_positions_by_zone(
        xbin_centers,
        layout=layout,
        rich_label=rich_label,
        poor_label=poor_label,
    )


def zone_component_names_for_layout(
    layout: CueZoneLayout,
    *,
    include_rich: bool = True,
    include_poor: bool = True,
) -> tuple[str, ...]:
    names: list[str] = []
    if include_rich:
        names.extend(f"rich_{i + 1}" for i in range(len(layout.rich_spans)))
    if include_poor:
        names.extend(f"poor_{i + 1}" for i in range(len(layout.poor_spans)))
    return tuple(names)


def label_xbin_centers_by_zone_component(
    xbin_centers: np.ndarray,
    *,
    layout: CueZoneLayout | None = None,
    condition_name: str | None = "",
    condition_family: str | None = "",
    rich_prefix: str = "rich",
    poor_prefix: str = "poor",
) -> np.ndarray:
    if layout is None:
        layout = cue_zone_layout_for_condition(
            condition_name,
            condition_family=condition_family,
        )
    return label_positions_by_zone_component(
        xbin_centers,
        layout=layout,
        rich_prefix=rich_prefix,
        poor_prefix=poor_prefix,
    )


def label_positions_by_zone_component(
    positions: np.ndarray,
    *,
    layout: CueZoneLayout,
    rich_prefix: str = "rich",
    poor_prefix: str = "poor",
) -> np.ndarray:
    values = np.asarray(positions, dtype=np.float64).ravel()

    max_index = max(len(layout.rich_spans), len(layout.poor_spans), 1)
    max_label_len = max(
        1,
        len(f"{rich_prefix}_{max_index}"),
        len(f"{poor_prefix}_{max_index}"),
    )
    labels = np.full(values.size, "", dtype=f"<U{max_label_len}")

    for i, (start, end) in enumerate(layout.poor_spans, start=1):
        labels[_span_mask(values, start, end, float(layout.track_end))] = f"{poor_prefix}_{i}"

    for i, (start, end) in enumerate(layout.rich_spans, start=1):
        labels[_span_mask(values, start, end, float(layout.track_end))] = f"{rich_prefix}_{i}"

    for start, end in layout.excluded_spans:
        labels[_span_mask(values, start, end, float(layout.track_end))] = ""
    return labels
