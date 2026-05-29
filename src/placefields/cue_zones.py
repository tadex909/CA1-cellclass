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


_CUE_LAYOUTS_BY_FAMILY: dict[str, dict[str, tuple[tuple[float, float], ...] | tuple[float, ...]]] = {
    "PO": {
        "rich_spans": ((13.0, 43.0), (81.0, 96.0)),
        "excluded_spans": ((0.0, 10.0),),
        "object_centers": (20.0, 36.0, 88.0),
    },
    "POM": {
        "rich_spans": ((13.0, 28.0), (57.0, 96.0)),
        "excluded_spans": ((0.0, 10.0),),
        "object_centers": (20.0, 64.0, 88.0),
    },
    "PNO": {
        "rich_spans": tuple(),
        "excluded_spans": ((0.0, 10.0),),
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
            "excluded_spans": ((0.0, 10.0),),
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


def label_xbin_centers_by_zone(
    xbin_centers: np.ndarray,
    *,
    layout: CueZoneLayout | None = None,
    condition_name: str | None = "",
    condition_family: str | None = "",
    rich_label: str = "rich",
    poor_label: str = "poor",
) -> np.ndarray:
    centers = np.asarray(xbin_centers, dtype=np.float64).ravel()
    if layout is None:
        layout = cue_zone_layout_for_condition(
            condition_name,
            condition_family=condition_family,
        )

    max_label_len = max(1, len(str(rich_label)), len(str(poor_label)))
    labels = np.full(centers.size, "", dtype=f"<U{max_label_len}")
    in_track = (
        np.isfinite(centers)
        & (centers >= float(layout.track_start))
        & (centers <= float(layout.track_end))
    )
    labels[in_track] = poor_label

    for start, end in layout.excluded_spans:
        in_excluded = np.isfinite(centers) & (centers >= float(start)) & (centers < float(end))
        labels[in_excluded] = ""

    for start, end in layout.rich_spans:
        in_span = np.isfinite(centers) & (centers >= float(start)) & (centers < float(end))
        labels[in_span] = rich_label
    return labels


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
    centers = np.asarray(xbin_centers, dtype=np.float64).ravel()
    if layout is None:
        layout = cue_zone_layout_for_condition(
            condition_name,
            condition_family=condition_family,
        )

    max_index = max(len(layout.rich_spans), len(layout.poor_spans), 1)
    max_label_len = max(
        1,
        len(f"{rich_prefix}_{max_index}"),
        len(f"{poor_prefix}_{max_index}"),
    )
    labels = np.full(centers.size, "", dtype=f"<U{max_label_len}")

    for i, (start, end) in enumerate(layout.poor_spans, start=1):
        in_span = np.isfinite(centers) & (centers >= float(start)) & (centers < float(end))
        labels[in_span] = f"{poor_prefix}_{i}"

    for i, (start, end) in enumerate(layout.rich_spans, start=1):
        in_span = np.isfinite(centers) & (centers >= float(start)) & (centers < float(end))
        labels[in_span] = f"{rich_prefix}_{i}"
    return labels
