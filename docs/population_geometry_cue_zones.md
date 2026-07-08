# Population Geometry Cue-Zone Definitions

This note records the canonical cue-rich and cue-poor spatial layouts used for
population-geometry summaries on the normalized `0-100` linear track.

## Conventions

- Track coordinates are in normalized physical track units from `0` to `100`.
- Cue-rich spans are recorded as continuous spatial intervals.
- Excluded spans, when present, are labeled as neither rich nor poor.
- Cue-poor spans are defined as the complement of cue-rich spans plus excluded spans on the track.
- For lag-profile analyses, bins are usually labeled by the zone containing the
  bin center.
- Component-wise labels are also available when the contiguous patches should be
  kept separate, for example `rich_1`, `rich_2`, `poor_1`, `poor_2`, and so on.

There is currently no canonical excluded span for the PO/POM/PNO layouts.

## PO Family

The `PO` family includes:

- `PO`
- `PO2`
- `PO3`
- `POnM`

Canonical cue-rich spans:

- `0` to `43`
- `81` to `100`

Canonical object centers:

- `20`
- `36`
- `88`

The resulting cue-poor spans are:

- `43` to `81`

Component labels for `PO` therefore map to:

- `rich_1`: `0` to `43`
- `rich_2`: `81` to `100`
- `poor_1`: `43` to `81`

## POM Family

The `POM` family includes:

- `POM`
- `POMb`

Relative to `PO`, the second object is moved to `100 - 36 = 64`.

Canonical cue-rich spans:

- `0` to `28`
- `57` to `100`

Canonical object centers:

- `20`
- `64`
- `88`

The resulting cue-poor spans are:

- `28` to `57`

Component labels for `POM` therefore map to:

- `rich_1`: `0` to `28`
- `rich_2`: `57` to `100`
- `poor_1`: `28` to `57`

## PNO

`PNO` has no objects and therefore no cue-rich spans.

- cue-rich spans: none
- excluded span: none
- cue-poor span: `0` to `100`

## Recommended Use In Comparisons

When comparing a cue session with a `PNO` session:

- use the cue-session layout to define homologous `rich` and `poor` windows
- apply those same windows to the `PNO` geometry for a fair spatial comparison

This allows summaries such as:

- `g_rich(Δx)`
- `g_poor(Δx)`
- `g_no(Δx)` on the homologous cue-defined windows
- `g_rich_1_within(Δx)` and `g_rich_2_within(Δx)` when the two cue-rich patches
  should be compared separately

## Code Source

The canonical programmatic definitions live in:

- `src/placefields/cue_zones.py`
