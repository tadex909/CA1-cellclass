# Population Geometry Cue-Zone Definitions

This note records the canonical cue-rich and cue-poor spatial layouts used for
population-geometry summaries on the normalized `0-100` linear track.

## Conventions

- Track coordinates are in normalized physical track units from `0` to `100`.
- The first `10` track units are excluded from zone-based rich/poor analyses.
- Cue-rich spans are recorded as continuous spatial intervals.
- Excluded spans are labeled as neither rich nor poor.
- Cue-poor spans are defined as the complement of cue-rich spans plus excluded spans on the track.
- For lag-profile analyses, bins are usually labeled by the zone containing the
  bin center.

Canonical excluded span for all condition families:

- `0` to `10`

## PO Family

The `PO` family includes:

- `PO`
- `PO2`
- `PO3`
- `POnM`

Canonical cue-rich spans:

- `13` to `43`
- `81` to `96`

Canonical object centers:

- `20`
- `36`
- `88`

The resulting cue-poor spans are:

- `10` to `13`
- `43` to `81`
- `96` to `100`

## POM Family

The `POM` family includes:

- `POM`
- `POMb`

Relative to `PO`, the second object is moved to `100 - 36 = 64`.

Canonical cue-rich spans:

- `13` to `28`
- `57` to `96`

Canonical object centers:

- `20`
- `64`
- `88`

The resulting cue-poor spans are:

- `10` to `13`
- `28` to `57`
- `96` to `100`

## PNO

`PNO` has no objects and therefore no cue-rich spans.

- cue-rich spans: none
- excluded span: `0` to `10`
- cue-poor span: `10` to `100`

## Recommended Use In Comparisons

When comparing a cue session with a `PNO` session:

- use the cue-session layout to define homologous `rich` and `poor` windows
- apply those same windows to the `PNO` geometry for a fair spatial comparison

This allows summaries such as:

- `g_rich(Δx)`
- `g_poor(Δx)`
- `g_no(Δx)` on the homologous cue-defined windows

## Code Source

The canonical programmatic definitions live in:

- `src/placefields/cue_zones.py`
