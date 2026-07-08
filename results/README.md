# Results Directory

`results/` is for generated analysis outputs and is ignored by git except for
this README. Do not rely on manual moves inside this folder as part of the
workflow; instead, use the canonical output roots and the maintenance utilities
below.

## Canonical Top-Level Layout

| Root | Purpose | Produced by |
| --- | --- | --- |
| `results/<AGE>/` | Age-group datasets for modeling, such as clean/all unit tables. | `scripts/pipelines/aggregate_by_age.py` |
| `results/<MOUSE>/` | Mouse-level aggregate outputs. | `scripts/pipelines/aggregate_mouse.py` |
| `results/model_selection*/` | GMM model-selection outputs. | `python -m models.fitting` |
| `results/type_u_comparison_valero_feats_3/` | Current GMM-vs-`type_u` comparison outputs and `cell_classification_table.csv`. | `python -m models.compare_type_u` |
| `results/type_u_comparison*/` | Historical comparison variants; keep for reference unless deliberately archiving/deleting old experiments. | older `python -m models.compare_type_u` runs |
| `results/stability*/` | Per-age stability outputs. | `python -m models.stability_analysis` |
| `results/feature_set_experiments*/` | Feature-set experiment outputs. | `python -m models.feature_set_experiments` |
| `results/ratemap/` | Saved ratemap packs. | `scripts/pipelines/build_ratemap_from_interim.py` |
| `results/placefield_null/` | Null place-field/SSI outputs. | `scripts/pipelines/build_placefield_null_from_interim.py` |
| `results/position_decoding/` | Bayesian position-decoding outputs. | `scripts/pipelines/build_bayesian_decoder_from_interim.py` |
| `results/population_geometry/` | Population-geometry outputs. | `scripts/pipelines/build_population_geometry_from_ratemap.py` |
| `results/tables/` | Registry and comparison tables. | `scripts/maintenance/*` and `scripts/reports/*` |
| `results/figures/` | Generated plots. | `scripts/reports/*` |

Only `results/README.md` should be maintained as a loose top-level file. Any
other top-level files are included in the registry and organized view as cleanup
candidates.

The current classification source for downstream place-field/SSI workflows is
`results/type_u_comparison_valero_feats_3/cell_classification_table.csv`.
Its `pred_type` column is the current GMM label; `u_type` is the legacy
`allcel__type_u` label normalized to `pyramidal`/`interneuron`.

## Non-Destructive Organization

Use the registry to index what is present:

```powershell
python scripts/maintenance/build_results_registry.py --results_root results --out_root results/tables/results_registry
```

Use the organized view when you want a cleaner browsing layout without moving
the original outputs:

```powershell
python scripts/maintenance/create_results_view.py --results_root results --view_root results/_organized
python scripts/maintenance/create_results_view.py --results_root results --view_root results/_organized --mode hardlink --execute
```

The first command writes only a dry-run plan. The second materializes the view
using hardlinks where possible.
