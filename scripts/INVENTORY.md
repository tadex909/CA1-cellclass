# Scripts Inventory

This inventory marks which scripts are canonical entry points and which scripts
only exist for compatibility. Use this file when deciding what can be deleted or
ignored.

## Status Labels

- `canonical-pipeline`: important production workflow script; use for new runs.
- `canonical-report`: downstream report/plot script; optional, but not duplicated.
- `maintenance`: repository/results utility; useful support script, not a data pipeline.
- `helper`: internal support code, not a user-facing command.

## Canonical Pipelines

| Script | Status | Keep? | Unique responsibility |
| --- | --- | --- | --- |
| `scripts/pipelines/build_type_u_comparison_from_raw.py` | `canonical-pipeline` | yes | One-command cellclass workflow from raw MAT files to `results/type_u_comparison_valero_feats_3`. |
| `scripts/pipelines/interim_to_processed.py` | `canonical-pipeline` | yes | Batch conversion from interim `*_allcel.npz` files to processed features, ACG arrays, and QC manifests. |
| `scripts/pipelines/aggregate_by_age.py` | `canonical-pipeline` | yes | Joins processed features with schedule metadata and writes age-group result tables. |
| `scripts/pipelines/aggregate_mouse.py` | `canonical-pipeline` | yes, secondary | Mouse-level aggregation for per-mouse inspection; not the same output as age-group aggregation. |
| `scripts/pipelines/build_ratemap_from_interim.py` | `canonical-pipeline` | yes | Builds saved ratemap packs from interim allcel+traj pairs. |
| `scripts/pipelines/build_placefield_null_from_interim.py` | `canonical-pipeline` | yes | Builds null place-field/SSI outputs and `ssi_classification.csv`. |
| `scripts/pipelines/build_bayesian_decoder_from_interim.py` | `canonical-pipeline` | yes | Builds Bayesian position-decoding outputs from interim allcel+traj pairs. |
| `scripts/pipelines/build_population_geometry_from_ratemap.py` | `canonical-pipeline` | yes | Builds population-geometry summaries from saved ratemap packs. |

## Canonical Reports

| Script | Status | Keep? | Unique responsibility |
| --- | --- | --- | --- |
| `scripts/reports/compare_placefield_methods.py` | `canonical-report` | yes | Compares two SSI/place-field classification CSV outputs. |
| `scripts/reports/compare_pf_npz_methods.py` | `canonical-report` | yes | Compares legacy MATLAB-style `*_pf.npz` and `*_pf_c.npz` place-field outputs. |
| `scripts/reports/plot_pf_both_methods_grouped_conditions.py` | `canonical-report` | yes | Plots cells with place fields in both compared methods, grouped by condition family. |
| `scripts/reports/plot_pf_method_disagreements.py` | `canonical-report` | yes | Plots ratemaps for rows where two place-field methods disagree. |
| `scripts/reports/plot_putative_pcells_from_ssi.py` | `canonical-report` | yes | Plots putative place cells from processed SSI classifications. |

## Maintenance Utilities

| Script | Status | Keep? | Unique responsibility |
| --- | --- | --- | --- |
| `scripts/maintenance/build_cell_classification_table.py` | `maintenance` | yes | Rebuilds compact `cell_classification_table.csv` from canonical comparison outputs or legacy classification files. |
| `scripts/maintenance/build_results_registry.py` | `maintenance` | yes | Builds inventory tables for result folders and run metadata. |
| `scripts/maintenance/build_ssi_session_cell_counts.py` | `maintenance` | yes | Summarizes SSI outputs by session using the canonical cell classification table. |
| `scripts/maintenance/check_npz.py` | `maintenance` | yes, diagnostic | Inspects NPZ keys/shapes for debugging data conversions. |
| `scripts/maintenance/create_results_view.py` | `maintenance` | yes | Creates an organized non-destructive view of `results/`. |

## Removed Compatibility Wrappers

These top-level files duplicated canonical scripts and have been removed. Use
the canonical targets instead.

| Removed wrapper | Use instead |
| --- | --- |
| `scripts/aggreggate_by_age.py` | `scripts/pipelines/aggregate_by_age.py` |
| `scripts/agreggate_mouse.py` | `scripts/pipelines/aggregate_mouse.py` |
| `scripts/build_placefield_null_from_interim.py` | `scripts/pipelines/build_placefield_null_from_interim.py` |
| `scripts/build_population_geometry_from_ratemap.py` | `scripts/pipelines/build_population_geometry_from_ratemap.py` |
| `scripts/build_ratemap_from_interim.py` | `scripts/pipelines/build_ratemap_from_interim.py` |
| `scripts/build_results_registry.py` | `scripts/maintenance/build_results_registry.py` |
| `scripts/check_npz.py` | `scripts/maintenance/check_npz.py` |
| `scripts/compare_pf_npz_methods.py` | `scripts/reports/compare_pf_npz_methods.py` |
| `scripts/compare_placefield_methods.py` | `scripts/reports/compare_placefield_methods.py` |
| `scripts/create_results_view.py` | `scripts/maintenance/create_results_view.py` |
| `scripts/interim_to_processed.py` | `scripts/pipelines/interim_to_processed.py` |
| `scripts/one_file_processing.py` | `python -m cellclass.pipeline` or `scripts/pipelines/interim_to_processed.py` |
| `scripts/legacy/one_file_processing.py` | `python -m cellclass.pipeline` or `scripts/pipelines/interim_to_processed.py` |
| `scripts/plot_pf_both_methods_grouped_conditions.py` | `scripts/reports/plot_pf_both_methods_grouped_conditions.py` |
| `scripts/plot_pf_method_disagreements.py` | `scripts/reports/plot_pf_method_disagreements.py` |
| `scripts/plot_putative_pcells_from_ssi.py` | `scripts/reports/plot_putative_pcells_from_ssi.py` |

## Deletion Candidates

No top-level compatibility wrappers remain. Do not delete the canonical
pipeline/report/maintenance scripts without first checking whether their outputs
are reproduced elsewhere.
