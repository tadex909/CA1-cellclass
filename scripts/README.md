# Scripts Overview

This folder is organized by intent:

- `scripts/pipelines/`
  - Data-production pipelines (interim -> processed, ratemap/null generation, aggregation).
- `scripts/reports/`
  - Plotting and cross-method comparison outputs.
- `scripts/maintenance/`
  - Repository/results maintenance utilities (registry, reorg view, NPZ checks).
- `scripts/legacy/`
  - Older scripts kept for compatibility/reference.

## Canonical Entry Points

Use these paths for new commands:

- `scripts/pipelines/interim_to_processed.py`
- `scripts/pipelines/aggregate_by_age.py`
- `scripts/pipelines/aggregate_mouse.py`
- `scripts/pipelines/build_ratemap_from_interim.py`
- `scripts/pipelines/build_placefield_null_from_interim.py`
- `scripts/reports/plot_putative_pcells_from_ssi.py`
- `scripts/reports/compare_placefield_methods.py`
- `scripts/reports/compare_pf_npz_methods.py`
- `scripts/reports/plot_pf_method_disagreements.py`
- `scripts/maintenance/build_results_registry.py`
- `scripts/maintenance/create_results_view.py`
- `scripts/maintenance/check_npz.py`

## Compatibility Wrappers

Legacy script paths still exist at `scripts/<old_name>.py` and forward execution to new paths.
They print a deprecation message and then execute the moved script.

These wrappers are temporary and can be removed once your workflows are updated.
