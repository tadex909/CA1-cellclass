# Scripts Overview

This folder is organized by intent. See `scripts/INVENTORY.md` for the full
status table and deletion candidates.

Run scripts from the repository root after installing the project once:

```powershell
python -m pip install -e .
```

- `scripts/pipelines/`
  - Data-production pipelines (interim -> processed, ratemap/null generation, aggregation).
- `scripts/reports/`
  - Plotting and cross-method comparison outputs.
- `scripts/maintenance/`
  - Repository/results maintenance utilities (registry, reorg view, NPZ checks).

## Status Rules

- Use files under `scripts/pipelines/`, `scripts/reports/`, and
  `scripts/maintenance/` for new commands.
- Do not add new top-level wrappers; add new commands to the appropriate
  subdirectory.

## Canonical Entry Points

Use these paths for new commands:

- `scripts/pipelines/interim_to_processed.py`
- `scripts/pipelines/build_type_u_comparison_from_raw.py`
- `scripts/pipelines/aggregate_by_age.py`
- `scripts/pipelines/aggregate_mouse.py`
- `scripts/pipelines/build_ratemap_from_interim.py`
- `scripts/pipelines/build_placefield_null_from_interim.py`
- `scripts/pipelines/build_bayesian_decoder_from_interim.py`
- `scripts/pipelines/build_population_geometry_from_ratemap.py`
- `scripts/reports/plot_putative_pcells_from_ssi.py`
- `scripts/reports/compare_placefield_methods.py`
- `scripts/reports/compare_pf_npz_methods.py`
- `scripts/reports/plot_pf_method_disagreements.py`
- `scripts/reports/plot_pf_both_methods_grouped_conditions.py`
- `scripts/maintenance/build_cell_classification_table.py`
- `scripts/maintenance/build_ssi_session_cell_counts.py`
- `scripts/maintenance/build_results_registry.py`
- `scripts/maintenance/create_results_view.py`
- `scripts/maintenance/check_npz.py`

## Removed Compatibility Wrappers

Top-level compatibility wrappers such as `scripts/interim_to_processed.py` and
`scripts/build_ratemap_from_interim.py` have been removed. Use the canonical
paths listed above.

The old `scripts/legacy/one_file_processing.py` shim has also been removed.
Use `python -m cellclass.pipeline` for one-file extraction or
`scripts/pipelines/interim_to_processed.py` for batch extraction.
