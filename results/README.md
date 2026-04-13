# Results Layout

This repository keeps `results/` as a working area for generated artifacts.

## Recommended Structure

- `results/datasets/...`
  - Aggregated datasets (age groups, per-mouse exports).
- `results/models/...`
  - Model selection, stability, and type-comparison outputs.
- `results/placefields/...`
  - Ratemaps, null bootstraps, and placefield-specific figures/tables.
- `results/reports/...`
  - Cross-method summaries and final reporting tables.

The current historical layout may still contain legacy top-level folders
(`ratemap`, `placefield_null`, `type_u_comparison_*`, etc.).

## Run-Scoped Outputs

Placefield build scripts now support run-scoped output directories:

- `scripts/pipelines/build_ratemap_from_interim.py`
- `scripts/pipelines/build_placefield_null_from_interim.py`

By default they create:

- `<out_root>/<run_id>/...`

and write:

- `run_config.json` in the run folder
- `LATEST_RUN.txt` in the base output folder

Use `--no_run_subdir` to keep legacy behavior.

## Registry And Discoverability

Build an inventory/registry of runs and outputs:

```powershell
python scripts/maintenance/build_results_registry.py --results_root results --out_root results/tables/results_registry
```

This writes:

- `results/tables/results_registry/run_registry.csv`
- `results/tables/results_registry/top_level_inventory.csv`
- `results/tables/results_registry/latest_pointers.csv`
- JSON mirrors plus `summary.json`

## Non-Destructive Reorganization View

Create an organized mirror of `results/` without deleting original files:

```powershell
# dry-run
python scripts/maintenance/create_results_view.py --results_root results --view_root results/_organized

# materialize using hardlinks when possible
python scripts/maintenance/create_results_view.py --results_root results --view_root results/_organized --mode hardlink --execute
```

The generated view keeps legacy paths untouched and provides a cleaner navigation tree.
