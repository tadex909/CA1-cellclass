# CA1-cellclass

End-to-end workflow for CA1 neuron feature extraction, unsupervised cell-type
classification, and downstream place-field/SSI analyses.

The current canonical cell-classification output is:

```text
results/type_u_comparison_valero_feats_3/
```

Downstream place-field and SSI scripts should use:

```text
results/type_u_comparison_valero_feats_3/cell_classification_table.csv
```

Older `results/type_u_comparison*` folders are historical experiment variants.

---

## Repository Layout

- `src/cellclass/`
  - Raw/interim conversion, feature extraction, shared configuration, and
    workflow validation.
- `src/models/`
  - GMM model selection, fixed `k=2` comparison to `allcel__type_u`,
    disagreement review, and stability analysis.
- `src/placefields/`
  - Ratemaps, place-field/null analyses, SSI helpers, Bayesian decoding, and
    population geometry utilities.
- `scripts/`
  - CLI entry points organized by intent:
    - `scripts/pipelines/`: data-producing workflows.
    - `scripts/reports/`: plots and comparison reports.
    - `scripts/maintenance/`: registries, checks, and non-destructive result views.
- `data/`
  - Local raw/interim/processed data artifacts.
- `results/`
  - Local generated analysis outputs. See `results/README.md`.

Compatibility wrappers at the top of `scripts/` were removed. Use the canonical
paths under `scripts/pipelines/`, `scripts/reports/`, and `scripts/maintenance/`.

---

## Setup

Run commands from the repository root.

```powershell
python -m pip install -e .
```

For development tools:

```powershell
python -m pip install -e .[dev]
```

If you have not installed the package, run commands with `PYTHONPATH=src`.

Run the test suite:

```powershell
python -m pytest tests
```

---

## Main Cell-Class Workflow

The one-command workflow from raw MAT files to the canonical comparison output is:

```powershell
python scripts/pipelines/build_type_u_comparison_from_raw.py --dry-run
python scripts/pipelines/build_type_u_comparison_from_raw.py
```

The default comparison uses the current Valero feature set:

```text
cv2, acg_peak_latency_ms, spk_duration_ms, spk_asymmetry, log_fr_hz_session
```

It writes to:

```text
results/type_u_comparison_valero_feats_3/
```

The same workflow can be run manually:

```powershell
python -m cellclass.mat_to_npz --mode ratemap --input data/raw --output data/interim --recursive

python scripts/pipelines/interim_to_processed.py `
  --interim_root data/interim `
  --processed_root data/processed `
  --pattern "*_allcel.npz" `
  --skip_existing

python scripts/pipelines/aggregate_by_age.py `
  --processed_root data/processed `
  --excel data/schedule.xlsx `
  --outdir results `
  --qc

python -m models.compare_type_u `
  --results_root results `
  --out_root results/type_u_comparison_valero_feats_3
```

---

## Optional Model Analyses

Model selection:

```powershell
python -m models.fitting `
  --results_root results `
  --out_root results/model_selection `
  --subset_mode leave_one_out `
  --k_min 2 `
  --k_max 6 `
  --n_init 5 `
  --min_cluster_size 3 `
  --min_units 30
```

Stability analysis:

```powershell
python -m models.stability_analysis `
  --results_root results `
  --out_root results/stability `
  --age_group P19-21 `
  --subset_mode leave_one_out `
  --seeds 0,1,2,3,4,5,6,7,8,9 `
  --n_init 10
```

Disagreement review/export:

```powershell
python -m models.export_disagreement_excel
python -m models.disagreement_review --top_n 500
```

---

## Place-Field And SSI Workflows

Build ratemaps from interim files:

```powershell
python scripts/pipelines/build_ratemap_from_interim.py `
  --interim_root data/interim `
  --out_root results/ratemap
```

Build null place-field/SSI outputs:

```powershell
python scripts/pipelines/build_placefield_null_from_interim.py `
  --interim_root data/interim `
  --out_root results/placefield_null
```

The null/SSI pipeline defaults to the canonical classification table:

```text
results/type_u_comparison_valero_feats_3/cell_classification_table.csv
```

In that table, `pred_type` is the current GMM classification and `u_type` is
the legacy `allcel__type_u` label normalized to `pyramidal`/`interneuron`.

Other available pipelines include:

- `scripts/pipelines/build_bayesian_decoder_from_interim.py`
- `scripts/pipelines/build_population_geometry_from_ratemap.py`

---

## Results Maintenance

Refresh the results registry:

```powershell
python scripts/maintenance/build_results_registry.py `
  --results_root results `
  --out_root results/tables/results_registry
```

Create or refresh a non-destructive organized view:

```powershell
python scripts/maintenance/create_results_view.py `
  --results_root results `
  --view_root results/_organized
```

Materialize that view only when needed:

```powershell
python scripts/maintenance/create_results_view.py `
  --results_root results `
  --view_root results/_organized `
  --mode hardlink `
  --execute
```

`results/_organized/` is derived and can be deleted/recreated. It is not a
canonical output.

---

## Conventions

- Shared feature lists, age groups, label conventions, and canonical result
  paths live in `src/cellclass/config.py`.
- Boundary validators live in `src/cellclass/validation.py`.
- Run model modules with `python -m models.<module>` after editable install.
- Do not add new top-level script wrappers.
- Generated outputs under `data/`, `results/`, and `path/` are local artifacts
  by default.
