# Cellclass and Models Workflow Audit

Scope: the cell-classification path from raw ratemap MAT files to the
`type_u_comparison` outputs.

## Canonical Path

1. Raw ratemap MAT files under `data/raw` are converted with
   `src/cellclass/mat_to_npz.py --mode ratemap`.
2. The conversion writes `data/interim/<mouse>/<date>/*_allcel.npz` files.
3. `scripts/pipelines/interim_to_processed.py` consumes only `*_allcel.npz`
   files through `cellclass.pipeline.extract_one` and writes:
   - `data/processed/<mouse>/features/*_features.parquet`
   - `data/processed/<mouse>/acg/*_acg_counts_*.npz`
   - `data/processed/<mouse>/qc/*_manifest.json`
4. `scripts/pipelines/aggregate_by_age.py` joins processed features to
   `data/schedule.xlsx` and writes age-group tables under `results/<AGE>/`.
5. `src/models/compare_type_u.py` fits the fixed two-cluster GMM and writes
   `results/type_u_comparison`.

The one-command entry point is:

```powershell
python scripts/pipelines/build_type_u_comparison_from_raw.py
```

Use `--dry-run` first to print the exact commands without writing outputs.

## Findings

- `interim_to_processed.py` previously walked every `*.npz` under
  `data/interim`, including `*_trajdata.npz`, `*_pf.npz`, and `*_pf_c.npz`.
  The cell-classification extractor expects allcel ratemap NPZs, so the
  script now defaults to `*_allcel.npz`.
- The main feature-extraction orchestration has been promoted to
  `src/cellclass/pipeline.py`. `scripts/legacy/one_file_processing.py` is now
  only a compatibility wrapper.
- Feature defaults, age groups, and `type_u` label conventions now live in
  `src/cellclass/config.py`. Aggregation and model scripts import those shared
  defaults instead of carrying independent copies.
- Boundary validators now live in `src/cellclass/validation.py` for interim
  allcel NPZ files, processed feature tables, and age-group model tables.
- `src/cellclass/io.py` and `src/cellclass/classify.py` are still empty or
  placeholder-like, so the package API is better but not yet fully organized.
- Some maintenance/reporting defaults still point at specific historical
  result roots such as `results/type_u_comparison_valero_feats_3`.

## Recommended Next Cleanup

1. Replace script-local `sys.path` hacks with either editable installs or
   `python -m` entry points.
2. Update maintenance/report defaults so canonical outputs are not tied to one
   historical experiment folder.
