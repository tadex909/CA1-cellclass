# Cellclass and Models Workflow Audit

Scope: the cell-classification path from raw ratemap MAT files to the
`type_u_comparison_valero_feats_3` outputs.

## Canonical Path

1. Raw ratemap MAT files under `data/raw` are converted with
   `python -m cellclass.mat_to_npz --mode ratemap`.
2. The conversion writes `data/interim/<mouse>/<date>/*_allcel.npz` files.
3. `scripts/pipelines/interim_to_processed.py` consumes only `*_allcel.npz`
   files through `cellclass.pipeline.extract_one` and writes:
   - `data/processed/<mouse>/features/*_features.parquet`
   - `data/processed/<mouse>/acg/*_acg_counts_*.npz`
   - `data/processed/<mouse>/qc/*_manifest.json`
4. `scripts/pipelines/aggregate_by_age.py` joins processed features to
   `data/schedule.xlsx` and writes age-group tables under `results/<AGE>/`.
5. `python -m models.compare_type_u` fits the fixed two-cluster GMM and writes
   `results/type_u_comparison_valero_feats_3`, including the compact
   `cell_classification_table.csv` consumed by downstream place-field/SSI
   scripts.
   The default comparison feature set is:
   `cv2, acg_peak_latency_ms, spk_duration_ms, spk_asymmetry, log_fr_hz_session`.

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
  `src/cellclass/pipeline.py`. The old
  `scripts/legacy/one_file_processing.py` compatibility wrapper has been
  removed.
- Feature defaults, age groups, and `type_u` label conventions now live in
  `src/cellclass/config.py`. Aggregation and model scripts import those shared
  defaults instead of carrying independent copies.
- Boundary validators now live in `src/cellclass/validation.py` for interim
  allcel NPZ files, processed feature tables, and age-group model tables.
- Script-local `sys.path` edits have been removed from package-consuming
  scripts. The supported setup is `python -m pip install -e .`, with model
  commands run as `python -m models.<module>`.
- Maintenance/report defaults now point at the current canonical
  `results/type_u_comparison_valero_feats_3` outputs.

## Recommended Next Cleanup

No remaining cleanup items from this workflow audit.
