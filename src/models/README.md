# Models Pipeline

This folder contains model-selection code for unsupervised cell-type clustering.

## Main Script

- `fitting.py`: systematic Gaussian Mixture Model (GMM) search across:
  - different numbers of clusters (`k`)
  - different feature subsets (including `n` vs `n-1` features)
- `compare_type_u.py`: fixed 2-cluster GMM on all features, then comparison against
  external labels `allcel__type_u` (`0=interneuron`, `1=pyramidal`).

It is designed to run on age-group datasets already produced in `results/<AGE>/<AGE>_clean_units.parquet`.

## What `fitting.py` Does

For each age group:

1. Loads clean units.
2. Builds candidate feature subsets (`full_only`, `leave_one_out`, or `all`).
3. Fits GMMs for each subset and each `k`.
4. Computes clustering quality metrics:
   - `BIC` (lower is better)
   - `AIC` (lower is better)
   - `silhouette` (higher is better)
   - `Calinski-Harabasz` (higher is better)
   - `Davies-Bouldin` (lower is better)
5. Rejects models with very small clusters using `--min_cluster_size`.
6. Selects the best valid model using ordered criteria:
   - lowest `BIC` (primary)
   - highest `silhouette` (tie-break 1)
   - lowest `Davies-Bouldin` (tie-break 2)
   - lowest `AIC` (tie-break 3)
   - highest `Calinski-Harabasz` (final tie-break)
7. Saves:
   - full evaluation table
   - ranked table
   - labels from the best model
   - summary JSON

## Basic Usage

Run from repository root:

```powershell
C:\Users\tadse\miniconda3\envs\odors\python.exe src/models/fitting.py `
  --results_root results `
  --out_root results/model_selection `
  --subset_mode leave_one_out `
  --k_min 2 --k_max 6 `
  --n_init 5 `
  --min_cluster_size 3
```

## Most Useful Arguments

- `--subset_mode`
  - `full_only`: only all features
  - `leave_one_out`: all features + every `n-1` subset
  - `all`: all subsets down to `--min_subset_size`
- `--k_min`, `--k_max`: cluster range
- `--min_cluster_size`: filters over-fragmented solutions
- `--min_units`: minimum rows required to evaluate a subset
- `--features`: comma-separated feature pool
- `--age_groups`: comma-separated age groups (default: all standard groups)
- `--no_log_fr`: disable `log10(fr_hz)` transform
- `--no_standardize`: disable z-scoring

## Outputs

Outputs are written under `results/model_selection` (or your chosen `--out_root`):

- `run_config.json`
- `/<AGE>/<AGE>_gmm_model_selection.csv`
- `/<AGE>/<AGE>_gmm_model_selection_ranked.csv`
- `/<AGE>/<AGE>_best_gmm_labels.parquet`
- `/<AGE>/<AGE>_gmm_summary.json`

## Notes

- `k=2` best models also get a `putative_type` mapping (`interneuron`/`pyramidal`) based on shorter `spk_duration_ms`.
- For small age groups, increase caution even with `min_cluster_size` filtering.

## Recommended Defaults (Current Dataset)

Given your current sample sizes (roughly `P25` smallest, `P17_18`/`P23_24` larger), these are sensible defaults:

- `--subset_mode leave_one_out`
  - good balance between speed and feature-ablation insight (`n` vs `n-1`)
- `--k_min 2 --k_max 6`
  - enough range to test simple vs richer structure
- `--min_cluster_size 3` (or `4` for stricter runs)
  - prevents tiny, unstable clusters
- `--min_units 30`
  - skips underpowered subset evaluations
- keep transforms enabled (default)
  - use `log10(fr_hz)` and z-scoring for stable optimization

Suggested baseline run:

```powershell
C:\Users\tadse\miniconda3\envs\odors\python.exe src/models/fitting.py `
  --results_root results `
  --out_root results/model_selection `
  --subset_mode leave_one_out `
  --k_min 2 --k_max 6 `
  --n_init 5 `
  --min_cluster_size 3 `
  --min_units 30
```

For a stricter stability pass:

```powershell
C:\Users\tadse\miniconda3\envs\odors\python.exe src/models/fitting.py `
  --results_root results `
  --out_root results/model_selection_strict `
  --subset_mode leave_one_out `
  --k_min 2 --k_max 5 `
  --n_init 10 `
  --min_cluster_size 4 `
  --min_units 35
```

## Compare With `type_u` Labels

Use `compare_type_u.py` when you want a direct benchmark of GMM-vs-`type_u` agreement.

What it does:

1. Uses all default features from `fitting.py`.
2. Fits a 2-cluster GMM (`k=2`) for each age group.
3. Maps clusters to classes using `spk_duration_ms`:
   - shorter duration cluster => `interneuron`
   - other cluster => `pyramidal`
4. Compares predicted class to `allcel__type_u` (`0=interneuron`, `1=pyramidal`).
5. Writes discrepancy summaries by neuron, session, and age group.

Example:

```powershell
C:\Users\tadse\miniconda3\envs\odors\python.exe src/models/compare_type_u.py `
  --results_root results `
  --out_root results/type_u_comparison `
  --n_init 10
```

Main outputs:

- `results/type_u_comparison/summary.json`
- `results/type_u_comparison/discrepancies_by_age.csv`
- `results/type_u_comparison/discrepancies_by_session.csv`
- `results/type_u_comparison/discrepant_neurons.parquet`
- `results/type_u_comparison/discrepant_neurons_top500.csv`
- per-age tables:
  - `results/type_u_comparison/<AGE>/<AGE>_gmm2_vs_type_u.parquet`
