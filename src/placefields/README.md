# Placefields

Current Python pipeline for ratemap construction (MATLAB-aligned preprocessing, Python-native data model).

This README reflects the current state of:

- `src/cellclass/mat_to_npz.py` (MAT -> interim NPZ extraction)
- `src/placefields/metrics.py` (core numerical helpers: smoothing/division)
- `src/placefields/matlab_compat.py` (MATLAB/Python index conversions)
- `src/placefields/trials.py` (trial metadata builders)
- `src/placefields/interim_io.py` (interim NPZ pairing/loading/saving helpers)
- `src/placefields/pipeline.py` (ratemap builder)
- `src/placefields/bootstrap.py` (MATLAB-like null-map simulation + empirical p-values)
- `src/placefields/ssi.py` (Spatial Selectivity Index + null-distribution p-values)
- `scripts/build_ratemap_from_interim.py` (session pairing + ratemap run)
- `scripts/build_placefield_null_from_interim.py` (session-level null bootstrap + p-values)

## Workflow

1. Convert MATLAB files to interim NPZ:
   - `--mode ratemap` for `*_Ratemap*.mat` -> `*_allcel.npz`
   - `--mode trajdata` for `*_TrajData*.mat` -> `*_trajdata.npz`
2. Keep paired files in the same folder:
   - `data/interim/<MOUSE>/<DATE>/<SESSION>_allcel.npz`
   - `data/interim/<MOUSE>/<DATE>/<SESSION>_trajdata.npz`
3. Build ratemaps by pairing files with the same session stem:
   - output: `results/ratemap/<MOUSE>/<DATE>/<SESSION>_rmap.npz`

## Data Extraction (MAT -> NPZ)

`src/cellclass/mat_to_npz.py`:

- `ratemap` mode exports `allcel__*` plus selected `allpf__ispf_cxu`.
- `trajdata` mode exports selected `Traj` fields by default:
  - `Cond, time, Wheel, VRtraj, condition, Speed, XSpeed, binSpX, BinSpW, WB, start, stop, tstart, tstop, endVR`
- You can override Traj export fields with `--traj-fields ...` or use `--traj-fields all`.

Example:

```powershell
C:\Users\tadse\miniconda3\envs\odors\python.exe src/cellclass/mat_to_npz.py `
  --mode trajdata `
  --input data/raw/trajdata `
  --output data/interim `
  --recursive
```

## Build Ratemaps From Paired Interim Files

`scripts/build_ratemap_from_interim.py`:

- finds `*_allcel.npz` and `*_trajdata.npz` pairs in `data/interim` (via `interim_io`)
- downsamples spike indices `25000 Hz -> 1000 Hz`
- reconstructs session trajectory from trial vectors
- builds ratemap tensors with `build_ratemap_from_trials(...)`

Example:

```powershell
C:\Users\tadse\miniconda3\envs\odors\python.exe scripts/build_ratemap_from_interim.py `
  --interim_root data/interim `
  --out_root results/ratemap `
  --min_speed 2.0
```

## Build Null Bootstrap Maps/P-values

`scripts/build_placefield_null_from_interim.py`:

- reuses the same interim allcel+traj pairs
- rebuilds ratemap tensors per session
- simulates null maps per selected cell with:
  - `random`
  - `poisson`
  - `circular_shift`
- saves:
  - `dwell_tx_x` (trial-level dwell time, seconds)
  - `dwell_cx_x` (condition-level dwell time, seconds)
  - `pfnull__fr_tx_ux` (observed unsmoothed trial maps for selected cells)
  - `pfnull__fr_s_tx_ux` (observed smoothed trial maps for selected cells)
  - `pfnull__fr_cx_ux` (observed unsmoothed condition maps for selected cells)
  - `pfnull__fr_s_cx_ux` (observed smoothed condition maps for selected cells)
  - `pfnull__pval_tx_ux` (trial-level p-values)
  - `pfnull__pval_cx_ux` (condition-level p-values)
  - `pfnull__ssi_obs_cu` (observed SSI per cell/condition)
  - `pfnull__ssi_pval_cu` (empirical p-value for SSI per cell/condition)
  - `occupancy_cx` (normalized occupancy per condition/bin used for SSI)
  - optional `pfnull__fr_s_txrep_uxtr` (full null maps, large)
  - optional `pfnull__ssi_null_cur` (full SSI null distributions, large)
  - `run_index.csv` (per-session run status)
  - `ssi_classification.csv` (one row per `session_id`/`cell_id`/`condition_1b` with `ssi_obs`, `p_value`, `SM`)
- classification controls:
  - `--sm_alpha` (default `0.05`, `SM=True` when `p_value < sm_alpha`)
  - `--classification_csv` (filename/path for the classification CSV)

Example:

```powershell
C:\Users\tadse\miniconda3\envs\odors\python.exe scripts/build_placefield_null_from_interim.py `
  --interim_root data/interim `
  --out_root results/placefield_null `
  --null_method random `
  --nb_rep 1000 `
  --save_ssi_null `
  --max_cells 25 `
  --min_speed 2.0
```

## Plot Putative Place Cells

`scripts/plot_putative_pcells_from_ssi.py`:

- reads `results/placefield_null/ssi_classification.csv`
- selects putative place cells with `p_value < sm_alpha`
- loads observed maps from each session `*_pfnull.npz` (so `smooth_sigma_bins` and `xbin_rem` match bootstrap)
- saves one figure per condition for each selected cell:
  - mean ratemap (unsmoothed + smoothed)
  - trial-by-trial smoothed heatmap (condition trials only)
  - SSI null histogram for that condition (with observed SSI marker, if available)
  - condition occupancy panel (dwell + normalized occupancy)
  - output path:
    - `results/figures/putative_pcells/smooth_<...>__p_<...>/<AGE>/<SESSION>/cell_<CELL_ID>/cond_<COND>.png`
- writes:
  - `.../plot_index.csv` in the same parameterized folder
- useful args:
  - `--sm_alpha` (default `0.05`, used for selection and plot labels)
  - `--pfnull_root` (default `results/placefield_null`)
  - `--schedule_xlsx` / `--schedule_sheet` (default `data/schedule.xlsx`, `VINCA`) to map session to age
  - `--unknown_age_label` (default `unmapped`) for sessions missing in schedule metadata
  - `--only_sm_conditions` (plot only condition(s) passing `p_value < sm_alpha`)
  - `--hist_bins` (SSI histogram bins)
  - `--heatmap_percentile` (color-scale upper percentile for trial heatmap)
  - `--no_param_subdir` (disable parameter-based output subdirectory)
- output folder naming:
  - smoothing tag is inferred automatically from `pfnull` metadata (`smooth_*` or `smooth_mixed`)
  - p-value tag comes from `--sm_alpha` (`p_*`)
- note:
  - if plotting reports missing keys like `pfnull__fr_s_tx_ux` / `dwell_cx_x`, re-run `build_placefield_null_from_interim.py` with the updated code so `*_pfnull.npz` includes trial maps and dwell arrays.
  - SSI histogram panel requires `pfnull__ssi_null_cur`; if missing, rerun bootstrap with `--save_ssi_null`.

Example:

```powershell
C:\Users\tadse\miniconda3\envs\odors\python.exe scripts/plot_putative_pcells_from_ssi.py `
  --classification_csv results/placefield_null/ssi_classification.csv `
  --pfnull_root results/placefield_null `
  --out_root results/figures/putative_pcells `
  --sm_alpha 0.01
```

## Bootstrapping Algorithm (Step-by-step)

This section describes exactly what the Python null-bootstrap pipeline does for one selected cell.

1. Build observed ratemap tensors with the same preprocessing used everywhere else.
2. For each lap `t`, define `valid_idx` as the intersection of:
   - lap time window,
   - finite position samples,
   - kept spatial range after `xbin_rem`,
   - optional speed mask (`speed >= min_speed`).
3. Keep only spikes that fall on `valid_idx` for that lap.
4. Let:
   - `N = len(valid_idx)` (number of valid samples in that lap),
   - `nbspk = number of kept spikes in that lap`.
5. Generate `nb_rep` null spike realizations for that lap:
   - `random`: sample `nbspk` events uniformly over the `N` valid samples.
   - `poisson`: sample independent Poisson counts per valid sample with mean `nbspk/N`.
   - `circular_shift`: convert kept spike absolute indices to local positions `0..N-1`, shift by `delta` with wrap, map back to absolute indices.
6. Convert each null realization to spatial spike counts (histogram over `xbin_edges`).
7. Smooth null spike counts along space (Gaussian, `smooth_sigma_bins`).
8. Divide by observed smoothed dwell for the same lap:
   - `fr_null = smoothed_counts / dwell_s_tx`.
9. Smooth the resulting null rate map again (same Gaussian).
10. Store `fr_s_txrep[t, x, rep]`.
11. Compute trial-level empirical p-values per bin:
   - `p_tx(t,x) = mean_rep[ fr_null(t,x,rep) >= fr_obs(t,x) ]`.
12. Compute condition-level empirical p-values per bin:
   - average null maps across laps in each condition first,
   - then compare null vs observed condition map bin-by-bin.
13. Compute occupancy per condition from trial dwell:
   - `p_i = occupancy_i / sum_j occupancy_j`.
14. Compute observed SSI per condition:
   - `SSI = sum_i p_i * (lambda_i/lambda_bar) * log2(lambda_i/lambda_bar)`,
   - with `lambda_bar = sum_i p_i * lambda_i`.
15. Compute null SSI distribution per condition:
   - for each repetition, first average null maps across laps of that condition,
   - then evaluate SSI with the same occupancy weights.
16. Compute SSI empirical p-value (one-sided, greater):
   - `k = count(SSI_null >= SSI_real)`,
   - `p = (k + 1) / (n_rep + 1)` (add-one correction).

Important:

- The same binning, masking, and smoothing settings are used for observed and null maps.
- For `circular_shift`, shifting is done on the lap-local valid-sample axis, not directly on absolute time.

## Indexing Convention

- Internal Python indexing is 0-based.
- MATLAB-style trial bounds (`start`, `stop`) are converted to Python slices:
  - start: inclusive
  - stop: exclusive

## Important Functions

`src/placefields/metrics.py`:

- `gaussian_kernel_1d(sigma_bins)`
- `smooth_last_axis(arr, sigma_bins)`
- `divide_with_nan(num, den)`

`src/placefields/interim_io.py`:

- `find_pairs(interim_root, sessions)`:
  pair discovery for `*_allcel.npz` + `*_trajdata.npz`.
- `load_allcel_spikes(allcel_path)`:
  loads spike times/cell ids from interim allcel files.
- `load_traj_fields(traj_path)`:
  loads Cond/WB/start/stop and decodes JSON trial vectors.
  For velocity masking, it uses `XSpeed` by default (falls back to `Speed`).
- `stitch_trial_series(start_1b, stop_1b, series_by_trial)`:
  reconstructs session-level vectors from per-trial arrays.
- `save_ratemap_pack(pack, out_path, meta)`:
  writes compressed ratemap outputs.

`src/placefields/matlab_compat.py`:

- `ifreq_swap(indices, f1_hz, f2_hz)`:
  MATLAB-compatible frequency index conversion
  (`floor((f2/f1)*(if1-1))+1`), used for spike index downsampling.
- `matlab_1b_to_python_0b(...)`

`src/placefields/trials.py`:

- `build_condway(cond, wb)`:
  maps `(Cond, WB)` to condition-direction labels (`W -> 2*i-1`, `B -> 2*i`).
- `build_trial_info_from_traj(...)`:
  builds per-trial metadata (`TrialInfo`).

`src/placefields/pipeline.py`:

- `normalize_x_to_100(position_x)`:
  mirrors MATLAB normalization `X_ds_n = X_ds_n*100/max(X_ds_n)`.
- `build_default_xbin(position_x)`:
  mirrors MATLAB default bin logic.
- `build_ratemap_from_trials(...)`:
  builds all ratemap tensors (`*_tx`, `*_cx`, smoothed/unsmoothed).

`src/placefields/bootstrap.py`:

- `NullBootstrapConfig(method='random'|'poisson'|'circular_shift', nb_rep=...)`
- `simulate_null_fr_s_txrep(...)`:
  one-cell null ratemap generator, Python counterpart of MATLAB `fct_placefield_simtrain`.
- `empirical_pval_tx(...)`:
  trial-level p-values `P(null >= observed)`.
- `empirical_pval_cx(...)`:
  condition-level p-values using condition-wise average of trial null maps.

`src/placefields/ssi.py`:

- `spatial_selectivity_index(rate_x, occupancy_x)`:
  computes
  `sum_i p_i * (lambda_i/lambda_bar) * log2(lambda_i/lambda_bar)`.
- `compute_condition_occupancy(dwell_tx_x, idcond_t, nb_cond)`:
  normalized occupancy per condition from trial dwell.
- `ssi_observed_by_condition(...)` and `ssi_null_by_condition(...)`:
  observed/null SSI per condition.
- `empirical_pvalue_from_null(...)`:
  empirical p-values from null SSI distribution.

## Ratemap Objects

### `RatemapConfig`

Main parameters:

- `freq_hz` (default `1000.0`)
- `smooth_sigma_bins` (default `10.0`)
- `xbin_rem` (default `0`)
- `nb_cond` (default `None`, inferred from trials)
- `min_speed` (default `2.0`)

### `RatemapPack`

Main fields:

- `cell_ids`
- `trial_info` (`list[TrialInfo]`)
- `idcond_t`
- `xbin_edges`, `xbin_centers`
- `nbspk_tx_ux`, `dwell_tx_x`, `fr_tx_ux`
- `nbspk_s_tx_ux`, `dwell_s_tx_x`, `fr_s_tx_ux`
- `nbspk_cx_ux`, `dwell_cx_x`, `fr_cx_ux`
- `nbspk_s_cx_ux`, `dwell_s_cx_x`, `fr_s_cx_ux`

Shape convention:

- trial-level (`*_tx_*`):
  - cell-dependent tensors: `(n_cells, n_trials, n_bins)`
  - dwell tensors: `(n_trials, n_bins)`
- condition-level (`*_cx_*`):
  - cell-dependent tensors: `(n_cells, n_cond, n_bins)`
  - dwell tensors: `(n_cond, n_bins)`

Notes:

- smoothing is applied to `nbspk` and `dwell`, then `fr = smoothed_nbspk / smoothed_dwell`
- low-speed filtering is controlled with `min_speed`
- condition-level `fr_cx`/`fr_s_cx` are unweighted means of trial-level rates (`nanmean` across trials)

## Place-Field Stage

Full place-field thresholding/stability filtering is not yet integrated in the new ratemap runner.
However, MATLAB-style null bootstrapping is now available in Python (`bootstrap.py`) and can be
used as a direct precursor to threshold/stability logic.

Minimal usage sketch (one cell):

```python
from placefields import NullBootstrapConfig, simulate_null_fr_s_txrep, empirical_pval_tx

cfg = NullBootstrapConfig(method="random", nb_rep=1000, seed=0)
fr_s_txrep = simulate_null_fr_s_txrep(
    position_x=x_session,
    spike_indices_cell_0b=spike_idx_cell,
    trials=pack.trial_info,
    xbin_edges=pack.xbin_edges,
    dwell_s_tx_x=pack.dwell_s_tx_x,
    smooth_sigma_bins=10.0,
    cfg=cfg,
    valid_sample_mask=sample_keep_mask,  # optional (e.g. speed-threshold mask)
)
pval_tx = empirical_pval_tx(fr_s_txrep, pack.fr_s_tx_ux[cell_u])
```

Null methods:

- `random`: redistributes observed spike count uniformly across valid samples in each lap.
- `poisson`: independent Poisson counts at each valid sample with mean `nbspk/N`.
- `circular_shift`: circularly shifts observed spike times within each lap to preserve temporal
  structure while breaking spike-position alignment.

How `circular_shift` is applied (important indexing detail):

1. Build `valid_idx` for one lap from the intersection of:
   - lap time window
   - finite/valid position samples
   - kept spatial range after `xbin_rem`
   - optional speed mask
2. Keep only spikes that occur at samples in `valid_idx`.
3. Convert each spike from absolute sample index to local index in `valid_idx` (`0..N-1`).
4. Shift locally with wrap: `local_shift = (local + delta) % N`.
5. Map back to absolute sample indices using `valid_idx[local_shift]`.

So the circular shift is not done directly on absolute time indices, but on the compressed
"valid-sample timeline" of each lap.

Example:

- `valid_idx = [2, 5, 6, 7, 8, 10, 12, 13]`
- observed spikes: `[5, 8, 10, 13]`
- local indices: `[1, 4, 5, 7]`
- with `delta=2` and wrap over `N=8`: `[3, 6, 7, 1]`
- shifted absolute spikes: `[7, 12, 13, 5]` (same set as sorted `[5, 7, 12, 13]`)
