# Rmap_G.m (Rate Maps + Place Fields)

This document explains the MATLAB pipeline in [Rmap_G.m](/c:/Users/tadse/OneDrive/Documenti/GitHub/CA1-cellclass/src/placefields/matlab/Rmap_G.m), focusing only on:

- `rmap` computation (rate maps)
- `pf` computation (place-field detection)

Everything after that (cell typing, active/silent labels, percentages, etc.) is intentionally out of scope here.

## 1) What this script is for

`Rmap_G.m` computes per-neuron spatial firing maps and place fields on a 1D trajectory, split by experimental condition and running direction.

Core output at this stage:

- `rmap(g)`: firing-rate maps for neuron `g` (trial-level and condition-level)
- `pf(g)`: place-field masks/statistics for neuron `g` (trial-level and condition-level)
- `[allrmap, allpf]`: reshaped population-level versions (via `fct_mapandfields_reshape`)

## 2) Data you must have before running

The script expects variables to already exist in MATLAB workspace (from your lab `.mat` files):

- From trajectory/session behavior:
  - `Traj` struct array with fields:
    - `start`, `stop` (trial boundaries, behavior sample index)
    - `Cond` (condition id)
    - `WB` (`'W'`/`'B'`, used as direction grouping)
  - `X_ds_n`: 1D position sampled at behavior frequency (script uses 1000 Hz convention)
- From ephys:
  - `allcel.id_cel`: unique cell ids
  - `allcel.id_spk`: cell id per spike
  - `allcel.itime_spk`: spike times in spike-sampling indices (25 kHz in this script)
- Session metadata:
  - `Namesession`
  - `OutputFolder`

The script header comment says to load `_TrajData.mat`, `_Phenosys.mat`, `_ePhy.mat` first.

## 3) Required code dependencies

Present in this repo:

- [fct_mapandfields.m](/c:/Users/tadse/OneDrive/Documenti/GitHub/CA1-cellclass/src/placefields/matlab/fct_mapandfields.m)
- [fct_rmap_wrap.m](/c:/Users/tadse/OneDrive/Documenti/GitHub/CA1-cellclass/src/placefields/matlab/fct_rmap_wrap.m)
- [fct_find_placefield.m](/c:/Users/tadse/OneDrive/Documenti/GitHub/CA1-cellclass/src/placefields/matlab/fct_find_placefield.m)
- [fct_ifreq_swap.m](/c:/Users/tadse/OneDrive/Documenti/GitHub/CA1-cellclass/src/placefields/matlab/fct_ifreq_swap.m)

Referenced but not included in this folder (you need them on MATLAB path):

- `fct_rmap`, `fct_rmap_features`
- `fct_placefield_simtrain`, `fct_placefield_threshold`, `fct_placefield_stability`
- `fct_placefield_peakorder`, `fct_placefield_keepoverlap`, `fct_placefield_features`
- `fct_mapandfields_reshape`
- (optional plotting) `fct_plot_rmap`
- utility/color scripts (`loadcolmathSFN`, `colormaprom.mat`)

## 4) Pipeline logic (rmap/pf stage)

## 4.1 Position normalization

`X_ds_n` is rescaled to `[0, 100]`:

```matlab
maxmaze = max(X_ds_n);
X_ds_n = X_ds_n * 100 / maxmaze;
```

Why: this standardizes spatial range across sessions so binning and field-size thresholds are comparable.

## 4.2 Build trial-to-condition-direction index (`condway`)

`Traj` trials are mapped to integer labels used by rate-map averaging:

- one label per `(condition, direction)` pair
- effectively two labels per condition (W/B)
- total labels `prm.rmap.nb_cond = nbsess * 2`

Why: `rmap` and `pf` are first computed per trial, then averaged/performed per condition-direction group.

## 4.3 Build trial boundaries (`indextr`)

For each trial `tt`, `indextr(tt,:) = [Traj(tt).start, Traj(tt).stop]`.

This becomes `prm.rmap.tbin` (time bins/trials used in map computation).

## 4.4 Parameters used for `rmap`

Main values set in `Rmap_G.m`:

- `prm.rmap.tbin = indextr`
- `prm.rmap.idcond_t = condway`
- `prm.rmap.freq = 1000` (behavior sampling frequency)
- `prm.rmap.xbin = 0:floor(maxmaze/100):maxmaze`
- `prm.rmap.xbin_rem = 10` (remove edge bins)
- `prm.rmap.ismooth = 10` (Gaussian smoothing half-width in bins)
- activity criteria copied into `prm.rmap`: `frmax`, `frmean`, `pctspklap`

Why these matter:

- `xbin`, `xbin_rem`, `ismooth` control spatial resolution and smoothness
- `tbin`/`idcond_t` define how trials are grouped
- `freq` ensures time/dwell computation is in correct units

## 4.5 Parameters used for `pf`

Main values set in `Rmap_G.m`:

- null model: `prm.pf.mtd='random'`, `prm.pf.nb_rep=1000`
- geometry constraints: `min_len=3`, `max_len=45`
- merge/extend rules: `max_btw_pf=3`, `max_ext_pf=5`
- significance: `pval_crit=0.01`, `pval_edge_pf=0.30`
- stability: `pctstab=0.4`, `stab=0.6`
- active-cell gating: `frmax`, `frmean`, `pctspklap` (copied from `rmap`)

Why: fields are detected by comparing observed maps to simulated maps, then filtered by significance, geometry, and lap stability.

## 4.6 Per-cell loop (main computation)

For each cell `g`:

1. Select spikes for that cell:
   - `idx_spk = allcel.itime_spk(allcel.id_spk == allcel.id_cel(g))`
2. Convert spike indices from 25 kHz to 1 kHz:
   - `idx_spk = fct_ifreq_swap(idx_spk, 25000, 1000)`
3. Compute maps and fields:
   - `[rmap(g), pf(g)] = fct_mapandfields(X_ds_n, idx_spk, prm)`

## 4.7 Inside `fct_mapandfields`

`fct_mapandfields` is just:

1. `rmap = fct_rmap_wrap(...)`
2. `pf   = fct_find_placefield(..., rmap, prm.pf)`

### `fct_rmap_wrap` does:

- calls `fct_rmap` to build per-trial matrices (`*_tx`):
  - spike counts (`nbspk_tx`)
  - dwell (`dwell_tx`)
  - firing rate (`fr_tx`)
  - smoothed versions (`*_s_tx`)
- averages trials into condition-direction matrices (`*_cx`) using `idcond_t`
- computes condition-wise ratemap features (`rmap.feat(c)`)

### `fct_find_placefield` does:

- builds simulated smoothed maps (`fr_s_txrep`) via `fct_placefield_simtrain`
- per trial:
  - computes binwise p-values vs null simulations
  - thresholds into `pf.ispf_tx`
- per condition-direction:
  - computes p-values on trial-averaged maps
  - thresholds + stability filtering + peak ordering
  - keeps lap fields overlapping condition field
  - extracts PF features into `pf.feat(c)`

## 4.8 Detailed place-field detection (bootstrap/null model)

This section expands the exact detection logic implemented in:

- [fct_placefield_simtrain.m](/c:/Users/tadse/OneDrive/Documenti/GitHub/CA1-cellclass/src/placefields/matlab/fct_placefield_simtrain.m)
- [fct_find_placefield.m](/c:/Users/tadse/OneDrive/Documenti/GitHub/CA1-cellclass/src/placefields/matlab/fct_find_placefield.m)

### Step A: Build null maps per lap (`fct_placefield_simtrain`)

For each trial/lap `t`, the function generates `nb_rep` surrogate maps under a no-place-coding null, while preserving trial-level spike load.

1. Preprocess position support:
   - Remove behavior samples listed in `prm.idx_rem` (`xpos(prm.idx_rem)=NaN`).
   - Remove edge spatial bins according to `xbin_rem`.
2. For each trial `t`:
   - Find valid time samples in that trial (`idxnonan`) that are non-NaN and inside spatial bin range.
   - Compute:
     - `N = number of valid samples`
     - `nbspk = total spikes on that trial` (`sum(rmap.nbspk_tx(t,:))`)
   - Edge case:
     - if `N == 0`, set `N=1`, `nbspk=0` so the simulation remains well-defined.
3. Simulate spike placement (`prm.mtd`):
   - `random`:
     - for each repetition, place `nbspk` spikes uniformly among the `N` valid samples
     - convert simulated spike positions to binned counts (`fct_hist`)
   - `poisson`:
     - for each repetition and valid sample, sample spike counts from `Poisson(nbspk/N)`
     - project counts to spatial bins (`fct_discretize` + `accumarray`)
4. Convert surrogate counts to surrogate smoothed rate maps:
   - smooth simulated counts (`fct_smoothgauss`)
   - divide by observed smoothed dwell (`rmap.dwell_s_tx(t,:)`)
   - smooth again
5. Store output:
   - `fr_s_txrep(t, x, rep)` = surrogate smoothed rate at bin `x`, repetition `rep`.

Interpretation:
- The null preserves trial duration/occupancy support and total spike count per trial, but destroys structured spike-position coupling expected from a real field.

### Step B: Trial-level p-values and binary fields (`fct_find_placefield`)

For each trial `t` and position bin `x`, compute a one-sided empirical p-value:

- `pval_pf_tx(t,x) = (# reps with fr_null(t,x,rep) >= fr_obs(t,x)) / nb_rep`

Then threshold p-values into a binary field mask:

- `ispf_tx(t,:) = fct_placefield_threshold(pval_pf_tx(t,:), prm)`

This threshold function applies:
- core significance (`pval_crit`)
- minimum/maximum field width (`min_len`, `max_len`)
- merge/extension rules (`max_btw_pf`, `max_ext_pf`, `pval_edge_pf`)

### Step C: Condition-level field from trial-averaged null

For each condition/direction `c`:

1. Select laps in that condition: `idx = (idcond_t == c)`.
2. Build condition-level null by averaging surrogate maps across selected laps:
   - `A(x,rep) = mean_t fr_s_txrep(t,x,rep)` over `t in idx`.
3. Compute condition-level p-value against observed condition map:
   - `pval_pf_cx(c,x) = (# reps with A(x,rep) >= fr_obs_c(c,x)) / nb_rep`
4. Threshold:
   - `ispf_cx(c,:) = fct_placefield_threshold(pval_pf_cx(c,:), prm)`

### Step D: Stability and consistency filters (condition level)

After significance thresholding, condition-level fields are refined:

1. Stability filter:
   - `ispf_cx(c,:) = fct_placefield_stability(rmap.fr_s_tx(idx,:), ispf_cx(c,:), prm)`
   - enforces consistency of trial maps with the condition field (`pctstab`, `stab` parameters).
2. Peak ordering:
   - `ispf_cx(c,:) = fct_placefield_peakorder(rmap.fr_s_cx(c,:), ispf_cx(c,:))`
   - orders fields by descending peak activity.
3. Trial-condition overlap constraint:
   - `pf.ispf_tx(idx,:) = fct_placefield_keepoverlap(pf.ispf_tx(idx,:), pf.ispf_cx(c,:))`
   - keeps only per-lap fields that overlap accepted condition-level fields.

### Step E: Feature extraction

For each condition `c`, final features are extracted from filtered trial/condition masks:

- `pf.feat(c) = fct_placefield_features(rmap.fr_s_tx(idx,:), pf.ispf_tx(idx,:), pf.ispf_cx(c,:))`

This produces the final per-field metrics used downstream (rate in/out, size, peak bin, lap variability, etc.).

### Notes on defaults currently used in code

`fct_find_placefield.m` contains updated defaults (different from the commented OLD block):

- `min_len = 2`
- `max_btw_pf = 3`
- `max_ext_pf = 3`
- `pval_crit = 0.05`
- `frmax = 1.5`
- `pctspklap = 0.45`

When reproducing analyses, always record the exact `pf.prm` saved with outputs, because these values strongly impact field detection rates.

## 5) Output structures at this stage

Per neuron (`rmap(g)`):

- `fr_tx`, `fr_s_tx`: trial-by-position maps
- `fr_cx`, `fr_s_cx`: condition-by-position maps
- plus `nbspk` and `dwell` counterparts

Per neuron (`pf(g)`):

- `pval_pf_tx`, `ispf_tx`: trial-level PF significance/mask
- `pval_pf_cx`, `ispf_cx`: condition-level PF significance/mask
- `feat(c)`: PF properties per condition-direction

Population reshape:

- `[allrmap, allpf] = fct_mapandfields_reshape(rmap, pf)`

Saved in final file:

- `<OutputFolder>/<Namesession>_Ratemap.mat` containing `allrmap`, `allpf`, `allcel`, `rmap`, `pf`

## 6) Practical checklist to reproduce only rmap/pf

1. Load required `.mat` files so `Traj`, `X_ds_n`, `allcel`, `Namesession`, `OutputFolder` exist.
2. Ensure missing helper functions (listed above) are on MATLAB path.
3. Confirm units:
   - `allcel.itime_spk` at 25 kHz
   - behavior indices (`Traj.start/stop`, `X_ds_n`) at 1 kHz
4. Run `Rmap_G.m`.
5. Inspect:
   - `rmap(g).fr_s_cx` for mean condition maps
   - `pf(g).ispf_cx` for detected condition-level fields
   - `allpf` after reshape for population analysis

## 7) Why this stage is needed

This stage converts raw trajectory + spikes into statistically filtered spatial tuning estimates:

- `rmap`: where each neuron fires along the track
- `pf`: whether elevated firing regions are significant/stable place fields

Without this stage, downstream labels like "spatially modulated", "uni/bidirectional", or condition comparisons are not grounded in a formal null-model detection step.

