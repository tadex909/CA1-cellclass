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

