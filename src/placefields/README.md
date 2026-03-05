# Placefields

Minimal scaffold for place-field analysis inside `CA1-cellclass`.

## Current Scope

- Bin-based occupancy and firing-rate maps.
- Gaussian smoothing on 1D rate maps.
- Skaggs-style spatial information (bits/spike).
- Simple place-field detection via thresholded contiguous bins.
- Session-level runner script: `scripts/placefield_analysis.py`.

## Expected Input Schema (`.npz`)

Each input session file must provide:

- `position_cm`: `(n_samples,)`
- `dt_s`: `(n_samples,)`
- `spike_positions_cm`: `(n_spikes,)`
- `spike_cell_ids`: `(n_spikes,)`
- `cell_ids`: `(n_cells,)`

## Example Run

```powershell
python scripts/placefield_analysis.py `
  --input_root data/placefields/interim `
  --out_root results/placefields `
  --bin_size_cm 2.0 `
  --min_occupancy_s 0.1 `
  --smooth_sigma_bins 1.5 `
  --field_threshold_ratio 0.2 `
  --min_field_bins 3
```

## Next Integration Step

Map your existing interim/processed structures to this input schema so place-field analysis
can run directly on your current sessions.
