# CA1-cellclass

Pipeline for CA1 neuron feature extraction and cell-type analysis across mouse age groups.

This repository has two main layers:

1. `src/cellclass`: data processing and feature extraction from raw/interim files.
2. `src/models`: unsupervised modeling, model-selection, stability analysis, and comparison with external labels.

---

## Project Structure

- `src/cellclass/`
  - Core library code for signal processing, ACG/waveform features, and IO utilities.
- `scripts/`
  - End-to-end orchestration scripts (interim -> processed, aggregation by mouse/age).
- `src/models/`
  - Modeling scripts:
    - `fitting.py`: systematic GMM model selection (`k`, feature subsets).
    - `compare_type_u.py`: fixed `k=2` clustering vs `allcel__type_u`.
    - `stability_analysis.py`: per-age-group robustness analysis across seeds/subsets.
- `data/`
  - Raw/interim/processed data.
- `results/`
  - Aggregated age-group datasets and modeling outputs.

---

## What `cellclass` Does

`cellclass` turns session-level spike/waveform data into processed per-neuron features.

Typical outputs include:

- firing/statistical features (`fr_hz`, `cv2`, `burst_index`)
- waveform features (`spk_duration_ms`, `spk_peaktrough_ms`, `spk_asymmetry`)
- ACG-based features (`refractory_*`, `acg_peak_latency_ms`)
- QC flags (`qc_*`)
- external prior label when present (`allcel__type_u`)

Main orchestration scripts:

1. `scripts/interim_to_processed.py`
   - Runs extraction on all interim `.npz` sessions.
   - Writes per-session files under `data/processed/<mouse>/...`.

2. `scripts/aggreggate_by_age.py`
   - Merges sessions with age metadata (`data/schedule.xlsx`).
   - Writes age-group datasets in `results/<AGE>/`.

---

## What `models` Does

`src/models` runs downstream unsupervised analyses on age-group datasets.

### 1) Model Selection (`fitting.py`)

- Tests multiple feature subsets and cluster counts.
- Computes metrics (`BIC`, `AIC`, `silhouette`, `Calinski-Harabasz`, `Davies-Bouldin`).
- Filters tiny clusters with `--min_cluster_size`.
- Selects best model using ordered rule:
  - lowest `BIC`
  - highest `silhouette`
  - lowest `Davies-Bouldin`
  - lowest `AIC`
  - highest `Calinski-Harabasz`

### 2) Label Comparison (`compare_type_u.py`)

- Runs fixed `k=2` GMM with all features.
- Maps clusters to `interneuron/pyramidal` via `spk_duration_ms`.
- Compares predictions with `allcel__type_u` (`0=interneuron`, `1=pyramidal`).
- Exports discrepancy summaries by neuron/session/age.

### 3) Stability (`stability_analysis.py`)

- Runs repeated `k=2` clustering for one age group over seeds/subsets.
- Quantifies run-to-run stability (pairwise ARI).
- Quantifies agreement with `type_u` (including flip-aware accuracy).
- Produces uncertainty ranking for neurons.

---

## Typical Workflow

From repository root:

1. Process interim sessions

```powershell
python scripts/interim_to_processed.py --interim_root data/interim --processed_root data/processed --skip_existing
```

2. Aggregate by age group

```powershell
python scripts/aggreggate_by_age.py --processed_root data/processed --excel data/schedule.xlsx --outdir results --qc
```

3. Run GMM model selection

```powershell
python src/models/fitting.py --results_root results --out_root results/model_selection --subset_mode leave_one_out --k_min 2 --k_max 6 --n_init 5 --min_cluster_size 3 --min_units 30
```

4. Compare fixed `k=2` model to `type_u`

```powershell
python src/models/compare_type_u.py --results_root results --out_root results/type_u_comparison --n_init 10
```

5. Run per-group stability analysis (example: `P23_24`)

```powershell
python src/models/stability_analysis.py --results_root results --out_root results/stability --age_group P23_24 --subset_mode leave_one_out --seeds 0,1,2,3,4,5,6,7,8,9 --n_init 10
```

---

## Key Outputs

- `results/<AGE>/...`
  - age-group clean/all units used for modeling.
- `results/model_selection/<AGE>/...`
  - full + ranked model selection tables and best model labels.
- `results/type_u_comparison/...`
  - agreement/discrepancy summaries vs `allcel__type_u`.
- `results/stability/<AGE>/...`
  - run-level stability metrics, pairwise ARI, uncertain neurons.

---

## Notes

- Cluster IDs are arbitrary; when comparing binary labels, use label-invariant metrics (e.g., ARI) or flip-aware accuracy.
- Low silhouette values indicate overlapping classes in current feature space, not necessarily a bug.
- For details on modeling options, see `src/models/README.md`.

