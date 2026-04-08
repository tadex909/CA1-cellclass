from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .metrics import smooth_last_axis
from .trials import TrialInfo


@dataclass(frozen=True)
class RatemapConfig:
    """
    Parameters for MATLAB-like ratemap construction.
    """
    freq_hz: float = 1000.0
    smooth_sigma_bins: float = 10.0
    xbin_rem: int = 0
    nb_cond: int | None = None
    min_speed: float | None = 2.0

    def validate(self) -> None:
        if self.freq_hz <= 0:
            raise ValueError("freq_hz must be > 0")
        if self.smooth_sigma_bins < 0:
            raise ValueError("smooth_sigma_bins must be >= 0")
        if self.xbin_rem < 0:
            raise ValueError("xbin_rem must be >= 0")
        if self.nb_cond is not None and self.nb_cond < 1:
            raise ValueError("nb_cond must be >= 1 when provided")
        if self.min_speed is not None and not np.isfinite(self.min_speed):
            raise ValueError("min_speed must be finite when provided")


@dataclass(frozen=True)
class RatemapPack:
    """
    Session ratemap tensors, analogous to MATLAB `rmap` outputs.

    Shapes
    ------
    - `nbspk_tx_ux`, `fr_tx_ux`, `nbspk_s_tx_ux`, `fr_s_tx_ux`: (n_cells, n_trials, n_bins)
    - `dwell_tx_x`, `dwell_s_tx_x`: (n_trials, n_bins)
    - `*_cx_ux`: (n_cells, n_cond, n_bins)
    - `dwell_cx_x`, `dwell_s_cx_x`: (n_cond, n_bins)
    """
    cell_ids: np.ndarray
    trial_info: list[TrialInfo]
    idcond_t: np.ndarray
    xbin_edges: np.ndarray
    xbin_centers: np.ndarray
    nb_cond: int
    nbspk_tx_ux: np.ndarray
    dwell_tx_x: np.ndarray
    fr_tx_ux: np.ndarray
    nbspk_s_tx_ux: np.ndarray
    dwell_s_tx_x: np.ndarray
    fr_s_tx_ux: np.ndarray
    nbspk_cx_ux: np.ndarray
    dwell_cx_x: np.ndarray
    fr_cx_ux: np.ndarray
    nbspk_s_cx_ux: np.ndarray
    dwell_s_cx_x: np.ndarray
    fr_s_cx_ux: np.ndarray

    def cell_index(self, cell_id: int) -> int:
        idx = np.where(self.cell_ids.astype(np.int64) == int(cell_id))[0]
        if idx.size == 0:
            raise KeyError(f"cell_id {cell_id} not found")
        return int(idx[0])

    def cell_view(self, *, cell_id: int | None = None, cell_index: int | None = None) -> dict[str, np.ndarray]:
        """
        Return one-cell view with MATLAB-like key names (`*_tx`, `*_cx`).
        """
        if cell_index is None:
            if cell_id is None:
                raise ValueError("Provide either cell_id or cell_index")
            cell_index = self.cell_index(cell_id)
        u = int(cell_index)
        return {
            "nbspk_tx": self.nbspk_tx_ux[u],
            "dwell_tx": self.dwell_tx_x,
            "fr_tx": self.fr_tx_ux[u],
            "nbspk_s_tx": self.nbspk_s_tx_ux[u],
            "dwell_s_tx": self.dwell_s_tx_x,
            "fr_s_tx": self.fr_s_tx_ux[u],
            "nbspk_cx": self.nbspk_cx_ux[u],
            "dwell_cx": self.dwell_cx_x,
            "fr_cx": self.fr_cx_ux[u],
            "nbspk_s_cx": self.nbspk_s_cx_ux[u],
            "dwell_s_cx": self.dwell_s_cx_x,
            "fr_s_cx": self.fr_s_cx_ux[u],
        }


def normalize_x_to_100(position_x: np.ndarray) -> np.ndarray:
    """
    Mirror the MATLAB step:
        X_ds_n = X_ds_n * 100 / max(X_ds_n)
    """
    x = np.asarray(position_x, dtype=np.float64).ravel()
    finite = np.isfinite(x)
    if not np.any(finite):
        raise ValueError("position_x has no finite values")
    mx = float(np.nanmax(x[finite]))
    if mx <= 0:
        raise ValueError("max(position_x) must be > 0 for normalization")
    out = x.copy()
    out[finite] = out[finite] * 100.0 / mx
    return out


def build_default_xbin(position_x: np.ndarray) -> np.ndarray:
    """
    MATLAB-like default:
        xbin = 0:floor(maxmaze/100):maxmaze
    """
    x = np.asarray(position_x, dtype=np.float64).ravel()
    finite = np.isfinite(x)
    if not np.any(finite):
        raise ValueError("position_x has no finite values")
    maxmaze = float(np.nanmax(x[finite]))
    step = int(np.floor(maxmaze / 100.0))
    step = max(1, step)
    edges = np.arange(0.0, maxmaze + step, step, dtype=np.float64)
    if edges[-1] < maxmaze:
        edges = np.append(edges, maxmaze)
    return edges


def build_ratemap_from_trials(
    *,
    position_x: np.ndarray,
    spike_indices_0b: np.ndarray,
    spike_cell_ids: np.ndarray,
    cell_ids: np.ndarray,
    trials: list[TrialInfo],
    xbin_edges: np.ndarray,
    cfg: RatemapConfig,
    speed: np.ndarray | None = None,
) -> RatemapPack:
    """
    Build ratemap tensors (`*_tx` and `*_cx`) in a MATLAB-compatible layout.

    Parameters
    ----------
    position_x
        1D position over time at behavior frequency (typically 1000 Hz).
    spike_indices_0b
        Spike sample indices in Python 0-based indexing at the same frequency as `position_x`.
    spike_cell_ids
        Cell id per spike, same length as `spike_indices_0b`.
    cell_ids
        Ordered unique cell ids.
    trials
        Per-trial descriptors (`start` inclusive, `stop` exclusive).
    xbin_edges
        Spatial bin edges.
    cfg
        Ratemap parameters.
    speed
        Optional speed vector aligned to `position_x`; when provided with `cfg.min_speed`,
        samples/spikes below threshold are removed.
    """
    cfg.validate()
    eps = np.finfo(np.float64).eps

    x = np.asarray(position_x, dtype=np.float64).ravel()
    n_samples = x.size
    if n_samples == 0:
        raise ValueError("position_x is empty")

    spk_idx = np.asarray(spike_indices_0b, dtype=np.int64).ravel()
    spk_cid = np.asarray(spike_cell_ids, dtype=np.int64).ravel()
    if spk_idx.size != spk_cid.size:
        raise ValueError(
            f"spike_indices_0b and spike_cell_ids size mismatch: {spk_idx.size} vs {spk_cid.size}"
        )

    cids = np.asarray(cell_ids, dtype=np.int64).ravel()
    if cids.size == 0:
        raise ValueError("cell_ids is empty")
    cell_to_u = {int(cid): i for i, cid in enumerate(cids)}

    edges_full = np.asarray(xbin_edges, dtype=np.float64).ravel()
    if edges_full.size < 2:
        raise ValueError("xbin_edges must contain at least 2 edges")
    if not np.all(np.diff(edges_full) > 0):
        raise ValueError("xbin_edges must be strictly increasing")
    n_bins_full = edges_full.size - 1
    if cfg.xbin_rem * 2 >= n_bins_full:
        raise ValueError(
            f"xbin_rem={cfg.xbin_rem} removes too many bins for n_bins={n_bins_full}"
        )

    if speed is not None:
        speed_arr = np.asarray(speed, dtype=np.float64).ravel()
        if speed_arr.size != n_samples:
            raise ValueError(
                f"speed length mismatch: expected {n_samples}, got {speed_arr.size}"
            )
    else:
        speed_arr = None

    sample_keep = np.isfinite(x)
    if (speed_arr is not None) and (cfg.min_speed is not None):
        sample_keep &= np.isfinite(speed_arr) & (speed_arr >= float(cfg.min_speed))

    n_cells = cids.size
    n_trials = len(trials)

    nbspk_tx_ux = np.zeros((n_cells, n_trials, n_bins_full), dtype=np.float64)
    dwell_tx_x = np.zeros((n_trials, n_bins_full), dtype=np.float64)

    for t, tr in enumerate(trials):
        s = max(0, int(tr.start_idx_0b))
        e = min(n_samples, int(tr.stop_idx_0b_exclusive))
        if e <= s:
            continue

        keep_t = sample_keep[s:e]
        x_t = x[s:e][keep_t]
        if x_t.size:
            occ_counts, _ = np.histogram(x_t, bins=edges_full)
            dwell_tx_x[t] = occ_counts.astype(np.float64) / float(cfg.freq_hz)

        in_trial = (spk_idx >= s) & (spk_idx < e)
        if not np.any(in_trial):
            continue

        spk_idx_t = spk_idx[in_trial]
        spk_cid_t = spk_cid[in_trial]

        valid_idx = (spk_idx_t >= 0) & (spk_idx_t < n_samples)
        if not np.any(valid_idx):
            continue
        spk_idx_t = spk_idx_t[valid_idx]
        spk_cid_t = spk_cid_t[valid_idx]

        keep_spk = sample_keep[spk_idx_t]
        if not np.any(keep_spk):
            continue
        spk_idx_t = spk_idx_t[keep_spk]
        spk_cid_t = spk_cid_t[keep_spk]

        spk_x_t = x[spk_idx_t]
        finite_spk_x = np.isfinite(spk_x_t)
        if not np.any(finite_spk_x):
            continue
        spk_x_t = spk_x_t[finite_spk_x]
        spk_cid_t = spk_cid_t[finite_spk_x]

        for cid in np.unique(spk_cid_t):
            u = cell_to_u.get(int(cid))
            if u is None:
                continue
            vals = spk_x_t[spk_cid_t == cid]
            if vals.size == 0:
                continue
            c, _ = np.histogram(vals, bins=edges_full)
            nbspk_tx_ux[u, t] = c.astype(np.float64)

    # MATLAB fct_rmap adds eps to dwell before computing firing rate.
    dwell_tx_x = dwell_tx_x + eps
    fr_tx_ux = nbspk_tx_ux / dwell_tx_x[None, :, :]

    nbspk_s_tx_ux = smooth_last_axis(nbspk_tx_ux, cfg.smooth_sigma_bins)
    dwell_s_tx_x = smooth_last_axis(dwell_tx_x, cfg.smooth_sigma_bins) + eps
    fr_s_tx_ux = nbspk_s_tx_ux / dwell_s_tx_x[None, :, :]

    if cfg.xbin_rem > 0:
        r = int(cfg.xbin_rem)
        bsl = slice(r, n_bins_full - r)
        edges = edges_full[r : (n_bins_full - r + 1)]
        nbspk_tx_ux = nbspk_tx_ux[:, :, bsl]
        dwell_tx_x = dwell_tx_x[:, bsl]
        fr_tx_ux = fr_tx_ux[:, :, bsl]
        nbspk_s_tx_ux = nbspk_s_tx_ux[:, :, bsl]
        dwell_s_tx_x = dwell_s_tx_x[:, bsl]
        fr_s_tx_ux = fr_s_tx_ux[:, :, bsl]
    else:
        edges = edges_full

    idcond_t = np.asarray([tr.condway for tr in trials], dtype=np.int64)
    if idcond_t.size == 0:
        nb_cond = int(cfg.nb_cond or 1)
    else:
        nb_cond = int(cfg.nb_cond or int(np.nanmax(idcond_t)))
    n_bins = edges.size - 1

    nbspk_cx_ux = np.full((n_cells, nb_cond, n_bins), np.nan, dtype=np.float64)
    dwell_cx_x = np.full((nb_cond, n_bins), np.nan, dtype=np.float64)
    fr_cx_ux = np.full((n_cells, nb_cond, n_bins), np.nan, dtype=np.float64)
    nbspk_s_cx_ux = np.full((n_cells, nb_cond, n_bins), np.nan, dtype=np.float64)
    dwell_s_cx_x = np.full((nb_cond, n_bins), np.nan, dtype=np.float64)
    fr_s_cx_ux = np.full((n_cells, nb_cond, n_bins), np.nan, dtype=np.float64)

    for c in range(1, nb_cond + 1):
        idx = idcond_t == c
        if not np.any(idx):
            continue
        nbspk_cx_ux[:, c - 1, :] = np.nanmean(nbspk_tx_ux[:, idx, :], axis=1)
        dwell_cx_x[c - 1, :] = np.nanmean(dwell_tx_x[idx, :], axis=0)
        nbspk_s_cx_ux[:, c - 1, :] = np.nanmean(nbspk_s_tx_ux[:, idx, :], axis=1)
        dwell_s_cx_x[c - 1, :] = np.nanmean(dwell_s_tx_x[idx, :], axis=0)
        fr_cx_ux[:, c - 1, :] = np.nanmean(fr_tx_ux[:, idx, :], axis=1)
        fr_s_cx_ux[:, c - 1, :] = np.nanmean(fr_s_tx_ux[:, idx, :], axis=1)

    centers = 0.5 * (edges[:-1] + edges[1:])
    return RatemapPack(
        cell_ids=cids.copy(),
        trial_info=list(trials),
        idcond_t=idcond_t,
        xbin_edges=edges,
        xbin_centers=centers,
        nb_cond=nb_cond,
        nbspk_tx_ux=nbspk_tx_ux,
        dwell_tx_x=dwell_tx_x,
        fr_tx_ux=fr_tx_ux,
        nbspk_s_tx_ux=nbspk_s_tx_ux,
        dwell_s_tx_x=dwell_s_tx_x,
        fr_s_tx_ux=fr_s_tx_ux,
        nbspk_cx_ux=nbspk_cx_ux,
        dwell_cx_x=dwell_cx_x,
        fr_cx_ux=fr_cx_ux,
        nbspk_s_cx_ux=nbspk_s_cx_ux,
        dwell_s_cx_x=dwell_s_cx_x,
        fr_s_cx_ux=fr_s_cx_ux,
    )
