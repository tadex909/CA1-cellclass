from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal, Mapping, Sequence

import numpy as np

from .metrics import smooth_last_axis
from .trials import TrialInfo, canonical_condition_name, condition_family_name


DecodeGroupBy = Literal["condway", "condition", "global"]


@dataclass(frozen=True)
class BayesianDecoderConfig:
    """
    Parameters for memoryless Bayesian position decoding.
    """

    freq_hz: float = 1000.0
    tau_s: float = 0.150
    bin_size_cm: float = 2.0
    min_speed: float | None = 2.0
    smooth_sigma_bins: float = 2.8
    xbin_rem: int = 0
    min_valid_window_fraction: float = 0.5
    rate_floor_hz: float = 1e-12
    decode_groupby: DecodeGroupBy = "condway"
    group_condition_families: bool = False
    train_all_laps: bool = False

    def validate(self) -> None:
        if not np.isfinite(self.freq_hz) or self.freq_hz <= 0:
            raise ValueError("freq_hz must be finite and > 0")
        if not np.isfinite(self.tau_s) or self.tau_s <= 0:
            raise ValueError("tau_s must be finite and > 0")
        if not np.isfinite(self.bin_size_cm) or self.bin_size_cm <= 0:
            raise ValueError("bin_size_cm must be finite and > 0")
        if self.min_speed is not None and not np.isfinite(self.min_speed):
            raise ValueError("min_speed must be finite when provided")
        if not np.isfinite(self.smooth_sigma_bins) or self.smooth_sigma_bins < 0:
            raise ValueError("smooth_sigma_bins must be finite and >= 0")
        if self.xbin_rem < 0:
            raise ValueError("xbin_rem must be >= 0")
        if not (
            np.isfinite(self.min_valid_window_fraction)
            and 0.0 < self.min_valid_window_fraction <= 1.0
        ):
            raise ValueError("min_valid_window_fraction must be in (0, 1]")
        if not np.isfinite(self.rate_floor_hz) or self.rate_floor_hz <= 0:
            raise ValueError("rate_floor_hz must be finite and > 0")
        if self.decode_groupby not in {"condway", "condition", "global"}:
            raise ValueError("decode_groupby must be one of {'condway','condition','global'}")


@dataclass(frozen=True)
class BayesianDecodingResult:
    cell_ids: np.ndarray
    idcond_t: np.ndarray
    xbin_edges: np.ndarray
    xbin_centers: np.ndarray
    posterior_wx: np.ndarray
    decoded_bin_w: np.ndarray
    decoded_x_w: np.ndarray
    actual_bin_w: np.ndarray
    actual_x_w: np.ndarray
    error_cm_w: np.ndarray
    prob_actual_w: np.ndarray
    trial_index_w: np.ndarray
    condway_w: np.ndarray
    train_group_w: np.ndarray
    train_group_label_w: np.ndarray
    window_start_w: np.ndarray
    window_stop_w: np.ndarray
    n_spikes_w: np.ndarray
    n_train_laps_w: np.ndarray

    def to_payload(self) -> dict[str, Any]:
        return {
            "cell_ids": self.cell_ids.astype(np.int64, copy=False),
            "idcond_t": self.idcond_t.astype(np.int64, copy=False),
            "xbin_edges": self.xbin_edges.astype(np.float64, copy=False),
            "xbin_centers": self.xbin_centers.astype(np.float64, copy=False),
            "decode__posterior_wx": self.posterior_wx.astype(np.float32, copy=False),
            "decode__decoded_bin_w": self.decoded_bin_w.astype(np.int64, copy=False),
            "decode__decoded_x_w": self.decoded_x_w.astype(np.float32, copy=False),
            "decode__actual_bin_w": self.actual_bin_w.astype(np.int64, copy=False),
            "decode__actual_x_w": self.actual_x_w.astype(np.float32, copy=False),
            "decode__error_cm_w": self.error_cm_w.astype(np.float32, copy=False),
            "decode__prob_actual_w": self.prob_actual_w.astype(np.float32, copy=False),
            "decode__trial_index_w": self.trial_index_w.astype(np.int64, copy=False),
            "decode__condway_w": self.condway_w.astype(np.int64, copy=False),
            "decode__train_group_w": self.train_group_w.astype(np.int64, copy=False),
            "decode__train_group_label_w": self.train_group_label_w.astype(np.str_, copy=False),
            "decode__window_start_w": self.window_start_w.astype(np.int64, copy=False),
            "decode__window_stop_w": self.window_stop_w.astype(np.int64, copy=False),
            "decode__n_spikes_w": self.n_spikes_w.astype(np.int64, copy=False),
            "decode__n_train_laps_w": self.n_train_laps_w.astype(np.int64, copy=False),
        }

    def summary_rows(self, *, session_id: str) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        if self.condway_w.size == 0:
            return rows

        for condway in sorted(set(self.condway_w.astype(np.int64).tolist())):
            idx = self.condway_w == int(condway)
            errors = self.error_cm_w[idx]
            probs = self.prob_actual_w[idx]
            rows.append(
                {
                    "session_id": str(session_id),
                    "condway": int(condway),
                    "train_group": _format_unique_ints(self.train_group_w[idx]),
                    "train_group_label": _format_unique_strings(self.train_group_label_w[idx]),
                    "n_windows": int(np.sum(idx)),
                    "n_cells_used": int(self.cell_ids.size),
                    "n_laps_decoded": int(np.unique(self.trial_index_w[idx]).size),
                    "median_error_cm": _nan_stat(errors, np.nanmedian),
                    "mean_error_cm": _nan_stat(errors, np.nanmean),
                    "mean_prob_actual": _nan_stat(probs, np.nanmean),
                }
            )
        return rows


def _nan_stat(arr: np.ndarray, fn: Any) -> float:
    x = np.asarray(arr, dtype=np.float64)
    if x.size == 0 or not np.any(np.isfinite(x)):
        return float("nan")
    return float(fn(x))


def _format_unique_ints(arr: np.ndarray) -> str:
    vals = sorted(set(np.asarray(arr, dtype=np.int64).ravel().tolist()))
    return ",".join(str(v) for v in vals)


def _format_unique_strings(arr: np.ndarray) -> str:
    vals = sorted(set(np.asarray(arr).astype(str).ravel().tolist()))
    return ",".join(v for v in vals if v)


def build_regular_xbin(position_x: np.ndarray, bin_size_cm: float = 2.0) -> np.ndarray:
    """
    Build regular spatial bin edges, defaulting to 2 cm bins over the track.
    """

    if not np.isfinite(bin_size_cm) or bin_size_cm <= 0:
        raise ValueError("bin_size_cm must be finite and > 0")
    x = np.asarray(position_x, dtype=np.float64).ravel()
    finite = x[np.isfinite(x)]
    if finite.size == 0:
        raise ValueError("position_x has no finite values")

    lo_raw = float(np.nanmin(finite))
    hi_raw = float(np.nanmax(finite))
    lo = 0.0 if lo_raw >= 0 else np.floor(lo_raw / bin_size_cm) * bin_size_cm
    hi = np.ceil(hi_raw / bin_size_cm) * bin_size_cm
    if hi <= lo:
        hi = lo + float(bin_size_cm)

    edges = np.arange(lo, hi + bin_size_cm * 0.5, bin_size_cm, dtype=np.float64)
    if edges.size < 2:
        edges = np.array([lo, lo + float(bin_size_cm)], dtype=np.float64)
    if edges[-1] < hi_raw:
        edges = np.append(edges, edges[-1] + float(bin_size_cm))
    return edges


def _bin_indices(x: np.ndarray, edges: np.ndarray) -> np.ndarray:
    vals = np.asarray(x, dtype=np.float64)
    bins = np.searchsorted(edges, vals, side="right") - 1
    bins = bins.astype(np.int64, copy=False)
    bins[vals == edges[-1]] = int(edges.size - 2)
    return bins


def _apply_xbin_rem(xbin_edges: np.ndarray, xbin_rem: int) -> np.ndarray:
    edges_full = np.asarray(xbin_edges, dtype=np.float64).ravel()
    if edges_full.size < 2:
        raise ValueError("xbin_edges must contain at least 2 edges")
    if not np.all(np.diff(edges_full) > 0):
        raise ValueError("xbin_edges must be strictly increasing")

    n_bins_full = int(edges_full.size - 1)
    if int(xbin_rem) * 2 >= n_bins_full:
        raise ValueError(
            f"xbin_rem={xbin_rem} removes too many bins for n_bins={n_bins_full}"
        )
    if xbin_rem <= 0:
        return edges_full.copy()

    r = int(xbin_rem)
    return edges_full[r : (n_bins_full - r + 1)].copy()


def _compute_trial_maps(
    *,
    position_x: np.ndarray,
    spike_indices_0b: np.ndarray,
    spike_cell_ids: np.ndarray,
    cell_ids: np.ndarray,
    trials: Sequence[TrialInfo],
    xbin_edges: np.ndarray,
    sample_keep: np.ndarray,
    freq_hz: float,
) -> tuple[np.ndarray, np.ndarray]:
    n_cells = int(cell_ids.size)
    n_trials = len(trials)
    n_bins = int(xbin_edges.size - 1)
    counts_utx = np.zeros((n_cells, n_trials, n_bins), dtype=np.float64)
    dwell_tx = np.zeros((n_trials, n_bins), dtype=np.float64)
    cell_to_u = {int(cid): u for u, cid in enumerate(cell_ids.astype(np.int64))}

    x = np.asarray(position_x, dtype=np.float64)
    spk_idx = np.asarray(spike_indices_0b, dtype=np.int64)
    spk_cid = np.asarray(spike_cell_ids, dtype=np.int64)

    for t, tr in enumerate(trials):
        s = max(0, int(tr.start_idx_0b))
        e = min(int(x.size), int(tr.stop_idx_0b_exclusive))
        if e <= s:
            continue

        keep_t = sample_keep[s:e]
        x_t = x[s:e][keep_t]
        if x_t.size:
            dwell_tx[t], _ = np.histogram(x_t, bins=xbin_edges)
            dwell_tx[t] /= float(freq_hz)

        in_trial = (spk_idx >= s) & (spk_idx < e)
        if not np.any(in_trial):
            continue
        idx_t = spk_idx[in_trial]
        cid_t = spk_cid[in_trial]
        keep_spk = sample_keep[idx_t]
        if not np.any(keep_spk):
            continue
        idx_t = idx_t[keep_spk]
        cid_t = cid_t[keep_spk]

        bins = _bin_indices(x[idx_t], xbin_edges)
        in_bins = (bins >= 0) & (bins < n_bins)
        if not np.any(in_bins):
            continue
        bins = bins[in_bins]
        cid_t = cid_t[in_bins]

        for cid in np.unique(cid_t):
            u = cell_to_u.get(int(cid))
            if u is None:
                continue
            counts_utx[u, t, :] = np.bincount(
                bins[cid_t == cid],
                minlength=n_bins,
            ).astype(np.float64, copy=False)

    return counts_utx, dwell_tx


def _pooled_tuning_rates(
    *,
    counts_utx: np.ndarray,
    dwell_tx: np.ndarray,
    train_idx: np.ndarray,
    smooth_sigma_bins: float,
    rate_floor_hz: float,
) -> tuple[np.ndarray, np.ndarray]:
    counts_ux = np.nansum(counts_utx[:, train_idx, :], axis=1)
    dwell_x = np.nansum(dwell_tx[train_idx, :], axis=0)

    counts_s_ux = smooth_last_axis(counts_ux, smooth_sigma_bins)
    dwell_s_x = smooth_last_axis(dwell_x, smooth_sigma_bins)
    valid_x = np.isfinite(dwell_s_x) & (dwell_s_x > 0)

    rates_ux = np.full_like(counts_s_ux, np.nan, dtype=np.float64)
    np.divide(
        counts_s_ux,
        dwell_s_x[None, :],
        out=rates_ux,
        where=valid_x[None, :],
    )
    rates_ux[:, valid_x] = np.maximum(rates_ux[:, valid_x], float(rate_floor_hz))
    rates_ux[:, ~valid_x] = np.nan
    return rates_ux, valid_x


def _posterior_from_counts(
    *,
    spike_counts_u: np.ndarray,
    rates_ux: np.ndarray,
    valid_x: np.ndarray,
    tau_s: float,
) -> np.ndarray | None:
    n_bins = int(rates_ux.shape[1])
    valid = np.asarray(valid_x, dtype=bool).ravel()
    if valid.size != n_bins or not np.any(valid):
        return None

    rates_v = rates_ux[:, valid]
    if not np.all(np.isfinite(rates_v)):
        return None

    n_valid = int(np.sum(valid))
    logp_v = (
        -np.log(float(n_valid))
        + spike_counts_u.astype(np.float64) @ np.log(rates_v)
        - float(tau_s) * np.sum(rates_v, axis=0)
    )
    m = float(np.max(logp_v))
    if not np.isfinite(m):
        return None

    p_v = np.exp(logp_v - m)
    total = float(np.sum(p_v))
    if not np.isfinite(total) or total <= 0:
        return None

    posterior = np.zeros(n_bins, dtype=np.float64)
    posterior[valid] = p_v / total
    return posterior


def _normalize_condition_name_map(
    condition_names_by_base: Mapping[int, str] | None,
) -> dict[int, str]:
    if condition_names_by_base is None:
        return {}
    return {
        int(k): canonical_condition_name(v)
        for k, v in condition_names_by_base.items()
        if canonical_condition_name(v)
    }


def _decode_group_info(
    trials: Sequence[TrialInfo],
    groupby: DecodeGroupBy,
    *,
    group_condition_families: bool,
    condition_names_by_base: Mapping[int, str] | None,
) -> tuple[np.ndarray, np.ndarray]:
    name_map = _normalize_condition_name_map(condition_names_by_base)

    if groupby == "global":
        return (
            np.ones(len(trials), dtype=np.int64),
            np.full(len(trials), "global", dtype=np.str_),
        )

    if not group_condition_families:
        if groupby == "condway":
            group_ids = np.asarray([tr.condway for tr in trials], dtype=np.int64)
            labels = []
            for tr in trials:
                name = name_map.get(int(tr.cond), f"cond{int(tr.cond)}")
                labels.append(f"{name} {str(tr.wb).upper()}")
            return group_ids, np.asarray(labels, dtype=np.str_)
        if groupby == "condition":
            group_ids = np.asarray([tr.cond for tr in trials], dtype=np.int64)
            labels = [name_map.get(int(tr.cond), f"cond{int(tr.cond)}") for tr in trials]
            return group_ids, np.asarray(labels, dtype=np.str_)
        raise ValueError(f"Unsupported decode_groupby: {groupby}")

    missing = sorted({int(tr.cond) for tr in trials if int(tr.cond) not in name_map})
    if missing:
        raise ValueError(
            "condition_names_by_base is required for family grouping and is missing "
            f"base condition ids: {missing}"
        )

    keys: list[tuple[str, ...]] = []
    labels: list[str] = []
    for tr in trials:
        family = condition_family_name(name_map[int(tr.cond)])
        direction = str(tr.wb).upper()
        if groupby == "condway":
            key = (family, direction)
            label = f"{family} {direction}"
        elif groupby == "condition":
            key = (family,)
            label = family
        else:
            raise ValueError(f"Unsupported decode_groupby: {groupby}")
        keys.append(key)
        labels.append(label)

    group_id_by_key = {key: i + 1 for i, key in enumerate(sorted(set(keys)))}
    group_ids = np.asarray([group_id_by_key[key] for key in keys], dtype=np.int64)
    return group_ids, np.asarray(labels, dtype=np.str_)


def decode_bayesian_position_from_trials(
    *,
    position_x: np.ndarray,
    spike_indices_0b: np.ndarray,
    spike_cell_ids: np.ndarray,
    cell_ids: np.ndarray,
    trials: Sequence[TrialInfo],
    xbin_edges: np.ndarray,
    cfg: BayesianDecoderConfig,
    speed: np.ndarray | None = None,
    condition_names_by_base: Mapping[int, str] | None = None,
) -> BayesianDecodingResult:
    """
    Decode position with a memoryless Bayesian decoder.

    By default, each lap is decoded using maps trained on other laps in the
    same decode group. Set ``cfg.train_all_laps=True`` to train on all laps in
    the same group, including the lap being decoded.
    """

    cfg.validate()
    x = np.asarray(position_x, dtype=np.float64).ravel()
    n_samples = int(x.size)
    if n_samples == 0:
        raise ValueError("position_x is empty")

    spk_idx = np.asarray(spike_indices_0b, dtype=np.int64).ravel()
    spk_cid = np.asarray(spike_cell_ids, dtype=np.int64).ravel()
    if spk_idx.size != spk_cid.size:
        raise ValueError(
            f"spike_indices_0b and spike_cell_ids size mismatch: {spk_idx.size} vs {spk_cid.size}"
        )

    cids_requested = np.asarray(cell_ids, dtype=np.int64).ravel()
    if cids_requested.size == 0:
        raise ValueError("cell_ids is empty")

    edges = _apply_xbin_rem(xbin_edges, cfg.xbin_rem)
    n_bins = int(edges.size - 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    if speed is None and cfg.min_speed is not None:
        raise ValueError("speed is required when min_speed is provided")
    if speed is not None:
        speed_arr = np.asarray(speed, dtype=np.float64).ravel()
        if speed_arr.size != n_samples:
            raise ValueError(f"speed length mismatch: expected {n_samples}, got {speed_arr.size}")
    else:
        speed_arr = None

    sample_keep = np.isfinite(x) & (x >= edges[0]) & (x <= edges[-1])
    if speed_arr is not None and cfg.min_speed is not None:
        sample_keep &= np.isfinite(speed_arr) & (speed_arr >= float(cfg.min_speed))

    valid_spikes = (spk_idx >= 0) & (spk_idx < n_samples)
    sample_valid_spikes = np.zeros_like(valid_spikes, dtype=bool)
    if np.any(valid_spikes):
        sample_valid_spikes[valid_spikes] = sample_keep[spk_idx[valid_spikes]]
    valid_spikes &= sample_valid_spikes
    valid_spikes &= np.isin(spk_cid, cids_requested)
    active_ids = np.unique(spk_cid[valid_spikes]).astype(np.int64, copy=False)
    keep_cells = np.isin(cids_requested, active_ids)
    cids = cids_requested[keep_cells]
    if cids.size == 0:
        raise ValueError("No active cells after applying position/speed filters")

    spk_idx_keep = spk_idx[valid_spikes]
    spk_cid_keep = spk_cid[valid_spikes]
    in_active = np.isin(spk_cid_keep, cids)
    spk_idx_keep = spk_idx_keep[in_active]
    spk_cid_keep = spk_cid_keep[in_active]

    if spk_idx_keep.size:
        order = np.argsort(spk_idx_keep, kind="stable")
        spk_idx_keep = spk_idx_keep[order]
        spk_cid_keep = spk_cid_keep[order]

    idcond_t = np.asarray([tr.condway for tr in trials], dtype=np.int64)
    train_group_t, train_group_label_t = _decode_group_info(
        trials,
        cfg.decode_groupby,
        group_condition_families=bool(cfg.group_condition_families),
        condition_names_by_base=condition_names_by_base,
    )
    counts_utx, dwell_tx = _compute_trial_maps(
        position_x=x,
        spike_indices_0b=spk_idx_keep,
        spike_cell_ids=spk_cid_keep,
        cell_ids=cids,
        trials=trials,
        xbin_edges=edges,
        sample_keep=sample_keep,
        freq_hz=float(cfg.freq_hz),
    )

    cell_to_u = {int(cid): u for u, cid in enumerate(cids.astype(np.int64))}
    tau_samples = max(1, int(round(float(cfg.tau_s) * float(cfg.freq_hz))))

    posterior_rows: list[np.ndarray] = []
    decoded_bin: list[int] = []
    decoded_x: list[float] = []
    actual_bin: list[int] = []
    actual_x: list[float] = []
    error_cm: list[float] = []
    prob_actual: list[float] = []
    trial_index: list[int] = []
    condway_w: list[int] = []
    train_group_w: list[int] = []
    train_group_label_w: list[str] = []
    window_start: list[int] = []
    window_stop: list[int] = []
    n_spikes_w: list[int] = []
    n_train_laps_w: list[int] = []

    for t, tr in enumerate(trials):
        same_group = train_group_t == int(train_group_t[t])
        if bool(cfg.train_all_laps):
            train_idx = np.where(same_group)[0]
        else:
            train_idx = np.where(same_group & (np.arange(train_group_t.size) != t))[0]
        if train_idx.size == 0:
            continue

        rates_ux, valid_x = _pooled_tuning_rates(
            counts_utx=counts_utx,
            dwell_tx=dwell_tx,
            train_idx=train_idx,
            smooth_sigma_bins=float(cfg.smooth_sigma_bins),
            rate_floor_hz=float(cfg.rate_floor_hz),
        )
        if not np.any(valid_x):
            continue

        s = max(0, int(tr.start_idx_0b))
        e = min(n_samples, int(tr.stop_idx_0b_exclusive))
        if e <= s:
            continue

        w0 = s
        while (w0 + tau_samples) <= e:
            w1 = w0 + tau_samples
            keep_w = sample_keep[w0:w1]
            valid_fraction = float(np.mean(keep_w)) if keep_w.size else 0.0
            if valid_fraction >= float(cfg.min_valid_window_fraction) and np.any(keep_w):
                actual_pos = float(np.nanmean(x[w0:w1][keep_w]))
                actual_b = int(_bin_indices(np.array([actual_pos], dtype=np.float64), edges)[0])
                if 0 <= actual_b < n_bins:
                    lo = int(np.searchsorted(spk_idx_keep, w0, side="left"))
                    hi = int(np.searchsorted(spk_idx_keep, w1, side="left"))
                    cids_w = spk_cid_keep[lo:hi]
                    counts_u = np.zeros(cids.size, dtype=np.float64)
                    for cid in cids_w:
                        u = cell_to_u.get(int(cid))
                        if u is not None:
                            counts_u[u] += 1.0

                    posterior = _posterior_from_counts(
                        spike_counts_u=counts_u,
                        rates_ux=rates_ux,
                        valid_x=valid_x,
                        tau_s=float(cfg.tau_s),
                    )
                    if posterior is not None:
                        dec_b = int(np.argmax(posterior))
                        dec_x = float(centers[dec_b])
                        posterior_rows.append(posterior)
                        decoded_bin.append(dec_b)
                        decoded_x.append(dec_x)
                        actual_bin.append(actual_b)
                        actual_x.append(actual_pos)
                        error_cm.append(abs(dec_x - actual_pos))
                        prob_actual.append(float(posterior[actual_b]))
                        trial_index.append(int(tr.trial_index))
                        condway_w.append(int(tr.condway))
                        train_group_w.append(int(train_group_t[t]))
                        train_group_label_w.append(str(train_group_label_t[t]))
                        window_start.append(int(w0))
                        window_stop.append(int(w1))
                        n_spikes_w.append(int(np.sum(counts_u)))
                        n_train_laps_w.append(int(train_idx.size))

            w0 += tau_samples

    if posterior_rows:
        posterior_wx = np.vstack(posterior_rows).astype(np.float64, copy=False)
    else:
        posterior_wx = np.zeros((0, n_bins), dtype=np.float64)

    return BayesianDecodingResult(
        cell_ids=cids.copy(),
        idcond_t=idcond_t.copy(),
        xbin_edges=edges,
        xbin_centers=centers,
        posterior_wx=posterior_wx,
        decoded_bin_w=np.asarray(decoded_bin, dtype=np.int64),
        decoded_x_w=np.asarray(decoded_x, dtype=np.float64),
        actual_bin_w=np.asarray(actual_bin, dtype=np.int64),
        actual_x_w=np.asarray(actual_x, dtype=np.float64),
        error_cm_w=np.asarray(error_cm, dtype=np.float64),
        prob_actual_w=np.asarray(prob_actual, dtype=np.float64),
        trial_index_w=np.asarray(trial_index, dtype=np.int64),
        condway_w=np.asarray(condway_w, dtype=np.int64),
        train_group_w=np.asarray(train_group_w, dtype=np.int64),
        train_group_label_w=np.asarray(train_group_label_w, dtype=np.str_),
        window_start_w=np.asarray(window_start, dtype=np.int64),
        window_stop_w=np.asarray(window_stop, dtype=np.int64),
        n_spikes_w=np.asarray(n_spikes_w, dtype=np.int64),
        n_train_laps_w=np.asarray(n_train_laps_w, dtype=np.int64),
    )
