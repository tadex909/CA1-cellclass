from pathlib import Path
from typing import Union, List, Dict, Tuple
from scipy.stats import median_abs_deviation
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import numpy as np
from scipy.io import loadmat
from scipy.signal import convolve
from scipy.signal.windows import gaussian
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np



class TSD:
    """Minimal Time-Series Data container used in your original code."""
    def __init__(self, tvec: np.ndarray, data: np.ndarray):
        self.tvec = np.asarray(tvec)
        self.data = np.asarray(data)
        self.usr = {}
        self.cfg = {'history': {'func_name': [], 'cfg': []}}


def zscore(signal: np.ndarray) -> np.ndarray:
    """Z-score by median / MAD across rows (keeps shape)."""
    signal = np.atleast_2d(signal)
    med = np.median(signal, axis=1)
    mad = median_abs_deviation(signal, axis=1)
    # avoid division by zero
    mad_safe = np.where(mad == 0, 1.0, mad)
    zscored_sig = ((signal.T - med) / mad_safe).T
    return zscored_sig


def gaussian_kernel(win_size_bins: int, sd_bins: float) -> np.ndarray:
    win = int(round(win_size_bins))
    if win < 1:
        win = 1
    gk = gaussian(win, std=sd_bins)
    s = gk.sum()
    return gk / s if s != 0 else np.ones(win) / win


def decode_argmax_posterior(posterior_data: np.ndarray) -> np.ndarray:
    """
    Decode argmax safely: posterior_data shape (n_time, n_states).
    Returns indices (nan when entire row is nan).
    """
    posterior_data = np.asarray(posterior_data)
    if posterior_data.ndim != 2:
        raise ValueError("posterior_data must be 2D (time x states).")
    n_time = posterior_data.shape[0]
    decoded = np.full(n_time, np.nan)
    valid_mask = ~np.all(np.isnan(posterior_data), axis=1)
    if np.any(valid_mask):
        decoded[valid_mask] = np.nanargmax(posterior_data[valid_mask, :], axis=1)
    return decoded


def compute_speed_sliding(tvec: np.ndarray, pos: np.ndarray, window_ms: float = 300) -> np.ndarray:
    """
    Compute speed (absolute derivative) then moving-average over window_ms.
    Handles NaNs in pos by linear interpolation over short gaps.
    """
    tvec = np.asarray(tvec)
    pos = np.asarray(pos)

    if len(tvec) < 2:
        return np.zeros_like(pos)

    # interpolate NaNs in pos (linear) for derivative computation
    pos_interp = pos.copy().astype(float)
    nans = np.isnan(pos_interp)
    if np.any(nans):
        good = ~nans
        if good.sum() < 2:
            # not enough points to interpolate
            pos_interp[nans] = 0.0
        else:
            pos_interp[nans] = np.interp(tvec[nans], tvec[good], pos_interp[good])

    # dt median (s)
    dt = np.median(np.diff(tvec))
    # instantaneous absolute velocity
    vel_inst = np.abs(np.gradient(pos_interp, tvec))
    # window in bins
    win_bins = max(1, int(round((window_ms / 1000.0) / dt)))
    kernel = np.ones(win_bins) / win_bins
    speed = convolve(vel_inst, kernel, mode='same')
    return speed


def tvec_to_edges(tvec: np.ndarray) -> np.ndarray:
    tvec = np.asarray(tvec)
    if len(tvec) < 2:
        dt = 1.0
    else:
        dt = np.median(np.diff(tvec))
    edges = np.concatenate((tvec - dt / 2.0, [tvec[-1] + dt / 2.0]))
    return edges


def compute_tuning_curves(spike_data: List[np.ndarray],
                          beh_tvec: np.ndarray, beh_values: np.ndarray,
                          xmin: float = 0, xmax: float = 145, n_bins: int = 50,
                          smooth_size: int = None, smooth_sd_bins: float = 2.0,
                          speed: np.ndarray = None, speed_tvec: np.ndarray = None,
                          speed_threshold: float = 2.0, speed_window_ms: float = 300,
                          arena_min: float = 0.0, arena_max: float = 145.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute tuning curves: firing rate per position bin.
    Excludes positions outside [arena_min, arena_max] and times with speed <= threshold.
    Returns (tuning_curves: n_cells x n_bins, bin_edges).
    """
    beh_tvec = np.asarray(beh_tvec)
    beh_values = np.asarray(beh_values)
    n_neurons = len(spike_data)
    bin_edges = np.linspace(xmin, xmax, n_bins + 1)

    # compute speed if not provided
    if speed is None:
        speed_used = compute_speed_sliding(beh_tvec, beh_values, window_ms=speed_window_ms)
    else:
        speed = np.asarray(speed)
        if speed_tvec is None or np.array_equal(np.asarray(speed_tvec), beh_tvec):
            speed_used = speed
        else:
            speed_used = np.interp(beh_tvec, np.asarray(speed_tvec), speed)

    # valid times (in-arena & moving)
    valid_mask = (beh_values >= arena_min) & (beh_values <= arena_max) & (speed_used > speed_threshold)

    # digitize position samples into position bins
    beh_digitized = np.clip(np.digitize(beh_values, bin_edges) - 1, 0, n_bins - 1)
    occupancy_counts = np.bincount(beh_digitized[valid_mask], minlength=n_bins).astype(float)

    if smooth_size is not None:
        gk = gaussian_kernel(smooth_size, smooth_sd_bins)
        occupancy_counts = convolve(occupancy_counts, gk, mode='same')

    # occupancy time in seconds: number_of_samples_per_bin * dt
    if len(beh_tvec) >= 2:
        dt = np.median(np.diff(beh_tvec))
    else:
        dt = 1.0
    occupancy_time = occupancy_counts * dt

    tuning_curves = np.full((n_neurons, n_bins), np.nan, dtype=float)
    tmax = beh_tvec[-1] if len(beh_tvec) else np.inf

    for n_idx in range(n_neurons):
        spike_times = np.asarray(spike_data[n_idx])
        if spike_times.size == 0:
            spike_counts = np.zeros(n_bins, dtype=float)
        else:
            # ignore spikes beyond behavior time range
            spike_times = spike_times[spike_times < tmax]
            if spike_times.size == 0:
                spike_counts = np.zeros(n_bins, dtype=float)
            else:
                # map each spike time to the nearest behavior sample index (digitize over time samples)
                time_bins = tvec_to_edges(beh_tvec)
                spike_time_bins = np.digitize(spike_times, time_bins) - 1
                valid_spike_idx = (spike_time_bins >= 0) & (spike_time_bins < len(beh_values))
                if not np.any(valid_spike_idx):
                    spike_counts = np.zeros(n_bins, dtype=float)
                else:
                    # keep only spikes that fall in valid_mask (speed & arena)
                    spk_bins = spike_time_bins[valid_spike_idx]
                    spk_mask = valid_mask[spk_bins]
                    spk_pos_bins = beh_digitized[spk_bins[spk_mask]]
                    spike_counts = np.bincount(spk_pos_bins, minlength=n_bins).astype(float)

        # smooth spike counts if requested
        if smooth_size is not None:
            gk = gaussian_kernel(smooth_size, smooth_sd_bins)
            spike_counts = convolve(spike_counts, gk, mode='same')

        with np.errstate(divide='ignore', invalid='ignore'):
            fr = np.divide(spike_counts, occupancy_time)
            fr[occupancy_time == 0] = np.nan
        tuning_curves[n_idx, :] = fr

    return tuning_curves, bin_edges


def make_q_from_s(l_spk_times: List[np.ndarray], t_vec: np.ndarray,
                  smooth_size: int = None, smooth_sd_bins: float = 2.0,
                  speed: np.ndarray = None, speed_tvec: np.ndarray = None,
                  speed_threshold: float = 2.0) -> TSD:
    """
    Build Q-matrix (n_cells x n_timebins) from spike time lists.
    If speed provided, Q contains spikes only for times where speed > threshold.
    Returns a TSD with tvec centers and data matrix.
    """
    func_name = 'make_q_from_s'
    t_vec = np.asarray(t_vec)
    if all(len(spikes) == 0 for spikes in l_spk_times):
        return TSD(np.array([]), np.array([[]]))

    edges = tvec_to_edges(t_vec)

    # prepare speed on t_vec
    speed_used = None
    if speed is not None:
        speed = np.asarray(speed)
        if speed_tvec is None or np.array_equal(np.asarray(speed_tvec), t_vec):
            speed_used = speed
        else:
            speed_used = np.interp(t_vec, np.asarray(speed_tvec), speed)

    gk = None
    if smooth_size is not None:
        gk = gaussian_kernel(smooth_size, smooth_sd_bins)

    n_cells = len(l_spk_times)
    n_time_bins = len(t_vec)
    q_mat = np.zeros((n_cells, n_time_bins), dtype=float)

    for icell, spk_times in enumerate(l_spk_times):
        if spk_times.size == 0:
            counts = np.zeros(n_time_bins, dtype=float)
        else:
            # Map spikes to Q time bins (histogram using edges)
            if speed_used is None:
                counts, _ = np.histogram(spk_times, bins=edges)
            else:
                spike_bins = np.digitize(spk_times, edges) - 1
                valid = (spike_bins >= 0) & (spike_bins < len(speed_used))
                if not np.any(valid):
                    counts = np.zeros(n_time_bins, dtype=float)
                else:
                    valid = (spike_bins >= 0) & (spike_bins < len(speed_used))

                    # Correction : on force spike_bins à rester entre 0 et len(speed_used)-1
                    spike_bins = np.clip(spike_bins, 0, len(speed_used) - 1)

                    keep_mask = valid & (speed_used[spike_bins] > speed_threshold)

                    if not np.any(keep_mask):
                        counts = np.zeros(n_time_bins, dtype=float)
                    else:
                        filtered_spikes = spk_times[keep_mask]
                        counts, _ = np.histogram(filtered_spikes, bins=edges)
        row = counts.astype(float)
        if gk is not None:
            row = convolve(row, gk, mode='same')
        # ensure length
        if row.shape[0] != n_time_bins:
            row = np.resize(row, (n_time_bins,))
        q_mat[icell, :] = row

    q_tsd = TSD(t_vec, q_mat)
    q_tsd.cfg['history']['func_name'].append(func_name)
    q_tsd.cfg['history']['cfg'].append({'smooth_size': smooth_size,
                                        'smooth_sd_bins': smooth_sd_bins,
                                        'speed_threshold': (speed_threshold if speed is not None else None)})
    return q_tsd


def decode_z(q_tsd: TSD, tuning_curve: np.ndarray, no_spikes_in_bin: str = 'nans',
             exclude_method: str = 'frate', n_min_neurons: int = 1, n_min_spikes: int = 1) -> TSD:
    """
    Bayesian decoding.
    q_tsd.data : (n_cells, n_time)
    tuning_curve : (n_cells, n_position_bins)
    Returns posterior as TSD with tvec (time centers) and data (n_time x n_pos)
    """
    q_data = np.asarray(q_tsd.data)
    tvec = np.asarray(q_tsd.tvec).copy()
    # transpose if q_data shape unexpected (ensure n_cells x n_time)
    if q_data.ndim == 1:
        q_data = q_data.reshape(1, -1)
    n_cells_q, n_time_q = q_data.shape

    if tuning_curve.ndim != 2:
        raise ValueError("tuning_curve must be 2D (n_cells x n_pos).")
    n_cells_tc, n_pos = tuning_curve.shape

    # If tuning curve and q_data have different cell dimension ordering, try to fix.
    if n_cells_q != n_cells_tc:
        # try transpose tuning_curve
        if n_cells_q == n_pos and n_cells_tc == n_time_q:
            tuning_curve = tuning_curve.T
            n_cells_tc, n_pos = tuning_curve.shape
        else:
            raise ValueError(f"Cell count mismatch: q has {n_cells_q} cells, tuning has {n_cells_tc} cells.")

    # ensure time dimension consistent: q_data columns == len(tvec)
    if n_time_q != len(tvec):
        # try convert tvec edges -> centers
        if len(tvec) == n_time_q + 1:
            tvec = (tvec[:-1] + tvec[1:]) / 2.0
        elif n_time_q == len(tvec) + 1:
            q_data = q_data[:, :len(tvec)]
            n_time_q = q_data.shape[1]
        else:
            raise ValueError(f"Time mismatch: q columns {n_time_q}, tvec {len(tvec)}")

    # bin duration
    t_bin_size = np.median(np.diff(tvec)) if len(tvec) > 1 else 1.0

    # uniform prior
    prior = np.ones(n_pos) / n_pos

    p = np.zeros((len(tvec), n_pos), dtype=float)
    eps = 1e-12

    # Precompute terms that depend on position only
    tc_safe = np.maximum(tuning_curve, eps)  # (n_cells, n_pos)
    sum_tc = np.sum(tc_safe, axis=0)  # sum over cells for each position (n_pos,)

    # For numerical stability compute in log-space where possible
    # q_data is counts s_i(t) for each cell i and time t
    # For each position x: log P(s|x) = sum_i s_i * log(f_i(x)) - T * sum_i f_i(x)
    log_tc = np.log(tc_safe)  # shape (n_cells, n_pos)

    # iterate time bins and compute posterior over positions
    for tt in range(len(tvec)):
        s_t = q_data[:, tt]  # (n_cells,)
        # compute sum_i s_i * log(f_i(x)) for each x: -> dot(s_t, log_tc) over cells
        # result shape (n_pos,)
        term1 = np.dot(s_t, log_tc)  # (n_pos,)
        term2 = -sum_tc * t_bin_size  # (n_pos,)
        logp_x = term1 + term2 + np.log(prior)
        # exponentiate safely (subtract max for stability)
        mx = np.nanmax(logp_x)
        if not np.isfinite(mx):
            p[tt, :] = np.nan
            continue
        rel = np.exp(logp_x - mx)
        denom = np.nansum(rel)
        if denom == 0:
            p[tt, :] = np.nan
        else:
            p[tt, :] = rel / denom

    # Exclusion rules: mark rows (time bins) with too few spikes as NaN or zeros
    if exclude_method == 'nNeurons':
        n_active_neurons = np.sum(q_data >= n_min_spikes, axis=0)
        toss_idx = n_active_neurons < n_min_neurons
    elif exclude_method == 'frate':
        total_spikes = np.sum(q_data, axis=0)
        toss_idx = total_spikes < n_min_spikes
    else:
        raise ValueError("Invalid exclude_method: choose 'nNeurons' or 'frate'.")

    if no_spikes_in_bin == 'zeros':
        p[toss_idx, :] = 0.0
    elif no_spikes_in_bin == 'nans':
        p[toss_idx, :] = np.nan
    else:
        raise ValueError("no_spikes_in_bin must be 'zeros' or 'nans'.")

    p_tsd = TSD(tvec, p)
    # store some diagnostic info
    p_tsd.usr['n_active_neurons'] = (np.sum(q_data >= n_min_spikes, axis=0) if q_data.size else np.array([]))
    p_tsd.cfg['history']['func_name'].append('decode_z')
    return p_tsd

def get_paths(base_path: Union[str, Path], session_name: str) -> Dict[str, Path]:
    base_path = Path(base_path)

    # dossier racine des ratemaps
    rmap_root = Path(r'//Epsztein-nas02//TEAM//Vinca//MATLAB//Data//Data_thèse_Vinca')

    phenosys_path = list(base_path.glob('*Phenosys.mat'))[0]
    traj_path = list(base_path.glob('*TrajData.mat'))[0]
    ephys_path = list(base_path.glob('*ePhy.mat'))[0]

    # --- ratemap spécifique à la session ---
    ratemap_path = rmap_root / f"{session_name}_Ratemap_final_thèse_CASE1.mat"

    if not ratemap_path.exists():
        raise FileNotFoundError(f"Ratemap introuvable : {ratemap_path}")

    return {
        'phenosys': phenosys_path,
        'traj': traj_path,
        'ephys': ephys_path,
        'ratemap': ratemap_path
    }

def load_cell_types(ratemap_path: Union[str, Path]) -> Dict[int, int]:
    """
    Retourne un dict {cell_id: cell_type}
    """
    rm = loadmat(ratemap_path)['allcel']

    cell_ids = np.squeeze(rm['id_cel'][0, 0]).astype(int)
    cell_types = np.squeeze(rm['type_u'][0, 0]).astype(int)
    print(cell_types)
    return dict(zip(cell_ids, cell_types))

def load_trajectory_blocks(traj_path: Union[str, Path]) -> Dict:
    traj_raw = loadmat(traj_path)
    traj_arr = traj_raw['Traj'].squeeze()
    block_ix = np.squeeze(np.hstack(traj_arr['Cond']))
    start_ix = np.squeeze(np.hstack(traj_arr['start']))
    stop_ix = np.squeeze(np.hstack(traj_arr['stop']))
    way = np.squeeze(np.hstack(traj_arr['icondway_tr']))
    forward = (way % 2 == 0)
    blocks = {}
    for c_bl_ix in np.unique(block_ix):
        blocks[c_bl_ix] = {}
        block_mask = block_ix == c_bl_ix
        for is_fwd in (0, 1):
            direction_mask = forward == is_fwd
            full_mask = block_mask & direction_mask
            blocks[c_bl_ix][is_fwd] = np.vstack((start_ix[full_mask], stop_ix[full_mask]))
    return blocks


def load_trajectory(phenosys_path: Union[str, Path], blocks: Dict) -> Dict[int, List[TSD]]:
    ph_raw = loadmat(phenosys_path)
    full_time = np.squeeze(ph_raw['t_ds'])
    full_time = full_time - full_time[0]
    x_raw = np.squeeze(ph_raw['X_ds_n'])
    print("Max de x_raw :", np.nanmax(x_raw))
    all_traj_dict = {}
    for block_ix, edges in blocks.items():
        for start_ix, stop_ix in zip(edges[0], edges[1]):
            c_raw_x = x_raw[start_ix:stop_ix]
            c_time = full_time[start_ix:stop_ix]
            prev_traj = all_traj_dict.setdefault(block_ix, [])
            prev_traj.append(TSD(c_time, c_raw_x))
            all_traj_dict[block_ix] = prev_traj
    return all_traj_dict


def load_spikes(ephys_path: Union[str, Path]) -> Tuple[Dict[int, np.ndarray], float]:
    ephys_raw = loadmat(ephys_path)['allcel']
    spike_times = np.squeeze(np.hstack(ephys_raw['time_spk'][0]))
    spike_ids = np.squeeze(ephys_raw['id_spk'][0, 0])
    cell_ids = np.squeeze(ephys_raw['id_cel'][0, 0])
    spikes = {}
    for c_cell in cell_ids:
        spikes[int(c_cell)] = spike_times[spike_ids == c_cell]
    return spikes, spike_times.max() if spike_times.size else 0.0


def decode_error(decoded_position: np.ndarray, real_position: np.ndarray,
                 n_bins: int = 100, xmin: float = 0, xmax: float = 145,
                 show_cm: bool = True):
    bin_edges = np.linspace(xmin, xmax, n_bins + 1)
    keep = ~np.isnan(decoded_position)
    decoded_keep = decoded_position[keep]
    real_keep = real_position[keep]
    real_digitized = np.clip(np.digitize(real_keep, bin_edges) - 1, 0, n_bins - 1)
    decoded_digitized = np.clip(np.digitize(decoded_keep, bin_edges) - 1, 0, n_bins - 1)
    err = np.abs(real_keep - decoded_keep)
    conf_mat = confusion_matrix(real_digitized, decoded_digitized, labels=np.arange(n_bins))
    if show_cm:
        plt.figure()
        ConfusionMatrixDisplay(conf_mat, display_labels=np.arange(n_bins)).plot()
    return keep, err, conf_mat


def run_decoding(data_dir: Union[str, Path], direction: int = 0, condition_ix: int = 1,
                 bin_dur_s: float = 0.15, plot_tuning_curves: bool = True,
                 plot_decoded_traj: bool = True, smooth_size: int = 11,
                 smooth_sd_bins: float = 2.0) -> Tuple[Dict[int, List[TSD]], np.ndarray, TSD, np.ndarray]:
    """
    High-level pipeline:
      - load paths, trajectories, spikes
      - compute tuning curves (averaged across trials)
      - compute Q from spikes on a regular time grid t_vec
      - decode posterior and return decoded positions (mapped to bin centers)
    """
    data_dir = Path(data_dir)
    all_paths = get_paths(data_dir)

    # blocks forward/backward
    bl_ix_both_dir = load_trajectory_blocks(all_paths['traj'])
    bl_ix_one_dir = {ix: c_traj[direction] for ix, c_traj in bl_ix_both_dir.items()}

    trajs_one_dir = load_trajectory(all_paths['phenosys'], bl_ix_one_dir)
    spikes, _ = load_spikes(all_paths['ephys'])

    # --- CHARGEMENT TYPES CELLULAIRES ---
    cell_types = load_cell_types(all_paths['ratemap'])
    PYRAMIDAL_CODE = 1
    print(cell_types)
    # --- FILTRAGE PYRAMIDALES ---
    spikes_pyr = {
        cid: spikes[cid]
        for cid in spikes
        if cid in cell_types and cell_types[cid] == PYRAMIDAL_CODE
    }

    spikes_list = list(spikes_pyr.values())

    print(f"{len(spikes)} cellules totales")
    print(f"{len(spikes_list)} cellules pyramidales")
    plt.figure(figsize=(12, 6))
    for i, spikes in enumerate(spikes_list):
        if len(spikes) > 0:
            plt.vlines(spikes[spikes < 1000], i, i+1, colors='blue', alpha=0.5)  # Focus sur 0-1000 s
    plt.ylim(0, len(spikes_list))
    plt.xlabel('Time (s)')
    plt.title('Raster plot des spikes (0-1000 s)')
    plt.show()
    if len(spike_list) < 10:
        print("Session ignorée : moins de 10 cellules pyramidales")
        raise SystemExit
    # pick first trajectory to get bin edges
    c_traj = trajs_one_dir[condition_ix][0]
    _, beh_edges = compute_tuning_curves(spikes_list, c_traj.tvec, c_traj.data,
                                         n_bins=100, smooth_size=smooth_size,
                                         smooth_sd_bins=smooth_sd_bins)
    # compute tuning curves per trial and average
    all_tcs = np.dstack([compute_tuning_curves(spikes_list, tr.tvec, tr.data,
                                               n_bins=100, smooth_size=smooth_size,
                                               smooth_sd_bins=smooth_sd_bins)[0]
                         for tr in trajs_one_dir[condition_ix]])
    avg_tuning_curves = np.nanmean(all_tcs, axis=2)  # shape n_cells x n_pos

    if plot_tuning_curves:
        fig, ax = plt.subplots()
        # compute centers for display
        centers = (beh_edges[:-1] + beh_edges[1:]) / 2.0
        # normalize per neuron for display
        min_tc = np.nanmin(avg_tuning_curves[:, 10:-10], axis=1)
        max_tc = np.nanmax(avg_tuning_curves[:, 10:-10], axis=1)
        denom = (max_tc - min_tc)
        denom[denom == 0] = 1.0
        normed_tc = ((avg_tuning_curves.T - min_tc) / denom).T
        normed_tc[:, :10] = np.nan
        normed_tc[:, -10:] = np.nan
        ax.imshow(normed_tc, aspect='auto', interpolation='none',
                  extent=(centers[0], centers[-1], 0, avg_tuning_curves.shape[0]),
                  origin='lower')
        ax.set_xlabel('Position (units)')
        ax.set_ylabel('Neuron index')

    # build global time vector covering all trials of the condition
    t_min = min([np.min(tr.tvec) for tr in trajs_one_dir[condition_ix]])
    t_max = max([np.max(tr.tvec) for tr in trajs_one_dir[condition_ix]])
    t_vec = np.arange(t_min, t_max, bin_dur_s)

    # build continuous position trace by concatenating and sorting
    all_times = np.hstack([tr.tvec for tr in trajs_one_dir[condition_ix]])
    all_pos = np.hstack([tr.data for tr in trajs_one_dir[condition_ix]])
    sort_idx = np.argsort(all_times)
    all_times = all_times[sort_idx]
    all_pos = all_pos[sort_idx]

    # interpolate positions on t_vec (NaN outside)
    pos_on_tvec = np.interp(t_vec, all_times, all_pos, left=np.nan, right=np.nan)

    # apply spatial mask
    valid_arena = (pos_on_tvec >= 0.0) & (pos_on_tvec <= 145.0)
    pos_on_tvec[~valid_arena] = np.nan

    # compute speed on t_vec
    speed_on_tvec = compute_speed_sliding(t_vec, pos_on_tvec, window_ms=300)

    # Q-matrix (TSD)
    q_tsd = make_q_from_s(spikes_list, t_vec=t_vec, smooth_size=10,
                          speed=speed_on_tvec, speed_tvec=t_vec,
                          speed_threshold=2.0)

    # decode posterior using averaged tuning curves
    posterior = decode_z(q_tsd, avg_tuning_curves, no_spikes_in_bin='nans',
                         exclude_method='frate', n_min_neurons=1, n_min_spikes=1)

    # argmax and map to bin centers
    pos_idx = decode_argmax_posterior(posterior.data)  # indices
    centers = (beh_edges[:-1] + beh_edges[1:]) / 2.0
    decoded_position = np.full_like(pos_idx, np.nan, dtype=float)
    valid = ~np.isnan(pos_idx)
    decoded_position[valid] = centers[pos_idx[valid].astype(np.intp)]

    if plot_decoded_traj:
        fig, ax = plt.subplots()
        for tr in trajs_one_dir[condition_ix]:
            ax.plot(tr.tvec, tr.data, c='C0', alpha=0.6)
        ax.plot(posterior.tvec, decoded_position, c='.2', lw=1)
        fig.suptitle(f'Decoding for {all_paths["ephys"].stem}')
        fig.tight_layout()

    return trajs_one_dir, avg_tuning_curves, posterior, decoded_position

if __name__ == "__main__":

    # --- PARAMÈTRES ---
    base_path = Path(r"//Epsztein-nas02//TEAM//Vinca//MATLAB//Data//Data_thèse_Vinca")
    condition_ix = 1
    bin_dur_s = 0.1
    n_bins = 100
    smooth_size = 11
    smooth_sd_bins = 2
    
    session_name = input("Entrez le nom de l'animal/session : ").strip()
    # --- CHARGEMENT DONNÉES ---
    paths = get_paths(base_path,session_name)
    blocks = load_trajectory_blocks(paths["traj"])
   
    spikes, _ = load_spikes(paths["ephys"])

    # --- TYPES CELLULAIRES ---
    cell_types = load_cell_types(paths["ratemap"])
    PYRAMIDAL_CODE = 1

    # --- FILTRAGE PYRAMIDALES ---
    spikes_pyr = {
        cid: spikes[cid]
        for cid in spikes
   #     if cid in cell_types and cell_types[cid] == PYRAMIDAL_CODE
        if cid in cell_types and cell_types[cid] == 1
    }

    spike_list = list(spikes_pyr.values())

    print(f"{len(spikes)} cellules totales")
    print(f"{len(spike_list)} cellules pyramidales")
    if len(spike_list) < 10:
        print("Session ignorée : moins de 10 cellules pyramidales")
        raise SystemExit
    # --- FIGURE COMPARATIVE ---
    fig, ax = plt.subplots(figsize=(12, 6))
 
    colors = {0: "dodgerblue", 1: "seagreen"}
    labels = {0: "Forward (real)", 1: "Backward (real)"}

    decoded_position = None
    posterior = None
    t_vec = None

    for direction in (0, 1):  # 0 = forward, 1 = backward
        blocks_one_dir = {ix: c_traj[direction] for ix, c_traj in blocks.items()}
        trajs = load_trajectory(paths["phenosys"], blocks_one_dir)

        # --- TUNING CURVES ---
        tuning_curves_all = np.dstack([
            compute_tuning_curves(spike_list, traj.tvec, traj.data,
                                  n_bins=n_bins,
                                  smooth_size=smooth_size,
                                  smooth_sd_bins=smooth_sd_bins)[0]
            for traj in trajs[condition_ix]
        ])
        avg_tuning_curves = np.nanmean(tuning_curves_all, axis=2)

        # --- TEMPS GLOBAL ---
        t_min = min([min(traj.tvec) for traj in trajs[condition_ix]])
        t_max = max([max(traj.tvec) for traj in trajs[condition_ix]])
        t_vec = np.arange(t_min, t_max, bin_dur_s)

        # Position continue
        all_times = np.hstack([traj.tvec for traj in trajs[condition_ix]])
        all_pos = np.hstack([traj.data for traj in trajs[condition_ix]])
        sort_idx = np.argsort(all_times)
        all_times, all_pos = all_times[sort_idx], all_pos[sort_idx]
        pos_on_tvec = np.interp(t_vec, all_times, all_pos, left=np.nan, right=np.nan)

        valid_arena = (pos_on_tvec >= 0.0) & (pos_on_tvec <= 145.0)
        pos_on_tvec[~valid_arena] = np.nan
        speed_on_tvec = compute_speed_sliding(t_vec, pos_on_tvec, window_ms=300)

        # --- Q-MATRIX & DECODING ---
        q_matrix = make_q_from_s(spike_list, t_vec=t_vec,
                                 smooth_size=10, speed=speed_on_tvec,
                                 speed_tvec=t_vec, speed_threshold=2.0)

        posterior = decode_z(q_matrix, avg_tuning_curves,
                             no_spikes_in_bin="nans",
                             exclude_method="frate",
                             n_min_neurons=1, n_min_spikes=1)


        decoded_bins = decode_argmax_posterior(posterior.data)
        decoded_position = np.full_like(decoded_bins, np.nan, dtype=float)
        valid_mask = ~np.isnan(decoded_bins)
        bin_edges = np.linspace(0, 145, n_bins + 1)
        decoded_position[valid_mask] = bin_edges[decoded_bins[valid_mask].astype(int)]

        # --- PLOT TRAJECTOIRES ---
        for traj in trajs[condition_ix]:
            ax.plot(traj.tvec, traj.data, color=colors[direction], alpha=0.3,
                    label=labels[direction] if traj is trajs[condition_ix][0] else "")

    plt.figure()
    plt.plot(t_vec, speed_on_tvec, label='Speed')
    plt.axhline(2.0, color='red', linestyle='--', label='Speed threshold')
    plt.xlabel('Time (s)')
    plt.ylabel('Speed (cm/s)')
    plt.legend()
    plt.show()

    # --- PLOT DECODED ---
    ax.plot(posterior.tvec, decoded_position, color="red", lw=1.5, label="Decoded position")
    print(decoded_position)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Position")
    ax.set_title("Forward vs Backward trajectories with decoded position")
    ax.legend()
    #plt.show()

    # --- créer une colormap type elife ---
    colors = [(1, 1, 1),(0.8, 0.3, 0)]  # orange clair -> blanc
    cmap_elife = mcolors.LinearSegmentedColormap.from_list("elife_orange", colors)

    n_bins = posterior.data.shape[1]
    bin_edges = np.linspace(0, 145, n_bins + 1)
    positions_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    import scipy.ndimage

    # --- Calcul de l'argmax par position ---
    posterior_data = posterior.data.copy()

    # normalisation (déjà fait, mais on peut recalculer juste pour être sûr)
    posterior_norm = posterior_data / np.nansum(posterior_data, axis=1, keepdims=True)
    posterior_norm[np.isnan(posterior_norm)] = 0

    # argmax (indice de la position avec proba max à chaque instant)
    argmax_idx = np.nanargmax(posterior_norm, axis=1)

    # smoothing très large pour avoir une ligne “globale”
    # ici sigma = 20 correspond à un lissage sur ~20 bins temporels
    smoothed_idx = scipy.ndimage.gaussian_filter1d(argmax_idx.astype(float), sigma=10)

    # conversion en positions réelles
    positions_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    smoothed_positions = positions_centers[smoothed_idx.astype(int)]

    # --- Superposition sur le heatmap existant ---
    fig, ax = plt.subplots(figsize=(12,6))

    im = ax.imshow(posterior_norm.T, origin='lower', aspect='auto',
                extent=[posterior.tvec[0], posterior.tvec[-1],
                        positions_centers[0], positions_centers[-1]],
                cmap=cmap_elife, interpolation='nearest')
    fig.colorbar(im, ax=ax, label='Posterior probability')

    # superposer toutes les trajectoires réelles
    colors_traj = {0: "dodgerblue", 1: "seagreen"}
    for direction in (0,1):
        blocks_one_dir = {ix: c_traj[direction] for ix, c_traj in blocks.items()}
        trajs = load_trajectory(paths["phenosys"], blocks_one_dir)
        for traj in trajs[condition_ix]:
            ax.plot(traj.tvec, traj.data, color=colors_traj[direction], lw=1.5, alpha=0.4, linestyle='--')

    # --- ligne pointillée smoothed ---
    #ax.plot(posterior.tvec, smoothed_positions, 'r--', lw=2, label='Decoded smoothed')

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Position")
    ax.set_title("Bayesian Decoding Posterior Heatmap with Smoothed Decoded Trajectory")
    ax.legend()
    plt.tight_layout()
    #plt.show()
    import numpy as np
    import matplotlib.pyplot as plt

    # Posterior original
    posterior_data = posterior.data.copy()

    # --- Normalisation par colonne (chaque temps) ---
    posterior_norm = posterior_data / np.nanmax(posterior_data, axis=1, keepdims=True)
    posterior_norm[np.isnan(posterior_norm)] = 0  # remplacer les NaN par 0

    # --- Argmax smoothed (pour ligne décodée) ---
    import scipy.ndimage

    argmax_idx = np.nanargmax(posterior_norm, axis=1)
    smoothed_idx = scipy.ndimage.gaussian_filter1d(argmax_idx.astype(float), sigma=2)
    smoothed_positions = positions_centers[smoothed_idx.astype(int)]

    # --- Plot heatmap normalisée ---
    fig, ax = plt.subplots(figsize=(12,6))

    im = ax.imshow(posterior_norm.T, origin='lower', aspect='auto',
                extent=[posterior.tvec[0], posterior.tvec[-1],
                        positions_centers[0], positions_centers[-1]],
                cmap=cmap_elife, interpolation='nearest')
    fig.colorbar(im, ax=ax, label='Normalized posterior probability')

    # --- Superposer trajectoires réelles en pointillés ---
    colors_traj = {0: "dodgerblue", 1: "seagreen"}
    for direction in (0,1):
        blocks_one_dir = {ix: c_traj[direction] for ix, c_traj in blocks.items()}
        trajs = load_trajectory(paths["phenosys"], blocks_one_dir)
        for traj in trajs[condition_ix]:
            ax.plot(traj.tvec, traj.data, color=colors_traj[direction],
                    lw=1.5, alpha=0.4, linestyle='--')

    # --- Ligne décodée smoothed ---
    #ax.plot(posterior.tvec, smoothed_positions, 'r--', lw=2, label='Decoded smoothed')
    # --- Ligne décodée smoothed remplacée par polyfit ---
    # Ignorer les NaN pour le fit
    from scipy.signal import savgol_filter

    # argmax des positions décodées
    argmax_idx = decode_argmax_posterior(posterior.data)
    valid_mask = ~np.isnan(argmax_idx)
    argmax_idx_clean = argmax_idx.copy()
    argmax_idx_clean[~valid_mask] = 0  # remplacer NaN par 0 pour le filtre

    # lissage Savitzky-Golay : fenêtre impair, polyorder 3
    smoothed_idx = savgol_filter(argmax_idx_clean.astype(float), window_length=31, polyorder=3)

    # clip pour éviter l'IndexError
    smoothed_idx = np.clip(smoothed_idx, 0, len(positions_centers)-1)

    # conversion en positions réelles
    smoothed_positions = positions_centers[smoothed_idx.astype(int)]

    # plot
    ax.plot(posterior.tvec, smoothed_positions, 'gray', lw=2, linestyle='--', label='Decoded smoothed')


    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Position")
    ax.set_title("Bayesian Decoding Posterior Heatmap (Normalized by time bin)")
    ax.legend()
    plt.tight_layout()
    #plt.show()
    print('hello')

    import numpy as np
    import matplotlib.pyplot as plt

    xmin = 0
    xmax = 145
    bin_size = 1.45  # cm
    bin_edges = np.arange(xmin, xmax + bin_size, bin_size)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    colors = {0: "dodgerblue", 1: "seagreen"}
    labels = {0: "Forward", 1: "Backward"}

    import matplotlib.pyplot as plt
    import numpy as np

    fig, axes = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

    for idx, direction in enumerate((0, 1)):
        condition = condition_ix
        blocks_one_dir = {ix: c_traj[direction] for ix, c_traj in blocks.items()}
        trajs_dir = load_trajectory(paths["phenosys"], blocks_one_dir)[condition]
        
        all_errors = []

        for traj in trajs_dir:
            # Interpoler la position décodée sur le temps de la trajectoire
            decoded_on_traj = np.interp(traj.tvec, posterior.tvec, decoded_position, left=np.nan, right=np.nan)
            # Erreur absolue
            err = np.abs(traj.data - decoded_on_traj)
            all_errors.append({'pos': traj.data, 'err': err})
        
        # Moyenne par bin
        bin_means = []
        for i in range(len(bin_edges)-1):
            bin_mask = [(e['pos'] >= bin_edges[i]) & (e['pos'] < bin_edges[i+1]) for e in all_errors]
            errors_in_bin = np.hstack([e['err'][mask] for e, mask in zip(all_errors, bin_mask)])
            bin_means.append(np.nanmean(errors_in_bin) if len(errors_in_bin) > 0 else np.nan)
        
        # Plot sur le subplot correspondant
        axes[idx].plot(bin_centers, bin_means, '-o', color=colors[direction])
        axes[idx].set_ylabel("Erreur abs (cm)")
        axes[idx].set_title(f"Erreur moyenne - {labels[direction]}")
        axes[idx].grid(True)

    axes[-1].set_xlabel("Position (cm)")
    plt.tight_layout()
    plt.show()

        #plt.show()

    import numpy as np
    import pickle
    from pathlib import Path

    # --- Demander le nom de la session/animal ---

    # --- Créer un dossier pour sauvegarder les données ---
    save_dir = Path("saved_sessions")
    save_dir.mkdir(exist_ok=True)

    # --- Exemple de données à sauvegarder (adapter selon ton code) ---
    # x_raw : positions brutes
    # all_traj_dict : dict bloc -> list of TSD
    # spikes : dict cellule -> spike_times
    # decoded_position : np.ndarray
    # posterior : TSD

    data_to_save = {
        "session_name": session_name,
        
        "all_traj_dict": trajs,
        "spikes": spikes,
        "decoded_position": decoded_position,
        "posterior_data": posterior.data if posterior is not None else None,
        "posterior_tvec": posterior.tvec if posterior is not None else None
    }

    # --- Sauvegarde au format .npz (numpy compressed) ---
    npz_path = save_dir / f"{session_name}_data.npz"
    np.savez_compressed(npz_path, **data_to_save)
    print(f"Données sauvegardées en .npz : {npz_path}")

    # --- Optionnel : sauvegarde au format pickle (plus flexible pour objets complexes) ---
    pkl_path = save_dir / f"{session_name}_data.pkl"
    with open(pkl_path, "wb") as f:
        pickle.dump(data_to_save, f)
    print(f"Données sauvegardées en .pkl : {pkl_path}")
