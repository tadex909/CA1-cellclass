from __future__ import annotations

import json
import runpy
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from placefields import (
    BayesianDecoderConfig,
    BayesianDecodingResult,
    TrialInfo,
    average_soft_confusions_by_condition,
    bayesian_decoding_accuracy,
    bayesian_decoding_result_to_frame,
    decoding_error_cm,
    decode_bayesian_position_from_trials,
    hard_decoding_accuracy,
    hard_decoding_correct,
    local_decoding_probability,
    local_decoding_probability_by_actual_bin,
    posterior_entropy_bits,
    posterior_expected_abs_error_cm,
    posterior_mean_x,
    posterior_peak_probability,
    posterior_std_cm,
    position_balanced_bayesian_decoding_accuracy,
    soft_confusion_matrix_for_windows,
    summarize_decoding_by_condition,
    summarize_decoding_by_trial,
    summarize_decoding_measures,
    summarize_soft_confusions,
)


def _trial(trial_index: int, start: int, stop: int) -> TrialInfo:
    return TrialInfo(
        trial_index=trial_index,
        cond=1,
        wb="W",
        condway=1,
        start_idx_0b=start,
        stop_idx_0b_exclusive=stop,
    )


def _cfg(
    *,
    min_speed: float | None = 2.0,
    decode_groupby: str = "condway",
    group_condition_families: bool = False,
    train_all_laps: bool = False,
) -> BayesianDecoderConfig:
    return BayesianDecoderConfig(
        freq_hz=10.0,
        tau_s=1.0,
        bin_size_cm=1.0,
        min_speed=min_speed,
        smooth_sigma_bins=0.0,
        xbin_rem=0,
        min_valid_window_fraction=0.5,
        rate_floor_hz=1e-12,
        decode_groupby=decode_groupby,
        group_condition_families=group_condition_families,
        train_all_laps=train_all_laps,
    )


def _manual_decoding_result() -> BayesianDecodingResult:
    posterior = np.array(
        [
            [0.7, 0.2, 0.1],
            [0.1, 0.6, 0.3],
            [0.2, 0.2, 0.6],
        ],
        dtype=np.float64,
    )
    actual_bin = np.array([0, 2, 1], dtype=np.int64)
    decoded_bin = np.array([0, 1, 2], dtype=np.int64)
    actual_x = np.array([0.0, 20.0, 10.0], dtype=np.float64)
    decoded_x = np.array([0.0, 10.0, 20.0], dtype=np.float64)
    return BayesianDecodingResult(
        cell_ids=np.array([10, 11], dtype=np.int64),
        idcond_t=np.array([1, 1, 1], dtype=np.int64),
        xbin_edges=np.array([-5.0, 5.0, 15.0, 25.0], dtype=np.float64),
        xbin_centers=np.array([0.0, 10.0, 20.0], dtype=np.float64),
        posterior_wx=posterior,
        decoded_bin_w=decoded_bin,
        decoded_x_w=decoded_x,
        actual_bin_w=actual_bin,
        actual_x_w=actual_x,
        error_cm_w=np.abs(decoded_x - actual_x),
        prob_actual_w=posterior[np.arange(posterior.shape[0]), actual_bin],
        trial_index_w=np.array([0, 0, 1], dtype=np.int64),
        condway_w=np.array([1, 1, 1], dtype=np.int64),
        train_group_w=np.array([1, 1, 1], dtype=np.int64),
        train_group_label_w=np.array(["PO W", "PO W", "PO W"], dtype=np.str_),
        window_start_w=np.array([0, 10, 20], dtype=np.int64),
        window_stop_w=np.array([10, 20, 30], dtype=np.int64),
        n_spikes_w=np.array([2, 3, 4], dtype=np.int64),
        n_train_laps_w=np.array([2, 2, 2], dtype=np.int64),
    )


def _position_imbalanced_decoding_result() -> BayesianDecodingResult:
    posterior = np.array(
        [
            [0.8, 0.1, 0.1],
            [0.6, 0.3, 0.1],
            [0.7, 0.2, 0.1],
            [0.3, 0.4, 0.3],
            [0.5, 0.3, 0.2],
        ],
        dtype=np.float64,
    )
    actual_bin = np.array([0, 0, 0, 1, 2], dtype=np.int64)
    decoded_bin = np.argmax(posterior, axis=1).astype(np.int64)
    xbin_centers = np.array([0.0, 10.0, 20.0], dtype=np.float64)
    actual_x = xbin_centers[actual_bin]
    decoded_x = xbin_centers[decoded_bin]
    return BayesianDecodingResult(
        cell_ids=np.array([10, 11], dtype=np.int64),
        idcond_t=np.array([1, 1, 1], dtype=np.int64),
        xbin_edges=np.array([-5.0, 5.0, 15.0, 25.0], dtype=np.float64),
        xbin_centers=xbin_centers,
        posterior_wx=posterior,
        decoded_bin_w=decoded_bin,
        decoded_x_w=decoded_x,
        actual_bin_w=actual_bin,
        actual_x_w=actual_x,
        error_cm_w=np.abs(decoded_x - actual_x),
        prob_actual_w=posterior[np.arange(posterior.shape[0]), actual_bin],
        trial_index_w=np.arange(posterior.shape[0], dtype=np.int64),
        condway_w=np.ones(posterior.shape[0], dtype=np.int64),
        train_group_w=np.ones(posterior.shape[0], dtype=np.int64),
        train_group_label_w=np.full(posterior.shape[0], "PO W", dtype=np.str_),
        window_start_w=np.arange(0, posterior.shape[0] * 10, 10, dtype=np.int64),
        window_stop_w=np.arange(10, (posterior.shape[0] + 1) * 10, 10, dtype=np.int64),
        n_spikes_w=np.ones(posterior.shape[0], dtype=np.int64),
        n_train_laps_w=np.full(posterior.shape[0], 2, dtype=np.int64),
    )


class BayesianDecodingTest(unittest.TestCase):
    def test_decoder_measure_functions_on_known_posteriors(self) -> None:
        result = _manual_decoding_result()

        np.testing.assert_allclose(decoding_error_cm(result), np.array([0.0, 10.0, 10.0]))
        np.testing.assert_array_equal(
            hard_decoding_correct(result),
            np.array([True, False, False]),
        )
        self.assertAlmostEqual(hard_decoding_accuracy(result), 1.0 / 3.0)
        self.assertAlmostEqual(hard_decoding_accuracy(result, tolerance_bins=1), 1.0)
        self.assertAlmostEqual(hard_decoding_accuracy(result, tolerance_cm=10.0), 1.0)

        np.testing.assert_allclose(
            local_decoding_probability(result),
            np.array([0.7, 0.3, 0.2]),
        )
        np.testing.assert_allclose(
            local_decoding_probability(result, radius_bins=1),
            np.array([0.9, 0.9, 1.0]),
        )
        self.assertAlmostEqual(bayesian_decoding_accuracy(result), 0.4)

        np.testing.assert_allclose(posterior_peak_probability(result), np.array([0.7, 0.6, 0.6]))
        expected_entropy = -np.sum(
            result.posterior_wx * np.log2(result.posterior_wx),
            axis=1,
        )
        np.testing.assert_allclose(posterior_entropy_bits(result), expected_entropy)
        np.testing.assert_allclose(
            posterior_entropy_bits(result, normalized=True),
            expected_entropy / np.log2(3.0),
        )

        np.testing.assert_allclose(posterior_mean_x(result), np.array([4.0, 12.0, 14.0]))
        np.testing.assert_allclose(posterior_std_cm(result), np.array([np.sqrt(44.0), 6.0, 8.0]))
        np.testing.assert_allclose(
            posterior_expected_abs_error_cm(result),
            np.array([4.0, 8.0, 8.0]),
        )

        summary = summarize_decoding_measures(result)
        self.assertEqual(summary["n_windows"], 3)
        self.assertAlmostEqual(summary["median_error_cm"], 10.0)
        self.assertAlmostEqual(summary["mean_error_cm"], 20.0 / 3.0)
        self.assertAlmostEqual(summary["hard_accuracy"], 1.0 / 3.0)
        self.assertAlmostEqual(summary["bayesian_accuracy"], 0.4)
        self.assertAlmostEqual(summary["mean_peak_probability"], 1.9 / 3.0)
        self.assertAlmostEqual(summary["mean_posterior_expected_abs_error_cm"], 20.0 / 3.0)

    def test_decoder_measure_selection_and_edge_cases(self) -> None:
        result = _manual_decoding_result()

        np.testing.assert_allclose(
            decoding_error_cm(result, window_idx=np.array([True, False, True])),
            np.array([0.0, 10.0]),
        )
        self.assertAlmostEqual(bayesian_decoding_accuracy(result, window_idx=[0, 2]), 0.45)

        empty_idx = np.array([], dtype=np.int64)
        self.assertEqual(decoding_error_cm(result, window_idx=empty_idx).size, 0)
        self.assertTrue(np.isnan(hard_decoding_accuracy(result, window_idx=empty_idx)))
        self.assertTrue(np.isnan(bayesian_decoding_accuracy(result, window_idx=empty_idx)))
        empty_summary = summarize_decoding_measures(result, window_idx=empty_idx)
        self.assertEqual(empty_summary["n_windows"], 0)
        self.assertTrue(np.isnan(empty_summary["mean_error_cm"]))
        self.assertTrue(np.isnan(empty_summary["bayesian_accuracy"]))

        with self.assertRaises(ValueError):
            hard_decoding_correct(result, tolerance_bins=-1)
        with self.assertRaises(ValueError):
            local_decoding_probability(result, radius_bins=-1)
        with self.assertRaises(ValueError):
            local_decoding_probability(result, radius_bins=1, radius_cm=10.0)

    def test_position_balanced_bayesian_accuracy_balances_actual_bins(self) -> None:
        result = _position_imbalanced_decoding_result()

        self.assertAlmostEqual(bayesian_decoding_accuracy(result), 0.54)
        self.assertAlmostEqual(
            position_balanced_bayesian_decoding_accuracy(result),
            (0.7 + 0.4 + 0.2) / 3.0,
        )

        x, probability_x, n_windows_x = local_decoding_probability_by_actual_bin(result)
        np.testing.assert_allclose(x, np.array([0.0, 10.0, 20.0]))
        np.testing.assert_allclose(probability_x, np.array([0.7, 0.4, 0.2]))
        np.testing.assert_array_equal(n_windows_x, np.array([3, 1, 1], dtype=np.int64))

        self.assertAlmostEqual(
            position_balanced_bayesian_decoding_accuracy(result, exclude_edge_bins=1),
            0.4,
        )
        self.assertAlmostEqual(
            position_balanced_bayesian_decoding_accuracy(result, include_bins=[0, 2]),
            0.45,
        )
        self.assertAlmostEqual(
            position_balanced_bayesian_decoding_accuracy(result, x_min_cm=5.0, x_max_cm=15.0),
            0.4,
        )
        self.assertAlmostEqual(
            position_balanced_bayesian_decoding_accuracy(result, min_windows_per_bin=2),
            0.7,
        )

        summary = summarize_decoding_measures(result)
        self.assertAlmostEqual(
            summary["position_balanced_bayesian_accuracy"],
            (0.7 + 0.4 + 0.2) / 3.0,
        )

        with self.assertRaises(ValueError):
            position_balanced_bayesian_decoding_accuracy(result, min_windows_per_bin=0)
        with self.assertRaises(ValueError):
            position_balanced_bayesian_decoding_accuracy(result, x_min_cm=15.0, x_max_cm=5.0)

    def test_known_two_bin_decoder_returns_normalized_posteriors(self) -> None:
        position_x = np.tile(
            np.array([0.5] * 10 + [1.5] * 10, dtype=np.float64),
            3,
        )
        speed = np.full(position_x.size, 3.0, dtype=np.float64)
        trials = [_trial(0, 0, 20), _trial(1, 20, 40), _trial(2, 40, 60)]
        spike_indices = np.array([2, 12, 22, 32, 42, 52], dtype=np.int64)
        spike_cell_ids = np.array([10, 11, 10, 11, 10, 11], dtype=np.int64)

        result = decode_bayesian_position_from_trials(
            position_x=position_x,
            spike_indices_0b=spike_indices,
            spike_cell_ids=spike_cell_ids,
            cell_ids=np.array([10, 11], dtype=np.int64),
            trials=trials,
            xbin_edges=np.array([0.0, 1.0, 2.0], dtype=np.float64),
            cfg=_cfg(),
            speed=speed,
        )

        self.assertEqual(result.posterior_wx.shape, (6, 2))
        np.testing.assert_allclose(np.sum(result.posterior_wx, axis=1), 1.0)
        np.testing.assert_array_equal(
            result.decoded_bin_w,
            np.array([0, 1, 0, 1, 0, 1], dtype=np.int64),
        )
        np.testing.assert_array_equal(result.actual_bin_w, result.decoded_bin_w)
        np.testing.assert_allclose(result.error_cm_w, 0.0)
        self.assertTrue(np.all(result.prob_actual_w > 0.5))

    def test_decoder_diagnostics_frames_summaries_and_soft_confusions(self) -> None:
        position_x = np.tile(
            np.array([0.5] * 10 + [1.5] * 10, dtype=np.float64),
            3,
        )
        speed = np.full(position_x.size, 3.0, dtype=np.float64)
        trials = [_trial(0, 0, 20), _trial(1, 20, 40), _trial(2, 40, 60)]
        result = decode_bayesian_position_from_trials(
            position_x=position_x,
            spike_indices_0b=np.array([2, 12, 22, 32, 42, 52], dtype=np.int64),
            spike_cell_ids=np.array([10, 11, 10, 11, 10, 11], dtype=np.int64),
            cell_ids=np.array([10, 11], dtype=np.int64),
            trials=trials,
            xbin_edges=np.array([0.0, 1.0, 2.0], dtype=np.float64),
            cfg=_cfg(),
            speed=speed,
        )

        frame = bayesian_decoding_result_to_frame(
            result,
            freq_hz=10.0,
            condition_names_by_base={1: "PO"},
        )
        self.assertEqual(len(frame), 6)
        self.assertEqual(set(frame["condition_label"]), {"PO W"})
        self.assertEqual(frame["window_mid_s"].iloc[0], 0.5)
        np.testing.assert_array_equal(frame["actual_bin"], frame["decoded_bin"])
        self.assertEqual(set(frame["cue_zone"]), {"rich", "poor"})
        self.assertNotIn("cue_layout", frame.columns)
        np.testing.assert_array_equal(
            frame["cue_zone_component"],
            np.array(["rich_1", "poor_1", "rich_1", "poor_1", "rich_1", "poor_1"]),
        )

        condition_summary = summarize_decoding_by_condition(frame)
        self.assertEqual(condition_summary["condition_label"].tolist(), ["PO W"])
        self.assertEqual(condition_summary["n_windows"].tolist(), [6])
        self.assertEqual(condition_summary["n_laps_decoded"].tolist(), [3])
        self.assertGreater(float(condition_summary.loc[0, "mean_prob_actual"]), 0.5)

        trial_summary = summarize_decoding_by_trial(frame)
        self.assertEqual(trial_summary["trial_index"].tolist(), [0, 1, 2])
        self.assertEqual(trial_summary["n_windows"].tolist(), [2, 2, 2])

        matrix, n_windows = soft_confusion_matrix_for_windows(
            result,
            np.arange(result.posterior_wx.shape[0]),
            direction="W",
        )
        self.assertEqual(n_windows, 6)
        self.assertEqual(matrix.shape, (2, 2))
        np.testing.assert_allclose(np.nansum(matrix, axis=1), np.array([1.0, 1.0]))
        self.assertGreater(matrix[0, 0], matrix[0, 1])
        self.assertGreater(matrix[1, 1], matrix[1, 0])

        average_confusions, trial_confusion = average_soft_confusions_by_condition(result, frame)
        self.assertEqual(set(average_confusions), {"PO W"})
        self.assertEqual(len(trial_confusion), 3)

        confusion_summary = summarize_soft_confusions(
            average_confusions,
            trial_confusion,
            result.xbin_centers,
            near_diagonal_cm=0.1,
        )
        self.assertEqual(confusion_summary["condition_label"].tolist(), ["PO W"])
        self.assertEqual(confusion_summary["n_trials"].tolist(), [3])
        self.assertGreater(float(confusion_summary.loc[0, "mean_diagonal_mass"]), 0.5)

    def test_leave_one_lap_out_training_excludes_held_out_lap(self) -> None:
        position_x = np.array([0.5] * 10 + [1.5] * 10, dtype=np.float64)
        speed = np.full(position_x.size, 3.0, dtype=np.float64)
        trials = [_trial(0, 0, 10), _trial(1, 10, 20)]

        result = decode_bayesian_position_from_trials(
            position_x=position_x,
            spike_indices_0b=np.array([2, 12], dtype=np.int64),
            spike_cell_ids=np.array([10, 10], dtype=np.int64),
            cell_ids=np.array([10], dtype=np.int64),
            trials=trials,
            xbin_edges=np.array([0.0, 1.0, 2.0], dtype=np.float64),
            cfg=_cfg(),
            speed=speed,
        )

        np.testing.assert_array_equal(result.actual_bin_w, np.array([0, 1], dtype=np.int64))
        np.testing.assert_array_equal(result.decoded_bin_w, np.array([1, 0], dtype=np.int64))
        np.testing.assert_array_equal(result.n_train_laps_w, np.array([1, 1], dtype=np.int64))

    def test_train_all_laps_includes_current_lap(self) -> None:
        position_x = np.array([0.5] * 10 + [1.5] * 10, dtype=np.float64)
        speed = np.full(position_x.size, 3.0, dtype=np.float64)
        trials = [_trial(0, 0, 10), _trial(1, 10, 20)]

        result = decode_bayesian_position_from_trials(
            position_x=position_x,
            spike_indices_0b=np.array([2, 12], dtype=np.int64),
            spike_cell_ids=np.array([10, 11], dtype=np.int64),
            cell_ids=np.array([10, 11], dtype=np.int64),
            trials=trials,
            xbin_edges=np.array([0.0, 1.0, 2.0], dtype=np.float64),
            cfg=_cfg(train_all_laps=True),
            speed=speed,
        )

        np.testing.assert_array_equal(result.actual_bin_w, np.array([0, 1], dtype=np.int64))
        np.testing.assert_array_equal(result.decoded_bin_w, np.array([0, 1], dtype=np.int64))
        np.testing.assert_array_equal(result.n_train_laps_w, np.array([2, 2], dtype=np.int64))

    def test_zero_spike_windows_keep_finite_normalized_posteriors(self) -> None:
        position_x = np.tile(
            np.array([0.5] * 10 + [1.5] * 10, dtype=np.float64),
            2,
        )
        speed = np.full(position_x.size, 3.0, dtype=np.float64)
        trials = [_trial(0, 0, 20), _trial(1, 20, 40)]

        result = decode_bayesian_position_from_trials(
            position_x=position_x,
            spike_indices_0b=np.array([2, 22], dtype=np.int64),
            spike_cell_ids=np.array([10, 10], dtype=np.int64),
            cell_ids=np.array([10], dtype=np.int64),
            trials=trials,
            xbin_edges=np.array([0.0, 1.0, 2.0], dtype=np.float64),
            cfg=_cfg(),
            speed=speed,
        )

        zero_spike = result.n_spikes_w == 0
        self.assertTrue(np.any(zero_spike))
        self.assertTrue(np.all(np.isfinite(result.posterior_wx[zero_spike])))
        np.testing.assert_allclose(np.sum(result.posterior_wx[zero_spike], axis=1), 1.0)

    def test_speed_filter_removes_low_speed_decode_windows(self) -> None:
        position_x = np.array([0.5] * 30, dtype=np.float64)
        speed = np.array([3.0] * 20 + [0.0] * 10, dtype=np.float64)
        trials = [_trial(0, 0, 10), _trial(1, 10, 20), _trial(2, 20, 30)]

        result = decode_bayesian_position_from_trials(
            position_x=position_x,
            spike_indices_0b=np.array([2, 12, 22], dtype=np.int64),
            spike_cell_ids=np.array([10, 10, 10], dtype=np.int64),
            cell_ids=np.array([10], dtype=np.int64),
            trials=trials,
            xbin_edges=np.array([0.0, 1.0, 2.0], dtype=np.float64),
            cfg=_cfg(),
            speed=speed,
        )

        np.testing.assert_array_equal(result.trial_index_w, np.array([0, 1], dtype=np.int64))
        self.assertNotIn(2, result.trial_index_w.tolist())

    def test_decode_groupby_condway_keeps_direction_groups_separate(self) -> None:
        position_x = np.array([0.5] * 10 + [1.5] * 10, dtype=np.float64)
        speed = np.full(position_x.size, 3.0, dtype=np.float64)
        trials = [
            TrialInfo(
                trial_index=0,
                cond=1,
                wb="W",
                condway=1,
                start_idx_0b=0,
                stop_idx_0b_exclusive=10,
            ),
            TrialInfo(
                trial_index=1,
                cond=1,
                wb="B",
                condway=2,
                start_idx_0b=10,
                stop_idx_0b_exclusive=20,
            ),
        ]

        kwargs = dict(
            position_x=position_x,
            spike_indices_0b=np.array([2, 12], dtype=np.int64),
            spike_cell_ids=np.array([10, 10], dtype=np.int64),
            cell_ids=np.array([10], dtype=np.int64),
            trials=trials,
            xbin_edges=np.array([0.0, 1.0, 2.0], dtype=np.float64),
            speed=speed,
        )

        by_condway = decode_bayesian_position_from_trials(cfg=_cfg(), **kwargs)
        self.assertEqual(by_condway.posterior_wx.shape[0], 0)

        global_decode = decode_bayesian_position_from_trials(
            cfg=_cfg(decode_groupby="global"),
            **kwargs,
        )
        self.assertEqual(global_decode.posterior_wx.shape[0], 2)
        np.testing.assert_array_equal(
            global_decode.n_train_laps_w,
            np.array([1, 1], dtype=np.int64),
        )

    def test_condition_family_grouping_pools_condition_variants_with_direction(self) -> None:
        position_x = np.array([0.5] * 10 + [1.5] * 10, dtype=np.float64)
        speed = np.full(position_x.size, 3.0, dtype=np.float64)
        trials = [
            TrialInfo(
                trial_index=0,
                cond=1,
                wb="W",
                condway=1,
                start_idx_0b=0,
                stop_idx_0b_exclusive=10,
            ),
            TrialInfo(
                trial_index=1,
                cond=2,
                wb="W",
                condway=3,
                start_idx_0b=10,
                stop_idx_0b_exclusive=20,
            ),
        ]
        kwargs = dict(
            position_x=position_x,
            spike_indices_0b=np.array([2, 12], dtype=np.int64),
            spike_cell_ids=np.array([10, 10], dtype=np.int64),
            cell_ids=np.array([10], dtype=np.int64),
            trials=trials,
            xbin_edges=np.array([0.0, 1.0, 2.0], dtype=np.float64),
            speed=speed,
        )

        exact = decode_bayesian_position_from_trials(cfg=_cfg(), **kwargs)
        self.assertEqual(exact.posterior_wx.shape[0], 0)

        family = decode_bayesian_position_from_trials(
            cfg=_cfg(group_condition_families=True),
            condition_names_by_base={1: "PO", 2: "PO2"},
            **kwargs,
        )
        self.assertEqual(family.posterior_wx.shape[0], 2)
        np.testing.assert_array_equal(
            family.condway_w,
            np.array([1, 3], dtype=np.int64),
        )
        self.assertEqual(set(family.train_group_label_w.astype(str)), {"PO W"})
        np.testing.assert_array_equal(
            family.n_train_laps_w,
            np.array([1, 1], dtype=np.int64),
        )

    def test_cli_smoke_builds_decoder_outputs(self) -> None:
        repo_root = Path(__file__).resolve().parents[1]
        script_path = repo_root / "scripts" / "pipelines" / "build_bayesian_decoder_from_interim.py"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            interim_root = tmp / "interim"
            out_root = tmp / "decode"
            session = "S1"
            session_dir = interim_root / "M1" / "2025-01-01"
            session_dir.mkdir(parents=True, exist_ok=True)

            np.savez_compressed(
                session_dir / f"{session}_allcel.npz",
                allcel__itime_spk=np.array([3, 13, 23, 33, 43, 53], dtype=np.int64),
                allcel__id_spk=np.array([10, 11, 10, 11, 10, 11], dtype=np.int64),
                allcel__id_cel=np.array([10, 11], dtype=np.int64),
            )
            vr = [[0.5] * 10 + [1.5] * 10 for _ in range(3)]
            speed = [[3.0] * 20 for _ in range(3)]
            np.savez_compressed(
                session_dir / f"{session}_trajdata.npz",
                traj__Cond=np.array([1, 1, 1], dtype=np.int64),
                traj__WB=np.array(["W", "W", "W"]),
                traj__condition=np.array(["PO", "PO", "PO"]),
                traj__start=np.array([1, 21, 41], dtype=np.int64),
                traj__stop=np.array([20, 40, 60], dtype=np.int64),
                traj__VRtraj__json=np.asarray(json.dumps(vr), dtype="S"),
                traj__XSpeed__json=np.asarray(json.dumps(speed), dtype="S"),
            )

            argv_prev = sys.argv[:]
            sys.argv = [
                str(script_path),
                "--interim_root",
                str(interim_root),
                "--out_root",
                str(out_root),
                "--no_run_subdir",
                "--spike_freq_hz",
                "10",
                "--behavior_freq_hz",
                "10",
                "--tau_s",
                "1.0",
                "--bin_size_cm",
                "1.0",
                "--smooth_sigma_bins",
                "0",
                "--min_speed",
                "2",
                "--no_normalize_x",
            ]
            try:
                runpy.run_path(str(script_path), run_name="__main__")
            finally:
                sys.argv = argv_prev

            decode_path = out_root / "M1" / "2025-01-01" / f"{session}_bayes_decode.npz"
            self.assertTrue(decode_path.exists())
            with np.load(decode_path, allow_pickle=False) as z:
                self.assertIn("decode__posterior_wx", z.files)
                self.assertIn("decode__decoded_x_w", z.files)
                self.assertIn("decode__prob_actual_w", z.files)
                self.assertIn("decode__train_group_w", z.files)
                self.assertIn("decode__train_group_label_w", z.files)
                self.assertIn("decode__condition_label_w", z.files)
                self.assertEqual(z["decode__posterior_wx"].shape, (6, 2))
                self.assertEqual(set(np.asarray(z["decode__condition_label_w"]).astype(str)), {"PO W"})
                self.assertEqual(set(np.asarray(z["decode__train_group_label_w"]).astype(str)), {"PO W"})

            summary_path = out_root / "decoding_summary.csv"
            self.assertTrue(summary_path.exists())
            summary = pd.read_csv(summary_path)
            self.assertEqual(summary["session_id"].tolist(), [session])
            self.assertEqual(summary["condition_label"].tolist(), ["PO W"])
            self.assertEqual(summary["train_group_label"].tolist(), ["PO W"])
            self.assertEqual(summary["group_condition_families"].tolist(), [False])
            self.assertEqual(int(summary.loc[0, "n_windows"]), 6)


if __name__ == "__main__":
    unittest.main()
