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
    TrialInfo,
    decode_bayesian_position_from_trials,
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


class BayesianDecodingTest(unittest.TestCase):
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
