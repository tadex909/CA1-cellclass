from __future__ import annotations

import unittest

import numpy as np

from placefields import RatemapConfig, TrialInfo, build_ratemap_from_trials


class PlacefieldsSmokeTest(unittest.TestCase):
    def test_build_ratemap_from_trials_returns_consistent_shapes(self) -> None:
        position_x = np.array([0.2, 0.8, 1.2, 1.8, 0.3, 0.9, 1.4, 1.9], dtype=np.float64)
        spike_indices_0b = np.array([1, 2, 5, 6], dtype=np.int64)
        spike_cell_ids = np.array([10, 10, 10, 10], dtype=np.int64)
        cell_ids = np.array([10], dtype=np.int64)
        trials = [
            TrialInfo(
                trial_index=0,
                cond=1,
                wb="W",
                condway=1,
                start_idx_0b=0,
                stop_idx_0b_exclusive=4,
            ),
            TrialInfo(
                trial_index=1,
                cond=1,
                wb="B",
                condway=2,
                start_idx_0b=4,
                stop_idx_0b_exclusive=8,
            ),
        ]
        xbin_edges = np.array([0.0, 1.0, 2.0], dtype=np.float64)
        cfg = RatemapConfig(
            freq_hz=1.0,
            smooth_sigma_bins=0.0,
            xbin_rem=0,
            nb_cond=2,
            min_speed=None,
        )

        pack = build_ratemap_from_trials(
            position_x=position_x,
            spike_indices_0b=spike_indices_0b,
            spike_cell_ids=spike_cell_ids,
            cell_ids=cell_ids,
            trials=trials,
            xbin_edges=xbin_edges,
            cfg=cfg,
            speed=None,
        )

        self.assertEqual(pack.nb_cond, 2)
        self.assertEqual(pack.nbspk_tx_ux.shape, (1, 2, 2))
        self.assertEqual(pack.fr_tx_ux.shape, (1, 2, 2))
        self.assertEqual(pack.fr_cx_ux.shape, (1, 2, 2))
        np.testing.assert_array_equal(pack.idcond_t, np.array([1, 2], dtype=np.int64))
        self.assertTrue(np.all(np.isfinite(pack.dwell_tx_x)))

        cell_view = pack.cell_view(cell_id=10)
        self.assertEqual(cell_view["fr_tx"].shape, (2, 2))
        self.assertEqual(cell_view["fr_cx"].shape, (2, 2))


if __name__ == "__main__":
    unittest.main()
