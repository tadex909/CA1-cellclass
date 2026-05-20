from __future__ import annotations

import unittest

import numpy as np

from cellclass.features import compute_cv2, waveform_features_from_bestswaveforms
from cellclass.processing import compute_acg


class CellclassSmokeTest(unittest.TestCase):
    def test_compute_acg_has_expected_shape(self) -> None:
        spike_times_s = np.array([0.01, 0.02, 0.04, 0.11, 0.13, 0.16], dtype=np.float64)
        spike_cluster_ids = np.array([1, 1, 1, 2, 2, 2], dtype=np.int64)

        result = compute_acg(
            spike_times_s=spike_times_s,
            spike_cluster_ids=spike_cluster_ids,
            cell_ids=np.array([1, 2], dtype=np.int64),
            bin_ms=10.0,
            window_ms=30.0,
            normalize="count",
        )

        self.assertEqual(result.acg.shape, (7, 2))
        np.testing.assert_array_equal(result.cell_ids, np.array([1, 2], dtype=np.int64))
        self.assertTrue(np.all(result.acg >= 0.0))

    def test_basic_feature_extractors_return_finite_values(self) -> None:
        spike_times_s = np.array([0.0, 1.0, 2.0, 0.0, 2.0, 5.0], dtype=np.float64)
        spike_cluster_ids = np.array([1, 1, 1, 2, 2, 2], dtype=np.int64)

        cell_ids, cv2 = compute_cv2(spike_times_s, spike_cluster_ids)
        np.testing.assert_array_equal(cell_ids, np.array([1, 2], dtype=np.int64))
        self.assertTrue(np.all(np.isfinite(cv2)))

        base_waveform = np.array([0.2, 0.0, -1.0, 0.5, 0.25], dtype=np.float64)
        bestswaveforms = np.stack(
            [base_waveform, base_waveform * 0.9, base_waveform * 1.1],
            axis=1,
        )[:, :, np.newaxis]

        waveform_result = waveform_features_from_bestswaveforms(
            bestswaveforms=bestswaveforms,
            fs_hz=25_000.0,
            trim=0,
        )

        self.assertEqual(waveform_result.spk_duration_ms.shape, (1,))
        self.assertTrue(bool(waveform_result.valid[0]))
        self.assertTrue(np.isfinite(waveform_result.spk_duration_ms[0]))
        self.assertTrue(np.isfinite(waveform_result.spk_peaktrough_ms[0]))
        self.assertTrue(np.isfinite(waveform_result.spk_asymmetry[0]))


if __name__ == "__main__":
    unittest.main()
