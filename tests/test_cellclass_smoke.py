from __future__ import annotations

import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from cellclass.config import (
    DEFAULT_AGE_GROUPS,
    DEFAULT_CELL_CLASSIFICATION_TABLE,
    DEFAULT_TYPE_U_COMPARISON_FEATURES,
    DEFAULT_TYPE_U_COMPARISON_N_INIT,
    DEFAULT_TYPE_U_COMPARISON_ROOT,
    age_group_from_age,
    normalize_type_u,
)
from cellclass.features import compute_cv2, waveform_features_from_bestswaveforms
from cellclass.pipeline import parse_session_id
from cellclass.processing import compute_acg


class CellclassSmokeTest(unittest.TestCase):
    def test_config_age_groups_and_type_u_normalization(self) -> None:
        self.assertEqual(DEFAULT_AGE_GROUPS, ("P16-18", "P19-21", "P22-24"))
        self.assertEqual(age_group_from_age(16), "P16-18")
        self.assertEqual(age_group_from_age(21), "P19-21")
        self.assertEqual(age_group_from_age(25), None)

        labels = normalize_type_u(pd.Series(["interneuron", "pyr", 1, "bad"]))
        self.assertEqual(labels.iloc[:3].tolist(), [0, 1, 1])
        self.assertTrue(pd.isna(labels.iloc[3]))

        self.assertEqual(DEFAULT_TYPE_U_COMPARISON_ROOT, "results/type_u_comparison_valero_feats_3")
        self.assertEqual(
            DEFAULT_CELL_CLASSIFICATION_TABLE,
            "results/type_u_comparison_valero_feats_3/cell_classification_table.csv",
        )
        self.assertEqual(
            DEFAULT_TYPE_U_COMPARISON_FEATURES,
            (
                "cv2",
                "acg_peak_latency_ms",
                "spk_duration_ms",
                "spk_asymmetry",
                "log_fr_hz_session",
            ),
        )
        self.assertEqual(DEFAULT_TYPE_U_COMPARISON_N_INIT, 15)

    def test_pipeline_parse_session_id(self) -> None:
        meta = parse_session_id(Path("VS57_2022-12-18_18-51-04_allcel.npz"))

        self.assertEqual(meta["mouse"], "VS57")
        self.assertEqual(meta["date"], "2022-12-18")
        self.assertEqual(meta["time"], "18-51-04")
        self.assertEqual(meta["session_id"], "VS57_2022-12-18_18-51-04")

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
