from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from cellclass.validation import (
    ValidationError,
    validate_age_group_table,
    validate_allcel_npz,
    validate_feature_table,
)


def _valid_allcel_payload(n_cells: int = 2) -> dict[str, np.ndarray]:
    n_spikes = 4
    return {
        "allcel__time_spk": np.array([0.1, 0.2, 0.3, 0.4], dtype=np.float64),
        "allcel__id_spk": np.array([10, 11, 10, 11], dtype=np.int64),
        "allcel__id_cel": np.array([10, 11], dtype=np.int64)[:n_cells],
        "allcel__type_u": np.array([1, 0], dtype=np.int64)[:n_cells],
        "allcel__burst_u": np.ones(n_cells, dtype=np.float64),
        "allcel__fr_u": np.ones(n_cells, dtype=np.float64),
        "allcel__duration_u": np.ones(n_cells, dtype=np.float64),
        "allcel__asymmetry_u": np.ones(n_cells, dtype=np.float64),
        "allcel__bestswaveforms": np.ones((5, 3, n_cells), dtype=np.float64),
        "allpf__ispf_cxu": np.zeros((10, 4, n_cells), dtype=np.float64),
    }


def _valid_feature_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "session_id": ["S1", "S1"],
            "mouse": ["M1", "M1"],
            "date": ["2025-01-01", "2025-01-01"],
            "time": ["12-00-00", "12-00-00"],
            "cell_id": [10, 11],
            "allcel__type_u": [1, 0],
            "n_spikes": [200, 250],
            "fr_hz": [2.0, 4.0],
            "fr_hz_session": [2.1, 4.1],
            "cv2": [0.5, 0.6],
            "refractory_ms_center": [2.0, 1.5],
            "refractory_ms_edge": [2.5, 1.8],
            "burst_index": [0.1, 0.2],
            "acg_peak_latency_ms": [5.0, 6.0],
            "spk_duration_ms": [0.4, 0.2],
            "spk_peaktrough_ms": [0.3, 0.2],
            "spk_asymmetry": [0.1, 0.2],
            "qc_min_spikes": [True, True],
            "qc_refractory": [True, True],
            "qc_waveform": [True, True],
        }
    )


class ValidationTest(unittest.TestCase):
    def test_validate_allcel_npz_accepts_expected_payload(self) -> None:
        validate_allcel_npz(_valid_allcel_payload(), source="synthetic allcel")

    def test_validate_allcel_npz_reports_missing_key(self) -> None:
        payload = _valid_allcel_payload()
        payload.pop("allpf__ispf_cxu")

        with self.assertRaisesRegex(ValidationError, "missing required key"):
            validate_allcel_npz(payload, source="bad allcel")

    def test_validate_allcel_npz_reports_per_cell_length_mismatch(self) -> None:
        payload = _valid_allcel_payload()
        payload["allcel__type_u"] = np.array([1, 0, 1], dtype=np.int64)

        with self.assertRaisesRegex(ValidationError, "allcel__type_u length"):
            validate_allcel_npz(payload, source="bad allcel")

    def test_validate_feature_table_reports_duplicate_unit_rows(self) -> None:
        df = _valid_feature_table()
        df.loc[1, "cell_id"] = 10

        with self.assertRaisesRegex(ValidationError, "duplicate"):
            validate_feature_table(df, source="bad features")

    def test_validate_age_group_table_checks_age_and_required_features(self) -> None:
        df = _valid_feature_table()
        df["unit_uid"] = df["session_id"] + "__cell" + df["cell_id"].astype(str)
        df["Age"] = [16, 16]
        df["age_group"] = ["P16-18", "P19-21"]

        with self.assertRaisesRegex(ValidationError, "mismatched"):
            validate_age_group_table(
                df,
                source="bad age table",
                age_group="P16-18",
                required_features=["fr_hz", "spk_duration_ms"],
            )

    def test_validate_age_group_table_reports_missing_feature(self) -> None:
        df = _valid_feature_table()
        df["unit_uid"] = df["session_id"] + "__cell" + df["cell_id"].astype(str)
        df["Age"] = [16, 16]
        df["age_group"] = ["P16-18", "P16-18"]

        with self.assertRaisesRegex(ValidationError, "missing required column"):
            validate_age_group_table(
                df,
                source="bad age table",
                age_group="P16-18",
                required_features=["missing_feature"],
            )


if __name__ == "__main__":
    unittest.main()
