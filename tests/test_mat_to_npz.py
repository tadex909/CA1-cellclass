from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.io import savemat

from cellclass.mat_to_npz import (
    DEFAULT_MATILDE_TRAJ_FIELDS,
    convert_one,
    matilde_traj_to_records,
)


def _matilde_structs() -> tuple[dict[str, object], dict[str, object]]:
    return (
        {
            "idtrack_tr": np.array([[1, 3], [4, 7], [8, 9], [10, 12]]),
            "way_tr": np.array([0, 1, 0, 1]),
            "icond_tr": np.array([2, 2, 1, 1]),
            "icondw_tr": np.array([1, 2, 3, 4]),
            "time_d": np.arange(12, dtype=np.float64) / 100.0,
            "p_x_ds": np.arange(12, dtype=np.float64) + 10.0,
            "v_x_ds2": np.arange(12, dtype=np.float64) + 2.0,
            "vel_tx": np.arange(8, dtype=np.float64).reshape(4, 2),
            "prm": {"freq_d": 100.0},
        },
        {
            "cond": np.array(["PO", "PNO"], dtype=object),
        },
    )


class MatildeTrajAdapterTest(unittest.TestCase):
    def test_records_preserve_blocks_directions_and_behavior_frequency(self) -> None:
        bhv, eprm = _matilde_structs()
        records, meta = matilde_traj_to_records(bhv, eprm)

        self.assertEqual([r["Cond"] for r in records], [1, 1, 2, 2])
        self.assertEqual([r["WB"] for r in records], ["W", "B", "W", "B"])
        self.assertEqual(
            [r["condition"] for r in records],
            ["PO", "PO", "PNO", "PNO"],
        )
        np.testing.assert_allclose(records[0]["time"], [0.0, 0.01, 0.02])
        np.testing.assert_allclose(records[0]["VRtraj"], [10.0, 11.0, 12.0])
        np.testing.assert_allclose(records[0]["XSpeed"], [2.0, 3.0, 4.0])
        np.testing.assert_allclose(records[0]["binSpX"], [0.0, 1.0])
        self.assertEqual(meta["behavior_freq_hz"], 100.0)

    def test_convert_one_writes_trajdata_npz_from_bhv0_mat(self) -> None:
        bhv, eprm = _matilde_structs()
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            mat_path = root / "M1_2025-01-01_12-00-00_Bhv0.mat"
            savemat(mat_path, {"bhv": bhv, "eprm": eprm})

            out_path = convert_one(
                mat_path,
                out_root=root / "interim",
                overwrite=False,
                include_spikes=False,
                mode="traj_matilde",
                traj_fields_req=DEFAULT_MATILDE_TRAJ_FIELDS,
            )

            self.assertEqual(out_path.name, "M1_2025-01-01_12-00-00_trajdata.npz")
            with np.load(out_path, allow_pickle=False) as z:
                self.assertIn("traj__VRtraj__json", z.files)
                self.assertIn("traj__XSpeed__json", z.files)
                np.testing.assert_array_equal(z["traj__Cond"], [1, 1, 2, 2])
                np.testing.assert_array_equal(z["traj__WB"], ["W", "B", "W", "B"])
                meta = json.loads(z["meta_json"].tobytes().decode("utf-8"))
                self.assertEqual(meta["mode"], "traj_matilde")
                self.assertEqual(meta["traj_matilde"]["behavior_freq_hz"], 100.0)


if __name__ == "__main__":
    unittest.main()
