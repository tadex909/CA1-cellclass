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
    complement_spans,
    compute_condition_component_displacement_profiles,
    compute_condition_zone_displacement_profiles,
    compute_displacement_profile,
    PopulationGeometryConfig,
    build_population_geometry_for_group,
    build_population_geometry_from_saved_ratemap,
    cue_zone_layout_for_condition,
    compute_normalization_scales,
    compute_zone_displacement_profiles,
    decode_condway,
    filter_cell_ids_by_pred_type,
    label_xbin_centers_by_zone,
    label_xbin_centers_by_zone_component,
    load_cell_classification_table,
    load_saved_ratemap_pack,
    rebin_trial_maps,
    zone_component_names_for_layout,
)


def _save_synthetic_rmap(rmap_path: Path, traj_path: Path) -> None:
    counts_lkn = np.array(
        [
            [[3.0, 0.0], [0.0, 2.0], [1.0, 0.0], [0.0, 1.0]],
            [[2.0, 0.0], [0.0, 3.0], [0.0, 1.0], [1.0, 0.0]],
        ],
        dtype=np.float64,
    )
    dwell_lk = np.ones((2, 4), dtype=np.float64)
    counts_uxk = np.transpose(counts_lkn, (2, 0, 1))
    fr_uxk = counts_uxk / dwell_lk[None, :, :]

    counts_cxk = np.nanmean(counts_uxk, axis=1, keepdims=True)
    dwell_cxk = np.nanmean(dwell_lk, axis=0, keepdims=True)
    fr_cxk = np.nanmean(fr_uxk, axis=1, keepdims=True)

    payload = {
        "meta_json": np.array(json.dumps({"source_traj_npz": str(traj_path)}), dtype=np.string_),
        "cell_ids": np.array([10, 11], dtype=np.int64),
        "idcond_t": np.array([1, 1], dtype=np.int64),
        "xbin_edges": np.array([0.0, 1.0, 2.0, 3.0, 4.0], dtype=np.float64),
        "xbin_centers": np.array([0.5, 1.5, 2.5, 3.5], dtype=np.float64),
        "rmap__nbspk_tx_ux": counts_uxk,
        "rmap__dwell_tx_x": dwell_lk,
        "rmap__fr_tx_ux": fr_uxk,
        "rmap__nbspk_s_tx_ux": counts_uxk,
        "rmap__dwell_s_tx_x": dwell_lk,
        "rmap__fr_s_tx_ux": fr_uxk,
        "rmap__nbspk_cx_ux": counts_cxk,
        "rmap__dwell_cx_x": dwell_cxk,
        "rmap__fr_cx_ux": fr_cxk,
        "rmap__nbspk_s_cx_ux": counts_cxk,
        "rmap__dwell_s_cx_x": dwell_cxk,
        "rmap__fr_s_cx_ux": fr_cxk,
    }
    rmap_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(rmap_path, **payload)


def _save_synthetic_traj(traj_path: Path) -> None:
    traj_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        traj_path,
        traj__Cond=np.array([1, 1], dtype=np.int64),
        traj__condition=np.array(["PO2", "PO2"]),
    )


def _save_synthetic_classification_csv(csv_path: Path) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            {"session_id": "S1", "cell_id": 10, "pred_type": "pyramidal", "age_group": "P16-18"},
            {
                "session_id": "S1",
                "cell_id": 11,
                "pred_type": "interneuron",
                "age_group": "P16-18",
                "gmm_p_interneuron": 0.91,
                "gmm_p_pyramidal": 0.09,
            },
            {
                "session_id": "S2",
                "cell_id": 20,
                "pred_type": "pyramidal",
                "age_group": "P19-21",
                "gmm_p_interneuron": 0.47,
                "gmm_p_pyramidal": 0.53,
            },
        ]
    ).to_csv(csv_path, index=False)


class PopulationGeometryTest(unittest.TestCase):
    def test_decode_condway(self) -> None:
        self.assertEqual(decode_condway(1), (1, "W"))
        self.assertEqual(decode_condway(2), (1, "B"))
        self.assertEqual(decode_condway(7), (4, "W"))
        with self.assertRaises(ValueError):
            decode_condway(0)

    def test_rebin_preserves_total_counts_and_dwell(self) -> None:
        counts = np.arange(2 * 4 * 3, dtype=np.float64).reshape(2, 4, 3)
        dwell = np.arange(1, 9, dtype=np.float64).reshape(2, 4)
        rebinned_counts, rebinned_dwell, edges, centers = rebin_trial_maps(
            counts_lkn=counts,
            dwell_lk=dwell,
            xbin_edges=np.array([0.0, 1.0, 2.0, 3.0, 4.0], dtype=np.float64),
            n_geom_bins=2,
        )

        np.testing.assert_allclose(np.sum(rebinned_counts, axis=1), np.sum(counts, axis=1))
        np.testing.assert_allclose(np.sum(rebinned_dwell, axis=1), np.sum(dwell, axis=1))
        np.testing.assert_allclose(edges, np.array([0.0, 2.0, 4.0], dtype=np.float64))
        np.testing.assert_allclose(centers, np.array([1.0, 3.0], dtype=np.float64))

    def test_cue_zone_layout_records_po_and_pom_definitions(self) -> None:
        po_layout = cue_zone_layout_for_condition("PO2")
        self.assertEqual(po_layout.condition_family, "PO")
        self.assertEqual(po_layout.rich_spans, ((13.0, 43.0), (81.0, 96.0)))
        self.assertEqual(po_layout.excluded_spans, ((0.0, 10.0),))
        self.assertEqual(po_layout.object_centers, (20.0, 36.0, 88.0))
        self.assertEqual(
            po_layout.poor_spans,
            ((10.0, 13.0), (43.0, 81.0), (96.0, 100.0)),
        )

        pom_layout = cue_zone_layout_for_condition("POM")
        self.assertEqual(pom_layout.condition_family, "POM")
        self.assertEqual(pom_layout.rich_spans, ((13.0, 28.0), (57.0, 96.0)))
        self.assertEqual(pom_layout.excluded_spans, ((0.0, 10.0),))
        self.assertEqual(pom_layout.object_centers, (20.0, 64.0, 88.0))
        self.assertEqual(
            pom_layout.poor_spans,
            ((10.0, 13.0), (28.0, 57.0), (96.0, 100.0)),
        )

        self.assertEqual(complement_spans(tuple()), ((0.0, 100.0),))

    def test_label_xbin_centers_by_zone_uses_bin_centers(self) -> None:
        centers = np.array([5.0, 20.0, 50.0, 90.0], dtype=np.float64)
        labels = label_xbin_centers_by_zone(centers, condition_name="PO3")
        np.testing.assert_array_equal(
            labels,
            np.array(["", "rich", "poor", "rich"], dtype=np.str_),
        )

    def test_label_xbin_centers_by_zone_component_splits_contiguous_spans(self) -> None:
        layout = cue_zone_layout_for_condition("PO2")
        self.assertEqual(
            zone_component_names_for_layout(layout),
            ("rich_1", "rich_2", "poor_1", "poor_2", "poor_3"),
        )

        centers = np.array([5.0, 12.0, 20.0, 50.0, 90.0, 97.0], dtype=np.float64)
        labels = label_xbin_centers_by_zone_component(centers, layout=layout)
        np.testing.assert_array_equal(
            labels,
            np.array(["", "poor_1", "rich_1", "poor_2", "rich_2", "poor_3"], dtype=np.str_),
        )

    def test_pno_layout_excludes_first_track_segment_from_zone_labels(self) -> None:
        layout = cue_zone_layout_for_condition("PNO")
        self.assertEqual(layout.excluded_spans, ((0.0, 10.0),))
        self.assertEqual(layout.poor_spans, ((10.0, 100.0),))
        labels = label_xbin_centers_by_zone(
            np.array([5.0, 12.0, 50.0], dtype=np.float64),
            layout=layout,
        )
        np.testing.assert_array_equal(
            labels,
            np.array(["", "poor", "poor"], dtype=np.str_),
        )

    def test_compute_displacement_profile_uses_positive_lags_only(self) -> None:
        sim = np.array(
            [
                [1.0, 0.9, 0.4, 0.1],
                [0.9, 1.0, 0.5, 0.3],
                [0.4, 0.5, 1.0, 0.8],
                [0.1, 0.3, 0.8, 1.0],
            ],
            dtype=np.float64,
        )
        valid_mask = np.array([True, True, False, True], dtype=bool)
        profile = compute_displacement_profile(
            sim,
            xbin_centers=np.array([0.0, 1.0, 2.0, 3.0], dtype=np.float64),
            valid_mask_k=valid_mask,
        )

        np.testing.assert_array_equal(profile.delta_bin, np.array([1, 2, 3], dtype=np.int64))
        np.testing.assert_allclose(profile.delta_x, np.array([1.0, 2.0, 3.0], dtype=np.float64))
        np.testing.assert_array_equal(profile.n_pairs, np.array([1, 1, 1], dtype=np.int64))
        np.testing.assert_allclose(profile.mean, np.array([0.9, 0.3, 0.1], dtype=np.float64))

    def test_zone_displacement_profiles_distinguish_anchor_and_within(self) -> None:
        sim = np.array(
            [
                [1.0, 0.9, 0.4, 0.1],
                [0.9, 1.0, 0.5, 0.3],
                [0.4, 0.5, 1.0, 0.8],
                [0.1, 0.3, 0.8, 1.0],
            ],
            dtype=np.float64,
        )
        labels = np.array(["rich", "rich", "poor", "poor"], dtype=np.str_)
        profiles = compute_zone_displacement_profiles(
            sim,
            zone_labels_k=labels,
            xbin_centers=np.array([0.0, 1.0, 2.0, 3.0], dtype=np.float64),
        )

        np.testing.assert_allclose(
            profiles.all_profile.mean,
            np.array([(0.9 + 0.5 + 0.8) / 3.0, 0.35, 0.1], dtype=np.float64),
        )
        np.testing.assert_array_equal(
            profiles.anchor_profiles["rich"].n_pairs,
            np.array([2, 2, 1], dtype=np.int64),
        )
        np.testing.assert_allclose(
            profiles.anchor_profiles["rich"].mean,
            np.array([0.7, 0.35, 0.1], dtype=np.float64),
        )
        np.testing.assert_array_equal(
            profiles.within_profiles["rich"].n_pairs,
            np.array([1, 0, 0], dtype=np.int64),
        )
        self.assertAlmostEqual(float(profiles.within_profiles["rich"].mean[0]), 0.9, places=6)
        self.assertTrue(np.isnan(profiles.within_profiles["rich"].mean[1]))
        self.assertTrue(np.isnan(profiles.within_profiles["rich"].mean[2]))
        np.testing.assert_allclose(
            profiles.within_profiles["poor"].mean[:1],
            np.array([0.8], dtype=np.float64),
        )

    def test_condition_zone_profiles_use_canonical_family_mapping(self) -> None:
        sim = np.array(
            [
                [1.0, 0.6, 0.2],
                [0.6, 1.0, 0.4],
                [0.2, 0.4, 1.0],
            ],
            dtype=np.float64,
        )
        centers = np.array([20.0, 60.0, 90.0], dtype=np.float64)
        profiles = compute_condition_zone_displacement_profiles(
            sim,
            xbin_centers=centers,
            condition_name="PO2",
        )

        self.assertIsNotNone(profiles.layout)
        self.assertEqual(profiles.layout.condition_family, "PO")
        np.testing.assert_array_equal(
            profiles.zone_labels_k,
            np.array(["rich", "poor", "rich"], dtype=np.str_),
        )

    def test_condition_component_profiles_return_separate_rich_components(self) -> None:
        sim = np.array(
            [
                [1.0, 0.7, 0.2, 0.1],
                [0.7, 1.0, 0.3, 0.2],
                [0.2, 0.3, 1.0, 0.8],
                [0.1, 0.2, 0.8, 1.0],
            ],
            dtype=np.float64,
        )
        centers = np.array([20.0, 22.0, 90.0, 92.0], dtype=np.float64)
        profiles = compute_condition_component_displacement_profiles(
            sim,
            xbin_centers=centers,
            condition_name="PO3",
            include_rich=True,
            include_poor=False,
        )

        self.assertIsNotNone(profiles.layout)
        np.testing.assert_array_equal(
            profiles.zone_labels_k,
            np.array(["rich_1", "rich_1", "rich_2", "rich_2"], dtype=np.str_),
        )
        self.assertEqual(set(profiles.within_profiles), {"rich_1", "rich_2"})
        np.testing.assert_array_equal(
            profiles.within_profiles["rich_1"].n_pairs,
            np.array([1, 0, 0], dtype=np.int64),
        )
        np.testing.assert_array_equal(
            profiles.within_profiles["rich_2"].n_pairs,
            np.array([1, 0, 0], dtype=np.int64),
        )
        self.assertAlmostEqual(float(profiles.within_profiles["rich_1"].mean[0]), 0.7, places=6)
        self.assertAlmostEqual(float(profiles.within_profiles["rich_2"].mean[0]), 0.8, places=6)

    def test_mean_geometry_uses_pooled_counts_over_pooled_occupancy(self) -> None:
        counts = np.array(
            [
                [[10.0, 0.0], [0.0, 10.0]],
                [[0.0, 1.0], [0.0, 0.0]],
            ],
            dtype=np.float64,
        )
        dwell = np.array(
            [
                [10.0, 10.0],
                [1.0, 10.0],
            ],
            dtype=np.float64,
        )
        cfg = PopulationGeometryConfig(
            n_geom_bins=2,
            min_occupancy_s=0.0,
            normalizations=("raw",),
            n_splits=10,
            max_exact_splits=128,
            seed=0,
        )
        result = build_population_geometry_for_group(
            counts_lkn=counts,
            dwell_lk=dwell,
            cfg=cfg,
            condway_1b=1,
        )

        expected = 1.0 / np.sqrt(101.0)
        self.assertAlmostEqual(float(result.G_nkk[0, 0, 1]), expected, places=6)
        self.assertNotAlmostEqual(float(result.G_nkk[0, 0, 1]), 1.0 / np.sqrt(2.0), places=3)

    def test_diagonal_is_one_only_for_valid_nonzero_bins(self) -> None:
        counts = np.array(
            [
                [[1.0], [0.0], [0.0]],
                [[1.0], [0.0], [0.0]],
            ],
            dtype=np.float64,
        )
        dwell = np.array(
            [
                [1.0, 1.0, 0.0],
                [1.0, 1.0, 0.0],
            ],
            dtype=np.float64,
        )
        cfg = PopulationGeometryConfig(
            n_geom_bins=3,
            min_occupancy_s=0.0,
            normalizations=("raw",),
            n_splits=10,
            max_exact_splits=128,
            seed=0,
        )
        result = build_population_geometry_for_group(
            counts_lkn=counts,
            dwell_lk=dwell,
            cfg=cfg,
            condway_1b=1,
        )

        diag = np.diag(result.G_nkk[0])
        self.assertEqual(diag[0], 1.0)
        self.assertTrue(np.isnan(diag[1]))
        self.assertTrue(np.isnan(diag[2]))

    def test_cross_validated_geometry_is_symmetric_and_diagonal_not_forced(self) -> None:
        counts = np.array(
            [
                [[10.0, 0.0], [0.0, 10.0]],
                [[10.0, 0.0], [0.0, 10.0]],
                [[0.0, 10.0], [10.0, 0.0]],
                [[0.0, 10.0], [10.0, 0.0]],
            ],
            dtype=np.float64,
        )
        dwell = np.ones((4, 2), dtype=np.float64)
        cfg = PopulationGeometryConfig(
            n_geom_bins=2,
            min_occupancy_s=0.0,
            normalizations=("raw",),
            n_splits=100,
            max_exact_splits=128,
            seed=0,
        )
        result = build_population_geometry_for_group(
            counts_lkn=counts,
            dwell_lk=dwell,
            cfg=cfg,
            condway_1b=1,
        )

        s_cv = result.S_cv_nkk[0]
        np.testing.assert_allclose(s_cv, s_cv.T, equal_nan=True)
        diag = np.diag(s_cv)
        self.assertTrue(np.any(np.isfinite(diag)))
        self.assertTrue(np.any(diag < 0.999))

    def test_low_occupancy_bins_propagate_to_nan(self) -> None:
        counts = np.array(
            [
                [[1.0], [2.0]],
                [[1.0], [2.0]],
            ],
            dtype=np.float64,
        )
        dwell = np.array(
            [
                [1.0, 0.01],
                [1.0, 0.01],
            ],
            dtype=np.float64,
        )
        cfg = PopulationGeometryConfig(
            n_geom_bins=2,
            min_occupancy_s=0.05,
            normalizations=("raw",),
            n_splits=10,
            max_exact_splits=128,
            seed=0,
        )
        result = build_population_geometry_for_group(
            counts_lkn=counts,
            dwell_lk=dwell,
            cfg=cfg,
            condway_1b=1,
        )

        self.assertFalse(bool(result.valid_mean_nk[0, 1]))
        self.assertFalse(bool(result.valid_cv_nk[0, 1]))
        self.assertTrue(np.all(np.isnan(result.G_nkk[0, 1, :])))
        self.assertTrue(np.all(np.isnan(result.G_nkk[0, :, 1])))
        self.assertTrue(np.all(np.isnan(result.S_cv_nkk[0, 1, :])))
        self.assertTrue(np.all(np.isnan(result.S_cv_nkk[0, :, 1])))

    def test_mean_rate_normalization_changes_geometry_and_handles_silent_cells(self) -> None:
        counts = np.array(
            [
                [[10.0, 1.0, 0.0], [10.0, 0.0, 0.0]],
                [[10.0, 1.0, 0.0], [10.0, 0.0, 0.0]],
            ],
            dtype=np.float64,
        )
        dwell = np.ones((2, 2), dtype=np.float64)
        cfg = PopulationGeometryConfig(
            n_geom_bins=2,
            min_occupancy_s=0.0,
            normalizations=("raw", "mean_rate"),
            n_splits=10,
            max_exact_splits=128,
            seed=0,
        )
        scales = compute_normalization_scales(
            counts_lkn=counts,
            dwell_lk=dwell,
            min_occupancy_s=0.0,
            normalizations=cfg.normalizations,
        )
        self.assertEqual(float(scales["mean_rate"][2]), 1.0)

        result = build_population_geometry_for_group(
            counts_lkn=counts,
            dwell_lk=dwell,
            cfg=cfg,
            condway_1b=1,
            normalization_scales=scales,
        )

        raw_val = float(result.G_nkk[0, 0, 1])
        norm_val = float(result.G_nkk[1, 0, 1])
        self.assertTrue(np.isfinite(norm_val))
        self.assertNotAlmostEqual(raw_val, norm_val, places=4)

    def test_single_lap_still_returns_mean_geometry_but_not_cv(self) -> None:
        counts = np.array([[[1.0, 0.0], [0.0, 1.0]]], dtype=np.float64)
        dwell = np.ones((1, 2), dtype=np.float64)
        cfg = PopulationGeometryConfig(
            n_geom_bins=2,
            min_occupancy_s=0.0,
            normalizations=("raw",),
            n_splits=10,
            max_exact_splits=128,
            seed=0,
        )
        result = build_population_geometry_for_group(
            counts_lkn=counts,
            dwell_lk=dwell,
            cfg=cfg,
            condway_1b=1,
        )

        self.assertEqual(result.cv_status, "too_few_laps_for_cv")
        self.assertEqual(result.n_splits, 0)
        self.assertTrue(np.all(np.isfinite(result.G_nkk[0])))
        self.assertTrue(np.all(np.isnan(result.S_cv_nkk[0])))

    def test_load_classification_table_and_filter_pred_type(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            root = tmp / "classifications"
            _save_synthetic_classification_csv(root / "P16-18" / "p16_18_classification_info.csv")
            table = load_cell_classification_table(root)

            self.assertEqual(
                table.columns.tolist(),
                [
                    "session_id",
                    "cell_id",
                    "pred_type",
                    "p_pred_type",
                    "Sure (P(pred_type) > 0.6)",
                    "age_group",
                ],
            )
            self.assertEqual(len(table), 3)
            self.assertAlmostEqual(
                float(table.loc[table["cell_id"] == 11, "p_pred_type"].iloc[0]),
                0.91,
            )
            self.assertAlmostEqual(
                float(table.loc[table["cell_id"] == 20, "p_pred_type"].iloc[0]),
                0.53,
            )
            self.assertTrue(bool(table.loc[table["cell_id"] == 11, "Sure (P(pred_type) > 0.6)"].iloc[0]))
            self.assertFalse(bool(table.loc[table["cell_id"] == 20, "Sure (P(pred_type) > 0.6)"].iloc[0]))

            selected = filter_cell_ids_by_pred_type(
                session_id="S1",
                cell_ids=np.array([10, 11, 12], dtype=np.int64),
                classification_table=table,
                pred_types=("pyramidal",),
            )
            np.testing.assert_array_equal(selected, np.array([10], dtype=np.int64))

    def test_saved_ratemap_builder_can_subset_selected_cells(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            traj_path = tmp / "interim" / "M1" / "2025-01-01" / "S1_trajdata.npz"
            rmap_path = tmp / "ratemap" / "M1" / "2025-01-01" / "S1_rmap.npz"
            _save_synthetic_traj(traj_path)
            _save_synthetic_rmap(rmap_path, traj_path)
            saved = load_saved_ratemap_pack(rmap_path)
            cfg = PopulationGeometryConfig(
                n_geom_bins=2,
                min_occupancy_s=0.0,
                normalizations=("raw", "mean_rate"),
                n_splits=10,
                max_exact_splits=128,
                seed=0,
            )

            result = build_population_geometry_from_saved_ratemap(
                saved=saved,
                cfg=cfg,
                condition_names_by_base={1: "PO"},
                selected_cell_ids=np.array([10], dtype=np.int64),
            )

            np.testing.assert_array_equal(result.cell_ids_u, np.array([10], dtype=np.int64))
            self.assertEqual(result.summary_rows(session_id="S1")[0]["n_cells_used"], 1)

    def test_cli_smoke_builds_geometry_outputs(self) -> None:
        repo_root = Path(__file__).resolve().parents[1]
        script_path = repo_root / "scripts" / "pipelines" / "build_population_geometry_from_ratemap.py"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            ratemap_root = tmp / "ratemap"
            out_root = tmp / "geom"
            traj_path = tmp / "interim" / "M1" / "2025-01-01" / "S1_trajdata.npz"
            rmap_path = ratemap_root / "M1" / "2025-01-01" / "S1_rmap.npz"
            _save_synthetic_traj(traj_path)
            _save_synthetic_rmap(rmap_path, traj_path)

            argv_prev = sys.argv[:]
            sys.argv = [
                str(script_path),
                "--ratemap_root",
                str(ratemap_root),
                "--out_root",
                str(out_root),
                "--no_run_subdir",
                "--n_geom_bins",
                "2",
                "--seed",
                "0",
            ]
            try:
                runpy.run_path(str(script_path), run_name="__main__")
            finally:
                sys.argv = argv_prev

            geom_files = sorted(out_root.rglob("*_geom.npz"))
            self.assertEqual(len(geom_files), 1)
            with np.load(geom_files[0], allow_pickle=False) as z:
                self.assertIn("geom__G_gnkk", z.files)
                self.assertIn("geom__S_cv_gnkk", z.files)
                self.assertIn("geom__condition_family_g", z.files)
                self.assertEqual(str(z["geom__condition_family_g"][0]), "PO")

            summary_csv = out_root / "geometry_summary.csv"
            self.assertTrue(summary_csv.exists())
            d = pd.read_csv(summary_csv)
            self.assertEqual(set(d["normalization"].tolist()), {"raw", "mean_rate"})
            self.assertEqual(set(d["condition_family"].tolist()), {"PO"})

    def test_cli_smoke_builds_geometry_outputs_with_pred_type_filter(self) -> None:
        repo_root = Path(__file__).resolve().parents[1]
        script_path = repo_root / "scripts" / "pipelines" / "build_population_geometry_from_ratemap.py"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            ratemap_root = tmp / "ratemap"
            out_root = tmp / "geom"
            traj_path = tmp / "interim" / "M1" / "2025-01-01" / "S1_trajdata.npz"
            rmap_path = ratemap_root / "M1" / "2025-01-01" / "S1_rmap.npz"
            classification_root = tmp / "classifications"
            _save_synthetic_traj(traj_path)
            _save_synthetic_rmap(rmap_path, traj_path)
            _save_synthetic_classification_csv(
                classification_root / "P16-18" / "p16_18_classification_info.csv"
            )

            argv_prev = sys.argv[:]
            sys.argv = [
                str(script_path),
                "--ratemap_root",
                str(ratemap_root),
                "--out_root",
                str(out_root),
                "--no_run_subdir",
                "--n_geom_bins",
                "2",
                "--seed",
                "0",
                "--cell_classification_source",
                str(classification_root),
                "--cell_pred_types",
                "pyramidal",
            ]
            try:
                runpy.run_path(str(script_path), run_name="__main__")
            finally:
                sys.argv = argv_prev

            geom_files = sorted(out_root.rglob("*_geom.npz"))
            self.assertEqual(len(geom_files), 1)
            with np.load(geom_files[0], allow_pickle=False) as z:
                np.testing.assert_array_equal(z["geom__cell_ids_u"], np.array([10], dtype=np.int64))
                meta = json.loads(z["meta_json"].tobytes().decode("utf-8"))
                self.assertEqual(meta["n_cells_total"], 2)
                self.assertEqual(meta["n_cells_used"], 1)
                self.assertEqual(meta["cell_pred_types"], ["pyramidal"])

            summary_csv = out_root / "geometry_summary.csv"
            d = pd.read_csv(summary_csv)
            self.assertEqual(set(d["n_cells_used"].tolist()), {1})


if __name__ == "__main__":
    unittest.main()
