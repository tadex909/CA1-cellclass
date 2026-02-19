from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.mixture import GaussianMixture

from fitting import DEFAULT_AGE_GROUPS, prepare_matrix


FEATURE_SETS: dict[str, list[str]] = {
    "all_features": [
        "fr_hz",
        "burst_index",
        "cv2",
        "spk_duration_ms",
        "spk_peaktrough_ms",
        "spk_asymmetry",
        "refractory_ms_edge",
        "acg_peak_latency_ms",
    ],
    "some_features": [
        "fr_hz",
        "cv2",
        "acg_peak_latency_ms",
        "spk_duration_ms",
        "spk_asymmetry",
    ],
    "few_features": [
        "fr_hz",
        "spk_duration_ms",
        "refractory_ms_edge",
    ],
    "valero_features": [
        "cv2",
        "acg_peak_latency_ms",
        "spk_duration_ms",
        "spk_asymmetry",
        "fr_hz",
    ],
}


def parse_csv_list(raw: str | None) -> list[str]:
    if not raw:
        return []
    return [x.strip() for x in raw.split(",") if x.strip()]


def parse_int_list(raw: str | None, default: list[int]) -> list[int]:
    if not raw:
        return default
    return [int(x.strip()) for x in raw.split(",") if x.strip()]


def fit_gmm_once(
    X: np.ndarray,
    k: int,
    *,
    seed: int,
    n_init: int,
    covariance_type: str,
) -> dict[str, Any]:
    gmm = GaussianMixture(
        n_components=k,
        random_state=seed,
        n_init=n_init,
        covariance_type=covariance_type,
    )
    labels = gmm.fit_predict(X)

    bic = float(gmm.bic(X))

    sil = float("nan")
    uniq = np.unique(labels)
    if uniq.size > 1 and X.shape[0] > uniq.size:
        try:
            sil = float(silhouette_score(X, labels))
        except Exception:
            sil = float("nan")

    return {
        "labels": labels,
        "bic": bic,
        "silhouette": sil,
    }


def pairwise_mean_ari(label_list: list[np.ndarray]) -> float:
    if len(label_list) < 2:
        return float("nan")
    vals: list[float] = []
    for a, b in itertools.combinations(label_list, 2):
        vals.append(float(adjusted_rand_score(a, b)))
    return float(np.mean(vals)) if vals else float("nan")


def rank01(series: pd.Series, ascending: bool) -> pd.Series:
    # Returns score in [0,1], where 1 is best.
    if series.empty:
        return series
    valid = series.notna()
    out = pd.Series(np.nan, index=series.index, dtype=float)
    if valid.sum() == 0:
        return out

    r = series[valid].rank(method="average", ascending=ascending)
    n = len(r)
    if n == 1:
        out.loc[valid] = 1.0
        return out

    # rank 1 = best for ascending=True; invert when ascending=False behavior already in rank.
    out.loc[valid] = (n - r) / (n - 1)
    return out


def add_general_score(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    d["score_bic"] = np.nan
    d["score_silhouette"] = np.nan
    d["score_stability_ari"] = np.nan

    for age_group, sub_idx in d.groupby("age_group").groups.items():
        sub = d.loc[sub_idx].copy()
        # Lower BIC is better.
        sub["score_bic"] = rank01(sub["bic_mean"], ascending=True)
        # Higher silhouette/stability is better.
        sub["score_silhouette"] = rank01(-sub["silhouette_mean"], ascending=True)
        sub["score_stability_ari"] = rank01(-sub["stability_ari"], ascending=True)

        d.loc[sub_idx, "score_bic"] = sub["score_bic"]
        d.loc[sub_idx, "score_silhouette"] = sub["score_silhouette"]
        d.loc[sub_idx, "score_stability_ari"] = sub["score_stability_ari"]

    # Weighted mean score (equal weights here).
    d["general_score"] = d[["score_bic", "score_silhouette", "score_stability_ari"]].mean(
        axis=1, skipna=True
    )
    return d


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Run fixed feature-set GMM experiments (k values + seeds) and summarize "
            "silhouette, BIC, stability ARI, and general score."
        )
    )
    ap.add_argument("--results_root", type=str, default="results")
    ap.add_argument("--out_root", type=str, default="results/feature_set_experiments")
    ap.add_argument("--age_groups", type=str, default=",".join(DEFAULT_AGE_GROUPS))
    ap.add_argument("--k_values", type=str, default="1,2,3")
    ap.add_argument("--seeds", type=str, default="0,1,2,3,4")
    ap.add_argument("--n_init", type=int, default=10)
    ap.add_argument(
        "--covariance_type",
        type=str,
        choices=["full", "tied", "diag", "spherical"],
        default="full",
    )
    ap.add_argument("--no_log_fr", action="store_true")
    ap.add_argument("--no_standardize", action="store_true")
    ap.add_argument("--min_units", type=int, default=10)
    args = ap.parse_args()

    results_root = Path(args.results_root)
    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    age_groups = parse_csv_list(args.age_groups) or DEFAULT_AGE_GROUPS
    k_values = parse_int_list(args.k_values, default=[1, 2, 3])
    seeds = parse_int_list(args.seeds, default=[0, 1, 2, 3, 4])

    rows: list[dict[str, Any]] = []

    for age_group in age_groups:
        in_path = results_root / age_group / f"{age_group}_clean_units.parquet"
        if not in_path.exists():
            print(f"{age_group}: missing {in_path}, skipping")
            continue

        df = pd.read_parquet(in_path)

        for set_name, feature_cols in FEATURE_SETS.items():
            pack = prepare_matrix(
                df,
                feature_cols,
                log_fr=not args.no_log_fr,
                standardize=not args.no_standardize,
            )

            if pack.rows_after < args.min_units:
                for k in k_values:
                    rows.append(
                        {
                            "age_group": age_group,
                            "feature_set": set_name,
                            "feature_cols": "+".join(feature_cols),
                            "k": int(k),
                            "status": "too_few_units",
                            "n_units": int(pack.rows_after),
                            "bic_mean": np.nan,
                            "bic_std": np.nan,
                            "silhouette_mean": np.nan,
                            "silhouette_std": np.nan,
                            "stability_ari": np.nan,
                        }
                    )
                continue

            X = pack.X

            for k in k_values:
                if k < 1 or k > X.shape[0]:
                    rows.append(
                        {
                            "age_group": age_group,
                            "feature_set": set_name,
                            "feature_cols": "+".join(feature_cols),
                            "k": int(k),
                            "status": "invalid_k",
                            "n_units": int(X.shape[0]),
                            "bic_mean": np.nan,
                            "bic_std": np.nan,
                            "silhouette_mean": np.nan,
                            "silhouette_std": np.nan,
                            "stability_ari": np.nan,
                        }
                    )
                    continue

                run_metrics: list[dict[str, Any]] = []
                labels_list: list[np.ndarray] = []

                for seed in seeds:
                    try:
                        m = fit_gmm_once(
                            X,
                            k,
                            seed=seed,
                            n_init=args.n_init,
                            covariance_type=args.covariance_type,
                        )
                        run_metrics.append(m)
                        labels_list.append(m["labels"])
                    except Exception:
                        pass

                if not run_metrics:
                    rows.append(
                        {
                            "age_group": age_group,
                            "feature_set": set_name,
                            "feature_cols": "+".join(feature_cols),
                            "k": int(k),
                            "status": "fit_failed",
                            "n_units": int(X.shape[0]),
                            "bic_mean": np.nan,
                            "bic_std": np.nan,
                            "silhouette_mean": np.nan,
                            "silhouette_std": np.nan,
                            "stability_ari": np.nan,
                        }
                    )
                    continue

                bic_vals = np.array([m["bic"] for m in run_metrics], dtype=float)
                sil_vals = np.array([m["silhouette"] for m in run_metrics], dtype=float)
                stab = pairwise_mean_ari(labels_list)

                rows.append(
                    {
                        "age_group": age_group,
                        "feature_set": set_name,
                        "feature_cols": "+".join(feature_cols),
                        "k": int(k),
                        "status": "ok",
                        "n_units": int(X.shape[0]),
                        "n_runs": int(len(run_metrics)),
                        "bic_mean": float(np.nanmean(bic_vals)),
                        "bic_std": float(np.nanstd(bic_vals)),
                        "silhouette_mean": float(np.nanmean(sil_vals)),
                        "silhouette_std": float(np.nanstd(sil_vals)),
                        "stability_ari": float(stab),
                    }
                )

    out = pd.DataFrame(rows)
    out = add_general_score(out)
    out = out.sort_values(["age_group", "general_score"], ascending=[True, False]).reset_index(drop=True)

    csv_path = out_root / "feature_set_k_experiments.csv"
    out.to_csv(csv_path, index=False)

    cfg = {
        "results_root": str(results_root),
        "out_root": str(out_root),
        "age_groups": age_groups,
        "feature_sets": FEATURE_SETS,
        "k_values": k_values,
        "seeds": seeds,
        "n_init": int(args.n_init),
        "covariance_type": args.covariance_type,
        "log_fr": bool(not args.no_log_fr),
        "standardize": bool(not args.no_standardize),
        "min_units": int(args.min_units),
        "output_csv": str(csv_path),
    }
    (out_root / "run_config.json").write_text(json.dumps(cfg, indent=2))

    print(f"Wrote {csv_path}")


if __name__ == "__main__":
    main()
