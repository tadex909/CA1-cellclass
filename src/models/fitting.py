from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import (
    calinski_harabasz_score,
    davies_bouldin_score,
    silhouette_score,
)
from sklearn.mixture import GaussianMixture


DEFAULT_AGE_GROUPS = ["P15_16", "P17_18", "P19_20", "P21_22", "P23_24", "P25"]
DEFAULT_FEATURES = [
    "fr_hz",
    "burst_index",
    "cv2",
    "spk_duration_ms",
    "spk_peaktrough_ms",
    "spk_asymmetry",
    "refractory_ms_edge",
    "acg_peak_latency_ms",
]


@dataclass
class MatrixPack:
    X: np.ndarray
    df_valid: pd.DataFrame
    used_cols: list[str]
    rows_before: int
    rows_after: int


def parse_csv_list(raw: str | None) -> list[str]:
    if not raw:
        return []
    return [x.strip() for x in raw.split(",") if x.strip()]


def iter_feature_subsets(
    features: list[str],
    mode: str,
    min_subset_size: int,
) -> list[list[str]]:
    if not features:
        return []

    if mode == "full_only":
        return [features]

    out: list[list[str]] = [features]

    if mode == "leave_one_out":
        if len(features) <= 1:
            return out
        out.extend([f for f in combinations(features, len(features) - 1)])
        return [list(x) for x in out]

    # mode == "all"
    for size in range(len(features) - 1, min_subset_size - 1, -1):
        out.extend(combinations(features, size))
    return [list(x) for x in out]


def prepare_matrix(
    df: pd.DataFrame,
    feature_cols: list[str],
    *,
    log_fr: bool,
    standardize: bool,
) -> MatrixPack:
    d = df.copy()
    used_cols: list[str] = []

    for col in feature_cols:
        if col == "fr_hz" and log_fr:
            d["log10_fr_hz"] = np.log10(d["fr_hz"].clip(lower=1e-6))
            used_cols.append("log10_fr_hz")
        elif col in d.columns:
            used_cols.append(col)

    if not used_cols:
        raise ValueError("No usable feature columns found in dataframe.")

    X = d[used_cols].to_numpy(dtype=np.float64)
    finite = np.isfinite(X).all(axis=1)
    d_valid = d.loc[finite].copy()
    X = X[finite]

    if standardize and X.size:
        mu = X.mean(axis=0)
        sd = X.std(axis=0, ddof=0)
        sd[sd == 0] = 1.0
        X = (X - mu) / sd

    return MatrixPack(
        X=X,
        df_valid=d_valid,
        used_cols=used_cols,
        rows_before=len(d),
        rows_after=len(d_valid),
    )


def safe_metric(metric_fn, X: np.ndarray, labels: np.ndarray) -> float:
    try:
        val = float(metric_fn(X, labels))
    except Exception:
        return float("nan")
    if not np.isfinite(val):
        return float("nan")
    return val


def evaluate_gmm(
    X: np.ndarray,
    k: int,
    *,
    random_state: int,
    n_init: int,
    covariance_type: str,
) -> dict[str, Any]:
    model = GaussianMixture(
        n_components=k,
        random_state=random_state,
        n_init=n_init,
        covariance_type=covariance_type,
    )
    labels = model.fit_predict(X)

    unique_labels, counts = np.unique(labels, return_counts=True)
    n_effective = int(unique_labels.size)
    min_cluster = int(counts.min()) if counts.size else 0
    max_cluster = int(counts.max()) if counts.size else 0

    sil = float("nan")
    cal = float("nan")
    dav = float("nan")
    if n_effective > 1 and X.shape[0] > n_effective:
        sil = safe_metric(silhouette_score, X, labels)
        cal = safe_metric(calinski_harabasz_score, X, labels)
        dav = safe_metric(davies_bouldin_score, X, labels)

    out = {
        "bic": float(model.bic(X)),
        "aic": float(model.aic(X)),
        "silhouette": sil,
        "calinski_harabasz": cal,
        "davies_bouldin": dav,
        "n_effective_clusters": n_effective,
        "min_cluster_size": min_cluster,
        "max_cluster_size": max_cluster,
        "labels": labels,
        "proba": model.predict_proba(X),
        "model": model,
    }
    return out


def rank_results(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    if d.empty:
        return d

    d["rank_bic"] = d["bic"].rank(method="average", ascending=True, na_option="bottom")
    d["rank_aic"] = d["aic"].rank(method="average", ascending=True, na_option="bottom")
    d["rank_silhouette"] = d["silhouette"].rank(
        method="average", ascending=False, na_option="bottom"
    )
    d["rank_calinski"] = d["calinski_harabasz"].rank(
        method="average", ascending=False, na_option="bottom"
    )
    d["rank_davies"] = d["davies_bouldin"].rank(
        method="average", ascending=True, na_option="bottom"
    )

    rank_cols = ["rank_bic", "rank_aic", "rank_silhouette", "rank_calinski", "rank_davies"]
    d["mean_rank"] = d[rank_cols].mean(axis=1, skipna=True)
    d = d.sort_values(["mean_rank", "bic", "aic"], ascending=[True, True, True]).reset_index(
        drop=True
    )
    return d


def save_best_labels(
    out_path: Path,
    matrix_pack: MatrixPack,
    labels: np.ndarray,
    proba: np.ndarray,
    *,
    age_group: str,
    subset_name: str,
    k: int,
) -> None:
    d = matrix_pack.df_valid.copy()
    d["gmm_cluster"] = labels
    d["gmm_pmax"] = proba.max(axis=1)
    d["age_group"] = age_group
    d["subset_name"] = subset_name
    d["k"] = int(k)

    # Optional mapping for 2-cluster results if spike duration exists.
    if k == 2 and "spk_duration_ms" in d.columns:
        cluster_mean_dur = (
            pd.DataFrame({"cluster": labels, "dur": d["spk_duration_ms"].to_numpy()})
            .groupby("cluster")["dur"]
            .mean()
        )
        interneuron_cluster = int(cluster_mean_dur.idxmin())
        d["putative_type"] = np.where(labels == interneuron_cluster, "interneuron", "pyramidal")

    keep_cols = [
        "session_id",
        "cell_id",
        "unit_uid",
        "allcel__type_u",
        "gmm_cluster",
        "gmm_pmax",
        "putative_type",
        "age_group",
        "subset_name",
        "k",
    ]
    keep_cols = [c for c in keep_cols if c in d.columns]
    d[keep_cols].to_parquet(out_path, index=False)


def evaluate_age_group(
    age_group: str,
    *,
    results_root: Path,
    out_root: Path,
    feature_pool: list[str],
    subset_mode: str,
    min_subset_size: int,
    k_values: list[int],
    min_units: int,
    log_fr: bool,
    standardize: bool,
    random_state: int,
    n_init: int,
    covariance_type: str,
    min_cluster_size: int,
) -> None:
    in_path = results_root / age_group / f"{age_group}_clean_units.parquet"
    if not in_path.exists():
        print(f"{age_group}: missing {in_path}, skipping")
        return

    df = pd.read_parquet(in_path)
    if "unit_uid" not in df.columns:
        df["unit_uid"] = df["session_id"].astype(str) + "__cell" + df["cell_id"].astype(str)

    out_dir = out_root / age_group
    out_dir.mkdir(parents=True, exist_ok=True)

    subsets = iter_feature_subsets(feature_pool, subset_mode, min_subset_size=min_subset_size)
    print(f"{age_group}: units={len(df)}, subsets={len(subsets)}, k_values={k_values}")

    rows: list[dict[str, Any]] = []
    pack_cache: dict[tuple[str, ...], MatrixPack] = {}

    for subset in subsets:
        subset_key = tuple(subset)
        subset_name = "+".join(subset)
        try:
            pack = prepare_matrix(df, subset, log_fr=log_fr, standardize=standardize)
        except Exception as exc:
            for k in k_values:
                rows.append(
                    {
                        "age_group": age_group,
                        "feature_subset": subset_name,
                        "k": int(k),
                        "status": "error_prepare",
                        "error": str(exc),
                        "n_rows_before": len(df),
                        "n_rows_after": 0,
                    }
                )
            continue

        pack_cache[subset_key] = pack
        if pack.rows_after < min_units:
            for k in k_values:
                rows.append(
                    {
                        "age_group": age_group,
                        "feature_subset": subset_name,
                        "k": int(k),
                        "status": "too_few_units",
                        "error": f"{pack.rows_after} < min_units {min_units}",
                        "n_rows_before": pack.rows_before,
                        "n_rows_after": pack.rows_after,
                        "used_cols": ",".join(pack.used_cols),
                    }
                )
            continue

        for k in k_values:
            if k < 2 or k >= pack.rows_after:
                rows.append(
                    {
                        "age_group": age_group,
                        "feature_subset": subset_name,
                        "k": int(k),
                        "status": "invalid_k",
                        "error": f"k={k} incompatible with n_rows_after={pack.rows_after}",
                        "n_rows_before": pack.rows_before,
                        "n_rows_after": pack.rows_after,
                        "used_cols": ",".join(pack.used_cols),
                    }
                )
                continue

            try:
                ev = evaluate_gmm(
                    pack.X,
                    k,
                    random_state=random_state,
                    n_init=n_init,
                    covariance_type=covariance_type,
                )
                meets_cluster_size = ev["min_cluster_size"] >= min_cluster_size
                status = "ok" if meets_cluster_size else "cluster_too_small"
                error = "" if meets_cluster_size else (
                    f"min_cluster_size={ev['min_cluster_size']} < required={min_cluster_size}"
                )
                rows.append(
                    {
                        "age_group": age_group,
                        "feature_subset": subset_name,
                        "k": int(k),
                        "status": status,
                        "error": error,
                        "n_rows_before": pack.rows_before,
                        "n_rows_after": pack.rows_after,
                        "used_cols": ",".join(pack.used_cols),
                        "bic": ev["bic"],
                        "aic": ev["aic"],
                        "silhouette": ev["silhouette"],
                        "calinski_harabasz": ev["calinski_harabasz"],
                        "davies_bouldin": ev["davies_bouldin"],
                        "n_effective_clusters": ev["n_effective_clusters"],
                        "min_cluster_size": ev["min_cluster_size"],
                        "max_cluster_size": ev["max_cluster_size"],
                    }
                )
            except Exception as exc:
                rows.append(
                    {
                        "age_group": age_group,
                        "feature_subset": subset_name,
                        "k": int(k),
                        "status": "error_fit",
                        "error": str(exc),
                        "n_rows_before": pack.rows_before,
                        "n_rows_after": pack.rows_after,
                        "used_cols": ",".join(pack.used_cols),
                    }
                )

    res = pd.DataFrame(rows)
    csv_path = out_dir / f"{age_group}_gmm_model_selection.csv"
    res.to_csv(csv_path, index=False)

    ok = res[res["status"] == "ok"].copy()
    if ok.empty:
        summary = {
            "age_group": age_group,
            "status": "no_valid_models",
            "n_rows_total": int(len(res)),
            "n_rows_ok": 0,
            "selection_csv": str(csv_path),
        }
        (out_dir / f"{age_group}_gmm_summary.json").write_text(json.dumps(summary, indent=2))
        print(f"{age_group}: no valid models")
        return

    ranked = rank_results(ok)
    ranked_path = out_dir / f"{age_group}_gmm_model_selection_ranked.csv"
    ranked.to_csv(ranked_path, index=False)

    best = ranked.iloc[0]
    best_subset = str(best["feature_subset"]).split("+")
    best_k = int(best["k"])
    best_pack = pack_cache[tuple(best_subset)]
    best_ev = evaluate_gmm(
        best_pack.X,
        best_k,
        random_state=random_state,
        n_init=n_init,
        covariance_type=covariance_type,
    )

    labels_path = out_dir / f"{age_group}_best_gmm_labels.parquet"
    save_best_labels(
        labels_path,
        best_pack,
        best_ev["labels"],
        best_ev["proba"],
        age_group=age_group,
        subset_name="+".join(best_subset),
        k=best_k,
    )

    summary = {
        "age_group": age_group,
        "status": "ok",
        "n_rows_total": int(len(res)),
        "n_rows_ok": int(len(ok)),
        "selection_csv": str(csv_path),
        "ranked_csv": str(ranked_path),
        "best_labels_parquet": str(labels_path),
        "best": {
            "feature_subset": str(best["feature_subset"]),
            "k": best_k,
            "bic": float(best["bic"]),
            "aic": float(best["aic"]),
            "silhouette": (
                None if (pd.isna(best["silhouette"]) or not np.isfinite(best["silhouette"])) else float(best["silhouette"])
            ),
            "calinski_harabasz": (
                None
                if (pd.isna(best["calinski_harabasz"]) or not np.isfinite(best["calinski_harabasz"]))
                else float(best["calinski_harabasz"])
            ),
            "davies_bouldin": (
                None
                if (pd.isna(best["davies_bouldin"]) or not np.isfinite(best["davies_bouldin"]))
                else float(best["davies_bouldin"])
            ),
            "mean_rank": float(best["mean_rank"]),
            "n_effective_clusters": int(best["n_effective_clusters"]),
            "min_cluster_size": int(best["min_cluster_size"]),
            "max_cluster_size": int(best["max_cluster_size"]),
            "n_rows_after": int(best["n_rows_after"]),
        },
    }
    (out_dir / f"{age_group}_gmm_summary.json").write_text(json.dumps(summary, indent=2))
    print(
        f"{age_group}: best subset={summary['best']['feature_subset']} "
        f"k={summary['best']['k']} mean_rank={summary['best']['mean_rank']:.3f}"
    )


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Systematic GMM model selection for age-group datasets: "
            "compare k values and feature subsets."
        )
    )
    ap.add_argument("--results_root", type=str, default="results")
    ap.add_argument("--out_root", type=str, default="results/model_selection")
    ap.add_argument("--age_groups", type=str, default=",".join(DEFAULT_AGE_GROUPS))
    ap.add_argument("--features", type=str, default=",".join(DEFAULT_FEATURES))
    ap.add_argument(
        "--subset_mode",
        type=str,
        choices=["full_only", "leave_one_out", "all"],
        default="leave_one_out",
        help=(
            "full_only: only all features; leave_one_out: all features and n-1 subsets; "
            "all: evaluate all subsets from n down to min_subset_size."
        ),
    )
    ap.add_argument("--min_subset_size", type=int, default=2)
    ap.add_argument("--k_min", type=int, default=2)
    ap.add_argument("--k_max", type=int, default=6)
    ap.add_argument("--min_units", type=int, default=30)
    ap.add_argument(
        "--min_cluster_size",
        type=int,
        default=3,
        help="Reject candidate models where the smallest cluster has fewer units than this.",
    )
    ap.add_argument("--no_log_fr", action="store_true")
    ap.add_argument("--no_standardize", action="store_true")
    ap.add_argument("--random_state", type=int, default=0)
    ap.add_argument("--n_init", type=int, default=10)
    ap.add_argument(
        "--covariance_type",
        type=str,
        choices=["full", "tied", "diag", "spherical"],
        default="full",
    )
    return ap


def main() -> None:
    ap = build_parser()
    args = ap.parse_args()

    results_root = Path(args.results_root)
    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    age_groups = parse_csv_list(args.age_groups) or DEFAULT_AGE_GROUPS
    feature_pool = parse_csv_list(args.features) or DEFAULT_FEATURES

    if args.k_min < 2:
        raise ValueError("--k_min must be >= 2")
    if args.k_max < args.k_min:
        raise ValueError("--k_max must be >= --k_min")
    if args.min_subset_size < 1:
        raise ValueError("--min_subset_size must be >= 1")
    if args.min_subset_size > len(feature_pool):
        raise ValueError("--min_subset_size cannot exceed number of features")
    if args.min_cluster_size < 1:
        raise ValueError("--min_cluster_size must be >= 1")

    k_values = list(range(args.k_min, args.k_max + 1))

    run_cfg = {
        "results_root": str(results_root),
        "out_root": str(out_root),
        "age_groups": age_groups,
        "features": feature_pool,
        "subset_mode": args.subset_mode,
        "min_subset_size": int(args.min_subset_size),
        "k_values": k_values,
        "min_units": int(args.min_units),
        "min_cluster_size": int(args.min_cluster_size),
        "log_fr": bool(not args.no_log_fr),
        "standardize": bool(not args.no_standardize),
        "random_state": int(args.random_state),
        "n_init": int(args.n_init),
        "covariance_type": args.covariance_type,
    }
    (out_root / "run_config.json").write_text(json.dumps(run_cfg, indent=2))

    for age_group in age_groups:
        evaluate_age_group(
            age_group,
            results_root=results_root,
            out_root=out_root,
            feature_pool=feature_pool,
            subset_mode=args.subset_mode,
            min_subset_size=args.min_subset_size,
            k_values=k_values,
            min_units=args.min_units,
            min_cluster_size=args.min_cluster_size,
            log_fr=not args.no_log_fr,
            standardize=not args.no_standardize,
            random_state=args.random_state,
            n_init=args.n_init,
            covariance_type=args.covariance_type,
        )

    print("Done.")


if __name__ == "__main__":
    main()
