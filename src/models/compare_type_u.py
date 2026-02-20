from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_samples

from fitting import DEFAULT_AGE_GROUPS, DEFAULT_FEATURES, evaluate_gmm, prepare_matrix


def normalize_type_u(series: pd.Series) -> pd.Series:
    """
    Map external labels to binary convention:
      0 = interneuron
      1 = pyramidal
    """
    s = series.copy()
    if pd.api.types.is_numeric_dtype(s):
        out = pd.to_numeric(s, errors="coerce")
        out = out.where(out.isin([0, 1]), np.nan)
        return out.astype("Int64")

    s = s.astype(str).str.strip().str.lower()
    mapping = {
        "0": 0,
        "interneuron": 0,
        "int": 0,
        "1": 1,
        "pyramidal": 1,
        "pyr": 1,
    }
    out = s.map(mapping)
    return out.astype("Int64")


def parse_csv_list(raw: str | None) -> list[str]:
    if not raw:
        return []
    return [x.strip() for x in raw.split(",") if x.strip()]


def compute_silhouette_samples_safe(X: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """
    Per-sample silhouette values, or NaN when undefined.
    """
    out = np.full(X.shape[0], np.nan, dtype=np.float64)
    uniq = np.unique(labels)
    if uniq.size < 2:
        return out
    if X.shape[0] <= uniq.size:
        return out
    try:
        out = silhouette_samples(X, labels)
    except Exception:
        return out
    return out


def compare_one_age_group(
    age_group: str,
    *,
    results_root: Path,
    out_root: Path,
    features: list[str],
    log_fr: bool,
    standardize: bool,
    random_state: int,
    n_init: int,
    covariance_type: str,
    min_units: int,
) -> pd.DataFrame:
    in_path = results_root / age_group / f"{age_group}_clean_units.parquet"
    if not in_path.exists():
        print(f"{age_group}: missing {in_path}, skipping")
        return pd.DataFrame()

    df = pd.read_parquet(in_path)
    if "allcel__type_u" not in df.columns:
        print(f"{age_group}: no allcel__type_u column, skipping")
        return pd.DataFrame()

    if "unit_uid" not in df.columns:
        df["unit_uid"] = df["session_id"].astype(str) + "__cell" + df["cell_id"].astype(str)

    pack = prepare_matrix(df, features, log_fr=log_fr, standardize=standardize)
    if pack.rows_after < min_units:
        print(f"{age_group}: too few finite units ({pack.rows_after} < {min_units}), skipping")
        return pd.DataFrame()

    ev = evaluate_gmm(
        pack.X,
        2,
        random_state=random_state,
        n_init=n_init,
        covariance_type=covariance_type,
    )
    labels = ev["labels"]
    proba = ev["proba"]
    sil_samples = compute_silhouette_samples_safe(pack.X, labels)

    d = pack.df_valid.copy()
    d["type_u_binary"] = normalize_type_u(d["allcel__type_u"])

    # Define cluster->cell type using spike duration:
    # lower mean spk_duration_ms => interneuron (0), other => pyramidal (1)
    if "spk_duration_ms" not in d.columns:
        raise ValueError("spk_duration_ms is required to map clusters to interneuron/pyramidal.")

    cstats = (
        pd.DataFrame({"cluster": labels, "dur": d["spk_duration_ms"].to_numpy()})
        .groupby("cluster")["dur"]
        .mean()
    )
    interneuron_cluster = int(cstats.idxmin())
    pyramidal_cluster = 1 - interneuron_cluster

    d["gmm_cluster"] = labels
    d["gmm_p_interneuron"] = proba[:, interneuron_cluster]
    d["gmm_p_pyramidal"] = proba[:, pyramidal_cluster]
    d["gmm_pmax"] = proba.max(axis=1)
    d["gmm_margin"] = np.abs(proba[:, interneuron_cluster] - proba[:, pyramidal_cluster])
    d["gmm_silhouette_sample"] = sil_samples

    d["pred_binary"] = np.where(labels == interneuron_cluster, 0, 1).astype("int64")
    d["pred_type"] = np.where(d["pred_binary"] == 0, "interneuron", "pyramidal")
    d["type_u_type"] = np.where(d["type_u_binary"] == 0, "interneuron", "pyramidal")

    valid_cmp = d["type_u_binary"].notna()
    d["comparable"] = valid_cmp
    d["match"] = np.where(valid_cmp, d["pred_binary"] == d["type_u_binary"], pd.NA)
    d["discrepancy"] = np.where(valid_cmp, ~(d["pred_binary"] == d["type_u_binary"]), pd.NA)

    out_dir = out_root / age_group
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{age_group}_gmm2_vs_type_u.parquet"
    d.to_parquet(out_path, index=False)

    n_cmp = int(valid_cmp.sum())
    n_match = int(((d["match"] == True) & valid_cmp).sum())
    rate = float(n_match / n_cmp) if n_cmp > 0 else float("nan")
    sil_mean = float(np.nanmean(d["gmm_silhouette_sample"])) if len(d) > 0 else float("nan")
    print(
        f"{age_group}: compared={n_cmp}, match={n_match}, "
        f"match_rate={rate:.3f}, silhouette_mean={sil_mean:.3f}"
    )
    return d


def build_summaries(df_all: pd.DataFrame, out_root: Path) -> None:
    if df_all.empty:
        (out_root / "summary.json").write_text(json.dumps({"status": "no_data"}, indent=2))
        return

    d = df_all.copy()
    d_cmp = d[d["comparable"] == True].copy()

    by_age = (
        d_cmp.groupby("age_group", dropna=False)
        .agg(
            n_compared=("comparable", "size"),
            n_match=("match", lambda x: int((x == True).sum())),
            n_discrepancy=("discrepancy", lambda x: int((x == True).sum())),
        )
        .reset_index()
    )
    by_age["match_rate"] = by_age["n_match"] / by_age["n_compared"]
    by_age["discrepancy_rate"] = by_age["n_discrepancy"] / by_age["n_compared"]
    by_age.to_csv(out_root / "discrepancies_by_age.csv", index=False)

    by_session = (
        d_cmp.groupby(["age_group", "session_id"], dropna=False)
        .agg(
            n_compared=("comparable", "size"),
            n_match=("match", lambda x: int((x == True).sum())),
            n_discrepancy=("discrepancy", lambda x: int((x == True).sum())),
            mean_confidence=("gmm_pmax", "mean"),
        )
        .reset_index()
    )
    by_session["match_rate"] = by_session["n_match"] / by_session["n_compared"]
    by_session["discrepancy_rate"] = by_session["n_discrepancy"] / by_session["n_compared"]
    by_session = by_session.sort_values(
        ["n_discrepancy", "discrepancy_rate"], ascending=[False, False]
    ).reset_index(drop=True)
    by_session.to_csv(out_root / "discrepancies_by_session.csv", index=False)

    by_cluster = (
        d.groupby(["age_group", "gmm_cluster"], dropna=False)
        .agg(
            n_units=("gmm_cluster", "size"),
            silhouette_mean=("gmm_silhouette_sample", "mean"),
            silhouette_median=("gmm_silhouette_sample", "median"),
            silhouette_min=("gmm_silhouette_sample", "min"),
            silhouette_max=("gmm_silhouette_sample", "max"),
        )
        .reset_index()
    )
    by_cluster.to_csv(out_root / "silhouette_by_age_cluster.csv", index=False)

    discrepant_neurons = d_cmp[d_cmp["discrepancy"] == True].copy()
    discrepant_neurons = discrepant_neurons.sort_values(
        ["gmm_pmax", "gmm_margin"], ascending=[False, False]
    ).reset_index(drop=True)
    keep_cols = [
        "age_group",
        "session_id",
        "cell_id",
        "unit_uid",
        "allcel__type_u",
        "type_u_binary",
        "type_u_type",
        "pred_binary",
        "pred_type",
        "gmm_cluster",
        "gmm_p_interneuron",
        "gmm_p_pyramidal",
        "gmm_pmax",
        "gmm_margin",
    ]
    keep_cols = [c for c in keep_cols if c in discrepant_neurons.columns]
    discrepant_neurons[keep_cols].to_parquet(out_root / "discrepant_neurons.parquet", index=False)
    discrepant_neurons[keep_cols].to_csv(
        out_root / "discrepant_neurons_all.csv", index=False
    )
    discrepant_neurons[keep_cols].head(500).to_csv(
        out_root / "discrepant_neurons_top500.csv", index=False
    )

    n_compared = int(len(d_cmp))
    n_match = int((d_cmp["match"] == True).sum())
    n_discrepancy = int((d_cmp["discrepancy"] == True).sum())
    overall_match_rate = float(n_match / n_compared) if n_compared > 0 else float("nan")
    overall_silhouette_mean = float(np.nanmean(d["gmm_silhouette_sample"]))
    overall_silhouette_median = float(np.nanmedian(d["gmm_silhouette_sample"]))

    worst_age = (
        by_age.sort_values(["n_discrepancy", "discrepancy_rate"], ascending=[False, False])
        .head(1)
        .to_dict(orient="records")
    )
    worst_session = by_session.head(1).to_dict(orient="records")

    summary: dict[str, Any] = {
        "status": "ok",
        "n_total_units_after_finite_filter": int(len(d)),
        "n_comparable_units": n_compared,
        "n_match": n_match,
        "n_discrepancy": n_discrepancy,
        "overall_match_rate": overall_match_rate,
        "overall_silhouette_mean": overall_silhouette_mean,
        "overall_silhouette_median": overall_silhouette_median,
        "age_group_with_most_discrepancies": worst_age[0] if worst_age else None,
        "session_with_most_discrepancies": worst_session[0] if worst_session else None,
        "outputs": {
            "discrepancies_by_age_csv": str(out_root / "discrepancies_by_age.csv"),
            "discrepancies_by_session_csv": str(out_root / "discrepancies_by_session.csv"),
            "silhouette_by_age_cluster_csv": str(out_root / "silhouette_by_age_cluster.csv"),
            "discrepant_neurons_parquet": str(out_root / "discrepant_neurons.parquet"),
            "discrepant_neurons_all_csv": str(out_root / "discrepant_neurons_all.csv"),
            "discrepant_neurons_top500_csv": str(out_root / "discrepant_neurons_top500.csv"),
        },
    }
    (out_root / "summary.json").write_text(json.dumps(summary, indent=2))


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Run fixed 2-cluster GMM (all features) for each age group and compare with allcel__type_u."
        )
    )
    ap.add_argument("--results_root", type=str, default="results")
    ap.add_argument("--out_root", type=str, default="results/type_u_comparison")
    ap.add_argument("--age_groups", type=str, default=",".join(DEFAULT_AGE_GROUPS))
    ap.add_argument("--features", type=str, default=",".join(DEFAULT_FEATURES))
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
    ap.add_argument("--min_units", type=int, default=10)
    args = ap.parse_args()

    results_root = Path(args.results_root)
    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    age_groups = parse_csv_list(args.age_groups) or DEFAULT_AGE_GROUPS
    features = parse_csv_list(args.features) or DEFAULT_FEATURES

    cfg = {
        "results_root": str(results_root),
        "out_root": str(out_root),
        "age_groups": age_groups,
        "features": features,
        "log_fr": bool(not args.no_log_fr),
        "standardize": bool(not args.no_standardize),
        "k": 2,
        "random_state": int(args.random_state),
        "n_init": int(args.n_init),
        "covariance_type": args.covariance_type,
        "min_units": int(args.min_units),
    }
    (out_root / "run_config.json").write_text(json.dumps(cfg, indent=2))

    parts: list[pd.DataFrame] = []
    for age_group in age_groups:
        d = compare_one_age_group(
            age_group,
            results_root=results_root,
            out_root=out_root,
            features=features,
            log_fr=not args.no_log_fr,
            standardize=not args.no_standardize,
            random_state=args.random_state,
            n_init=args.n_init,
            covariance_type=args.covariance_type,
            min_units=args.min_units,
        )
        if not d.empty:
            parts.append(d)

    df_all = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    if not df_all.empty:
        df_all.to_parquet(out_root / "all_age_groups_gmm2_vs_type_u.parquet", index=False)
    build_summaries(df_all, out_root)
    print("Done.")


if __name__ == "__main__":
    main()
