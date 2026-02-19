from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

from fitting import DEFAULT_FEATURES, evaluate_gmm, prepare_matrix


def parse_csv_list(raw: str | None) -> list[str]:
    if not raw:
        return []
    return [x.strip() for x in raw.split(",") if x.strip()]


def parse_int_list(raw: str | None, default: list[int]) -> list[int]:
    if not raw:
        return default
    return [int(x.strip()) for x in raw.split(",") if x.strip()]


def normalize_type_u(series: pd.Series) -> pd.Series:
    """
    Standardize type_u to:
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


def iter_feature_subsets(features: list[str], mode: str, min_subset_size: int) -> list[list[str]]:
    if mode == "full_only":
        return [features]

    out: list[list[str]] = [features]
    if mode == "leave_one_out":
        if len(features) > 1:
            out.extend([list(x) for x in combinations(features, len(features) - 1)])
        return out

    # mode == all
    for size in range(len(features) - 1, min_subset_size - 1, -1):
        out.extend([list(x) for x in combinations(features, size)])
    return out


def evaluate_run(
    df: pd.DataFrame,
    *,
    subset: list[str],
    seed: int,
    log_fr: bool,
    standardize: bool,
    n_init: int,
    covariance_type: str,
    min_units: int,
) -> tuple[dict[str, Any], pd.DataFrame]:
    pack = prepare_matrix(df, subset, log_fr=log_fr, standardize=standardize)
    if pack.rows_after < min_units:
        return (
            {
                "status": "too_few_units",
                "error": f"{pack.rows_after} < min_units {min_units}",
                "n_rows_after": pack.rows_after,
                "used_cols": "+".join(pack.used_cols),
            },
            pd.DataFrame(),
        )

    ev = evaluate_gmm(
        pack.X,
        2,
        random_state=seed,
        n_init=n_init,
        covariance_type=covariance_type,
    )

    labels = ev["labels"]
    proba = ev["proba"]
    d = pack.df_valid.copy()

    if "spk_duration_ms" not in d.columns:
        return (
            {
                "status": "error",
                "error": "spk_duration_ms missing; cannot map clusters to classes",
                "n_rows_after": pack.rows_after,
                "used_cols": "+".join(pack.used_cols),
            },
            pd.DataFrame(),
        )

    cstats = (
        pd.DataFrame({"cluster": labels, "dur": d["spk_duration_ms"].to_numpy()})
        .groupby("cluster")["dur"]
        .mean()
    )
    interneuron_cluster = int(cstats.idxmin())
    pyramidal_cluster = 1 - interneuron_cluster

    d["type_u_binary"] = normalize_type_u(d["allcel__type_u"]) if "allcel__type_u" in d.columns else pd.Series([pd.NA] * len(d), index=d.index, dtype="Int64")
    d["pred_binary"] = np.where(labels == interneuron_cluster, 0, 1).astype("int64")
    d["gmm_cluster"] = labels.astype("int64")
    d["gmm_p_interneuron"] = proba[:, interneuron_cluster]
    d["gmm_p_pyramidal"] = proba[:, pyramidal_cluster]
    d["gmm_pmax"] = proba.max(axis=1)
    d["gmm_margin"] = np.abs(proba[:, interneuron_cluster] - proba[:, pyramidal_cluster])

    valid = d["type_u_binary"].notna()
    n_cmp = int(valid.sum())
    acc_raw = float((d.loc[valid, "pred_binary"] == d.loc[valid, "type_u_binary"]).mean()) if n_cmp > 0 else float("nan")
    acc_flip = float(((1 - d.loc[valid, "pred_binary"]) == d.loc[valid, "type_u_binary"]).mean()) if n_cmp > 0 else float("nan")
    acc_best = float(max(acc_raw, acc_flip)) if n_cmp > 0 else float("nan")
    ari = (
        float(adjusted_rand_score(d.loc[valid, "type_u_binary"].astype(int), d.loc[valid, "pred_binary"].astype(int)))
        if n_cmp > 0
        else float("nan")
    )

    metrics = {
        "status": "ok",
        "error": "",
        "n_rows_after": pack.rows_after,
        "used_cols": "+".join(pack.used_cols),
        "bic": ev["bic"],
        "aic": ev["aic"],
        "silhouette": ev["silhouette"],
        "davies_bouldin": ev["davies_bouldin"],
        "calinski_harabasz": ev["calinski_harabasz"],
        "min_cluster_size": ev["min_cluster_size"],
        "max_cluster_size": ev["max_cluster_size"],
        "n_compared_type_u": n_cmp,
        "acc_raw": acc_raw,
        "acc_flip": acc_flip,
        "acc_best_label_permutation": acc_best,
        "ari_vs_type_u": ari,
    }
    return metrics, d


def build_pairwise_ari(pred_df: pd.DataFrame) -> pd.DataFrame:
    runs = pred_df["run_id"].drop_duplicates().tolist()
    per_run = {rid: pred_df[pred_df["run_id"] == rid][["unit_uid", "pred_binary"]] for rid in runs}
    rows: list[dict[str, Any]] = []
    for i, r1 in enumerate(runs):
        for r2 in runs[i + 1 :]:
            a = per_run[r1].rename(columns={"pred_binary": "p1"})
            b = per_run[r2].rename(columns={"pred_binary": "p2"})
            m = a.merge(b, on="unit_uid", how="inner")
            n = len(m)
            if n < 2:
                ari = float("nan")
            else:
                ari = float(adjusted_rand_score(m["p1"].to_numpy(), m["p2"].to_numpy()))
            rows.append({"run_1": r1, "run_2": r2, "n_overlap_units": n, "ari": ari})
    return pd.DataFrame(rows)


def build_unit_uncertainty(pred_df: pd.DataFrame, top_n: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    d = pred_df.copy()
    # instability in [0,1], max at split votes 50/50
    unit = (
        d.groupby(["unit_uid", "age_group", "session_id", "cell_id"], dropna=False)
        .agg(
            n_runs=("pred_binary", "size"),
            pred_fraction_pyramidal=("pred_binary", "mean"),
            mean_pmax=("gmm_pmax", "mean"),
            mean_margin=("gmm_margin", "mean"),
            std_p_pyramidal=("gmm_p_pyramidal", "std"),
            type_u_binary=("type_u_binary", "first"),
        )
        .reset_index()
    )
    unit["std_p_pyramidal"] = unit["std_p_pyramidal"].fillna(0.0)
    unit["instability_vote"] = 1.0 - np.abs(2.0 * unit["pred_fraction_pyramidal"] - 1.0)
    unit["posterior_uncertainty"] = 2.0 * (1.0 - unit["mean_pmax"])  # maps [0.5,1]->[1,0]
    unit["posterior_uncertainty"] = unit["posterior_uncertainty"].clip(lower=0.0, upper=1.0)
    unit["uncertainty_score"] = 0.5 * unit["instability_vote"] + 0.5 * unit["posterior_uncertainty"]

    top = unit.sort_values(
        ["uncertainty_score", "instability_vote", "posterior_uncertainty"],
        ascending=[False, False, False],
    ).head(top_n)
    return unit, top


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Stability analysis for one age group: fixed k=2 GMM, multiple feature subsets and seeds."
        )
    )
    ap.add_argument("--results_root", type=str, default="results")
    ap.add_argument("--out_root", type=str, default="results/stability")
    ap.add_argument("--age_group", type=str, required=True)
    ap.add_argument("--features", type=str, default=",".join(DEFAULT_FEATURES))
    ap.add_argument(
        "--subset_mode",
        type=str,
        choices=["full_only", "leave_one_out", "all"],
        default="leave_one_out",
    )
    ap.add_argument("--min_subset_size", type=int, default=2)
    ap.add_argument("--seeds", type=str, default="0,1,2,3,4,5,6,7,8,9")
    ap.add_argument("--n_init", type=int, default=10)
    ap.add_argument(
        "--covariance_type",
        type=str,
        choices=["full", "tied", "diag", "spherical"],
        default="full",
    )
    ap.add_argument("--no_log_fr", action="store_true")
    ap.add_argument("--no_standardize", action="store_true")
    ap.add_argument("--min_units", type=int, default=20)
    ap.add_argument("--top_n_uncertain", type=int, default=100)
    args = ap.parse_args()

    results_root = Path(args.results_root)
    in_path = results_root / args.age_group / f"{args.age_group}_clean_units.parquet"
    if not in_path.exists():
        raise FileNotFoundError(f"Missing file: {in_path}")

    df = pd.read_parquet(in_path)
    if "unit_uid" not in df.columns:
        df["unit_uid"] = df["session_id"].astype(str) + "__cell" + df["cell_id"].astype(str)

    features = parse_csv_list(args.features) or DEFAULT_FEATURES
    seeds = parse_int_list(args.seeds, default=[0])
    subsets = iter_feature_subsets(features, args.subset_mode, min_subset_size=args.min_subset_size)

    out_dir = Path(args.out_root) / args.age_group
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg = {
        "age_group": args.age_group,
        "input": str(in_path),
        "features": features,
        "subset_mode": args.subset_mode,
        "min_subset_size": args.min_subset_size,
        "seeds": seeds,
        "k": 2,
        "n_init": args.n_init,
        "covariance_type": args.covariance_type,
        "log_fr": bool(not args.no_log_fr),
        "standardize": bool(not args.no_standardize),
        "min_units": args.min_units,
        "top_n_uncertain": args.top_n_uncertain,
    }
    (out_dir / "run_config.json").write_text(json.dumps(cfg, indent=2))

    metrics_rows: list[dict[str, Any]] = []
    pred_rows: list[pd.DataFrame] = []
    run_idx = 0

    for subset in subsets:
        subset_name = "+".join(subset)
        for seed in seeds:
            run_idx += 1
            run_id = f"run_{run_idx:04d}"
            m, d_run = evaluate_run(
                df,
                subset=subset,
                seed=seed,
                log_fr=not args.no_log_fr,
                standardize=not args.no_standardize,
                n_init=args.n_init,
                covariance_type=args.covariance_type,
                min_units=args.min_units,
            )
            m.update(
                {
                    "run_id": run_id,
                    "subset_name": subset_name,
                    "seed": int(seed),
                    "n_features_requested": int(len(subset)),
                }
            )
            metrics_rows.append(m)

            if m["status"] == "ok" and not d_run.empty:
                keep_cols = [
                    "unit_uid",
                    "session_id",
                    "cell_id",
                    "type_u_binary",
                    "pred_binary",
                    "gmm_cluster",
                    "gmm_pmax",
                    "gmm_margin",
                    "gmm_p_pyramidal",
                ]
                keep_cols = [c for c in keep_cols if c in d_run.columns]
                tmp = d_run[keep_cols].copy()
                tmp["run_id"] = run_id
                tmp["subset_name"] = subset_name
                tmp["seed"] = int(seed)
                tmp["age_group"] = args.age_group
                pred_rows.append(tmp)

    metrics_df = pd.DataFrame(metrics_rows)
    metrics_df.to_csv(out_dir / "run_metrics.csv", index=False)

    ok_df = metrics_df[metrics_df["status"] == "ok"].copy()
    if ok_df.empty:
        summary = {"status": "no_valid_runs", "n_total_runs": int(len(metrics_df)), "n_ok_runs": 0}
        (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
        print("No valid runs.")
        return

    pred_df = pd.concat(pred_rows, ignore_index=True) if pred_rows else pd.DataFrame()
    pred_df.to_parquet(out_dir / "run_predictions.parquet", index=False)

    pairwise = build_pairwise_ari(pred_df)
    pairwise.to_csv(out_dir / "pairwise_run_ari.csv", index=False)

    unit_unc, top_unc = build_unit_uncertainty(pred_df, top_n=args.top_n_uncertain)
    unit_unc.to_csv(out_dir / "unit_uncertainty.csv", index=False)
    top_unc.to_csv(out_dir / "top_uncertain_neurons.csv", index=False)

    session_stability = (
        pred_df.groupby("session_id", dropna=False)
        .agg(
            n_units=("unit_uid", lambda x: int(pd.Series(x).nunique())),
            n_rows=("unit_uid", "size"),
            mean_pmax=("gmm_pmax", "mean"),
            mean_margin=("gmm_margin", "mean"),
            mean_pred_pyramidal=("pred_binary", "mean"),
        )
        .reset_index()
    )
    session_stability.to_csv(out_dir / "session_stability.csv", index=False)

    summary = {
        "status": "ok",
        "age_group": args.age_group,
        "n_total_runs": int(len(metrics_df)),
        "n_ok_runs": int(len(ok_df)),
        "best_run_by_bic": (
            ok_df.sort_values(["bic", "silhouette", "davies_bouldin"], ascending=[True, False, True])
            .head(1)
            .to_dict(orient="records")[0]
        ),
        "agreement_vs_type_u": {
            "acc_raw_mean": float(ok_df["acc_raw"].mean(skipna=True)),
            "acc_best_label_permutation_mean": float(ok_df["acc_best_label_permutation"].mean(skipna=True)),
            "ari_mean": float(ok_df["ari_vs_type_u"].mean(skipna=True)),
        },
        "run_to_run_stability": {
            "pairwise_ari_mean": float(pairwise["ari"].mean(skipna=True)),
            "pairwise_ari_median": float(pairwise["ari"].median(skipna=True)),
            "pairwise_ari_min": float(pairwise["ari"].min(skipna=True)),
            "pairwise_ari_max": float(pairwise["ari"].max(skipna=True)),
        },
        "uncertainty_definition": {
            "instability_vote": "1 - abs(2*f_pyramidal - 1), where f_pyramidal is run-wise vote fraction",
            "posterior_uncertainty": "2*(1-mean_pmax) for binary GMM",
            "uncertainty_score": "0.5*instability_vote + 0.5*posterior_uncertainty",
        },
        "outputs": {
            "run_metrics_csv": str(out_dir / "run_metrics.csv"),
            "run_predictions_parquet": str(out_dir / "run_predictions.parquet"),
            "pairwise_run_ari_csv": str(out_dir / "pairwise_run_ari.csv"),
            "unit_uncertainty_csv": str(out_dir / "unit_uncertainty.csv"),
            "top_uncertain_neurons_csv": str(out_dir / "top_uncertain_neurons.csv"),
            "session_stability_csv": str(out_dir / "session_stability.csv"),
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print("Done.")


if __name__ == "__main__":
    main()

