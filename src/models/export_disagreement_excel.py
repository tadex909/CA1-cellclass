from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def parse_csv_list(raw: str | None) -> list[str]:
    if not raw:
        return []
    return [x.strip() for x in raw.split(",") if x.strip()]


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


def load_all_rows(comparison_root: Path) -> pd.DataFrame:
    all_path = comparison_root / "all_age_groups_gmm2_vs_type_u.parquet"
    if all_path.exists():
        return pd.read_parquet(all_path)

    parts: list[pd.DataFrame] = []
    for age_dir in sorted([p for p in comparison_root.iterdir() if p.is_dir()]):
        p = age_dir / f"{age_dir.name}_gmm2_vs_type_u.parquet"
        if not p.exists():
            continue
        d = pd.read_parquet(p)
        if "age_group" not in d.columns:
            d["age_group"] = age_dir.name
        parts.append(d)

    if not parts:
        raise FileNotFoundError(
            f"No *_gmm2_vs_type_u.parquet files found under {comparison_root}"
        )
    return pd.concat(parts, ignore_index=True)


def ensure_discrepancy_columns(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()

    if "type_u_type" not in d.columns:
        if "type_u_binary" not in d.columns and "allcel__type_u" in d.columns:
            d["type_u_binary"] = normalize_type_u(d["allcel__type_u"])
        if "type_u_binary" in d.columns:
            d["type_u_type"] = np.where(
                d["type_u_binary"] == 0, "interneuron", "pyramidal"
            )

    if "pred_type" not in d.columns and "pred_binary" in d.columns:
        d["pred_type"] = np.where(d["pred_binary"] == 0, "interneuron", "pyramidal")

    if "discrepancy" not in d.columns:
        if {"pred_binary", "type_u_binary"}.issubset(d.columns):
            valid = d["type_u_binary"].notna()
            d["discrepancy"] = np.where(
                valid, d["pred_binary"] != d["type_u_binary"], pd.NA
            )
        elif {"pred_type", "type_u_type"}.issubset(d.columns):
            valid = d["type_u_type"].notna() & d["pred_type"].notna()
            d["discrepancy"] = np.where(valid, d["pred_type"] != d["type_u_type"], pd.NA)
        else:
            raise ValueError(
                "Could not build discrepancy column. Missing comparison labels."
            )

    if "comparable" not in d.columns:
        if "type_u_binary" in d.columns:
            d["comparable"] = d["type_u_binary"].notna()
        else:
            d["comparable"] = d["type_u_type"].notna()

    return d


def add_posterior_probability(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()

    if "gmm_pmax" in d.columns:
        d["posterior_probability"] = d["gmm_pmax"].astype(float)
        return d

    p_int = d.get("gmm_p_interneuron")
    p_pyr = d.get("gmm_p_pyramidal")
    pred = d.get("pred_type")
    if p_int is not None and p_pyr is not None and pred is not None:
        d["posterior_probability"] = np.where(
            pred.astype(str).str.lower().eq("interneuron"),
            pd.to_numeric(p_int, errors="coerce"),
            pd.to_numeric(p_pyr, errors="coerce"),
        )
    else:
        d["posterior_probability"] = np.nan
    return d


def make_disagreement_sheet(
    d_cmp: pd.DataFrame,
    *,
    features_used: list[str],
) -> pd.DataFrame:
    d_dis = d_cmp[d_cmp["discrepancy"] == True].copy()  # noqa: E712
    d_dis = add_posterior_probability(d_dis)

    d_dis["vincas_classification"] = d_dis.get("type_u_type", pd.NA)
    d_dis["gmm_classification"] = d_dis.get("pred_type", pd.NA)
    d_dis["features_used_for_gmm"] = ", ".join(features_used)

    preferred_cols = [
        "age_group",
        "session_id",
        "cell_id",
        "unit_uid",
        "vincas_classification",
        "gmm_classification",
        "posterior_probability",
        "allcel__sm_u_any",
    ]
    keep_cols = [c for c in preferred_cols if c in d_dis.columns]
    d_dis = d_dis[keep_cols]
    d_dis = d_dis.sort_values(
        ["age_group", "posterior_probability"], ascending=[True, False]
    ).reset_index(drop=True)
    return d_dis


def make_overall_summary(
    d_cmp: pd.DataFrame,
    d_dis: pd.DataFrame,
    *,
    features_used: list[str],
    comparison_root: Path,
) -> pd.DataFrame:
    n_total = int(len(d_cmp))
    n_dis = int(len(d_dis))
    frac = float(n_dis / n_total) if n_total > 0 else float("nan")

    rows = [
        {"metric": "comparison_root", "value": str(comparison_root)},
        {"metric": "n_total_compared_neurons", "value": n_total},
        {"metric": "n_disagreeing_neurons", "value": n_dis},
        {"metric": "disagreement_fraction", "value": frac},
        {"metric": "features_used_for_gmm", "value": ", ".join(features_used)},
    ]
    return pd.DataFrame(rows)


def make_by_age_summary(d_cmp: pd.DataFrame) -> pd.DataFrame:
    by_age_total = (
        d_cmp.groupby("age_group", dropna=False)
        .size()
        .rename("n_total_compared")
        .reset_index()
    )
    by_age_dis = (
        d_cmp[d_cmp["discrepancy"] == True]  # noqa: E712
        .groupby("age_group", dropna=False)
        .size()
        .rename("n_disagreeing")
        .reset_index()
    )
    out = by_age_total.merge(by_age_dis, on="age_group", how="left")
    out["n_disagreeing"] = out["n_disagreeing"].fillna(0).astype(int)
    out["disagreement_fraction"] = out["n_disagreeing"] / out["n_total_compared"]
    return out.sort_values("age_group").reset_index(drop=True)


def make_classification_counts(d_cmp: pd.DataFrame) -> pd.DataFrame:
    d = d_cmp.copy()
    d["vincas_classification"] = d.get("type_u_type", pd.NA)
    d["gmm_classification"] = d.get("pred_type", pd.NA)

    grouped = (
        d.groupby(["vincas_classification", "gmm_classification"], dropna=False)
        .size()
        .rename("n_neurons")
        .reset_index()
        .sort_values("n_neurons", ascending=False)
        .reset_index(drop=True)
    )
    return grouped


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Export an Excel report of disagreements between Vinca's type_u and GMM labels."
        )
    )
    ap.add_argument(
        "--comparison_root",
        type=str,
        default="results/type_u_comparison_valero_feats_3",
    )
    ap.add_argument(
        "--out_xlsx",
        type=str,
        default="",
        help="Default: <comparison_root>/disagreement_report.xlsx",
    )
    ap.add_argument(
        "--age_groups",
        type=str,
        default="",
        help="Optional comma-separated filter, e.g. P16-18,P19-21",
    )
    args = ap.parse_args()

    comparison_root = Path(args.comparison_root)
    if not comparison_root.exists():
        raise FileNotFoundError(f"Missing comparison root: {comparison_root}")

    out_xlsx = (
        Path(args.out_xlsx)
        if args.out_xlsx
        else comparison_root / "disagreement_report.xlsx"
    )
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)

    run_config_path = comparison_root / "run_config.json"
    features_used: list[str] = []
    if run_config_path.exists():
        with run_config_path.open("r", encoding="utf-8") as f:
            cfg = json.load(f)
        features_used = [str(x) for x in cfg.get("features", [])]

    d = load_all_rows(comparison_root)
    d = ensure_discrepancy_columns(d)

    ages = set(parse_csv_list(args.age_groups))
    if ages:
        d = d[d["age_group"].isin(ages)].copy()

    d_cmp = d[d["comparable"] == True].copy()  # noqa: E712
    d_dis = make_disagreement_sheet(d_cmp, features_used=features_used)
    overall = make_overall_summary(
        d_cmp,
        d_dis,
        features_used=features_used,
        comparison_root=comparison_root,
    )
    by_age = make_by_age_summary(d_cmp)
    counts = make_classification_counts(d_cmp)
    features_df = pd.DataFrame(
        {"feature_index": range(1, len(features_used) + 1), "feature": features_used}
    )

    with pd.ExcelWriter(out_xlsx) as writer:
        d_dis.to_excel(writer, sheet_name="disagreements", index=False)
        overall.to_excel(writer, sheet_name="overall_summary", index=False)
        by_age.to_excel(writer, sheet_name="by_age_group", index=False)
        counts.to_excel(writer, sheet_name="classification_counts", index=False)
        features_df.to_excel(writer, sheet_name="features_used", index=False)

    print(f"Wrote {out_xlsx}")
    print(f"Compared neurons: {len(d_cmp)}")
    print(f"Disagreeing neurons: {len(d_dis)}")


if __name__ == "__main__":
    main()
