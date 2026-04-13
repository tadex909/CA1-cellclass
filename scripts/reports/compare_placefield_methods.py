from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

KEY_COLS = ["session_id", "cell_id", "condition_1b"]


def normalize_session_key(value: object) -> str:
    return " ".join(str(value).strip().split())


def sanitize_token(value: str) -> str:
    out = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip())
    return out.strip("_") or "x"


def format_float_token(name: str, value: float) -> str:
    txt = f"{float(value):g}".replace("-", "m").replace(".", "p")
    return f"{name}_{txt}"


def resolve_classification_csv(path_like: str | Path) -> Path:
    p = Path(path_like)
    if p.is_file():
        return p
    cand = p / "ssi_classification.csv"
    if cand.exists():
        return cand
    raise FileNotFoundError(
        f"Could not resolve classification CSV from {p}. "
        f"Expected file or directory containing {cand.name}."
    )


def infer_name(path_like: str | Path) -> str:
    p = Path(path_like)
    if p.is_file():
        return p.parent.name or p.stem
    return p.name or "method"


def load_classification(csv_path: Path) -> pd.DataFrame:
    d = pd.read_csv(csv_path)
    required = ["session_id", "cell_id", "condition_1b", "p_value"]
    missing = [c for c in required if c not in d.columns]
    if missing:
        raise KeyError(f"Missing required columns in {csv_path}: {missing}")

    d = d.copy()
    d["session_id"] = d["session_id"].astype(str)
    d["cell_id"] = pd.to_numeric(d["cell_id"], errors="coerce").astype("Int64")
    d["condition_1b"] = pd.to_numeric(d["condition_1b"], errors="coerce").astype("Int64")
    d["p_value"] = pd.to_numeric(d["p_value"], errors="coerce")
    if "ssi_obs" in d.columns:
        d["ssi_obs"] = pd.to_numeric(d["ssi_obs"], errors="coerce")
    else:
        d["ssi_obs"] = np.nan

    d = d.dropna(subset=["session_id", "cell_id", "condition_1b"]).copy()
    d["cell_id"] = d["cell_id"].astype(np.int64)
    d["condition_1b"] = d["condition_1b"].astype(np.int64)

    keep = ["session_id", "cell_id", "condition_1b", "p_value", "ssi_obs"]
    for col in ["SM", "sm_alpha", "source"]:
        if col in d.columns:
            keep.append(col)
    d = (
        d[keep]
        .sort_values(KEY_COLS, kind="stable")
        .drop_duplicates(subset=KEY_COLS, keep="first")
        .reset_index(drop=True)
    )
    return d


def add_sm_columns(d: pd.DataFrame, alpha: float) -> pd.DataFrame:
    out = d.copy()
    out["is_sm_alpha"] = np.isfinite(out["p_value"]) & (out["p_value"] < float(alpha))
    return out


def confusion_counts(sm_a: np.ndarray, sm_b: np.ndarray) -> dict[str, int]:
    a = np.asarray(sm_a, dtype=bool)
    b = np.asarray(sm_b, dtype=bool)
    return {
        "both_not_sm": int((~a & ~b).sum()),
        "a_only_sm": int((a & ~b).sum()),
        "b_only_sm": int((~a & b).sum()),
        "both_sm": int((a & b).sum()),
    }


def load_age_map(schedule_xlsx: Path, sheet_name: str) -> dict[str, int]:
    if not schedule_xlsx.exists():
        return {}
    d = pd.read_excel(schedule_xlsx, sheet_name=sheet_name)
    required = {"SessionName", "Age"}
    if not required.issubset(d.columns):
        return {}
    d = d[["SessionName", "Age"]].copy()
    d["session_key"] = d["SessionName"].map(normalize_session_key)
    d["Age"] = pd.to_numeric(d["Age"], errors="coerce")
    d = d.dropna(subset=["session_key", "Age"])
    d["Age"] = d["Age"].astype(np.int64)
    d = d.drop_duplicates(subset=["session_key"], keep="first")
    return {str(k): int(v) for k, v in zip(d["session_key"], d["Age"])}


def plot_confusion_matrix(
    conf: dict[str, int],
    label_a: str,
    label_b: str,
    title: str,
    out_png: Path,
) -> None:
    mat = np.array(
        [
            [conf["both_not_sm"], conf["b_only_sm"]],
            [conf["a_only_sm"], conf["both_sm"]],
        ],
        dtype=np.int64,
    )
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    im = ax.imshow(mat, cmap="Blues")
    ax.set_xticks([0, 1], labels=[f"{label_b}: not SM", f"{label_b}: SM"])
    ax.set_yticks([0, 1], labels=[f"{label_a}: not SM", f"{label_a}: SM"])
    ax.set_title(title)
    for i in range(2):
        for j in range(2):
            ax.text(j, i, str(int(mat[i, j])), ha="center", va="center", color="black")
    fig.colorbar(im, ax=ax, fraction=0.05, pad=0.04)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_pvalue_scatter(merged: pd.DataFrame, label_a: str, label_b: str, out_png: Path) -> None:
    d = merged.copy()
    d = d[np.isfinite(d["p_a"]) & np.isfinite(d["p_b"])].copy()
    if d.empty:
        return
    min_pos = np.nanmin(np.concatenate([d["p_a"].to_numpy(), d["p_b"].to_numpy()]))
    eps = float(max(min_pos / 2.0, 1e-6)) if np.isfinite(min_pos) and min_pos > 0 else 1e-6
    x = np.clip(d["p_a"].to_numpy(dtype=np.float64), eps, 1.0)
    y = np.clip(d["p_b"].to_numpy(dtype=np.float64), eps, 1.0)
    disc = d["discordant_sm"].to_numpy(dtype=bool)

    fig, ax = plt.subplots(figsize=(6.0, 5.8))
    ax.scatter(x[~disc], y[~disc], s=8, alpha=0.35, label="agree")
    ax.scatter(x[disc], y[disc], s=12, alpha=0.75, label="disagree")
    ax.plot([eps, 1.0], [eps, 1.0], "k--", linewidth=1.0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(f"{label_a} p-value")
    ax.set_ylabel(f"{label_b} p-value")
    ax.set_title("P-value Comparison (Row-level)")
    ax.grid(alpha=0.25, which="both")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_delta_hist(merged: pd.DataFrame, out_png: Path) -> None:
    d = merged.copy()
    d = d[np.isfinite(d["delta_p"])].copy()
    if d.empty:
        return
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.hist(d["delta_p"].to_numpy(dtype=np.float64), bins=80, color="#9ecae1", edgecolor="#4a6fa5")
    ax.axvline(0.0, color="k", linestyle="--", linewidth=1.0)
    ax.set_xlabel("delta_p = p_b - p_a")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of P-value Differences")
    ax.grid(alpha=0.2, axis="y")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_age_disagreement(age_summary: pd.DataFrame, out_png: Path) -> None:
    if age_summary.empty:
        return
    d = age_summary.sort_values("age").copy()
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.bar(d["age"].astype(str), d["discordant_rate_row"], color="#74c476", edgecolor="#2b8cbe")
    ax.set_xlabel("Age")
    ax.set_ylabel("Row-level disagreement rate")
    ax.set_title("Disagreement Rate by Age")
    ax.set_ylim(0.0, max(0.05, float(d["discordant_rate_row"].max()) * 1.2))
    ax.grid(alpha=0.2, axis="y")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Compare two placefield SSI classification datasets (e.g., circular-shift vs random). "
            "Outputs overlap/disagreement tables plus summary plots."
        )
    )
    ap.add_argument("--a", type=str, required=True, help="Method A: classification CSV or its parent folder.")
    ap.add_argument("--b", type=str, required=True, help="Method B: classification CSV or its parent folder.")
    ap.add_argument("--name_a", type=str, default="", help="Display name for method A.")
    ap.add_argument("--name_b", type=str, default="", help="Display name for method B.")
    ap.add_argument("--alpha", type=float, default=0.05, help="SM threshold: p_value < alpha.")
    ap.add_argument("--out_root", type=str, default="results/tables/placefield_method_compare")
    ap.add_argument("--schedule_xlsx", type=str, default="data/schedule.xlsx")
    ap.add_argument("--schedule_sheet", type=str, default="VINCA")
    return ap


def main() -> None:
    args = build_parser().parse_args()
    if not (0.0 < float(args.alpha) < 1.0):
        raise ValueError("--alpha must be strictly between 0 and 1.")

    csv_a = resolve_classification_csv(args.a)
    csv_b = resolve_classification_csv(args.b)
    name_a = args.name_a.strip() or infer_name(args.a)
    name_b = args.name_b.strip() or infer_name(args.b)

    out_root = Path(args.out_root)
    run_tag = "__".join(
        [
            f"{sanitize_token(name_a)}_vs_{sanitize_token(name_b)}",
            format_float_token("alpha", float(args.alpha)),
        ]
    )
    out_dir = out_root / run_tag
    out_dir.mkdir(parents=True, exist_ok=True)

    age_map = load_age_map(Path(args.schedule_xlsx), str(args.schedule_sheet))

    a = add_sm_columns(load_classification(csv_a), float(args.alpha))
    b = add_sm_columns(load_classification(csv_b), float(args.alpha))

    a2 = a.rename(
        columns={
            "p_value": "p_a",
            "ssi_obs": "ssi_obs_a",
            "is_sm_alpha": "sm_a",
            "SM": "sm_src_a",
            "sm_alpha": "sm_alpha_src_a",
            "source": "source_a",
        }
    )
    b2 = b.rename(
        columns={
            "p_value": "p_b",
            "ssi_obs": "ssi_obs_b",
            "is_sm_alpha": "sm_b",
            "SM": "sm_src_b",
            "sm_alpha": "sm_alpha_src_b",
            "source": "source_b",
        }
    )

    merged = a2.merge(b2, on=KEY_COLS, how="inner")
    merged["delta_p"] = merged["p_b"] - merged["p_a"]
    merged["abs_delta_p"] = np.abs(merged["delta_p"])
    merged["discordant_sm"] = merged["sm_a"] != merged["sm_b"]
    merged["session_key"] = merged["session_id"].map(normalize_session_key)
    merged["age"] = merged["session_key"].map(age_map)

    only_a = a2.merge(b2[KEY_COLS], on=KEY_COLS, how="left", indicator=True)
    only_a = only_a[only_a["_merge"] == "left_only"].drop(columns=["_merge"])
    only_b = b2.merge(a2[KEY_COLS], on=KEY_COLS, how="left", indicator=True)
    only_b = only_b[only_b["_merge"] == "left_only"].drop(columns=["_merge"])

    merged.to_csv(out_dir / "overlap_rows.csv", index=False)
    only_a.to_csv(out_dir / "only_in_a_rows.csv", index=False)
    only_b.to_csv(out_dir / "only_in_b_rows.csv", index=False)
    merged[merged["discordant_sm"]].to_csv(out_dir / "disagreement_rows.csv", index=False)

    # Session-level summary (overlap rows only).
    sess = (
        merged.groupby("session_id", dropna=False)
        .agg(
            n_rows=("sm_a", "size"),
            n_discordant=("discordant_sm", "sum"),
            discordant_rate=("discordant_sm", "mean"),
            n_sm_a=("sm_a", "sum"),
            n_sm_b=("sm_b", "sum"),
        )
        .reset_index()
    )
    sess["session_key"] = sess["session_id"].map(normalize_session_key)
    sess["age"] = sess["session_key"].map(age_map)
    sess.to_csv(out_dir / "session_disagreement_summary.csv", index=False)

    # Unit-level (session/cell) summary from full A and B.
    ua = (
        a2.groupby(["session_id", "cell_id"], dropna=False)
        .agg(any_sm_a=("sm_a", "any"), n_cond_a=("condition_1b", "nunique"))
        .reset_index()
    )
    ub = (
        b2.groupby(["session_id", "cell_id"], dropna=False)
        .agg(any_sm_b=("sm_b", "any"), n_cond_b=("condition_1b", "nunique"))
        .reset_index()
    )
    unit = ua.merge(ub, on=["session_id", "cell_id"], how="outer")
    unit["present_a"] = unit["any_sm_a"].notna()
    unit["present_b"] = unit["any_sm_b"].notna()
    unit["any_sm_a"] = unit["any_sm_a"].astype("boolean").fillna(False).astype(bool)
    unit["any_sm_b"] = unit["any_sm_b"].astype("boolean").fillna(False).astype(bool)
    unit["discordant_sm"] = unit["any_sm_a"] != unit["any_sm_b"]
    unit["session_key"] = unit["session_id"].map(normalize_session_key)
    unit["age"] = unit["session_key"].map(age_map)
    unit.to_csv(out_dir / "unit_summary.csv", index=False)
    unit[(unit["present_a"] & unit["present_b"] & unit["discordant_sm"])].to_csv(
        out_dir / "disagreement_units.csv", index=False
    )

    age_summary = pd.DataFrame()
    if age_map:
        age_summary = (
            merged.dropna(subset=["age"])
            .groupby("age", dropna=False)
            .agg(
                n_rows=("sm_a", "size"),
                n_discordant_rows=("discordant_sm", "sum"),
                discordant_rate_row=("discordant_sm", "mean"),
                n_sm_a=("sm_a", "sum"),
                n_sm_b=("sm_b", "sum"),
            )
            .reset_index()
            .sort_values("age", kind="stable")
        )
        age_summary.to_csv(out_dir / "age_disagreement_summary.csv", index=False)

    conf_row = confusion_counts(merged["sm_a"].to_numpy(bool), merged["sm_b"].to_numpy(bool))
    unit_overlap = unit[unit["present_a"] & unit["present_b"]].copy()
    conf_unit = confusion_counts(
        unit_overlap["any_sm_a"].to_numpy(bool),
        unit_overlap["any_sm_b"].to_numpy(bool),
    )

    row_dis = int(merged["discordant_sm"].sum())
    row_n = int(len(merged))
    unit_dis = int(unit_overlap["discordant_sm"].sum())
    unit_n = int(len(unit_overlap))

    summary: dict[str, Any] = {
        "inputs": {
            "a_csv": str(csv_a),
            "b_csv": str(csv_b),
            "name_a": name_a,
            "name_b": name_b,
            "alpha": float(args.alpha),
        },
        "method_a": {
            "n_rows": int(len(a2)),
            "n_units": int(a2[["session_id", "cell_id"]].drop_duplicates().shape[0]),
            "n_sm_rows": int(a2["sm_a"].sum()),
            "sm_row_rate": float(a2["sm_a"].mean()) if len(a2) else float("nan"),
        },
        "method_b": {
            "n_rows": int(len(b2)),
            "n_units": int(b2[["session_id", "cell_id"]].drop_duplicates().shape[0]),
            "n_sm_rows": int(b2["sm_b"].sum()),
            "sm_row_rate": float(b2["sm_b"].mean()) if len(b2) else float("nan"),
        },
        "comparison_rows": {
            "n_overlap": row_n,
            "n_only_a": int(len(only_a)),
            "n_only_b": int(len(only_b)),
            "n_discordant": row_dis,
            "discordant_rate": (float(row_dis) / float(row_n)) if row_n else float("nan"),
            "confusion": conf_row,
        },
        "comparison_units": {
            "n_overlap_units": unit_n,
            "n_discordant_units": unit_dis,
            "discordant_rate_units": (float(unit_dis) / float(unit_n)) if unit_n else float("nan"),
            "confusion": conf_unit,
            "n_only_a_units": int((unit["present_a"] & ~unit["present_b"]).sum()),
            "n_only_b_units": int((~unit["present_a"] & unit["present_b"]).sum()),
        },
    }
    with open(out_dir / "summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    # Plots
    plot_confusion_matrix(
        conf_row,
        label_a=name_a,
        label_b=name_b,
        title="Row-level SM Confusion",
        out_png=out_dir / "confusion_rows.png",
    )
    plot_confusion_matrix(
        conf_unit,
        label_a=name_a,
        label_b=name_b,
        title="Unit-level SM Confusion (any condition)",
        out_png=out_dir / "confusion_units.png",
    )
    plot_pvalue_scatter(merged, label_a=name_a, label_b=name_b, out_png=out_dir / "pvalue_scatter.png")
    plot_delta_hist(merged, out_png=out_dir / "pvalue_delta_hist.png")
    if not age_summary.empty:
        plot_age_disagreement(age_summary, out_png=out_dir / "age_disagreement_rate.png")

    print(f"Wrote comparison output: {out_dir}")
    print(
        f"Row-level overlap/disagreement: {row_n} / {row_dis} "
        f"({(100.0 * row_dis / row_n) if row_n else float('nan'):.2f}%)"
    )
    print(
        f"Unit-level overlap/disagreement: {unit_n} / {unit_dis} "
        f"({(100.0 * unit_dis / unit_n) if unit_n else float('nan'):.2f}%)"
    )
    print(f"Rows only in A: {len(only_a)} | Rows only in B: {len(only_b)}")


if __name__ == "__main__":
    main()
