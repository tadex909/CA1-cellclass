from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def sanitize_token(value: str) -> str:
    out = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip())
    return out.strip("_") or "x"


def session_from_filename(path: Path, suffix: str) -> str | None:
    name = path.name
    if not name.endswith(suffix):
        return None
    return name[: -len(suffix)]


def list_by_session(interim_root: Path, suffix: str) -> tuple[dict[str, Path], list[dict[str, str]]]:
    mapping: dict[str, Path] = {}
    duplicates: list[dict[str, str]] = []
    patt = f"*{suffix}"
    for p in sorted(interim_root.rglob(patt)):
        sid = session_from_filename(p, suffix)
        if sid is None:
            continue
        if sid in mapping:
            duplicates.append(
                {
                    "session_id": sid,
                    "kept": str(mapping[sid]),
                    "ignored": str(p),
                    "reason": "duplicate_session_file",
                }
            )
            continue
        mapping[sid] = p
    return mapping, duplicates


def parse_meta_json(z: np.lib.npyio.NpzFile) -> dict[str, Any]:
    if "meta_json" not in z.files:
        return {}
    try:
        raw = z["meta_json"]
        if hasattr(raw, "item"):
            raw = raw.item()
        if isinstance(raw, bytes):
            txt = raw.decode("utf-8", errors="replace")
        else:
            txt = str(raw)
        val = json.loads(txt)
        return val if isinstance(val, dict) else {}
    except Exception:
        return {}


def infer_allcel_path(pf_path: Path) -> Path:
    name = pf_path.name
    if name.endswith("_pf_c.npz"):
        return pf_path.with_name(name.replace("_pf_c.npz", "_allcel.npz"))
    if name.endswith("_pf.npz"):
        return pf_path.with_name(name.replace("_pf.npz", "_allcel.npz"))
    return pf_path.with_name(f"{pf_path.stem}_allcel.npz")


def try_load_cell_ids(pf_path: Path, expected_n_cells: int) -> np.ndarray:
    allcel_path = infer_allcel_path(pf_path)
    if allcel_path.exists():
        try:
            z = np.load(allcel_path, allow_pickle=True)
            if "allcel__id_cel" in z.files:
                ids = np.asarray(z["allcel__id_cel"]).reshape(-1)
                if ids.size == expected_n_cells:
                    return ids.astype(np.int64, copy=False)
        except Exception:
            pass
    return np.arange(1, int(expected_n_cells) + 1, dtype=np.int64)


def orient_pf_ispf(arr_raw: np.ndarray, n_records_hint: int | None = None) -> np.ndarray:
    arr = np.asarray(arr_raw)
    if arr.ndim != 3:
        raise ValueError(f"Expected 3D pf__ispf_cx, got shape {arr.shape}.")

    # Expected final orientation: [cell, condition, bin]
    if n_records_hint is not None:
        matches = [ax for ax, dim in enumerate(arr.shape) if dim == int(n_records_hint)]
        if len(matches) == 1:
            arr = np.moveaxis(arr, matches[0], 0)

    # If still looks like [condition, bin, cell], rotate to cell-first.
    # Typical data has bins as the largest axis, and cells are not always largest.
    # If axis 0 is large and axis 2 is small, this check will not trigger.
    # We only swap when axis 0 is clearly not the cell axis.
    if arr.shape[0] > max(arr.shape[1], arr.shape[2]) and arr.shape[2] < arr.shape[0]:
        arr = np.moveaxis(arr, 2, 0)

    # For remaining axes, bins usually >= conditions.
    if arr.shape[1] > arr.shape[2]:
        arr = np.swapaxes(arr, 1, 2)

    return arr


def load_pf_has_matrix(npz_path: Path) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    z = np.load(npz_path, allow_pickle=True)
    if "pf__ispf_cx" not in z.files:
        raise KeyError(f"{npz_path} is missing key 'pf__ispf_cx'.")

    meta = parse_meta_json(z)
    n_records_hint: int | None = None
    if "pf_n_records" in meta:
        try:
            n_records_hint = int(meta["pf_n_records"])
        except Exception:
            n_records_hint = None

    arr_raw = np.asarray(z["pf__ispf_cx"])
    arr = orient_pf_ispf(arr_raw, n_records_hint=n_records_hint)
    has_pf = np.any(arr > 0, axis=2)
    cell_ids = try_load_cell_ids(npz_path, expected_n_cells=has_pf.shape[0])

    info = {
        "path": str(npz_path),
        "shape_raw": tuple(int(x) for x in arr_raw.shape),
        "shape_oriented": tuple(int(x) for x in arr.shape),
        "n_cells": int(has_pf.shape[0]),
        "n_conditions": int(has_pf.shape[1]),
    }
    return has_pf.astype(bool), cell_ids, info


def to_long_df(session_id: str, has_pf: np.ndarray, cell_ids: np.ndarray, col_name: str) -> pd.DataFrame:
    n_cells, n_cond = has_pf.shape
    cell_idx = np.repeat(np.arange(1, n_cells + 1, dtype=np.int64), n_cond)
    cond_idx = np.tile(np.arange(1, n_cond + 1, dtype=np.int64), n_cells)
    cell_id_long = np.repeat(np.asarray(cell_ids).reshape(-1), n_cond)
    out = pd.DataFrame(
        {
            "session_id": str(session_id),
            "cell_index_1b": cell_idx,
            "cell_id": cell_id_long,
            "condition_1b": cond_idx,
            col_name: has_pf.reshape(-1).astype(bool),
        }
    )
    out["cell_id"] = pd.to_numeric(out["cell_id"], errors="coerce").astype("Int64")
    return out


def confusion_counts(a: np.ndarray, b: np.ndarray) -> dict[str, int]:
    aa = np.asarray(a, dtype=bool)
    bb = np.asarray(b, dtype=bool)
    return {
        "both_no_pf": int((~aa & ~bb).sum()),
        "a_only_pf": int((aa & ~bb).sum()),
        "b_only_pf": int((~aa & bb).sum()),
        "both_pf": int((aa & bb).sum()),
    }


def plot_confusion(conf: dict[str, int], label_a: str, label_b: str, title: str, out_png: Path) -> None:
    mat = np.array(
        [
            [conf["both_no_pf"], conf["b_only_pf"]],
            [conf["a_only_pf"], conf["both_pf"]],
        ],
        dtype=np.int64,
    )
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    im = ax.imshow(mat, cmap="Blues")
    ax.set_xticks([0, 1], labels=[f"{label_b}: no PF", f"{label_b}: PF"])
    ax.set_yticks([0, 1], labels=[f"{label_a}: no PF", f"{label_a}: PF"])
    ax.set_title(title)
    for i in range(2):
        for j in range(2):
            ax.text(j, i, str(int(mat[i, j])), ha="center", va="center", color="black")
    fig.colorbar(im, ax=ax, fraction=0.05, pad=0.04)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_disagreement_by_condition(cond_summary: pd.DataFrame, out_png: Path) -> None:
    if cond_summary.empty:
        return
    d = cond_summary.sort_values("condition_1b", kind="stable")
    fig, ax = plt.subplots(figsize=(8.6, 4.2))
    ax.bar(d["condition_1b"].astype(str), d["disagreement_rate"], color="#9ecae1", edgecolor="#4a6fa5")
    ax.set_xlabel("Condition (1-based)")
    ax.set_ylabel("Disagreement rate")
    ax.set_title("Row-level PF Disagreement by Condition")
    ax.set_ylim(0.0, max(0.05, float(d["disagreement_rate"].max()) * 1.2))
    ax.grid(alpha=0.2, axis="y")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def write_missing_txt(missing: pd.DataFrame, out_path: Path, suffix_a: str, suffix_b: str) -> None:
    lines: list[str] = []
    n_total = int(len(missing))
    n_missing_a = int((missing["missing"] == "missing_a").sum()) if n_total else 0
    n_missing_b = int((missing["missing"] == "missing_b").sum()) if n_total else 0
    n_missing_both = int((missing["missing"] == "missing_both").sum()) if n_total else 0

    lines.append(f"Missing pair report ({n_total} sessions with missing files)")
    lines.append(f"- missing_a ({suffix_a}): {n_missing_a}")
    lines.append(f"- missing_b ({suffix_b}): {n_missing_b}")
    lines.append(f"- missing_both: {n_missing_both}")
    lines.append("")

    if n_total:
        for row in missing.itertuples(index=False):
            lines.append(f"{row.session_id}: {row.missing}")
            lines.append(f"  path_a: {row.path_a}")
            lines.append(f"  path_b: {row.path_b}")
            lines.append("")

    out_path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Compare session-wise placefield calls from interim PF files "
            "(e.g., *_pf.npz vs *_pf_c.npz)."
        )
    )
    ap.add_argument("--interim_root", type=str, default="data/interim")
    ap.add_argument("--suffix_a", type=str, default="_pf.npz", help="Filename suffix for method A.")
    ap.add_argument("--suffix_b", type=str, default="_pf_c.npz", help="Filename suffix for method B.")
    ap.add_argument("--name_a", type=str, default="random_poisson")
    ap.add_argument("--name_b", type=str, default="circular_shift")
    ap.add_argument("--out_root", type=str, default="results/tables/pf_npz_compare")
    ap.add_argument("--run_tag", type=str, default="")
    return ap


def main() -> None:
    args = build_parser().parse_args()

    interim_root = Path(args.interim_root)
    if not interim_root.exists():
        raise FileNotFoundError(f"interim_root does not exist: {interim_root}")

    suffix_a = str(args.suffix_a)
    suffix_b = str(args.suffix_b)
    name_a = str(args.name_a).strip() or sanitize_token(suffix_a)
    name_b = str(args.name_b).strip() or sanitize_token(suffix_b)

    auto_tag = f"{sanitize_token(name_a)}_vs_{sanitize_token(name_b)}"
    run_tag = sanitize_token(args.run_tag) if str(args.run_tag).strip() else auto_tag

    out_dir = Path(args.out_root) / run_tag
    out_dir.mkdir(parents=True, exist_ok=True)

    by_a, dup_a = list_by_session(interim_root, suffix_a)
    by_b, dup_b = list_by_session(interim_root, suffix_b)

    if dup_a:
        pd.DataFrame(dup_a).to_csv(out_dir / "duplicates_a.csv", index=False)
    if dup_b:
        pd.DataFrame(dup_b).to_csv(out_dir / "duplicates_b.csv", index=False)

    all_sessions = sorted(set(by_a) | set(by_b))
    pair_rows: list[dict[str, Any]] = []
    missing_rows: list[dict[str, Any]] = []

    for sid in all_sessions:
        pa = by_a.get(sid)
        pb = by_b.get(sid)
        has_a = pa is not None
        has_b = pb is not None

        if has_a and has_b:
            pair_rows.append(
                {
                    "session_id": sid,
                    "path_a": str(pa),
                    "path_b": str(pb),
                }
            )
        else:
            if (not has_a) and (not has_b):
                miss = "missing_both"
            elif not has_a:
                miss = "missing_a"
            else:
                miss = "missing_b"
            missing_rows.append(
                {
                    "session_id": sid,
                    "missing": miss,
                    "path_a": str(pa) if pa is not None else "",
                    "path_b": str(pb) if pb is not None else "",
                }
            )

    paired = pd.DataFrame(pair_rows)
    missing = pd.DataFrame(missing_rows)

    paired.to_csv(out_dir / "paired_sessions.csv", index=False)
    missing.to_csv(out_dir / "missing_pairs.csv", index=False)
    write_missing_txt(missing, out_dir / "missing_pairs.txt", suffix_a=suffix_a, suffix_b=suffix_b)

    long_a_list: list[pd.DataFrame] = []
    long_b_list: list[pd.DataFrame] = []
    shape_rows: list[dict[str, Any]] = []
    load_err_rows: list[dict[str, Any]] = []

    for row in paired.itertuples(index=False):
        sid = str(row.session_id)
        pa = Path(row.path_a)
        pb = Path(row.path_b)
        try:
            has_a, ids_a, info_a = load_pf_has_matrix(pa)
            has_b, ids_b, info_b = load_pf_has_matrix(pb)
        except Exception as exc:
            load_err_rows.append(
                {
                    "session_id": sid,
                    "path_a": str(pa),
                    "path_b": str(pb),
                    "error": str(exc),
                }
            )
            continue

        shape_rows.append(
            {
                "session_id": sid,
                "path_a": str(pa),
                "path_b": str(pb),
                "shape_a_raw": str(info_a["shape_raw"]),
                "shape_b_raw": str(info_b["shape_raw"]),
                "shape_a_oriented": str(info_a["shape_oriented"]),
                "shape_b_oriented": str(info_b["shape_oriented"]),
            }
        )

        dfa = to_long_df(sid, has_a, ids_a, col_name="has_pf_a")
        dfb = to_long_df(sid, has_b, ids_b, col_name="has_pf_b")
        long_a_list.append(dfa)
        long_b_list.append(dfb)

    shape_df = pd.DataFrame(shape_rows)
    shape_df.to_csv(out_dir / "shape_alignment.csv", index=False)

    if load_err_rows:
        pd.DataFrame(load_err_rows).to_csv(out_dir / "load_errors.csv", index=False)

    if not long_a_list or not long_b_list:
        summary = {
            "inputs": {
                "interim_root": str(interim_root),
                "suffix_a": suffix_a,
                "suffix_b": suffix_b,
                "name_a": name_a,
                "name_b": name_b,
            },
            "pairs": {
                "n_sessions_total": int(len(all_sessions)),
                "n_sessions_paired": int(len(paired)),
                "n_sessions_missing": int(len(missing)),
            },
            "comparison": "No paired sessions could be loaded.",
        }
        with open(out_dir / "summary.json", "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2)
        print(f"Wrote output: {out_dir}")
        print("No valid paired files to compare.")
        return

    a_long = pd.concat(long_a_list, ignore_index=True)
    b_long = pd.concat(long_b_list, ignore_index=True)

    key = ["session_id", "cell_index_1b", "condition_1b"]
    merged = a_long.merge(b_long, on=key, how="inner", suffixes=("_a", "_b"))
    merged["discordant"] = merged["has_pf_a"] != merged["has_pf_b"]

    only_a = a_long.merge(b_long[key], on=key, how="left", indicator=True)
    only_a = only_a[only_a["_merge"] == "left_only"].drop(columns=["_merge"])
    only_b = b_long.merge(a_long[key], on=key, how="left", indicator=True)
    only_b = only_b[only_b["_merge"] == "left_only"].drop(columns=["_merge"])

    # Keep one cell_id column for convenience.
    merged["cell_id"] = merged["cell_id_a"]
    mismatch_cell_id = (
        merged["cell_id_a"].notna()
        & merged["cell_id_b"].notna()
        & (merged["cell_id_a"].astype("Int64") != merged["cell_id_b"].astype("Int64"))
    )
    merged.loc[mismatch_cell_id, "cell_id"] = pd.NA

    overlap_cols = [
        "session_id",
        "cell_index_1b",
        "cell_id",
        "condition_1b",
        "has_pf_a",
        "has_pf_b",
        "discordant",
        "cell_id_a",
        "cell_id_b",
    ]
    merged[overlap_cols].to_csv(out_dir / "overlap_rows.csv", index=False)
    merged[merged["discordant"]][overlap_cols].to_csv(out_dir / "disagreement_rows.csv", index=False)
    only_a.to_csv(out_dir / "only_in_a_rows.csv", index=False)
    only_b.to_csv(out_dir / "only_in_b_rows.csv", index=False)

    sess = (
        merged.groupby("session_id", dropna=False)
        .agg(
            n_rows=("discordant", "size"),
            n_discordant=("discordant", "sum"),
            disagreement_rate=("discordant", "mean"),
            n_pf_a=("has_pf_a", "sum"),
            n_pf_b=("has_pf_b", "sum"),
        )
        .reset_index()
        .sort_values("disagreement_rate", ascending=False, kind="stable")
    )
    sess.to_csv(out_dir / "session_disagreement_summary.csv", index=False)

    cond = (
        merged.groupby("condition_1b", dropna=False)
        .agg(
            n_rows=("discordant", "size"),
            n_discordant=("discordant", "sum"),
            disagreement_rate=("discordant", "mean"),
            n_pf_a=("has_pf_a", "sum"),
            n_pf_b=("has_pf_b", "sum"),
        )
        .reset_index()
        .sort_values("condition_1b", kind="stable")
    )
    cond.to_csv(out_dir / "condition_disagreement_summary.csv", index=False)

    ua = (
        a_long.groupby(["session_id", "cell_index_1b"], dropna=False)
        .agg(any_pf_a=("has_pf_a", "any"), cell_id_a=("cell_id", "first"))
        .reset_index()
    )
    ub = (
        b_long.groupby(["session_id", "cell_index_1b"], dropna=False)
        .agg(any_pf_b=("has_pf_b", "any"), cell_id_b=("cell_id", "first"))
        .reset_index()
    )
    unit = ua.merge(ub, on=["session_id", "cell_index_1b"], how="outer")
    unit["present_a"] = unit["any_pf_a"].notna()
    unit["present_b"] = unit["any_pf_b"].notna()
    unit["any_pf_a"] = unit["any_pf_a"].astype("boolean").fillna(False).astype(bool)
    unit["any_pf_b"] = unit["any_pf_b"].astype("boolean").fillna(False).astype(bool)
    unit["discordant"] = unit["any_pf_a"] != unit["any_pf_b"]
    unit["cell_id"] = unit["cell_id_a"]
    mismatch_u_cell = (
        unit["cell_id_a"].notna()
        & unit["cell_id_b"].notna()
        & (unit["cell_id_a"].astype("Int64") != unit["cell_id_b"].astype("Int64"))
    )
    unit.loc[mismatch_u_cell, "cell_id"] = pd.NA
    unit.to_csv(out_dir / "unit_summary.csv", index=False)
    unit[(unit["present_a"] & unit["present_b"] & unit["discordant"])].to_csv(
        out_dir / "disagreement_units.csv", index=False
    )

    conf_row = confusion_counts(merged["has_pf_a"].to_numpy(bool), merged["has_pf_b"].to_numpy(bool))
    unit_overlap = unit[unit["present_a"] & unit["present_b"]].copy()
    conf_unit = confusion_counts(unit_overlap["any_pf_a"].to_numpy(bool), unit_overlap["any_pf_b"].to_numpy(bool))

    plot_confusion(
        conf_row,
        label_a=name_a,
        label_b=name_b,
        title="Row-level PF Confusion",
        out_png=out_dir / "confusion_rows.png",
    )
    plot_confusion(
        conf_unit,
        label_a=name_a,
        label_b=name_b,
        title="Unit-level PF Confusion (any condition)",
        out_png=out_dir / "confusion_units.png",
    )
    plot_disagreement_by_condition(cond, out_png=out_dir / "disagreement_by_condition.png")

    row_dis = int(merged["discordant"].sum())
    row_n = int(len(merged))
    unit_dis = int(unit_overlap["discordant"].sum())
    unit_n = int(len(unit_overlap))

    summary = {
        "inputs": {
            "interim_root": str(interim_root),
            "suffix_a": suffix_a,
            "suffix_b": suffix_b,
            "name_a": name_a,
            "name_b": name_b,
        },
        "pairs": {
            "n_sessions_total": int(len(all_sessions)),
            "n_sessions_paired": int(len(paired)),
            "n_sessions_missing": int(len(missing)),
            "n_duplicates_a": int(len(dup_a)),
            "n_duplicates_b": int(len(dup_b)),
            "n_load_errors": int(len(load_err_rows)),
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

    print(f"Wrote comparison output: {out_dir}")
    print(
        f"Paired sessions: {len(paired)} | Missing sessions: {len(missing)} | "
        f"Load errors: {len(load_err_rows)}"
    )
    print(
        f"Row-level overlap/disagreement: {row_n} / {row_dis} "
        f"({(100.0 * row_dis / row_n) if row_n else float('nan'):.2f}%)"
    )
    print(
        f"Unit-level overlap/disagreement: {unit_n} / {unit_dis} "
        f"({(100.0 * unit_dis / unit_n) if unit_n else float('nan'):.2f}%)"
    )


if __name__ == "__main__":
    main()
