from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_meta_json(raw: np.ndarray) -> dict[str, Any]:
    try:
        txt = raw.tobytes().decode("utf-8", errors="ignore")
        obj = json.loads(txt)
        if isinstance(obj, dict):
            return obj
    except Exception:
        pass
    return {}


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Plot ratemaps for disagreement rows between two PF methods "
            "(one figure per session/cell/condition row)."
        )
    )
    ap.add_argument(
        "--compare_dir",
        type=str,
        default="results/tables/pf_npz_compare/random_poisson_vs_circular_shift",
        help="Folder produced by compare_pf_npz_methods.py.",
    )
    ap.add_argument(
        "--disagreement_csv",
        type=str,
        default="",
        help="Optional explicit disagreement_rows.csv path. Overrides --compare_dir when provided.",
    )
    ap.add_argument(
        "--ratemap_root",
        type=str,
        default="results/ratemap",
        help="Root containing <SESSION>_rmap.npz files.",
    )
    ap.add_argument(
        "--out_root",
        type=str,
        default="results/figures/pf_method_disagreement_maps",
    )
    ap.add_argument(
        "--max_rows",
        type=int,
        default=0,
        help="Optional cap on number of disagreement rows to plot (0 = all).",
    )
    ap.add_argument(
        "--max_units",
        type=int,
        default=0,
        help="Deprecated alias for --max_rows.",
    )
    ap.add_argument(
        "--heatmap_percentile",
        type=float,
        default=99.0,
        help="Upper percentile for heatmap color scaling.",
    )
    ap.add_argument("--dpi", type=int, default=140)
    ap.add_argument("--overwrite", action="store_true")
    return ap


def resolve_disagreement_csv(compare_dir: Path, disagreement_csv: str) -> Path:
    if str(disagreement_csv).strip():
        p = Path(disagreement_csv)
        if not p.exists():
            raise FileNotFoundError(f"Missing --disagreement_csv: {p}")
        return p
    p = compare_dir / "disagreement_rows.csv"
    if not p.exists():
        raise FileNotFoundError(
            f"Could not find disagreement_rows.csv in {compare_dir}. "
            "Run compare_pf_npz_methods.py first or pass --disagreement_csv."
        )
    return p


def find_rmap_file(ratemap_root: Path, session_id: str) -> Path | None:
    cands = sorted(ratemap_root.rglob(f"{session_id}_rmap.npz"))
    if not cands:
        return None
    if len(cands) > 1:
        print(f"WARN: multiple rmap files for {session_id}; using {cands[-1]}")
    return cands[-1]


def prepare_disagreement_rows(d_rows: pd.DataFrame) -> pd.DataFrame:
    d = d_rows.copy()
    required = ["session_id", "cell_index_1b", "condition_1b", "has_pf_a", "has_pf_b"]
    missing = [c for c in required if c not in d.columns]
    if missing:
        raise KeyError(f"Missing required columns in disagreement CSV: {missing}")

    d["session_id"] = d["session_id"].astype(str)
    d["cell_index_1b"] = pd.to_numeric(d["cell_index_1b"], errors="coerce").astype("Int64")
    d["cell_id"] = pd.to_numeric(d.get("cell_id", pd.Series(dtype=float)), errors="coerce")
    d["condition_1b"] = pd.to_numeric(d["condition_1b"], errors="coerce").astype("Int64")
    d["has_pf_a"] = d["has_pf_a"].astype(bool)
    d["has_pf_b"] = d["has_pf_b"].astype(bool)

    d = d.dropna(subset=["session_id", "cell_index_1b", "condition_1b"]).copy()
    d["cell_index_1b"] = d["cell_index_1b"].astype(np.int64)
    d["condition_1b"] = d["condition_1b"].astype(np.int64)

    def row_class(row: pd.Series) -> str:
        a = bool(row["has_pf_a"])
        b = bool(row["has_pf_b"])
        if a and (not b):
            return "pf_a_no_pf_b"
        if (not a) and b:
            return "no_pf_a_pf_b"
        return "other"

    d["row_class"] = d.apply(row_class, axis=1)

    key = ["session_id", "cell_index_1b", "condition_1b", "has_pf_a", "has_pf_b"]
    d = d.sort_values(key, kind="stable").drop_duplicates(subset=key, keep="first").reset_index(drop=True)
    return d


def map_pf_cond_to_rmap_cond(cond_pf_1b: int, n_cond_rmap: int) -> tuple[int, str]:
    if n_cond_rmap <= 0:
        raise ValueError("n_cond_rmap must be >= 1")

    c = int(cond_pf_1b)
    if 1 <= c <= int(n_cond_rmap):
        return c, "direct"

    if int(n_cond_rmap) == 2:
        # Common case: PF table has 10 condition slots (5 pairs x direction),
        # while ratemap condition axis is directional (1/2).
        mapped = 1 if (c % 2 == 1) else 2
        return mapped, "odd_even_to_1_2"

    mapped = ((c - 1) % int(n_cond_rmap)) + 1
    return mapped, "modulo_wrap"


def cell_axis_index(cell_ids: np.ndarray, cell_id: float, cell_index_1b: int) -> tuple[int | None, str]:
    if np.isfinite(cell_id):
        cid = int(cell_id)
        hits = np.where(cell_ids == cid)[0]
        if hits.size > 0:
            return int(hits[0]), "cell_id"

    idx0 = int(cell_index_1b) - 1
    if 0 <= idx0 < int(cell_ids.size):
        return idx0, "cell_index_1b"
    return None, "not_found"


def load_rmap_payload(rmap_path: Path) -> dict[str, Any]:
    with np.load(rmap_path, allow_pickle=False) as z:
        req = [
            "cell_ids",
            "idcond_t",
            "xbin_centers",
            "rmap__fr_tx_ux",
            "rmap__fr_s_tx_ux",
            "rmap__fr_cx_ux",
            "rmap__fr_s_cx_ux",
        ]
        miss = [k for k in req if k not in z.files]
        if miss:
            raise KeyError(f"missing keys in {rmap_path}: {miss}")

        return {
            "cell_ids": z["cell_ids"].astype(np.int64, copy=False),
            "idcond_t": z["idcond_t"].astype(np.int64, copy=False),
            "x": z["xbin_centers"].astype(np.float64, copy=False),
            "fr_tx": z["rmap__fr_tx_ux"].astype(np.float64, copy=False),
            "fr_s_tx": z["rmap__fr_s_tx_ux"].astype(np.float64, copy=False),
            "fr_cx": z["rmap__fr_cx_ux"].astype(np.float64, copy=False),
            "fr_s_cx": z["rmap__fr_s_cx_ux"].astype(np.float64, copy=False),
            "meta": parse_meta_json(z["meta_json"]) if "meta_json" in z.files else {},
        }


def save_row_plot(
    *,
    out_png: Path,
    session_id: str,
    row_class: str,
    condition_pf_1b: int,
    condition_rmap_1b: int,
    cond_mapping: str,
    cell_index_1b: int,
    cell_id: float,
    matched_by: str,
    x: np.ndarray,
    idcond_t: np.ndarray,
    fr_tx_u: np.ndarray,
    fr_s_tx_u: np.ndarray,
    fr_cx_u: np.ndarray,
    fr_s_cx_u: np.ndarray,
    rmap_path: Path,
    meta: dict[str, Any],
    heatmap_percentile: float,
    dpi: int,
) -> tuple[int, int]:
    c0 = int(condition_rmap_1b) - 1
    n_cond = int(fr_cx_u.shape[0])
    if c0 < 0 or c0 >= n_cond:
        raise IndexError(
            f"Mapped condition {condition_rmap_1b} is out of range for n_cond_rmap={n_cond}."
        )

    # Align direction: even conditions are mirrored so all maps are shown in the same orientation.
    reflect_even = (int(condition_pf_1b) % 2) == 0

    y_raw = np.asarray(fr_cx_u[c0], dtype=np.float64)
    y_sm = np.asarray(fr_s_cx_u[c0], dtype=np.float64)
    if reflect_even:
        y_raw = y_raw[::-1]
        y_sm = y_sm[::-1]
    y_stack = np.concatenate([y_raw.ravel(), y_sm.ravel()])
    y_stack = y_stack[np.isfinite(y_stack)]
    y_upper = float(np.max(y_stack)) if y_stack.size else 1.0
    if (not np.isfinite(y_upper)) or (y_upper <= 0):
        y_upper = 1.0

    idx_cond = np.where(np.asarray(idcond_t, dtype=np.int64) == int(condition_rmap_1b))[0]
    heat = np.asarray(fr_s_tx_u[idx_cond, :], dtype=np.float64) if idx_cond.size else np.empty((0, x.size), float)
    if reflect_even and idx_cond.size:
        heat = heat[:, ::-1]

    fig, axes = plt.subplots(2, 1, figsize=(9.6, 7.2), squeeze=False)
    ax_heat = axes[0, 0]
    ax_mean = axes[1, 0]

    if idx_cond.size == 0:
        ax_heat.text(
            0.5,
            0.5,
            f"No trials found for rmap condition {condition_rmap_1b}",
            ha="center",
            va="center",
            transform=ax_heat.transAxes,
        )
        ax_heat.set_axis_off()
    else:
        heat_fin = heat[np.isfinite(heat)]
        vmax = np.nan
        if heat_fin.size:
            vmax = float(np.nanpercentile(heat_fin, float(heatmap_percentile)))
        if (not np.isfinite(vmax)) or (vmax <= 0):
            vmax = float(np.nanmax(heat_fin)) if heat_fin.size else 1.0
        if (not np.isfinite(vmax)) or (vmax <= 0):
            vmax = 1.0

        im = ax_heat.imshow(
            heat,
            aspect="auto",
            origin="lower",
            interpolation="nearest",
            extent=[float(x[0]), float(x[-1]), 1, int(idx_cond.size)],
            vmin=0.0,
            vmax=float(vmax),
            cmap="viridis",
        )
        ax_heat.set_title(
            f"Smoothed trial ratemaps | rmap cond {condition_rmap_1b} | n_trials={int(idx_cond.size)}",
            fontsize=10,
        )
        ax_heat.set_xlabel("Position")
        ax_heat.set_ylabel("Trial")
        cbar = fig.colorbar(im, ax=ax_heat, fraction=0.046, pad=0.04)
        cbar.set_label("FR (smoothed)")

    ax_mean.plot(x, y_raw, color="0.65", linewidth=1.6, label="mean raw")
    ax_mean.plot(x, y_sm, color="#1f77b4", linewidth=2.0, label="mean smoothed")
    ax_mean.set_title(f"Condition mean ratemap | rmap cond {condition_rmap_1b}", fontsize=10)
    ax_mean.set_xlabel("Position")
    ax_mean.set_ylabel("FR")
    ax_mean.set_ylim(bottom=0.0, top=y_upper * 1.05)
    ax_mean.grid(alpha=0.25)
    ax_mean.legend(loc="best", fontsize=8)

    cell_txt = f"{int(cell_id)}" if np.isfinite(cell_id) else "NA"
    suptitle = (
        f"{session_id} | class={row_class} | cell_id={cell_txt} | cell_idx={cell_index_1b} "
        f"(matched by {matched_by})\n"
        f"PF cond={condition_pf_1b} -> rmap cond={condition_rmap_1b} [{cond_mapping}] | reflected_even={reflect_even}"
    )
    fig.suptitle(suptitle, fontsize=10, y=0.99)

    smooth = meta.get("smooth_sigma_bins")
    xrem = meta.get("xbin_rem")
    info = f"rmap={rmap_path.name}"
    if smooth is not None:
        info += f" | smooth_sigma_bins={smooth}"
    if xrem is not None:
        info += f" | xbin_rem={xrem}"
    fig.text(0.01, 0.01, info, fontsize=8, ha="left", va="bottom")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(out_png, dpi=int(dpi), bbox_inches="tight")
    plt.close(fig)

    return int(idx_cond.size), int(n_cond)


def main() -> None:
    args = build_parser().parse_args()

    compare_dir = Path(args.compare_dir)
    ratemap_root = Path(args.ratemap_root)
    out_root = Path(args.out_root)

    if not ratemap_root.exists():
        raise FileNotFoundError(f"Missing ratemap root: {ratemap_root}")

    dis_csv = resolve_disagreement_csv(compare_dir, args.disagreement_csv)
    d_rows = pd.read_csv(dis_csv)
    if d_rows.empty:
        print(f"No disagreement rows in: {dis_csv}")
        return

    d = prepare_disagreement_rows(d_rows)
    if d.empty:
        print(f"No valid disagreement rows derived from: {dis_csv}")
        return

    max_rows = int(args.max_rows) if int(args.max_rows) > 0 else int(args.max_units)
    if max_rows > 0:
        d = d.head(max_rows).copy()

    run_tag = compare_dir.name if compare_dir.exists() else "manual_compare"
    out_dir = out_root / run_tag
    out_dir.mkdir(parents=True, exist_ok=True)

    # Cache ratemap payload per session.
    session_rmap_path: dict[str, Path | None] = {}
    session_payload: dict[str, dict[str, Any]] = {}

    for sid in d["session_id"].astype(str).drop_duplicates():
        p = find_rmap_file(ratemap_root, sid)
        session_rmap_path[sid] = p
        if p is not None:
            try:
                session_payload[sid] = load_rmap_payload(p)
            except Exception as exc:
                print(f"WARN: failed to load rmap for {sid}: {exc}")

    plot_rows: list[dict[str, Any]] = []
    missing_rows: list[dict[str, Any]] = []

    n_total = int(len(d))
    for i, row in enumerate(d.itertuples(index=False), start=1):
        sid = str(row.session_id)
        cidx = int(row.cell_index_1b)
        cid = float(row.cell_id) if pd.notna(row.cell_id) else float("nan")
        cond_pf_1b = int(row.condition_1b)
        row_class = str(row.row_class)

        rmap_path = session_rmap_path.get(sid)
        payload = session_payload.get(sid)

        if (rmap_path is None) or (payload is None):
            missing_rows.append(
                {
                    "session_id": sid,
                    "cell_index_1b": cidx,
                    "cell_id": cid,
                    "condition_pf_1b": cond_pf_1b,
                    "row_class": row_class,
                    "reason": "missing_or_unreadable_rmap",
                    "rmap_path": "" if rmap_path is None else str(rmap_path),
                }
            )
            continue

        try:
            cell_ids = payload["cell_ids"]
            idx_u, matched_by = cell_axis_index(cell_ids, cid, cidx)
            if idx_u is None:
                missing_rows.append(
                    {
                        "session_id": sid,
                        "cell_index_1b": cidx,
                        "cell_id": cid,
                        "condition_pf_1b": cond_pf_1b,
                        "row_class": row_class,
                        "reason": "cell_not_found_in_rmap",
                        "rmap_path": str(rmap_path),
                    }
                )
                continue

            fr_cx_u = payload["fr_cx"][idx_u]
            fr_s_cx_u = payload["fr_s_cx"][idx_u]
            fr_tx_u = payload["fr_tx"][idx_u]
            fr_s_tx_u = payload["fr_s_tx"][idx_u]
            n_cond_rmap = int(fr_cx_u.shape[0])
            reflect_even = (int(cond_pf_1b) % 2) == 0

            cond_rmap_1b, cond_mapping = map_pf_cond_to_rmap_cond(cond_pf_1b, n_cond_rmap)

            cell_label = f"{int(cid)}" if np.isfinite(cid) else f"idx_{cidx:04d}"
            out_png = (
                out_dir
                / row_class
                / sid
                / f"cell_{cell_label}"
                / f"cond_pf_{cond_pf_1b:02d}__rmap_{cond_rmap_1b:02d}.png"
            )
            if out_png.exists() and not bool(args.overwrite):
                plot_rows.append(
                    {
                        "session_id": sid,
                        "cell_index_1b": cidx,
                        "cell_id": cid,
                        "condition_pf_1b": cond_pf_1b,
                        "condition_rmap_1b": cond_rmap_1b,
                        "row_class": row_class,
                        "status": "exists_skip",
                        "out_png": str(out_png),
                        "rmap_path": str(rmap_path),
                        "matched_by": matched_by,
                        "cond_mapping": cond_mapping,
                        "reflected_even": bool(reflect_even),
                    }
                )
                continue

            n_trials, n_cond = save_row_plot(
                out_png=out_png,
                session_id=sid,
                row_class=row_class,
                condition_pf_1b=cond_pf_1b,
                condition_rmap_1b=cond_rmap_1b,
                cond_mapping=cond_mapping,
                cell_index_1b=cidx,
                cell_id=cid,
                matched_by=matched_by,
                x=payload["x"],
                idcond_t=payload["idcond_t"],
                fr_tx_u=fr_tx_u,
                fr_s_tx_u=fr_s_tx_u,
                fr_cx_u=fr_cx_u,
                fr_s_cx_u=fr_s_cx_u,
                rmap_path=rmap_path,
                meta=payload["meta"],
                heatmap_percentile=float(args.heatmap_percentile),
                dpi=int(args.dpi),
            )

            plot_rows.append(
                {
                    "session_id": sid,
                    "cell_index_1b": cidx,
                    "cell_id": cid,
                    "condition_pf_1b": cond_pf_1b,
                    "condition_rmap_1b": cond_rmap_1b,
                    "row_class": row_class,
                    "status": "plotted",
                    "out_png": str(out_png),
                    "rmap_path": str(rmap_path),
                    "matched_by": matched_by,
                    "cond_mapping": cond_mapping,
                    "reflected_even": bool(reflect_even),
                    "n_trials_rmap_condition": int(n_trials),
                    "n_conditions_in_rmap": int(n_cond),
                }
            )

        except Exception as exc:
            missing_rows.append(
                {
                    "session_id": sid,
                    "cell_index_1b": cidx,
                    "cell_id": cid,
                    "condition_pf_1b": cond_pf_1b,
                    "row_class": row_class,
                    "reason": "plot_error",
                    "rmap_path": str(rmap_path),
                    "error": str(exc),
                }
            )

        if (i % 25 == 0) or (i == n_total):
            print(f"progress: {i}/{n_total} rows")

    plotted_cols = [
        "session_id",
        "cell_index_1b",
        "cell_id",
        "condition_pf_1b",
        "condition_rmap_1b",
        "row_class",
        "status",
        "out_png",
        "rmap_path",
        "matched_by",
        "cond_mapping",
        "reflected_even",
        "n_trials_rmap_condition",
        "n_conditions_in_rmap",
    ]
    missing_cols = [
        "session_id",
        "cell_index_1b",
        "cell_id",
        "condition_pf_1b",
        "row_class",
        "reason",
        "rmap_path",
        "error",
    ]

    plotted = pd.DataFrame(plot_rows, columns=plotted_cols)
    missing = pd.DataFrame(missing_rows, columns=missing_cols)
    plotted.to_csv(out_dir / "plotted_rows.csv", index=False)
    missing.to_csv(out_dir / "plot_missing_or_errors.csv", index=False)

    summary = {
        "inputs": {
            "compare_dir": str(compare_dir),
            "disagreement_csv": str(dis_csv),
            "ratemap_root": str(ratemap_root),
            "max_rows": int(max_rows),
            "heatmap_percentile": float(args.heatmap_percentile),
        },
        "counts": {
            "n_rows_requested": int(n_total),
            "n_plotted": int((plotted["status"] == "plotted").sum()) if not plotted.empty else 0,
            "n_exists_skip": int((plotted["status"] == "exists_skip").sum()) if not plotted.empty else 0,
            "n_missing_or_errors": int(len(missing)),
            "n_class_pf_a_no_pf_b": int((d["row_class"] == "pf_a_no_pf_b").sum()),
            "n_class_no_pf_a_pf_b": int((d["row_class"] == "no_pf_a_pf_b").sum()),
            "n_class_other": int((d["row_class"] == "other").sum()),
        },
        "out_dir": str(out_dir),
    }
    with open(out_dir / "summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"Wrote disagreement ratemap plots to: {out_dir}")
    print(
        f"requested={summary['counts']['n_rows_requested']} | "
        f"plotted={summary['counts']['n_plotted']} | "
        f"skipped={summary['counts']['n_exists_skip']} | "
        f"missing_or_errors={summary['counts']['n_missing_or_errors']}"
    )


if __name__ == "__main__":
    main()
