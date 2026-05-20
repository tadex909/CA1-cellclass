from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
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


@dataclass(frozen=True)
class CueLayout:
    name: str
    cue_rich: tuple[float, float]
    cue_poor: tuple[float, float]
    object_zone: tuple[float, float]
    object_centers: tuple[float, ...]
    moved_object_span: tuple[float, float] | None = None
    note: str | None = None


@dataclass(frozen=True)
class DirectionView:
    direction_label: str
    pf_cond_1b: int
    base_condition_1b: int
    condition_name: str
    condition_rmap_1b: int
    cond_mapping: str
    reflect_even: bool
    mean_raw: np.ndarray
    mean_sm: np.ndarray
    heat: np.ndarray
    n_trials: int
    state: dict[str, Any] | None


def canonical_condition_name(name: str | None) -> str:
    if name is None:
        return ""
    return str(name).strip().upper()


def cue_layout_for_condition(name: str | None) -> CueLayout | None:
    cname = canonical_condition_name(name)
    base = {
        "cue_rich": (13.0, 43.0),
        "cue_poor": (43.0, 81.0),
        "object_zone": (81.0, 96.0),
    }
    if cname in {"PO", "PO2", "PO3", "PONM"}:
        return CueLayout(
            name=cname,
            cue_rich=base["cue_rich"],
            cue_poor=base["cue_poor"],
            object_zone=base["object_zone"],
            object_centers=(20.0, 36.0, 88.0),
        )
    if cname in {"POM", "POMB"}:
        return CueLayout(
            name=cname,
            cue_rich=base["cue_rich"],
            cue_poor=base["cue_poor"],
            object_zone=base["object_zone"],
            object_centers=(20.0, 64.0, 88.0),
            moved_object_span=(57.0, 72.0),
            note="Moved second object",
        )
    if cname == "PNO":
        return CueLayout(
            name=cname,
            cue_rich=base["cue_rich"],
            cue_poor=base["cue_poor"],
            object_zone=base["object_zone"],
            object_centers=tuple(),
            note="No objects",
        )
    return None


def load_compare_labels(compare_dir: Path) -> tuple[str, str]:
    summary_path = compare_dir / "summary.json"
    if not summary_path.exists():
        return "A", "B"
    try:
        with open(summary_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        inputs = data.get("inputs", {}) if isinstance(data, dict) else {}
        label_a = str(inputs.get("name_a", "A")).strip() or "A"
        label_b = str(inputs.get("name_b", "B")).strip() or "B"
        return label_a, label_b
    except Exception:
        return "A", "B"


def load_overlap_lookup(compare_dir: Path) -> dict[tuple[str, int, int], dict[str, Any]]:
    overlap_path = compare_dir / "overlap_rows.csv"
    if not overlap_path.exists():
        return {}
    d = pd.read_csv(overlap_path)
    required = ["session_id", "cell_index_1b", "condition_1b", "has_pf_a", "has_pf_b"]
    if any(col not in d.columns for col in required):
        return {}
    out: dict[tuple[str, int, int], dict[str, Any]] = {}
    for row in d.itertuples(index=False):
        try:
            key = (str(row.session_id), int(row.cell_index_1b), int(row.condition_1b))
        except Exception:
            continue
        out[key] = {
            "has_pf_a": bool(row.has_pf_a),
            "has_pf_b": bool(row.has_pf_b),
            "discordant": bool(getattr(row, "discordant", False)),
        }
    return out


def load_condition_name_map(traj_path: Path) -> dict[int, str]:
    if not traj_path.exists():
        return {}
    with np.load(traj_path, allow_pickle=False) as z:
        if "traj__Cond" not in z.files:
            return {}
        cond = np.asarray(z["traj__Cond"]).astype(np.int64, copy=False).ravel()
        if "traj__condition" in z.files:
            names = np.asarray(z["traj__condition"]).astype(str).ravel()
        elif "traj__condition__json" in z.files:
            raw = z["traj__condition__json"].tobytes().decode("utf-8", errors="ignore")
            names = np.asarray(json.loads(raw), dtype=object).astype(str).ravel()
        else:
            return {}
    if names.size != cond.size:
        return {}
    out: dict[int, str] = {}
    for c, name in zip(cond, names):
        c_int = int(c)
        name_str = str(name).strip()
        if c_int not in out and name_str:
            out[c_int] = name_str
    return out


def format_pf_state(state: dict[str, Any] | None, label_a: str, label_b: str) -> str:
    if state is None:
        return "status unavailable"
    txt_a = "PF" if bool(state.get("has_pf_a", False)) else "no PF"
    txt_b = "PF" if bool(state.get("has_pf_b", False)) else "no PF"
    return f"{label_a}={txt_a} | {label_b}={txt_b}"


def add_cue_overlays(
    ax: plt.Axes,
    layout: CueLayout | None,
    *,
    add_labels: bool,
) -> None:
    if layout is None:
        return
    ax.axvspan(layout.cue_rich[0], layout.cue_rich[1], color="#f4c095", alpha=0.12, zorder=0)
    ax.axvspan(layout.cue_poor[0], layout.cue_poor[1], color="#d7ebba", alpha=0.10, zorder=0)
    ax.axvspan(layout.object_zone[0], layout.object_zone[1], color="#a8cbe6", alpha=0.12, zorder=0)
    if layout.moved_object_span is not None:
        ax.axvspan(
            layout.moved_object_span[0],
            layout.moved_object_span[1],
            color="#f7a1a1",
            alpha=0.16,
            zorder=0.1,
        )
    for center in layout.object_centers:
        ax.axvline(center, color="#8b1e3f", linestyle="--", linewidth=1.1, alpha=0.85, zorder=3)
    if add_labels:
        ax.text(
            np.mean(layout.cue_rich),
            0.98,
            "cue-rich",
            ha="center",
            va="top",
            fontsize=8,
            transform=ax.get_xaxis_transform(),
            color="#7a4b00",
        )
        ax.text(
            np.mean(layout.cue_poor),
            0.98,
            "cue-poor",
            ha="center",
            va="top",
            fontsize=8,
            transform=ax.get_xaxis_transform(),
            color="#426a1f",
        )
        ax.text(
            np.mean(layout.object_zone),
            0.98,
            "object zone",
            ha="center",
            va="top",
            fontsize=8,
            transform=ax.get_xaxis_transform(),
            color="#2b5b84",
        )
        if layout.note:
            ax.text(
                0.99,
                0.02,
                layout.note,
                ha="right",
                va="bottom",
                fontsize=8,
                transform=ax.transAxes,
                color="#6b1e1e",
            )


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


def resolve_existing_path(raw_path: Any, *, working_root: Path) -> Path | None:
    if raw_path is None:
        return None
    txt = str(raw_path).strip()
    if not txt:
        return None
    p = Path(txt)
    if p.exists():
        return p
    if not p.is_absolute():
        cand = working_root / p
        if cand.exists():
            return cand
    return None


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


def build_direction_view(
    *,
    session_id: str,
    cell_index_1b: int,
    pf_cond_1b: int,
    x: np.ndarray,
    idcond_t: np.ndarray,
    fr_s_tx_u: np.ndarray,
    fr_cx_u: np.ndarray,
    fr_s_cx_u: np.ndarray,
    n_cond_rmap: int,
    condition_name_map: dict[int, str],
    overlap_lookup: dict[tuple[str, int, int], dict[str, Any]],
) -> DirectionView:
    cond_rmap_1b, cond_mapping = map_pf_cond_to_rmap_cond(pf_cond_1b, n_cond_rmap)
    c0 = int(cond_rmap_1b) - 1
    if c0 < 0 or c0 >= int(fr_cx_u.shape[0]):
        raise IndexError(
            f"Mapped condition {cond_rmap_1b} is out of range for n_cond_rmap={int(fr_cx_u.shape[0])}."
        )

    reflect_even = (int(pf_cond_1b) % 2) == 0
    direction_label = "Backward (B)" if reflect_even else "Forward (W)"
    base_condition_1b = (int(pf_cond_1b) + 1) // 2
    condition_name = condition_name_map.get(base_condition_1b, f"cond_{base_condition_1b}")

    mean_raw = np.asarray(fr_cx_u[c0], dtype=np.float64)
    mean_sm = np.asarray(fr_s_cx_u[c0], dtype=np.float64)
    if reflect_even:
        mean_raw = mean_raw[::-1]
        mean_sm = mean_sm[::-1]

    idx_cond = np.where(np.asarray(idcond_t, dtype=np.int64) == int(cond_rmap_1b))[0]
    heat = np.asarray(fr_s_tx_u[idx_cond, :], dtype=np.float64) if idx_cond.size else np.empty((0, x.size), float)
    if reflect_even and idx_cond.size:
        heat = heat[:, ::-1]

    state = overlap_lookup.get((session_id, int(cell_index_1b), int(pf_cond_1b)))
    return DirectionView(
        direction_label=direction_label,
        pf_cond_1b=int(pf_cond_1b),
        base_condition_1b=base_condition_1b,
        condition_name=str(condition_name),
        condition_rmap_1b=int(cond_rmap_1b),
        cond_mapping=str(cond_mapping),
        reflect_even=bool(reflect_even),
        mean_raw=mean_raw,
        mean_sm=mean_sm,
        heat=heat,
        n_trials=int(idx_cond.size),
        state=state,
    )


def save_row_plot(
    *,
    out_png: Path,
    label_a: str,
    label_b: str,
    session_id: str,
    row_class: str,
    condition_pf_1b: int,
    cell_index_1b: int,
    cell_id: float,
    matched_by: str,
    x: np.ndarray,
    idcond_t: np.ndarray,
    fr_s_tx_u: np.ndarray,
    fr_cx_u: np.ndarray,
    fr_s_cx_u: np.ndarray,
    rmap_path: Path,
    meta: dict[str, Any],
    condition_name_map: dict[int, str],
    overlap_lookup: dict[tuple[str, int, int], dict[str, Any]],
    heatmap_percentile: float,
    dpi: int,
) -> tuple[int, int]:
    n_cond = int(fr_cx_u.shape[0])
    base_condition_1b = (int(condition_pf_1b) + 1) // 2
    forward_pf_cond_1b = (2 * base_condition_1b) - 1
    backward_pf_cond_1b = 2 * base_condition_1b
    condition_name = condition_name_map.get(base_condition_1b, f"cond_{base_condition_1b}")
    cue_layout = cue_layout_for_condition(condition_name)

    direction_views = [
        build_direction_view(
            session_id=session_id,
            cell_index_1b=cell_index_1b,
            pf_cond_1b=forward_pf_cond_1b,
            x=x,
            idcond_t=idcond_t,
            fr_s_tx_u=fr_s_tx_u,
            fr_cx_u=fr_cx_u,
            fr_s_cx_u=fr_s_cx_u,
            n_cond_rmap=n_cond,
            condition_name_map=condition_name_map,
            overlap_lookup=overlap_lookup,
        ),
        build_direction_view(
            session_id=session_id,
            cell_index_1b=cell_index_1b,
            pf_cond_1b=backward_pf_cond_1b,
            x=x,
            idcond_t=idcond_t,
            fr_s_tx_u=fr_s_tx_u,
            fr_cx_u=fr_cx_u,
            fr_s_cx_u=fr_s_cx_u,
            n_cond_rmap=n_cond,
            condition_name_map=condition_name_map,
            overlap_lookup=overlap_lookup,
        ),
    ]
    focus_pf_cond_1b = int(condition_pf_1b)

    y_stack = np.concatenate(
        [view.mean_raw.ravel() for view in direction_views] + [view.mean_sm.ravel() for view in direction_views]
    )
    y_stack = y_stack[np.isfinite(y_stack)]
    y_upper = float(np.max(y_stack)) if y_stack.size else 1.0
    if (not np.isfinite(y_upper)) or (y_upper <= 0):
        y_upper = 1.0

    heat_parts = [view.heat[np.isfinite(view.heat)] for view in direction_views if view.heat.size > 0]
    heat_stack = np.concatenate(heat_parts) if heat_parts else np.array([], dtype=np.float64)
    heat_vmax = np.nan
    if heat_stack.size:
        heat_vmax = float(np.nanpercentile(heat_stack, float(heatmap_percentile)))
    if (not np.isfinite(heat_vmax)) or (heat_vmax <= 0):
        heat_vmax = float(np.nanmax(heat_stack)) if heat_stack.size else 1.0
    if (not np.isfinite(heat_vmax)) or (heat_vmax <= 0):
        heat_vmax = 1.0

    fig, axes = plt.subplots(2, 2, figsize=(13.4, 7.8), squeeze=False, sharex="col")
    total_trials = 0

    for col, view in enumerate(direction_views):
        ax_heat = axes[0, col]
        ax_mean = axes[1, col]
        total_trials += int(view.n_trials)
        is_focus = int(view.pf_cond_1b) == focus_pf_cond_1b
        title_prefix = "FOCUS" if is_focus else "Paired"

        if view.n_trials == 0:
            ax_heat.text(
                0.5,
                0.5,
                f"No trials found for rmap condition {view.condition_rmap_1b}",
                ha="center",
                va="center",
                transform=ax_heat.transAxes,
            )
            ax_heat.set_axis_off()
        else:
            im = ax_heat.imshow(
                view.heat,
                aspect="auto",
                origin="lower",
                interpolation="nearest",
                extent=[float(x[0]), float(x[-1]), 1, int(view.n_trials)],
                vmin=0.0,
                vmax=float(heat_vmax),
                cmap="viridis",
            )
            add_cue_overlays(ax_heat, cue_layout, add_labels=False)
            ax_heat.set_title(
                (
                    f"{title_prefix} {view.direction_label} | {view.condition_name} | "
                    f"rmap cond {view.condition_rmap_1b} | n_trials={view.n_trials}"
                ),
                fontsize=10,
            )
            ax_heat.set_xlabel("Position")
            ax_heat.set_ylabel("Trial")
            fig.colorbar(im, ax=ax_heat, fraction=0.046, pad=0.04).set_label("FR (smoothed)")

        add_cue_overlays(ax_mean, cue_layout, add_labels=True)
        ax_mean.plot(x, view.mean_raw, color="0.65", linewidth=1.5, label="mean raw")
        ax_mean.plot(x, view.mean_sm, color="#1f77b4", linewidth=2.0, label="mean smoothed")
        ax_mean.set_title(
            (
                f"{view.direction_label} | PF cond {view.pf_cond_1b} -> "
                f"rmap cond {view.condition_rmap_1b} [{view.cond_mapping}]"
            ),
            fontsize=10,
        )
        ax_mean.set_xlabel("Position")
        ax_mean.set_ylabel("FR")
        ax_mean.set_ylim(bottom=0.0, top=y_upper * 1.05)
        ax_mean.grid(alpha=0.25)
        ax_mean.legend(loc="upper right", fontsize=8)
        ax_mean.text(
            0.01,
            0.98,
            format_pf_state(view.state, label_a, label_b),
            ha="left",
            va="top",
            fontsize=8,
            transform=ax_mean.transAxes,
            bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none", "pad": 2.0},
        )
        if is_focus:
            for spine in ax_heat.spines.values():
                spine.set_linewidth(2.0)
                spine.set_edgecolor("#c0392b")
            for spine in ax_mean.spines.values():
                spine.set_linewidth(2.0)
                spine.set_edgecolor("#c0392b")

    cell_txt = f"{int(cell_id)}" if np.isfinite(cell_id) else "NA"
    suptitle = (
        f"{session_id} | class={row_class} | cell_id={cell_txt} | cell_idx={cell_index_1b} "
        f"(matched by {matched_by})\n"
        f"base condition={base_condition_1b} ({condition_name}) | focus PF cond={focus_pf_cond_1b}"
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

    return int(total_trials), int(n_cond)


def main() -> None:
    args = build_parser().parse_args()

    compare_dir = Path(args.compare_dir)
    ratemap_root = Path(args.ratemap_root)
    out_root = Path(args.out_root)
    working_root = Path.cwd()
    label_a, label_b = load_compare_labels(compare_dir)
    overlap_lookup = load_overlap_lookup(compare_dir)

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
    session_traj_path: dict[str, Path | None] = {}
    session_condition_name_map: dict[str, dict[int, str]] = {}

    for sid in d["session_id"].astype(str).drop_duplicates():
        p = find_rmap_file(ratemap_root, sid)
        session_rmap_path[sid] = p
        if p is not None:
            try:
                payload = load_rmap_payload(p)
                session_payload[sid] = payload
                traj_path = resolve_existing_path(payload.get("meta", {}).get("source_traj_npz"), working_root=working_root)
                session_traj_path[sid] = traj_path
                session_condition_name_map[sid] = load_condition_name_map(traj_path) if traj_path is not None else {}
            except Exception as exc:
                print(f"WARN: failed to load rmap for {sid}: {exc}")
                session_traj_path[sid] = None
                session_condition_name_map[sid] = {}
        else:
            session_traj_path[sid] = None
            session_condition_name_map[sid] = {}

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
        traj_path = session_traj_path.get(sid)
        condition_name_map = session_condition_name_map.get(sid, {})
        base_condition_1b = (int(cond_pf_1b) + 1) // 2
        condition_name = condition_name_map.get(base_condition_1b, f"cond_{base_condition_1b}")
        focus_direction = "backward" if (int(cond_pf_1b) % 2 == 0) else "forward"

        if (rmap_path is None) or (payload is None):
            missing_rows.append(
                {
                    "session_id": sid,
                    "cell_index_1b": cidx,
                    "cell_id": cid,
                    "condition_pf_1b": cond_pf_1b,
                    "base_condition_1b": base_condition_1b,
                    "condition_name": condition_name,
                    "row_class": row_class,
                    "reason": "missing_or_unreadable_rmap",
                    "rmap_path": "" if rmap_path is None else str(rmap_path),
                    "traj_path": "" if traj_path is None else str(traj_path),
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
                        "base_condition_1b": base_condition_1b,
                        "condition_name": condition_name,
                        "row_class": row_class,
                        "reason": "cell_not_found_in_rmap",
                        "rmap_path": str(rmap_path),
                        "traj_path": "" if traj_path is None else str(traj_path),
                    }
                )
                continue

            fr_cx_u = payload["fr_cx"][idx_u]
            fr_s_cx_u = payload["fr_s_cx"][idx_u]
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
                        "base_condition_1b": base_condition_1b,
                        "condition_name": condition_name,
                        "focus_direction": focus_direction,
                        "condition_rmap_1b": cond_rmap_1b,
                        "row_class": row_class,
                        "status": "exists_skip",
                        "out_png": str(out_png),
                        "rmap_path": str(rmap_path),
                        "traj_path": "" if traj_path is None else str(traj_path),
                        "matched_by": matched_by,
                        "cond_mapping": cond_mapping,
                        "reflected_even": bool(reflect_even),
                    }
                )
                continue

            n_trials, n_cond = save_row_plot(
                out_png=out_png,
                label_a=label_a,
                label_b=label_b,
                session_id=sid,
                row_class=row_class,
                condition_pf_1b=cond_pf_1b,
                cell_index_1b=cidx,
                cell_id=cid,
                matched_by=matched_by,
                x=payload["x"],
                idcond_t=payload["idcond_t"],
                fr_s_tx_u=fr_s_tx_u,
                fr_cx_u=fr_cx_u,
                fr_s_cx_u=fr_s_cx_u,
                rmap_path=rmap_path,
                meta=payload["meta"],
                condition_name_map=condition_name_map,
                overlap_lookup=overlap_lookup,
                heatmap_percentile=float(args.heatmap_percentile),
                dpi=int(args.dpi),
            )

            plot_rows.append(
                {
                    "session_id": sid,
                    "cell_index_1b": cidx,
                    "cell_id": cid,
                    "condition_pf_1b": cond_pf_1b,
                    "base_condition_1b": base_condition_1b,
                    "condition_name": condition_name,
                    "focus_direction": focus_direction,
                    "condition_rmap_1b": cond_rmap_1b,
                    "row_class": row_class,
                    "status": "plotted",
                    "out_png": str(out_png),
                    "rmap_path": str(rmap_path),
                    "traj_path": "" if traj_path is None else str(traj_path),
                    "matched_by": matched_by,
                    "cond_mapping": cond_mapping,
                    "reflected_even": bool(reflect_even),
                    "n_trials_plotted_total": int(n_trials),
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
                    "base_condition_1b": base_condition_1b,
                    "condition_name": condition_name,
                    "row_class": row_class,
                    "reason": "plot_error",
                    "rmap_path": str(rmap_path),
                    "traj_path": "" if traj_path is None else str(traj_path),
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
        "base_condition_1b",
        "condition_name",
        "focus_direction",
        "condition_rmap_1b",
        "row_class",
        "status",
        "out_png",
        "rmap_path",
        "traj_path",
        "matched_by",
        "cond_mapping",
        "reflected_even",
        "n_trials_plotted_total",
        "n_conditions_in_rmap",
    ]
    missing_cols = [
        "session_id",
        "cell_index_1b",
        "cell_id",
        "condition_pf_1b",
        "base_condition_1b",
        "condition_name",
        "row_class",
        "reason",
        "rmap_path",
        "traj_path",
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
            "label_a": label_a,
            "label_b": label_b,
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
