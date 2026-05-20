from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PLOT_X_MIN = 0.0
PLOT_X_MAX = 100.0
HEATMAP_CMAP = "magma"


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
class TrialFamily:
    key: str
    title: str
    member_conditions: tuple[str, ...]
    cue_layout: CueLayout


@dataclass(frozen=True)
class DirectionPlotData:
    direction_key: str
    direction_label: str
    trial_indices_0b: np.ndarray
    original_condition_counts: tuple[tuple[str, int], ...]
    condition_blocks: tuple[tuple[str, int, int], ...]
    heat: np.ndarray
    mean_raw: np.ndarray
    mean_sm: np.ndarray


PRETTY_CONDITION_NAMES: dict[str, str] = {
    "PO": "PO",
    "PO2": "PO2",
    "PO3": "PO3",
    "PONM": "POnM",
    "PNO": "PNO",
    "POM": "POM",
    "POMB": "POMb",
    "PO_NOTREE": "PO_NOTREE",
}

STANDARD_CUE_LAYOUT = CueLayout(
    name="object_standard",
    cue_rich=(13.0, 43.0),
    cue_poor=(43.0, 81.0),
    object_zone=(81.0, 96.0),
    object_centers=(20.0, 36.0, 88.0),
)

PNO_CUE_LAYOUT = CueLayout(
    name="pno",
    cue_rich=(13.0, 43.0),
    cue_poor=(43.0, 81.0),
    object_zone=(81.0, 96.0),
    object_centers=tuple(),
    note="No objects",
)

POM_CUE_LAYOUT = CueLayout(
    name="pom_family",
    cue_rich=(13.0, 28.0),
    cue_poor=(28.0, 57.0),
    object_zone=(57.0, 96.0),
    object_centers=(20.0, 64.0, 88.0),
    moved_object_span=None,
    note=None,
)

TRIAL_FAMILIES: tuple[TrialFamily, ...] = (
    TrialFamily(
        key="object_standard",
        title="PO / PO2 / PO3 / POnM",
        member_conditions=("PO", "PO2", "PO3", "PONM"),
        cue_layout=STANDARD_CUE_LAYOUT,
    ),
    TrialFamily(
        key="pno",
        title="PNO",
        member_conditions=("PNO",),
        cue_layout=PNO_CUE_LAYOUT,
    ),
    TrialFamily(
        key="pom_family",
        title="POM / POMb",
        member_conditions=("POM", "POMB"),
        cue_layout=POM_CUE_LAYOUT,
    ),
)

IGNORED_CONDITIONS = {"PO_NOTREE"}


def parse_meta_json(raw: np.ndarray) -> dict[str, Any]:
    try:
        txt = raw.tobytes().decode("utf-8", errors="ignore")
        obj = json.loads(txt)
        if isinstance(obj, dict):
            return obj
    except Exception:
        pass
    return {}


def canonical_condition_name(name: str | None) -> str:
    if name is None:
        return ""
    return str(name).strip().upper()


def pretty_condition_name(name: str | None) -> str:
    cname = canonical_condition_name(name)
    return PRETTY_CONDITION_NAMES.get(cname, cname or "unknown")


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


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Plot grouped-condition ratemaps for cells that have place fields in both "
            "PF comparison methods."
        )
    )
    ap.add_argument(
        "--compare_dir",
        type=str,
        default="results/tables/pf_npz_compare/random_poisson_vs_circular_shift",
        help="Folder produced by compare_pf_npz_methods.py.",
    )
    ap.add_argument(
        "--unit_summary_csv",
        type=str,
        default="",
        help="Optional explicit unit_summary.csv path. Overrides --compare_dir when provided.",
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
        default="results/figures/pf_both_methods_grouped_conditions",
    )
    ap.add_argument(
        "--max_units",
        type=int,
        default=0,
        help="Optional cap on number of selected cells to plot (0 = all).",
    )
    ap.add_argument(
        "--heatmap_percentile",
        type=float,
        default=99.0,
        help="Upper percentile for heatmap color scaling within each grouped-family figure.",
    )
    ap.add_argument("--dpi", type=int, default=140)
    ap.add_argument("--overwrite", action="store_true")
    return ap


def resolve_unit_summary_csv(compare_dir: Path, unit_summary_csv: str) -> Path:
    if str(unit_summary_csv).strip():
        p = Path(unit_summary_csv)
        if not p.exists():
            raise FileNotFoundError(f"Missing --unit_summary_csv: {p}")
        return p
    p = compare_dir / "unit_summary.csv"
    if not p.exists():
        raise FileNotFoundError(
            f"Could not find unit_summary.csv in {compare_dir}. "
            "Run compare_pf_npz_methods.py first or pass --unit_summary_csv."
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
            "xbin_centers",
            "xbin_edges",
            "idcond_t",
            "rmap__fr_tx_ux",
            "rmap__fr_s_tx_ux",
        ]
        miss = [k for k in req if k not in z.files]
        if miss:
            raise KeyError(f"missing keys in {rmap_path}: {miss}")

        return {
            "cell_ids": z["cell_ids"].astype(np.int64, copy=False),
            "x": z["xbin_centers"].astype(np.float64, copy=False),
            "x_edges": z["xbin_edges"].astype(np.float64, copy=False),
            "idcond_t": z["idcond_t"].astype(np.int64, copy=False),
            "fr_tx": z["rmap__fr_tx_ux"].astype(np.float64, copy=False),
            "fr_s_tx": z["rmap__fr_s_tx_ux"].astype(np.float64, copy=False),
            "meta": parse_meta_json(z["meta_json"]) if "meta_json" in z.files else {},
        }


def load_trial_table(traj_path: Path) -> pd.DataFrame:
    if not traj_path.exists():
        raise FileNotFoundError(f"Missing trajdata file: {traj_path}")

    with np.load(traj_path, allow_pickle=False) as z:
        if "traj__WB" not in z.files:
            raise KeyError(f"{traj_path} missing key traj__WB")

        if "traj__condition" in z.files:
            cond_raw = np.asarray(z["traj__condition"]).astype(str).ravel()
        elif "traj__condition__json" in z.files:
            raw = z["traj__condition__json"].tobytes().decode("utf-8", errors="ignore")
            cond_raw = np.asarray(json.loads(raw), dtype=object).astype(str).ravel()
        else:
            raise KeyError(f"{traj_path} missing traj__condition / traj__condition__json")

        wb_raw = np.asarray(z["traj__WB"]).astype(str).ravel()

    if cond_raw.size != wb_raw.size:
        raise ValueError(
            f"traj condition/WB size mismatch in {traj_path}: {cond_raw.size} vs {wb_raw.size}"
        )

    cond_clean = np.asarray([str(v).strip() for v in cond_raw], dtype=object)
    wb_clean = np.asarray([str(v).strip().upper() for v in wb_raw], dtype=object)
    bad_wb = sorted(set(wb_clean[~np.isin(wb_clean, ["W", "B"])].tolist()))
    if bad_wb:
        raise ValueError(f"Unsupported WB values in {traj_path}: {bad_wb}")

    d = pd.DataFrame(
        {
            "trial_index_0b": np.arange(cond_clean.size, dtype=np.int64),
            "condition_raw": cond_clean.astype(str),
            "condition_canon": [canonical_condition_name(v) for v in cond_clean.tolist()],
            "condition_pretty": [pretty_condition_name(v) for v in cond_clean.tolist()],
            "wb": wb_clean.astype(str),
        }
    )
    return d


def load_selected_units(unit_summary_csv: Path) -> pd.DataFrame:
    d = pd.read_csv(unit_summary_csv)
    required = [
        "session_id",
        "cell_index_1b",
        "any_pf_a",
        "any_pf_b",
        "present_a",
        "present_b",
    ]
    missing = [c for c in required if c not in d.columns]
    if missing:
        raise KeyError(f"Missing required columns in {unit_summary_csv}: {missing}")

    d = d.copy()
    d["session_id"] = d["session_id"].astype(str)
    d["cell_index_1b"] = pd.to_numeric(d["cell_index_1b"], errors="coerce").astype("Int64")
    d["any_pf_a"] = d["any_pf_a"].astype(bool)
    d["any_pf_b"] = d["any_pf_b"].astype(bool)
    d["present_a"] = d["present_a"].astype(bool)
    d["present_b"] = d["present_b"].astype(bool)

    if "cell_id" in d.columns:
        cell_id_series = d["cell_id"]
    elif "cell_id_a" in d.columns:
        cell_id_series = d["cell_id_a"]
    elif "cell_id_b" in d.columns:
        cell_id_series = d["cell_id_b"]
    else:
        cell_id_series = pd.Series(dtype=float)

    d["cell_id"] = pd.to_numeric(cell_id_series, errors="coerce")
    d = d.dropna(subset=["session_id", "cell_index_1b"]).copy()
    d["cell_index_1b"] = d["cell_index_1b"].astype(np.int64)

    keep = d["present_a"] & d["present_b"] & d["any_pf_a"] & d["any_pf_b"]
    out = (
        d.loc[keep, ["session_id", "cell_index_1b", "cell_id", "any_pf_a", "any_pf_b"]]
        .drop_duplicates(subset=["session_id", "cell_index_1b"], keep="first")
        .sort_values(["session_id", "cell_index_1b"], kind="stable")
        .reset_index(drop=True)
    )
    return out


def add_cue_overlays(
    ax: plt.Axes,
    layout: CueLayout,
    *,
    add_labels: bool,
) -> None:
    ax.axvspan(layout.cue_rich[0], layout.cue_rich[1], color="#f4c095", alpha=0.22, zorder=1)
    ax.axvspan(layout.cue_poor[0], layout.cue_poor[1], color="#d7ebba", alpha=0.18, zorder=1)
    ax.axvspan(layout.object_zone[0], layout.object_zone[1], color="#a8cbe6", alpha=0.22, zorder=1)
    if layout.moved_object_span is not None:
        ax.axvspan(
            layout.moved_object_span[0],
            layout.moved_object_span[1],
            color="#f7a1a1",
            alpha=0.28,
            zorder=0.1,
        )
    for center in layout.object_centers:
        ax.axvline(center, color="#8b1e3f", linestyle="--", linewidth=1.1, alpha=0.9, zorder=3)
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


def format_condition_counts(counts: tuple[tuple[str, int], ...]) -> str:
    if not counts:
        return "No trials"
    return " | ".join(f"{name}={count}" for name, count in counts)


def add_condition_block_annotations(
    ax: plt.Axes,
    *,
    x_min: float,
    x_max: float,
    condition_blocks: tuple[tuple[str, int, int], ...],
) -> None:
    if not condition_blocks:
        return

    x_text = float(x_min) + (float(x_max) - float(x_min)) * 0.012
    for i, (name, start_row_1b, end_row_1b) in enumerate(condition_blocks):
        if i > 0:
            y_sep = float(start_row_1b) - 0.5
            ax.axhline(y_sep, color="white", linewidth=1.2, alpha=0.9, zorder=4)
        y_mid = 0.5 * (float(start_row_1b) + float(end_row_1b))
        ax.text(
            x_text,
            y_mid,
            str(name),
            ha="left",
            va="center",
            fontsize=7.5,
            color="black",
            bbox={"facecolor": "white", "alpha": 0.72, "edgecolor": "none", "pad": 1.5},
            zorder=5,
        )


def select_family_trials(trial_table: pd.DataFrame, family: TrialFamily, direction_key: str) -> pd.DataFrame:
    order = {name: i for i, name in enumerate(family.member_conditions)}
    keep = trial_table["condition_canon"].isin(family.member_conditions) & (trial_table["wb"] == direction_key)
    out = trial_table.loc[keep].copy()
    if out.empty:
        return out
    out = out[~out["condition_canon"].isin(IGNORED_CONDITIONS)].copy()
    if out.empty:
        return out
    out["family_rank"] = out["condition_canon"].map(order).astype(np.int64)
    out = out.sort_values(["family_rank", "trial_index_0b"], kind="stable").reset_index(drop=True)
    return out


def build_direction_plot_data(
    *,
    family: TrialFamily,
    direction_key: str,
    trial_table: pd.DataFrame,
    fr_tx_u: np.ndarray,
    fr_s_tx_u: np.ndarray,
) -> DirectionPlotData:
    subset = select_family_trials(trial_table, family, direction_key)
    n_bins = int(fr_tx_u.shape[1])
    direction_label = "Forward (W)" if direction_key == "W" else "Backward (B)"

    if subset.empty:
        return DirectionPlotData(
            direction_key=direction_key,
            direction_label=direction_label,
            trial_indices_0b=np.empty((0,), dtype=np.int64),
            original_condition_counts=tuple(),
            condition_blocks=tuple(),
            heat=np.empty((0, n_bins), dtype=np.float64),
            mean_raw=np.full(n_bins, np.nan, dtype=np.float64),
            mean_sm=np.full(n_bins, np.nan, dtype=np.float64),
        )

    idx = subset["trial_index_0b"].to_numpy(dtype=np.int64, copy=True)
    heat = np.asarray(fr_s_tx_u[idx, :], dtype=np.float64)
    mean_raw = np.asarray(np.nanmean(fr_tx_u[idx, :], axis=0), dtype=np.float64)
    mean_sm = np.asarray(np.nanmean(fr_s_tx_u[idx, :], axis=0), dtype=np.float64)

    counts: list[tuple[str, int]] = []
    blocks: list[tuple[str, int, int]] = []
    row_start = 1
    for cond in family.member_conditions:
        n = int((subset["condition_canon"] == cond).sum())
        if n > 0:
            pretty = pretty_condition_name(cond)
            counts.append((pretty, n))
            blocks.append((pretty, row_start, row_start + n - 1))
            row_start += n

    if direction_key == "B":
        heat = heat[:, ::-1]
        mean_raw = mean_raw[::-1]
        mean_sm = mean_sm[::-1]

    return DirectionPlotData(
        direction_key=direction_key,
        direction_label=direction_label,
        trial_indices_0b=idx,
        original_condition_counts=tuple(counts),
        condition_blocks=tuple(blocks),
        heat=heat,
        mean_raw=mean_raw,
        mean_sm=mean_sm,
    )


def save_family_plot(
    *,
    out_png: Path,
    family: TrialFamily,
    label_a: str,
    label_b: str,
    session_id: str,
    cell_index_1b: int,
    cell_id: float,
    matched_by: str,
    x: np.ndarray,
    x_edges: np.ndarray,
    fr_tx_u: np.ndarray,
    fr_s_tx_u: np.ndarray,
    trial_table: pd.DataFrame,
    rmap_path: Path,
    traj_path: Path,
    meta: dict[str, Any],
    heatmap_percentile: float,
    dpi: int,
) -> tuple[int, int]:
    directions = [
        build_direction_plot_data(
            family=family,
            direction_key="W",
            trial_table=trial_table,
            fr_tx_u=fr_tx_u,
            fr_s_tx_u=fr_s_tx_u,
        ),
        build_direction_plot_data(
            family=family,
            direction_key="B",
            trial_table=trial_table,
            fr_tx_u=fr_tx_u,
            fr_s_tx_u=fr_s_tx_u,
        ),
    ]

    mean_parts = [
        arr[np.isfinite(arr)]
        for view in directions
        for arr in (view.mean_raw, view.mean_sm)
        if np.isfinite(arr).any()
    ]
    mean_stack = np.concatenate(mean_parts) if mean_parts else np.array([], dtype=np.float64)
    y_upper = float(np.nanmax(mean_stack)) if mean_stack.size else 1.0
    if (not np.isfinite(y_upper)) or (y_upper <= 0):
        y_upper = 1.0

    heat_parts = [view.heat[np.isfinite(view.heat)] for view in directions if view.heat.size > 0]
    heat_stack = np.concatenate(heat_parts) if heat_parts else np.array([], dtype=np.float64)
    heat_vmax = np.nan
    if heat_stack.size:
        heat_vmax = float(np.nanpercentile(heat_stack, float(heatmap_percentile)))
    if (not np.isfinite(heat_vmax)) or (heat_vmax <= 0):
        heat_vmax = float(np.nanmax(heat_stack)) if heat_stack.size else 1.0
    if (not np.isfinite(heat_vmax)) or (heat_vmax <= 0):
        heat_vmax = 1.0

    fig = plt.figure(figsize=(13.6, 7.9))
    gs = fig.add_gridspec(
        2,
        4,
        width_ratios=[1.0, 0.026, 1.0, 0.026],
        height_ratios=[1.0, 1.0],
        wspace=0.08,
        hspace=0.22,
    )
    ax_heat_left = fig.add_subplot(gs[0, 0])
    ax_heat_right = fig.add_subplot(gs[0, 2])
    ax_mean_left = fig.add_subplot(gs[1, 0], sharex=ax_heat_left)
    ax_mean_right = fig.add_subplot(gs[1, 2], sharex=ax_heat_right)
    cax_left = fig.add_subplot(gs[0, 1])
    cax_right = fig.add_subplot(gs[0, 3])
    axes = ((ax_heat_left, ax_heat_right), (ax_mean_left, ax_mean_right))
    caxes = (cax_left, cax_right)
    total_trials = 0
    for col, view in enumerate(directions):
        ax_heat = axes[0][col]
        ax_mean = axes[1][col]
        cax = caxes[col]
        total_trials += int(view.trial_indices_0b.size)

        if view.trial_indices_0b.size == 0:
            ax_heat.set_xlim(float(x[0]), float(x[-1]))
            ax_heat.set_ylim(0.0, 1.0)
            ax_heat.set_title(f"{view.direction_label} | n_trials=0", fontsize=10)
            ax_heat.set_xlabel("Position")
            ax_heat.set_ylabel("Trial")
            ax_heat.text(
                0.5,
                0.5,
                f"No {view.direction_key.lower()} trials in {family.title}",
                ha="center",
                va="center",
                transform=ax_heat.transAxes,
            )
            ax_heat.set_yticks([])
        else:
            im = ax_heat.imshow(
                view.heat,
                aspect="auto",
                origin="lower",
                interpolation="nearest",
                extent=[float(x_edges[0]), float(x_edges[-1]), 0.5, float(int(view.trial_indices_0b.size) + 0.5)],
                vmin=0.0,
                vmax=float(heat_vmax),
                cmap=HEATMAP_CMAP,
            )
            ax_heat.set_title(
                f"{view.direction_label} | n_trials={int(view.trial_indices_0b.size)}",
                fontsize=10,
            )
            ax_heat.set_xlabel("Position")
            ax_heat.set_ylabel("Trial")
        add_cue_overlays(ax_heat, family.cue_layout, add_labels=False)
        ax_heat.set_xlim(PLOT_X_MIN, PLOT_X_MAX)
        if view.trial_indices_0b.size > 0:
            add_condition_block_annotations(
                ax_heat,
                x_min=float(x[0]),
                x_max=float(x[-1]),
                condition_blocks=view.condition_blocks,
            )
            fig.colorbar(im, cax=cax).set_label("FR (smoothed)")
        else:
            cax.set_visible(False)
        ax_heat.text(
            0.01,
            0.98,
            format_condition_counts(view.original_condition_counts),
            ha="left",
            va="top",
            fontsize=8,
            transform=ax_heat.transAxes,
            bbox={"facecolor": "white", "alpha": 0.78, "edgecolor": "none", "pad": 2.0},
        )

        add_cue_overlays(ax_mean, family.cue_layout, add_labels=True)
        if view.trial_indices_0b.size > 0:
            ax_mean.plot(x, view.mean_raw, color="0.65", linewidth=1.5, label="mean raw")
            ax_mean.plot(x, view.mean_sm, color="#1f77b4", linewidth=2.0, label="mean smoothed")
            ax_mean.legend(loc="upper right", fontsize=8)
        else:
            ax_mean.text(
                0.5,
                0.5,
                f"No {view.direction_key.lower()} trials",
                ha="center",
                va="center",
                transform=ax_mean.transAxes,
            )
        ax_mean.set_title(f"{view.direction_label} mean FR", fontsize=10)
        ax_mean.set_xlabel("Position")
        ax_mean.set_ylabel("FR")
        ax_mean.set_xlim(PLOT_X_MIN, PLOT_X_MAX)
        ax_mean.set_ylim(bottom=0.0, top=y_upper * 1.05)
        ax_mean.grid(alpha=0.25)
        ax_mean.text(
            0.01,
            0.98,
            format_condition_counts(view.original_condition_counts),
            ha="left",
            va="top",
            fontsize=8,
            transform=ax_mean.transAxes,
            bbox={"facecolor": "white", "alpha": 0.78, "edgecolor": "none", "pad": 2.0},
        )

    cell_txt = f"{int(cell_id)}" if np.isfinite(cell_id) else "NA"
    included = ", ".join(pretty_condition_name(c) for c in family.member_conditions)
    fig.suptitle(
        (
            f"{session_id} | cell_id={cell_txt} | cell_idx={cell_index_1b} "
            f"(matched by {matched_by})\n"
            f"family={family.title} | grouped labels={included} | selected because PF in both: {label_a} + {label_b}"
        ),
        fontsize=10,
        y=0.99,
    )

    smooth = meta.get("smooth_sigma_bins")
    xrem = meta.get("xbin_rem")
    info = f"rmap={rmap_path.name} | traj={traj_path.name}"
    if smooth is not None:
        info += f" | smooth_sigma_bins={smooth}"
    if xrem is not None:
        info += f" | xbin_rem={xrem}"
    fig.text(0.01, 0.01, info, fontsize=8, ha="left", va="bottom")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(left=0.06, right=0.97, bottom=0.10, top=0.88)
    fig.savefig(out_png, dpi=int(dpi), bbox_inches="tight")
    plt.close(fig)
    return int(directions[0].trial_indices_0b.size), int(directions[1].trial_indices_0b.size)


def main() -> None:
    args = build_parser().parse_args()

    compare_dir = Path(args.compare_dir)
    ratemap_root = Path(args.ratemap_root)
    out_root = Path(args.out_root)
    working_root = Path.cwd()
    label_a, label_b = load_compare_labels(compare_dir)

    if not ratemap_root.exists():
        raise FileNotFoundError(f"Missing ratemap root: {ratemap_root}")

    unit_summary_csv = resolve_unit_summary_csv(compare_dir, args.unit_summary_csv)
    units = load_selected_units(unit_summary_csv)
    if units.empty:
        print(f"No cells with any_pf_a && any_pf_b in: {unit_summary_csv}")
        return

    if int(args.max_units) > 0:
        units = units.head(int(args.max_units)).copy()

    run_tag = compare_dir.name if compare_dir.exists() else "manual_compare"
    out_dir = out_root / run_tag
    out_dir.mkdir(parents=True, exist_ok=True)

    session_rmap_path: dict[str, Path | None] = {}
    session_payload: dict[str, dict[str, Any]] = {}
    session_traj_path: dict[str, Path | None] = {}
    session_trial_table: dict[str, pd.DataFrame] = {}

    for sid in units["session_id"].astype(str).drop_duplicates():
        rmap_path = find_rmap_file(ratemap_root, sid)
        session_rmap_path[sid] = rmap_path
        if rmap_path is None:
            session_traj_path[sid] = None
            continue
        try:
            payload = load_rmap_payload(rmap_path)
            session_payload[sid] = payload
            traj_path = resolve_existing_path(payload.get("meta", {}).get("source_traj_npz"), working_root=working_root)
            session_traj_path[sid] = traj_path
            if traj_path is None:
                print(f"WARN: missing source_traj_npz for {sid}")
                continue
            trial_table = load_trial_table(traj_path)
            n_trials_rmap = int(payload["fr_tx"].shape[1])
            if int(len(trial_table)) != n_trials_rmap:
                raise ValueError(
                    f"trial count mismatch for {sid}: traj={len(trial_table)} vs rmap={n_trials_rmap}"
                )
            session_trial_table[sid] = trial_table
        except Exception as exc:
            print(f"WARN: failed to prepare session {sid}: {exc}")

    plot_rows: list[dict[str, Any]] = []
    missing_rows: list[dict[str, Any]] = []

    n_units = int(len(units))
    n_requested_families = 0
    for i, row in enumerate(units.itertuples(index=False), start=1):
        sid = str(row.session_id)
        cidx = int(row.cell_index_1b)
        cid = float(row.cell_id) if pd.notna(row.cell_id) else float("nan")
        rmap_path = session_rmap_path.get(sid)
        payload = session_payload.get(sid)
        traj_path = session_traj_path.get(sid)
        trial_table = session_trial_table.get(sid)

        if (rmap_path is None) or (payload is None) or (traj_path is None) or (trial_table is None):
            missing_rows.append(
                {
                    "session_id": sid,
                    "cell_index_1b": cidx,
                    "cell_id": cid,
                    "reason": "missing_session_payload",
                    "rmap_path": "" if rmap_path is None else str(rmap_path),
                    "traj_path": "" if traj_path is None else str(traj_path),
                }
            )
            continue

        try:
            idx_u, matched_by = cell_axis_index(payload["cell_ids"], cid, cidx)
            if idx_u is None:
                missing_rows.append(
                    {
                        "session_id": sid,
                        "cell_index_1b": cidx,
                        "cell_id": cid,
                        "reason": "cell_not_found_in_rmap",
                        "rmap_path": str(rmap_path),
                        "traj_path": str(traj_path),
                    }
                )
                continue

            fr_tx_u = payload["fr_tx"][idx_u]
            fr_s_tx_u = payload["fr_s_tx"][idx_u]
            x = payload["x"]

            cell_label = f"{int(cid)}" if np.isfinite(cid) else f"idx_{cidx:04d}"
            unit_out_dir = out_dir / sid / f"cell_{cell_label}"
            for family in TRIAL_FAMILIES:
                fam_mask = trial_table["condition_canon"].isin(family.member_conditions)
                fam_mask &= ~trial_table["condition_canon"].isin(IGNORED_CONDITIONS)
                if not bool(fam_mask.any()):
                    plot_rows.append(
                        {
                            "session_id": sid,
                            "cell_index_1b": cidx,
                            "cell_id": cid,
                            "family_key": family.key,
                            "family_title": family.title,
                            "status": "skip_no_trials_in_session",
                            "out_png": "",
                            "rmap_path": str(rmap_path),
                            "traj_path": str(traj_path),
                            "matched_by": matched_by,
                        }
                    )
                    continue

                n_requested_families += 1
                out_png = unit_out_dir / f"{family.key}.png"
                if out_png.exists() and not bool(args.overwrite):
                    plot_rows.append(
                        {
                            "session_id": sid,
                            "cell_index_1b": cidx,
                            "cell_id": cid,
                            "family_key": family.key,
                            "family_title": family.title,
                            "status": "exists_skip",
                            "out_png": str(out_png),
                            "rmap_path": str(rmap_path),
                            "traj_path": str(traj_path),
                            "matched_by": matched_by,
                        }
                    )
                    continue

                n_forward, n_backward = save_family_plot(
                    out_png=out_png,
                    family=family,
                    label_a=label_a,
                    label_b=label_b,
                    session_id=sid,
                    cell_index_1b=cidx,
                    cell_id=cid,
                    matched_by=matched_by,
                    x=x,
                    x_edges=payload["x_edges"],
                    fr_tx_u=fr_tx_u,
                    fr_s_tx_u=fr_s_tx_u,
                    trial_table=trial_table,
                    rmap_path=rmap_path,
                    traj_path=traj_path,
                    meta=payload["meta"],
                    heatmap_percentile=float(args.heatmap_percentile),
                    dpi=int(args.dpi),
                )
                plot_rows.append(
                    {
                        "session_id": sid,
                        "cell_index_1b": cidx,
                        "cell_id": cid,
                        "family_key": family.key,
                        "family_title": family.title,
                        "status": "plotted",
                        "out_png": str(out_png),
                        "rmap_path": str(rmap_path),
                        "traj_path": str(traj_path),
                        "matched_by": matched_by,
                        "n_trials_forward": int(n_forward),
                        "n_trials_backward": int(n_backward),
                    }
                )

        except Exception as exc:
            missing_rows.append(
                {
                    "session_id": sid,
                    "cell_index_1b": cidx,
                    "cell_id": cid,
                    "reason": "plot_error",
                    "rmap_path": str(rmap_path),
                    "traj_path": str(traj_path),
                    "error": str(exc),
                }
            )

        if (i % 25 == 0) or (i == n_units):
            print(f"progress: {i}/{n_units} selected cells")

    selected_out = units.copy()
    selected_out.to_csv(out_dir / "selected_units.csv", index=False)

    plotted_cols = [
        "session_id",
        "cell_index_1b",
        "cell_id",
        "family_key",
        "family_title",
        "status",
        "out_png",
        "rmap_path",
        "traj_path",
        "matched_by",
        "n_trials_forward",
        "n_trials_backward",
    ]
    missing_cols = [
        "session_id",
        "cell_index_1b",
        "cell_id",
        "reason",
        "rmap_path",
        "traj_path",
        "error",
    ]

    plotted = pd.DataFrame(plot_rows, columns=plotted_cols)
    missing = pd.DataFrame(missing_rows, columns=missing_cols)
    plotted.to_csv(out_dir / "plotted_families.csv", index=False)
    missing.to_csv(out_dir / "plot_missing_or_errors.csv", index=False)

    summary = {
        "inputs": {
            "compare_dir": str(compare_dir),
            "unit_summary_csv": str(unit_summary_csv),
            "ratemap_root": str(ratemap_root),
            "label_a": label_a,
            "label_b": label_b,
            "max_units": int(args.max_units),
            "heatmap_percentile": float(args.heatmap_percentile),
        },
        "counts": {
            "n_selected_units": int(len(units)),
            "n_requested_families": int(n_requested_families),
            "n_plotted": int((plotted["status"] == "plotted").sum()) if not plotted.empty else 0,
            "n_exists_skip": int((plotted["status"] == "exists_skip").sum()) if not plotted.empty else 0,
            "n_skip_no_trials_in_session": int((plotted["status"] == "skip_no_trials_in_session").sum())
            if not plotted.empty
            else 0,
            "n_missing_or_errors": int(len(missing)),
        },
        "families": [
            {
                "key": fam.key,
                "title": fam.title,
                "member_conditions": [pretty_condition_name(v) for v in fam.member_conditions],
            }
            for fam in TRIAL_FAMILIES
        ],
        "out_dir": str(out_dir),
    }
    with open(out_dir / "summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"Wrote grouped-condition figures to: {out_dir}")
    print(
        f"selected_units={summary['counts']['n_selected_units']} | "
        f"requested_families={summary['counts']['n_requested_families']} | "
        f"plotted={summary['counts']['n_plotted']} | "
        f"skipped={summary['counts']['n_exists_skip']} | "
        f"missing_or_errors={summary['counts']['n_missing_or_errors']}"
    )


if __name__ == "__main__":
    main()
