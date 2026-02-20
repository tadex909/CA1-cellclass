from __future__ import annotations

import argparse
from functools import lru_cache
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_csv_list(raw: str | None) -> list[str]:
    if not raw:
        return []
    return [x.strip() for x in raw.split(",") if x.strip()]


@lru_cache(maxsize=1024)
def find_interim_npz(interim_root: str, session_id: str) -> Optional[Path]:
    root = Path(interim_root)
    hits = sorted(root.rglob(f"{session_id}_*.npz"))
    return hits[0] if hits else None


@lru_cache(maxsize=1024)
def load_waveform_arrays(npz_path: str) -> tuple[np.ndarray, np.ndarray]:
    with np.load(npz_path, allow_pickle=False) as z:
        cell_ids = np.asarray(z["allcel__id_cel"]).astype(np.int64).ravel()
        bestw = np.asarray(z["allcel__bestwaveform"]).astype(np.float64)
    return cell_ids, bestw


def extract_typical_waveform(
    npz_path: Path,
    cell_id: int,
    fs_hz: float,
) -> Optional[tuple[np.ndarray, np.ndarray]]:
    """
    Return (time_ms, waveform) for one neuron from allcel__bestwaveform.
    """
    try:
        cell_ids, bestw = load_waveform_arrays(str(npz_path))
    except Exception:
        return None

    idx = np.where(cell_ids == int(cell_id))[0]
    if idx.size == 0:
        return None

    cell_idx = int(idx[0])
    wf: Optional[np.ndarray] = None

    # Common expected shapes:
    # (n_samples, n_cells) or (n_cells, n_samples)
    if bestw.ndim == 2:
        if bestw.shape[1] == cell_ids.size:
            wf = bestw[:, cell_idx]
        elif bestw.shape[0] == cell_ids.size:
            wf = bestw[cell_idx, :]
    # Fallback if an unexpected 3D layout appears: pick by cell axis where possible.
    elif bestw.ndim == 3:
        if bestw.shape[2] == cell_ids.size:
            wf = np.asarray(bestw[:, :, cell_idx]).reshape(-1)
        elif bestw.shape[0] == cell_ids.size:
            wf = np.asarray(bestw[cell_idx, :, :]).reshape(-1)

    if wf is None:
        return None

    wf = np.asarray(wf, dtype=np.float64).ravel()
    if wf.size == 0:
        return None

    t_ms = np.arange(wf.size, dtype=np.float64) * (1000.0 / fs_hz)
    return t_ms, wf


def build_population_cache(comparison_root: Path) -> dict[str, pd.DataFrame]:
    out: dict[str, pd.DataFrame] = {}
    for age_dir in sorted([p for p in comparison_root.iterdir() if p.is_dir()]):
        age = age_dir.name
        p = age_dir / f"{age}_gmm2_vs_type_u.parquet"
        if p.exists():
            out[age] = pd.read_parquet(p)
    return out


def get_neuron_xy(
    pop: Optional[pd.DataFrame],
    row: pd.Series,
) -> tuple[float, float]:
    """
    Get (spk_duration_ms, fr_hz) for the highlighted neuron.
    Prefer explicit values in row, fallback to lookup in population table.
    """
    x = row.get("spk_duration_ms", np.nan)
    y = row.get("fr_hz", np.nan)
    if pd.notna(x) and pd.notna(y):
        return float(x), float(y)

    if pop is None:
        return float("nan"), float("nan")

    m = pop[
        (pop["session_id"].astype(str) == str(row["session_id"]))
        & (pop["cell_id"].astype(int) == int(row["cell_id"]))
    ]
    if m.empty:
        return float("nan"), float("nan")
    return float(m.iloc[0]["spk_duration_ms"]), float(m.iloc[0]["fr_hz"])


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Create per-neuron review artifacts for disagreements between type_u and GMM labels."
        )
    )
    ap.add_argument("--comparison_root", type=str, default="results/type_u_comparison")
    ap.add_argument("--interim_root", type=str, default="data/interim")
    ap.add_argument("--out_root", type=str, default="results/type_u_comparison/review")
    ap.add_argument("--ages", type=str, default="")
    ap.add_argument("--top_n", type=int, default=0, help="0 means all discrepant neurons")
    ap.add_argument("--waveform_fs_hz", type=float, default=25000.0)
    args = ap.parse_args()

    comparison_root = Path(args.comparison_root)
    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    discrepant_path = comparison_root / "discrepant_neurons.parquet"
    if not discrepant_path.exists():
        raise FileNotFoundError(f"Missing {discrepant_path}")

    d = pd.read_parquet(discrepant_path)
    if d.empty:
        print("No discrepant neurons found.")
        return

    ages = set(parse_csv_list(args.ages))
    if ages:
        d = d[d["age_group"].isin(ages)].copy()
    if args.top_n > 0:
        d = d.head(args.top_n).copy()

    pop_by_age = build_population_cache(comparison_root)

    review_rows: list[dict] = []
    for _, row in d.iterrows():
        age = str(row["age_group"])
        session_id = str(row["session_id"])
        cell_id = int(row["cell_id"])
        unit_uid = str(row.get("unit_uid", f"{session_id}__cell{cell_id}"))

        age_plot_dir = out_root / "plots" / age
        age_plot_dir.mkdir(parents=True, exist_ok=True)
        plot_path = age_plot_dir / f"{session_id}__cell{cell_id}.png"

        npz_path = find_interim_npz(args.interim_root, session_id)
        wf = None
        if npz_path is not None:
            wf = extract_typical_waveform(npz_path, cell_id, args.waveform_fs_hz)

        fig, axes = plt.subplots(1, 2, figsize=(10, 4))

        # Left: waveform
        ax0 = axes[0]
        if wf is None:
            ax0.text(0.5, 0.5, "Waveform not found", ha="center", va="center")
            ax0.set_axis_off()
        else:
            t_ms, y = wf
            ax0.plot(t_ms, y, lw=1.5)
            ax0.set_xlabel("time (ms)")
            ax0.set_ylabel("amplitude (a.u.)")
            ax0.set_title("Typical waveform")

        # Right: FR vs duration context
        ax1 = axes[1]
        pop = pop_by_age.get(age)
        if pop is not None and {"spk_duration_ms", "fr_hz"}.issubset(pop.columns):
            ax1.scatter(
                pop["spk_duration_ms"],
                pop["fr_hz"],
                s=8,
                alpha=0.25,
                color="gray",
                zorder=1,
            )
        x_hi, y_hi = get_neuron_xy(pop, row)
        ax1.scatter(
            [x_hi],
            [y_hi],
            s=80,
            alpha=1.0,
            color="red",
            edgecolors="black",
            linewidths=0.7,
            marker="o",
            zorder=5,
            label="discrepant neuron",
        )
        ax1.set_xlabel("spk_duration_ms")
        ax1.set_ylabel("fr_hz")
        ax1.set_title("Feature context")
        ax1.legend(frameon=False, fontsize=8)

        fig.suptitle(
            f"{unit_uid} | type_u={row.get('type_u_type','?')} | pred={row.get('pred_type','?')}",
            fontsize=10,
        )
        fig.tight_layout()
        fig.savefig(plot_path, dpi=140)
        plt.close(fig)

        review_rows.append(
            {
                "age_group": age,
                "session_id": session_id,
                "cell_id": cell_id,
                "unit_uid": unit_uid,
                "type_u_type": row.get("type_u_type", pd.NA),
                "pred_type": row.get("pred_type", pd.NA),
                "allcel__type_u": row.get("allcel__type_u", pd.NA),
                "gmm_cluster": row.get("gmm_cluster", pd.NA),
                "gmm_pmax": row.get("gmm_pmax", pd.NA),
                "gmm_margin": row.get("gmm_margin", pd.NA),
                "fr_hz": row.get("fr_hz", pd.NA),
                "spk_duration_ms": row.get("spk_duration_ms", pd.NA),
                "spk_peaktrough_ms": row.get("spk_peaktrough_ms", pd.NA),
                "spk_asymmetry": row.get("spk_asymmetry", pd.NA),
                "cv2": row.get("cv2", pd.NA),
                "burst_index": row.get("burst_index", pd.NA),
                "interim_npz_found": bool(npz_path is not None),
                "interim_npz_path": str(npz_path) if npz_path is not None else "",
                "review_plot_png": str(plot_path),
            }
        )

    review_df = pd.DataFrame(review_rows)
    review_df.to_csv(out_root / "discrepant_neurons_review.csv", index=False)
    review_df.to_parquet(out_root / "discrepant_neurons_review.parquet", index=False)
    print(f"Wrote {len(review_df)} rows to {out_root / 'discrepant_neurons_review.csv'}")


if __name__ == "__main__":
    main()
