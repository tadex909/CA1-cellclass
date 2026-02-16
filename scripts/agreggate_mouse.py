from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd


DEFAULT_FEATURES = [
    "fr_hz",
    "burst_index",
    "cv2",
    "spk_duration_ms",
    "spk_peaktrough_ms",
    "spk_asymmetry",
    "refractory_ms_center",
    "acg_peak_latency_ms",
]


def load_mouse_features(mouse_dir: Path) -> pd.DataFrame:
    feat_dir = mouse_dir / "features"
    if not feat_dir.is_dir():
        raise FileNotFoundError(f"Missing features dir: {feat_dir}")

    files = sorted(feat_dir.glob("*_features.parquet"))
    if not files:
        raise FileNotFoundError(f"No *_features.parquet found in {feat_dir}")

    dfs = []
    for fp in files:
        df = pd.read_parquet(fp)
        df["features_file"] = fp.name
        dfs.append(df)

    out = pd.concat(dfs, ignore_index=True)
    out["unit_uid"] = out["session_id"].astype(str) + "__cell" + out["cell_id"].astype(str)
    return out


def apply_qc_filters(df: pd.DataFrame, *, strict: bool) -> pd.DataFrame:
    qc_cols = [c for c in df.columns if c.startswith("qc_")]
    if not qc_cols:
        return df.copy()

    if strict:
        mask = df[qc_cols].fillna(False).all(axis=1)
    else:
        # moderate defaults: require these if present; else fall back to any qc_*
        keep_cols = [c for c in ["qc_min_spikes", "qc_waveform", "qc_refractory"] if c in df.columns]
        mask = df[keep_cols].fillna(False).all(axis=1) if keep_cols else df[qc_cols].fillna(False).any(axis=1)

    return df.loc[mask].copy()


def build_ml_matrix(
    df: pd.DataFrame,
    feature_cols: list[str],
    *,
    log_fr: bool = True,
    standardize: bool = True,
) -> tuple[np.ndarray, list[str], pd.DataFrame, dict]:
    """
    Returns:
      X: (n_units, n_features)
      feature_names
      df_out: df filtered to finite rows, with transformed columns added
      meta: dict with transform info
    """
    df2 = df.copy()

    # Optional transform: log10(fr_hz)
    used = []
    for col in feature_cols:
        if col == "fr_hz" and log_fr:
            new = "log10_fr_hz"
            df2[new] = np.log10(df2["fr_hz"].clip(lower=1e-6))
            used.append(new)
        else:
            used.append(col)

    # Keep only columns that exist
    used = [c for c in used if c in df2.columns]
    if not used:
        raise ValueError("No usable feature columns found in dataframe.")

    X = df2[used].to_numpy(dtype=np.float64)
    finite_mask = np.isfinite(X).all(axis=1)
    df3 = df2.loc[finite_mask].copy()
    X = X[finite_mask]

    meta = {
        "input_feature_cols": feature_cols,
        "used_feature_cols": used,
        "log_fr": log_fr,
        "standardize": standardize,
        "n_units_before_finite_filter": int(len(df2)),
        "n_units_after_finite_filter": int(len(df3)),
    }

    if standardize:
        mu = X.mean(axis=0)
        sd = X.std(axis=0, ddof=0)
        sd[sd == 0] = 1.0
        Xz = (X - mu) / sd
        meta["standardize_mu"] = mu.tolist()
        meta["standardize_sd"] = sd.tolist()
        return Xz, used, df3, meta

    return X, used, df3, meta


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mouse", type=str, required=True, help="Mouse code, e.g. VS67")
    ap.add_argument("--processed_root", type=str, default="data/processed")
    ap.add_argument("--outdir", type=str, default="results", help="Where to write aggregated outputs")

    ap.add_argument("--qc", action="store_true", help="Apply QC filtering for CLEAN output")
    ap.add_argument("--qc_strict", action="store_true", help="Require all qc_* True (if --qc)")

    ap.add_argument("--features", type=str, default=",".join(DEFAULT_FEATURES),
                    help="Comma-separated feature columns to include in ML matrix")
    ap.add_argument("--no_log_fr", action="store_true", help="Do not log-transform fr_hz")
    ap.add_argument("--no_standardize", action="store_true", help="Do not z-score features")

    args = ap.parse_args()

    processed_root = Path(args.processed_root)
    mouse_dir = processed_root / args.mouse
    outdir = Path(args.outdir) / args.mouse
    outdir.mkdir(parents=True, exist_ok=True)

    df_all = load_mouse_features(mouse_dir)
    all_path = outdir / f"{args.mouse}_all_units.parquet"
    df_all.to_parquet(all_path, index=False)

    # CLEAN
    df_clean = df_all
    if args.qc:
        df_clean = apply_qc_filters(df_clean, strict=args.qc_strict)

    clean_path = outdir / f"{args.mouse}_clean_units.parquet"
    df_clean.to_parquet(clean_path, index=False)

    # ML matrix
    feature_cols = [s.strip() for s in args.features.split(",") if s.strip()]
    X, used_cols, df_ml, meta = build_ml_matrix(
        df_clean,
        feature_cols,
        log_fr=not args.no_log_fr,
        standardize=not args.no_standardize,
    )

    ml_parquet = outdir / f"{args.mouse}_ml_units.parquet"
    df_ml.to_parquet(ml_parquet, index=False)

    ml_npz = outdir / f"{args.mouse}_X.npz"
    np.savez_compressed(
        ml_npz,
        X=X,
        feature_names=np.array(used_cols, dtype=object),
        unit_uid=df_ml["unit_uid"].to_numpy(dtype=object),
        session_id=df_ml["session_id"].to_numpy(dtype=object),
        cell_id=df_ml["cell_id"].to_numpy(dtype=np.int64),
    )

    meta_path = outdir / f"{args.mouse}_ml_meta.json"
    meta_path.write_text(json.dumps(meta, indent=2))

    print("Wrote:")
    print(" ", all_path, f"(rows={len(df_all)})")
    print(" ", clean_path, f"(rows={len(df_clean)})")
    print(" ", ml_parquet, f"(rows={len(df_ml)})")
    print(" ", ml_npz)
    print(" ", meta_path)


if __name__ == "__main__":
    main()
