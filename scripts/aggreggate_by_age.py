from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional, List, Dict, Tuple

import numpy as np
import pandas as pd


DEFAULT_FEATURES = [
    "fr_hz",
    "burst_index",
    "cv2",
    "spk_duration_ms",
    "spk_peaktrough_ms",
    "spk_asymmetry",
    "refractory_ms_edge",
    "acg_peak_latency_ms",
    "type_u"
]


def normalize_session_name(value: object) -> str:
    """Normalize session identifiers to make schedule matching robust."""
    return " ".join(str(value).strip().split())


def age_group_from_age(age: int) -> Optional[str]:
    if age in (15, 16):
        return "P15_16"
    if age in (17, 18):
        return "P17_18"
    if age in (19, 20):
        return "P19_20"
    if age in (21, 22):
        return "P21_22"
    if age in (23, 24):
        return "P23_24"
    if age == 25:
        return "P25"
    return None


def load_session_metadata(excel_path: Path) -> pd.DataFrame:
    meta = pd.read_excel(excel_path, sheet_name='VINCA')
    # Expect columns like: SessionName, Age
    # Normalize
    if "SessionName" not in meta.columns or "Age" not in meta.columns:
        raise ValueError(f"Excel must contain columns 'SessionName' and 'Age'. Got: {list(meta.columns)}")

    meta = meta.copy()
    meta["SessionName"] = meta["SessionName"].astype(str)
    meta["session_key"] = meta["SessionName"].map(normalize_session_name)
    meta["Age"] = pd.to_numeric(meta["Age"], errors="coerce").astype("Int64")
    meta["age_group"] = meta["Age"].apply(lambda x: age_group_from_age(int(x)) if pd.notna(x) else None)
    return meta[["SessionName", "session_key", "Age", "age_group"]]


def iter_feature_parquets(processed_root: Path) -> List[Path]:
    # data/processed/<mouse>/features/*_features.parquet
    return sorted(processed_root.glob("*/features/*_features.parquet"))


def load_all_features(processed_root: Path) -> pd.DataFrame:
    files = iter_feature_parquets(processed_root)
    if not files:
        raise FileNotFoundError(f"No feature parquet files found under {processed_root}/<mouse>/features/")

    dfs = []
    for fp in files:
        df = pd.read_parquet(fp)
        df["features_file"] = fp.name
        df["mouse_dir"] = fp.parent.parent.name  # <mouse>
        dfs.append(df)

    out = pd.concat(dfs, ignore_index=True)

    # Ensure unit_uid exists
    if "unit_uid" not in out.columns:
        out["unit_uid"] = out["session_id"].astype(str) + "__cell" + out["cell_id"].astype(str)

    return out


def apply_qc_filters(df: pd.DataFrame, *, strict: bool = False) -> pd.DataFrame:
    qc_cols = [c for c in df.columns if c.startswith("qc_")]
    if not qc_cols:
        return df.copy()

    if strict:
        mask = df[qc_cols].fillna(False).all(axis=1)
    else:
        keep_cols = [c for c in ["qc_min_spikes", "qc_waveform", "qc_refractory"] if c in df.columns]
        mask = df[keep_cols].fillna(False).all(axis=1) if keep_cols else df[qc_cols].fillna(False).any(axis=1)

    return df.loc[mask].copy()


def build_ml_matrix(
    df: pd.DataFrame,
    feature_cols: List[str],
    *,
    log_fr: bool = True,
    standardize: bool = True,
) -> Tuple[np.ndarray, List[str], pd.DataFrame, Dict]:
    df2 = df.copy()

    used_cols: List[str] = []
    for col in feature_cols:
        if col == "fr_hz" and log_fr:
            df2["log10_fr_hz"] = np.log10(df2["fr_hz"].clip(lower=1e-6))
            used_cols.append("log10_fr_hz")
        else:
            if col in df2.columns:
                used_cols.append(col)

    if not used_cols:
        raise ValueError("No usable feature columns found after checking dataframe columns.")

    X = df2[used_cols].to_numpy(dtype=np.float64)
    finite = np.isfinite(X).all(axis=1)
    df3 = df2.loc[finite].copy()
    X = X[finite]

    meta = {
        "input_feature_cols": feature_cols,
        "used_feature_cols": used_cols,
        "log_fr": log_fr,
        "standardize": standardize,
        "n_units_before_finite_filter": int(len(df2)),
        "n_units_after_finite_filter": int(len(df3)),
    }

    if standardize:
        mu = X.mean(axis=0)
        sd = X.std(axis=0, ddof=0)
        sd[sd == 0] = 1.0
        X = (X - mu) / sd
        meta["standardize_mu"] = mu.tolist()
        meta["standardize_sd"] = sd.tolist()

    return X, used_cols, df3, meta


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--processed_root", type=str, default="data/processed")
    ap.add_argument("--excel", type=str, required=False,default = "data/schedule.xlsx",help="Excel file mapping SessionName -> Age")
    #ap.add_argument("--excel", type=str, required=True, help="Excel file mapping SessionName -> Age")
    ap.add_argument("--outdir", type=str, default="results")

    ap.add_argument("--qc", action="store_true", help="Write clean_units with QC filtering")
    ap.add_argument("--qc_strict", action="store_true")

    ap.add_argument("--features", type=str, default=",".join(DEFAULT_FEATURES))
    ap.add_argument("--no_log_fr", action="store_true")
    ap.add_argument("--no_standardize", action="store_true")

    ap.add_argument("--drop_unmapped", action="store_true",
                    help="Drop sessions not found in Excel (otherwise keep with age_group=NaN but they won't be written).")
    args = ap.parse_args()

    processed_root = Path(args.processed_root)
    outdir = Path(args.outdir)

    meta = load_session_metadata(Path(args.excel))
    df = load_all_features(processed_root)
    df["session_key"] = df["session_id"].map(normalize_session_name)

    # Merge session metadata
    df = df.merge(meta, on="session_key", how="left")
    if args.drop_unmapped:
        df = df[df["age_group"].notna()].copy()

    # Report mapping
    n_mapped = int(df["age_group"].notna().sum())
    print(f"Loaded {len(df)} units total; mapped to age groups: {n_mapped}")

    # Only write supported age groups
    groups = ["P15_16", "P17_18", "P19_20", "P21_22", "P23_24", "P25"]
    feature_cols = [s.strip() for s in args.features.split(",") if s.strip()]

    for g in groups:
        df_g = df[df["age_group"] == g].copy()
        if df_g.empty:
            print(f"{g}: no units, skipping")
            continue

        gdir = outdir / g
        gdir.mkdir(parents=True, exist_ok=True)

        all_path = gdir / f"{g}_all_units.parquet"
        df_g.to_parquet(all_path, index=False)

        df_clean = df_g
        if args.qc:
            df_clean = apply_qc_filters(df_clean, strict=args.qc_strict)

        clean_path = gdir / f"{g}_clean_units.parquet"
        df_clean.to_parquet(clean_path, index=False)

        # ML matrix from clean
        X, used_cols, df_ml, ml_meta = build_ml_matrix(
            df_clean,
            feature_cols=feature_cols,
            log_fr=not args.no_log_fr,
            standardize=not args.no_standardize,
        )

        ml_path = gdir / f"{g}_ml_units.parquet"
        df_ml.to_parquet(ml_path, index=False)

        npz_path = gdir / f"{g}_X.npz"
        np.savez_compressed(
            npz_path,
            X=X,
            feature_names=np.array(used_cols, dtype=object),
            unit_uid=df_ml["unit_uid"].to_numpy(dtype=object),
            session_id=df_ml["session_id"].to_numpy(dtype=object),
            cell_id=df_ml["cell_id"].to_numpy(dtype=np.int64),
            mouse=df_ml["mouse"].to_numpy(dtype=object) if "mouse" in df_ml.columns else df_ml["mouse_dir"].to_numpy(dtype=object),
            age=df_ml["Age"].to_numpy(),
        )

        ml_meta.update({
            "age_group": g,
            "n_all_units": int(len(df_g)),
            "n_clean_units": int(len(df_clean)),
            "n_ml_units": int(len(df_ml)),
        })

        meta_path = gdir / f"{g}_ml_meta.json"
        meta_path.write_text(json.dumps(ml_meta, indent=2))

        print(f"{g}: wrote {all_path.name}, {clean_path.name}, {ml_path.name}, {npz_path.name}")

    print("Done.")


if __name__ == "__main__":
    main()


