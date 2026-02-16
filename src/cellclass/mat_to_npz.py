#!/usr/bin/env python3
"""
Convert MATLAB *_Ratemap*.mat files to compressed NPZ files containing the `allcel` struct.

Typical input filename:
  VS57_2022-12-18_18-51-04_Ratemap_final_thèse.mat

Output:
  data/interim/VS57/2022-12-18/VS57_2022-12-18_18-51-04_allcel.npz

The NPZ contains:
  - meta (JSON string): mouse, date, time, session_name, source_path, source_mtime, etc.
  - allcel__<field> arrays for each numeric field in allcel (Spikes excluded by default)
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Tuple

import numpy as np


# -------------------------
# Filename parsing
# -------------------------

@dataclass(frozen=True)
class SessionInfo:
    mouse: str
    date: str
    time: str
    stem: str  # original filename stem (no extension)

    @property
    def session_name(self) -> str:
        # Keep what MATLAB uses a lot: MOUSE_DATE_TIME
        if self.time:
            return f"{self.mouse}_{self.date}_{self.time}"
        return f"{self.mouse}_{self.date}"


FILENAME_RE = re.compile(
    r"^(?P<mouse>[^_]+)_(?P<date>\d{4}-\d{2}-\d{2})_(?P<time>\d{2}-\d{2}-\d{2})_.*$"
)


def parse_session_info(mat_path: Path) -> SessionInfo:
    stem = mat_path.stem
    m = FILENAME_RE.match(stem)
    if not m:
        # fallback: take first token as mouse, and try to find date/time anywhere
        parts = stem.split("_")
        mouse = parts[0] if parts else "UNKNOWN"
        date = next((p for p in parts if re.fullmatch(r"\d{4}-\d{2}-\d{2}", p)), "UNKNOWNDATE")
        time = next((p for p in parts if re.fullmatch(r"\d{2}-\d{2}-\d{2}", p)), "")
        return SessionInfo(mouse=mouse, date=date, time=time, stem=stem)
    return SessionInfo(mouse=m.group("mouse"), date=m.group("date"), time=m.group("time"), stem=stem)


def safe_filename(s: str) -> str:
    # Keep letters/numbers/_/-
    s = s.strip()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^A-Za-z0-9_\-]+", "-", s)
    s = re.sub(r"-{2,}", "-", s).strip("-")
    return s


# -------------------------
# MATLAB loading helpers
# -------------------------

def load_mat_file(mat_path: Path) -> Dict[str, Any]:
    """
    Loads MATLAB .mat into a Python dict.
    Supports v7.2 (scipy.io.loadmat) and optionally v7.3 (mat73).
    """
    try:
        from scipy.io import loadmat
        from scipy.io.matlab import mat_struct  # type: ignore

        md = loadmat(mat_path, squeeze_me=True, struct_as_record=False)
        # Remove scipy metadata keys
        md = {k: v for k, v in md.items() if not k.startswith("__")}

        # Convert mat_struct recursively
        def _convert(x: Any) -> Any:
            if isinstance(x, mat_struct):
                return {fn: _convert(getattr(x, fn)) for fn in x._fieldnames}
            if isinstance(x, np.ndarray) and x.dtype == object:
                # convert each element
                return np.array([_convert(e) for e in x.ravel()], dtype=object).reshape(x.shape)
            return x

        return {k: _convert(v) for k, v in md.items()}

    except NotImplementedError:
        # Likely v7.3
        pass
    except Exception as e:
        # Could still be v7.3 or corrupted
        # We'll try v7.3 loader next, otherwise re-raise
        last_err = e
    else:
        last_err = None  # type: ignore

    # Try mat73 (optional dependency)
    try:
        import mat73  # type: ignore
        md = mat73.loadmat(str(mat_path))
        return md
    except Exception:
        if last_err is not None:
            raise RuntimeError(
                f"Failed to read {mat_path.name}. If this is a MATLAB v7.3 file, "
                f"install 'mat73' (pip install mat73). Original error: {last_err}"
            )
        raise RuntimeError(
            f"Failed to read {mat_path.name}. If this is a MATLAB v7.3 file, "
            f"install 'mat73' (pip install mat73)."
        )


def extract_allcel(mat_dict: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extract allcel/allcell struct as a plain dict.
    """
    if "allcel" in mat_dict and isinstance(mat_dict["allcel"], dict):
        return mat_dict["allcel"]
    if "allcell" in mat_dict and isinstance(mat_dict["allcell"], dict):
        return mat_dict["allcell"]

    # Some files might store it nested (rare); provide helpful message
    keys = ", ".join(sorted(mat_dict.keys()))
    raise KeyError(f"Could not find 'allcel' (or 'allcell') in MAT file. Top-level keys: {keys}")


# Better: just treat any field starting with id_ or containing "channel" as int-ish if it looks integer
def maybe_cast_int(name: str, arr: Any) -> Any:
    if not isinstance(arr, np.ndarray):
        return arr
    if arr.dtype.kind in ("i", "u"):
        return arr
    # If floats but all close to integers -> cast
    if arr.dtype.kind == "f":
        if np.all(np.isfinite(arr)) and np.allclose(arr, np.round(arr)):
            return np.round(arr).astype(np.int64)
    return arr


def flatten_cell(x: Any) -> Any:
    # Convert MATLAB column vectors like (N,1) into (N,)
    if isinstance(x, np.ndarray) and x.ndim == 2 and 1 in x.shape:
        return x.reshape(-1)
    return x

# -------------------------
# Conversion
# -------------------------

def convert_one(mat_path: Path, out_root: Path, overwrite: bool, include_spikes: bool) -> Path:
    if "Ratemap" not in mat_path.name:
        raise ValueError(f"File does not contain 'Ratemap' in name: {mat_path.name}")

    info = parse_session_info(mat_path)

    # Output path: data/interim/MOUSE/DATE/<session>_allcel.npz
    out_dir = out_root / safe_filename(info.mouse) / safe_filename(info.date)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{safe_filename(info.session_name)}_allcel.npz"

    if out_path.exists() and not overwrite:
        return out_path

    md = load_mat_file(mat_path)
    allcel = extract_allcel(md)

    # Build NPZ payload
    payload: Dict[str, Any] = {}

    meta = {
        "mouse": info.mouse,
        "date": info.date,
        "time": info.time,
        "session_name": info.session_name,
        "source_path": str(mat_path.resolve()),
        "source_mtime": mat_path.stat().st_mtime,
        "mat_keys": sorted(md.keys()),
        "allcel_fields": sorted(allcel.keys()) if isinstance(allcel, dict) else None,
        "include_spikes": include_spikes,
    }
    payload["meta_json"] = np.array(json.dumps(meta), dtype=np.string_)

    for k, v in allcel.items():
        if (not include_spikes) and (k.lower() == "spikes"):
            continue

        # Many entries are numeric arrays already; some may be nested dicts (structs)
        if isinstance(v, dict):
            # Store nested struct as JSON (small) or skip if large
            payload[f"allcel__{k}__json"] = np.array(json.dumps(v, default=str), dtype=np.string_)
            continue

        if isinstance(v, np.ndarray):
            v = flatten_cell(v)
            v = maybe_cast_int(k, v)
            payload[f"allcel__{k}"] = v
        else:
            # store scalars / strings as JSON-safe
            payload[f"allcel__{k}__json"] = np.array(json.dumps(v, default=str), dtype=np.string_)

    np.savez_compressed(out_path, **payload)
    return out_path


def iter_mat_files(input_root: Path, recursive: bool) -> Iterable[Path]:
    if input_root.is_file():
        yield input_root
        return
    if recursive:
        yield from input_root.rglob("*.mat")
    else:
        yield from input_root.glob("*.mat")


def main() -> None:
    
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=str, default="Epsztein-nas02/TEAM/Tadeo/data/raw", help="Input .mat file or directory (can be network path).")
    ap.add_argument("--output", type=str, default="data/interim", help="Output directory for .npz files.")
    ap.add_argument("--recursive", action="store_true", help="Recursively search for .mat files.")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing .npz files.")
    ap.add_argument("--include-spikes", action="store_true", help="Include allcel.Spikes if present.")
    ap.add_argument("--dry-run", action="store_true", help="Only print what would be converted.")
    args = ap.parse_args()

    raw_dir = Path(r"\\Epsztein-nas02\TEAM\Tadeo\data\raw")
    print(raw_dir.exists())

    in_path = raw_dir
    out_root = Path(args.output)
    out_root.mkdir(parents=True, exist_ok=True)

    mat_files = list(iter_mat_files(in_path, recursive=args.recursive))
    # Filter to those containing "Ratemap" in the name
    mat_files = [p for p in mat_files if ("Ratemap" in p.name and p.suffix.lower() == ".mat")]

    if not mat_files:
        raise SystemExit(f"No '*Ratemap*.mat' files found in {in_path}")

    print(f"Found {len(mat_files)} file(s).")
    for p in mat_files:
        info = parse_session_info(p)
        out_dir = out_root / safe_filename(info.mouse) / safe_filename(info.date)
        out_path = out_dir / f"{safe_filename(info.session_name)}_allcel.npz"
        if args.dry_run:
            print(f"[DRY] {p} -> {out_path}")
            continue

        try:
            outp = convert_one(p, out_root=out_root, overwrite=args.overwrite, include_spikes=args.include_spikes)
            print(f"[OK]  {p.name} -> {outp}")
        except Exception as e:
            print(f"[ERR] {p.name}: {e}")


if __name__ == "__main__":
    main()