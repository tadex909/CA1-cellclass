#!/usr/bin/env python3
"""
Inspect a converted NPZ file produced by mat_to_npz and verify key arrays exist.

Usage:
  python scripts/inspect_npz.py --file data/interim/VS57/2022-12-18/VS57_2022-12-18_18-51-04_allcel.npz
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Tuple

import numpy as np


def summarize(name: str, arr: np.ndarray) -> str:
    s = f"{name}: shape={arr.shape}, dtype={arr.dtype}"
    if arr.dtype.kind in ("i", "u", "f") and arr.size:
        try:
            mn = float(np.nanmin(arr.astype(float)))
            mx = float(np.nanmax(arr.astype(float)))
            s += f", min={mn:g}, max={mx:g}"
        except Exception:
            pass
    return s


def get_required_keys() -> Tuple[str, ...]:
    return (
        "allcel__id_cel",
        "allcel__itime_spk",
        "allcel__id_spk",
        "allcel__bestwaveform",
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    #ap.add_argument("--file", required=True, type=str, help="Path to one .npz file")
    ap.add_argument("--file", required=False, default="data/interim/VS57/2022-12-18/VS57_2022-12-18_18-51-04_allcel.npz", type=str, help="Path to one .npz file")
    ap.add_argument("--list-keys", action="store_true", help="List all keys in the npz")
    args = ap.parse_args()

    npz_path = Path(args.file)
    if not npz_path.exists():
        raise SystemExit(f"File not found: {npz_path}")

    data = np.load(npz_path, allow_pickle=False)

    print(f"Loaded: {npz_path}")
    print(f"Keys: {len(data.files)}")

    if args.list_keys:
        for k in sorted(data.files):
            print("  -", k)

    # Meta
    if "meta_json" in data.files:
        raw = data["meta_json"]
        # meta_json stored as np.string_ scalar/array
        meta_str = raw.tobytes().decode("utf-8") if raw.dtype.kind in ("S", "a") else str(raw)
        try:
            meta = json.loads(meta_str)
        except Exception:
            meta = {"_raw": meta_str}
        print("\nmeta_json:")
        for k in sorted(meta.keys()):
            print(f"  {k}: {meta[k]}")
    else:
        print("\nWARNING: missing meta_json")

    # Required arrays
    print("\nRequired fields check:")
    ok = True
    arrays: Dict[str, np.ndarray] = {}
    for k in get_required_keys():
        if k not in data.files:
            print(f"  ✗ missing: {k}")
            ok = False
        else:
            arr = data[k]
            arrays[k] = arr
            print(f"  ✓ {summarize(k, arr)}")

    # Optional but useful waveform arrays
    for k in ("allcel__waveform", "allcel__meanwaveform", "allcel__bestswaveforms"):
        if k in data.files:
            print(f"  ✓ {summarize(k, data[k])}")

    # Consistency checks (only if present)
    print("\nConsistency checks:")
    if "allcel__id_spk" in arrays and "allcel__itime_spk" in arrays:
        n1 = arrays["allcel__id_spk"].shape[0]
        n2 = arrays["allcel__itime_spk"].shape[0]
        if n1 != n2:
            print(f"  ✗ spike length mismatch: id_spk={n1} vs itime_spk={n2}")
            ok = False
        else:
            print(f"  ✓ spikes length consistent: {n1}")

    if "allcel__id_cel" in arrays and "allcel__bestwaveform" in arrays:
        n_cells = arrays["allcel__id_cel"].size
        bw = arrays["allcel__bestwaveform"]
        # Expect bestwaveform ~ (51, n_cells) or (n_cells, 51) depending on how it was saved
        if bw.ndim != 2:
            print(f"  ✗ bestwaveform unexpected ndim={bw.ndim} (expected 2)")
            ok = False
        else:
            if n_cells in bw.shape:
                print(f"  ✓ bestwaveform matches n_cells={n_cells} (shape={bw.shape})")
            else:
                print(f"  ✗ bestwaveform shape {bw.shape} does not include n_cells={n_cells}")
                ok = False

    # Quick content sanity
    if "allcel__itime_spk" in arrays:
        it = arrays["allcel__itime_spk"]
        if it.dtype.kind not in ("i", "u"):
            print(f"  ! itime_spk dtype is {it.dtype} (often uint64 is expected)")
        if it.size and np.any(np.diff(it[: min(it.size, 1000)]) < 0):
            print("  ! itime_spk is not non-decreasing (first 1000 samples) — check conversion/order")
        else:
            print("  ✓ itime_spk looks non-decreasing (first chunk)")

    print("\nRESULT:", "OK" if ok else "FAILED")
    if not ok:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
