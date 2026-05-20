#!/usr/bin/env python3
"""
Export all top-level MATLAB structs from one *_Ratemap*.mat file to separate NPZ files.

Example:
  python src/cellclass/mat_to_npz_structs.py \
    --input data/raw/VS57_2022-12-18_18-51-04_Ratemap.mat \
    --output data/interim/VS57_2022-12-18_18-51-04_structs
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

import numpy as np


def safe_name(name: str) -> str:
    out = name.strip()
    out = re.sub(r"\s+", "_", out)
    out = re.sub(r"[^A-Za-z0-9_\-]+", "-", out)
    out = re.sub(r"-{2,}", "-", out).strip("-")
    return out or "unnamed"


def load_mat_file(mat_path: Path) -> dict[str, Any]:
    """
    Load .mat content into native Python structures.
    Supports:
      - MATLAB v7.2 and earlier via scipy.io.loadmat
      - MATLAB v7.3 via mat73 (if installed)
    """
    try:
        from scipy.io import loadmat
        from scipy.io.matlab import mat_struct  # type: ignore

        raw = loadmat(mat_path, squeeze_me=True, struct_as_record=False)
        raw = {k: v for k, v in raw.items() if not k.startswith("__")}

        def _convert(x: Any) -> Any:
            if isinstance(x, mat_struct):
                return {fn: _convert(getattr(x, fn)) for fn in x._fieldnames}
            if isinstance(x, np.ndarray) and x.dtype == object:
                return np.array([_convert(v) for v in x.ravel()], dtype=object).reshape(x.shape)
            return x

        return {k: _convert(v) for k, v in raw.items()}

    except NotImplementedError:
        pass
    except Exception as e:
        last_err = e
    else:
        last_err = None  # type: ignore

    try:
        import mat73  # type: ignore

        out = mat73.loadmat(str(mat_path))
        if not isinstance(out, dict):
            raise TypeError("mat73.loadmat did not return a dictionary.")
        return out
    except Exception as e:
        if last_err is not None:
            raise RuntimeError(
                f"Failed to read {mat_path.name}. If this is a MATLAB v7.3 file, "
                f"install 'mat73' (pip install mat73). Original error: {last_err}"
            ) from e
        raise RuntimeError(
            f"Failed to read {mat_path.name}. If this is a MATLAB v7.3 file, "
            f"install 'mat73' (pip install mat73)."
        ) from e


def to_jsonable(x: Any) -> Any:
    if isinstance(x, dict):
        return {k: to_jsonable(v) for k, v in x.items()}
    if isinstance(x, np.ndarray):
        if x.dtype == object:
            return [to_jsonable(v) for v in x.ravel()]
        return x.tolist()
    if isinstance(x, (list, tuple)):
        return [to_jsonable(v) for v in x]
    if isinstance(x, (np.integer, np.floating, np.bool_)):
        return x.item()
    return x


def flatten_vector(x: np.ndarray) -> np.ndarray:
    if x.ndim == 2 and 1 in x.shape:
        return x.reshape(-1)
    return x


def maybe_cast_int(arr: np.ndarray) -> np.ndarray:
    if arr.dtype.kind in ("i", "u"):
        return arr
    if arr.dtype.kind == "f":
        if np.all(np.isfinite(arr)) and np.allclose(arr, np.round(arr)):
            return np.round(arr).astype(np.int64)
    return arr


def is_struct_like(x: Any) -> bool:
    if isinstance(x, dict):
        return True
    if isinstance(x, (list, tuple)):
        return any(isinstance(v, dict) for v in x)
    if isinstance(x, np.ndarray) and x.dtype == object:
        return any(isinstance(v, dict) for v in x.ravel())
    return False


def struct_to_records(obj: Any, struct_name: str) -> list[dict[str, Any]]:
    if isinstance(obj, dict):
        return [obj]
    if isinstance(obj, (list, tuple)):
        recs = [x for x in obj if isinstance(x, dict)]
        if recs:
            return recs
    if isinstance(obj, np.ndarray) and obj.dtype == object:
        recs = [x for x in obj.ravel() if isinstance(x, dict)]
        if recs:
            return recs
    raise TypeError(
        f"Unexpected format for '{struct_name}'. Expected dict or array/list of dict records."
    )


def pack_single_value(value: Any) -> np.ndarray | None:
    if isinstance(value, np.ndarray):
        arr = flatten_vector(value)
    elif isinstance(value, (str, bytes, np.str_)):
        arr = np.asarray(value)
    elif np.isscalar(value):
        arr = np.asarray(value)
    elif isinstance(value, (list, tuple)):
        arr = np.asarray(value)
    else:
        return None

    if arr.dtype == object:
        return None
    return maybe_cast_int(arr)


def pack_multi_values(values: list[Any]) -> np.ndarray | None:
    arrays: list[np.ndarray] = []
    for value in values:
        packed = pack_single_value(value)
        if packed is None:
            return None
        arrays.append(packed)

    if not arrays:
        return None

    if len(arrays) == 1:
        return arrays[0]

    if all(a.ndim == 0 for a in arrays):
        raw = [a.item() for a in arrays]
        if all(isinstance(x, str) for x in raw):
            return np.asarray(raw, dtype=str)
        out = np.asarray(raw)
        if out.dtype == object:
            return None
        return maybe_cast_int(out)

    shape0 = arrays[0].shape
    if all(a.shape == shape0 for a in arrays):
        try:
            out = np.stack(arrays, axis=0)
        except Exception:
            return None
        if out.dtype != object:
            return maybe_cast_int(out)

    return None


def export_struct(
    mat_path: Path,
    struct_name: str,
    struct_obj: Any,
    out_dir: Path,
    overwrite: bool,
) -> Path:
    records = struct_to_records(struct_obj, struct_name=struct_name)
    fields = sorted({k for rec in records for k in rec.keys()})

    payload: dict[str, Any] = {}
    meta = {
        "source_path": str(mat_path.resolve()),
        "source_name": mat_path.name,
        "struct_name": struct_name,
        "n_records": len(records),
        "fields": fields,
    }
    payload["meta_json"] = np.array(json.dumps(meta), dtype=np.string_)

    for field in fields:
        vals = [rec.get(field, np.nan) for rec in records]
        packed = pack_multi_values(vals)
        key = f"{struct_name}__{field}"
        if packed is not None:
            payload[key] = packed
        else:
            payload[f"{key}__json"] = np.array(
                json.dumps([to_jsonable(v) for v in vals], default=str),
                dtype=np.string_,
            )

    out_name = f"{mat_path.stem}_{safe_name(struct_name)}.npz"
    out_path = out_dir / out_name
    if out_path.exists() and not overwrite:
        raise FileExistsError(f"File exists (use --overwrite): {out_path}")
    np.savez_compressed(out_path, **payload)
    return out_path


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Export all top-level structs from one *_Ratemap*.mat file into separate .npz files."
        )
    )
    ap.add_argument("--input", required=True, type=str, help="Path to one *_Ratemap*.mat file.")
    ap.add_argument(
        "--output",
        required=True,
        type=str,
        help="Output directory where .npz files will be written.",
    )
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing output files.")
    args = ap.parse_args()

    mat_path = Path(args.input)
    if not mat_path.exists() or not mat_path.is_file():
        raise SystemExit(f"Input file not found: {mat_path}")
    if mat_path.suffix.lower() != ".mat":
        raise SystemExit(f"Input must be a .mat file: {mat_path}")
    if "ratemap" not in mat_path.name.lower():
        raise SystemExit(f"Expected '*_Ratemap*.mat' file, got: {mat_path.name}")

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    mat_dict = load_mat_file(mat_path)
    struct_items = [(k, v) for k, v in mat_dict.items() if is_struct_like(v)]
    skipped = [k for k, v in mat_dict.items() if not is_struct_like(v)]

    if not struct_items:
        keys = ", ".join(sorted(mat_dict.keys()))
        raise SystemExit(
            f"No top-level struct-like variables found in {mat_path.name}. Keys: {keys}"
        )

    print(f"Found {len(struct_items)} top-level struct(s) in {mat_path.name}.")
    exported: list[Path] = []
    for struct_name, struct_obj in struct_items:
        out_path = export_struct(
            mat_path=mat_path,
            struct_name=struct_name,
            struct_obj=struct_obj,
            out_dir=out_dir,
            overwrite=args.overwrite,
        )
        exported.append(out_path)
        print(f"[OK]  {struct_name} -> {out_path.name}")

    if skipped:
        print(f"[INFO] Skipped non-struct keys: {', '.join(sorted(skipped))}")
    print(f"Done. Wrote {len(exported)} file(s) to: {out_dir.resolve()}")


if __name__ == "__main__":
    main()
