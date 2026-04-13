#!/usr/bin/env python3
"""
Convert MATLAB files to compressed NPZ files.

Supported modes:
1) ratemap  -> expects *_Ratemap*.mat and exports `allcel` (+ selected `allpf`)
               and, when present, also exports `pf` to a sibling *_pf.npz.
2) trajdata -> expects *_TrajData*.mat and exports only `Traj` fields.
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable

import numpy as np

DEFAULT_TRAJ_FIELDS = [
    "Cond",
    "time",
    "Wheel",
    "VRtraj",
    "condition",
    "Speed",
    "XSpeed",
    "binSpX",
    "BinSpW",
    "WB",
    "start",
    "stop",
    "tstart",
    "tstop",
    "endVR",
]


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
        parts = stem.split("_")
        mouse = parts[0] if parts else "UNKNOWN"
        date = next((p for p in parts if re.fullmatch(r"\d{4}-\d{2}-\d{2}", p)), "UNKNOWNDATE")
        time = next((p for p in parts if re.fullmatch(r"\d{2}-\d{2}-\d{2}", p)), "")
        return SessionInfo(mouse=mouse, date=date, time=time, stem=stem)
    return SessionInfo(mouse=m.group("mouse"), date=m.group("date"), time=m.group("time"), stem=stem)


def safe_filename(s: str) -> str:
    s = s.strip()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^A-Za-z0-9_\-]+", "-", s)
    s = re.sub(r"-{2,}", "-", s).strip("-")
    return s


def parse_csv_list(raw: str | None) -> list[str]:
    if not raw:
        return []
    return [x.strip() for x in raw.split(",") if x.strip()]


# -------------------------
# MATLAB loading helpers
# -------------------------

def load_mat_file(mat_path: Path, variable_names: list[str] | None = None) -> Dict[str, Any]:
    """
    Loads MATLAB .mat into a Python dict.
    Supports v7.2 (scipy.io.loadmat) and optionally v7.3 (mat73).
    """
    try:
        from scipy.io import loadmat
        from scipy.io.matlab import mat_struct  # type: ignore

        load_kwargs: Dict[str, Any] = {
            "squeeze_me": True,
            "struct_as_record": False,
        }
        if variable_names:
            load_kwargs["variable_names"] = variable_names
        md = loadmat(mat_path, **load_kwargs)
        md = {k: v for k, v in md.items() if not k.startswith("__")}

        def _convert(x: Any) -> Any:
            if isinstance(x, mat_struct):
                return {fn: _convert(getattr(x, fn)) for fn in x._fieldnames}
            if isinstance(x, np.ndarray) and x.dtype == object:
                return np.array([_convert(e) for e in x.ravel()], dtype=object).reshape(x.shape)
            return x

        return {k: _convert(v) for k, v in md.items()}

    except NotImplementedError:
        pass
    except Exception as e:
        last_err = e
    else:
        last_err = None  # type: ignore

    try:
        import mat73  # type: ignore

        # mat73 may not support selective variable loading across versions.
        md = mat73.loadmat(str(mat_path))
        if variable_names:
            keep = set(variable_names)
            md = {k: v for k, v in md.items() if k in keep}
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
    if "allcel" in mat_dict and isinstance(mat_dict["allcel"], dict):
        return mat_dict["allcel"]
    if "allcell" in mat_dict and isinstance(mat_dict["allcell"], dict):
        return mat_dict["allcell"]
    keys = ", ".join(sorted(mat_dict.keys()))
    raise KeyError(f"Could not find 'allcel' (or 'allcell') in MAT file. Top-level keys: {keys}")


def extract_optional_struct(mat_dict: Dict[str, Any], names: tuple[str, ...]) -> Dict[str, Any] | None:
    for n in names:
        if n in mat_dict and isinstance(mat_dict[n], dict):
            return mat_dict[n]
    return None


def extract_optional_value(mat_dict: Dict[str, Any], names: tuple[str, ...]) -> Any | None:
    for n in names:
        if n in mat_dict:
            return mat_dict[n]
    return None


def extract_traj(mat_dict: Dict[str, Any]) -> Any:
    if "Traj" in mat_dict:
        return mat_dict["Traj"]
    keys = ", ".join(sorted(mat_dict.keys()))
    raise KeyError(f"Could not find 'Traj' in MAT file. Top-level keys: {keys}")


def maybe_cast_int(_name: str, arr: Any) -> Any:
    if not isinstance(arr, np.ndarray):
        return arr
    if arr.dtype.kind in ("i", "u"):
        return arr
    if arr.dtype.kind == "f":
        if np.all(np.isfinite(arr)) and np.allclose(arr, np.round(arr)):
            return np.round(arr).astype(np.int64)
    return arr


def flatten_cell(x: Any) -> Any:
    if isinstance(x, np.ndarray) and x.ndim == 2 and 1 in x.shape:
        return x.reshape(-1)
    return x


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


def traj_to_records(traj_obj: Any) -> list[dict[str, Any]]:
    if isinstance(traj_obj, dict):
        return [traj_obj]
    if isinstance(traj_obj, (list, tuple)):
        out = [x for x in traj_obj if isinstance(x, dict)]
        if out:
            return out
    if isinstance(traj_obj, np.ndarray) and traj_obj.dtype == object:
        out = [x for x in traj_obj.ravel() if isinstance(x, dict)]
        if out:
            return out
    raise TypeError("Unexpected Traj format. Expected dict or array/list of dict records.")


def struct_to_records(struct_obj: Any, struct_name: str = "struct") -> list[dict[str, Any]]:
    if isinstance(struct_obj, dict):
        return [struct_obj]
    if isinstance(struct_obj, (list, tuple)):
        out = [x for x in struct_obj if isinstance(x, dict)]
        if out:
            return out
    if isinstance(struct_obj, np.ndarray) and struct_obj.dtype == object:
        out = [x for x in struct_obj.ravel() if isinstance(x, dict)]
        if out:
            return out
    raise TypeError(
        f"Unexpected {struct_name} format. Expected dict or array/list of dict records."
    )


def pack_traj_field(field_name: str, values: list[Any]) -> tuple[np.ndarray | None, bool]:
    arrs: list[np.ndarray] = []
    for v in values:
        if isinstance(v, np.ndarray):
            a = flatten_cell(v)
        else:
            a = np.asarray(v)
        arrs.append(a)

    if all(a.ndim == 0 for a in arrs):
        raw = [a.item() for a in arrs]
        if all(isinstance(x, str) for x in raw):
            return np.asarray(raw, dtype=str), False
        out = np.asarray(raw)
        if out.dtype != object:
            return maybe_cast_int(field_name, out), False
        return None, True

    shape0 = arrs[0].shape
    if all(a.shape == shape0 for a in arrs):
        try:
            out = np.stack(arrs, axis=0)
        except Exception:
            return None, True
        if out.dtype != object:
            return maybe_cast_int(field_name, out), False

    return None, True


# -------------------------
# Conversion
# -------------------------

def build_out_path(info: SessionInfo, out_root: Path, mode: str) -> Path:
    out_dir = out_root / safe_filename(info.mouse) / safe_filename(info.date)
    out_dir.mkdir(parents=True, exist_ok=True)
    if mode == "ratemap":
        return out_dir / f"{safe_filename(info.session_name)}_allcel.npz"
    return out_dir / f"{safe_filename(info.session_name)}_trajdata.npz"


def build_pf_out_path(info: SessionInfo, out_root: Path, pf_tag: str | None = None) -> Path:
    out_dir = out_root / safe_filename(info.mouse) / safe_filename(info.date)
    out_dir.mkdir(parents=True, exist_ok=True)
    base = f"{safe_filename(info.session_name)}_pf"
    if pf_tag:
        tag = safe_filename(pf_tag).strip("_-")
        if tag:
            base = f"{base}_{tag}"
    return out_dir / f"{base}.npz"


def convert_one(
    mat_path: Path,
    out_root: Path,
    overwrite: bool,
    include_spikes: bool,
    mode: str,
    traj_fields_req: list[str] | None = None,
    pf_tag: str | None = None,
) -> Path:
    marker = "Ratemap" if mode == "ratemap" else "TrajData"
    if marker.lower() not in mat_path.name.lower():
        raise ValueError(f"File does not contain '{marker}' in name: {mat_path.name}")

    info = parse_session_info(mat_path)
    out_path = build_out_path(info, out_root, mode=mode)
    skip_main_write = out_path.exists() and not overwrite
    if skip_main_write:
        if mode != "ratemap":
            return out_path
        pf_out_path = build_pf_out_path(info, out_root=out_root, pf_tag=pf_tag)
        if pf_out_path.exists():
            return out_path

    if mode == "ratemap":
        # Load only what we need from potentially large/complex MAT files.
        md = load_mat_file(mat_path, variable_names=["allcel", "allcell", "allpf", "pf"])
    else:
        md = load_mat_file(mat_path, variable_names=["Traj"])
    payload: Dict[str, Any] = {}

    meta: Dict[str, Any] = {
        "mouse": info.mouse,
        "date": info.date,
        "time": info.time,
        "session_name": info.session_name,
        "mode": mode,
        "source_path": str(mat_path.resolve()),
        "source_mtime": mat_path.stat().st_mtime,
        "mat_keys": sorted(md.keys()),
    }

    if mode == "ratemap":
        allcel = extract_allcel(md)
        allpf = extract_optional_struct(md, ("allpf",))
        pf_raw = extract_optional_value(md, ("pf",))
        meta["allcel_fields"] = sorted(allcel.keys())
        meta["allpf_fields"] = sorted(allpf.keys()) if isinstance(allpf, dict) else None
        meta["pf_present_in_mat"] = pf_raw is not None
        meta["include_spikes"] = include_spikes

        for k, v in allcel.items():
            if (not include_spikes) and (k.lower() == "spikes"):
                continue
            if isinstance(v, dict):
                payload[f"allcel__{k}__json"] = np.array(json.dumps(v, default=str), dtype=np.string_)
                continue
            if isinstance(v, np.ndarray):
                vv = maybe_cast_int(k, flatten_cell(v))
                payload[f"allcel__{k}"] = vv
            else:
                payload[f"allcel__{k}__json"] = np.array(json.dumps(v, default=str), dtype=np.string_)

        if isinstance(allpf, dict) and "ispf_cxu" in allpf:
            v = allpf["ispf_cxu"]
            if isinstance(v, np.ndarray):
                payload["allpf__ispf_cxu"] = maybe_cast_int("ispf_cxu", v)
            else:
                payload["allpf__ispf_cxu__json"] = np.array(json.dumps(v, default=str), dtype=np.string_)

        pf_written = False
        pf_out_path = build_pf_out_path(info, out_root=out_root, pf_tag=pf_tag)
        if pf_raw is not None:
            try:
                pf_records = struct_to_records(pf_raw, struct_name="pf")
                pf_fields = sorted({k for rec in pf_records for k in rec.keys()})
                pf_meta: Dict[str, Any] = {
                    "mouse": info.mouse,
                    "date": info.date,
                    "time": info.time,
                    "session_name": info.session_name,
                    "mode": "pf",
                    "source_path": str(mat_path.resolve()),
                    "source_mtime": mat_path.stat().st_mtime,
                    "mat_keys": sorted(md.keys()),
                    "pf_n_records": len(pf_records),
                    "pf_fields": pf_fields,
                    "pf_tag": (pf_tag or ""),
                }
                pf_payload: Dict[str, Any] = {
                    "meta_json": np.array(json.dumps(pf_meta), dtype=np.string_)
                }

                for field in pf_fields:
                    vals = [rec.get(field, np.nan) for rec in pf_records]
                    packed, needs_json = pack_traj_field(field, vals)
                    if (not needs_json) and packed is not None:
                        pf_payload[f"pf__{field}"] = packed
                    else:
                        pf_payload[f"pf__{field}__json"] = np.array(
                            json.dumps([to_jsonable(v) for v in vals], default=str),
                            dtype=np.string_,
                        )

                np.savez_compressed(pf_out_path, **pf_payload)
                pf_written = True
                meta["pf_export_path"] = str(pf_out_path.resolve())
            except Exception as e:
                meta["pf_export_error"] = str(e)
                pf_written = False

        meta["pf_exported"] = pf_written
        payload["meta_json"] = np.array(json.dumps(meta), dtype=np.string_)

    else:
        traj = extract_traj(md)
        records = traj_to_records(traj)
        fields_all = sorted({k for rec in records for k in rec.keys()})
        if traj_fields_req:
            lookup = {f.lower(): f for f in fields_all}
            selected: list[str] = []
            missing: list[str] = []
            for req in traj_fields_req:
                hit = lookup.get(req.lower())
                if hit is None:
                    missing.append(req)
                elif hit not in selected:
                    selected.append(hit)
            fields = selected
            meta["traj_fields_requested"] = traj_fields_req
            meta["traj_fields_missing"] = missing
        else:
            fields = fields_all

        meta["traj_n_trials"] = len(records)
        meta["traj_fields"] = fields
        payload["meta_json"] = np.array(json.dumps(meta), dtype=np.string_)

        for field in fields:
            vals = [rec.get(field, np.nan) for rec in records]
            packed, needs_json = pack_traj_field(field, vals)
            if (not needs_json) and packed is not None:
                payload[f"traj__{field}"] = packed
            else:
                payload[f"traj__{field}__json"] = np.array(
                    json.dumps([to_jsonable(v) for v in vals], default=str),
                    dtype=np.string_,
                )

    if not skip_main_write:
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
    ap.add_argument(
        "--mode",
        type=str,
        choices=["ratemap", "trajdata"],
        default="ratemap",
        help="Choose conversion mode.",
    )
    ap.add_argument(
        "--input",
        type=str,
        default="",
        help=(
            "Input .mat file or directory. "
            "Default: data/raw (ratemap) or data/raw/trajdata (trajdata)."
        ),
    )
    ap.add_argument(
        "--output",
        type=str,
        default="",
        help=(
            "Output directory. "
            "Default: data/interim for both modes."
        ),
    )
    ap.add_argument("--recursive", action="store_true", help="Recursively search for .mat files.")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing .npz files.")
    ap.add_argument("--include-spikes", action="store_true", help="Include allcel.Spikes if present.")
    ap.add_argument(
        "--traj-fields",
        type=str,
        default="",
        help=(
            "Comma-separated Traj fields to export in trajdata mode. "
            "Default uses a compact whitelist. Use 'all' to export every Traj field."
        ),
    )
    ap.add_argument(
        "--pf-tag",
        type=str,
        default="",
        help=(
            "Optional suffix tag for PF export file in ratemap mode. "
            "Example: --pf-tag c writes <SESSION>_pf_c.npz."
        ),
    )
    ap.add_argument("--dry-run", action="store_true", help="Only print what would be converted.")
    args = ap.parse_args()

    if args.input:
        in_path = Path(args.input)
    else:
        in_path = Path("data/raw/trajdata") if args.mode == "trajdata" else Path("data/raw")

    if args.output:
        out_root = Path(args.output)
    else:
        out_root = Path("data/interim")
    out_root.mkdir(parents=True, exist_ok=True)

    mat_files = list(iter_mat_files(in_path, recursive=args.recursive))
    marker = "ratemap" if args.mode == "ratemap" else "trajdata"
    mat_files = [p for p in mat_files if (marker in p.name.lower() and p.suffix.lower() == ".mat")]

    if not mat_files:
        raise SystemExit(f"No '*{marker}*.mat' files found in {in_path}")

    traj_fields_req: list[str] | None = None
    if args.mode == "trajdata":
        if args.traj_fields.strip().lower() == "all":
            traj_fields_req = None
        else:
            traj_fields_req = parse_csv_list(args.traj_fields) or DEFAULT_TRAJ_FIELDS

    print(f"Found {len(mat_files)} file(s).")
    pf_tag = args.pf_tag.strip() or None
    for p in mat_files:
        info = parse_session_info(p)
        out_path = build_out_path(info, out_root=out_root, mode=args.mode)
        if args.dry_run:
            print(f"[DRY] {p} -> {out_path}")
            if args.mode == "ratemap":
                print(f"[DRY] {p} -> {build_pf_out_path(info, out_root=out_root, pf_tag=pf_tag)}")
            continue

        try:
            outp = convert_one(
                p,
                out_root=out_root,
                overwrite=args.overwrite,
                include_spikes=args.include_spikes,
                mode=args.mode,
                traj_fields_req=traj_fields_req,
                pf_tag=pf_tag,
            )
            print(f"[OK]  {p.name} -> {outp}")
        except Exception as e:
            print(f"[ERR] {p.name}: {e}")


if __name__ == "__main__":
    main()
