from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd


AGE_GROUP_RE = re.compile(r"^P\d+(?:[-_]\d+)?$")
MOUSE_RE = re.compile(r"^(?:VS|V)\d+$")


def parse_csv_list(raw: str | None) -> list[str]:
    if not raw:
        return []
    return [x.strip() for x in raw.split(",") if x.strip()]


def dir_stats(path: Path) -> tuple[int, int]:
    n_files = 0
    n_bytes = 0
    for p in path.rglob("*"):
        if p.is_file():
            n_files += 1
            try:
                n_bytes += p.stat().st_size
            except OSError:
                continue
    return n_files, n_bytes


def classify_top_level(name: str) -> tuple[str, str]:
    if name in {"ratemap", "figures"} or name.startswith("placefield_"):
        return "placefields", "session_outputs"
    if name in {"tables"}:
        return "reports", "tables"
    if name in {"model_selection", "stability", "feature_set_experiments"} or name.startswith(
        "type_u_comparison"
    ):
        return "models", "analysis_outputs"
    if AGE_GROUP_RE.match(name):
        return "datasets", "age_group"
    if MOUSE_RE.match(name):
        return "datasets", "mouse"
    return "misc", "other"


def safe_read_json(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def collect_top_level_inventory(results_root: Path, skip_names: set[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for child in sorted(results_root.iterdir(), key=lambda p: p.name.lower()):
        if not child.is_dir():
            continue
        if child.name.startswith("."):
            continue
        if child.name in skip_names:
            continue
        category, subcategory = classify_top_level(child.name)
        n_files, n_bytes = dir_stats(child)
        rows.append(
            {
                "name": child.name,
                "path": str(child),
                "category": category,
                "subcategory": subcategory,
                "n_files": int(n_files),
                "size_mb": float(n_bytes / 1e6),
            }
        )
    return pd.DataFrame(rows)


def read_status_counts(run_index_csv: Path) -> dict[str, int]:
    out = {"ok": 0, "error": 0, "skip_exists": 0, "other": 0, "rows": 0}
    if not run_index_csv.exists():
        return out
    try:
        d = pd.read_csv(run_index_csv)
    except Exception:
        return out
    out["rows"] = int(len(d))
    if "status" not in d.columns:
        return out
    vc = d["status"].astype(str).value_counts(dropna=False).to_dict()
    out["ok"] = int(vc.get("ok", 0))
    out["error"] = int(vc.get("error", 0))
    out["skip_exists"] = int(vc.get("skip_exists", 0))
    used = {"ok", "error", "skip_exists"}
    out["other"] = int(sum(v for k, v in vc.items() if k not in used))
    return out


def collect_run_registry(results_root: Path, skip_names: set[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for cfg_path in sorted(results_root.rglob("run_config.json")):
        if any(name in skip_names for name in cfg_path.parts):
            continue
        run_dir = cfg_path.parent
        cfg = safe_read_json(cfg_path)
        script = str(cfg.get("script") or cfg.get("runner") or "")
        if not script:
            script = str(run_dir.relative_to(results_root))
        args = cfg.get("args") if isinstance(cfg.get("args"), dict) else {}
        run_id = str(cfg.get("run_id") or run_dir.name)
        created_at = cfg.get("created_at_utc")

        run_index_csv = run_dir / "run_index.csv"
        status = read_status_counts(run_index_csv)
        class_csv = run_dir / "ssi_classification.csv"
        if isinstance(cfg.get("classification_csv_resolved"), str):
            class_csv = Path(cfg["classification_csv_resolved"])

        n_files, n_bytes = dir_stats(run_dir)
        rows.append(
            {
                "run_id": run_id,
                "run_dir": str(run_dir),
                "script": script,
                "created_at_utc": created_at,
                "git_head": cfg.get("git_head"),
                "null_method": args.get("null_method"),
                "smooth_sigma_bins": args.get("smooth_sigma_bins"),
                "xbin_rem": args.get("xbin_rem"),
                "sm_alpha": args.get("sm_alpha"),
                "nb_rep": args.get("nb_rep"),
                "run_subdir_enabled": cfg.get("run_subdir_enabled"),
                "has_run_index": bool(run_index_csv.exists()),
                "run_index_csv": str(run_index_csv),
                "n_rows_index": int(status["rows"]),
                "n_ok": int(status["ok"]),
                "n_error": int(status["error"]),
                "n_skip_exists": int(status["skip_exists"]),
                "n_other_status": int(status["other"]),
                "has_ssi_classification": bool(class_csv.exists()),
                "ssi_classification_csv": str(class_csv),
                "n_files": int(n_files),
                "size_mb": float(n_bytes / 1e6),
            }
        )
    return pd.DataFrame(rows)


def collect_latest_pointers(results_root: Path, skip_names: set[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for p in sorted(results_root.rglob("LATEST_RUN.txt")):
        if any(name in skip_names for name in p.parts):
            continue
        try:
            run_id = p.read_text(encoding="utf-8").strip()
        except Exception:
            run_id = ""
        base = p.parent
        resolved = base / run_id if run_id else base
        rows.append(
            {
                "latest_pointer": str(p),
                "base_dir": str(base),
                "run_id": run_id,
                "resolved_dir": str(resolved),
                "resolved_exists": bool(resolved.exists()),
            }
        )
    return pd.DataFrame(rows)


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Build an inventory/registry of results outputs and analysis runs."
    )
    ap.add_argument("--results_root", type=str, default="results")
    ap.add_argument("--out_root", type=str, default="results/tables/results_registry")
    ap.add_argument(
        "--skip_dirs",
        type=str,
        default="_organized,_archive",
        help="Comma-separated top-level names to skip when scanning.",
    )
    return ap


def main() -> None:
    args = build_parser().parse_args()
    results_root = Path(args.results_root)
    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    skip_names = set(parse_csv_list(args.skip_dirs))

    top = collect_top_level_inventory(results_root, skip_names)
    runs = collect_run_registry(results_root, skip_names)
    latest = collect_latest_pointers(results_root, skip_names)

    if not top.empty:
        top = top.sort_values(["category", "subcategory", "name"], kind="stable").reset_index(drop=True)
    if not runs.empty:
        runs = runs.sort_values(["script", "created_at_utc", "run_dir"], kind="stable").reset_index(drop=True)
    if not latest.empty:
        latest = latest.sort_values(["base_dir"], kind="stable").reset_index(drop=True)

    top.to_csv(out_root / "top_level_inventory.csv", index=False)
    runs.to_csv(out_root / "run_registry.csv", index=False)
    latest.to_csv(out_root / "latest_pointers.csv", index=False)

    write_json(
        out_root / "top_level_inventory.json",
        [] if top.empty else top.to_dict(orient="records"),
    )
    write_json(
        out_root / "run_registry.json",
        [] if runs.empty else runs.to_dict(orient="records"),
    )
    write_json(
        out_root / "latest_pointers.json",
        [] if latest.empty else latest.to_dict(orient="records"),
    )

    summary = {
        "results_root": str(results_root),
        "out_root": str(out_root),
        "n_top_level_dirs": int(len(top)),
        "n_run_configs": int(len(runs)),
        "n_latest_pointers": int(len(latest)),
        "top_level_category_counts": (
            {} if top.empty else {str(k): int(v) for k, v in top["category"].value_counts().to_dict().items()}
        ),
        "scripts_with_runs": (
            {} if runs.empty else {str(k): int(v) for k, v in runs["script"].value_counts().to_dict().items()}
        ),
    }
    write_json(out_root / "summary.json", summary)

    print(f"Wrote registry to: {out_root}")
    print(f"Top-level directories indexed: {summary['n_top_level_dirs']}")
    print(f"Runs indexed (run_config.json): {summary['n_run_configs']}")


if __name__ == "__main__":
    main()
