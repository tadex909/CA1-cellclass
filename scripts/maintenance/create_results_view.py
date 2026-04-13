from __future__ import annotations

import argparse
import json
import os
import re
import shutil
from pathlib import Path
from typing import Any

import pandas as pd


AGE_GROUP_RE = re.compile(r"^P\d+(?:[-_]\d+)?$")
MOUSE_RE = re.compile(r"^(?:VS|V)\d+$")


def parse_csv_list(raw: str | None) -> list[str]:
    if not raw:
        return []
    return [x.strip() for x in raw.split(",") if x.strip()]


def classify_top_level(name: str) -> tuple[str, str]:
    if name in {"ratemap", "figures"} or name.startswith("placefield_"):
        return "placefields", name
    if name in {"tables"}:
        return "reports", name
    if name in {"model_selection", "stability", "feature_set_experiments"} or name.startswith(
        "type_u_comparison"
    ):
        return "models", name
    if AGE_GROUP_RE.match(name):
        return "datasets/age_groups", name
    if MOUSE_RE.match(name):
        return "datasets/mice", name
    return "misc", name


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


def copytree_with_mode(src: Path, dst: Path, mode: str) -> str:
    if mode not in {"hardlink", "copy"}:
        raise ValueError(f"Unsupported mode: {mode}")

    def _copy_fn(src_file: str, dst_file: str) -> str:
        if mode == "hardlink":
            try:
                os.link(src_file, dst_file)
                return dst_file
            except OSError:
                return shutil.copy2(src_file, dst_file)
        return shutil.copy2(src_file, dst_file)

    shutil.copytree(src, dst, copy_function=_copy_fn)
    return "hardlinked" if mode == "hardlink" else "copied"


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Create an organized non-destructive view of results/. "
            "By default this performs a dry-run plan only."
        )
    )
    ap.add_argument("--results_root", type=str, default="results")
    ap.add_argument("--view_root", type=str, default="results/_organized")
    ap.add_argument(
        "--mode",
        type=str,
        default="hardlink",
        choices=["hardlink", "copy"],
        help=(
            "How files are materialized in the organized view. "
            "'hardlink' avoids data duplication where supported."
        ),
    )
    ap.add_argument("--execute", action="store_true", help="Apply the plan.")
    ap.add_argument("--overwrite", action="store_true", help="Replace existing destination directories.")
    ap.add_argument(
        "--skip_dirs",
        type=str,
        default="_organized,_archive",
        help="Comma-separated top-level names in results/ to skip.",
    )
    return ap


def main() -> None:
    args = build_parser().parse_args()
    results_root = Path(args.results_root)
    view_root = Path(args.view_root)
    skip_names = set(parse_csv_list(args.skip_dirs))

    if not results_root.exists():
        raise FileNotFoundError(f"results_root does not exist: {results_root}")

    plan_rows: list[dict[str, Any]] = []
    children = sorted(results_root.iterdir(), key=lambda p: p.name.lower())
    for src in children:
        if not src.is_dir():
            continue
        if src.name.startswith("."):
            continue
        if src.name in skip_names:
            continue
        # Avoid recursive self-inclusion if view_root is under results_root.
        if src.resolve() == view_root.resolve():
            continue

        category_root, leaf_name = classify_top_level(src.name)
        dst = view_root / category_root / leaf_name
        n_files, n_bytes = dir_stats(src)

        action = "planned"
        status = "ok"
        if args.execute:
            dst.parent.mkdir(parents=True, exist_ok=True)
            if dst.exists():
                if not args.overwrite:
                    action = "skip_exists"
                    status = "skipped"
                else:
                    shutil.rmtree(dst)
                    action = copytree_with_mode(src, dst, args.mode)
            else:
                action = copytree_with_mode(src, dst, args.mode)

        plan_rows.append(
            {
                "source_dir": str(src),
                "dest_dir": str(dst),
                "category": category_root,
                "name": leaf_name,
                "mode": args.mode,
                "n_files": int(n_files),
                "size_mb": float(n_bytes / 1e6),
                "action": action,
                "status": status,
            }
        )

    view_root.mkdir(parents=True, exist_ok=True)
    plan_df = pd.DataFrame(plan_rows).sort_values(["category", "name"], kind="stable")
    plan_df.to_csv(view_root / "reorg_plan.csv", index=False)
    write_json(
        view_root / "reorg_plan.json",
        [] if plan_df.empty else plan_df.to_dict(orient="records"),
    )

    summary = {
        "results_root": str(results_root),
        "view_root": str(view_root),
        "mode": args.mode,
        "execute": bool(args.execute),
        "overwrite": bool(args.overwrite),
        "n_entries": int(len(plan_df)),
        "n_skipped": int((plan_df["status"] != "ok").sum()) if not plan_df.empty else 0,
        "total_files": int(plan_df["n_files"].sum()) if not plan_df.empty else 0,
        "total_size_mb": float(plan_df["size_mb"].sum()) if not plan_df.empty else 0.0,
    }
    write_json(view_root / "reorg_summary.json", summary)

    if args.execute:
        print(f"Created organized view: {view_root}")
    else:
        print("Dry-run only. No files were copied/linked.")
        print("Re-run with --execute to materialize the organized view.")
    print(f"Wrote plan: {view_root / 'reorg_plan.csv'}")


if __name__ == "__main__":
    main()

