from __future__ import annotations

import sys
from pathlib import Path


def _ensure_src_on_path() -> None:
    here = Path(__file__).resolve().parent
    root = here
    while root != root.parent and not (root / "src" / "cellclass").is_dir():
        root = root.parent
    src_dir = root / "src"
    if not (src_dir / "cellclass").is_dir():
        raise RuntimeError(f"Could not find src/cellclass starting from {here}")
    sys.path.insert(0, str(src_dir))


def main() -> None:
    _ensure_src_on_path()
    from cellclass.pipeline import main as pipeline_main

    print(
        "[compat] `scripts/legacy/one_file_processing.py` is deprecated. "
        "Use `python -m cellclass.pipeline` or "
        "`python scripts/pipelines/interim_to_processed.py`.",
        file=sys.stderr,
    )
    pipeline_main()


if __name__ == "__main__":
    main()
