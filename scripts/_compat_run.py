from __future__ import annotations

import runpy
import sys
from pathlib import Path


def run_moved(relative_target: str, *, legacy_name: str) -> None:
    here = Path(__file__).resolve().parent
    target = (here / relative_target).resolve()
    if not target.exists():
        raise FileNotFoundError(f"Moved script target not found: {target}")
    print(
        f"[compat] `scripts/{legacy_name}` is deprecated. "
        f"Use `python scripts/{relative_target}`.",
        file=sys.stderr,
    )
    runpy.run_path(str(target), run_name="__main__")

