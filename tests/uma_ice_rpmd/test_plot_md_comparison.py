"""Smoke test: plot_md_comparison runs when pipeline CSVs exist."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent


def test_plot_md_comparison_runs() -> None:
    out = _SCRIPT_DIR / "pipeline_out" / "_test_md_order_plot.png"
    r = subprocess.run(
        [sys.executable, str(_SCRIPT_DIR / "plot_md_comparison.py"), "-o", str(out)],
        cwd=str(_SCRIPT_DIR),
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert r.returncode == 0, r.stderr + r.stdout
    assert out.is_file()
    out.unlink(missing_ok=True)
