"""CLI smoke: i-PI orchestrator documents UMA minimization flags (OpenMM parity)."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

_SCRIPT = Path(__file__).resolve().parent / "run_ipi_lammps_uma_rpmd.py"


def test_help_lists_ipi_minimize_options() -> None:
    out = subprocess.run(
        [sys.executable, str(_SCRIPT), "--help"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert "--no-ipi-minimize" in out.stdout
    assert "--ipi-min-maxiter" in out.stdout
    assert "--ipi-min-maxeval" in out.stdout
    assert "--order-thermo-prefer-ring-temperature" in out.stdout

    order_script = Path(__file__).resolve().parent / "ipi_order_from_traj.py"
    out2 = subprocess.run(
        [sys.executable, str(order_script), "--help"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert "--thermo-prefer-ring-temperature" in out2.stdout
