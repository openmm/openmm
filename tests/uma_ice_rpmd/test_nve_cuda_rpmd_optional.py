"""Opt-in CUDA vs Reference RPMD NVE drift (microcanonical, no thermostat).

Run when GPU + CPU reference are acceptable::

  export OPENMM_PLUGIN_DIR="${CONDA_PREFIX}/lib/plugins"
  export OPENMM_NVE_CUDA_SMOKE=1
  pytest tests/uma_ice_rpmd/test_nve_cuda_rpmd_optional.py -q

Uses 2×2×2 ice (64 molecules) to match ``ipi/init.xyz`` when present, 4 beads,
and short production. Parses ``NVE_diag`` lines from ``test_uma_ice_rpmd`` output.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

_SCRIPT_DIR = Path(__file__).resolve().parent


def _require_pythonforce() -> None:
    import openmm

    if not hasattr(openmm, "PythonForce"):
        pytest.skip(
            "OpenMM without PythonForce in this Python (use a project/ML-enabled OpenMM build)."
        )

_NVE_DIAG_RE = re.compile(
    r"NVE_diag:\s*sum_bead\(PE\+KE\)=([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)"
)


def _parse_nve_diag_energies_kjmol(text: str) -> list[float]:
    return [float(m.group(1)) for m in _NVE_DIAG_RE.finditer(text)]


def test_parse_nve_diag_regex() -> None:
    text = (
        "header\n"
        "         NVE_diag: sum_bead(PE+KE)=-12845000.123456 kJ/mol (approximate conserved quantity)\n"
    )
    vals = _parse_nve_diag_energies_kjmol(text)
    assert vals == [-12845000.123456]


def _drift_rate_kj_per_mol_ps(energies: list[float], n_molecules: int, prod_ps: float) -> float:
    if len(energies) < 2:
        raise ValueError("need at least two NVE_diag samples")
    de = abs(energies[-1] - energies[0])
    return float(de) / max(float(n_molecules), 1.0) / max(prod_ps, 1e-15)


def _base_cmd_common(
    *,
    out_csv: Path,
    platform: str,
    prod_ps: float,
    precision: str | None,
) -> list[str]:
    cmd = [
        sys.executable,
        str(_SCRIPT_DIR / "run_openmm_rpmd_reference.py"),
        "--prod",
        str(prod_ps),
        "--dt",
        "0.1",
        "--beads",
        "4",
        "--nx",
        "2",
        "--ny",
        "2",
        "--nz",
        "2",
        "--platform",
        platform,
        "--rpmd-thermostat",
        "none",
        "--equil",
        "0",
        "--order-every",
        "4",
        "--order-csv",
        str(out_csv),
        "--report-every-steps",
        "10",
    ]
    if precision:
        cmd.extend(["--precision", precision])
    return cmd


@pytest.mark.skipif(
    os.environ.get("OPENMM_NVE_CUDA_SMOKE") != "1",
    reason="Set OPENMM_NVE_CUDA_SMOKE=1 to run (needs fairchem, OPENMM_PLUGIN_DIR for CUDA).",
)
def test_cuda_rpmd_nve_matches_reference_drift_order() -> None:
    _require_pythonforce()
    if not os.environ.get("OPENMM_PLUGIN_DIR", "").strip():
        pytest.skip("OPENMM_PLUGIN_DIR not set")

    prod_ps = 0.004
    n_mol = 64
    env = os.environ.copy()
    env.setdefault("PYTHONUNBUFFERED", "1")
    env.setdefault("OPENMMML_UMA_RPMD_CHUNK", "4")

    out_dir = _SCRIPT_DIR / "pipeline_out"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_cuda = out_dir / "nve_cuda_pytest_smoke.csv"
    out_ref = out_dir / "nve_reference_pytest_smoke.csv"

    cmd_ref = _base_cmd_common(out_csv=out_ref, platform="reference", prod_ps=prod_ps, precision=None)
    cmd_ref.extend(["--ml-device", "cpu"])
    proc_ref = subprocess.run(
        cmd_ref,
        cwd=str(_SCRIPT_DIR),
        env=env,
        capture_output=True,
        text=True,
        timeout=7200,
    )
    if proc_ref.returncode != 0:
        pytest.fail(
            "Reference NVE run failed:\n"
            + proc_ref.stdout
            + "\n"
            + proc_ref.stderr
        )
    ref_text = proc_ref.stdout + "\n" + proc_ref.stderr
    e_ref = _parse_nve_diag_energies_kjmol(ref_text)
    assert len(e_ref) >= 2, f"expected NVE_diag lines in reference output; got:\n{ref_text[-4000:]}"
    rate_ref = _drift_rate_kj_per_mol_ps(e_ref, n_mol, prod_ps)

    cmd_cuda = _base_cmd_common(
        out_csv=out_cuda, platform="cuda", prod_ps=prod_ps, precision="double"
    )
    proc_cuda = subprocess.run(
        cmd_cuda,
        cwd=str(_SCRIPT_DIR),
        env=env,
        capture_output=True,
        text=True,
        timeout=7200,
    )
    if proc_cuda.returncode != 0:
        pytest.fail(
            "CUDA NVE run failed:\n" + proc_cuda.stdout + "\n" + proc_cuda.stderr
        )
    cuda_text = proc_cuda.stdout + "\n" + proc_cuda.stderr
    e_cuda = _parse_nve_diag_energies_kjmol(cuda_text)
    assert len(e_cuda) >= 2, f"expected NVE_diag lines in CUDA output; got:\n{cuda_text[-4000:]}"
    rate_cuda = _drift_rate_kj_per_mol_ps(e_cuda, n_mol, prod_ps)

    # Reference sets the scale; CUDA should not show orders-of-magnitude worse drift
    # from fixed-point / GPU integration (loose bound for heterogeneous stacks).
    ceiling = max(250.0, 20.0 * max(rate_ref, 1e-9))
    assert rate_cuda <= ceiling, (
        f"CUDA NVE drift rate {rate_cuda:.4f} kJ/(mol·ps) exceeds "
        f"max({250.0}, 20*ref) = {ceiling:.4f} (reference rate {rate_ref:.4f})"
    )
