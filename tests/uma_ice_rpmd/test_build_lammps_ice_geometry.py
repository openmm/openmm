"""Geometry sanity for LAMMPS ice data built with -n (no CIF truncation)."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

_SCRIPT_DIR = Path(__file__).resolve().parent
_LAMMPS_DIR = _SCRIPT_DIR / "lammps"


def _min_oo_angstrom_from_data(path: Path) -> float:
    from benchmark_same_structure_openmm_vs_lammps import parse_lammps_data

    pos, cell, z = parse_lammps_data(path)
    n_mol = z.shape[0] // 3
    box = np.array([[cell[0, 0], 0.0, 0.0], [0.0, cell[1, 1], 0.0], [0.0, 0.0, cell[2, 2]]])

    def mic(dr: np.ndarray) -> np.ndarray:
        out = dr.copy()
        for k in range(3):
            l_ = box[k, k]
            out[k] -= np.round(out[k] / l_) * l_
        return out

    o_idx = [3 * i for i in range(n_mol)]
    min_d = 1e9
    for ii, i in enumerate(o_idx):
        for j in o_idx[ii + 1 :]:
            d = np.linalg.norm(mic(pos[j] - pos[i]))
            min_d = min(min_d, d)
    return float(min_d)


@pytest.mark.parametrize("n_mol", [32, 64])
def test_build_lammps_ice_min_oo_sane(n_mol: int, tmp_path: Path) -> None:
    """Built ice must have O–O ≥ 2.4 Å (MIC); truncated CIF gives ~1.7 Å."""
    out = tmp_path / f"data_{n_mol}"
    subprocess.run(
        [sys.executable, str(_LAMMPS_DIR / "build_lammps_ice_data.py"), "-n", str(n_mol), "-o", str(out)],
        cwd=str(_SCRIPT_DIR),
        check=True,
    )
    d_min = _min_oo_angstrom_from_data(out)
    # Truncated CIF (first-N molecules) gives ~1.7 Å; valid ice Ih ~2.76 Å; small ice Ic
    # supercells can show ~2.2 Å shortest O–O while remaining equilibratable.
    assert d_min >= 2.15, f"min O–O {d_min:.3f} Å too small (likely bad truncation / clash)"
