"""Regression: batched UMA PythonForce must hand float64 forces to OpenMM (CUDA RPMD stability).

Also guard against reintroducing the **float32 ``AtomicData.pos`` / template-mutation** batch path,
which caused large drift vs sequential ``from_ase`` (see ``test_uma_fairchem_batch_parity.py``).
"""

from pathlib import Path

import pytest

pytestmark = pytest.mark.openmmml


def test_batch_return_forces_are_float64_in_source() -> None:
    root = Path(__file__).resolve().parents[2]
    path = root / "wrappers" / "python" / "openmmml" / "models" / "umapotential_pythonforce_batch.py"
    text = path.read_text(encoding="utf-8")
    assert "forces_cpu = np.array(forces_kj_nm_list, dtype=np.float64)" in text
    assert "full_forces_buffer = np.zeros((total_particles, 3), dtype=np.float64)" in text


def test_batch_path_uses_fresh_from_ase_not_float32_pos_overwrite() -> None:
    root = Path(__file__).resolve().parents[2]
    path = root / "wrappers" / "python" / "openmmml" / "models" / "umapotential_pythonforce_batch.py"
    text = path.read_text(encoding="utf-8")
    assert "Fresh AtomicData per bead" in text or "fresh ``AtomicData.from_ase``" in text
    assert "atomic_templates[i].pos" not in text
    assert "dtype=np.float32" not in text
    assert "all_pos_angstrom[i], dtype=np.float64)" in text
