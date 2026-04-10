"""FairChem-only parity: batched ``predict`` vs single-graph ``predict`` for UMA.

UMA weights are typically float32; we expect batch vs single to match within **float32
noise** (~1e-5–1e-4 eV/Å on forces) when inputs are built the same way as production
(``r_edges=False``, fresh ``AtomicData.from_ase`` per graph, float64 NumPy positions).

Run GPU checks optionally::

  OPENMMML_FAIRCHEM_PARITY_DEVICE=cuda pytest tests/openmmml/test_uma_fairchem_batch_parity.py -q
"""

from __future__ import annotations

import os

import numpy as np
import pytest
import torch

pytest.importorskip("fairchem.core")
from fairchem.core.calculate import pretrained_mlip  # noqa: E402
from fairchem.core.datasets.atomic_data import AtomicData, atomicdata_list_to_batch  # noqa: E402
from ase import Atoms  # noqa: E402

pytestmark = [pytest.mark.openmmml, pytest.mark.fairchem, pytest.mark.slow]


def _parity_device() -> str:
    d = os.environ.get("OPENMMML_FAIRCHEM_PARITY_DEVICE", "cpu").strip().lower()
    if d == "cuda" and not torch.cuda.is_available():
        pytest.skip("OPENMMML_FAIRCHEM_PARITY_DEVICE=cuda but CUDA not available")
    return d


def _small_water_system(*, rng: np.random.Generator) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Three waters in Å, orthorhombic box."""
    symbols = (["O", "H", "H"] * 3)[:]
    lx, ly, lz = 12.0, 12.0, 12.0
    cell = np.diag([lx, ly, lz]).astype(np.float64)
    pos = rng.standard_normal((9, 3)) * 0.5 + np.array([[2.0, 2.0, 2.0]] * 9)
    pos = pos % np.array([lx, ly, lz])
    return symbols, pos, cell


def _from_ase_omol(
    symbols: list[str],
    pos_angstrom: np.ndarray,
    cell_3x3: np.ndarray | None,
    *,
    task_name: str = "omol",
    charge: int = 0,
    spin: int = 1,
) -> AtomicData:
    pbc = cell_3x3 is not None
    atoms = Atoms(
        symbols=symbols,
        positions=np.ascontiguousarray(pos_angstrom, dtype=np.float64),
        pbc=pbc,
    )
    atoms.info["charge"] = charge
    atoms.info["spin"] = spin
    if cell_3x3 is not None:
        atoms.set_cell(np.ascontiguousarray(cell_3x3, dtype=np.float64))
    return AtomicData.from_ase(
        atoms,
        task_name=task_name,
        r_edges=False,
        r_data_keys=["spin", "charge"],
    )


@pytest.fixture(scope="module")
def predict_unit_uma():
    device = _parity_device()
    return pretrained_mlip.get_predict_unit("uma-s-1p1", device=device), device


def test_batched_identical_graphs_match_single_forces(predict_unit_uma) -> None:
    """N fresh ``from_ase`` graphs with identical positions → batch forces match single."""
    predict_unit, device = predict_unit_uma
    rng = np.random.default_rng(0)
    symbols, pos, cell = _small_water_system(rng=rng)
    d0 = _from_ase_omol(symbols, pos, cell)
    d0.sid = ["bead-0"]
    if device == "cuda":
        d0 = d0.to(device)

    pred_s = predict_unit.predict(d0)
    f_s = pred_s["forces"].detach().cpu().numpy()

    n_copy = 4
    data_list = []
    for i in range(n_copy):
        di = _from_ase_omol(symbols, pos, cell)
        di.sid = [f"bead-{i}"]
        if device == "cuda":
            di = di.to(device)
        data_list.append(di)

    batch = atomicdata_list_to_batch(data_list)
    if device == "cuda":
        batch = batch.to(device)
    pred_b = predict_unit.predict(batch)
    f_b = pred_b["forces"].detach().cpu().numpy()
    bidx = batch.batch.cpu().numpy()

    for i in range(n_copy):
        fi = f_b[bidx == i]
        np.testing.assert_allclose(
            fi,
            f_s,
            rtol=5e-4,
            atol=5e-5,
            err_msg=f"batched graph {i} forces differ from single-graph predict",
        )


def test_batched_vs_loop_sequential_predictions(predict_unit_uma) -> None:
    """``atomicdata_list_to_batch`` vs N× ``predict`` on distinct small perturbations."""
    predict_unit, device = predict_unit_uma
    rng = np.random.default_rng(1)
    symbols, pos0, cell = _small_water_system(rng=rng)

    n_graphs = 3
    singles: list[AtomicData] = []
    for k in range(n_graphs):
        pos = pos0 + 0.02 * rng.standard_normal(pos0.shape)
        dk = _from_ase_omol(symbols, pos, cell)
        dk.sid = [f"g-{k}"]
        if device == "cuda":
            dk = dk.to(device)
        singles.append(dk)

    f_loop = []
    e_loop = []
    for dk in singles:
        pk = predict_unit.predict(dk)
        f_loop.append(pk["forces"].detach().cpu().numpy())
        e_loop.append(float(pk["energy"].detach().cpu().item()))

    batch = atomicdata_list_to_batch(singles)
    if device == "cuda":
        batch = batch.to(device)
    pred_b = predict_unit.predict(batch)
    f_b = pred_b["forces"].detach().cpu().numpy()
    bidx = batch.batch.cpu().numpy()
    e_b = pred_b["energy"].detach().cpu().numpy()

    for k in range(n_graphs):
        fk = f_b[bidx == k]
        np.testing.assert_allclose(
            fk,
            f_loop[k],
            rtol=5e-4,
            atol=5e-5,
            err_msg=f"loop vs batch forces differ for graph {k}",
        )
    np.testing.assert_allclose(
        np.sum(e_b),
        float(np.sum(e_loop)),
        rtol=1e-4,
        atol=1e-5,
        err_msg="batched total energy vs sum of singles",
    )


def test_float32_pos_tensor_increases_drift_vs_float64_template(predict_unit_uma) -> None:
    """Characterization: mutating ``AtomicData.pos`` to float32 can diverge from fresh from_ase."""
    predict_unit, device = predict_unit_uma
    rng = np.random.default_rng(2)
    symbols, pos, cell = _small_water_system(rng=rng)
    d_ref = _from_ase_omol(symbols, pos, cell)
    if device == "cuda":
        d_ref = d_ref.to(device)
    f_ref = predict_unit.predict(d_ref)["forces"].detach().cpu().numpy()

    d_mut = _from_ase_omol(symbols, pos, cell)
    if device == "cuda":
        d_mut = d_mut.to(device)
    # Only downgrade positions to float32 (production batch path did this); leave cell as from_ase.
    d_mut.pos = torch.from_numpy(np.ascontiguousarray(pos, dtype=np.float32)).to(
        device=d_mut.pos.device,
        dtype=torch.float32,
    )

    f_mut = predict_unit.predict(d_mut)["forces"].detach().cpu().numpy()
    max_abs = float(np.max(np.abs(f_mut - f_ref)))
    # With float32 pos, drift should exceed tight float64 noise (order 1e-6 eV/Å or larger).
    assert max_abs > 1e-7, (
        "expected float32 .pos mutation to measurably change forces vs float64 from_ase; "
        f"max |ΔF|={max_abs:.3e} eV/Å"
    )
