#!/usr/bin/env python3
"""
Diagnose why the i-PI + Python UMA driver produces catastrophically wrong
RPMD dynamics (O-H bonds break from 0.96 to 3.9 A in 100 steps).

Three force-evaluation paths on the same initial ice structure:

  Path A  --  Exact i-PI driver path (template_ase.copy(), get_properties)
  Path B  --  Direct ASE reference   (get_potential_energy / get_forces)
  Path C  --  OpenMM UMA PythonForce (from compare_forces_ice.py helpers)

Also tests 32 sequential calculator calls with perturbed positions to rule
out stale-cache or GPU-state corruption.

Usage:
  cd tests/uma_ice_rpmd
  python diagnose_ipi_forces.py
"""
from __future__ import annotations

import os
import re
import sys
import time
from pathlib import Path

import numpy as np

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

if os.getenv("OPENMM_PLUGIN_DIR") is None:
    for _base in [os.getenv("CONDA_PREFIX"), os.path.expanduser("~/miniconda3")]:
        if _base:
            _plugins = os.path.join(_base, "lib", "plugins")
            if os.path.isdir(_plugins):
                os.environ["OPENMM_PLUGIN_DIR"] = _plugins
                break

INIT_XYZ = _SCRIPT_DIR / "ipi" / "init.xyz"
MODEL = "uma-s-1p1"

_CELL_RE = re.compile(
    r"CELL\(abcABC\):\s+"
    r"([\d.eE+\-]+)\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)"
    r"\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)"
)


def _read_ipi_cell_angstrom(xyz_path: str | Path) -> np.ndarray:
    with open(xyz_path) as fh:
        fh.readline()
        comment = fh.readline()
    m = _CELL_RE.search(comment)
    if m is None:
        raise ValueError(f"No CELL(abcABC) in {xyz_path}")
    a, b, c = float(m.group(1)), float(m.group(2)), float(m.group(3))
    alpha = np.deg2rad(float(m.group(4)))
    beta = np.deg2rad(float(m.group(5)))
    gamma = np.deg2rad(float(m.group(6)))
    cx = c * np.cos(beta)
    cy = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    cz = np.sqrt(max(c**2 - cx**2 - cy**2, 0.0))
    return np.array([
        [a, 0.0, 0.0],
        [b * np.cos(gamma), b * np.sin(gamma), 0.0],
        [cx, cy, cz],
    ])


def _load_init_structure():
    """Return (positions_ang, cell_ang_3x3, symbols) from ipi/init.xyz."""
    from ase.io import read as ase_read

    atoms = ase_read(str(INIT_XYZ))
    cell = _read_ipi_cell_angstrom(INIT_XYZ)
    atoms.set_cell(cell)
    atoms.set_pbc(True)
    return np.array(atoms.positions), cell, list(atoms.get_chemical_symbols())


def _make_calculator():
    """Build a fresh FAIRChemCalculator (shared by paths A/B)."""
    try:
        from fairchem.core import FAIRChemCalculator
        from fairchem.core.calculate import pretrained_mlip
    except ImportError:
        from fairchem.core.calculate import pretrained_mlip
        from fairchem.core.calculate.ase_calculator import FAIRChemCalculator

    import torch

    dev = "cuda" if torch.cuda.is_available() else "cpu"
    predictor = pretrained_mlip.get_predict_unit(MODEL, device=dev)
    return FAIRChemCalculator(predictor, task_name="omol")


def path_a_ipi_driver(pos_ang, cell_ang, symbols, calc):
    """Reproduce compute_structure: template.copy(), set pos/cell, get_properties."""
    from ase import Atoms

    template = Atoms(symbols=symbols, positions=pos_ang)
    template.set_cell(cell_ang)
    template.set_pbc(True)
    template.info["charge"] = 0
    template.info["spin"] = 1

    structure = template.copy()
    structure.positions[:] = pos_ang
    structure.cell[:] = cell_ang
    structure.pbc = [True, True, True]
    structure.info["charge"] = 0
    structure.info["spin"] = 1
    structure.calc = calc

    props = structure.get_properties(["energy", "forces"])
    e = float(props["energy"])
    f = np.asarray(props["forces"], dtype=np.float64).copy()
    return e, f


def path_b_direct_ase(pos_ang, cell_ang, symbols, calc):
    """Standard ASE path: get_potential_energy / get_forces."""
    from ase import Atoms
    from ase.geometry import wrap_positions

    wrapped = wrap_positions(pos_ang, cell_ang, pbc=[True, True, True])
    atoms = Atoms(symbols=symbols, positions=wrapped, cell=cell_ang, pbc=True)
    atoms.info["charge"] = 0
    atoms.info["spin"] = 1
    atoms.calc = calc

    e = float(atoms.get_potential_energy())
    f = np.asarray(atoms.get_forces(), dtype=np.float64).copy()
    return e, f


def _report_forces(label, e, f):
    max_f = np.max(np.abs(f))
    rms_f = np.sqrt(np.mean(np.sum(f**2, axis=1)))
    has_nan = np.any(np.isnan(f))
    has_inf = np.any(np.isinf(f))
    print(f"  [{label}]  E = {e:.6f} eV   max|F| = {max_f:.4f} eV/A"
          f"   RMS|F| = {rms_f:.4f} eV/A   NaN={has_nan}  Inf={has_inf}")


def _compare(label_a, ea, fa, label_b, eb, fb):
    dE = ea - eb
    df = fa - fb
    max_df = np.max(np.abs(df))
    rms_df = np.sqrt(np.mean(np.sum(df**2, axis=1)))
    per_atom = np.max(np.abs(df), axis=1)
    worst = int(np.argmax(per_atom))
    print(f"  {label_a} vs {label_b}:  ΔE = {dE:+.6f} eV  max|ΔF| = {max_df:.6f} eV/A"
          f"  RMS|ΔF| = {rms_df:.6f}  worst atom {worst}")


def test_sequential_calls(pos_ang, cell_ang, symbols, calc, n_calls=32):
    """Call the calculator n_calls times with small random perturbations.

    Returns True if all energies differ and no NaN/Inf detected.
    """
    from ase import Atoms

    rng = np.random.default_rng(42)
    energies = []
    force_norms = []
    ok = True

    for i in range(n_calls):
        perturbed = pos_ang + rng.normal(scale=0.01, size=pos_ang.shape)
        atoms = Atoms(symbols=symbols, positions=perturbed, cell=cell_ang, pbc=True)
        atoms.info["charge"] = 0
        atoms.info["spin"] = 1
        atoms.calc = calc

        props = atoms.get_properties(["energy", "forces"])
        e = float(props["energy"])
        f = np.asarray(props["forces"], dtype=np.float64)
        max_f = np.max(np.abs(f))

        if np.any(np.isnan(f)) or np.any(np.isinf(f)):
            print(f"    call {i}: E={e:.6f}  max|F|={max_f:.4f}  *** NaN/Inf ***")
            ok = False
        energies.append(e)
        force_norms.append(max_f)

    energies = np.array(energies)
    force_norms = np.array(force_norms)

    n_unique = len(np.unique(np.round(energies, 8)))
    all_different = n_unique == n_calls
    e_spread = energies.max() - energies.min()
    f_spread = force_norms.max() - force_norms.min()

    print(f"  {n_calls} sequential calls: "
          f"E range = {e_spread:.6f} eV   "
          f"max|F| range = [{force_norms.min():.4f}, {force_norms.max():.4f}] eV/A   "
          f"unique E = {n_unique}/{n_calls}")

    if not all_different:
        print("  *** WARNING: some energies are identical -> possible stale cache ***")
        ok = False

    return ok


def main():
    print("=" * 70)
    print("i-PI UMA driver force diagnostic")
    print("=" * 70)

    pos_ang, cell_ang, symbols = _load_init_structure()
    n_atoms = len(pos_ang)
    print(f"Structure: {INIT_XYZ}  ({n_atoms} atoms)")
    print(f"Cell diagonal: {np.diag(cell_ang)}")
    print()

    print("[1] Building FAIRChemCalculator ...")
    calc = _make_calculator()
    print()

    print("[2] Path A: i-PI driver path (template.copy + get_properties) ...")
    t0 = time.time()
    ea, fa = path_a_ipi_driver(pos_ang, cell_ang, symbols, calc)
    dt_a = time.time() - t0
    _report_forces("Path A", ea, fa)
    print(f"      time = {dt_a:.2f}s")
    print()

    print("[3] Path B: direct ASE (get_potential_energy + get_forces) ...")
    t0 = time.time()
    eb, fb = path_b_direct_ase(pos_ang, cell_ang, symbols, calc)
    dt_b = time.time() - t0
    _report_forces("Path B", eb, fb)
    print(f"      time = {dt_b:.2f}s")
    print()

    print("[4] Comparison A vs B ...")
    _compare("Path A", ea, fa, "Path B", eb, fb)
    print()

    print("[5] Path C: OpenMM UMA PythonForce ...")
    try:
        import test_uma_ice_rpmd as _ice_mod
        from openmm import app, unit as omm_unit, Vec3, Context, Platform, VerletIntegrator
        from openmmml import MLPotential

        topology, positions_nm, box_vectors = _ice_mod._load_structure_from_file(
            str(INIT_XYZ), max_molecules=None
        )
        pos_nm_arr = np.array([[p[0], p[1], p[2]] for p in positions_nm])
        from ase.geometry import wrap_positions
        cell_omm = cell_ang
        pos_omm_ang = pos_nm_arr * 10.0
        wrapped_ang = wrap_positions(pos_omm_ang, cell_omm, pbc=[True, True, True])
        wrapped_nm = [(wrapped_ang[i, 0] * 0.1, wrapped_ang[i, 1] * 0.1,
                       wrapped_ang[i, 2] * 0.1) for i in range(n_atoms)]

        platform = Platform.getPlatformByName("CUDA") if _has_cuda_omm() else Platform.getPlatformByName("CPU")
        ml_dev = "cuda" if platform.getName() == "CUDA" else "cpu"
        potential = MLPotential("uma-s-1p1-pythonforce-batch")
        system = potential.createSystem(
            topology,
            task_name="omol",
            charge=0,
            spin=1,
            device=ml_dev,
            use_atom_wrap_for_lammps_parity=True,
        )
        integrator = VerletIntegrator(1.0 * omm_unit.femtoseconds)
        context = Context(system, integrator, platform)
        context.setPositions([Vec3(p[0], p[1], p[2]) * omm_unit.nanometer for p in wrapped_nm])
        context.setPeriodicBoxVectors(*box_vectors)

        state = context.getState(getEnergy=True, getForces=True)
        ec_kj = state.getPotentialEnergy().value_in_unit(omm_unit.kilojoules_per_mole)
        fc_kj_nm = state.getForces(asNumpy=True).value_in_unit(
            omm_unit.kilojoules_per_mole / omm_unit.nanometer)

        EV_TO_KJ = 96.4853
        ec_ev = ec_kj / EV_TO_KJ
        fc_ev_ang = np.asarray(fc_kj_nm) / (EV_TO_KJ * 10.0)

        _report_forces("Path C", ec_ev, fc_ev_ang)
        _compare("Path A", ea, fa, "Path C", ec_ev, fc_ev_ang)
        _compare("Path B", eb, fb, "Path C", ec_ev, fc_ev_ang)

    except Exception as exc:
        print(f"  (Path C skipped: {exc})")
    print()

    print("[6] Sequential call test (32 calls, small perturbations) ...")
    seq_ok = test_sequential_calls(pos_ang, cell_ang, symbols, calc, n_calls=32)
    print()

    print("=" * 70)
    if np.max(np.abs(fa - fb)) < 0.01 and np.max(np.abs(fa)) > 0.1:
        print("CONCLUSION: Forces match between Path A and B and are non-trivial.")
        print("            The i-PI driver force computation is likely correct.")
        print("            Issue may be in i-PI dynamics / thermostat / communication.")
    elif np.max(np.abs(fa)) < 0.001:
        print("CONCLUSION: Path A forces are near-zero! Calculator returning empty results.")
    elif np.max(np.abs(fa - fb)) > 0.1:
        print("CONCLUSION: Forces DIFFER between Path A (get_properties) and B (get_forces)!")
        print("            The ASE get_properties path is returning wrong results.")
    print("=" * 70)


def _has_cuda_omm():
    try:
        from openmm import Platform
        Platform.getPlatformByName("CUDA")
        return True
    except Exception:
        return False


if __name__ == "__main__":
    main()
