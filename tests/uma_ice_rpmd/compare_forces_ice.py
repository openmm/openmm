#!/usr/bin/env python3
"""
Compare single-point forces and energy: OpenMM-UMA vs ASE-FAIRChem on ice structure.

UMA has no native LAMMPS interface. ASE + FAIRChemCalculator is the canonical reference
(per FAIRChem docs and tests/uma_ice_rpmd/README.md). This script helps debug force
discrepancies that could cause ice to disintegrate during OpenMM simulation.

Usage:
  python compare_forces_ice.py --input ice.cif
  python compare_forces_ice.py --input ice.cif --molecules 16   # subset for fast debug

Output: Energy diff, max |ΔF|, RMS |ΔF|, per-atom stats. Fails if discrepancy too large.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

if os.getenv('OPENMM_PLUGIN_DIR') is None:
    for _base in [os.getenv('CONDA_PREFIX'), os.path.expanduser('~/miniconda3')]:
        if _base:
            _plugins = os.path.join(_base, 'lib', 'plugins')
            if os.path.isdir(_plugins):
                os.environ['OPENMM_PLUGIN_DIR'] = _plugins
                break

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

import numpy as np
import test_uma_ice_rpmd as _ice_mod
from openmm import app, unit, Vec3, Context, Platform, VerletIntegrator
from openmmml import MLPotential


# eV -> kJ/mol
EV_TO_KJ = 96.4853
# eV/Å -> kJ/(mol·nm): 1 eV/Å = 96.4853 / 0.1 = 964.853
EV_ANG_TO_KJ_NM = EV_TO_KJ * 10.0

MODEL_ASE = "uma-s-1p1"
MODEL_OMM = "uma-s-1p1-pythonforce-batch"


def get_ase_forces(positions_angstrom: np.ndarray, cell_angstrom: np.ndarray) -> tuple[float, np.ndarray]:
    """Single-point energy and forces via ASE + FAIRChem (reference)."""
    from ase import Atoms
    from ase.geometry import wrap_positions
    try:
        from fairchem.core import pretrained_mlip, FAIRChemCalculator
    except ImportError:
        from fairchem.core.calculate import pretrained_mlip
        from fairchem.core.calculate.ase_calculator import FAIRChemCalculator

    import torch
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    predictor = pretrained_mlip.get_predict_unit(MODEL_ASE, device=dev)
    calc = FAIRChemCalculator(predictor, task_name="omol")

    n_atoms = len(positions_angstrom)
    symbols = (["O", "H", "H"] * (n_atoms // 3))[:n_atoms]
    # Wrap positions into primary cell (matches UMA batch path and AtomicData.from_ase)
    wrapped = wrap_positions(positions_angstrom, cell_angstrom, pbc=[True, True, True])
    atoms = Atoms(symbols=symbols, positions=wrapped)
    atoms.set_cell(cell_angstrom)
    atoms.set_pbc(True)
    # Match OpenMM UMA PythonForce (charge=0, spin=1)
    atoms.info["charge"] = 0
    atoms.info["spin"] = 1
    atoms.calc = calc

    e_ev = atoms.get_potential_energy()
    f_ev_ang = atoms.get_forces()
    return float(e_ev), np.asarray(f_ev_ang)


def get_openmm_forces(topology: app.Topology, positions_nm: list, box_vectors: tuple) -> tuple[float, np.ndarray]:
    """Single-point energy and forces via OpenMM + UMA PythonForce."""
    platform = Platform.getPlatformByName("CUDA") if _has_cuda() else Platform.getPlatformByName("CPU")
    ml_dev = "cuda" if platform.getName() == "CUDA" else "cpu"

    potential = MLPotential(MODEL_OMM)
    system = potential.createSystem(topology, task_name="omol", charge=0, spin=1, device=ml_dev)

    integrator = VerletIntegrator(1.0 * unit.femtoseconds)
    context = Context(system, integrator, platform)
    context.setPositions([Vec3(p[0], p[1], p[2]) * unit.nanometer for p in positions_nm])
    context.setPeriodicBoxVectors(*box_vectors)

    state = context.getState(getEnergy=True, getForces=True)
    e_kj = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    f_kj_nm = state.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
    return e_kj, np.asarray(f_kj_nm)


def _has_cuda() -> bool:
    try:
        Platform.getPlatformByName("CUDA")
        return True
    except Exception:
        return False


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Compare OpenMM-UMA vs ASE-FAIRChem forces on ice")
    parser.add_argument("--input", "-i", type=Path, default=_SCRIPT_DIR / "ice.cif",
                        help="Input structure (CIF/PDB/XYZ)")
    parser.add_argument("--molecules", "-n", type=int, default=None,
                        help="Trim to N molecules (for fast debug)")
    parser.add_argument("--tol-energy", type=float, default=10.0,
                        help="Max allowed |ΔE| kJ/mol (default 10.0)")
    parser.add_argument("--tol-force", type=float, default=150.0,
                        help="Max allowed max|ΔF| kJ/(mol·nm) (default 150)")
    parser.add_argument("--output", "-o", type=Path, default=None,
                        help="Output PNG for force error vs particle index (default: force_error_vs_particle.png)")
    args = parser.parse_args()

    print("=" * 70)
    print("OpenMM-UMA vs ASE-FAIRChem single-point force comparison (ice)")
    print("=" * 70)
    print(f"Model: {MODEL_ASE} (ASE) vs {MODEL_OMM} (OpenMM)")
    print(f"Input: {args.input}")
    if args.molecules:
        print(f"Trimmed to {args.molecules} molecules")
    print()

    if not args.input.exists():
        print(f"Error: {args.input} not found")
        sys.exit(1)

    topology, positions_nm, box_vectors = _ice_mod._load_structure_from_file(
        str(args.input), max_molecules=args.molecules
    )
    n_atoms = topology.getNumAtoms()

    positions_nm_arr = np.array([[p[0], p[1], p[2]] for p in positions_nm])
    positions_angstrom = positions_nm_arr * 10.0

    box_matrix = np.array([[box_vectors[i][j].value_in_unit(unit.nanometer) for j in range(3)]
                           for i in range(3)])
    cell_angstrom = box_matrix * 10.0

    # Wrap positions once so both ASE and OpenMM see identical coordinates.
    # UMA batch path wraps; single-copy does not. Mismatch would cause force error.
    from ase.geometry import wrap_positions
    wrapped_angstrom = wrap_positions(positions_angstrom, cell_angstrom, pbc=[True, True, True])
    wrapped_nm = [(wrapped_angstrom[i, 0] * 0.1, wrapped_angstrom[i, 1] * 0.1, wrapped_angstrom[i, 2] * 0.1)
                  for i in range(n_atoms)]

    print("[1] ASE + FAIRChem (reference)...")
    e_ase_ev, f_ase_ev_ang = get_ase_forces(wrapped_angstrom, cell_angstrom)
    e_ase_kj = e_ase_ev * EV_TO_KJ
    f_ase_kj_nm = f_ase_ev_ang * EV_ANG_TO_KJ_NM

    print("[2] OpenMM + UMA...")
    e_omm_kj, f_omm_kj_nm = get_openmm_forces(topology, wrapped_nm, box_vectors)

    print()
    print("--- Results ---")
    print(f"  N_atoms: {n_atoms}")

    dE = e_omm_kj - e_ase_kj
    print(f"  Energy:")
    print(f"    ASE:   {e_ase_kj:.4f} kJ/mol")
    print(f"    OpenMM: {e_omm_kj:.4f} kJ/mol")
    print(f"    ΔE = {dE:+.4f} kJ/mol")

    df = f_omm_kj_nm - f_ase_kj_nm
    max_df = np.max(np.abs(df))
    rms_df = np.sqrt(np.mean(np.sum(df**2, axis=1)))
    print(f"\n  Forces:")
    print(f"    max |ΔF| = {max_df:.4f} kJ/(mol·nm)")
    print(f"    RMS |ΔF| = {rms_df:.4f} kJ/(mol·nm)")

    per_atom_max = np.max(np.abs(df), axis=1)
    worst_atom = np.argmax(per_atom_max)
    print(f"\n  Worst atom (index {worst_atom}):")
    print(f"    ΔF = [{df[worst_atom, 0]:+.4f}, {df[worst_atom, 1]:+.4f}, {df[worst_atom, 2]:+.4f}] kJ/(mol·nm)")
    print(f"    ASE F   = [{f_ase_kj_nm[worst_atom, 0]:.4f}, {f_ase_kj_nm[worst_atom, 1]:.4f}, {f_ase_kj_nm[worst_atom, 2]:.4f}]")
    print(f"    OMM F   = [{f_omm_kj_nm[worst_atom, 0]:.4f}, {f_omm_kj_nm[worst_atom, 1]:.4f}, {f_omm_kj_nm[worst_atom, 2]:.4f}]")

    # Plot force error vs particle index
    png_path = args.output or Path("force_error_vs_particle.png")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.bar(np.arange(n_atoms), per_atom_max, color="steelblue", edgecolor="navy", alpha=0.8)
        ax.set_xlabel("Particle index")
        ax.set_ylabel(r"Force error $|\Delta F|$ (kJ/(mol·nm))")
        ax.set_title(f"OpenMM vs ASE-FAIRChem force error per particle (max component)\n{args.input.name}")
        ax.axhline(y=max_df, color="red", linestyle="--", alpha=0.7, label=f"max = {max_df:.2f}")
        ax.legend()
        ax.set_xlim(-0.5, n_atoms - 0.5)
        plt.tight_layout()
        plt.savefig(png_path, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"\n  Plot saved: {png_path}")
    except ImportError:
        print(f"\n  (matplotlib not available, skipping PNG plot)")

    ok = abs(dE) <= args.tol_energy and max_df <= args.tol_force
    if ok:
        print(f"\n✓ PASS: Discrepancy within tolerance")
    else:
        print(f"\n✗ FAIL: Discrepancy exceeds tolerance (--tol-energy {args.tol_energy}, --tol-force {args.tol_force})")
        sys.exit(1)


if __name__ == "__main__":
    main()
