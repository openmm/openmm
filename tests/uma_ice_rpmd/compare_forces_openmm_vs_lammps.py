#!/usr/bin/env python3
"""
Compare single-point forces and energy: OpenMM-UMA vs LAMMPS-style UMA (same model, same structure).

Uses the same LAMMPS data file and replicates how each engine feeds positions to the model:
- OpenMM: molecular wrap (_wrap_water_molecules_per_bead) then batch path.
- LAMMPS: atom wrap (wrap_positions) then atomic_data_from_lammps_data (mirrors fix external).

This validates force parity between the two codes; discrepancies (wrap, units, order) can
explain different order-parameter evolution in the pipeline.

Usage:
  cd tests/uma_ice_rpmd
  python compare_forces_openmm_vs_lammps.py --data lammps/data.ice_uma
  python compare_forces_openmm_vs_lammps.py --data lammps/data.ice_uma --molecules 16
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

if os.getenv("OPENMM_PLUGIN_DIR") is None:
    for _base in [os.getenv("CONDA_PREFIX"), os.path.expanduser("~/miniconda3")]:
        if _base:
            _plugins = os.path.join(_base, "lib", "plugins")
            if os.path.isdir(_plugins):
                os.environ["OPENMM_PLUGIN_DIR"] = _plugins
                break

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

import numpy as np
from openmm import unit, Vec3, Context, Platform, VerletIntegrator
from openmm.app import Topology, element
from openmmml import MLPotential

from benchmark_same_structure_openmm_vs_lammps import parse_lammps_data
from run_openmm_ice_lammps_match import water_topology

# eV -> kJ/mol; eV/Å -> kJ/(mol·nm)
EV_TO_KJ = 96.4853
EV_ANG_TO_KJ_NM = EV_TO_KJ * 10.0

MODEL_OMM = "uma-s-1p1-pythonforce-batch"
MODEL_LMP = "uma-s-1p1"


def _has_cuda() -> bool:
    try:
        Platform.getPlatformByName("CUDA")
        return True
    except Exception:
        return False


def get_openmm_forces(
    topology: Topology,
    positions_nm: list[tuple[float, float, float]],
    box_vectors: tuple,
) -> tuple[float, np.ndarray]:
    """Single-point energy and forces via OpenMM + UMA (same as pipeline)."""
    platform = Platform.getPlatformByName("CUDA") if _has_cuda() else Platform.getPlatformByName("CPU")
    ml_dev = "cuda" if platform.getName() == "CUDA" else "cpu"
    potential = MLPotential(MODEL_OMM)
    system = potential.createSystem(
        topology,
        task_name="omol",
        charge=0,
        spin=1,
        device=ml_dev,
        use_atom_wrap_for_lammps_parity=True,
    )
    integrator = VerletIntegrator(1.0 * unit.femtoseconds)
    context = Context(system, integrator, platform)
    context.setPositions([Vec3(p[0], p[1], p[2]) * unit.nanometer for p in positions_nm])
    context.setPeriodicBoxVectors(*box_vectors)
    state = context.getState(getEnergy=True, getForces=True)
    e_kj = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    f_kj_nm = state.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
    return float(e_kj), np.asarray(f_kj_nm)


def get_lammps_style_forces(
    pos_ang: np.ndarray,
    cell: np.ndarray,
    atomic_numbers: np.ndarray,
    device: str = "cuda",
) -> tuple[float, np.ndarray]:
    """Single-point energy and forces using LAMMPS-style path (atom wrap + atomic_data_from_lammps_data)."""
    from ase.geometry import wrap_positions
    try:
        from fairchem.core import pretrained_mlip
    except ImportError:
        from fairchem.core.calculate import pretrained_mlip
    from fairchem.core.datasets.atomic_data import AtomicData
    from fairchem.lammps.lammps_fc import atomic_data_from_lammps_data, restricted_cell_from_lammps_box

    n = pos_ang.shape[0]
    periodicity = np.array([True, True, True])
    task_name = "omol"
    lx, ly, lz = cell[0, 0], cell[1, 1], cell[2, 2]
    xy, xz, yz = cell[1, 0], cell[2, 0], cell[2, 1]
    boxlo = np.array([0.0, 0.0, 0.0])
    boxhi = np.array([lx, ly, lz])
    cell_lmp = restricted_cell_from_lammps_box(boxlo, boxhi, xy, yz, xz)

    x_wrapped = wrap_positions(pos_ang, cell=cell, pbc=periodicity, eps=0)
    ad = atomic_data_from_lammps_data(
        x_wrapped,
        atomic_numbers,
        n,
        cell_lmp,
        periodicity,
        task_name,
        charge=0,
        spin=1,
    )
    predictor = pretrained_mlip.get_predict_unit(MODEL_LMP, device=device)
    pred = predictor.predict(ad)
    energy_ev = float(pred["energy"].detach().cpu().numpy().flat[0])
    forces_ev_ang = pred["forces"].detach().cpu().numpy()
    if forces_ev_ang.shape[0] != n:
        # Batched output: select first system
        batch = getattr(ad, "batch", None)
        if batch is not None:
            mask = batch.cpu().numpy() == 0
            forces_ev_ang = forces_ev_ang[mask]
        elif hasattr(pred["forces"], "batch"):
            mask = pred["forces"].batch.cpu().numpy() == 0
            forces_ev_ang = forces_ev_ang[mask]
    assert forces_ev_ang.shape[0] == n, f"forces shape {forces_ev_ang.shape} vs n={n}"
    energy_kj = energy_ev * EV_TO_KJ
    forces_kj_nm = forces_ev_ang * EV_ANG_TO_KJ_NM
    return energy_kj, forces_kj_nm


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(description="Compare OpenMM-UMA vs LAMMPS-style UMA forces on ice (same data file)")
    ap.add_argument("--data", type=Path, default=_SCRIPT_DIR / "lammps" / "data.ice_uma")
    ap.add_argument("--molecules", "-n", type=int, default=None, help="Use first N molecules only")
    ap.add_argument("--tol-energy", type=float, default=50.0, help="Max allowed |ΔE| kJ/mol")
    ap.add_argument("--tol-force", type=float, default=200.0, help="Max allowed max|ΔF| kJ/(mol·nm)")
    ap.add_argument("--device", default="cuda", choices=["cuda", "cpu"])
    args = ap.parse_args()

    if not args.data.is_file():
        print(f"Missing {args.data}. Run: python lammps/build_lammps_ice_data.py -n 128", file=sys.stderr)
        sys.exit(1)

    pos_ang, cell, Z = parse_lammps_data(args.data)
    n = Z.shape[0]
    if args.molecules is not None:
        n = min(n, args.molecules * 3)
        pos_ang = pos_ang[:n]
        Z = Z[:n]
    assert n % 3 == 0
    n_mol = n // 3

    lx, ly, lz = cell[0, 0], cell[1, 1], cell[2, 2]
    xy, xz, yz = cell[1, 0], cell[2, 0], cell[2, 1]
    a = Vec3(lx / 10, 0, 0) * unit.nanometer
    b = Vec3(xy / 10, ly / 10, 0) * unit.nanometer
    c = Vec3(xz / 10, yz / 10, lz / 10) * unit.nanometer
    box = (a, b, c)
    positions_nm = [tuple(float(x) for x in pos_ang[i] / 10.0) for i in range(n)]

    # OpenMM uses molecular wrap internally; we pass same coords as LAMMPS data (already consistent with data file)
    topology = water_topology(n_mol)
    print("OpenMM + UMA (pythonforce-batch)...")
    e_omm, f_omm = get_openmm_forces(topology, positions_nm, box)

    print("LAMMPS-style (atom wrap + atomic_data_from_lammps_data)...")
    e_lmp, f_lmp = get_lammps_style_forces(pos_ang, cell, Z, device=args.device)

    print()
    print("--- Results ---")
    print(f"  N_atoms: {n}")

    dE = e_omm - e_lmp
    print(f"  Energy:")
    print(f"    OpenMM: {e_omm:.4f} kJ/mol")
    print(f"    LAMMPS: {e_lmp:.4f} kJ/mol")
    print(f"    ΔE = {dE:+.4f} kJ/mol")

    df = f_omm - f_lmp
    max_df = np.max(np.abs(df))
    rms_df = np.sqrt(np.mean(np.sum(df**2, axis=1)))
    print(f"\n  Forces:")
    print(f"    max |ΔF| = {max_df:.4f} kJ/(mol·nm)")
    print(f"    RMS |ΔF| = {rms_df:.4f} kJ/(mol·nm)")

    per_atom_max = np.max(np.abs(df), axis=1)
    worst = int(np.argmax(per_atom_max))
    print(f"\n  Worst atom (index {worst}):")
    print(f"    ΔF = [{df[worst, 0]:+.4f}, {df[worst, 1]:+.4f}, {df[worst, 2]:+.4f}] kJ/(mol·nm)")
    print(f"    OpenMM F = [{f_omm[worst, 0]:.4f}, {f_omm[worst, 1]:.4f}, {f_omm[worst, 2]:.4f}]")
    print(f"    LAMMPS F = [{f_lmp[worst, 0]:.4f}, {f_lmp[worst, 1]:.4f}, {f_lmp[worst, 2]:.4f}]")

    ok = abs(dE) <= args.tol_energy and max_df <= args.tol_force
    if ok:
        print(f"\nPASS: Discrepancy within tolerance (--tol-energy {args.tol_energy}, --tol-force {args.tol_force})")
    else:
        print(f"\nFAIL: Discrepancy exceeds tolerance")
        sys.exit(1)


if __name__ == "__main__":
    main()
