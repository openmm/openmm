#!/usr/bin/env python3
"""
Ring-polymer (RPMD) check: ML energy/forces on each bead vs sequential Verlet reference.

For each bead copy ``b``, OpenMM reports forces from the same UMA potential evaluated at
that copy's positions (batched RPMD path). This script compares those to a reference
obtained by ``N`` independent Verlet ``Context.getState`` calls with the same positions.

The RPMD ring springs live in the integrator, not in the Context's Force list; both paths
use the full System force (typically UMA only), so energies/forces should match bead-wise.

Usage:
  cd tests/uma_ice_rpmd
  python compare_rpmd_bead_forces_openmm.py --beads 8
  python compare_rpmd_bead_forces_openmm.py --beads 4 --structure ice --data lammps/data.ice_uma --molecules 16

Defaults use identical bead positions (``--sigma-perturb-nm 0``) so each copy should match
sequential Verlet within tight tolerance. With random displacements between beads, the batched
RPMD path can differ from Verlet by a few kJ/(mol·nm) on small systems; use e.g.
``--tol-force 2.5 --tol-energy 2`` or match ``compare_forces_openmm_vs_lammps.py`` for ice
(``--tol-force 200 --tol-energy 50``).
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
_root = _SCRIPT_DIR.parent.parent
for _plugdir in (os.path.join(_root, "build"), os.path.join(_root, "build", "lib", "plugins")):
    if os.path.isdir(_plugdir) and any(
        f.startswith("libOpenMMRPMD") for f in os.listdir(_plugdir)
    ):
        os.environ["OPENMM_PLUGIN_DIR"] = _plugdir
        break

if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

import numpy as np
from openmm import Context, Platform, RPMDIntegrator, VerletIntegrator, Vec3, unit

try:
    import torch

    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
except ImportError:
    torch = None  # type: ignore[misc, assignment]

from openmmml import MLPotential
from openmm.app import Topology, element

from benchmark_same_structure_openmm_vs_lammps import parse_lammps_data
from run_openmm_ice_lammps_match import water_topology

MODEL = "uma-s-1p1-pythonforce-batch"


def _platform_ml_device(openmm_platform: str) -> tuple[Platform, str]:
    if openmm_platform == "cuda":
        return Platform.getPlatformByName("CUDA"), "cuda"
    if openmm_platform == "cpu":
        for name in ("CPU", "Reference"):
            try:
                return Platform.getPlatformByName(name), "cpu"
            except Exception:
                continue
        raise RuntimeError("No CPU/Reference platform")
    # auto
    try:
        return Platform.getPlatformByName("CUDA"), "cuda"
    except Exception:
        for name in ("CPU", "Reference"):
            try:
                return Platform.getPlatformByName(name), "cpu"
            except Exception:
                continue
    raise RuntimeError("No OpenMM platform available")


def _water_topology_positions_box(n_molecules: int = 1):
    from ase.build import molecule

    waters = [molecule("H2O") for _ in range(n_molecules)]
    atoms = waters[0]
    for w in waters[1:]:
        atoms += w
    cell_size = 12.0
    atoms.set_cell([cell_size, cell_size, cell_size])
    atoms.set_pbc(True)
    atoms.center()

    topology = Topology()
    for i in range(n_molecules):
        chain = topology.addChain()
        res = topology.addResidue("HOH", chain)
        topology.addAtom("O", element.oxygen, res)
        topology.addAtom("H", element.hydrogen, res)
        topology.addAtom("H", element.hydrogen, res)

    cell_nm = np.eye(3) * (cell_size * 0.1)
    box_vectors = (
        Vec3(cell_nm[0, 0], 0, 0) * unit.nanometer,
        Vec3(0, cell_nm[1, 1], 0) * unit.nanometer,
        Vec3(0, 0, cell_nm[2, 2]) * unit.nanometer,
    )
    topology.setPeriodicBoxVectors(box_vectors)
    pos_nm = (atoms.get_positions() * 0.1).tolist()
    return topology, pos_nm, box_vectors


def _ice_topology_positions_box(data_path: Path, n_molecules: int | None):
    pos_ang, cell, Z = parse_lammps_data(data_path)
    n = Z.shape[0]
    if n_molecules is not None:
        n = min(n, n_molecules * 3)
        pos_ang = pos_ang[:n]
    assert n % 3 == 0
    n_mol = n // 3
    lx, ly, lz = cell[0, 0], cell[1, 1], cell[2, 2]
    xy, xz, yz = cell[1, 0], cell[2, 0], cell[2, 1]
    a = Vec3(lx / 10, 0, 0) * unit.nanometer
    b = Vec3(xy / 10, ly / 10, 0) * unit.nanometer
    c = Vec3(xz / 10, yz / 10, lz / 10) * unit.nanometer
    box = (a, b, c)
    positions_nm = [tuple(float(x) for x in pos_ang[i] / 10.0) for i in range(n)]
    return water_topology(n_mol), positions_nm, box


def _vec_list(pos_nm: list[tuple[float, float, float]]):
    return [Vec3(float(p[0]), float(p[1]), float(p[2])) * unit.nanometer for p in pos_nm]


def _perturb_beads(
    base_nm: list[tuple[float, float, float]],
    n_beads: int,
    rng: np.random.Generator,
    sigma_nm: float,
) -> list[list[tuple[float, float, float]]]:
    base = np.asarray(base_nm, dtype=np.float64)
    out: list[list[tuple[float, float, float]]] = []
    for _ in range(n_beads):
        noise = rng.normal(scale=sigma_nm, size=base.shape)
        p = base + noise
        out.append([tuple(float(x) for x in row) for row in p])
    return out


def reference_bead_forces(
    system,
    platform: Platform,
    bead_positions: list[list[tuple[float, float, float]]],
    box_vectors,
) -> tuple[list[float], list[np.ndarray]]:
    integrator = VerletIntegrator(1.0 * unit.femtoseconds)
    ctx = Context(system, integrator, platform)
    ctx.setPeriodicBoxVectors(*box_vectors)
    energies: list[float] = []
    forces: list[np.ndarray] = []
    for pos in bead_positions:
        ctx.setPositions(_vec_list(pos))
        st = ctx.getState(getEnergy=True, getForces=True)
        energies.append(
            float(st.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole))
        )
        forces.append(
            np.asarray(
                st.getForces(asNumpy=True).value_in_unit(
                    unit.kilojoules_per_mole / unit.nanometer
                )
            )
        )
    return energies, forces


def rpmd_bead_forces(
    system,
    platform: Platform,
    bead_positions: list[list[tuple[float, float, float]]],
    box_vectors,
    temperature_K: float,
    friction: float,
) -> tuple[list[float], list[np.ndarray]]:
    n_beads = len(bead_positions)
    integrator = RPMDIntegrator(
        n_beads,
        temperature_K * unit.kelvin,
        friction / unit.picosecond,
        0.5 * unit.femtoseconds,
    )
    integrator.setThermostatType(RPMDIntegrator.Pile)
    ctx = Context(system, integrator, platform)
    ctx.setPeriodicBoxVectors(*box_vectors)
    ctx.setPositions(_vec_list(bead_positions[0]))
    for b in range(n_beads):
        integrator.setPositions(b, _vec_list(bead_positions[b]))

    energies: list[float] = []
    forces: list[np.ndarray] = []
    for b in range(n_beads):
        st = integrator.getState(
            b,
            getEnergy=True,
            getForces=True,
            enforcePeriodicBox=False,
        )
        energies.append(
            float(st.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole))
        )
        forces.append(
            np.asarray(
                st.getForces(asNumpy=True).value_in_unit(
                    unit.kilojoules_per_mole / unit.nanometer
                )
            )
        )
    return energies, forces


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(
        description="Compare per-bead UMA forces: RPMDIntegrator vs sequential Verlet"
    )
    ap.add_argument("--beads", type=int, default=8, help="Number of RPMD copies")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument(
        "--sigma-perturb-nm",
        type=float,
        default=0.0,
        help="RMS displacement per atom between beads (nm); 0 = all beads identical (strictest)",
    )
    ap.add_argument(
        "--structure",
        choices=("water", "ice"),
        default="water",
        help="water: single H2O in box; ice: LAMMPS data file",
    )
    ap.add_argument("--data", type=Path, default=_SCRIPT_DIR / "lammps" / "data.ice_uma")
    ap.add_argument(
        "--molecules",
        type=int,
        default=None,
        help="Ice only: use first N water molecules",
    )
    ap.add_argument(
        "--openmm-platform",
        default="auto",
        choices=["auto", "cuda", "cpu"],
    )
    ap.add_argument(
        "--tol-energy",
        type=float,
        default=0.05,
        help="Max |ΔE| per bead (kJ/mol); loosen if using --sigma-perturb-nm > 0",
    )
    ap.add_argument(
        "--tol-force",
        type=float,
        default=0.2,
        help="Max max|ΔF| per bead (kJ/(mol·nm)); try ~2.5 (water) or ~200 (ice) if perturbed",
    )
    args = ap.parse_args()

    if args.beads < 2:
        print("Use --beads >= 2 for a ring-polymer copy test (bead 0 alone is covered by test_force_layout).", file=sys.stderr)
        sys.exit(2)

    rng = np.random.default_rng(args.seed)
    if args.structure == "water":
        topology, base_nm, box = _water_topology_positions_box(1)
    else:
        if not args.data.is_file():
            print(f"Missing {args.data}", file=sys.stderr)
            sys.exit(1)
        topology, base_nm, box = _ice_topology_positions_box(args.data, args.molecules)

    bead_pos = _perturb_beads(base_nm, args.beads, rng, args.sigma_perturb_nm)

    platform, ml_dev = _platform_ml_device(args.openmm_platform)
    print(f"OpenMM platform: {platform.getName()}, UMA device: {ml_dev}")
    potential = MLPotential(MODEL)
    system = potential.createSystem(
        topology,
        task_name="omol",
        charge=0,
        spin=1,
        device=ml_dev,
        use_atom_wrap_for_lammps_parity=True,
    )

    print("Sequential Verlet (reference, one Context evaluation per bead)...")
    e_ref, f_ref = reference_bead_forces(system, platform, bead_pos, box)

    print("RPMDIntegrator (batched path, getState per bead)...")
    e_rpmd, f_rpmd = rpmd_bead_forces(
        system, platform, bead_pos, box, temperature_K=243.0, friction=1.0
    )

    n_atoms = len(base_nm)
    print()
    print("--- Per-bead ---")
    worst_bead = -1
    worst_f = 0.0
    all_ok = True
    for b in range(args.beads):
        dE = e_rpmd[b] - e_ref[b]
        df = f_rpmd[b] - f_ref[b]
        max_df = float(np.max(np.abs(df)))
        rms_df = float(np.sqrt(np.mean(np.sum(df**2, axis=1))))
        ok_b = abs(dE) <= args.tol_energy and max_df <= args.tol_force
        all_ok = all_ok and ok_b
        if max_df > worst_f:
            worst_f = max_df
            worst_bead = b
        status = "ok" if ok_b else "FAIL"
        print(
            f"  bead {b}: ΔE = {dE:+.6f} kJ/mol  max|ΔF| = {max_df:.6f}  RMS|ΔF| = {rms_df:.6f}  [{status}]"
        )

    print()
    print(f"N_atoms = {n_atoms}, N_beads = {args.beads}, structure = {args.structure}")
    if worst_bead >= 0:
        print(f"Worst bead by max|ΔF|: {worst_bead} (max |ΔF| = {worst_f:.6f} kJ/(mol·nm))")

    if all_ok:
        print(
            f"\nPASS: All beads within tolerance (|ΔE| <= {args.tol_energy}, max|ΔF| <= {args.tol_force})"
        )
    else:
        print(
            f"\nFAIL: At least one bead exceeds tolerance. "
            f"For large ice systems try e.g. --tol-force 200 --tol-energy 50",
            file=sys.stderr,
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
