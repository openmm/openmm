#!/usr/bin/env python3
"""
Phase 4: Validate GSD→OpenMM mapping by building the same system as run_simulation
when --initial-gsd and --no-cavity are used, then reporting potential energy.

If PE is reasonable (~-4e3 kJ/mol for 500 particles), the build path is correct and
the bug is elsewhere (e.g. how run_simulation uses the system after build).
If PE is huge (~3e7), the bug is in _build_system_from_gsd or force field/topology.

Usage (from examples/cavity/dimer_system/):
  python validate_gsd_openmm_energy.py [gsd_path] [frame]
  Default: init-0.gsd, frame 0
"""

from __future__ import annotations

import sys
from pathlib import Path

script_dir = Path(__file__).resolve().parent
if str(script_dir) not in sys.path:
    sys.path.insert(0, str(script_dir))

import run_simulation
from openmm import unit
from openmm import openmm


def main():
    gsd_path = sys.argv[1] if len(sys.argv) > 1 else "init-0.gsd"
    frame = int(sys.argv[2]) if len(sys.argv) > 2 else 0

    system, positions, topology, cavity_index, n_mol, box_nm = run_simulation._build_system_from_gsd(
        gsd_path, frame, gsd_in_nm=False, no_cavity=True, ff_dir=script_dir
    )
    n_particles = system.getNumParticles()
    assert n_particles == len(positions), f"positions length {len(positions)} != system particles {n_particles}"

    integrator = openmm.VerletIntegrator(0.001 * unit.picosecond)
    platform = openmm.Platform.getPlatformByName("CPU")
    context = openmm.Context(system, integrator, platform)
    context.setPositions(positions)
    context.setPeriodicBoxVectors(
        (box_nm, 0, 0) * unit.nanometer,
        (0, box_nm, 0) * unit.nanometer,
        (0, 0, box_nm) * unit.nanometer,
    )
    state = context.getState(getEnergy=True)
    pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    print("Phase 4: GSD→OpenMM energy validation")
    print(f"  GSD: {gsd_path}, frame {frame}")
    print(f"  Particles: {n_particles} (no cavity)")
    print(f"  Potential energy: {pe:.4f} kJ/mol")
    if abs(pe) > 1e6:
        print("  → HUGE PE: bug is in _build_system_from_gsd or force field/topology mapping.")
    elif -1e5 < pe < 0 or 0 < pe < 1e5:
        print("  → Reasonable PE: build path is correct; bug is elsewhere in run_simulation.")
    else:
        print("  → PE outside typical range; interpret manually.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
