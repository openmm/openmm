#!/usr/bin/env python3
"""
Dump per-particle forces from cav-hoomd for a GSD snapshot (reference for OpenMM comparison).

This script must be run in an environment where cav-hoomd (hoomd + hoomd.cavitymd) is
installed. It loads the GSD, builds a HOOMD simulation with the same dimer + cavity
setup, runs one step so forces are computed, and writes forces (N, 3) in Hartree/Bohr
to NPZ for use with compare_openmm_cav_hoomd_forces.py.

Expected NPZ format for compare_openmm_cav_hoomd_forces.py --hoomd-forces-npz:
  forces: (N, 3) array, dtype float64, in Hartree/Bohr (cav-hoomd native units).
  Particle order is GSD particle (tag) order so compare_openmm_cav_hoomd_forces.py
  can align HOOMD forces to OpenMM bond order via gsd_idx.

Usage (with cav-hoomd on PYTHONPATH):
  python dump_hoomd_forces_from_gsd.py init-0.gsd -o hoomd_forces.npz --lambda 0.001
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def get_hoomd_net_forces_in_tag_order(simulation) -> "np.ndarray":
    """Return net forces (N_global, 3) in sim force units, in GSD particle (tag) order.
    Uses state.cpu_local_snapshot.net_force + rtag so compare_openmm_cav_hoomd_forces.py
    can align via gsd_idx. F_global[tag] = net_force[rtag[tag]].
    """
    import numpy as np
    with simulation.state.cpu_local_snapshot as snap:
        net_force = np.array(snap.particles.net_force, copy=True)[:, :3]
        rtag = np.array(snap.particles.rtag, copy=True)
    return np.asarray(net_force)[np.asarray(rtag, dtype=np.intp)]


def main():
    import numpy as np

    parser = argparse.ArgumentParser(description="Dump cav-hoomd per-particle forces to NPZ (reference for OpenMM).")
    parser.add_argument("gsd", type=str, help="GSD file (dimer + cavity)")
    parser.add_argument("-o", "--output", type=str, default="hoomd_forces.npz", help="Output NPZ path")
    parser.add_argument("--lambda", type=float, default=0.001, dest="lambda_au", help="Coupling λ in a.u.")
    parser.add_argument("--frequency", type=float, default=1560.0, help="Cavity frequency cm^-1")
    parser.add_argument("--frame", type=int, default=0, help="GSD frame index")
    args = parser.parse_args()

    try:
        import hoomd
    except ImportError:
        print("HOOMD-blue is required. Install cav-hoomd or hoomd.", file=sys.stderr)
        sys.exit(1)
    try:
        from hoomd.cavitymd.simulation import CavityMDSimulation
    except ImportError:
        print("cav-hoomd (hoomd.cavitymd) is required. Run this script with cav-hoomd installed.", file=sys.stderr)
        sys.exit(1)

    gsd_path = Path(args.gsd).resolve()
    if not gsd_path.exists():
        print(f"GSD not found: {gsd_path}", file=sys.stderr)
        sys.exit(1)

    base_kwargs = dict(
        job_dir=str(gsd_path.parent),
        replica=0,
        freq=args.frequency,
        lambda_coupling=args.lambda_au,
        incavity=True,
        runtime_ps=0.002,
        input_gsd=str(gsd_path),
        frame=args.frame,
        name="force_dump",
        error_tolerance=0.0,
        temperature=100.0,
        molecular_thermostat="none",
        cavity_thermostat="none",
        finite_q=False,
        dt_fs=1.0,
        device="CPU",
        log_level="WARNING",
    )
    try:
        sim = CavityMDSimulation(**base_kwargs, enable_fkt=False)
    except TypeError:
        sim = CavityMDSimulation(**base_kwargs)
    # sim.sim is only created inside sim.run(); run minimally then read forces
    sim.run()
    hoomd_sim = getattr(sim, "sim", None) or getattr(sim, "sim_obj", None)
    if hoomd_sim is not None and hasattr(hoomd_sim, "sim"):
        hoomd_sim = hoomd_sim.sim
    if hoomd_sim is None:
        print("CavityMDSimulation does not expose the HOOMD simulation after run().", file=sys.stderr)
        sys.exit(1)
    F = get_hoomd_net_forces_in_tag_order(hoomd_sim)
    np.savez(args.output, forces=F)
    print(f"Wrote {args.output} (forces shape {F.shape}, Hartree/Bohr, GSD particle/tag order)")


if __name__ == "__main__":
    main()
