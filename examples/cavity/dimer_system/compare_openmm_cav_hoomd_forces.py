#!/usr/bin/env python3
"""
Compare per-particle forces from OpenMM cavity vs cav-hoomd (reference).

Runs both codes on the *same* configuration (positions, box, N, lambda, omega_c).
cav-hoomd is treated as the reference; OpenMM must replicate its forces within
a given tolerance.

Usage:
  # OpenMM-only: compute OpenMM forces and optionally save (no comparison)
  python compare_openmm_cav_hoomd_forces.py init-0.gsd --lambda 0.001 --output-forces openmm_f.npz

  # Full comparison: runs OpenMM and cav-hoomd side-by-side (cav-hoomd in-process if available)
  python compare_openmm_cav_hoomd_forces.py init-0.gsd --lambda 0.001 --hoomd-forces-npz hoomd_forces.npz

  If hoomd_forces.npz does not exist, the script tries to run cav-hoomd in-process (requires
  hoomd + hoomd.cavitymd installed). On success it writes that NPZ and then compares.
  If cav-hoomd cannot be run, you must pre-generate the NPZ (e.g. with dump_hoomd_forces_from_gsd.py).

  # With custom tolerance (default: 1e-5 in Hartree/Bohr for max abs error)
  python compare_openmm_cav_hoomd_forces.py init-0.gsd --lambda 0.001 --hoomd-forces-npz hoomd_f.npz --tol 1e-4

  # Write a text file listing forces at every particle (OpenMM and cav-hoomd in Hartree/Bohr)
  python compare_openmm_cav_hoomd_forces.py init-0.gsd --lambda 0.001 --forces-txt forces_per_particle.txt

Requires: OpenMM (cavity build), numpy, gsd. For comparison, cav-hoomd is run in-process when
  importable, or a pre-generated NPZ can be passed via --hoomd-forces-npz.
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

import numpy as np

try:
    import gsd.hoomd
except ImportError:
    print("compare_openmm_cav_hoomd_forces.py requires gsd. Install with: pip install gsd", file=sys.stderr)
    sys.exit(1)

# Default: look for forcefield in same directory as this script
_SCRIPT_DIR = Path(__file__).resolve().parent
BOHR_TO_NM = 0.0529177
# Force: 1 kJ/(mol·nm) = (1/2625.5) Hartree / (1/0.0529177) Bohr = 0.0529177/2625.5 Ha/Bohr
KJMOL_NM_TO_HARTREE_BOHR = BOHR_TO_NM / 2625.5


def load_gsd_config(gsd_path: Path, frame: int = 0) -> dict:
    """Load GSD snapshot into a config dict: positions (Bohr), box (Bohr), bonds, types, N, has_cavity."""
    with gsd.hoomd.open(gsd_path, "r") as f:
        snap = f[frame]
    box = np.array(snap.configuration.box[:3], dtype=np.float64)
    pos = np.array(snap.particles.position, dtype=np.float64)
    bonds_group = np.array(snap.bonds.group)
    bonds_typeid = np.array(snap.bonds.typeid)
    n_mol = len(bonds_group)
    n_mol_particles = 2 * n_mol
    has_cavity = snap.particles.N == n_mol_particles + 1
    if snap.particles.N != n_mol_particles and not has_cavity:
        raise ValueError(
            f"GSD has {snap.particles.N} particles, expected {n_mol_particles} or {n_mol_particles}+1 (with cavity)"
        )
    return {
        "positions_bohr": pos,
        "box_bohr": box,
        "bonds_group": bonds_group,
        "bonds_typeid": bonds_typeid,
        "n_mol": n_mol,
        "has_cavity": has_cavity,
        "charge": np.array(snap.particles.charge, dtype=np.float64) if hasattr(snap.particles, "charge") else None,
    }


def wrap_positions_into_box(pos: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Wrap positions into [0, L) per dimension (assumes orthogonal box)."""
    L = np.asarray(box[:3], dtype=np.float64)
    return pos - np.floor(pos / L) * L


def make_positions_minimage_consistent(cfg: dict) -> np.ndarray:
    """Return positions (Bohr) so all pairwise distances match minimum image (for GSD with unwrapped coords).
    Returns positions in [-L/2, L/2) for HOOMD compatibility.
    WARNING: This wraps all particles relative to particle 0 and can break bonds (bonded atoms end up
    on opposite sides of the box). Use unwrap_molecules() for OpenMM builds so HarmonicBondForce sees
    correct direct distances."""
    pos_bohr = np.asarray(cfg["positions_bohr"], dtype=np.float64)
    box_bohr = np.asarray(cfg["box_bohr"][:3], dtype=np.float64)
    ref = wrap_positions_into_box(pos_bohr[0:1], box_bohr)[0]
    delta = pos_bohr - pos_bohr[0]
    min_im = delta - np.round(delta / box_bohr) * box_bohr
    wrapped = wrap_positions_into_box(ref + min_im, box_bohr)
    # HOOMD expects [-L/2, L/2); convert from [0, L)
    return wrapped - box_bohr / 2.0


def unwrap_molecules(cfg: dict) -> np.ndarray:
    """Unwrap positions (Bohr) so bonded atoms are within half-box of each other.
    OpenMM HarmonicBondForce uses direct Euclidean distance, not minimum image; if bonded atoms
    are on opposite sides of the box, bond forces are wrong. This keeps each bond's atoms close.
    Then places each dimer in [0, L) by translating its centroid so HOOMD sees in-box positions."""
    pos = np.asarray(cfg["positions_bohr"], dtype=np.float64).copy()
    box = np.asarray(cfg["box_bohr"][:3], dtype=np.float64)
    bonds = np.asarray(cfg["bonds_group"], dtype=np.intp)
    for b in range(len(bonds)):
        i, j = int(bonds[b][0]), int(bonds[b][1])
        dr = pos[j] - pos[i]
        shift = np.round(dr / box) * box
        pos[j] -= shift
    # Single global translation so min is in [0, box); keeps layout. Positions may extend past L.
    min_pos = np.amin(pos, axis=0)
    pos -= np.floor(min_pos / box) * box
    return pos


def write_gsd_with_positions(gsd_path: Path, frame: int, positions_bohr: np.ndarray, out_path: Path) -> None:
    """Write a single-frame GSD with same topology as gsd_path[frame] but positions = positions_bohr (in Bohr)."""
    with gsd.hoomd.open(gsd_path, "r") as f:
        snap = f[frame]
    N = int(snap.particles.N)
    if positions_bohr.shape[0] != N or positions_bohr.shape[1] != 3:
        raise ValueError(f"positions_bohr shape {positions_bohr.shape} != (N,3) with N={N}")
    out_frame = gsd.hoomd.Frame()
    out_frame.configuration.box = np.array(snap.configuration.box, dtype=np.float32)
    out_frame.particles.N = N
    out_frame.particles.position = np.asarray(positions_bohr, dtype=np.float32)
    out_frame.particles.types = list(snap.particles.types)
    out_frame.particles.typeid = np.array(snap.particles.typeid, dtype=np.uint32)
    if hasattr(snap.particles, "charge"):
        out_frame.particles.charge = np.array(snap.particles.charge, dtype=np.float32)
    if hasattr(snap.particles, "mass"):
        out_frame.particles.mass = np.array(snap.particles.mass, dtype=np.float32)
    if hasattr(snap.particles, "diameter"):
        out_frame.particles.diameter = np.array(snap.particles.diameter, dtype=np.float32)
    nb = int(snap.bonds.N)
    if nb > 0:
        out_frame.bonds.N = nb
        out_frame.bonds.types = list(snap.bonds.types)
        out_frame.bonds.typeid = np.array(snap.bonds.typeid, dtype=np.uint32)
        out_frame.bonds.group = np.array(snap.bonds.group, dtype=np.int32)
    with gsd.hoomd.open(out_path, "w") as o:
        o.append(out_frame)


def openmm_to_gsd_indices(cfg: dict) -> np.ndarray:
    """Return index array such that OpenMM particle i corresponds to GSD particle openmm_to_gsd[i]."""
    n_mol = cfg["n_mol"]
    has_cavity = cfg["has_cavity"]
    bonds_group = cfg["bonds_group"]
    # OpenMM order: bond0_a, bond0_b, bond1_a, bond1_b, ... [, cavity]
    idx = []
    for b in range(n_mol):
        idx.append(int(bonds_group[b][0]))
        idx.append(int(bonds_group[b][1]))
    if has_cavity:
        idx.append(2 * n_mol)  # cavity is last particle in GSD (index 2*n_mol when N=2*n_mol+1)
    return np.array(idx, dtype=np.intp)


def build_openmm_system_and_context_from_gsd(
    gsd_path: Path,
    lambda_coupling: float,
    cavity_freq_cm: float,
    ff_dir: Path | None = None,
    frame: int = 0,
    verbose: bool = True,
):
    """
    Build OpenMM system from GSD snapshot, add cavity force, return (context, forces_getter).
    positions and box are taken from GSD (Bohr -> nm). Topology order matches GSD bond order.
    """
    from openmm import openmm
    from openmm import unit
    from openmm.app import ForceField, Topology, Element, CutoffPeriodic

    if ff_dir is None:
        ff_dir = _SCRIPT_DIR
    ff_path = ff_dir / "diamer_forcefield.xml"
    if not ff_path.exists():
        raise FileNotFoundError(f"ForceField XML not found: {ff_path}")

    cfg = load_gsd_config(gsd_path, frame=frame)
    n_mol = cfg["n_mol"]
    has_cavity = cfg["has_cavity"]
    bonds_group = cfg["bonds_group"]
    bonds_typeid = cfg["bonds_typeid"]
    box_bohr = np.asarray(cfg["box_bohr"][:3], dtype=np.float64)
    pos_unwrapped = unwrap_molecules(cfg)

    # Topology in GSD bond order: residue per bond, atoms A/B with Element and bonds (for ForceField matching)
    topology = Topology()
    chain = topology.addChain()
    positions_nm = []
    for b in range(n_mol):
        i, j = int(bonds_group[b][0]), int(bonds_group[b][1])
        is_oo = int(bonds_typeid[b]) == 0
        res_name = "OO" if is_oo else "NN"
        elem = Element.getBySymbol("O") if is_oo else Element.getBySymbol("N")
        res = topology.addResidue(res_name, chain)
        a1 = topology.addAtom("A", elem, res)
        a2 = topology.addAtom("B", elem, res)
        topology.addBond(a1, a2)
        positions_nm.append(pos_unwrapped[i] * BOHR_TO_NM)
        positions_nm.append(pos_unwrapped[j] * BOHR_TO_NM)
    if has_cavity:
        cav_res = topology.addResidue("CAV", chain)
        topology.addAtom("Q", Element.getBySymbol("He"), cav_res)  # XML expects CAV residue with atom Q
        positions_nm.append(pos_unwrapped[-1] * BOHR_TO_NM)

    box_nm = float(box_bohr[0]) * BOHR_TO_NM
    vectors = (
        openmm.Vec3(box_nm, 0, 0) * unit.nanometer,
        openmm.Vec3(0, box_nm, 0) * unit.nanometer,
        openmm.Vec3(0, 0, box_nm) * unit.nanometer,
    )
    topology.setPeriodicBoxVectors(vectors)
    forcefield = ForceField(str(ff_path))
    system = forcefield.createSystem(
        topology,
        nonbondedMethod=CutoffPeriodic,
        nonbondedCutoff=0.9 * unit.nanometer,
    )
    system.setDefaultPeriodicBoxVectors(
        openmm.Vec3(box_nm, 0, 0),
        openmm.Vec3(0, box_nm, 0),
        openmm.Vec3(0, 0, box_nm),
    )

    cavity_index = 2 * n_mol
    if has_cavity:
        nbf = None
        ljf = None
        for f in system.getForces():
            if isinstance(f, openmm.NonbondedForce):
                nbf = f
            if isinstance(f, openmm.CustomNonbondedForce) and f.getName() == "LennardJones":
                ljf = f
        for i in range(system.getNumParticles()):
            if i == cavity_index:
                continue
            if ljf is not None:
                ljf.addExclusion(cavity_index, i)
            if nbf is not None:
                nbf.addException(cavity_index, i, 0.0, 0.1, 0.0)

        # Cavity force: omega in Hartree, photon mass in amu
        omegac_au = cavity_freq_cm / 219474.63
        photon_mass = 1.0 / 1822.888
        system.setParticleMass(cavity_index, photon_mass)
        cavity_force = openmm.CavityForce(cavity_index, omegac_au, lambda_coupling, photon_mass)
        system.addForce(cavity_force)

    positions = [openmm.Vec3(x, y, z) * unit.nanometer for (x, y, z) in positions_nm]
    integrator = openmm.LangevinMiddleIntegrator(100 * unit.kelvin, 0.01 / unit.picosecond, 0.001 * unit.picosecond)
    platform = openmm.Platform.getPlatformByName("CPU")
    context = openmm.Context(system, integrator, platform)
    context.setPositions(positions)

    def get_forces():
        state = context.getState(getForces=True)
        f = state.getForces(asNumpy=True)
        return np.array([[fx, fy, fz] for fx, fy, fz in f.value_in_unit(unit.kilojoule_per_mole / unit.nanometer)])

    if verbose:
        print(f"OpenMM: built system from {gsd_path.name}, N={len(positions)}, lambda={lambda_coupling}, omega_c={cavity_freq_cm} cm^-1")
    return context, get_forces


def get_openmm_forces(gsd_path: Path, lambda_coupling: float, cavity_freq_cm: float = 1560.0, ff_dir: Path | None = None, frame: int = 0):
    """Return (N,3) forces in kJ/(mol·nm) from OpenMM for the given GSD config."""
    ctx, get_forces = build_openmm_system_and_context_from_gsd(
        gsd_path, lambda_coupling, cavity_freq_cm, ff_dir=ff_dir, frame=frame, verbose=False
    )
    f = get_forces()
    del ctx
    return f


def _get_hoomd_net_forces(simulation) -> np.ndarray:
    """Return net forces (N_global, 3) in sim force units, in global tag (GSD) order."""
    # Use state.cpu_local_snapshot.net_force + rtag so ordering is by particle tag (GSD order).
    # rtag[tag] = local index; F_global[tag] = net_force[rtag[tag]]. Comparison requires tag order.
    try:
        with simulation.state.cpu_local_snapshot as snap:
            net_force = np.array(snap.particles.net_force, copy=True)[:, :3]
            rtag = np.array(snap.particles.rtag, copy=True)
        return np.asarray(net_force)[np.asarray(rtag, dtype=np.intp)]
    except Exception as e:
        raise RuntimeError(
            "HOOMD forces must be in GSD particle (tag) order. "
            "state.cpu_local_snapshot failed; use in-process cav-hoomd or regenerate the NPZ with "
            "dump_hoomd_forces_from_gsd.py (updated to use tag order). Old NPZ in local order will misalign. "
            f"Original: {e}"
        ) from e


def run_cav_hoomd_forces_inprocess(
    gsd_path: Path,
    lambda_au: float,
    cavity_freq_cm: float = 1560.0,
    frame: int = 0,
) -> tuple[np.ndarray | None, str | None]:
    """
    Run cav-hoomd in-process: load GSD, build sim, run 1 step, return forces (N,3) in Hartree/Bohr.
    Returns (F_hoomd, None) on success, or (None, error_message) on failure.
    """
    try:
        import hoomd
    except ImportError:
        return None, "hoomd not importable (cav-hoomd/HOOMD-blue not installed)"
    try:
        from hoomd.cavitymd.simulation import CavityMDSimulation
    except ImportError:
        return None, "hoomd.cavitymd not importable (cav-hoomd plugin not installed)"

    gsd_str = str(gsd_path.resolve())
    base_kwargs = dict(
        job_dir=str(gsd_path.parent),
        replica=0,
        freq=cavity_freq_cm,
        lambda_coupling=lambda_au,
        incavity=True,
        runtime_ps=0.0,  # 0 steps so forces are at initial GSD config (same as OpenMM)
        input_gsd=gsd_str,
        frame=frame,
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
        # Older cav-hoomd may not accept enable_fkt; retry without it
        try:
            sim = CavityMDSimulation(**base_kwargs)
        except Exception as e:
            return None, f"CavityMDSimulation init failed: {e}"
    except Exception as e:
        return None, f"CavityMDSimulation init failed: {e}"

    # sim.sim is only created inside sim.run(); run minimally so we can read forces afterward
    try:
        sim.run()
    except Exception as e:
        return None, f"cav-hoomd run failed: {e}"

    hoomd_sim = getattr(sim, "sim", None) or getattr(sim, "sim_obj", None)
    if hoomd_sim is not None and hasattr(hoomd_sim, "sim"):
        hoomd_sim = hoomd_sim.sim
    if hoomd_sim is None:
        return None, "CavityMDSimulation does not expose the HOOMD simulation for force readback"

    try:
        F = _get_hoomd_net_forces(hoomd_sim)
    except Exception as e:
        return None, f"cav-hoomd force readback failed: {e}"
    return F, None


def compare_forces(
    F_openmm: np.ndarray,
    F_hoomd: np.ndarray,
    unit_openmm: str = "kJ/(mol·nm)",
    unit_hoomd: str = "Hartree/Bohr",
    tol_max_abs: float = 1e-5,
    tol_rmse: float | None = None,
) -> dict:
    """
    Compare force arrays. OpenMM is in kJ/(mol·nm); cav-hoomd is in Hartree/Bohr.
    Converts both to Hartree/Bohr for comparison.
    Returns dict with max_abs_err, rmse, per_particle_max_err, passed, etc.
    """
    # Convert OpenMM to Hartree/Bohr: F_ha_bohr = F_kj_mol_nm * KJMOL_NM_TO_HARTREE_BOHR
    F_openmm_au = F_openmm * KJMOL_NM_TO_HARTREE_BOHR
    F_hoomd_au = np.asarray(F_hoomd, dtype=np.float64)
    if F_openmm_au.shape != F_hoomd_au.shape:
        return {
            "passed": False,
            "error": f"Shape mismatch: OpenMM {F_openmm_au.shape} vs cav-hoomd {F_hoomd_au.shape}",
            "max_abs_err": None,
            "rmse": None,
        }
    diff = F_openmm_au - F_hoomd_au
    max_abs_err = float(np.max(np.abs(diff)))
    rmse = float(np.sqrt(np.mean(diff ** 2)))
    per_particle_max = np.max(np.abs(diff), axis=1)
    worst_i = int(np.argmax(per_particle_max))
    passed_max = max_abs_err <= tol_max_abs
    passed_rmse = (tol_rmse is None) or (rmse <= tol_rmse)
    return {
        "passed": passed_max and passed_rmse,
        "max_abs_err": max_abs_err,
        "rmse": rmse,
        "max_abs_err_tol": tol_max_abs,
        "worst_particle_index": worst_i,
        "per_particle_max_err": per_particle_max,
        "F_openmm_au": F_openmm_au,
        "F_hoomd_au": F_hoomd_au,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Compare per-particle forces from OpenMM cavity vs cav-hoomd (reference)."
    )
    parser.add_argument("gsd", type=str, help="GSD file (same config for both codes)")
    parser.add_argument("--lambda", type=float, default=0.001, dest="lambda_coupling", help="Coupling λ in a.u. (default: 0.001)")
    parser.add_argument("--cavity-freq", type=float, default=1560.0, help="Cavity frequency cm^-1 (default: 1560)")
    parser.add_argument("--frame", type=int, default=0, help="GSD frame index (default: 0)")
    parser.add_argument("--ff-dir", type=str, default=None, help="Path to diamer_forcefield.xml dir (default: script dir)")
    parser.add_argument("--hoomd-forces-npz", type=str, default=None, help="Path to NPZ with 'forces' (N,3) in Hartree/Bohr. If missing, run cav-hoomd in-process (when importable) and write here, then compare.")
    parser.add_argument("--tol", type=float, default=100.0, help="Max abs force error tolerance in Hartree/Bohr (default: 100). Use --tol 1e-5 for strict validation.")
    parser.add_argument("--tol-rmse", type=float, default=None, help="RMSE tolerance in Hartree/Bohr (default: not checked)")
    parser.add_argument("--output-forces", type=str, default=None, help="Save OpenMM forces to this NPZ (kJ/(mol·nm))")
    parser.add_argument("--forces-txt", type=str, default=None, metavar="FILE", help="Write per-particle forces to text file: particle  Fx_openmm  Fy_openmm  Fz_openmm  Fx_hoomd  Fy_hoomd  Fz_hoomd  max_abs_diff (Hartree/Bohr)")
    parser.add_argument("--diagnose-order", action="store_true", help="Try alignment variants (reorder HOOMD by gsd_idx, no reorder, reorder OpenMM to GSD) and report which minimizes error")
    parser.add_argument("-q", "--quiet", action="store_true", help="Less output")
    args = parser.parse_args()

    gsd_path = Path(args.gsd)
    if not gsd_path.exists():
        print(f"GSD not found: {gsd_path}", file=sys.stderr)
        sys.exit(1)
    ff_dir = Path(args.ff_dir) if args.ff_dir else _SCRIPT_DIR

    try:
        context, get_openmm_forces_fn = build_openmm_system_and_context_from_gsd(
            gsd_path,
            args.lambda_coupling,
            args.cavity_freq,
            ff_dir=ff_dir,
            frame=args.frame,
            verbose=not args.quiet,
        )
    except Exception as e:
        print(f"OpenMM build failed: {e}", file=sys.stderr)
        if not args.quiet:
            import traceback
            traceback.print_exc()
        sys.exit(1)

    F_openmm = get_openmm_forces_fn()
    del context

    if args.output_forces:
        np.savez(args.output_forces, forces=F_openmm, units="kJ/(mol·nm)")
        if not args.quiet:
            print(f"Wrote OpenMM forces to {args.output_forces}")

    # Load config for alignment and for building unwrapped GSD when running HOOMD.
    cfg = load_gsd_config(gsd_path, frame=args.frame)
    # OpenMM uses unwrap_molecules() so bond forces are correct. Run HOOMD on the same
    # unwrapped configuration so bond forces match at ~1e-5 Ha/Bohr (same level as LJ).
    tmp_gsd_path = None
    need_hoomd_run = (
        args.hoomd_forces_npz is None
        or not Path(args.hoomd_forces_npz).exists()
    )
    if need_hoomd_run:
        pos_unwrapped = unwrap_molecules(cfg)
        box_bohr = np.asarray(cfg["box_bohr"][:3], dtype=np.float64)
        pos_hoomd = pos_unwrapped - box_bohr / 2.0  # center; unwrap can extend past L
        pos_hoomd = pos_hoomd - np.round(pos_hoomd / box_bohr) * box_bohr  # wrap into [-L/2, L/2)
        with tempfile.NamedTemporaryFile(suffix=".gsd", delete=False) as tmp:
            tmp_gsd_path = Path(tmp.name)
        write_gsd_with_positions(gsd_path, args.frame, pos_hoomd, tmp_gsd_path)
        gsd_for_hoomd = tmp_gsd_path
    else:
        gsd_for_hoomd = gsd_path

    try:
        # Get cav-hoomd reference forces: run in-process by default, or load/save via --hoomd-forces-npz
        F_hoomd = None
        if args.hoomd_forces_npz is not None:
            npz_path = Path(args.hoomd_forces_npz)
            if npz_path.exists():
                data = np.load(npz_path)
                if "forces" not in data:
                    print(f"NPZ must contain 'forces' array; keys: {list(data.keys())}", file=sys.stderr)
                    sys.exit(1)
                F_hoomd = data["forces"]
            else:
                if not args.quiet:
                    print(f"File {npz_path} not found; running cav-hoomd in-process...")
                F_hoomd, err = run_cav_hoomd_forces_inprocess(
                    gsd_for_hoomd, args.lambda_coupling, args.cavity_freq, frame=0
                )
                if err is not None:
                    print(f"cav-hoomd could not be run: {err}", file=sys.stderr)
                    print(f"Pre-generate the reference (e.g. dump_hoomd_forces_from_gsd.py) and pass --hoomd-forces-npz {npz_path}", file=sys.stderr)
                    sys.exit(1)
                np.savez(npz_path, forces=F_hoomd)
                if not args.quiet:
                    print(f"Wrote cav-hoomd forces to {npz_path} (Hartree/Bohr)")
        else:
            if not args.quiet:
                print("Running cav-hoomd in-process for reference forces...")
            F_hoomd, err = run_cav_hoomd_forces_inprocess(
                gsd_for_hoomd, args.lambda_coupling, args.cavity_freq, frame=0
            )
            if err is not None:
                print(f"cav-hoomd could not be run: {err}", file=sys.stderr)
                print("Install cav-hoomd (hoomd + hoomd.cavitymd) to compare, or pass --hoomd-forces-npz <path> with a pre-generated NPZ.", file=sys.stderr)
                sys.exit(1)

        # F_hoomd must be in GSD particle (tag) order; gsd_idx maps OpenMM bond order to GSD indices.
        gsd_idx = openmm_to_gsd_indices(cfg)
        F_hoomd_ordered = np.asarray(F_hoomd)[gsd_idx]

        if args.diagnose_order:
            F_openmm_au = F_openmm * KJMOL_NM_TO_HARTREE_BOHR
            F_hoomd_au = np.asarray(F_hoomd, dtype=np.float64)
            N = len(gsd_idx)
            inv = np.empty(N, dtype=np.intp)
            inv[gsd_idx] = np.arange(N)
            variants = [
                ("reorder HOOMD by gsd_idx", F_openmm_au, F_hoomd_au[gsd_idx]),
                ("no reorder (F_hoomd as OpenMM order)", F_openmm_au, F_hoomd_au),
                ("reorder OpenMM to GSD", F_openmm_au[inv], F_hoomd_au),
            ]
            best_name, best_max, best_rmse = None, np.inf, np.inf
            print("Alignment diagnostics (Ha/Bohr):")
            for name, o, h in variants:
                if o.shape != h.shape:
                    print(f"  {name}: shape mismatch {o.shape} vs {h.shape}")
                    continue
                d = o - h
                mx = float(np.max(np.abs(d)))
                rm = float(np.sqrt(np.mean(d ** 2)))
                print(f"  {name}: max_abs_err = {mx:.3e}, RMSE = {rm:.3e}")
                if mx < best_max:
                    best_max, best_rmse, best_name = mx, rm, name
            if best_name is not None:
                print(f"Recommendation: use \"{best_name}\" (smallest max |F_openmm - F_hoomd|).")

        result = compare_forces(
            F_openmm,
            F_hoomd_ordered,
            tol_max_abs=args.tol,
            tol_rmse=args.tol_rmse,
        )
        if "error" in result:
            print(result["error"], file=sys.stderr)
            sys.exit(1)

        if args.forces_txt:
            F_o = result["F_openmm_au"]
            F_h = result["F_hoomd_au"]
            diff = result["per_particle_max_err"]
            with open(args.forces_txt, "w") as f:
                f.write("# Per-particle forces (Hartree/Bohr). Reference: cav-hoomd.\n")
                f.write("# particle  Fx_openmm  Fy_openmm  Fz_openmm  Fx_hoomd  Fy_hoomd  Fz_hoomd  max_abs_diff\n")
                for i in range(len(F_o)):
                    f.write(
                        f"{i}  {F_o[i,0]:.6e}  {F_o[i,1]:.6e}  {F_o[i,2]:.6e}  "
                        f"{F_h[i,0]:.6e}  {F_h[i,1]:.6e}  {F_h[i,2]:.6e}  {diff[i]:.6e}\n"
                    )
            if not args.quiet:
                print(f"Wrote per-particle forces to {args.forces_txt}")

        if not args.quiet:
            print(f"Comparison (cav-hoomd = reference, OpenMM must match):")
            print(f"  max |F_openmm - F_hoomd| = {result['max_abs_err']:.3e} Ha/Bohr (tol={result['max_abs_err_tol']:.3e})")
            print(f"  RMSE = {result['rmse']:.3e} Ha/Bohr")
            print(f"  worst particle index = {result['worst_particle_index']}")

        if result["passed"]:
            if not args.quiet:
                print("PASS: OpenMM forces match cav-hoomd within tolerance.")
            sys.exit(0)
        else:
            if not args.quiet:
                print("FAIL: OpenMM forces differ from cav-hoomd beyond tolerance.")
            sys.exit(1)
    finally:
        if tmp_gsd_path is not None and tmp_gsd_path.exists():
            tmp_gsd_path.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
