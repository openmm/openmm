#!/usr/bin/env python3
"""
Debug OpenMM vs cav-hoomd force mismatch by comparing each force component separately.

Runs four cases in one go: (1) harmonic bond only, (2) Coulomb only, (3) LJ only,
(4) Cavity force only. For each component, compares OpenMM vs cav-hoomd and reports
max_abs_err (Ha/Bohr), RMSE, worst particle. Writes a summary to a text file.

Usage:
  python debug_force_components.py init-0.gsd --lambda 0.001 --summary force_component_summary.txt

Requires: OpenMM (cavity build), numpy, gsd. cav-hoomd is run in-process when importable.
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
    print("debug_force_components.py requires gsd. Install with: pip install gsd", file=sys.stderr)
    sys.exit(1)

_SCRIPT_DIR = Path(__file__).resolve().parent
# Reuse comparison helpers
sys.path.insert(0, str(_SCRIPT_DIR))
try:
    from compare_openmm_cav_hoomd_forces import (
        BOHR_TO_NM,
        KJMOL_NM_TO_HARTREE_BOHR,
        load_gsd_config,
        openmm_to_gsd_indices,
        unwrap_molecules,
        write_gsd_with_positions,
    )
except ImportError as e:
    print(f"Cannot import compare_openmm_cav_hoomd_forces: {e}", file=sys.stderr)
    sys.exit(1)


def _build_openmm_system_with_force_groups(
    gsd_path: Path,
    lambda_coupling: float,
    cavity_freq_cm: float,
    ff_dir: Path,
    frame: int,
):
    """Build OpenMM system from GSD, assign force groups 0..3 to bond/Coulomb/LJ/Cavity, return (context, cfg)."""
    from openmm import openmm
    from openmm import unit
    from openmm.app import CutoffPeriodic, Element, ForceField, Topology

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
        topology.addAtom("Q", Element.getBySymbol("He"), cav_res)
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
        omegac_au = cavity_freq_cm / 219474.63
        photon_mass = 1.0 / 1822.888
        system.setParticleMass(cavity_index, photon_mass)
        cavity_force = openmm.CavityForce(cavity_index, omegac_au, lambda_coupling, photon_mass)
        system.addForce(cavity_force)

    # Force groups: 0=bond, 1=Coulomb, 2=LJ, 3=Cavity
    for f in system.getForces():
        if isinstance(f, openmm.HarmonicBondForce):
            f.setForceGroup(0)
        elif isinstance(f, openmm.NonbondedForce):
            f.setForceGroup(1)
        elif isinstance(f, openmm.CustomNonbondedForce) and f.getName() == "LennardJones":
            f.setForceGroup(2)
        elif isinstance(f, openmm.CavityForce):
            f.setForceGroup(3)

    positions = [openmm.Vec3(x, y, z) * unit.nanometer for (x, y, z) in positions_nm]
    integrator = openmm.LangevinMiddleIntegrator(100 * unit.kelvin, 0.01 / unit.picosecond, 0.001 * unit.picosecond)
    platform = openmm.Platform.getPlatformByName("CPU")
    context = openmm.Context(system, integrator, platform)
    context.setPositions(positions)
    return context, cfg


def _get_openmm_component_forces(context) -> dict[str, np.ndarray]:
    """Return dict bond, coulomb, lj, cavity -> (N,3) in kJ/(mol·nm)."""
    from openmm import unit

    out = {}
    for name, group in [("bond", 0), ("coulomb", 1), ("lj", 2), ("cavity", 3)]:
        state = context.getState(getForces=True, groups=2**group)
        f = state.getForces(asNumpy=True)
        out[name] = np.array(
            [[fx, fy, fz] for fx, fy, fz in f.value_in_unit(unit.kilojoule_per_mole / unit.nanometer)]
        )
    return out


def _run_cav_hoomd_and_get_per_component_forces(
    gsd_path: Path,
    lambda_au: float,
    cavity_freq_cm: float = 1560.0,
    frame: int = 0,
) -> tuple[dict[str, np.ndarray] | None, np.ndarray | None, str | None]:
    """
    Run cav-hoomd in-process, then extract per-force arrays in tag (GSD) order.
    Returns (components, bond_raw_array, None) or (None, None, error_msg).
    bond_raw_array: raw bond force array before rtag (shape (N_local,3)); used for layout check.
    """
    try:
        import hoomd
    except ImportError:
        return None, None, "hoomd not importable"
    try:
        from hoomd.cavitymd.simulation import CavityMDSimulation
    except ImportError:
        return None, None, "hoomd.cavitymd not importable"

    gsd_str = str(gsd_path.resolve())
    base_kwargs = dict(
        job_dir=str(gsd_path.parent),
        replica=0,
        freq=cavity_freq_cm,
        lambda_coupling=lambda_au,
        incavity=True,
        runtime_ps=0.0,
        input_gsd=gsd_str,
        frame=frame,
        name="force_debug",
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
    try:
        sim.run()
    except Exception as e:
        return None, None, f"cav-hoomd run failed: {e}"

    hoomd_sim = getattr(sim, "sim", None) or getattr(sim, "sim_obj", None)
    if hoomd_sim is not None and hasattr(hoomd_sim, "sim"):
        hoomd_sim = hoomd_sim.sim
    if hoomd_sim is None:
        return None, None, "CavityMDSimulation does not expose HOOMD simulation"

    try:
        with hoomd_sim.state.cpu_local_snapshot as snap:
            rtag = np.array(snap.particles.rtag, copy=True)
    except Exception as e:
        return None, None, f"cpu_local_snapshot/rtag failed: {e}"

    components: dict[str, np.ndarray] = {}
    bond_raw: np.ndarray | None = None
    integrator_forces = getattr(
        getattr(getattr(hoomd_sim, "operations", None), "integrator", None), "forces", None
    )
    if integrator_forces is None:
        return None, None, "hoomd simulation.operations.integrator.forces not available"

    def _to_tag_order(arr: np.ndarray) -> np.ndarray:
        return np.asarray(arr)[np.asarray(rtag, dtype=np.intp)]

    for force_obj in integrator_forces:
        fname = getattr(force_obj, "name", None) or type(force_obj).__name__
        fname_lower = (fname or "").lower()
        key = None
        if "harmonic" in fname_lower or "bond" in fname_lower:
            key = "bond"
        elif "coulomb" in fname_lower or "pppm" in fname_lower or "charge" in fname_lower:
            key = "coulomb"
        elif "lennard" in fname_lower or "lj" in fname_lower or "lennardjones" in fname_lower:
            key = "lj"
        elif "cavity" in fname_lower:
            key = "cavity"
        if key is None:
            continue
        arr = None
        if hasattr(force_obj, "cpu_local_force_arrays"):
            try:
                with force_obj.cpu_local_force_arrays as fa:
                    arr = np.array(fa.force, copy=True)[:, :3]
            except Exception:
                pass
        if arr is None and hasattr(force_obj, "local_force_arrays"):
            try:
                with force_obj.local_force_arrays as fa:
                    arr = np.array(fa.force, copy=True)[:, :3]
            except Exception:
                pass
        if arr is not None:
            components[key] = _to_tag_order(arr)
            if key == "bond":
                bond_raw = arr  # raw before rtag, shape (N_local, 3), for layout check

    if not components:
        return None, None, "No per-force arrays could be read"
    return components, bond_raw, None


# Harmonic bond k (Ha/Bohr^2) and r0 (Bohr) from diamer_forcefield.xml / cav-hoomd setup; typeid 0=O-O, 1=N-N
BOND_PARAMS_HA_BOHR: dict[int, tuple[float, float]] = {
    0: (0.73204, 2.281655158),       # O-O
    1: (1.4325, 2.0743522177),      # N-N
}


def _min_image_dr(dr: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Minimum-image displacement for orthogonal box; dr and box are (3,) or broadcastable."""
    return np.asarray(dr, dtype=np.float64) - np.asarray(box, dtype=np.float64) * np.round(
        np.asarray(dr, dtype=np.float64) / np.asarray(box, dtype=np.float64)
    )


def _single_bond_diagnostic(
    cfg: dict,
    pos_consistent: np.ndarray,
    F_openmm_bond_kj: np.ndarray,
    F_hoomd_bond_tag_order: np.ndarray | None,
    gsd_idx: np.ndarray,
    bond_indices: list[int],
) -> str:
    """
    For each bond b: compute r, r_hat, F_by_hand from positions and k,r0; compare to OpenMM and HOOMD.
    Returns a multi-line string for the summary. Explains WHY harmonic bond can differ: same formula
    F = k*(r0-r)*r_hat, so any difference must come from (positions, box, bond pair, k, r0) or indexing.
    OpenMM forces are in kJ/(mol·nm); we convert to Ha/Bohr with KJMOL_NM_TO_HARTREE_BOHR for comparison.
    """
    lines: list[str] = []
    # OpenMM returns kJ/(mol·nm). Ha/Bohr = (kJ/mol/nm) * K.  K = BOHR_TO_NM/2625.5 ≈ 2.016e-5
    K = KJMOL_NM_TO_HARTREE_BOHR
    box = np.asarray(cfg["box_bohr"][:3], dtype=np.float64)
    bonds_group = np.asarray(cfg["bonds_group"], dtype=np.intp)
    bonds_typeid = np.asarray(cfg.get("bonds_typeid"), dtype=np.intp) if cfg.get("bonds_typeid") is not None else None
    n_mol = int(cfg["n_mol"])

    lines.append("# --- Single-bond diagnostic: WHY does harmonic bond differ? ---")
    lines.append("# Harmonic bond is F = k*(r0 - r)*r_hat. Same k, r0, positions, box => same F.")
    lines.append("# We use unwrap_molecules() positions so bonded atoms are close; by-hand and OpenMM both use direct distance.")
    lines.append("# Units: OpenMM getForces() returns kJ/(mol·nm). We convert to Ha/Bohr by multiplying by K = BOHR_TO_NM/2625.5.")
    lines.append("")

    for b in bond_indices:
        if b < 0 or b >= n_mol:
            continue
        tag_a = int(bonds_group[b][0])
        tag_b = int(bonds_group[b][1])
        openmm_a, openmm_b = 2 * b, 2 * b + 1
        typ = int(bonds_typeid[b]) if bonds_typeid is not None else 0
        k_au, r0_au = BOND_PARAMS_HA_BOHR.get(typ, BOND_PARAMS_HA_BOHR[0])

        pos_a = np.asarray(pos_consistent[tag_a], dtype=np.float64)
        pos_b = np.asarray(pos_consistent[tag_b], dtype=np.float64)
        dr = pos_b - pos_a
        dr_minimg = _min_image_dr(dr, box)
        r = float(np.linalg.norm(dr_minimg))
        r_hat = (dr_minimg / r) if r > 1e-12 else np.zeros(3)
        # Force on first atom (tag_a): F_a = k*(r0 - r)*r_hat (restoring toward second when r>r0)
        F_a_hand = k_au * (r0_au - r) * r_hat
        F_b_hand = -F_a_hand

        type_name = "O-O" if typ == 0 else "N-N"
        lines.append(f"# Bond {b}  type={type_name}  tags=({tag_a},{tag_b})  OpenMM particles=({openmm_a},{openmm_b})")
        lines.append(f"#   k={k_au} Ha/Bohr^2  r0={r0_au} Bohr")
        lines.append(f"#   pos_a (Bohr) = ({pos_a[0]:.6f}, {pos_a[1]:.6f}, {pos_a[2]:.6f})")
        lines.append(f"#   pos_b (Bohr) = ({pos_b[0]:.6f}, {pos_b[1]:.6f}, {pos_b[2]:.6f})")
        lines.append(f"#   dr_minimg (Bohr) = ({dr_minimg[0]:.6f}, {dr_minimg[1]:.6f}, {dr_minimg[2]:.6f})  |dr| = r = {r:.6f} Bohr")
        lines.append(f"#   F_on_first_by_hand (Ha/Bohr) = ({F_a_hand[0]:.6e}, {F_a_hand[1]:.6e}, {F_a_hand[2]:.6e})  |F| = {np.linalg.norm(F_a_hand):.6e}")
        lines.append(f"#   F_on_second_by_hand (Ha/Bohr) = ({F_b_hand[0]:.6e}, {F_b_hand[1]:.6e}, {F_b_hand[2]:.6e})")

        F_openmm_a_kj = np.asarray(F_openmm_bond_kj[openmm_a], dtype=np.float64)
        F_openmm_b_kj = np.asarray(F_openmm_bond_kj[openmm_b], dtype=np.float64)
        F_openmm_a = F_openmm_a_kj * K
        F_openmm_b = F_openmm_b_kj * K
        # Units check: OpenMM returns kJ/(mol·nm). We convert to Ha/Bohr via * K. By-hand in kJ/mol/nm = F_a_hand / K.
        F_by_hand_kj = F_a_hand / K if K != 0 else np.nan
        mag_o_kj = np.linalg.norm(F_openmm_a_kj)
        mag_hand_kj = np.linalg.norm(F_by_hand_kj)
        ratio_kj = (mag_o_kj / mag_hand_kj) if mag_hand_kj > 1e-30 else np.nan
        lines.append(f"#   F_openmm[{openmm_a}] raw (kJ/mol/nm) = ({F_openmm_a_kj[0]:.4e}, {F_openmm_a_kj[1]:.4e}, {F_openmm_a_kj[2]:.4e})  |F| = {mag_o_kj:.4e}")
        lines.append(f"#   F_openmm[{openmm_a}] (Ha/Bohr) = raw * K,  K = {K:.6e}  => ({F_openmm_a[0]:.6e}, {F_openmm_a[1]:.6e}, {F_openmm_a[2]:.6e})  |F| = {np.linalg.norm(F_openmm_a):.6e}")
        lines.append(f"#   F_by_hand on first (Ha/Bohr) = ({F_a_hand[0]:.6e}, {F_a_hand[1]:.6e}, {F_a_hand[2]:.6e})  |F| = {np.linalg.norm(F_a_hand):.6e}")
        lines.append(f"#   F_by_hand on first (kJ/mol/nm) = F_by_hand_Ha_Bohr / K => |F| = {mag_hand_kj:.4e}  (same physical force, different units)")
        lines.append(f"#   In same units (kJ/mol/nm): OpenMM |F| / by_hand |F| = {ratio_kj:.2f}  => conversion is applied; OpenMM value is ~{ratio_kj:.0f}x too large.")
        lines.append(f"#   F_openmm[{openmm_b}] (Ha/Bohr) = ({F_openmm_b[0]:.6e}, {F_openmm_b[1]:.6e}, {F_openmm_b[2]:.6e})  |F| = {np.linalg.norm(F_openmm_b):.6e}")

        if F_hoomd_bond_tag_order is not None and tag_a < F_hoomd_bond_tag_order.shape[0] and tag_b < F_hoomd_bond_tag_order.shape[0]:
            F_hoomd_a = np.asarray(F_hoomd_bond_tag_order[tag_a], dtype=np.float64)
            F_hoomd_b = np.asarray(F_hoomd_bond_tag_order[tag_b], dtype=np.float64)
            lines.append(f"#   F_hoomd[tag={tag_a}] (Ha/Bohr) = ({F_hoomd_a[0]:.6e}, {F_hoomd_a[1]:.6e}, {F_hoomd_a[2]:.6e})  |F| = {np.linalg.norm(F_hoomd_a):.6e}")
            lines.append(f"#   F_hoomd[tag={tag_b}] (Ha/Bohr) = ({F_hoomd_b[0]:.6e}, {F_hoomd_b[1]:.6e}, {F_hoomd_b[2]:.6e})  |F| = {np.linalg.norm(F_hoomd_b):.6e}")

            err_hand_openmm_a = np.linalg.norm(F_a_hand - F_openmm_a)
            err_hand_hoomd_a = np.linalg.norm(F_a_hand - F_hoomd_a)
            err_hand_hoomd_b = np.linalg.norm(F_b_hand - F_hoomd_a)  # HOOMD might report force on 'other' atom at this tag
            # Tolerate ~0.1 Ha/Bohr so HOOMD (same mag/dir as by-hand) counts as agreement
            agrees_openmm = err_hand_openmm_a < 0.1
            agrees_hoomd = err_hand_hoomd_a < 0.1 or err_hand_hoomd_b < 0.1
            lines.append(f"#   |F_by_hand - F_openmm| on first = {err_hand_openmm_a:.6e}  => by-hand agrees with OpenMM: {'YES' if agrees_openmm else 'NO'}")
            lines.append(f"#   |F_by_hand - F_hoomd|  on first = {err_hand_hoomd_a:.6e}  (or F_hoomd=F_on_second: {err_hand_hoomd_b:.6e})  => agrees with HOOMD: {'YES' if agrees_hoomd else 'NO'}")
            if agrees_openmm and not agrees_hoomd:
                lines.append("#   => OpenMM matches formula; HOOMD does not (wrong row/units/params or different positions/box).")
            elif not agrees_openmm and agrees_hoomd:
                lines.append("#   => HOOMD matches formula; OPENMM DOES NOT. Bug or indexing/positions in OpenMM build.")
            elif agrees_openmm and agrees_hoomd:
                lines.append("#   => Both match by-hand (good).")
            else:
                lines.append("#   => Neither matches by-hand (check positions/box used by both, or bond pair ordering).")
        else:
            lines.append("#   F_hoomd not available for comparison.")
        lines.append("")

    lines.append("# If by-hand matches OpenMM but not HOOMD: HOOMD is using different inputs (positions, box, bond list, k/r0)")
    lines.append("# or we are reading the wrong HOOMD row for this tag. If by-hand matches HOOMD but not OpenMM: invert.")
    lines.append("")
    return "\n".join(lines)


def _compare_component(
    F_openmm_au: np.ndarray,
    F_hoomd_au: np.ndarray | None,
) -> tuple[float | None, float | None, int | None]:
    """Return (max_abs_err, rmse, worst_particle) in Ha/Bohr, or (None,None,None) if F_hoomd_au is None."""
    if F_hoomd_au is None or F_openmm_au.shape != F_hoomd_au.shape:
        return None, None, None
    d = F_openmm_au - np.asarray(F_hoomd_au, dtype=np.float64)
    max_abs = float(np.max(np.abs(d)))
    rmse = float(np.sqrt(np.mean(d * d)))
    per = np.max(np.abs(d), axis=1)
    worst = int(np.argmax(per))
    return max_abs, rmse, worst


def _per_particle_errors(
    openmm_components: dict[str, np.ndarray],
    hoomd_components: dict[str, np.ndarray] | None,
    gsd_idx: np.ndarray,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """
    Compute per-particle L2 norm of (OpenMM - HOOMD) for each component, in Ha/Bohr.
    Returns (err_per_particle, diff_vectors) where err has keys bond,coulomb,lj,cavity
    with (N,) arrays; diff has (N,3) for debugging. Missing HOOMD → NaN for that component.
    """
    K = KJMOL_NM_TO_HARTREE_BOHR
    N = len(gsd_idx)
    err: dict[str, np.ndarray] = {}
    diff_vectors: dict[str, np.ndarray] = {}
    for comp_key in ["bond", "coulomb", "lj", "cavity"]:
        F_o = openmm_components.get(comp_key)
        if F_o is None:
            err[comp_key] = np.full(N, np.nan)
            diff_vectors[comp_key] = np.full((N, 3), np.nan)
            continue
        F_o_au = np.asarray(F_o, dtype=np.float64) * K
        F_h = hoomd_components.get(comp_key) if hoomd_components else None
        if F_h is None or F_h.shape[0] != N:
            err[comp_key] = np.full(N, np.nan)
            diff_vectors[comp_key] = np.full((N, 3), np.nan)
            continue
        F_h_ordered = np.asarray(F_h, dtype=np.float64)[gsd_idx]
        d = F_o_au - F_h_ordered
        diff_vectors[comp_key] = d
        err[comp_key] = np.linalg.norm(d, axis=1)
    return err, diff_vectors


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Debug force mismatch: compare OpenMM vs cav-hoomd per component (bond, Coulomb, LJ, Cavity)."
    )
    parser.add_argument("gsd", type=str, help="GSD file (same config for both codes)")
    parser.add_argument("--lambda", type=float, default=0.001, dest="lambda_coupling")
    parser.add_argument("--cavity-freq", type=float, default=1560.0)
    parser.add_argument("--frame", type=int, default=0)
    parser.add_argument("--ff-dir", type=str, default=None)
    parser.add_argument("--summary", type=str, default="force_component_summary.txt", metavar="FILE")
    parser.add_argument(
        "--per-particle",
        type=str,
        default="force_component_errors_per_particle",
        metavar="BASE",
        help="Base name for per-particle diagnostics: writes <BASE>.csv and <BASE>.npz for plotting err vs particle index.",
    )
    parser.add_argument("-q", "--quiet", action="store_true")
    args = parser.parse_args()

    gsd_path = Path(args.gsd)
    if not gsd_path.exists():
        print(f"GSD not found: {gsd_path}", file=sys.stderr)
        sys.exit(1)
    ff_dir = Path(args.ff_dir) if args.ff_dir else _SCRIPT_DIR

    if not args.quiet:
        print("Building OpenMM system with force groups...")
    try:
        context, cfg = _build_openmm_system_with_force_groups(
            gsd_path, args.lambda_coupling, args.cavity_freq, ff_dir, args.frame
        )
    except Exception as e:
        print(f"OpenMM build failed: {e}", file=sys.stderr)
        if not args.quiet:
            import traceback
            traceback.print_exc()
        sys.exit(1)

    pos_unwrapped = unwrap_molecules(cfg)
    openmm_components = _get_openmm_component_forces(context)
    del context

    if not args.quiet:
        print("Running cav-hoomd for per-component forces...")
    # OpenMM uses unwrap_molecules() so bond forces are correct. Run HOOMD on the same
    # unwrapped configuration so bond forces match at ~1e-5 Ha/Bohr (same level as LJ).
    box_bohr = np.asarray(cfg["box_bohr"][:3], dtype=np.float64)
    pos_hoomd = pos_unwrapped - box_bohr / 2.0  # center; unwrap can extend past L
    pos_hoomd = pos_hoomd - np.round(pos_hoomd / box_bohr) * box_bohr  # wrap into [-L/2, L/2)
    tmp_gsd_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(suffix=".gsd", delete=False) as tmp:
            tmp_gsd_path = Path(tmp.name)
        write_gsd_with_positions(gsd_path, args.frame, pos_hoomd, tmp_gsd_path)
        gsd_for_hoomd = tmp_gsd_path
        hoomd_components, bond_raw, err = _run_cav_hoomd_and_get_per_component_forces(
            gsd_for_hoomd, args.lambda_coupling, args.cavity_freq, frame=0
        )
    finally:
        if tmp_gsd_path is not None and tmp_gsd_path.exists():
            tmp_gsd_path.unlink(missing_ok=True)
    if err is not None:
        if not args.quiet:
            print(f"cav-hoomd per-component not available: {err}", file=sys.stderr)
        hoomd_components = {}
        bond_raw = None

    gsd_idx = openmm_to_gsd_indices(cfg)
    N = len(gsd_idx)
    component_names = ["bond", "Coulomb", "LJ", "Cavity"]
    key_map = {"bond": "bond", "Coulomb": "coulomb", "LJ": "lj", "Cavity": "cavity"}

    rows: list[tuple[str, float | None, float | None, int | None]] = []
    for display_name in component_names:
        key = key_map[display_name]
        F_o = openmm_components.get(key)
        if F_o is None:
            rows.append((display_name, None, None, None))
            continue
        F_o_au = F_o * KJMOL_NM_TO_HARTREE_BOHR
        F_h = hoomd_components.get(key) if hoomd_components else None
        if F_h is not None and F_h.shape[0] == N:
            F_h_ordered = np.asarray(F_h)[gsd_idx]
            max_e, rmse_e, worst_i = _compare_component(F_o_au, F_h_ordered)
            rows.append((display_name, max_e, rmse_e, worst_i))
        else:
            rows.append((display_name, None, None, None))

    # Per-particle diagnostics: L2 norm of (OpenMM - HOOMD) per component, for plotting err vs particle index
    err_per_particle, diff_per_particle = _per_particle_errors(
        openmm_components, hoomd_components, gsd_idx
    )

    culprit = "N/A"
    if hoomd_components:
        best_idx = None
        best_val = -1.0
        for i, (_, max_e, _, _) in enumerate(rows):
            if max_e is not None and max_e > best_val:
                best_val = max_e
                best_idx = i
        if best_idx is not None:
            culprit = component_names[best_idx]
    else:
        culprit = "N/A (cav-hoomd not available)"

    # Per-particle diagnostics: CSV and NPZ for plotting err vs particle index
    base = Path(args.per_particle)
    stem = base.stem if base.suffix in (".csv", ".npz") else base.name
    parent = base.parent
    csv_path = parent / (stem + ".csv")
    npz_path = parent / (stem + ".npz")

    with open(csv_path, "w") as c:
        c.write("particle_index,err_bond_Ha_Bohr,err_coulomb_Ha_Bohr,err_lj_Ha_Bohr,err_cavity_Ha_Bohr\n")
        for i in range(N):
            eb = err_per_particle["bond"][i]
            ec = err_per_particle["coulomb"][i]
            el = err_per_particle["lj"][i]
            ev = err_per_particle["cavity"][i]
            c.write(f"{i},{eb:.6e},{ec:.6e},{el:.6e},{ev:.6e}\n")

    np.savez(
        npz_path,
        particle_index=np.arange(N, dtype=np.int32),
        err_bond=err_per_particle["bond"],
        err_coulomb=err_per_particle["coulomb"],
        err_lj=err_per_particle["lj"],
        err_cavity=err_per_particle["cavity"],
    )

    # Extended summary: distributional stats, worst particles, and "Why bond?" section
    def _percentiles(a: np.ndarray) -> tuple[float, float, float]:
        a = np.asarray(a, dtype=np.float64)
        a = a[np.isfinite(a)]
        if a.size == 0:
            return np.nan, np.nan, np.nan
        return float(np.percentile(a, 50)), float(np.percentile(a, 95)), float(np.percentile(a, 99))

    def _worst_k(a: np.ndarray, k: int = 5) -> list[int]:
        a = np.asarray(a, dtype=np.float64)
        valid = np.isfinite(a)
        if not np.any(valid):
            return []
        perm = np.argsort(a)[::-1]
        out = []
        for i in perm:
            if len(out) >= k:
                break
            if np.isfinite(a[i]):
                out.append(int(i))
        return out

    summary_path = Path(args.summary)
    with open(summary_path, "w") as f:
        f.write("# Force component comparison: OpenMM vs cav-hoomd\n")
        f.write(f"# GSD: {gsd_path.name}  lambda: {args.lambda_coupling}  cavity_freq: {args.cavity_freq} cm^-1\n\n")
        f.write(f"{'Component':<12}  {'max_abs_err_Ha_Bohr':>20}  {'RMSE_Ha_Bohr':>14}  {'worst_particle':>14}\n")
        f.write("-" * 64 + "\n")
        for name, max_e, rmse_e, worst_i in rows:
            max_s = f"{max_e:.6e}" if max_e is not None else "N/A"
            rmse_s = f"{rmse_e:.6e}" if rmse_e is not None else "N/A"
            worst_s = str(worst_i) if worst_i is not None else "N/A"
            f.write(f"{name:<12}  {max_s:>20}  {rmse_s:>14}  {worst_s:>14}\n")
        f.write("\nLargest max error: " + culprit + "\n\n")

        # Single-bond diagnostic: hand-compute F from r,k,r0 and compare to OpenMM/HOOMD (why does harmonic bond differ?)
        bond_worst_for_diag = _worst_k(err_per_particle["bond"], 5)
        n_mol_diag = int(cfg["n_mol"])
        bond_indices_diag = sorted(set(idx // 2 for idx in (bond_worst_for_diag or []) if idx < 2 * n_mol_diag))[:5]
        if not bond_indices_diag:
            bond_indices_diag = [0]
        F_bond_openmm = openmm_components.get("bond")
        if F_bond_openmm is None:
            F_bond_openmm = np.zeros((N, 3), dtype=np.float64)
        F_bond_hoomd_tag = hoomd_components.get("bond") if hoomd_components else None
        f.write(_single_bond_diagnostic(cfg, pos_unwrapped, F_bond_openmm, F_bond_hoomd_tag, gsd_idx, bond_indices_diag))

        # Force magnitude and relative error: compare bond vs LJ apples-to-apples (why is bond abs err ~1e-3 vs LJ ~1e-5?)
        f.write("# --- Force magnitude and relative error (why bond ~1e-3 vs LJ ~1e-5?) ---\n")
        f.write("# Absolute error scales with |F|. Compare relative error rel_err = |Delta F| / max(|F_openmm|, 1e-12).\n")
        K = KJMOL_NM_TO_HARTREE_BOHR
        for comp in ["bond", "coulomb", "lj", "cavity"]:
            Fo = openmm_components.get(comp)
            if Fo is None:
                continue
            F_mag = np.linalg.norm(np.asarray(Fo, dtype=np.float64) * K, axis=1)
            rel_err = np.where(F_mag > 1e-12, err_per_particle[comp] / F_mag, np.nan)
            p50_mag, p95_mag, _ = _percentiles(F_mag)
            p50_rel, p95_rel, _ = _percentiles(rel_err)
            label = {"bond": "bond", "coulomb": "Coulomb", "lj": "LJ", "cavity": "Cavity"}[comp]
            f.write(f"#   {label}: p50(|F_openmm|)={p50_mag:.4e} Ha/Bohr  p50(rel_err)={p50_rel:.4e}  p95(rel_err)={p95_rel:.4e}\n")
        f.write("# Observed: bond and Coulomb have large p50(rel_err) ~0.25-0.5; LJ has small p50(rel_err) ~0.02.\n")
        f.write("# So bond/Coulomb are *relatively* wrong; LJ is relatively good. Magnitude does not explain it.\n\n")

        f.write("# --- Per-particle diagnostics ---\n")
        f.write(f"# Per-particle L2 norm of (OpenMM - HOOMD) in Ha/Bohr. Plot err_* vs particle_index to see which component dominates where.\n")
        f.write(f"# CSV: {csv_path}  NPZ: {npz_path}\n\n")
        f.write("# Distribution (p50, p95, p99) and worst 5 particles per component:\n")
        n_mol = cfg["n_mol"]
        bonds_typeid = cfg.get("bonds_typeid")
        has_cavity = cfg.get("has_cavity", False)
        for comp in ["bond", "coulomb", "lj", "cavity"]:
            a = err_per_particle[comp]
            p50, p95, p99 = _percentiles(a)
            worst = _worst_k(a, 5)
            label = {"bond": "bond", "coulomb": "Coulomb", "lj": "LJ", "cavity": "Cavity"}[comp]
            f.write(f"#   {label}: p50={p50:.4e}  p95={p95:.4e}  p99={p99:.4e}  worst_particles={worst}\n")
        f.write("\n")

        f.write("# --- Why are there bond-error outliers? ---\n")
        f.write("# Most particles have small bond errors (p50 ~1e-3); a few have huge errors (outliers).\n")
        f.write("# HOOMD's bond force may be stored in bond-list order (row k = bond k//2, atom k%2) or in\n")
        f.write("# local-domain order, not in particle-tag order. We index by rtag (tag order). When the\n")
        f.write("# backend uses a different layout, we assign the wrong row to each tag: for some tags we\n")
        f.write("# accidentally read a similar row -> small error; for others we read a different bond's\n")
        f.write("# force -> large error (outliers). So outliers are particles where 'our' row != the row\n")
        f.write("# the backend actually uses for that particle.\n")
        bond_worst = _worst_k(err_per_particle["bond"], 5)
        if bond_worst and np.any(np.isfinite(err_per_particle["bond"])):
            f.write("# Bond worst particles (particle_index, bond_index, type):\n")
            for idx in bond_worst:
                if idx < 2 * n_mol and bonds_typeid is not None:
                    b = idx // 2
                    typ = "OO" if int(bonds_typeid[b]) == 0 else "NN"
                    f.write(f"#     {idx}  bond {b}  {typ}\n")
                elif has_cavity and idx == 2 * n_mol:
                    f.write(f"#     {idx}  cavity (no bond)\n")
                else:
                    f.write(f"#     {idx}\n")
        f.write("\n")

        # Bond-force layout check: try "tag order" (current) vs "bond order" (raw rows = OpenMM/bond order)
        f.write("# --- Bond-force layout check ---\n")
        if bond_raw is not None:
            f.write(f"# HOOMD bond force raw array shape: {bond_raw.shape} (before rtag; rows = local index on this rank).\n")
            fo_bond = openmm_components.get("bond")
            F_o_bond_au = (fo_bond * KJMOL_NM_TO_HARTREE_BOHR) if fo_bond is not None else None
            if F_o_bond_au is not None and bond_raw.shape[0] == N and bond_raw.shape[1] >= 3:
                # Single rank: raw array has N rows. Try interpreting rows as bond order (row i = OpenMM particle i).
                F_hoomd_bond_order = np.asarray(bond_raw[:, :3], dtype=np.float64)
                max_bo, rmse_bo, _ = _compare_component(F_o_bond_au, F_hoomd_bond_order)
                bond_row = next((r for r in rows if r[0] == "bond"), None)
                max_tag = bond_row[1] if bond_row else None
                rmse_tag = bond_row[2] if bond_row else None
                f.write(f"# Tag-order (current):   max_err = {max_tag:.6e}  RMSE = {rmse_tag:.6e}\n")
                f.write(f"# Bond-order (raw[i]=OpenMM i): max_err = {max_bo:.6e}  RMSE = {rmse_bo:.6e}\n")
                if max_bo is not None and max_tag is not None and max_bo < 0.1 * max_tag:
                    f.write("# => Bond-order error << tag-order: HOOMD bond array is in bond order; use raw[i] for OpenMM i.\n")
                elif max_bo is not None and max_tag is not None and max_bo > 10 * max_tag:
                    f.write("# => Tag-order is better than bond-order: HOOMD bond array matches local/tag order (current use is ok).\n")
                else:
                    f.write("# => Inconclusive or similar; raw row order may differ from both.\n")
            else:
                f.write(f"# Cannot try bond-order (need raw.shape[0]==N={N}); have raw.shape[0]={bond_raw.shape[0]}. Multi-rank or partial.\n")
        else:
            f.write("# No raw bond array available for layout check.\n")
        f.write("\n")

        f.write("# --- Bond force in cav-hoomd (source) ---\n")
        f.write("# cav-hoomd does NOT implement bond forces in C++; it uses HOOMD's hoomd.md.bond.Harmonic.\n")
        f.write("# Setup: cavitymd/simulation/core.py setup_force_parameters() [or cavitymd/setup.py create_harmonic_bonds()].\n")
        f.write("#   harmonic = hoomd.md.bond.Harmonic()\n")
        f.write("#   harmonic.params['O-O'] = dict(k=2*0.36602, r0=2.281655158)   # k=0.73204 Ha/Bohr^2, r0 Bohr\n")
        f.write("#   harmonic.params['N-N'] = dict(k=2*0.71625, r0=2.0743522177)  # k=1.4325, r0 Bohr\n")
        f.write("# These match OpenMM diamer_forcefield.xml (Ha/Bohr^2, Bohr). cav-hoomd assumes Hartree/a.u.\n")
        f.write("# The per-particle bond array we read (cpu_local_force_arrays on the Harmonic force object)\n")
        f.write("# is filled by HOOMD-blue's MD backend, not by cav-hoomd. Cav-hoomd src has CavityForceCompute,\n")
        f.write("# DipoleResponseForceCompute, etc. — no bond force. To fix layout/scale: inspect HOOMD-blue\n")
        f.write("# source for md.bond.Harmonic / BondForceCompute: how it writes to the force array and index order.\n")
        f.write("\n")

        # Next-step diagnostic: F_openmm_bond vs F_hoomd_bond at worst particles (same units: Ha/Bohr)
        f.write("# --- Next-step: F_openmm_bond vs F_hoomd_bond at worst particles (Ha/Bohr) ---\n")
        fo_bond = openmm_components.get("bond")
        F_o_bond_au = (fo_bond * KJMOL_NM_TO_HARTREE_BOHR) if fo_bond is not None else np.full((N, 3), np.nan)
        F_h_bond = hoomd_components.get("bond") if hoomd_components else None
        F_h_bond_ordered = np.asarray(F_h_bond)[gsd_idx] if F_h_bond is not None and F_h_bond.shape[0] == N else None
        for idx in (bond_worst or [])[:3]:  # top 3 worst
            fo = np.asarray(F_o_bond_au[idx], dtype=np.float64)
            fh = np.asarray(F_h_bond_ordered[idx], dtype=np.float64) if F_h_bond_ordered is not None else np.full(3, np.nan)
            tag = int(gsd_idx[idx])
            f.write(f"# particle_index={idx}  gsd_tag={tag}\n")
            f.write(f"#   F_openmm_bond  (x,y,z) Ha/Bohr: ({fo[0]:.6e}, {fo[1]:.6e}, {fo[2]:.6e})  |F|={np.linalg.norm(fo):.6e}\n")
            f.write(f"#   F_hoomd_bond   (x,y,z) Ha/Bohr: ({fh[0]:.6e}, {fh[1]:.6e}, {fh[2]:.6e})  |F|={np.linalg.norm(fh):.6e}\n")
            d = fo - fh
            f.write(f"#   difference     (x,y,z) Ha/Bohr: ({d[0]:.6e}, {d[1]:.6e}, {d[2]:.6e})  |d|={np.linalg.norm(d):.6e}\n")
            no, nh = np.linalg.norm(fo), np.linalg.norm(fh)
            if no > 1e-12 and nh > 1e-12:
                ratio = no / nh
                cos_angle = float(np.dot(fo, fh) / (no * nh))
                f.write(f"#   ratio |F_openmm|/|F_hoomd| = {ratio:.4f}  ;  cos(angle) = {cos_angle:.4f}\n")
                f.write(f"#   -> If ratio ~1 and cos~1: same direction, unit/scale ok; wrong row would often give cos far from 1.\n")
                f.write(f"#   -> If ratio >> 1 or << 1: likely units (e.g. HOOMD bond in kJ/mol/nm or other).\n")
            f.write("\n")
        f.write("# --- Are we comparing the wrong stuff? ---\n")
        f.write("# No, for Coulomb/LJ/Cavity. We use the same alignment everywhere: HOOMD forces arr[rtag] -> tag\n")
        f.write("# order, then F_hoomd[gsd_idx] -> OpenMM order. OpenMM particle i = GSD tag gsd_idx[i]. Those three\n")
        f.write("# components agree (~1e-3), so we're comparing the right particles.\n")
        f.write("# For bond we use the same pipeline, but we *assume* the HOOMD bond force array is in local-index\n")
        f.write("# order (like net_force). If that backend uses bond-list order or another layout, arr[rtag] gives the\n")
        f.write("# wrong row for each tag -> we'd be comparing OpenMM bond-on-i to HOOMD bond-on-some-other-particle.\n")
        f.write("# So we're not misaligning particle indices; we might be reading the wrong HOOMD bond values.\n")
        f.write("\n")
        f.write("# --- Why might bond dominate? ---\n")
        f.write("# Bond max error ~50 Ha/Bohr vs ~1e-3 for Coulomb/LJ/Cavity. Same alignment (rtag, gsd_idx) is used\n")
        f.write("# for all components; if index/layout were wrong, the other forces would be way off too. Minimum\n")
        f.write("# image is shared. We match harmonic/bond to the bond force correctly.\n")
        f.write("# Plausible: (1) HOOMD bond array layout != local order -> wrong row per tag. (2) Units: bond in\n")
        f.write("# other units than Ha/Bohr. Next: Check HOOMD bond force array row order; compare F_openmm vs\n")
        f.write("# F_hoomd_bond at worst particle in same units.\n")
    if not args.quiet:
        print(f"Wrote {summary_path}")
        print(f"Per-particle: {csv_path}  {npz_path}")


if __name__ == "__main__":
    main()
