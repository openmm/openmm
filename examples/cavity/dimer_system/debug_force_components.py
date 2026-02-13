#!/usr/bin/env python3
"""
Debug OpenMM vs cav-hoomd force mismatch by comparing each force component separately.

Runs four cases in one go: (1) harmonic bond only, (2) Coulomb only, (3) LJ only,
(4) Cavity force only. For each component, compares OpenMM vs cav-hoomd and reports
max_abs_err (Ha/Bohr), RMSE, worst particle. Writes a summary to a text file.

Outputs:
  - Summary TXT (component stats, worst-particle diagnostics)
  - Per-particle error CSV/NPZ (err_bond, err_coulomb, err_lj, err_cavity)
  - Raw forces TXT/NPZ (Fx,Fy,Fz per particle per component from both codes)
  - PNGs: force_component_errors.png, _cdf.png, _scatter.png, _magnitude.png

Usage:
  python debug_force_components.py init-0.gsd --lambda 0.001 --summary force_component_summary.txt
  python debug_force_components.py init-0.gsd --no-raw-forces --no-plots  # minimal output

Requires: OpenMM (cavity build), numpy, gsd. cav-hoomd is run in-process when importable.
Matplotlib optional (for PNG generation).
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
        wrap_positions_into_box,
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
    cutoff_nm: float = 1.0,
):
    """Build OpenMM system from GSD, assign force groups 0..3 to bond/Coulomb/LJ/Cavity, return (context, cfg)."""
    from openmm import openmm
    from openmm import unit
    from openmm.app import CutoffPeriodic,  Ewald, PME, Element, ForceField, Topology

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
        #nonbondedMethod=CutoffPeriodic,

        nonbondedMethod=PME,
        nonbondedCutoff=cutoff_nm * unit.nanometer,
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
) -> tuple[
    dict[str, np.ndarray] | None,
    np.ndarray | None,
    np.ndarray | None,
    str | None,
]:
    """
    Run cav-hoomd in-process, then extract per-force arrays in tag (GSD) order.
    Returns (components, bond_raw_array, pos_tag_order, None) or (None, None, None, error_msg).
    - components: per-force arrays indexed by particle tag (same order as GSD).
    - bond_raw_array: raw bond force before rtag (shape (N_local,3)); for layout check.
    - pos_tag_order: positions from HOOMD state (after step) in tag order, same length units as
      box (Bohr if GSD was written in Bohr). None if unavailable.
    Force array layout: HOOMD stores forces in local particle index order; rtag[tag] = local index,
    so tag-order force = arr[rtag]. This matches HOOMD-blue PotentialBond (writes to h_force.data[idx]).
    """
    try:
        import hoomd
    except ImportError:
        return None, None, None, "hoomd not importable"
    try:
        from hoomd.cavitymd.simulation import CavityMDSimulation
    except ImportError:
        return None, None, None, "hoomd.cavitymd not importable"

    gsd_str = str(gsd_path.resolve())
    # runtime_ps=0 → 0 steps so forces are at initial GSD config (matches run_cav_hoomd_forces_inprocess
    # and OpenMM). Previously used 1 step; that made F_hoomd differ from optimization's reference.
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
        return None, None, None, f"cav-hoomd run failed: {e}"

    hoomd_sim = getattr(sim, "sim", None) or getattr(sim, "sim_obj", None)
    if hoomd_sim is not None and hasattr(hoomd_sim, "sim"):
        hoomd_sim = hoomd_sim.sim
    if hoomd_sim is None:
        return None, None, None, "CavityMDSimulation does not expose HOOMD simulation"

    pos_tag_order: np.ndarray | None = None
    try:
        with hoomd_sim.state.cpu_local_snapshot as snap:
            rtag = np.array(snap.particles.rtag, copy=True)
            # Positions in local order; tag order = pos[rtag] so diagnostic uses same config HOOMD used
            pos_local = np.array(snap.particles.position, copy=True)
            if pos_local.ndim == 2 and pos_local.shape[1] >= 3:
                pos_tag_order = np.asarray(pos_local, dtype=np.float64)[np.asarray(rtag, dtype=np.intp), :3]
    except Exception as e:
        return None, None, None, f"cpu_local_snapshot/rtag failed: {e}"

    components: dict[str, np.ndarray] = {}
    bond_raw: np.ndarray | None = None
    integrator_forces = getattr(
        getattr(getattr(hoomd_sim, "operations", None), "integrator", None), "forces", None
    )
    if integrator_forces is None:
        return None, None, pos_tag_order, "hoomd simulation.operations.integrator.forces not available"

    def _to_tag_order(arr: np.ndarray) -> np.ndarray:
        return np.asarray(arr)[np.asarray(rtag, dtype=np.intp)]

    for force_obj in integrator_forces:
        fname = getattr(force_obj, "name", None) or type(force_obj).__name__
        fname_lower = (fname or "").lower()
        key = None
        if "harmonic" in fname_lower or "bond" in fname_lower:
            key = "bond"
        # PPPM: short-range = hoomd.md.pair.Ewald, long-range = hoomd.md.long_range.pppm.Coulomb; sum both
        elif (
            "coulomb" in fname_lower
            or "pppm" in fname_lower
            or "charge" in fname_lower
            or "ewald" in fname_lower
        ):
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
            tag_order_arr = _to_tag_order(arr)
            if key == "coulomb":
                if "coulomb" not in components:
                    components["coulomb"] = tag_order_arr.copy()
                else:
                    components["coulomb"] = components["coulomb"] + tag_order_arr
            else:
                components[key] = tag_order_arr
            if key == "bond":
                bond_raw = arr  # raw before rtag, shape (N_local, 3), for layout check

    if not components:
        return None, None, pos_tag_order, "No per-force arrays could be read"
    return components, bond_raw, pos_tag_order, None


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
    pos_hoomd: np.ndarray | None = None,
) -> str:
    """
    For each bond b: compute r, r_hat, F_by_hand from positions and k,r0; compare to OpenMM and HOOMD.
    When pos_hoomd is provided, also output (r0-r) from OpenMM vs HOOMD geometry and hand-F from each
    geometry vs actual forces, to isolate geometry vs formula/units.
    OpenMM forces are in kJ/(mol·nm); we convert to Ha/Bohr with KJMOL_NM_TO_HARTREE_BOHR for comparison.
    """
    lines: list[str] = []
    K = KJMOL_NM_TO_HARTREE_BOHR
    box = np.asarray(cfg["box_bohr"][:3], dtype=np.float64)
    bonds_group = np.asarray(cfg["bonds_group"], dtype=np.intp)
    bonds_typeid = np.asarray(cfg.get("bonds_typeid"), dtype=np.intp) if cfg.get("bonds_typeid") is not None else None
    n_mol = int(cfg["n_mol"])

    lines.append("# --- Single-bond diagnostic: WHY does harmonic bond differ? ---")
    lines.append("# Harmonic bond is F = k*(r0 - r)*r_hat. Same k, r0, positions, box => same F.")
    lines.append("")

    for b in bond_indices:
        if b < 0 or b >= n_mol:
            continue
        tag_a = int(bonds_group[b][0])
        tag_b = int(bonds_group[b][1])
        openmm_a, openmm_b = 2 * b, 2 * b + 1
        typ = int(bonds_typeid[b]) if bonds_typeid is not None else 0
        k_au, r0_au = BOND_PARAMS_HA_BOHR.get(typ, BOND_PARAMS_HA_BOHR[0])

        # --- OpenMM geometry (pos_consistent = pos_unwrapped) ---
        pos_a_openmm = np.asarray(pos_consistent[tag_a], dtype=np.float64)
        pos_b_openmm = np.asarray(pos_consistent[tag_b], dtype=np.float64)
        dr_openmm = pos_b_openmm - pos_a_openmm
        dr_openmm = _min_image_dr(dr_openmm, box)
        r_openmm = float(np.linalg.norm(dr_openmm))
        r_hat_openmm = (dr_openmm / r_openmm) if r_openmm > 1e-12 else np.zeros(3)
        r0_r_openmm = r0_au - r_openmm
        F_hand_openmm = k_au * r0_r_openmm * r_hat_openmm

        # --- (r0-r) and hand-force isolation block ---
        lines.append("# --- (r0-r) and hand-force isolation ---")
        lines.append(f"# Bond {b}  type={'O-O' if typ == 0 else 'N-N'}  tags=({tag_a},{tag_b})  r0={r0_au} Bohr  k={k_au} Ha/Bohr^2")
        lines.append("# OpenMM geometry (pos_unwrapped):")
        lines.append(f"#   r_openmm = {r_openmm:.6e} Bohr   (r0-r)_openmm = {r0_r_openmm:.6e} Bohr")

        if pos_hoomd is not None and tag_a < pos_hoomd.shape[0] and tag_b < pos_hoomd.shape[0]:
            pos_a_hoomd = np.asarray(pos_hoomd[tag_a], dtype=np.float64)
            pos_b_hoomd = np.asarray(pos_hoomd[tag_b], dtype=np.float64)
            dr_hoomd = pos_b_hoomd - pos_a_hoomd
            dr_hoomd = _min_image_dr(dr_hoomd, box)
            r_hoomd = float(np.linalg.norm(dr_hoomd))
            r_hat_hoomd = (dr_hoomd / r_hoomd) if r_hoomd > 1e-12 else np.zeros(3)
            r0_r_hoomd = r0_au - r_hoomd
            F_hand_hoomd = k_au * r0_r_hoomd * r_hat_hoomd

            lines.append("# HOOMD geometry (pos_hoomd):")
            lines.append(f"#   r_hoomd = {r_hoomd:.6e} Bohr   (r0-r)_hoomd = {r0_r_hoomd:.6e} Bohr")
            diff_geom = r0_r_openmm - r0_r_hoomd
            rel_geom = (diff_geom / r0_r_openmm) if abs(r0_r_openmm) > 1e-12 else float("nan")
            lines.append(f"#   (r0-r)_openmm vs (r0-r)_hoomd:  diff = {diff_geom:.6e} Bohr   rel_diff = {rel_geom:.6e}")

        F_openmm_a_kj = np.asarray(F_openmm_bond_kj[openmm_a], dtype=np.float64)
        F_openmm_b_kj = np.asarray(F_openmm_bond_kj[openmm_b], dtype=np.float64)
        F_openmm_a = F_openmm_a_kj * K
        F_openmm_b = F_openmm_b_kj * K
        mag_F_openmm = float(np.linalg.norm(F_openmm_a))
        r0_r_from_F_openmm = (mag_F_openmm / k_au) if k_au > 1e-30 else float("nan")
        lines.append("# (r0-r) implied by forces (|F|/k):")
        lines.append(f"#   (r0-r)_from_F_openmm = {r0_r_from_F_openmm:.6e} Bohr")

        if F_hoomd_bond_tag_order is not None and tag_a < F_hoomd_bond_tag_order.shape[0] and tag_b < F_hoomd_bond_tag_order.shape[0]:
            F_hoomd_a = np.asarray(F_hoomd_bond_tag_order[tag_a], dtype=np.float64)
            F_hoomd_b = np.asarray(F_hoomd_bond_tag_order[tag_b], dtype=np.float64)
            mag_F_hoomd = float(np.linalg.norm(F_hoomd_a))
            r0_r_from_F_hoomd = (mag_F_hoomd / k_au) if k_au > 1e-30 else float("nan")
            lines.append(f"#   (r0-r)_from_F_hoomd = {r0_r_from_F_hoomd:.6e} Bohr")

        lines.append("# Hand F from OpenMM geometry vs F_openmm (on first atom):")
        err_hand_vs_openmm = float(np.linalg.norm(F_hand_openmm - F_openmm_a))
        lines.append(f"#   |F_hand_openmm - F_openmm_a| = {err_hand_vs_openmm:.6e} Ha/Bohr")

        if (
            pos_hoomd is not None
            and F_hoomd_bond_tag_order is not None
            and tag_a < pos_hoomd.shape[0]
            and tag_b < pos_hoomd.shape[0]
            and tag_a < F_hoomd_bond_tag_order.shape[0]
        ):
            F_hoomd_a = np.asarray(F_hoomd_bond_tag_order[tag_a], dtype=np.float64)
            err_hand_vs_hoomd = float(np.linalg.norm(F_hand_hoomd - F_hoomd_a))
            lines.append("# Hand F from HOOMD geometry vs F_hoomd (on first atom):")
            lines.append(f"#   |F_hand_hoomd - F_hoomd_a| = {err_hand_vs_hoomd:.6e} Ha/Bohr")
            # Flag when HOOMD force disagrees with the same geometry OpenMM uses
            if err_hand_vs_hoomd > 1e-6:
                lines.append("#   => HOOMD bond force inconsistent with shared geometry (same r, r0): likely wrong")
                lines.append("#      force array index/order (bond-list vs tag) or different units; OpenMM matches F=k(r0-r)r_hat.")
        lines.append("")

        # --- Legacy single-bond block (pos, dr, by-hand, F_openmm, F_hoomd, agrees) ---
        type_name = "O-O" if typ == 0 else "N-N"
        lines.append(f"# Bond {b}  type={type_name}  tags=({tag_a},{tag_b})  OpenMM particles=({openmm_a},{openmm_b})")
        lines.append(f"#   k={k_au} Ha/Bohr^2  r0={r0_au} Bohr")
        lines.append(f"#   pos_a (Bohr) = ({pos_a_openmm[0]:.6f}, {pos_a_openmm[1]:.6f}, {pos_a_openmm[2]:.6f})")
        lines.append(f"#   pos_b (Bohr) = ({pos_b_openmm[0]:.6f}, {pos_b_openmm[1]:.6f}, {pos_b_openmm[2]:.6f})")
        lines.append(f"#   dr_minimg (Bohr) = ({dr_openmm[0]:.6f}, {dr_openmm[1]:.6f}, {dr_openmm[2]:.6f})  |dr| = r = {r_openmm:.6f} Bohr")
        F_a_hand = F_hand_openmm
        F_b_hand = -F_a_hand
        lines.append(f"#   F_on_first_by_hand (Ha/Bohr) = ({F_a_hand[0]:.6e}, {F_a_hand[1]:.6e}, {F_a_hand[2]:.6e})  |F| = {np.linalg.norm(F_a_hand):.6e}")
        lines.append(f"#   F_on_second_by_hand (Ha/Bohr) = ({F_b_hand[0]:.6e}, {F_b_hand[1]:.6e}, {F_b_hand[2]:.6e})")
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
            err_hand_hoomd_b = np.linalg.norm(F_b_hand - F_hoomd_a)
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

    lines.append("# If (r0-r)_openmm != (r0-r)_hoomd: geometry/r0/positions differ. If hand F != actual F for a code: formula/units/convention in that code.")
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


def _write_raw_forces_txt(
    path: Path,
    openmm_components: dict[str, np.ndarray],
    hoomd_components: dict[str, np.ndarray] | None,
    gsd_idx: np.ndarray,
    err_per_particle: dict[str, np.ndarray],
) -> None:
    """Write per-particle, per-component raw forces (Ha/Bohr) in block-per-component format."""
    K = KJMOL_NM_TO_HARTREE_BOHR
    N = len(gsd_idx)
    comp_labels = [("bond", "Bond"), ("coulomb", "Coulomb"), ("lj", "LJ"), ("cavity", "Cavity")]
    with open(path, "w") as f:
        f.write("# Per-particle, per-component raw forces (Ha/Bohr). Reference: cav-hoomd.\n")
        for comp_key, display_name in comp_labels:
            F_o = openmm_components.get(comp_key)
            if F_o is None:
                continue
            F_o_au = np.asarray(F_o, dtype=np.float64) * K
            F_h = hoomd_components.get(comp_key) if hoomd_components else None
            F_h_ordered = np.asarray(F_h, dtype=np.float64)[gsd_idx] if (F_h is not None and F_h.shape[0] == N) else np.full((N, 3), np.nan)
            err = err_per_particle.get(comp_key, np.full(N, np.nan))
            f.write(f"\n# {display_name} forces (Ha/Bohr)\n")
            f.write("# particle  Fx_openmm  Fy_openmm  Fz_openmm  Fx_hoomd  Fy_hoomd  Fz_hoomd  err_L2\n")
            for i in range(N):
                ox, oy, oz = F_o_au[i, 0], F_o_au[i, 1], F_o_au[i, 2]
                hx, hy, hz = F_h_ordered[i, 0], F_h_ordered[i, 1], F_h_ordered[i, 2]
                e = err[i] if np.isfinite(err[i]) else np.nan
                f.write(f"{i}  {ox:.6e}  {oy:.6e}  {oz:.6e}  {hx:.6e}  {hy:.6e}  {hz:.6e}  {e:.6e}\n")


def _write_raw_forces_npz(
    path: Path,
    openmm_components: dict[str, np.ndarray],
    hoomd_components: dict[str, np.ndarray] | None,
    gsd_idx: np.ndarray,
    err_per_particle: dict[str, np.ndarray],
) -> None:
    """Save raw force arrays to NPZ (Ha/Bohr). HOOMD reordered to OpenMM particle order."""
    K = KJMOL_NM_TO_HARTREE_BOHR
    N = len(gsd_idx)
    save_dict = {"particle_index": np.arange(N, dtype=np.int32)}
    for comp in ["bond", "coulomb", "lj", "cavity"]:
        F_o = openmm_components.get(comp)
        if F_o is not None:
            save_dict[f"F_openmm_{comp}"] = np.asarray(F_o, dtype=np.float64) * K
        F_h = hoomd_components.get(comp) if hoomd_components else None
        if F_h is not None and F_h.shape[0] == N:
            save_dict[f"F_hoomd_{comp}"] = np.asarray(F_h, dtype=np.float64)[gsd_idx]
        save_dict[f"err_{comp}"] = err_per_particle.get(comp, np.full(N, np.nan))
    np.savez(path, **save_dict)


def _generate_plots(
    plots_dir: Path,
    err_per_particle: dict[str, np.ndarray],
    openmm_components: dict[str, np.ndarray],
    hoomd_components: dict[str, np.ndarray] | None,
    gsd_idx: np.ndarray,
    dpi: int = 150,
) -> list[Path]:
    """Generate PNGs: errors, CDF, scatter, magnitude. Returns list of written paths."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return []
    K = KJMOL_NM_TO_HARTREE_BOHR
    N = len(gsd_idx)
    x = np.arange(N, dtype=np.int32)
    components = [
        ("Bond", "bond", "C0"),
        ("Coulomb", "coulomb", "C1"),
        ("LJ", "lj", "C2"),
        ("Cavity", "cavity", "C3"),
    ]
    written: list[Path] = []
    ymin_log = 1e-20

    # 1. force_component_errors.png: err vs particle index
    fig1, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
    axes = axes.flatten()
    for ax, (label, key, color) in zip(axes, components):
        arr = err_per_particle.get(key, np.full(N, np.nan))
        arr = np.asarray(arr, dtype=np.float64)
        valid = np.isfinite(arr)
        if np.any(valid):
            y_plot = np.maximum(arr, ymin_log)
            ax.plot(x, y_plot, color=color, marker="o", markersize=3, linestyle="-", linewidth=1)
        ax.set_ylabel("|OpenMM − HOOMD| (Ha/Bohr)")
        ax.set_title(label)
        ax.set_yscale("log")
        ax.set_ylim(bottom=ymin_log)
        ax.grid(True, alpha=0.3)
    axes[2].set_xlabel("Particle index (OpenMM order)")
    axes[3].set_xlabel("Particle index (OpenMM order)")
    fig1.suptitle("Force error by component (OpenMM vs cav-hoomd)", fontsize=12)
    for ax in axes:
        ax.set_xlim(x.min() - 0.5, x.max() + 0.5)
    plt.tight_layout()
    p1 = plots_dir / "force_component_errors.png"
    fig1.savefig(p1, dpi=dpi, bbox_inches="tight")
    plt.close(fig1)
    written.append(p1)

    # 2. force_component_errors_cdf.png
    fig2, ax = plt.subplots(1, 1, figsize=(6, 4))
    for label, key, color in components:
        a = np.asarray(err_per_particle.get(key, np.full(N, np.nan)), dtype=np.float64)
        a = a[np.isfinite(a)]
        if a.size == 0:
            continue
        a = np.sort(a)
        cdf = np.arange(1, a.size + 1, dtype=np.float64) / a.size
        ax.plot(a, cdf, label=label, color=color)
    ax.set_xscale("log")
    ax.set_xlabel("Error (Ha/Bohr)")
    ax.set_ylabel("CDF")
    ax.set_title("Distribution of per-particle error by component")
    ax.legend(loc="lower right")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    p2 = plots_dir / "force_component_errors_cdf.png"
    fig2.savefig(p2, dpi=dpi, bbox_inches="tight")
    plt.close(fig2)
    written.append(p2)

    # 3. force_component_scatter.png: F_openmm vs F_hoomd (magnitude) per component
    fig3, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=False, sharey=False)
    axes = axes.flatten()
    for ax, (label, key, color) in zip(axes, components):
        F_o = openmm_components.get(key)
        if F_o is None:
            ax.set_title(label)
            ax.set_visible(True)
            continue
        F_o_au = np.asarray(F_o, dtype=np.float64) * K
        F_h = hoomd_components.get(key) if hoomd_components else None
        F_h_ordered = np.asarray(F_h, dtype=np.float64)[gsd_idx] if (F_h is not None and F_h.shape[0] == N) else np.full((N, 3), np.nan)
        mag_o = np.linalg.norm(F_o_au, axis=1)
        mag_h = np.linalg.norm(F_h_ordered, axis=1)
        valid = np.isfinite(mag_o) & np.isfinite(mag_h)
        if np.any(valid):
            ax.scatter(mag_o[valid], mag_h[valid], c=color, alpha=0.5, s=10)
        lim_max = max(np.nanmax(mag_o), np.nanmax(mag_h), 1e-20)
        ax.plot([0, lim_max], [0, lim_max], "k--", alpha=0.5, label="y=x")
        ax.set_xlabel("|F_openmm| (Ha/Bohr)")
        ax.set_ylabel("|F_hoomd| (Ha/Bohr)")
        ax.set_title(label)
        ax.set_xscale("symlog", linthresh=1e-10)
        ax.set_yscale("symlog", linthresh=1e-10)
        ax.grid(True, alpha=0.3)
    fig3.suptitle("Force magnitude: OpenMM vs cav-hoomd", fontsize=12)
    plt.tight_layout()
    p3 = plots_dir / "force_component_scatter.png"
    fig3.savefig(p3, dpi=dpi, bbox_inches="tight")
    plt.close(fig3)
    written.append(p3)

    # 4. force_component_magnitude.png: |F_openmm| vs |F_hoomd| per particle (same as scatter but linear for small values)
    fig4, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=False, sharey=False)
    axes = axes.flatten()
    for ax, (label, key, color) in zip(axes, components):
        F_o = openmm_components.get(key)
        if F_o is None:
            ax.set_title(label)
            continue
        F_o_au = np.asarray(F_o, dtype=np.float64) * K
        F_h = hoomd_components.get(key) if hoomd_components else None
        F_h_ordered = np.asarray(F_h, dtype=np.float64)[gsd_idx] if (F_h is not None and F_h.shape[0] == N) else np.full((N, 3), np.nan)
        mag_o = np.linalg.norm(F_o_au, axis=1)
        mag_h = np.linalg.norm(F_h_ordered, axis=1)
        valid = np.isfinite(mag_o) & np.isfinite(mag_h)
        if np.any(valid):
            ax.plot(x[valid], mag_o[valid], color=color, alpha=0.7, label="OpenMM", linewidth=1)
            ax.plot(x[valid], mag_h[valid], color=color, alpha=0.5, linestyle="--", label="cav-hoomd", linewidth=1)
        ax.set_xlabel("Particle index")
        ax.set_ylabel("|F| (Ha/Bohr)")
        ax.set_title(label)
        ax.set_yscale("log")
        ax.set_ylim(bottom=ymin_log)
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(True, alpha=0.3)
    fig4.suptitle("Force magnitude by particle (OpenMM vs cav-hoomd)", fontsize=12)
    plt.tight_layout()
    p4 = plots_dir / "force_component_magnitude.png"
    fig4.savefig(p4, dpi=dpi, bbox_inches="tight")
    plt.close(fig4)
    written.append(p4)

    return written


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Debug force mismatch: compare OpenMM vs cav-hoomd per component (bond, Coulomb, LJ, Cavity)."
    )
    parser.add_argument("gsd", type=str, help="GSD file (same config for both codes)")
    parser.add_argument("--lambda", type=float, default=0.001, dest="lambda_coupling")
    parser.add_argument("--cavity-freq", type=float, default=1560.0)
    parser.add_argument("--frame", type=int, default=0)
    parser.add_argument("--ff-dir", type=str, default=None)
    parser.add_argument("--cutoff", type=float, default=1.0, help="OpenMM nonbonded cutoff (nm)")
    parser.add_argument("--summary", type=str, default="force_component_summary.txt", metavar="FILE")
    parser.add_argument(
        "--per-particle",
        type=str,
        default="force_component_errors_per_particle",
        metavar="BASE",
        help="Base name for per-particle diagnostics: writes <BASE>.csv and <BASE>.npz for plotting err vs particle index.",
    )
    parser.add_argument("--raw-forces-txt", type=str, default="force_components_raw.txt", help="TXT with raw F_openmm, F_hoomd per particle per component")
    parser.add_argument("--raw-forces-npz", type=str, default="force_components_raw.npz", help="NPZ with raw force arrays")
    parser.add_argument("--hoomd-components-npz", type=str, default=None, help="Pre-computed cav-hoomd per-component forces (skip cav-hoomd run; use same reference as optimization)")
    parser.add_argument("--no-raw-forces", action="store_true", help="Skip raw force TXT/NPZ output")
    parser.add_argument("--plots-dir", type=str, default=".", help="Directory for PNG outputs")
    parser.add_argument("--no-plots", action="store_true", help="Skip PNG generation")
    parser.add_argument("--dpi", type=int, default=150, help="DPI for PNG outputs")
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
            gsd_path, args.lambda_coupling, args.cavity_freq, ff_dir, args.frame, cutoff_nm=args.cutoff
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

    hoomd_components: dict[str, np.ndarray] | None = None
    bond_raw: np.ndarray | None = None
    pos_hoomd_actual: np.ndarray | None = None

    if args.hoomd_components_npz and Path(args.hoomd_components_npz).exists():
        if not args.quiet:
            print(f"Loading cav-hoomd per-component forces from {args.hoomd_components_npz}...")
        data = np.load(args.hoomd_components_npz)
        hoomd_components = {}
        for key in ["bond", "coulomb", "lj", "cavity"]:
            arr = data.get(f"F_hoomd_{key}")
            if arr is not None:
                hoomd_components[key] = np.asarray(arr, dtype=np.float64)
    else:
        if not args.quiet:
            print("Running cav-hoomd for per-component forces...")
        # OpenMM uses unwrap_molecules() so bond forces are correct. Run HOOMD on the same
        # configuration. HOOMD expects positions in [-L/2, L/2); write temp GSD in that convention
        # so cav-hoomd does not raise "Not all particles were found inside the given box".
        box_bohr = np.asarray(cfg["box_bohr"][:3], dtype=np.float64)
        pos_hoomd = wrap_positions_into_box(pos_unwrapped, box_bohr)  # [0, L)
        pos_hoomd = pos_hoomd - box_bohr / 2.0  # [-L/2, L/2) for HOOMD
        tmp_gsd_path = _SCRIPT_DIR / "_force_debug_temp.gsd"
        try:
            write_gsd_with_positions(gsd_path, args.frame, pos_hoomd, tmp_gsd_path)
            gsd_for_hoomd = tmp_gsd_path
            hoomd_components, bond_raw, pos_hoomd_actual, err = _run_cav_hoomd_and_get_per_component_forces(
                gsd_for_hoomd, args.lambda_coupling, args.cavity_freq, frame=0
            )
        finally:
            if tmp_gsd_path.exists():
                tmp_gsd_path.unlink(missing_ok=True)
        if err is not None:
            if not args.quiet:
                print(f"cav-hoomd per-component not available: {err}", file=sys.stderr)
            hoomd_components = {}
            bond_raw = None
            pos_hoomd_actual = None

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

    # Raw force output (TXT and NPZ)
    raw_txt_path: Path | None = None
    raw_npz_path: Path | None = None
    plot_paths: list[Path] = []
    if not args.no_raw_forces:
        raw_txt_path = Path(args.raw_forces_txt)
        raw_npz_path = Path(args.raw_forces_npz)
        _write_raw_forces_txt(raw_txt_path, openmm_components, hoomd_components, gsd_idx, err_per_particle)
        _write_raw_forces_npz(raw_npz_path, openmm_components, hoomd_components, gsd_idx, err_per_particle)
        if not args.quiet:
            print(f"Raw forces: {raw_txt_path}  {raw_npz_path}")

    # Integrated PNG generation
    if not args.no_plots:
        plots_dir = Path(args.plots_dir)
        plots_dir.mkdir(parents=True, exist_ok=True)
        plot_paths = _generate_plots(
            plots_dir, err_per_particle, openmm_components, hoomd_components, gsd_idx, dpi=args.dpi
        )
        if plot_paths and not args.quiet:
            print(f"Plots: {'  '.join(str(p) for p in plot_paths)}")
        elif not plot_paths and not args.quiet:
            print("Plots skipped (matplotlib not available)", file=sys.stderr)

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
        f.write(f"# GSD: {gsd_path.name}  lambda: {args.lambda_coupling}  cavity_freq: {args.cavity_freq} cm^-1")
        f.write(f"  OpenMM_cutoff: {args.cutoff} nm\n\n")
        if raw_txt_path is not None:
            f.write(f"# Raw forces file: {raw_txt_path}\n")
        if raw_npz_path is not None:
            f.write(f"# Raw forces NPZ: {raw_npz_path}\n")
        if plot_paths:
            f.write("# Plots:\n")
            for p in plot_paths:
                f.write(f"#   {p}\n")
        f.write("\n")
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
        box_bohr_diag = np.asarray(cfg["box_bohr"][:3], dtype=np.float64)
        pos_hoomd_diag = pos_unwrapped - box_bohr_diag / 2.0
        pos_hoomd_diag = pos_hoomd_diag - np.round(pos_hoomd_diag / box_bohr_diag) * box_bohr_diag
        # Use HOOMD state positions (after 1 step) when available so hand F vs F_hoomd is like-with-like
        pos_for_hoomd_cmp = (
            pos_hoomd_actual
            if pos_hoomd_actual is not None and pos_hoomd_actual.shape[0] == N
            else pos_hoomd_diag
        )
        if pos_hoomd_actual is not None and pos_hoomd_actual.shape[0] == N:
            f.write("# Single-bond diagnostic uses HOOMD state positions (after 1 step) for hand vs F_hoomd so geometry matches.\n")
        f.write(_single_bond_diagnostic(cfg, pos_unwrapped, F_bond_openmm, F_bond_hoomd_tag, gsd_idx, bond_indices_diag, pos_hoomd=pos_for_hoomd_cmp))

        # Force magnitude and relative error: compare bond vs LJ apples-to-apples (why is bond abs err ~1e-3 vs LJ ~1e-5?)
        f.write("# --- Force magnitude and relative error (why bond ~1e-3 vs LJ ~1e-5?) ---\n")
        f.write("# Absolute error scales with |F|. Compare relative error rel_err = |Delta F| / max(|F_openmm|, 1e-12).\n")
        K = KJMOL_NM_TO_HARTREE_BOHR
        for comp in ["bond", "coulomb", "lj", "cavity"]:
            Fo = openmm_components.get(comp)
            if Fo is None:
                continue
            F_mag = np.linalg.norm(np.asarray(Fo, dtype=np.float64) * K, axis=1)
            # Only divide where F_mag > threshold to avoid RuntimeWarning (divide by zero / invalid)
            rel_err = np.full_like(F_mag, np.nan, dtype=np.float64)
            mask = F_mag > 1e-12
            np.divide(err_per_particle[comp], F_mag, out=rel_err, where=mask)
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
        f.write("# HOOMD bond force is in local index order (see Correct HOOMD->OpenMM mapping); we use arr[rtag].\n")
        f.write("# Outliers can still appear if: (1) OpenMM config (initial) vs HOOMD config (after 1 step) differ;\n")
        f.write("# (2) hand vs F_hoomd used different positions (we now use HOOMD state positions when available).\n")
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

        # Coulomb: F_openmm vs F_hoomd at worst particles (sign check)
        f.write("# --- Coulomb: F_openmm vs F_hoomd at worst particles (sign check) ---\n")
        if (
            hoomd_components
            and openmm_components.get("coulomb") is not None
            and hoomd_components.get("coulomb") is not None
        ):
            coulomb_worst = _worst_k(err_per_particle["coulomb"], 5)
            F_o_coulomb_au = np.asarray(openmm_components["coulomb"], dtype=np.float64) * KJMOL_NM_TO_HARTREE_BOHR
            F_h_coulomb = hoomd_components["coulomb"]
            F_h_coulomb_ordered = (
                np.asarray(F_h_coulomb)[gsd_idx]
                if F_h_coulomb.shape[0] == N
                else np.full((N, 3), np.nan)
            )
            for idx in (coulomb_worst or [])[:5]:
                fo = np.asarray(F_o_coulomb_au[idx], dtype=np.float64)
                fh = np.asarray(F_h_coulomb_ordered[idx], dtype=np.float64)
                tag = int(gsd_idx[idx])
                f.write(f"# particle_index={idx}  gsd_tag={tag}\n")
                f.write(f"#   F_openmm_coulomb (x,y,z) Ha/Bohr: ({fo[0]:.6e}, {fo[1]:.6e}, {fo[2]:.6e})  |F|={np.linalg.norm(fo):.6e}\n")
                f.write(f"#   F_hoomd_coulomb  (x,y,z) Ha/Bohr: ({fh[0]:.6e}, {fh[1]:.6e}, {fh[2]:.6e})  |F|={np.linalg.norm(fh):.6e}\n")
                d = fo - fh
                f.write(f"#   difference       (x,y,z) Ha/Bohr: ({d[0]:.6e}, {d[1]:.6e}, {d[2]:.6e})  |d|={np.linalg.norm(d):.6e}\n")
                no, nh = np.linalg.norm(fo), np.linalg.norm(fh)
                if no > 1e-12 and nh > 1e-12:
                    ratio = no / nh
                    cos_angle = float(np.dot(fo, fh) / (no * nh))
                    f.write(f"#   ratio |F_openmm|/|F_hoomd| = {ratio:.4f}  ;  cos(angle) = {cos_angle:.4f}\n")
                    f.write("#   -> If cos(angle) ~ -1: forces antiparallel (sign/convention check).\n")
                f.write("\n")
        else:
            f.write("# (Coulomb data missing: skipped.)\n\n")
        f.write("# --- Correct HOOMD -> OpenMM mapping ---\n")
        f.write("# HOOMD force arrays (cpu_local_force_arrays) are in LOCAL particle index order (same as net_force).\n")
        f.write("# ParticleLocalAccessBase: rtag[tag] = local index; so force on particle with tag = arr[rtag[tag]].\n")
        f.write("# Therefore tag-order force: F_tag[tag] = arr[rtag]. We then align to OpenMM: OpenMM particle i has\n")
        f.write("# GSD tag gsd_idx[i], so F_hoomd_in_openmm_order[i] = F_tag[gsd_idx[i]] = arr[rtag][gsd_idx[i]].\n")
        f.write("# HOOMD-blue PotentialBond writes to h_force.data[idx_a], h_force.data[idx_b] with idx = rtag[tag],\n")
        f.write("# so bond force is also in local index order; arr[rtag] is the correct mapping for bond as well.\n")
        f.write("# Single-bond diagnostic: when HOOMD state positions are available, we use them for hand F vs F_hoomd\n")
        f.write("# so the comparison is like-with-like (same geometry HOOMD used; cav-hoomd runs 1 step before we read).\n")
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
