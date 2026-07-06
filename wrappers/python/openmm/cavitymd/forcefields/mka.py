"""Modified Kob–Andersen (mKA) dimer force field (A–A / B–B mixture)."""

from __future__ import annotations

import numpy as np
import openmm
from openmm import unit

from openmm.cavitymd.constants import Units

# Masses (amu) — A ≈ O, B ≈ N in the cav-hoomd mKA model
MASS_A = 16.0
MASS_B = 14.0
PHOTON_MASS_AMU = 1.0

# Harmonic bonds (atomic units)
K_AA_AU = 0.73204
R0_AA_AU = 2.281655158
K_BB_AU = 1.4325
R0_BB_AU = 2.0743522177

# Lennard-Jones (atomic units)
EPS_AA_AU = 0.00016685201
SIG_AA_AU = 6.230426584
EPS_BB_AU = 0.000083426
SIG_BB_AU = 5.48277488
EPS_AB_AU = 0.00025027802
SIG_AB_AU = 4.9832074319

CHARGE_MAG = 0.3
OMEGA_C_CM1 = 1560.0

_REF_N = 250
_REF_L_AU = 40.0
_RCUT_AU = 15.0
_FRACTION_AA = 0.8

B2NM = Units.BOHR_TO_NM
H2K = Units.HARTREE_TO_KJMOL


def box_au_for_num_molecules(num_molecules: int) -> float:
    """Cubic box edge (Bohr) at constant density (reference N=250, L=40 Bohr)."""
    return _REF_L_AU * (num_molecules / _REF_N) ** (1.0 / 3.0)


def _to_omm_bond(k_au: float, r0_au: float) -> tuple[float, float]:
    return k_au * H2K / B2NM**2, r0_au * B2NM


def _to_omm_lj(eps_au: float, sig_au: float) -> tuple[float, float]:
    return eps_au * H2K, sig_au * B2NM


def build_mka_system(
    num_molecules: int,
    seed: int = 42,
    sample_bonds_at_T: float | None = None,
    fraction_AA: float = _FRACTION_AA,
) -> tuple[openmm.System, list, int]:
    """Build periodic mKA dimer box with PME Coulomb and KA LJ table."""
    rng = np.random.default_rng(seed)
    box_nm = box_au_for_num_molecules(num_molecules) * B2NM
    if box_nm < 2.0 * _RCUT_AU * B2NM:
        raise ValueError(
            f"Box too small for PME cutoff ({_RCUT_AU} Bohr): "
            f"need N large enough (>=106 typical), got edge {box_nm / B2NM:.1f} Bohr"
        )

    k_aa, r0_aa = _to_omm_bond(K_AA_AU, R0_AA_AU)
    k_bb, r0_bb = _to_omm_bond(K_BB_AU, R0_BB_AU)
    eps_aa, sig_aa = _to_omm_lj(EPS_AA_AU, SIG_AA_AU)
    eps_bb, sig_bb = _to_omm_lj(EPS_BB_AU, SIG_BB_AU)
    eps_ab, sig_ab = _to_omm_lj(EPS_AB_AU, SIG_AB_AU)

    system = openmm.System()
    bond_force = openmm.HarmonicBondForce()
    nonbonded = openmm.NonbondedForce()
    nonbonded.setNonbondedMethod(openmm.NonbondedForce.PME)
    nonbonded.setCutoffDistance(_RCUT_AU * B2NM)

    num_aa = int(fraction_AA * num_molecules)
    side = int(np.ceil(num_molecules ** (1.0 / 3.0)))
    spacing = box_nm / side
    positions: list = []
    a_indices: list[tuple[int, float]] = []
    b_indices: list[tuple[int, float]] = []

    mol = 0
    for i in range(side):
        for j in range(side):
            for k in range(side):
                if mol >= num_molecules:
                    break
                is_aa = mol < num_aa
                mass = MASS_A if is_aa else MASS_B
                r0 = r0_aa if is_aa else r0_bb
                k_bond = k_aa if is_aa else k_bb
                sigma = sig_aa if is_aa else sig_bb
                epsilon = eps_aa if is_aa else eps_bb

                center = np.array([(i + 0.5) * spacing, (j + 0.5) * spacing, (k + 0.5) * spacing])
                theta = rng.random() * 2.0 * np.pi
                phi = np.arccos(2.0 * rng.random() - 1.0)
                direction = np.array([
                    np.sin(phi) * np.cos(theta),
                    np.sin(phi) * np.sin(theta),
                    np.cos(phi),
                ])
                if sample_bonds_at_T is not None and sample_bonds_at_T > 0:
                    kT = Units.kelvin_to_kT_kjmol(sample_bonds_at_T)
                    r0_eff = float(rng.normal(r0, np.sqrt(kT / k_bond)))
                else:
                    r0_eff = r0

                r1 = center - 0.5 * r0_eff * direction
                r2 = center + 0.5 * r0_eff * direction
                idx1 = system.addParticle(mass)
                idx2 = system.addParticle(mass)
                positions.append(openmm.Vec3(*r1) * unit.nanometer)
                positions.append(openmm.Vec3(*r2) * unit.nanometer)
                bond_force.addBond(idx1, idx2, r0, k_bond)
                q1, q2 = -CHARGE_MAG, +CHARGE_MAG
                nonbonded.addParticle(q1, sigma, epsilon)
                nonbonded.addParticle(q2, sigma, epsilon)
                nonbonded.addException(idx1, idx2, 0.0, 1.0, 0.0)
                if is_aa:
                    a_indices.extend([(idx1, q1), (idx2, q2)])
                else:
                    b_indices.extend([(idx1, q1), (idx2, q2)])
                mol += 1

    for idx_a, qa in a_indices:
        for idx_b, qb in b_indices:
            nonbonded.addException(idx_a, idx_b, qa * qb, sig_ab, eps_ab)

    system.addForce(bond_force)
    system.addForce(nonbonded)
    system.setDefaultPeriodicBoxVectors(
        openmm.Vec3(box_nm, 0, 0),
        openmm.Vec3(0, box_nm, 0),
        openmm.Vec3(0, 0, box_nm),
    )
    return system, positions, system.getNumParticles()


def add_cavity_particle(system: openmm.System, positions: list) -> int:
    """Append photon; register q=0 ghost in NonbondedForce."""
    idx = system.addParticle(PHOTON_MASS_AMU)
    positions.append(openmm.Vec3(0, 0, 0) * unit.nanometer)
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            force.addParticle(0.0, 0.1, 0.0)
    return idx
