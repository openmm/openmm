"""Shared helpers for the mKA cavity MD tutorial and validation scripts."""

from __future__ import annotations

import numpy as np
import openmm
from openmm import unit

from openmm.cavitymd.constants import Units
from openmm.cavitymd.forcefields.mka import (
    CHARGE_MAG,
    K_AA_AU,
    MASS_A,
    OMEGA_C_CM1,
    PHOTON_MASS_AMU,
    R0_AA_AU,
)
from openmm.cavitymd.thermostats import (
    DEFAULT_LANGEVIN_FRICTION_PS,
    create_langevin_integrator,
)

B2NM = Units.BOHR_TO_NM
H2K = Units.HARTREE_TO_KJMOL
K_AA_OMM = K_AA_AU * H2K / B2NM**2
R0_AA_OMM = R0_AA_AU * B2NM


def select_platform(prefer_cuda: bool = True) -> openmm.Platform:
    """Return CUDA (mixed precision) when available, otherwise CPU/Reference."""
    if prefer_cuda:
        try:
            platform = openmm.Platform.getPlatformByName("CUDA")
            platform.setPropertyDefaultValue("Precision", "mixed")
            return platform
        except Exception:
            pass
    return openmm.Platform.getPlatformByName("CPU")


def build_single_aa_dimer_system(
    lambda_coupling: float,
    omegac_au: float | None = None,
) -> tuple[openmm.System, openmm.CavityParticleDisplacer, list]:
    """Build a single A-A dimer + photon system (no integrator/thermostat)."""
    if omegac_au is None:
        omegac_au = Units.cm1_to_au(OMEGA_C_CM1)

    system = openmm.System()
    system.addParticle(MASS_A)
    system.addParticle(MASS_A)
    system.addParticle(PHOTON_MASS_AMU)

    bond_force = openmm.HarmonicBondForce()
    bond_force.addBond(0, 1, R0_AA_OMM, K_AA_OMM)
    system.addForce(bond_force)

    cavity_idx = 2
    cavity_force = openmm.CavityForce(
        cavity_idx, omegac_au, lambda_coupling, PHOTON_MASS_AMU
    )
    system.addForce(cavity_force)

    displacer = openmm.CavityParticleDisplacer(
        cavity_idx, omegac_au, PHOTON_MASS_AMU
    )
    displacer.setSwitchOnStep(2**31 - 1)
    system.addForce(displacer)

    half_r0 = R0_AA_OMM / 2.0
    positions = [
        openmm.Vec3(-half_r0, 0, 0) * unit.nanometer,
        openmm.Vec3(+half_r0, 0, 0) * unit.nanometer,
        openmm.Vec3(0, 0, 0) * unit.nanometer,
    ]
    return system, displacer, positions


def create_context(
    system: openmm.System,
    dt_fs: float,
    temperature_K: float,
    seed: int,
    platform_name: str | None = None,
    *,
    use_langevin: bool = False,
    friction_ps_inv: float = DEFAULT_LANGEVIN_FRICTION_PS,
) -> openmm.Context:
    """Create a Context with thermal velocities on the selected platform."""
    dt_ps = dt_fs * 1e-3
    if use_langevin:
        integrator = create_langevin_integrator(
            temperature_K, dt_ps, friction_ps_inv=friction_ps_inv
        )
    else:
        integrator = openmm.VerletIntegrator(dt_ps)
    if platform_name is None:
        platform = select_platform()
    else:
        platform = openmm.Platform.getPlatformByName(platform_name)
    context = openmm.Context(system, integrator, platform)
    context.setVelocitiesToTemperature(temperature_K, seed)
    return context


def molecular_kinetic_energy(
    state: openmm.State,
    molecular_indices: list[int],
    system: openmm.System,
) -> float:
    """Kinetic energy of selected particles in kJ/mol."""
    vel = state.getVelocities(asNumpy=True).value_in_unit(
        unit.nanometer / unit.picosecond
    )
    ke = 0.0
    for idx in molecular_indices:
        mass = system.getParticleMass(idx).value_in_unit(unit.dalton)
        ke += 0.5 * mass * float(np.sum(vel[idx] ** 2))
    return ke


def molecular_kinetic_temperature(
    state: openmm.State,
    system: openmm.System,
    molecular_indices: list[int] | None = None,
    subtract_com: bool = False,
) -> float:
    """Molecular kinetic temperature from translational DOFs.

    Use subtract_com=False (3N DOF) with LangevinMiddleIntegrator, which
    thermostats each particle independently.  Use subtract_com=True (3N-3)
    when a BussiThermostat with setSubtractCMMotion(True) enforces the bath.
    """
    if molecular_indices is None:
        molecular_indices = [0, 1]
    ke = molecular_kinetic_energy(state, molecular_indices, system)
    dof = 3 * len(molecular_indices)
    if subtract_com:
        dof -= 3
    return 2.0 * ke / (dof * Units.KB_KJMOL_PER_K)


def system_kinetic_temperature(
    state: openmm.State,
    system: openmm.System,
) -> float:
    """Total kinetic temperature (3 DOF per particle) from an OpenMM State."""
    ke = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
    dof = 3 * system.getNumParticles()
    return 2.0 * ke / (dof * Units.KB_KJMOL_PER_K)


def photon_kinetic_temperature(
    state: openmm.State,
    system: openmm.System,
    photon_index: int = 2,
    dof: int = 2,
    in_plane_components: tuple[int, int] = (1, 2),
) -> float:
    """Photon kinetic temperature from selected translational DOFs.

    With LangevinMiddleIntegrator all particles (including the photon) are
    thermostatted at T_bath; use dof=3 for the canonical kinetic temperature.
    Use dof=2 for the cavity plane (y, z when the dimer lies along x) as a
    supplementary in-plane diagnostic.
    """
    vel = state.getVelocities(asNumpy=True).value_in_unit(
        unit.nanometer / unit.picosecond
    )
    mass = system.getParticleMass(photon_index).value_in_unit(unit.dalton)
    if dof == 2:
        ke = 0.5 * mass * float(
            vel[photon_index, in_plane_components[0]] ** 2
            + vel[photon_index, in_plane_components[1]] ** 2
        )
    else:
        ke = 0.5 * mass * float(np.sum(vel[photon_index] ** 2))
    return 2.0 * ke / (dof * Units.KB_KJMOL_PER_K)


def collect_dipole_trajectory(
    context: openmm.Context,
    integrator: openmm.Integrator,
    n_steps: int,
    sample_stride: int = 1,
) -> np.ndarray:
    """Collect molecular dipole d(t) for the single A-A dimer."""
    charges = np.array([-CHARGE_MAG, +CHARGE_MAG])
    dipoles = []
    for step in range(n_steps):
        integrator.step(1)
        if step % sample_stride != 0:
            continue
        state = context.getState(getPositions=True)
        pos = state.getPositions(asNumpy=True)
        dipoles.append(charges[0] * pos[0] + charges[1] * pos[1])
    return np.asarray(dipoles)


def dipole_spectrum_cm1(
    dipoles: np.ndarray,
    dt_fs: float,
    component: int = 0,
    min_freq_cm1: float = 500.0,
    max_freq_cm1: float = 3000.0,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Compute dipole power spectrum (direct FFT) and return peak frequency in cm^-1."""
    d_x = dipoles[:, component].copy()
    n_sig = len(d_x)
    d_x -= d_x.mean()

    dt_s = dt_fs * 1e-15
    freqs_hz = np.fft.rfftfreq(n_sig, d=dt_s)
    freqs_cm1 = freqs_hz / 3e10
    spectrum = np.abs(np.fft.rfft(d_x)) ** 2

    mask = (freqs_cm1 > min_freq_cm1) & (freqs_cm1 < max_freq_cm1)
    peak_idx = int(np.argmax(spectrum[mask]))
    peak_cm1 = float(freqs_cm1[mask][peak_idx])
    return freqs_cm1, spectrum, peak_cm1


def run_nvt_single_dimer(
    lambda_coupling: float = 0.01,
    temperature_K: float = 100.0,
    dt_fs: float = 1.0,
    n_steps: int = 5000,
    seed: int = 42,
    sample_stride: int = 1,
    platform_name: str | None = None,
    friction_ps_inv: float = DEFAULT_LANGEVIN_FRICTION_PS,
) -> dict:
    """Run tutorial Section 2 (Langevin NVT) and return diagnostics."""
    omegac_au = Units.cm1_to_au(OMEGA_C_CM1)
    system, displacer, positions = build_single_aa_dimer_system(
        lambda_coupling=lambda_coupling,
        omegac_au=omegac_au,
    )

    dt_ps = dt_fs * 1e-3
    integrator = create_langevin_integrator(
        temperature_K, dt_ps, friction_ps_inv=friction_ps_inv
    )
    platform = (
        select_platform()
        if platform_name is None
        else openmm.Platform.getPlatformByName(platform_name)
    )
    context = openmm.Context(system, integrator, platform)
    context.setPositions(positions)
    context.setVelocitiesToTemperature(temperature_K, seed)
    displacer.displaceToEquilibrium(context, lambda_coupling)

    system_temperatures = []
    molecular_temperatures = []
    photon_temperatures = []
    dipoles = []
    charges = np.array([-CHARGE_MAG, +CHARGE_MAG])

    for step in range(n_steps):
        integrator.step(1)
        if step % sample_stride != 0:
            continue
        state = context.getState(
            getPositions=True, getVelocities=True, getEnergy=True
        )
        pos = state.getPositions(asNumpy=True)
        dipoles.append(charges[0] * pos[0] + charges[1] * pos[1])
        system_temperatures.append(system_kinetic_temperature(state, system))
        molecular_temperatures.append(
            molecular_kinetic_temperature(state, system, subtract_com=False)
        )
        photon_temperatures.append(
            photon_kinetic_temperature(state, system, dof=3)
        )

    dipoles_arr = np.asarray(dipoles)
    freqs, spectrum, peak_cm1 = dipole_spectrum_cm1(dipoles_arr, dt_fs)

    return {
        "mean_system_temperature_K": float(np.mean(system_temperatures)),
        "mean_temperature_K": float(np.mean(molecular_temperatures)),
        "mean_photon_temperature_K": float(np.mean(photon_temperatures)),
        "peak_frequency_cm1": peak_cm1,
        "omega_c_cm1": OMEGA_C_CM1,
        "dipoles": dipoles_arr,
        "system_temperatures": np.asarray(system_temperatures),
        "temperatures": np.asarray(molecular_temperatures),
        "photon_temperatures": np.asarray(photon_temperatures),
        "freqs_cm1": freqs,
        "spectrum": spectrum,
    }
