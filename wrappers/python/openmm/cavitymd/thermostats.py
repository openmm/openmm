from typing import List, Optional
import numpy as np
import math

from .constants import Units


class DualThermostat:
    """Separate thermostat control for molecular and cavity particles.

    Molecular particles: BussiThermostat (CSVR) with dynamic temperature via
    context.setParameter("BussiTemperature", T_kelvin).

    Cavity particle: NVE by default (equilibrates via coupling). Optional
    Langevin-style friction applied as post-hoc velocity rescaling.

    Parameters
    ----------
    context : openmm.Context
    system : openmm.System
    cavity_particle_index : int or None
    cavity_friction_ps_inv : float
        Langevin friction for cavity particle (1/ps). 0 = NVE.
    cavity_temperature_K : float
        Initial cavity bath temperature (only if friction > 0).
    """

    def __init__(
        self,
        context,
        system,
        cavity_particle_index: Optional[int] = None,
        cavity_friction_ps_inv: float = 0.0,
        cavity_temperature_K: float = 100.0,
    ) -> None:
        self._context = context
        self._system = system
        self._cav_idx = cavity_particle_index
        self.cavity_friction = float(cavity_friction_ps_inv)
        self.cavity_temperature_K = float(cavity_temperature_K)

    def set_molecular_temperature(self, T_kelvin: float) -> None:
        """Update Bussi thermostat target temperature (K)."""
        self._context.setParameter("BussiTemperature", T_kelvin)

    def get_molecular_temperature(self) -> float:
        """Read current Bussi thermostat temperature (K)."""
        return self._context.getParameter("BussiTemperature")

    def set_cavity_temperature(self, T_kelvin: float) -> None:
        self.cavity_temperature_K = float(T_kelvin)

    def apply_cavity_thermostat_step(self, dt_ps: float) -> None:
        """Apply one Langevin velocity-rescaling step to the cavity particle.

        v_new = c1*v_old + c2*R, where c1 = exp(-gamma*dt),
        c2 = sqrt(kT/m * (1 - c1^2)), R ~ N(0,1).
        """
        if self._cav_idx is None or self.cavity_friction <= 0.0:
            return

        from openmm import unit

        state = self._context.getState(getVelocities=True)
        vel_all = state.getVelocities(asNumpy=True).value_in_unit(
            unit.nanometer / unit.picosecond
        )

        cav_mass_amu = self._system.getParticleMass(self._cav_idx).value_in_unit(unit.dalton)

        gamma = self.cavity_friction
        c1 = math.exp(-gamma * dt_ps)
        kT_kjmol = Units.kelvin_to_kT_kjmol(self.cavity_temperature_K)
        # kT/m in (nm/ps)^2: kT [kJ/mol] / m [amu] * (1 amu*nm^2/ps^2 per kJ/mol)
        # 1 kJ/mol = 1 amu * nm^2 / ps^2  (OpenMM internal unit consistency)
        kT_over_m = kT_kjmol / cav_mass_amu
        c2 = math.sqrt(kT_over_m * (1.0 - c1 * c1))

        v_old = np.array(vel_all[self._cav_idx], dtype=np.float64)
        v_new = c1 * v_old + c2 * np.random.randn(3)

        vel_all[self._cav_idx] = v_new
        self._context.setVelocities(vel_all)

    @staticmethod
    def setup_bussi_for_system(
        system,
        molecular_indices: List[int],
        temperature_K: float,
        tau_ps: float = 1.0,
        random_number_seed: int = 0,
    ) -> None:
        """Add BussiThermostat to a System for molecular particles only.

        Must be called before Context creation.
        """
        import openmm

        bussi = openmm.BussiThermostat(temperature_K, tau_ps)
        bussi.setApplyToAllParticles(False)
        bussi.setSubtractCMMotion(True)  # DOF = 3N-3, matches cav-hoomd
        bussi.setRandomNumberSeed(int(random_number_seed))
        for idx in molecular_indices:
            bussi.addParticle(idx)
        system.addForce(bussi)
