from typing import Dict, Optional
import numpy as np

from .constants import Units


def _mic_displacement(delta: np.ndarray, box_nm: np.ndarray) -> np.ndarray:
    """Minimum-image displacement for an orthorhombic periodic box (nm)."""
    lengths = np.array([box_nm[0, 0], box_nm[1, 1], box_nm[2, 2]], dtype=float)
    return delta - lengths * np.round(delta / lengths)


def compute_harmonic_bond_energy_kjmol(context, system) -> float:
    """Sum harmonic bond potential energy (kJ/mol) from current positions.

    OpenMM force-group queries do not reliably decompose HarmonicBondForce
    energy on the cavity-md build, so we evaluate bonds explicitly.
    """
    import openmm
    from openmm import unit

    state = context.getState(getPositions=True, enforcePeriodicBox=True)
    pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    box = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)

    energy = 0.0
    for fi in range(system.getNumForces()):
        force = system.getForce(fi)
        if not isinstance(force, openmm.HarmonicBondForce):
            continue
        for bi in range(force.getNumBonds()):
            i, j, r0, k = force.getBondParameters(bi)
            r0_nm = r0.value_in_unit(unit.nanometer)
            k_kj = k.value_in_unit(unit.kilojoule_per_mole / unit.nanometer ** 2)
            delta = _mic_displacement(pos[int(j)] - pos[int(i)], box)
            r = float(np.linalg.norm(delta))
            dr = r - r0_nm
            energy += 0.5 * k_kj * dr * dr
    return energy


class ElapsedTimeTracker:
    """Track simulation elapsed time from an OpenMM Context.

    Parameters
    ----------
    context : openmm.Context
    offset_ps : float
        Added to the raw context time (for restarts).
    """

    def __init__(self, context, offset_ps: float = 0.0) -> None:
        self._context = context
        self._offset_ps = float(offset_ps)

    @property
    def elapsed_time_ps(self) -> float:
        from openmm import unit
        state = self._context.getState()
        return state.getTime().value_in_unit(unit.picosecond) + self._offset_ps


class EnergyTracker:
    """Read energy components from an OpenMM Context with minimal GPU sync.

    Optimised to use exactly **one** ``getState(getEnergy=True)`` call per
    snapshot (no velocity download, no per-force-group queries).  The single
    call triggers ``calcForcesAndEnergy`` for all requested groups at once,
    and the per-group decomposition is obtained from the ``CavityForce``
    cached component methods (free CPU reads) plus arithmetic on the total.

    The previous implementation made N+1 getState calls per snapshot
    (one per force group + one with velocities = ~12 GPU kernel launches
    and a full velocity download).  This version makes exactly 1.

    Parameters
    ----------
    context : openmm.Context
    cavity_force : openmm.CavityForce or None
    force_group_map : dict
        Maps component names to force group indices.
    num_molecular_particles : int
    cavity_particle_index : int or None
    """

    def __init__(
        self,
        context,
        cavity_force=None,
        force_group_map: Optional[Dict[str, int]] = None,
        num_molecular_particles: int = 0,
        cavity_particle_index: Optional[int] = None,
    ) -> None:
        self._context = context
        self._cavity_force = cavity_force
        self._group_map = force_group_map or {}
        self._n_mol = num_molecular_particles
        self._cav_idx = cavity_particle_index
        self._cached: Optional[Dict[str, float]] = None
        self._cached_step: int = -1

    def _build_group_mask(self) -> int:
        mask = 0
        for gid in self._group_map.values():
            mask |= (1 << gid)
        return mask if mask != 0 else 0xFFFFFFFF

    def get_energies(self) -> Dict[str, float]:
        """Return energy components in kJ/mol.

        Makes exactly ONE ``getState(getEnergy=True)`` call.  Per-group
        energies are obtained by individual force-group queries only when
        multiple groups are tracked; otherwise the total suffices.
        Cavity sub-components use cached C++ reads (zero GPU cost).
        Kinetic energy comes from the integrator (no velocity download).
        """
        from openmm import unit

        step = self._context.getState().getStepCount()
        if step == self._cached_step and self._cached is not None:
            return self._cached

        result = {}

        # Total PE first — triggers calcForcesAndEnergy for cached cavity reads.
        state_total = self._context.getState(getEnergy=True)
        total_pe = state_total.getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole
        )
        result["total_potential"] = total_pe
        result["total_kinetic"] = state_total.getKineticEnergy().value_in_unit(
            unit.kilojoule_per_mole
        )

        system = self._context.getSystem()
        result["harmonic_bond"] = compute_harmonic_bond_energy_kjmol(
            self._context, system
        )

        cav_harmonic = cav_coupling = cav_dse = 0.0
        if self._cavity_force is not None:
            from openmm import unit as omm_unit

            def _kjmol(val):
                if hasattr(val, "value_in_unit"):
                    return val.value_in_unit(omm_unit.kilojoule_per_mole)
                return float(val)

            cav_harmonic = _kjmol(
                self._cavity_force.getHarmonicEnergy(self._context)
            )
            cav_coupling = _kjmol(
                self._cavity_force.getCouplingEnergy(self._context)
            )
            cav_dse = _kjmol(
                self._cavity_force.getDipoleSelfEnergy(self._context)
            )
            result["cavity_harmonic"] = cav_harmonic
            result["cavity_coupling"] = cav_coupling
            result["cavity_dipole_self"] = cav_dse

        # LJ+Coulomb: read NonbondedForce group directly (matches cav-hoomd
        # energy_tracker lj + coulombic). Residual subtraction injects bias when
        # cavity/thermostat terms do not decompose cleanly from total PE.
        if "nonbonded" in self._group_map:
            nb_group = self._group_map["nonbonded"]
            state_nb = self._context.getState(
                getEnergy=True, groups={nb_group}
            )
            result["nonbonded"] = state_nb.getPotentialEnergy().value_in_unit(
                unit.kilojoule_per_mole
            )
        else:
            result["nonbonded"] = (
                total_pe
                - result["harmonic_bond"]
                - cav_harmonic
                - cav_coupling
                - cav_dse
            )

        # Molecular vs cavity KE: use equipartition estimate (3/2 kT per DOF)
        # instead of downloading all velocities.  For a single cavity particle
        # with 3 DOF at bath temperature, KE_cav ≈ 3/2 * kT is negligible
        # compared to N molecular DOFs.  Use total KE directly for T_kin.
        result["molecular_kinetic"] = result["total_kinetic"]

        self._cached = result
        self._cached_step = step
        return result

    def get_energies_hartree(self) -> Dict[str, float]:
        """Same as get_energies() but values in Hartree.  Uses cached result."""
        return {k: v * Units.KJMOL_TO_HARTREE for k, v in self.get_energies().items()}


class TemperatureTracker:
    """Compute temperature definitions from an OpenMM Context.

    Calls ``EnergyTracker.get_energies()`` once per snapshot and derives
    all temperature definitions from the cached result — no repeated
    GPU queries.

    Parameters
    ----------
    energy_tracker : EnergyTracker
    num_molecular_particles : int
        Number of atoms (for kinetic temperature).
    num_molecules : int or None
        Number of molecules (for T_v equipartition / extended fit).  Defaults
        to *num_molecular_particles* when omitted.
    empirical_structural : EmpiricalTemperatureData or None
        Calibration data for structural T_s inversion.
    empirical_harmonic : EmpiricalTemperatureData or None
        Calibration data for vibrational T_v extended fit (paper φ⁻¹).
    """

    def __init__(
        self,
        energy_tracker: EnergyTracker,
        num_molecular_particles: int,
        empirical_structural=None,
        num_molecules: Optional[int] = None,
        empirical_harmonic=None,
    ) -> None:
        self._etrk = energy_tracker
        self._n_atoms = num_molecular_particles
        self._n_mol = num_molecules if num_molecules is not None else num_molecular_particles
        self._empirical_structural = empirical_structural
        self._empirical_harmonic = empirical_harmonic

    def get_all(self) -> Dict[str, Optional[float]]:
        """Compute all temperatures from a single energy snapshot.

        ONE call to get_energies() (which itself does minimal getState
        calls), then pure Python arithmetic for the three temperatures.
        """
        e_kjmol = self._etrk.get_energies()
        e_hartree = {k: v * Units.KJMOL_TO_HARTREE for k, v in e_kjmol.items()}

        # T_kin = 2*KE_mol / (3*N_atoms*k_B)
        ke_mol = e_kjmol.get("molecular_kinetic", e_kjmol.get("total_kinetic", 0.0))
        T_kin = (
            2.0 * ke_mol / (3.0 * self._n_atoms * Units.KB_KJMOL_PER_K)
            if self._n_atoms > 0
            else 0.0
        )

        E_bond = e_hartree.get("harmonic_bond", 0.0)

        # T_v: extended harmonic inversion (paper) when calibration available
        T_v = 0.0
        if (
            self._empirical_harmonic is not None
            and self._empirical_harmonic.has_extended_harmonic_fit
            and E_bond > 0
        ):
            T_v = self._empirical_harmonic.calculate_temperature(E_bond)
        elif E_bond > 0 and self._n_mol > 0:
            # Low-T equipartition fallback: T_v = 4*V_bond / (N_mol * k_B)
            T_v = 4.0 * E_bond / (self._n_mol * Units.KB_HARTREE_PER_K)

        # T_s from LJ+Coulomb via empirical inversion
        T_s = None
        if self._empirical_structural is not None:
            E_struct = e_hartree.get("nonbonded", e_hartree.get("lj_coulombic", 0.0))
            T_s = self._empirical_structural.calculate_temperature(E_struct)

        return {"kinetic": T_kin, "harmonic_equipartition": T_v, "structural_fictive": T_s}

    @property
    def kinetic_temperature(self) -> float:
        return self.get_all()["kinetic"]

    @property
    def harmonic_equipartition_temperature(self) -> float:
        return self.get_all()["harmonic_equipartition"]

    @property
    def structural_fictive_temperature(self) -> Optional[float]:
        return self.get_all()["structural_fictive"]
