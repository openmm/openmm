from typing import Optional, Dict, TextIO
import time as wall_time
from pathlib import Path

from .constants import Units
from .trackers import ElapsedTimeTracker, EnergyTracker, TemperatureTracker
from .thermostats import DualThermostat
from .variants import CouplingVariant
from .feedback import EmpiricalTemperatureFeedback, GradientDescentFeedback
from .controllers import DiffEqController, SimpleSetpointController, PIDControl


def setup_gpu_square_wave(cavity_force, amplitude: float, period_ps: float,
                          duty_cycle: float = 0.5, start_time_ps: float = 0.0,
                          stop_time_ps: float = -1.0) -> None:
    """Configure GPU-side square-wave modulation on a CavityForce.

    Once set, the CUDA/OpenCL kernel evaluates lambda from time_ps each
    step with zero host round-trips.  Call updateParametersInContext()
    once after this to push the modulation params to the device.

    Parameters
    ----------
    cavity_force : openmm.CavityForce
    amplitude : float
        Lambda value when coupling is "on".
    period_ps : float
        Full period (ps).
    duty_cycle : float
        Fraction of period that is "on" [0,1]. Default 0.5.
    start_time_ps : float
        When modulation activates (ps). Default 0.
    stop_time_ps : float
        When modulation deactivates (-1 = never). Default -1.
    """
    import openmm
    cavity_force.setCouplingModulation(
        openmm.CavityForce.ModulationSquareWave,
        amplitude, period_ps, duty_cycle, start_time_ps, stop_time_ps, 1.0,
    )


def setup_gpu_step(cavity_force, amplitude: float, start_time_ps: float,
                   stop_time_ps: float = -1.0) -> None:
    """Configure GPU-side step-function modulation."""
    import openmm
    cavity_force.setCouplingModulation(
        openmm.CavityForce.ModulationStep,
        amplitude, 0.0, 0.5, start_time_ps, stop_time_ps, 1.0,
    )


def setup_gpu_decaying_step(cavity_force, amplitude: float,
                            start_time_ps: float, decay_tau_ps: float,
                            stop_time_ps: float = -1.0) -> None:
    """Configure GPU-side decaying-step modulation."""
    import openmm
    cavity_force.setCouplingModulation(
        openmm.CavityForce.ModulationDecayingStep,
        amplitude, 0.0, 0.5, start_time_ps, stop_time_ps, decay_tau_ps,
    )


def setup_gpu_adaptive_square_wave(
    cavity_force,
    target_coupling: float,
    target_temperature_K: float,
    period_ps: float,
    duty_cycle: float = 0.5,
    start_time_ps: float = 0.0,
    stop_time_ps: float = -1.0,
    min_amplitude: float = 1e-8,
    max_amplitude: float = 0.1,
) -> None:
    """Configure GPU-side adaptive square-wave modulation.

    The kernel updates amplitude once per period on the GPU:
    ``g_next = target_coupling * sqrt(target_temperature_K / T_bath)``
    where T_bath is read from the BussiThermostat context parameter.
    No Python round-trip needed.

    Call ``cavity_force.updateParametersInContext(context)`` once after
    this to push the params to the device.
    """
    cavity_force.setAdaptiveSquareWaveModulation(
        target_coupling, target_temperature_K,
        period_ps, duty_cycle, start_time_ps, stop_time_ps,
        min_amplitude, max_amplitude,
    )


def setup_gpu_decaying_square_wave(
    cavity_force,
    initial_amplitude: float,
    period_ps: float,
    duty_cycle: float = 0.5,
    decay_rate_per_period: float = 0.0,
    start_time_ps: float = 0.0,
    stop_time_ps: float = -1.0,
    minimum_amplitude: float = 1e-8,
) -> None:
    """Configure GPU-side decaying square-wave modulation."""
    cavity_force.setDecayingSquareWaveModulation(
        initial_amplitude, period_ps, duty_cycle, decay_rate_per_period,
        start_time_ps, stop_time_ps, minimum_amplitude,
    )


def setup_gpu_sinusoid(
    cavity_force,
    amplitude: float,
    period_ps: float,
    phase_offset: float = 0.0,
    start_time_ps: float = 0.0,
    stop_time_ps: float = -1.0,
) -> None:
    """Configure GPU-side sinusoidal coupling modulation."""
    cavity_force.setSinusoidModulation(
        amplitude, period_ps, phase_offset, start_time_ps, stop_time_ps,
    )


def setup_gpu_exponential_wave(
    cavity_force,
    amplitude: float,
    period_ps: float,
    decay_tau_ps: float,
    start_time_ps: float = 0.0,
    stop_time_ps: float = -1.0,
) -> None:
    """Configure GPU-side exponential-wave coupling modulation."""
    cavity_force.setExponentialWaveModulation(
        amplitude, period_ps, decay_tau_ps, start_time_ps, stop_time_ps,
    )


def setup_multimode_adaptive_square_wave(
    force,
    period_ps: float,
    duty_cycle: float,
    mode_params,
    start_time_ps: float = 0.0,
    stop_time_ps: float = -1.0,
) -> None:
    """Configure per-mode adaptive square-wave on a MultiModeCavityForce.

    All modes share timing (period, duty cycle, start/stop). Each mode
    adapts its amplitude independently on the GPU:
    ``g_next_n = g_target_n * sqrt(T_target_n / T_bath)``

    Parameters
    ----------
    force : openmm.MultiModeCavityForce
    period_ps : float
        Square-wave period (ps), shared by all modes.
    duty_cycle : float
        Fraction ON [0,1], shared by all modes.
    mode_params : list of (g_target, T_target_K, min_amp, max_amp) tuples
        One entry per mode, in order (mode 0, mode 1, ...).
    start_time_ps, stop_time_ps : float
        Activation window.
    """
    force.setAdaptiveSquareWaveModulation(period_ps, duty_cycle,
                                          start_time_ps, stop_time_ps)
    for i, (g, t, mn, mx) in enumerate(mode_params):
        force.setModeModulationParams(i, g, t, mn, mx)


def configure_dipole_self_energy(force, include: bool = True) -> None:
    """Enable or disable the dipole self-energy (self-polarization) term.

    Works with ``CavityForce`` and ``MultiModeCavityForce``. When disabled,
    the DSE energy and its force contribution via the displaced coordinate
    are omitted; molecular forces reduce to pure bilinear coupling.
    """
    force.setIncludeDipoleSelfEnergy(include)


def assign_force_groups(system, include_dipole_self_energy: bool = True) -> Dict[str, int]:
    """Assign force groups for energy decomposition.

    Must be called before Context creation.  Returns a map
    ``{component_name: group_id}`` for use by EnergyTracker.
    """
    import openmm

    group_map: Dict[str, int] = {}
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        if isinstance(force, openmm.HarmonicBondForce):
            force.setForceGroup(1)
            group_map["harmonic_bond"] = 1
        elif isinstance(force, openmm.NonbondedForce):
            force.setForceGroup(0)
            group_map["nonbonded"] = 0
        elif isinstance(force, openmm.CavityForce):
            force.setIncludeDipoleSelfEnergy(include_dipole_self_energy)
            force.setForceGroup(2)
            group_map["cavity"] = 2
        elif isinstance(force, openmm.MultiModeCavityForce):
            force.setIncludeDipoleSelfEnergy(include_dipole_self_energy)
            force.setForceGroup(2)
            group_map["cavity"] = 2
        elif isinstance(force, openmm.BussiThermostat):
            force.setForceGroup(3)
            group_map["bussi"] = 3
        elif isinstance(force, openmm.HarmonicAngleForce):
            force.setForceGroup(4)
            group_map["harmonic_angle"] = 4
        elif isinstance(force, openmm.PeriodicTorsionForce):
            force.setForceGroup(5)
            group_map["torsion"] = 5
        elif isinstance(force, openmm.CMMotionRemover):
            force.setForceGroup(31)
        elif isinstance(force, openmm.CavityParticleDisplacer):
            force.setForceGroup(6)
    return group_map


class CavityMDSimulation:
    """C2F simulation orchestrator for OpenMM cavity-MD.

    Implements a macro-step control loop::

        for each macro-step:
            integrator.step(steps_per_macro)
            evaluate coupling variant -> push to CavityForce
            read energies / temperatures
            run feedback controller
            apply new bath T via DualThermostat
            optional cavity thermostat step
            log state

    Parameters
    ----------
    context : openmm.Context
    cavity_force : openmm.CavityForce or None
    system : openmm.System
    integrator : openmm.Integrator
    variant : CouplingVariant or None
    feedback : EmpiricalTemperatureFeedback or GradientDescentFeedback or None
    thermostat : DualThermostat or None
    energy_tracker : EnergyTracker
    temperature_tracker : TemperatureTracker or None
    runtime_ps : float
    macro_step_ps : float
        Control-loop granularity (default 0.01 ps = 10 fs).
    output_interval_ps : float
        Logging interval (default 1.0 ps).
    output_prefix : str
        Prefix for output CSV files.
    """

    def __init__(
        self,
        context,
        cavity_force,
        system,
        integrator,
        *,
        variant: Optional[CouplingVariant] = None,
        feedback=None,
        thermostat: Optional[DualThermostat] = None,
        energy_tracker: Optional[EnergyTracker] = None,
        temperature_tracker: Optional[TemperatureTracker] = None,
        runtime_ps: float = 100.0,
        macro_step_ps: float = 0.01,
        output_interval_ps: float = 1.0,
        output_prefix: str = "c2f",
    ) -> None:
        self._context = context
        self._cavity_force = cavity_force
        self._system = system
        self._integrator = integrator
        self._variant = variant
        self._feedback = feedback
        self._thermostat = thermostat
        self._etrk = energy_tracker
        self._ttrk = temperature_tracker
        self.runtime_ps = float(runtime_ps)
        self.macro_step_ps = float(macro_step_ps)
        self.output_interval_ps = float(output_interval_ps)
        self._output_prefix = output_prefix

        self._csv_file: Optional[TextIO] = None
        self._next_output_time = 0.0
        self._current_lambda = 0.0

    def run(self) -> None:
        """Execute the complete C2F simulation."""
        from openmm import unit

        dt_ps = self._integrator.getStepSize().value_in_unit(unit.picosecond)
        steps_per_macro = max(1, round(self.macro_step_ps / dt_ps))
        actual_macro_ps = steps_per_macro * dt_ps
        total_macros = int(self.runtime_ps / actual_macro_ps)

        csv_path = f"{self._output_prefix}_energies.csv"
        self._csv_file = open(csv_path, "w")
        self._write_csv_header()

        t0 = wall_time.time()
        try:
            for i in range(total_macros):
                self._integrator.step(steps_per_macro)
                time_ps = self._context.getState().getTime().value_in_unit(unit.picosecond)
                self._macro_step(time_ps, actual_macro_ps)
        finally:
            if self._csv_file is not None:
                self._csv_file.close()

        elapsed = wall_time.time() - t0
        print(f"Simulation complete: {self.runtime_ps:.1f} ps in {elapsed:.1f} s wall time")

    def _macro_step(self, time_ps: float, dt_macro_ps: float) -> None:
        new_lambda = self._current_lambda
        if self._variant is not None:
            new_lambda = self._variant.evaluate(time_ps)
            if new_lambda != self._current_lambda:
                self._cavity_force.setLambdaCoupling(new_lambda)
                self._cavity_force.updateParametersInContext(self._context)
                self._current_lambda = new_lambda

        new_bath_T = None
        if self._feedback is not None and self._thermostat is not None:
            current_T = self._thermostat.get_molecular_temperature()
            if isinstance(
                self._feedback,
                (
                    GradientDescentFeedback,
                    DiffEqController,
                    SimpleSetpointController,
                    PIDControl,
                ),
            ):
                new_bath_T = self._feedback.step(time_ps, current_T)
            elif isinstance(self._feedback, EmpiricalTemperatureFeedback):
                new_bath_T = self._feedback.step(time_ps)
            if new_bath_T is not None:
                self._thermostat.set_molecular_temperature(new_bath_T)

        if self._thermostat is not None and self._thermostat.cavity_friction > 0:
            self._thermostat.apply_cavity_thermostat_step(dt_macro_ps)

        if time_ps >= self._next_output_time:
            self._log_state(time_ps, new_lambda, new_bath_T)
            self._next_output_time += self.output_interval_ps

    def _write_csv_header(self) -> None:
        cols = [
            "time_ps", "lambda", "T_bath_K", "T_kinetic_K",
            "T_v_fictive_K", "T_s_fictive_K",
            "E_bond_kjmol", "E_nonbonded_kjmol",
            "E_cav_harmonic_kjmol", "E_cav_coupling_kjmol", "E_cav_dse_kjmol",
            "E_total_PE_kjmol", "E_total_KE_kjmol",
        ]
        self._csv_file.write(",".join(cols) + "\n")

    def _log_state(self, time_ps: float, lambda_val: float, bath_T: float = None) -> None:
        energies = self._etrk.get_energies() if self._etrk else {}
        temps = self._ttrk.get_all() if self._ttrk else {}

        T_bath = bath_T if bath_T is not None else (
            self._thermostat.get_molecular_temperature() if self._thermostat else 0.0
        )
        T_kin = temps.get("kinetic", 0.0)
        T_v = temps.get("harmonic_equipartition", 0.0)
        T_s = temps.get("structural_fictive")
        T_s_str = f"{T_s:.4f}" if T_s is not None else ""

        row = [
            f"{time_ps:.6f}",
            f"{lambda_val:.8e}",
            f"{T_bath:.4f}",
            f"{T_kin:.4f}",
            f"{T_v:.4f}",
            T_s_str,
            f"{energies.get('harmonic_bond', 0.0):.6f}",
            f"{energies.get('nonbonded', 0.0):.6f}",
            f"{energies.get('cavity_harmonic', 0.0):.6f}",
            f"{energies.get('cavity_coupling', 0.0):.6f}",
            f"{energies.get('cavity_dipole_self', 0.0):.6f}",
            f"{energies.get('total_potential', 0.0):.6f}",
            f"{energies.get('total_kinetic', 0.0):.6f}",
        ]
        self._csv_file.write(",".join(row) + "\n")
        self._csv_file.flush()
