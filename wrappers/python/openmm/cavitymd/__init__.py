"""C2F (Cavity Configurational Feedback) protocol for OpenMM cavity-MD.

Provides coupling modulation, fictive temperature computation, feedback
controllers, and a simulation orchestrator for cavity molecular dynamics
with the C2F cooling protocol.
"""

from .constants import Units
from .variants import (
    CouplingVariant,
    ConstantVariant,
    StepVariant,
    SquareWaveVariant,
    AdaptiveSquareWaveVariant,
    ExponentialDecayVariant,
    DecayingSquareWaveVariant,
    SinusoidVariant,
    ExponentialWaveVariant,
)
from .empirical import EmpiricalTemperatureData
from .calibration import (
    run_nvt_energy_calibration,
    run_legacy_equilibrium_calibration,
    validate_calibration_file,
    crosscheck_calibration_against_reference,
    check_calibration_sanity,
    BUSSI_TAU_PS,
)
from .trackers import ElapsedTimeTracker, EnergyTracker, TemperatureTracker, compute_harmonic_bond_energy_kjmol
from .feedback import EmpiricalTemperatureFeedback, GradientDescentFeedback
from .controllers import DiffEqController, SimpleSetpointController, PIDControl
from .analysis import RelaxationTimeModel, ToolNarayanaswamy
from .thermostats import (
    DEFAULT_LANGEVIN_FRICTION_PS,
    DualThermostat,
    create_langevin_integrator,
)
from .adaptive import (
    EPS_STAR_NM,
    F0,
    TAU_RAMP_PS,
    DT_MAX_PS,
    epsilon_tolerance,
    square_wave_on,
    lambda_transition,
    create_adaptive_integrator,
    advance_to_time,
)
from .simulation import (
    CavityMDSimulation,
    assign_force_groups,
    setup_gpu_square_wave,
    setup_gpu_step,
    setup_gpu_decaying_step,
    setup_gpu_adaptive_square_wave,
    setup_gpu_decaying_square_wave,
    setup_gpu_sinusoid,
    setup_gpu_exponential_wave,
    setup_multimode_adaptive_square_wave,
)

__all__ = [
    "Units",
    "CouplingVariant",
    "ConstantVariant",
    "StepVariant",
    "SquareWaveVariant",
    "AdaptiveSquareWaveVariant",
    "ExponentialDecayVariant",
    "DecayingSquareWaveVariant",
    "SinusoidVariant",
    "ExponentialWaveVariant",
    "EmpiricalTemperatureData",
    "run_nvt_energy_calibration",
    "run_legacy_equilibrium_calibration",
    "validate_calibration_file",
    "crosscheck_calibration_against_reference",
    "check_calibration_sanity",
    "BUSSI_TAU_PS",
    "ElapsedTimeTracker",
    "EnergyTracker",
    "TemperatureTracker",
    "compute_harmonic_bond_energy_kjmol",
    "EmpiricalTemperatureFeedback",
    "GradientDescentFeedback",
    "DiffEqController",
    "SimpleSetpointController",
    "PIDControl",
    "RelaxationTimeModel",
    "ToolNarayanaswamy",
    "DEFAULT_LANGEVIN_FRICTION_PS",
    "DualThermostat",
    "create_langevin_integrator",
    "EPS_STAR_NM",
    "F0",
    "TAU_RAMP_PS",
    "DT_MAX_PS",
    "epsilon_tolerance",
    "square_wave_on",
    "lambda_transition",
    "create_adaptive_integrator",
    "advance_to_time",
    "CavityMDSimulation",
    "assign_force_groups",
    "setup_gpu_square_wave",
    "setup_gpu_step",
    "setup_gpu_decaying_step",
    "setup_gpu_adaptive_square_wave",
    "setup_gpu_decaying_square_wave",
    "setup_gpu_sinusoid",
    "setup_gpu_exponential_wave",
    "setup_multimode_adaptive_square_wave",
]
