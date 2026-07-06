"""SI-accurate adaptive timestepping helpers for C2F cavity-MD.

Implements Eq. 3.15–3.16 from the paper SI:
  ε* = 5.0 nm, f₀ = 10⁻³, τ* = 50 ps, dt cap ≈ 1.5 fs.
Error tolerance resets on every square-wave λ edge (rising and falling).
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional

# SI Section 3.2 constants (NOT dimer reference values 1.0 / 1e-5 nm)
EPS_STAR_NM = 5.0
F0 = 1e-3
TAU_RAMP_PS = 50.0
DT_MAX_PS = 0.0015  # ~1.5 fs stable maximum


def epsilon_tolerance(time_ps: float, ramp_start_ps: Optional[float]) -> float:
    """SI Eq. 3.16: strict at t₀, relax to EPS_STAR over TAU."""
    if ramp_start_ps is None or time_ps < ramp_start_ps:
        return EPS_STAR_NM
    t_since = time_ps - ramp_start_ps
    return EPS_STAR_NM * (1.0 - (1.0 - F0) * math.exp(-t_since / TAU_RAMP_PS))


def square_wave_on(
    time_ps: float,
    start_ps: float,
    period_ps: float,
    duty: float,
) -> bool:
    """Return True when analytic square-wave λ(t) is ON."""
    if time_ps < start_ps or period_ps <= 0.0:
        return False
    phase = ((time_ps - start_ps) / period_ps) % 1.0
    return phase < duty


def lambda_transition(prev_on: bool, curr_on: bool) -> bool:
    """True on every rising or falling square-wave edge."""
    return prev_on != curr_on


def create_adaptive_integrator(dt_max_ps: float = DT_MAX_PS):
    """Create a VariableVerletIntegrator with SI default relaxed tolerance."""
    import openmm
    from openmm import unit

    integrator = openmm.VariableVerletIntegrator(EPS_STAR_NM)
    integrator.setMaximumStepSize(dt_max_ps * unit.picosecond)
    return integrator


def advance_to_time(
    context,
    integrator,
    thermostat,
    displacer,
    lambda_coupling: float,
    target_time_ps: float,
    coupling_start_ps: float,
    period_ps: float,
    duty: float,
    state: Dict[str, Any],
    log_dt: bool = False,
) -> List[float]:
    """Integrate from current time to *target_time_ps* with adaptive tolerance.

    Updates *state* keys ``ramp_t0`` and ``prev_lambda_on`` in place.
    Applies cavity Langevin after each internal MD step.

    Parameters
    ----------
    state : dict
        Must contain ``ramp_t0`` (Optional[float]) and ``prev_lambda_on`` (bool).

    Returns
    -------
    list of float
        Step sizes (ps) when *log_dt* is True, else empty list.
    """
    from openmm import unit

    logged_dts: List[float] = []
    target_time_ps = float(target_time_ps)

    while True:
        time_ps = context.getState().getTime().value_in_unit(unit.picosecond)
        if time_ps >= target_time_ps - 1e-15:
            break

        curr_on = square_wave_on(time_ps, coupling_start_ps, period_ps, duty)
        if lambda_transition(state["prev_lambda_on"], curr_on):
            state["ramp_t0"] = time_ps
            if curr_on and displacer is not None and lambda_coupling > 0.0:
                displacer.displaceToEquilibrium(context, lambda_coupling)

        integrator.setErrorTolerance(
            epsilon_tolerance(time_ps, state["ramp_t0"])
        )
        state["prev_lambda_on"] = curr_on

        integrator.step(1)
        dt_ps = integrator.getStepSize().value_in_unit(unit.picosecond)
        thermostat.apply_cavity_thermostat_step(dt_ps)
        if log_dt:
            logged_dts.append(dt_ps)

    return logged_dts
