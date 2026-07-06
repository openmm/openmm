"""Feedback controllers for C2F cavity-MD (ported from cav-hoomd)."""

from typing import Optional
import numpy as np


def _measure_temperature(temperature_tracker, method: str) -> Optional[float]:
    """Read a temperature signal from TemperatureTracker."""
    if method == "kinetic":
        return temperature_tracker.kinetic_temperature
    if method in ("harmonic_equipartition", "harmonic", "harmonic_fictive"):
        return temperature_tracker.harmonic_equipartition_temperature
    if method in ("lj_coulombic", "structural", "structural_fictive"):
        return temperature_tracker.structural_fictive_temperature
    return None


def _apply_clamps_and_rate_limit(
    raw: float,
    current_bath_T: float,
    T_min: float,
    T_max: Optional[float],
    rate_limit: Optional[float],
    dt_ps: float,
) -> float:
    clamped = max(T_min, raw)
    if T_max is not None:
        clamped = min(T_max, clamped)
    if rate_limit is not None:
        max_change = rate_limit * dt_ps
        delta = clamped - current_bath_T
        if abs(delta) > max_change:
            clamped = current_bath_T + max_change * (1.0 if delta > 0 else -1.0)
    return clamped


class _BaseController:
    """Shared activation window and update-interval logic."""

    def __init__(
        self,
        temperature_tracker,
        update_interval_ps: float = 0.1,
        T_min: float = 0.0,
        T_max: Optional[float] = None,
        rate_limit_K_per_ps: Optional[float] = None,
        turn_on_time_ps: float = 0.0,
        turn_off_time_ps: Optional[float] = None,
    ) -> None:
        self._ttrk = temperature_tracker
        self.update_interval_ps = float(update_interval_ps)
        self.T_min = float(T_min)
        self.T_max = float(T_max) if T_max is not None else None
        self.rate_limit = (
            float(rate_limit_K_per_ps) if rate_limit_K_per_ps is not None else None
        )
        self.turn_on_time_ps = float(turn_on_time_ps)
        self.turn_off_time_ps = (
            float(turn_off_time_ps) if turn_off_time_ps is not None else None
        )
        self._last_update_time = -1e10
        self._is_active = False

    def _check_active(self, time_ps: float) -> bool:
        active = (
            time_ps >= self.turn_on_time_ps
            and (self.turn_off_time_ps is None or time_ps < self.turn_off_time_ps)
        )
        self._is_active = active
        return active

    def _should_update(self, time_ps: float) -> Optional[float]:
        if not self._check_active(time_ps):
            return None
        if time_ps - self._last_update_time < self.update_interval_ps:
            return None
        dt_ps = (
            time_ps - self._last_update_time
            if self._last_update_time > 0
            else self.update_interval_ps
        )
        self._last_update_time = time_ps
        return dt_ps

    @property
    def is_active(self) -> bool:
        return self._is_active


class DiffEqController(_BaseController):
    """First-order tracking: dT_bath/dt = -(T_bath - T_signal) / tau.

    With ``temperature_method='structural'`` (default), T_signal is the
    structural fictive temperature T_s — the C2F bath-tracking step from
    the paper.
    """

    def __init__(
        self,
        temperature_tracker,
        time_constant_ps: float = 5.0,
        temperature_method: str = "structural",
        relaxation_model=None,
        time_constant_auto: bool = False,
        update_interval_ps: float = 5.0,
        T_min: float = 0.0,
        T_max: Optional[float] = None,
        rate_limit_K_per_ps: Optional[float] = None,
        turn_on_time_ps: float = 0.0,
        turn_off_time_ps: Optional[float] = None,
    ) -> None:
        super().__init__(
            temperature_tracker,
            update_interval_ps=update_interval_ps,
            T_min=T_min,
            T_max=T_max,
            rate_limit_K_per_ps=rate_limit_K_per_ps,
            turn_on_time_ps=turn_on_time_ps,
            turn_off_time_ps=turn_off_time_ps,
        )
        self.time_constant_ps = float(time_constant_ps)
        self.temperature_method = temperature_method
        self.relaxation_model = relaxation_model
        self.time_constant_auto = bool(time_constant_auto)

    def _effective_tau(self, bath_T: float) -> float:
        if self.time_constant_auto and self.relaxation_model is not None:
            return self.relaxation_model.get_relaxation_time(bath_T)
        return self.time_constant_ps

    def step(
        self,
        time_ps: float,
        current_bath_T: float,
        signal_override: Optional[float] = None,
        force: bool = False,
        dt_ps_override: Optional[float] = None,
    ) -> Optional[float]:
        if force:
            if not self._check_active(time_ps):
                return None
            if dt_ps_override is not None:
                dt_ps = float(dt_ps_override)
            else:
                dt_ps = self.update_interval_ps
            self._last_update_time = time_ps
        else:
            dt_ps = self._should_update(time_ps)
            if dt_ps is None:
                return None

        if signal_override is not None:
            T_signal = signal_override
        else:
            T_signal = _measure_temperature(self._ttrk, self.temperature_method)
        if T_signal is None:
            return None

        tau = self._effective_tau(current_bath_T)
        if tau <= 0:
            return None

        # Exponential integrator avoids Euler overshoot when dt >> tau (e.g. λ-off window).
        alpha = np.exp(-dt_ps / tau)
        raw = T_signal + (current_bath_T - T_signal) * alpha
        return _apply_clamps_and_rate_limit(
            raw, current_bath_T, self.T_min, self.T_max, self.rate_limit, dt_ps
        )


class SimpleSetpointController(_BaseController):
    """Capture a fixed setpoint at turn-on, then drive kinetic T toward it.

    dT_bath/dt = -(T_kinetic - T_setpoint) / tau
    """

    def __init__(
        self,
        temperature_tracker,
        signal_method: str = "structural",
        time_constant_ps: float = 5.0,
        update_interval_ps: float = 5.0,
        T_min: float = 0.0,
        T_max: Optional[float] = None,
        rate_limit_K_per_ps: Optional[float] = None,
        turn_on_time_ps: float = 0.0,
        turn_off_time_ps: Optional[float] = None,
    ) -> None:
        super().__init__(
            temperature_tracker,
            update_interval_ps=update_interval_ps,
            T_min=T_min,
            T_max=T_max,
            rate_limit_K_per_ps=rate_limit_K_per_ps,
            turn_on_time_ps=turn_on_time_ps,
            turn_off_time_ps=turn_off_time_ps,
        )
        self.signal_method = signal_method
        self.time_constant_ps = float(time_constant_ps)
        self.setpoint_temperature: Optional[float] = None
        self._setpoint_captured = False

    def step(self, time_ps: float, current_bath_T: float) -> Optional[float]:
        if not self._check_active(time_ps):
            return None

        if not self._setpoint_captured:
            sp = _measure_temperature(self._ttrk, self.signal_method)
            if sp is None:
                return None
            self.setpoint_temperature = sp
            self._setpoint_captured = True

        dt_ps = self._should_update(time_ps)
        if dt_ps is None:
            return None

        T_kin = _measure_temperature(self._ttrk, "kinetic")
        if T_kin is None or self.setpoint_temperature is None:
            return None

        error = T_kin - self.setpoint_temperature
        derivative = -error / self.time_constant_ps
        raw = current_bath_T + derivative * dt_ps
        return _apply_clamps_and_rate_limit(
            raw, current_bath_T, self.T_min, self.T_max, self.rate_limit, dt_ps
        )


class PIDControl(_BaseController):
    """Classical PID bath-temperature controller.

    u(t) = feedforward + Kp*(e + (1/Ti)*integral(e) + Td*de/dt)
    where e = T_setpoint - T_signal.
    """

    def __init__(
        self,
        temperature_tracker,
        signal_choice: str = "structural",
        target_temperature: float = 100.0,
        Kp: float = 0.5,
        Ti: float = 20.0,
        Td: float = 5.0,
        self_loop: bool = False,
        update_interval_ps: float = 0.1,
        T_min: float = 0.1,
        T_max: Optional[float] = None,
        rate_limit_K_per_ps: Optional[float] = None,
        turn_on_time_ps: float = 0.0,
        turn_off_time_ps: Optional[float] = None,
        derivative_filter_N: float = 10.0,
        enable_anti_windup: bool = True,
    ) -> None:
        super().__init__(
            temperature_tracker,
            update_interval_ps=update_interval_ps,
            T_min=T_min,
            T_max=T_max,
            rate_limit_K_per_ps=rate_limit_K_per_ps,
            turn_on_time_ps=turn_on_time_ps,
            turn_off_time_ps=turn_off_time_ps,
        )
        self.signal_choice = signal_choice
        self.target_temperature = float(target_temperature)
        self.self_loop = bool(self_loop)
        self.Kp = float(Kp)
        self.Ti = float(Ti)
        self.Td = float(Td)
        self.derivative_filter_N = float(derivative_filter_N)
        self.enable_anti_windup = bool(enable_anti_windup)

        self._integral = 0.0
        self._derivative_filtered = 0.0
        self._error_prev = 0.0
        self._bath_prev: Optional[float] = None

    def step(self, time_ps: float, current_bath_T: float) -> Optional[float]:
        dt_ps = self._should_update(time_ps)
        if dt_ps is None:
            return None

        T_signal = _measure_temperature(self._ttrk, self.signal_choice)
        if T_signal is None:
            return None

        if self.self_loop:
            if self._bath_prev is None:
                self._bath_prev = current_bath_T
            T_setpoint = self._bath_prev
        else:
            T_setpoint = self.target_temperature

        error = T_setpoint - T_signal
        P = self.Kp * error

        if self.Ti > 0:
            self._integral += dt_ps * (error + self._error_prev) / 2.0
            I = self.Kp * (1.0 / self.Ti) * self._integral
        else:
            I = 0.0

        if self.Td > 0:
            tau_f = self.Td / self.derivative_filter_N
            alpha = tau_f / (tau_f + dt_ps)
            deriv_raw = self.Kp * self.Td * (error - self._error_prev) / dt_ps
            self._derivative_filtered = (
                alpha * self._derivative_filtered + (1.0 - alpha) * deriv_raw
            )
            D = self._derivative_filtered
        else:
            D = 0.0

        if self.self_loop:
            u_unsat = self._bath_prev + (P + I + D)
        else:
            u_unsat = T_setpoint + (P + I + D)

        u_sat = _apply_clamps_and_rate_limit(
            u_unsat, current_bath_T, self.T_min, self.T_max, self.rate_limit, dt_ps
        )

        if self.enable_anti_windup and self.Ti > 0:
            self._integral += (dt_ps / self.Ti) * (u_sat - u_unsat)

        self._error_prev = error
        if self.self_loop:
            self._bath_prev = u_sat
        return u_sat
