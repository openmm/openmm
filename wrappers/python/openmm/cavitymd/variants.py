from abc import ABC, abstractmethod
from typing import Optional, List, Dict, Callable
import math


class CouplingVariant(ABC):
    """Base class for time-varying coupling modulation.

    Subclasses implement evaluate(time_ps) -> float, returning the
    dimensionless lambda coupling value at a given simulation time.
    The simulation loop pushes this into CavityForce via
    setLambdaCoupling() + updateParametersInContext().
    """

    @abstractmethod
    def evaluate(self, time_ps: float) -> float:
        ...

    @property
    @abstractmethod
    def target_value(self) -> float:
        ...


class ConstantVariant(CouplingVariant):
    """Time-independent coupling: g(t) = value for all t."""

    def __init__(self, value: float) -> None:
        self._value = float(value)

    def evaluate(self, time_ps: float) -> float:
        return self._value

    @property
    def target_value(self) -> float:
        return self._value


class StepVariant(CouplingVariant):
    r"""Step function with optional exponential decay.

    Without decay::

        g(t) = 0                    if t < t_switch
        g(t) = g_target             if t >= t_switch

    With decay::

        g(t) = g_target * exp(-(t - t_switch) / tau_decay)   if t >= t_switch

    Parameters
    ----------
    target_value : float
        Lambda coupling value after switching.
    switch_time_ps : float
        Time (ps) at which coupling activates.
    decay_time_constant_ps : float or None
        Exponential decay time constant (ps). None = no decay.
    turn_off_time_ps : float or None
        Time (ps) at which coupling returns to zero.
    """

    def __init__(
        self,
        target_value: float,
        switch_time_ps: float,
        decay_time_constant_ps: Optional[float] = None,
        turn_off_time_ps: Optional[float] = None,
    ) -> None:
        self._target = float(target_value)
        self.switch_time_ps = float(switch_time_ps)
        self.decay_time_constant_ps = (
            float(decay_time_constant_ps) if decay_time_constant_ps is not None else None
        )
        self.turn_off_time_ps = (
            float(turn_off_time_ps) if turn_off_time_ps is not None else None
        )

    def evaluate(self, time_ps: float) -> float:
        if self.turn_off_time_ps is not None and time_ps >= self.turn_off_time_ps:
            return 0.0
        if time_ps < self.switch_time_ps:
            return 0.0
        if self.decay_time_constant_ps is not None:
            dt = time_ps - self.switch_time_ps
            return self._target * math.exp(-dt / self.decay_time_constant_ps)
        return self._target

    @property
    def target_value(self) -> float:
        return self._target


class SquareWaveVariant(CouplingVariant):
    r"""Periodic on/off pulsing with configurable duty cycle.

    During the active window [start_time, stop_time)::

        phase = 2*pi*(t - start_time) / period + phase_offset
        phase_mod = phase mod 2*pi
        g(t) = amplitude   if phase_mod < duty_cycle * 2*pi
        g(t) = 0           otherwise

    Parameters
    ----------
    amplitude : float
        Lambda value when "high".
    period_ps : float
        Full period in ps.
    duty_cycle : float
        Fraction of period that is "high" (0 to 1). Default 0.5.
    start_time_ps : float
        When pulsing begins (ps). Default 0.
    stop_time_ps : float or None
        When pulsing ends (None = indefinite).
    """

    def __init__(
        self,
        amplitude: float,
        period_ps: float,
        duty_cycle: float = 0.5,
        start_time_ps: float = 0.0,
        stop_time_ps: Optional[float] = None,
    ) -> None:
        if not 0.0 <= duty_cycle <= 1.0:
            raise ValueError(f"duty_cycle must be in [0,1], got {duty_cycle}")
        self._amplitude = float(amplitude)
        self.period_ps = float(period_ps)
        self.duty_cycle = float(duty_cycle)
        self.start_time_ps = float(start_time_ps)
        self.stop_time_ps = (
            float(stop_time_ps) if stop_time_ps is not None else None
        )
        self._angular_freq = 2.0 * math.pi / self.period_ps
        self._high_phase = self.duty_cycle * 2.0 * math.pi
        self.is_high = False

    def evaluate(self, time_ps: float) -> float:
        if time_ps < self.start_time_ps:
            self.is_high = False
            return 0.0
        if self.stop_time_ps is not None and time_ps >= self.stop_time_ps:
            self.is_high = False
            return 0.0

        dt = time_ps - self.start_time_ps
        phase_mod = (self._angular_freq * dt) % (2.0 * math.pi)

        if phase_mod < self._high_phase:
            self.is_high = True
            return self._amplitude
        else:
            self.is_high = False
            return 0.0

    @property
    def target_value(self) -> float:
        return self._amplitude

    @property
    def current_period(self) -> int:
        return 0  # updated externally if needed


class AdaptiveSquareWaveVariant(CouplingVariant):
    r"""Square wave with amplitude adapted by temperature feedback.

    At each high-to-low transition the amplitude is updated::

        g_next = g_target * sqrt(T_target / T_bath)

    where T_bath is read via the supplied callable.

    Parameters
    ----------
    target_coupling : float
        Reference coupling g_target.
    target_temperature_K : float
        Desired molecular temperature (K).
    period_ps : float
        Square wave period (ps).
    get_bath_temperature : callable
        () -> float returning current bath temperature in K.
    duty_cycle : float
        Fraction of period that is "high". Default 0.5.
    start_time_ps : float
        When pulsing begins.
    stop_time_ps : float or None
        When pulsing ends.
    min_amplitude, max_amplitude : float
        Safety clamps on the adaptive amplitude.
    """

    def __init__(
        self,
        target_coupling: float,
        target_temperature_K: float,
        period_ps: float,
        get_bath_temperature: Callable[[], float],
        duty_cycle: float = 0.5,
        start_time_ps: float = 0.0,
        stop_time_ps: Optional[float] = None,
        min_amplitude: float = 1e-8,
        max_amplitude: float = 1e-1,
    ) -> None:
        if not 0.0 <= duty_cycle <= 1.0:
            raise ValueError(f"duty_cycle must be in [0,1], got {duty_cycle}")
        if target_temperature_K <= 0.0:
            raise ValueError("target_temperature_K must be positive")

        self._target_coupling = float(target_coupling)
        self.target_temperature_K = float(target_temperature_K)
        self.period_ps = float(period_ps)
        self._get_bath_temperature = get_bath_temperature
        self.duty_cycle = float(duty_cycle)
        self.start_time_ps = float(start_time_ps)
        self.stop_time_ps = (
            float(stop_time_ps) if stop_time_ps is not None else None
        )
        self.min_amplitude = float(min_amplitude)
        self.max_amplitude = float(max_amplitude)

        self._angular_freq = 2.0 * math.pi / self.period_ps
        self._high_phase = self.duty_cycle * 2.0 * math.pi

        self.current_amplitude = float(target_coupling)
        self.is_high = False
        self._was_high = False
        self._amplitude_history: List[Dict] = []

    def evaluate(self, time_ps: float) -> float:
        if time_ps < self.start_time_ps:
            self.is_high = False
            self._was_high = False
            return 0.0
        if self.stop_time_ps is not None and time_ps >= self.stop_time_ps:
            self.is_high = False
            self._was_high = False
            return 0.0

        dt = time_ps - self.start_time_ps
        phase_mod = (self._angular_freq * dt) % (2.0 * math.pi)

        self._was_high, self.is_high = self.is_high, phase_mod < self._high_phase

        if self._was_high and not self.is_high:
            self._update_amplitude(time_ps)

        return self.current_amplitude if self.is_high else 0.0

    def _update_amplitude(self, time_ps: float) -> None:
        try:
            T_bath = self._get_bath_temperature()
            if T_bath is None or T_bath <= 0.0:
                return
            ratio = self.target_temperature_K / T_bath
            new_amp = self._target_coupling * math.sqrt(ratio)
            new_amp = max(self.min_amplitude, min(self.max_amplitude, new_amp))
            self._amplitude_history.append(
                {
                    "time_ps": time_ps,
                    "T_bath": T_bath,
                    "T_target": self.target_temperature_K,
                    "old_amplitude": self.current_amplitude,
                    "new_amplitude": new_amp,
                }
            )
            self.current_amplitude = new_amp
        except Exception:
            pass

    @property
    def target_value(self) -> float:
        return self._target_coupling

    @property
    def amplitude_history(self) -> List[Dict]:
        return list(self._amplitude_history)


class ExponentialDecayVariant(CouplingVariant):
    """Exponential decay from an initial amplitude after turn-on."""

    def __init__(
        self,
        amplitude: float,
        decay_time_constant_ps: float,
        turn_on_time_ps: float = 0.0,
        turn_off_time_ps: Optional[float] = None,
    ) -> None:
        self._amplitude = float(amplitude)
        self.decay_time_constant_ps = float(decay_time_constant_ps)
        self.turn_on_time_ps = float(turn_on_time_ps)
        self.turn_off_time_ps = (
            float(turn_off_time_ps) if turn_off_time_ps is not None else None
        )

    def evaluate(self, time_ps: float) -> float:
        if time_ps < self.turn_on_time_ps:
            return 0.0
        if self.turn_off_time_ps is not None and time_ps >= self.turn_off_time_ps:
            return 0.0
        dt = time_ps - self.turn_on_time_ps
        return self._amplitude * math.exp(-dt / self.decay_time_constant_ps)

    @property
    def target_value(self) -> float:
        return self._amplitude


class DecayingSquareWaveVariant(CouplingVariant):
    """Square wave with amplitude A_n = A0 * (1-r)^n after each period."""

    def __init__(
        self,
        initial_amplitude: float,
        period_ps: float,
        decay_rate_per_period: float = 0.0,
        duty_cycle: float = 0.5,
        start_time_ps: float = 0.0,
        stop_time_ps: Optional[float] = None,
        minimum_amplitude: float = 1e-8,
    ) -> None:
        if not 0.0 <= duty_cycle <= 1.0:
            raise ValueError(f"duty_cycle must be in [0,1], got {duty_cycle}")
        if not 0.0 <= decay_rate_per_period <= 1.0:
            raise ValueError(
                f"decay_rate_per_period must be in [0,1], got {decay_rate_per_period}"
            )
        self._initial_amplitude = float(initial_amplitude)
        self.current_amplitude = float(initial_amplitude)
        self.period_ps = float(period_ps)
        self.decay_rate_per_period = float(decay_rate_per_period)
        self.duty_cycle = float(duty_cycle)
        self.start_time_ps = float(start_time_ps)
        self.stop_time_ps = (
            float(stop_time_ps) if stop_time_ps is not None else None
        )
        self.minimum_amplitude = float(minimum_amplitude)
        self._completed_periods = 0

    def evaluate(self, time_ps: float) -> float:
        if time_ps < self.start_time_ps:
            return 0.0
        if self.stop_time_ps is not None and time_ps >= self.stop_time_ps:
            return 0.0

        dt = time_ps - self.start_time_ps
        period_index = int(dt / self.period_ps)
        if period_index > self._completed_periods:
            for _ in range(period_index - self._completed_periods):
                self.current_amplitude *= 1.0 - self.decay_rate_per_period
            self._completed_periods = period_index

        if self.current_amplitude < self.minimum_amplitude:
            return 0.0

        phase = (dt / self.period_ps) % 1.0
        return self.current_amplitude if phase < self.duty_cycle else 0.0

    @property
    def target_value(self) -> float:
        return self._initial_amplitude


class SinusoidVariant(CouplingVariant):
    """Sinusoidal coupling: A * (1 + sin(2*pi*dt/period + phase)) / 2."""

    def __init__(
        self,
        amplitude: float,
        period_ps: float,
        phase_offset: float = 0.0,
        start_time_ps: float = 0.0,
        stop_time_ps: Optional[float] = None,
    ) -> None:
        self._amplitude = float(amplitude)
        self.period_ps = float(period_ps)
        self.phase_offset = float(phase_offset)
        self.start_time_ps = float(start_time_ps)
        self.stop_time_ps = (
            float(stop_time_ps) if stop_time_ps is not None else None
        )

    def evaluate(self, time_ps: float) -> float:
        if time_ps < self.start_time_ps:
            return 0.0
        if self.stop_time_ps is not None and time_ps >= self.stop_time_ps:
            return 0.0
        dt = time_ps - self.start_time_ps
        phase = 2.0 * math.pi * dt / self.period_ps + self.phase_offset
        return self._amplitude * (1.0 + math.sin(phase)) / 2.0

    @property
    def target_value(self) -> float:
        return self._amplitude


class ExponentialWaveVariant(CouplingVariant):
    """Periodic exponential pulses: A * exp(-t_in_period / tau) each period."""

    def __init__(
        self,
        amplitude: float,
        period_ps: float,
        decay_tau_ps: float,
        start_time_ps: float = 0.0,
        stop_time_ps: Optional[float] = None,
    ) -> None:
        self._amplitude = float(amplitude)
        self.period_ps = float(period_ps)
        self.decay_tau_ps = float(decay_tau_ps)
        self.start_time_ps = float(start_time_ps)
        self.stop_time_ps = (
            float(stop_time_ps) if stop_time_ps is not None else None
        )

    def evaluate(self, time_ps: float) -> float:
        if time_ps < self.start_time_ps:
            return 0.0
        if self.stop_time_ps is not None and time_ps >= self.stop_time_ps:
            return 0.0
        dt = time_ps - self.start_time_ps
        t_in_period = dt % self.period_ps
        return self._amplitude * math.exp(-t_in_period / self.decay_tau_ps)

    @property
    def target_value(self) -> float:
        return self._amplitude
