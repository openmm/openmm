from typing import Optional, List, Tuple
import numpy as np

from .constants import Units


class EmpiricalTemperatureFeedback:
    """Windowed-average feedback: measure energy -> T_fictive -> set bath T.

    Control loop (called each macro-step by the orchestrator):
    1. Read LJ+Coulomb energy via EnergyTracker.
    2. Convert to structural fictive temperature via EmpiricalTemperatureData.
    3. Accumulate in sliding window (default 5 ps).
    4. Every update_interval_ps, return clamped mean as new bath T.

    Parameters
    ----------
    empirical_data : EmpiricalTemperatureData
    energy_tracker : EnergyTracker
    update_interval_ps : float
        How often to emit a new target temperature. Default 5.0.
    averaging_window_ps : float
        Width of sliding measurement window. Default 5.0.
    T_min, T_max : float
        Temperature clamps (K).
    switch_time_ps : float or None
        Time when feedback activates (skip before coupling turns on).
    turn_off_time_ps : float or None
        Time to deactivate feedback.
    energy_component : str
        Which energy to read from tracker. Default 'nonbonded'.
    """

    def __init__(
        self,
        empirical_data,
        energy_tracker,
        update_interval_ps: float = 5.0,
        averaging_window_ps: float = 5.0,
        T_min: float = 0.0,
        T_max: float = 500.0,
        switch_time_ps: Optional[float] = None,
        turn_off_time_ps: Optional[float] = None,
        energy_component: str = "nonbonded",
    ) -> None:
        self._empirical = empirical_data
        self._etrk = energy_tracker
        self.update_interval_ps = float(update_interval_ps)
        self.averaging_window_ps = float(averaging_window_ps)
        self.T_min = float(T_min)
        self.T_max = float(T_max)
        self.switch_time_ps = float(switch_time_ps) if switch_time_ps is not None else None
        self.turn_off_time_ps = float(turn_off_time_ps) if turn_off_time_ps is not None else None
        self._energy_key = energy_component

        self._measurements: List[Tuple[float, float]] = []
        self._last_update_time = -1e10
        self._last_applied_T: Optional[float] = None
        self._active = True

    def step(self, time_ps: float) -> Optional[float]:
        """Execute one feedback step. Returns new bath T if update applied."""
        if self.switch_time_ps is not None and time_ps < self.switch_time_ps:
            return None
        if self.turn_off_time_ps is not None and time_ps >= self.turn_off_time_ps:
            self._active = False
            return None
        if not self._active:
            return None

        energies_h = self._etrk.get_energies_hartree()
        E = energies_h.get(self._energy_key, 0.0)
        T_inst = self._empirical.calculate_temperature(E)

        self._measurements.append((time_ps, T_inst))

        cutoff = time_ps - self.averaging_window_ps
        self._measurements = [(t, v) for t, v in self._measurements if t >= cutoff]

        if time_ps - self._last_update_time < self.update_interval_ps:
            return None

        self._last_update_time = time_ps
        temps = [v for _, v in self._measurements]
        avg_T = float(np.mean(temps)) if temps else T_inst
        target_T = max(self.T_min, min(self.T_max, avg_T))
        self._last_applied_T = target_T
        return target_T

    @property
    def is_active(self) -> bool:
        return self._active

    @property
    def last_applied_temperature(self) -> Optional[float]:
        return self._last_applied_T


class GradientDescentFeedback:
    """Gradient descent controller for bath temperature tracking.

    Objective: J = 0.5 * (T_eff - T_target)^2
    where T_eff = (T_measured + T_bath) / 2.

    Update rule::

        T_bath_new = T_bath_old - alpha * 0.5 * (T_eff - T_target)
        alpha = update_interval_ps / tau_time_constant

    Parameters
    ----------
    temperature_method : str
        'kinetic', 'harmonic_equipartition', 'lj_coulombic', 'harmonic'.
    time_constant_ps : float
        Controls convergence speed.
    target_temperature_K : float
    temperature_tracker : TemperatureTracker
    update_interval_ps : float
        How often to apply GD update. Default 0.1.
    T_min, T_max : float
        Bath temperature clamps (K).
    rate_limit_K_per_ps : float or None
        Maximum dT/dt.
    turn_on_time_ps, turn_off_time_ps : float or None
        Activation window.
    """

    def __init__(
        self,
        temperature_method: str,
        time_constant_ps: float,
        target_temperature_K: float,
        temperature_tracker,
        update_interval_ps: float = 0.1,
        T_min: float = 0.0,
        T_max: Optional[float] = None,
        rate_limit_K_per_ps: Optional[float] = None,
        turn_on_time_ps: float = 0.0,
        turn_off_time_ps: Optional[float] = None,
    ) -> None:
        self.temperature_method = temperature_method
        self.time_constant_ps = float(time_constant_ps)
        self.target_temperature_K = float(target_temperature_K)
        self._ttrk = temperature_tracker
        self.update_interval_ps = float(update_interval_ps)
        self.T_min = float(T_min)
        self.T_max = float(T_max) if T_max is not None else None
        self.rate_limit = float(rate_limit_K_per_ps) if rate_limit_K_per_ps is not None else None
        self.turn_on_time_ps = float(turn_on_time_ps)
        self.turn_off_time_ps = float(turn_off_time_ps) if turn_off_time_ps is not None else None

        self.alpha = self.update_interval_ps / self.time_constant_ps
        self._last_update_time = -1e10
        self._is_active = False

    def step(self, time_ps: float, current_bath_T: float) -> Optional[float]:
        """Execute one GD step. Returns new bath T or None."""
        should_be_active = (
            time_ps >= self.turn_on_time_ps
            and (self.turn_off_time_ps is None or time_ps < self.turn_off_time_ps)
        )
        if not should_be_active:
            self._is_active = False
            return None
        self._is_active = True

        if time_ps - self._last_update_time < self.update_interval_ps:
            return None
        dt_ps = time_ps - self._last_update_time if self._last_update_time > 0 else self.update_interval_ps
        self._last_update_time = time_ps

        T_meas = self._measure_temperature()
        if T_meas is None:
            return None

        T_eff = (T_meas + current_bath_T) / 2.0
        error = T_eff - self.target_temperature_K
        gradient = 0.5 * error

        raw = current_bath_T - self.alpha * gradient

        clamped = max(self.T_min, raw)
        if self.T_max is not None:
            clamped = min(self.T_max, clamped)

        if self.rate_limit is not None:
            max_change = self.rate_limit * dt_ps
            delta = clamped - current_bath_T
            if abs(delta) > max_change:
                clamped = current_bath_T + max_change * (1.0 if delta > 0 else -1.0)

        return clamped

    def _measure_temperature(self) -> Optional[float]:
        if self.temperature_method == "kinetic":
            return self._ttrk.kinetic_temperature
        elif self.temperature_method == "harmonic_equipartition":
            return self._ttrk.harmonic_equipartition_temperature
        elif self.temperature_method in ("lj_coulombic", "structural"):
            return self._ttrk.structural_fictive_temperature
        elif self.temperature_method == "harmonic":
            return self._ttrk.harmonic_equipartition_temperature
        return None

    @property
    def is_active(self) -> bool:
        return self._is_active
