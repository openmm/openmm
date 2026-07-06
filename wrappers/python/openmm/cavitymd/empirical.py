from pathlib import Path
from typing import Optional, Dict
import numpy as np

from .constants import Units


class EmpiricalTemperatureData:
    """Energy-to-temperature inversion using empirical calibration data.

    Loads equilibrium energy-vs-temperature data and fits analytical models
    whose inversions give fictive temperatures from instantaneous energies.

    For harmonic (vibrational) energy::

        phi(T) = a*T / (1 + b*T)
        T_v = E / (a - b*E)

    For LJ+Coulomb (structural) energy (Rosenfeld-Tarazona)::

        theta(T) = theta_0 + A*T^(3/5) / (1 + C*T^(3/5))
        T_s = numerical inversion via fsolve

    When no calibration data is available, falls back to direct equipartition::

        T_v = 4*V_bond / (N * k_B)

    Parameters
    ----------
    data_file_path : str or None
        Whitespace-separated file with columns: temperature, lj_hartree,
        coulombic_hartree, harmonic_hartree. None for direct-equipartition mode.
    energy_component : str
        'lj_coulombic', 'harmonic', or 'total_PE'.
    use_direct_harmonic : bool
        If True and component='harmonic', skip fitting; use T = 4E/(N*kB).
    """

    def __init__(
        self,
        data_file_path: Optional[str] = None,
        energy_component: str = "lj_coulombic",
        use_direct_harmonic: bool = False,
    ) -> None:
        self.data_file_path = Path(data_file_path) if data_file_path else None
        self.energy_component = energy_component
        self.use_direct_harmonic = use_direct_harmonic

        self.temperatures: Optional[np.ndarray] = None
        self.energies: Optional[np.ndarray] = None

        self.has_extended_harmonic_fit = False
        self.has_extended_t35_fit = False
        self.extended_harmonic_fit: Dict[str, float] = {}
        self.extended_t35_fit: Dict[str, float] = {}

        if self.data_file_path is not None:
            self._load_data()
            if not self.use_direct_harmonic or self.energy_component != "harmonic":
                if self.energy_component == "harmonic":
                    self.fit_extended_harmonic_function()
                else:
                    self.fit_extended_t35_function()

    def _load_data(self) -> None:
        if not self.data_file_path.exists():
            raise FileNotFoundError(f"Empirical data file not found: {self.data_file_path}")

        import pandas as pd

        data = pd.read_csv(self.data_file_path, sep=r"\s+", comment="#")

        if "temperature" not in data.columns:
            raise ValueError("temperature column not found in empirical data")
        self.temperatures = data["temperature"].values

        if self.energy_component == "lj_coulombic":
            if "lj_hartree" in data.columns and "coulombic_hartree" in data.columns:
                self.energies = data["lj_hartree"].values + data["coulombic_hartree"].values
            else:
                raise ValueError("lj_hartree and coulombic_hartree columns required")
        elif self.energy_component == "harmonic":
            if "harmonic_hartree" in data.columns:
                self.energies = data["harmonic_hartree"].values
            else:
                raise ValueError("harmonic_hartree column required")
        elif self.energy_component == "total_PE":
            if "total_potential_energy_hartree" in data.columns:
                self.energies = data["total_potential_energy_hartree"].values
            else:
                raise ValueError("total_potential_energy_hartree column required")
        else:
            raise ValueError(f"Unknown energy component: {self.energy_component}")

    def fit_extended_harmonic_function(self) -> None:
        """Fit E = a*T / (1 + b*T) to harmonic energy data."""
        from scipy.optimize import curve_fit

        def model(T, a, b):
            return a * T / (1 + b * T)

        a0 = self.energies.max() / self.temperatures.max()
        popt, _ = curve_fit(
            model, self.temperatures, self.energies,
            p0=[a0, 1e-3], bounds=([0, 0], [np.inf, np.inf]),
        )
        predicted = model(self.temperatures, *popt)
        ss_res = np.sum((self.energies - predicted) ** 2)
        ss_tot = np.sum((self.energies - np.mean(self.energies)) ** 2)
        self.extended_harmonic_fit = {"a": popt[0], "b": popt[1], "r2": 1 - ss_res / ss_tot}
        self.has_extended_harmonic_fit = True

    def fit_extended_t35_function(self) -> None:
        """Fit E = E0 + A*T^(3/5) / (1 + C*T^(3/5)) to LJ+Coulomb data."""
        from scipy.optimize import curve_fit

        def model(T, e0, a, c):
            t_pow = T ** (3.0 / 5.0)
            return e0 + a * t_pow / (1 + c * t_pow)

        e0_guess = self.energies.min()
        a_guess = (self.energies.max() - e0_guess) / (self.temperatures.max() ** 0.6)
        popt, _ = curve_fit(
            model, self.temperatures, self.energies,
            p0=[e0_guess, a_guess, 1e-3],
            bounds=(
                [self.energies.min() - 10, -np.inf, 0],
                [self.energies.max() + 10, np.inf, np.inf],
            ),
        )
        predicted = model(self.temperatures, *popt)
        ss_res = np.sum((self.energies - predicted) ** 2)
        ss_tot = np.sum((self.energies - np.mean(self.energies)) ** 2)
        self.extended_t35_fit = {"e0": popt[0], "a": popt[1], "c": popt[2], "r2": 1 - ss_res / ss_tot}
        self.has_extended_t35_fit = True

    def calculate_temperature(
        self, energy_hartree: float, num_particles: Optional[int] = None
    ) -> float:
        """Invert energy -> temperature (K) using fitted model.

        Parameters
        ----------
        energy_hartree : float
            Instantaneous potential energy in Hartree.
        num_particles : int or None
            Required for direct harmonic calculation.
        """
        if self.use_direct_harmonic and self.energy_component == "harmonic":
            if num_particles is None:
                raise ValueError("num_particles required for direct harmonic calculation")
            if energy_hartree <= 0:
                return 0.0
            return 4.0 * energy_hartree / (num_particles * Units.KB_HARTREE_PER_K)

        if self.has_extended_harmonic_fit and self.energy_component == "harmonic":
            a = self.extended_harmonic_fit["a"]
            b = self.extended_harmonic_fit["b"]
            denom = a - b * energy_hartree
            if energy_hartree <= 0 or denom <= 0:
                return 0.0
            return max(energy_hartree / denom, 0.0)

        if self.has_extended_t35_fit and self.energy_component != "harmonic":
            from scipy.optimize import brentq

            e0 = self.extended_t35_fit["e0"]
            a = self.extended_t35_fit["a"]
            c = self.extended_t35_fit["c"]
            if energy_hartree <= e0:
                return 0.0

            def residual(T):
                t_pow = T ** (3.0 / 5.0)
                return e0 + a * t_pow / (1 + c * t_pow) - energy_hartree

            try:
                T_lo, T_hi = 0.1, 1e5
                if residual(T_lo) * residual(T_hi) > 0:
                    T_hi = max(((energy_hartree - e0) / a) ** (5.0 / 3.0) * 10, 1e4)
                return max(brentq(residual, T_lo, T_hi), 0.0)
            except Exception:
                pass

        if self.temperatures is not None and self.energies is not None:
            return max(float(np.interp(energy_hartree, self.energies, self.temperatures)), 0.0)

        return 0.0

    @staticmethod
    def direct_harmonic_temperature(energy_hartree: float, num_particles: int) -> float:
        """T_v = 4*V_bond / (N * k_B) — no calibration data needed."""
        if energy_hartree <= 0 or num_particles <= 0:
            return 0.0
        return 4.0 * energy_hartree / (num_particles * Units.KB_HARTREE_PER_K)

    @property
    def fit_params(self) -> Dict[str, Dict]:
        return {
            "extended_harmonic": self.extended_harmonic_fit,
            "extended_t35": self.extended_t35_fit,
        }

    def check_inversion_at_calibration_points(self, tol_K: float = 15.0) -> bool:
        """Invert fitted curves at calibration energies; warn if T_recovered differs."""
        if self.temperatures is None or self.energies is None:
            return True

        issues = []
        for T_cal, E_cal in zip(self.temperatures, self.energies):
            T_inv = self.calculate_temperature(float(E_cal))
            if T_inv <= 0:
                issues.append(f"T_cal={T_cal:.1f} K: inversion returned {T_inv:.1f} K")
                continue
            if abs(T_inv - T_cal) > tol_K:
                issues.append(
                    f"T_cal={T_cal:.1f} K: inverted T={T_inv:.1f} K (Δ={T_inv - T_cal:+.1f} K)"
                )

        if issues:
            label = self.energy_component
            print(f"WARNING: {label} inversion check failed at {len(issues)} points:")
            for issue in issues[:8]:
                print(f"  {issue}")
            return False

        print(
            f"  {self.energy_component} inversion check passed "
            f"({len(self.temperatures)} calibration points, tol={tol_K:.0f} K)"
        )
        return True
