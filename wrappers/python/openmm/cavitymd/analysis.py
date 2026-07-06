"""Relaxation-time and Tool-Narayanaswamy analysis for C2F cavity-MD."""

from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np

k_B = 8.617333e-5  # Boltzmann constant in eV/K


class RelaxationTimeModel:
    """Temperature-dependent tau_s(T) from dual-regime (Arrhenius + parabolic) fit."""

    def __init__(self, data_file_path: Optional[str] = None) -> None:
        self.data_file_path = data_file_path
        self.is_fitted = False
        self.T_onset: Optional[float] = None
        self.fit_results: Dict = {}

        if data_file_path and Path(data_file_path).exists():
            self._load_and_fit_data()

    def _load_and_fit_data(self) -> None:
        try:
            data = np.loadtxt(self.data_file_path, skiprows=3, usecols=(0, 2))
            temperatures = data[:, 0]
            relaxation_times = data[:, 1]
            self.T_onset, self.fit_results = self._find_onset_temperature(
                temperatures, relaxation_times
            )
            self.is_fitted = not np.isnan(self.T_onset)
        except Exception as exc:
            print(f"RelaxationTimeModel: failed to load {self.data_file_path}: {exc}")

    @staticmethod
    def _fit_arrhenius(beta: np.ndarray, log_tau: np.ndarray) -> Tuple[float, float, float]:
        from scipy.optimize import curve_fit

        def model(b, Ea, beta_0, ln_tau_0):
            return ln_tau_0 + Ea * (b - beta_0)

        coeffs = np.polyfit(beta, log_tau, 1)
        beta_0_guess = float(np.mean(beta))
        ln_tau_0_guess = float(coeffs[1] + coeffs[0] * beta_0_guess)
        popt, _ = curve_fit(
            model, beta, log_tau,
            p0=[coeffs[0], beta_0_guess, ln_tau_0_guess],
        )
        return float(popt[0]), float(popt[1]), float(popt[2])

    @staticmethod
    def _fit_parabolic(
        beta: np.ndarray, log_tau: np.ndarray, beta_0: float, arr: Tuple[float, float]
    ) -> Tuple[float, float, float]:
        from scipy.optimize import curve_fit

        Ea_arr, ln_tau_0_arr = arr

        def model(b, Ea, J):
            d = b - beta_0
            return ln_tau_0_arr + Ea * d + J ** 2 * d ** 2

        delta = beta - beta_0
        residuals = log_tau - (ln_tau_0_arr + Ea_arr * delta)
        J_guess = float(np.mean(residuals / (delta ** 2 + 1e-10)))
        popt, _ = curve_fit(model, beta, log_tau, p0=[Ea_arr, max(J_guess, 1e-6)])
        return float(popt[0]), float(popt[1]), ln_tau_0_arr

    def _find_onset_temperature(
        self, temperatures: np.ndarray, relaxation_times: np.ndarray
    ) -> Tuple[float, Dict]:
        valid = ~np.isnan(relaxation_times)
        T = temperatures[valid]
        tau = relaxation_times[valid]
        if len(T) < 6:
            return np.nan, {}

        order = np.argsort(T)
        T = T[order]
        tau = tau[order]
        beta = 1.0 / (k_B * T)
        log_tau = np.log(tau)

        min_pts = 5
        best_r2 = -np.inf
        best_idx = None
        best_fits: Dict = {}

        for split in range(min_pts, len(T) - min_pts + 1):
            low_b, low_lt = beta[:split], log_tau[:split]
            high_b, high_lt = beta[split:], log_tau[split:]
            try:
                Ea, beta_0, ln_tau_0 = self._fit_arrhenius(high_b, high_lt)
                pred_h = ln_tau_0 + Ea * (high_b - beta_0)
                r2_h = 1 - np.sum((high_lt - pred_h) ** 2) / (
                    np.sum((high_lt - np.mean(high_lt)) ** 2) + 1e-10
                )
                Ea_p, J_p, ln_tau_0_p = self._fit_parabolic(
                    low_b, low_lt, beta_0, (Ea, ln_tau_0)
                )
                d = low_b - beta_0
                pred_l = ln_tau_0_p + Ea_p * d + J_p ** 2 * d ** 2
                r2_l = 1 - np.sum((low_lt - pred_l) ** 2) / (
                    np.sum((low_lt - np.mean(low_lt)) ** 2) + 1e-10
                )
                n_h, n_l = len(high_b), len(low_b)
                total_r2 = (n_h * r2_h + n_l * r2_l) / (n_h + n_l)
                if J_p > 0 and total_r2 > best_r2:
                    best_r2 = total_r2
                    best_idx = split
                    best_fits = {
                        "beta_0": beta_0,
                        "arrhenius": {"Ea": Ea, "ln_tau_0": ln_tau_0, "r2": r2_h},
                        "parabolic": {"Ea": Ea_p, "J": J_p, "ln_tau_0": ln_tau_0_p, "r2": r2_l},
                    }
            except Exception:
                continue

        if best_idx is None:
            return np.nan, {}
        return float(T[best_idx]), best_fits

    def get_relaxation_time(self, temperature_K: float) -> float:
        if temperature_K <= 0:
            return 1.0
        if not self.is_fitted:
            return 100.0 * np.exp(0.1 / (k_B * temperature_K))

        beta = 1.0 / (k_B * temperature_K)
        beta_0 = self.fit_results["beta_0"]
        if temperature_K > self.T_onset:
            p = self.fit_results["arrhenius"]
            ln_tau = p["ln_tau_0"] + p["Ea"] * (beta - beta_0)
        else:
            p = self.fit_results["parabolic"]
            d = beta - beta_0
            ln_tau = p["ln_tau_0"] + p["Ea"] * d + p["J"] ** 2 * d ** 2
        return float(np.exp(ln_tau))


class ToolNarayanaswamy:
    """Material-time reconstruction and TN-model integration."""

    def __init__(
        self,
        relaxation_model: Optional[RelaxationTimeModel] = None,
        beta: float = 0.55,
        smoothness_alpha: float = 1.0,
    ) -> None:
        self.relaxation_model = relaxation_model
        self.beta = float(beta)
        self.smoothness_alpha = float(smoothness_alpha)

    @staticmethod
    def stretched_exponential(h: np.ndarray, beta: float = 0.55) -> np.ndarray:
        return np.exp(-np.power(np.maximum(h, 0.0), beta))

    @staticmethod
    def _second_difference_matrix(n: int) -> np.ndarray:
        D = np.zeros((max(n - 2, 0), n))
        for i in range(n - 2):
            D[i, i] = 1.0
            D[i, i + 1] = -2.0
            D[i, i + 2] = 1.0
        return D

    def reconstruct_material_time(
        self,
        waiting_times_ps: np.ndarray,
        relaxation_times_ps: np.ndarray,
        time_grid_ps: Optional[np.ndarray] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Regularized LS reconstruction of h(t) from MTTI constraints."""
        tw = np.asarray(waiting_times_ps, dtype=float)
        tau = np.asarray(relaxation_times_ps, dtype=float)
        if time_grid_ps is None:
            t_max = float(np.max(tw + tau) * 1.1)
            time_grid_ps = np.linspace(0.0, t_max, max(len(tw) * 4, 20))

        t_grid = np.asarray(time_grid_ps, dtype=float)
        n = len(t_grid)
        A = np.zeros((len(tw), n))
        b = np.ones(len(tw))

        for r, (t_w, t_r) in enumerate(zip(tw, tau)):
            t_lo, t_hi = t_w, t_w + t_r
            for m in range(n):
                t_m = t_grid[m]
                t_lo_m = t_grid[m - 1] if m > 0 else t_grid[m]
                t_hi_m = t_grid[m + 1] if m + 1 < n else t_grid[m]
                A[r, m] += self._hat(t_hi, t_lo, t_hi_m, t_lo_m)
                A[r, m] -= self._hat(t_lo, t_lo, t_hi_m, t_lo_m)

        D = self._second_difference_matrix(n)
        if D.size == 0:
            h = np.linalg.lstsq(A, b, rcond=None)[0]
        else:
            lhs = A.T @ A + self.smoothness_alpha * (D.T @ D)
            rhs = A.T @ b
            h = np.linalg.solve(lhs, rhs)
        h = np.maximum.accumulate(np.maximum(h, 0.0))
        return t_grid, h

    @staticmethod
    def _hat(t: float, t_lo: float, t_hi: float, t_lo_m: float) -> float:
        if t_hi <= t_lo_m or t <= t_lo_m:
            return 0.0
        if t_lo >= t_hi:
            return 0.0
        denom = t_hi - t_lo_m
        if denom <= 0:
            return 0.0
        t_clip = min(t, t_hi)
        if t_clip <= t_lo_m:
            return 0.0
        if t_clip >= t_hi:
            return 1.0
        return (t_clip - t_lo_m) / denom

    def integrate_tn(
        self,
        times_ps: np.ndarray,
        structural_temperatures_K: np.ndarray,
    ) -> np.ndarray:
        """Integrate dh/dt = 1/tau_s,eq[T_s(t)] by trapezoidal rule."""
        t = np.asarray(times_ps, dtype=float)
        T_s = np.asarray(structural_temperatures_K, dtype=float)
        if len(t) < 2:
            return np.zeros_like(t)

        h = np.zeros(len(t))
        for i in range(1, len(t)):
            dt = t[i] - t[i - 1]
            if self.relaxation_model is not None:
                tau_prev = self.relaxation_model.get_relaxation_time(T_s[i - 1])
                tau_curr = self.relaxation_model.get_relaxation_time(T_s[i])
            else:
                tau_prev = tau_curr = 100.0
            rate = 0.5 * (1.0 / max(tau_prev, 1e-12) + 1.0 / max(tau_curr, 1e-12))
            h[i] = h[i - 1] + dt * rate
        return h

    def collapse_isf(
        self,
        correlation_times_ps: np.ndarray,
        waiting_times_ps: np.ndarray,
        h_at_times: np.ndarray,
        time_grid_ps: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Reparameterize ISF onto material time using interpolated h(t)."""
        t_corr = np.asarray(correlation_times_ps, dtype=float)
        t_w = np.asarray(waiting_times_ps, dtype=float)
        h_grid = np.interp(t_w + t_corr, time_grid_ps, h_at_times)
        h_w = np.interp(t_w, time_grid_ps, h_at_times)
        h_diff = np.maximum(h_grid - h_w, 0.0)
        phi = self.stretched_exponential(h_diff, self.beta)
        return h_diff, phi
