"""
F(k,t) intermediate scattering function tracker for OpenMM cavity simulations.

Mirrors cav-hoomd's FieldAutocorrelationTracker: ρ_k(t) = Σ_j exp(i k·r_j),
F(k,t) = Re⟨ρ_k(t) ρ*_k(t_0)⟩ with k-averaging over a Fibonacci sphere,
multiple reference times t_w, and per-reference output files.

Positions in nm; wavevectors in nm⁻¹. Caller passes molecular positions only
(exclude cavity particle). Use wrapped positions (cav-hoomd parity) so that
at 0 K F(k,t) stays constant; unwrapped positions can drift from round-off
and cause spurious decorrelation.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple

import numpy as np


def fibonacci_sphere(samples: int) -> np.ndarray:
    """
    Generate unit vectors on a sphere using the Fibonacci spiral.

    Parameters
    ----------
    samples : int
        Number of points on the sphere.

    Returns
    -------
    np.ndarray
        Shape (samples, 3), each row a unit vector [x, y, z].
    """
    if samples < 1:
        raise ValueError("samples must be at least 1")
    phi = np.pi * (3.0 - np.sqrt(5.0))  # Golden angle in radians
    points = []
    n = max(1, samples - 1)
    for i in range(samples):
        y = 1 - (i / float(n)) * 2 if n > 0 else 0.0
        radius = np.sqrt(1 - y * y) if abs(y) < 1 else 0.0
        theta = phi * i
        x = np.cos(theta) * radius
        z = np.sin(theta) * radius
        points.append([x, y, z])
    return np.array(points, dtype=np.float64)


def compute_rhok(
    positions_nm: np.ndarray,
    wavevectors_nm_inv: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute density modes ρ_k = Σ_j exp(i k·r_j) for all wavevectors.

    Parameters
    ----------
    positions_nm : np.ndarray
        Shape (N, 3), positions in nm (molecular particles only).
    wavevectors_nm_inv : np.ndarray
        Shape (M, 3), wavevectors in nm⁻¹.

    Returns
    -------
    rhok_real, rhok_imag : np.ndarray
        Each shape (M,), real and imaginary parts of ρ_k.
    """
    # k_dot_r: (M, N)
    k_dot_r = np.dot(wavevectors_nm_inv, positions_nm.T)
    rhok_complex = np.sum(np.exp(1j * k_dot_r), axis=1)
    return np.real(rhok_complex), np.imag(rhok_complex)


def compute_fkt(
    rhok0_real: np.ndarray,
    rhok0_imag: np.ndarray,
    rhok_real: np.ndarray,
    rhok_imag: np.ndarray,
) -> float:
    """
    Compute F(k,t) = Re⟨ρ_k(t) ρ*_k(t_0)⟩ averaged over k-vectors.

    ρ_k(t)·ρ*_k(t_0) has real part: real(t)·real(t_0) + imag(t)·imag(t_0).

    Parameters
    ----------
    rhok0_real, rhok0_imag : np.ndarray
        Real and imaginary parts of ρ_k at reference time.
    rhok_real, rhok_imag : np.ndarray
        Real and imaginary parts of ρ_k at current time.

    Returns
    -------
    float
        k-averaged F(k,t).
    """
    correlations = rhok_real * rhok0_real + rhok_imag * rhok0_imag
    return float(np.mean(correlations))


class FKTTracker:
    """
    Tracks F(k,t) during a simulation with multiple reference times t_w.

    Writes one file per reference: {output_prefix}_fkt_ref_{idx:03d}.txt
    with header "# Reference time: {t_w} ps" and columns lag_time_ps, F(k,t).
    """

    def __init__(
        self,
        kmag_nm_inv: float,
        num_wavevectors: int,
        reference_interval_ps: float,
        max_references: int,
        output_period_ps: float,
        output_prefix: str,
        log_time_spacing: bool = False,
        min_log_ps: Optional[float] = None,
        max_log_ps: Optional[float] = None,
        log_num_points: int = 50,
    ) -> None:
        """
        Initialize the F(k,t) tracker.

        Parameters
        ----------
        kmag_nm_inv : float
            |k| in nm⁻¹ (e.g. 113.4 for 6.0 a.u. in paper).
        num_wavevectors : int
            Number of k-vectors on Fibonacci sphere.
        reference_interval_ps : float
            Interval between new reference times t_w (ps).
        max_references : int
            Maximum number of reference frames to keep (FIFO).
        output_period_ps : float
            How often to write F(k,t) points in linear mode (ps).
        output_prefix : str
            Prefix for output files (e.g. "fkt_lambda0.07").
        log_time_spacing : bool
            If True, use log-spaced target lags (requires min_log_ps, max_log_ps).
        min_log_ps, max_log_ps : float, optional
            Log-spacing range in ps. Required if log_time_spacing=True.
        log_num_points : int
            Number of log-spaced points. Used if log_time_spacing=True.
        """
        self.kmag_nm_inv = kmag_nm_inv
        self.num_wavevectors = num_wavevectors
        self.reference_interval_ps = reference_interval_ps
        self.max_references = max_references
        self.output_period_ps = output_period_ps
        self.output_prefix = str(output_prefix)
        self.log_time_spacing = log_time_spacing
        self.min_log_ps = min_log_ps
        self.max_log_ps = max_log_ps
        self.log_num_points = log_num_points

        if log_time_spacing:
            if min_log_ps is None or max_log_ps is None:
                raise ValueError(
                    "min_log_ps and max_log_ps required when log_time_spacing=True"
                )
            if min_log_ps >= max_log_ps or min_log_ps <= 0:
                raise ValueError("Need 0 < min_log_ps < max_log_ps")
            self.log_time_points_ps = np.logspace(
                np.log10(min_log_ps),
                np.log10(max_log_ps),
                log_num_points,
                dtype=np.float64,
            )
            self._written_log_points: dict[int, set[int]] = {}
        else:
            self.log_time_points_ps = np.array([])
            self._written_log_points = {}

        # Wavevectors: unit vectors * kmag, shape (num_wavevectors, 3)
        unit_vectors = fibonacci_sphere(num_wavevectors)
        self.wavevectors = unit_vectors * kmag_nm_inv

        self.references: list[dict] = []
        self.last_reference_time_ps: Optional[float] = None
        self.last_output_time_ps: Optional[float] = None
        self._next_ref_file_idx: int = 0  # Monotonic; stable file name after FIFO drop

        # Ensure output directory exists if prefix contains a path
        prefix_path = Path(self.output_prefix)
        if prefix_path.parent != Path(".") and not prefix_path.parent.exists():
            prefix_path.parent.mkdir(parents=True, exist_ok=True)

    def _should_add_reference(self, current_time_ps: float) -> bool:
        if self.last_reference_time_ps is None:
            return True
        return (current_time_ps - self.last_reference_time_ps) >= self.reference_interval_ps

    def _should_output(self, current_time_ps: float) -> bool:
        if self.last_output_time_ps is None:
            return True
        return (current_time_ps - self.last_output_time_ps) >= self.output_period_ps

    def _write_row(
        self, ref_file_idx: int, lag_time_ps: float, fkt_value: float
    ) -> None:
        """Append one row to the file for reference ref_file_idx."""
        path = f"{self.output_prefix}_fkt_ref_{ref_file_idx:03d}.txt"
        with open(path, "a", encoding="utf-8") as f:
            f.write(f"{lag_time_ps:.6f}\t{fkt_value:.8f}\n")

    def _create_ref_file(self, ref_file_idx: int, time_ps: float) -> None:
        """Create and write header for a new reference file."""
        path = f"{self.output_prefix}_fkt_ref_{ref_file_idx:03d}.txt"
        with open(path, "w", encoding="utf-8") as f:
            f.write("# F(k,t) correlation function\n")
            f.write(f"# Reference time: {time_ps:.3f} ps\n")
            if self.log_time_spacing:
                f.write("# Using logarithmic time spacing\n")
                f.write(
                    f"# Time range: {self.min_log_ps:.3f} - {self.max_log_ps:.1f} ps\n"
                )
            f.write("# lag_time_ps\tF(k,t)\n")

    def _output_for_ref(
        self,
        ref: dict,
        current_time_ps: float,
        rhok_real: np.ndarray,
        rhok_imag: np.ndarray,
    ) -> None:
        """Compute F(k,t) for this reference and write if needed."""
        file_idx = ref["ref_file_idx"]
        lag_time_ps = current_time_ps - ref["time_ps"]
        fkt_value = compute_fkt(
            ref["rhok_real"], ref["rhok_imag"], rhok_real, rhok_imag
        )

        if self.log_time_spacing:
            if file_idx not in self._written_log_points:
                self._written_log_points[file_idx] = set()
            time_diffs = np.abs(self.log_time_points_ps - lag_time_ps)
            closest_idx = int(np.argmin(time_diffs))
            closest_time = self.log_time_points_ps[closest_idx]
            min_diff = time_diffs[closest_idx]
            if min_diff < 0.005 and closest_idx not in self._written_log_points[file_idx]:
                self._written_log_points[file_idx].add(closest_idx)
                self._write_row(file_idx, closest_time, fkt_value)
        else:
            self._write_row(file_idx, lag_time_ps, fkt_value)

    def update(self, current_time_ps: float, positions_molecular_nm: np.ndarray) -> None:
        """
        Update the tracker with current time and molecular positions.

        Adds a new reference if the interval has elapsed; writes F(k,t) for
        all references when the output period has elapsed (linear) or when
        a log-spaced target is hit.

        Parameters
        ----------
        current_time_ps : float
            Current simulation time in ps.
        positions_molecular_nm : np.ndarray
            Shape (N, 3), molecular positions in nm (cavity excluded).
        """
        positions = np.asarray(positions_molecular_nm, dtype=np.float64)
        if positions.ndim != 2 or positions.shape[1] != 3:
            raise ValueError("positions_molecular_nm must be (N, 3)")

        rhok_real, rhok_imag = compute_rhok(positions, self.wavevectors)

        if self._should_add_reference(current_time_ps):
            file_idx = self._next_ref_file_idx
            self._next_ref_file_idx += 1
            ref = {
                "time_ps": current_time_ps,
                "rhok_real": rhok_real.copy(),
                "rhok_imag": rhok_imag.copy(),
                "ref_file_idx": file_idx,
            }
            self.references.append(ref)
            self.last_reference_time_ps = current_time_ps
            self._create_ref_file(file_idx, current_time_ps)
            if len(self.references) > self.max_references:
                self.references.pop(0)

        if self.log_time_spacing:
            for ref in self.references:
                self._output_for_ref(ref, current_time_ps, rhok_real, rhok_imag)
        else:
            if self._should_output(current_time_ps):
                for ref in self.references:
                    self._output_for_ref(ref, current_time_ps, rhok_real, rhok_imag)
                self.last_output_time_ps = current_time_ps

    def finalize(self) -> None:
        """
        Flush or close any resources. Optional hook for log-spaced post-processing.

        Currently a no-op; linear output is written inline in update().
        """
        pass
