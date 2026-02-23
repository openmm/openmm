#!/usr/bin/env python3
"""
Compute IR spectrum from CACE dipole trajectory.
Adapted from compute_ir_spectrum_rpmd.py
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.ndimage import gaussian_filter1d
import sys
from pathlib import Path

def compute_ir_spectrum(npz_file, cutoff_ps=0.0, max_freq_cm=6000, output_prefix=None):
    """
    Compute IR spectrum from dipole trajectory.
    """
    print("=" * 60)
    print("IR Spectrum Calculation from CACE Dipoles")
    print("=" * 60)
    
    # Load data
    print(f"\nLoading data from: {npz_file}")
    data = np.load(npz_file, allow_pickle=True)
    time_ps = data['time_ps']
    dipole_eA = data['dipole_eA']
    metadata = data['metadata'].item()
    
    print(f"  Data points: {len(time_ps)}")
    print(f"  Time range: {time_ps[0]:.2f} - {time_ps[-1]:.2f} ps")
    print(f"  Molecules: {metadata.get('num_molecules', 'unknown')}")
    
    # Warn about short simulation times
    sim_time = time_ps[-1] - time_ps[0]
    if sim_time < 10.0:
        print(f"\n⚠ WARNING: Simulation time ({sim_time:.2f} ps) is very short for IR spectra.")
        print(f"  Recommended: at least 10-50 ps for smooth, converged spectra.")
        print(f"  Current spectrum may be noisy. Consider running a longer simulation.")
    
    # Apply time cutoff
    idx_start = np.where(time_ps >= cutoff_ps)[0]
    idx_start = idx_start[0] if len(idx_start) > 0 else 0
    
    time_cut = time_ps[idx_start:]
    dipole_cut = dipole_eA[idx_start:]
    
    # Timestep
    dt_ps = time_cut[1] - time_cut[0] if len(time_cut) > 1 else metadata.get('dt_ps', 0.0005)
    dt_fs = dt_ps * 1000.0
    
    # Compute autocorrelation
    print("\nComputing dipole autocorrelation function...")
    n_points = len(dipole_cut)
    
    # Compute total dipole autocorrelation (sum over x, y, z components)
    dipole_total = np.sum(dipole_cut**2, axis=1)
    dipole_mean = dipole_total.mean()
    dipole_centered = dipole_total - dipole_mean
    
    # Use scipy's correlate for better autocorrelation
    # Full autocorrelation
    autocorr_full = signal.correlate(dipole_centered, dipole_centered, mode='full', method='fft')
    autocorr = autocorr_full[len(autocorr_full)//2:]
    autocorr = autocorr[:n_points]
    
    # Normalize by the first value (variance)
    if autocorr[0] > 0:
        autocorr_norm = autocorr / autocorr[0]
    else:
        autocorr_norm = autocorr
    
    # Apply window function BEFORE normalization for spectrum
    # Use a longer window to reduce high-frequency noise
    window = signal.windows.hann(n_points)
    # Apply exponential decay to emphasize early times
    decay = np.exp(-time_cut / (time_cut[-1] * 0.5))
    window_weighted = window * decay
    autocorr_windowed = autocorr * window_weighted
    
    # Zero-pad for better frequency resolution
    n_pad = 8 * n_points  # More zero-padding for smoother spectrum
    autocorr_padded = np.zeros(n_pad)
    autocorr_padded[:n_points] = autocorr_windowed
    
    # FFT to frequency domain
    print("Computing FFT...")
    spectrum_raw = np.abs(np.fft.rfft(autocorr_padded))
    freqs_ps = np.fft.rfftfreq(n_pad, d=dt_ps)
    
    # Convert ps^-1 to cm^-1
    # 1 ps^-1 = 33.356 cm^-1
    freqs_cm = freqs_ps * 33.356
    
    # Apply smoothing to reduce noise
    # Use a Gaussian filter with appropriate width
    # Smooth with sigma proportional to frequency resolution
    freq_resolution = freqs_cm[1] - freqs_cm[0]
    sigma_smooth = max(5.0, freq_resolution * 2)  # At least 5 cm^-1 smoothing
    spectrum_smooth = gaussian_filter1d(spectrum_raw, sigma=sigma_smooth / freq_resolution)
    
    # Units and intensity
    # The IR intensity is proportional to ω * sinh(βω/2) * I(ω) or just ω^2 * I(ω)
    # We'll use a simple ω^2 * I(ω) scaling for visualization
    intensity = freqs_cm**2 * spectrum_smooth
    
    # Plotting
    if output_prefix is None:
        output_prefix = Path(npz_file).stem
        
    fig = plt.figure(figsize=(12, 10))
    
    # 1. Dipole moment time series
    ax1 = plt.subplot(3, 1, 1)
    time_rel = time_cut - time_cut[0]
    ax1.plot(time_rel, dipole_cut[:, 0], label='μₓ', alpha=0.7, linewidth=0.5)
    ax1.plot(time_rel, dipole_cut[:, 1], label='μᵧ', alpha=0.7, linewidth=0.5)
    ax1.plot(time_rel, dipole_cut[:, 2], label='μᵤ', alpha=0.7, linewidth=0.5)
    # Total magnitude
    dipole_mag = np.sqrt(np.sum(dipole_cut**2, axis=1))
    ax1.plot(time_rel, dipole_mag, 'k-', label='|μ|', alpha=0.9, linewidth=1.0)
    ax1.set_title(f"Dipole Moment Time Series ({metadata['model']})")
    ax1.set_xlabel("Time (ps)")
    ax1.set_ylabel("Dipole Moment (e·Å)")
    ax1.legend(loc='upper right', ncol=4)
    ax1.grid(True, alpha=0.3)
    
    # 2. Dipole autocorrelation
    ax2 = plt.subplot(3, 1, 2)
    ax2.plot(time_rel, autocorr_norm)
    ax2.set_title("Dipole Autocorrelation Function")
    ax2.set_xlabel("Time (ps)")
    ax2.set_ylabel("Normalized ACF")
    ax2.grid(True, alpha=0.3)
    
    # 3. IR spectrum
    ax3 = plt.subplot(3, 1, 3)
    mask = freqs_cm <= max_freq_cm
    ax3.plot(freqs_cm[mask], intensity[mask] / intensity[mask].max())
    ax3.set_title("IR Spectrum (Scaled Intensity)")
    ax3.set_xlabel("Frequency (cm⁻¹)")
    ax3.set_ylabel("Intensity (arb. units)")
    ax3.grid(True, alpha=0.3)
    
    # Mark characteristic peaks
    peaks = [1600, 3400]  # Bending and Stretching
    for p in peaks:
        ax3.axvline(x=p, color='r', linestyle='--', alpha=0.3)
        
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_spectrum.png", dpi=300)
    print(f"\nSaved spectrum plot: {output_prefix}_spectrum.png")
    
    # Save spectrum data
    np.savez(f"{output_prefix}_spectrum.npz", freqs_cm=freqs_cm, intensity=intensity)
    print(f"Saved spectrum data: {output_prefix}_spectrum.npz")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python compute_ir_spectrum_cace.py <dipole_file.npz>")
        sys.exit(1)
        
    compute_ir_spectrum(sys.argv[1])
