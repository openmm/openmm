#!/usr/bin/env python3
"""
Real-time IR Spectrum Calculator
=================================
Reads dipole moments from the real-time file and computes IR spectrum
while the simulation is still running.

This implementation matches tests/water_system/analyze_spectrum.py for consistency:
- Uses FFT-based autocorrelation (faster than direct method)
- Subtracts mean dipole before ACF (removes DC component)
- Normalizes ACF by number of samples
- Computes spectrum as ω² × |FT[ACF]| (not ω × |FT[ACF]|²)
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import time
from pathlib import Path

def compute_acf(dipole_vector_timeseries, subtract_mean=True):
    """
    Compute the autocorrelation function of the dipole vector using FFT method.
    
    For IR spectroscopy, we need the VECTOR dot product autocorrelation:
    C(t) = <μ(0)·μ(t)>
    
    This implementation matches tests/water_system/analyze_spectrum.py
    
    Parameters
    ----------
    dipole_vector_timeseries : ndarray (N, 3)
        Dipole moment trajectory (e·Å)
    subtract_mean : bool
        Whether to subtract mean (removes DC component)
        
    Returns
    -------
    C_mumu : ndarray (N,)
        Autocorrelation function
    """
    # Process each component
    dipole = dipole_vector_timeseries.copy()
    
    if subtract_mean:
        dipole = dipole - dipole.mean(axis=0)
    
    N = len(dipole)
    C_mumu = np.zeros(N)
    
    # Compute autocorrelation using FFT (fast method)
    # Sum over all 3 components (x, y, z)
    for alpha in range(3):
        mu_alpha = dipole[:, alpha]
        
        # FFT-based autocorrelation
        # Pad to avoid circular convolution artifacts
        nfft = 2 * N
        mu_fft = np.fft.fft(mu_alpha, n=nfft)
        acf = np.fft.ifft(mu_fft * np.conj(mu_fft))
        C_mumu += acf[:N].real
    
    # Normalize by number of samples
    norm = np.arange(N, 0, -1)
    C_mumu /= norm
    
    return C_mumu

def compute_ir_spectrum_realtime(dipole_file, output_prefix="realtime_ir", 
                                 max_freq_cm=4500, window_ps=None):
    """
    Compute IR spectrum from real-time dipole data.
    
    Args:
        dipole_file: Path to real-time dipole data file
        output_prefix: Prefix for output files
        max_freq_cm: Maximum frequency to plot (cm^-1)
        window_ps: If specified, only use last N picoseconds of data
    """
    
    # Read dipole data
    print(f"Reading dipole data from {dipole_file}...")
    try:
        data = np.loadtxt(dipole_file, comments='#')
        times_ps = data[:, 0]
        dipoles_eA = data[:, 1:4]  # 3D vector
    except Exception as e:
        print(f"Error reading file: {e}")
        return None
    
    n_points = len(times_ps)
    if n_points < 100:
        print(f"Not enough data points yet ({n_points}), need at least 100")
        return None
    
    # Apply time window if specified
    if window_ps is not None:
        max_time = times_ps[-1]
        min_time = max(0, max_time - window_ps)
        mask = times_ps >= min_time
        times_ps = times_ps[mask]
        dipoles_eA = dipoles_eA[mask]
        n_points = len(times_ps)
        print(f"  Using time window: {min_time:.2f} - {max_time:.2f} ps ({n_points} points)")
    
    dt_ps = times_ps[1] - times_ps[0]
    total_time_ps = times_ps[-1] - times_ps[0]
    
    print(f"  Data points: {n_points}")
    print(f"  Time step: {dt_ps*1000:.3f} fs")
    print(f"  Total time: {total_time_ps:.3f} ps")
    print(f"  Nyquist frequency: {1/(2*dt_ps*1e-12*2.998e10):.1f} cm⁻¹")
    
    # Compute dipole vector autocorrelation function
    print("Computing dipole autocorrelation function...")
    acf = compute_acf(dipoles_eA, subtract_mean=True)
    
    # Apply Hann window to reduce spectral leakage
    window = np.hanning(len(acf))
    acf_windowed = acf * window
    
    # Compute power spectrum via FFT
    print("Computing FFT...")
    fft_result = np.fft.rfft(acf_windowed)
    
    # Frequency axis
    freqs_hz = np.fft.rfftfreq(len(acf_windowed), d=dt_ps * 1e-12)
    freqs_cm = freqs_hz / (2.998e10)  # Convert Hz to cm^-1
    
    # IR intensity is proportional to ω² * |FT[ACF]|
    # Matches water_system
    omega = 2 * np.pi * freqs_hz
    ir_intensity = omega**2 * np.abs(fft_result)
    
    # Normalize
    ir_intensity /= np.max(ir_intensity[freqs_cm > 100])  # Normalize excluding low freq
    
    # Restrict to plotting range
    mask = (freqs_cm >= 0) & (freqs_cm <= max_freq_cm)
    freqs_plot = freqs_cm[mask]
    intensity_plot = ir_intensity[mask]
    
    # Save spectrum data
    spectrum_file = f"{output_prefix}_spectrum.npz"
    np.savez(spectrum_file,
             frequency_cm=freqs_plot,
             intensity=intensity_plot,
             acf=acf,
             time_ps=times_ps[:len(acf)],
             metadata={'dt_ps': dt_ps, 'n_points': n_points, 'total_time_ps': total_time_ps})
    
    # Plot
    fig, axes = plt.subplots(3, 1, figsize=(12, 14))
    
    # Dipole moment time series
    ax0 = axes[0]
    dipole_magnitude = np.linalg.norm(dipoles_eA, axis=1)
    time_array = times_ps[:len(dipole_magnitude)]
    
    # Plot magnitude
    ax0.plot(time_array, dipole_magnitude, 'g-', linewidth=1.0, alpha=0.8, label='|μ|')
    
    # Also plot components with reduced opacity
    ax0.plot(time_array, dipoles_eA[:, 0], 'r-', linewidth=0.5, alpha=0.4, label='μx')
    ax0.plot(time_array, dipoles_eA[:, 1], 'b-', linewidth=0.5, alpha=0.4, label='μy')
    ax0.plot(time_array, dipoles_eA[:, 2], 'm-', linewidth=0.5, alpha=0.4, label='μz')
    
    ax0.set_xlabel('Time (ps)', fontsize=12)
    ax0.set_ylabel('Dipole Moment (e·Å)', fontsize=12)
    ax0.set_title(f'Dipole Moment Time Series ({n_points} points, {total_time_ps:.2f} ps)', fontsize=14)
    ax0.grid(True, alpha=0.3)
    ax0.legend(loc='upper right', fontsize=10)
    
    # Calculate and show statistics
    mean_mag = np.mean(dipole_magnitude)
    std_mag = np.std(dipole_magnitude)
    ax0.axhline(mean_mag, color='k', linestyle='--', alpha=0.5, linewidth=1)
    ax0.text(0.02, 0.98, f'Mean: {mean_mag:.2f} e·Å\nStd: {std_mag:.2f} e·Å', 
             transform=ax0.transAxes, verticalalignment='top', 
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5), fontsize=10)
    
    # ACF plot
    ax1 = axes[1]
    time_acf_ps = np.arange(len(acf)) * dt_ps
    ax1.plot(time_acf_ps, acf / acf[0], 'b-', linewidth=1.5, label='Dipole ACF')
    ax1.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax1.set_xlabel('Time (ps)', fontsize=12)
    ax1.set_ylabel('Normalized ACF', fontsize=12)
    ax1.set_title(f'Dipole Autocorrelation Function', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # IR spectrum plot
    ax2 = axes[2]
    ax2.plot(freqs_plot, intensity_plot, 'r-', linewidth=1.5)
    ax2.set_xlabel('Frequency (cm⁻¹)', fontsize=12)
    ax2.set_ylabel('IR Intensity (normalized)', fontsize=12)
    ax2.set_title(f'Infrared Spectrum (Real-time)', fontsize=14)
    ax2.set_xlim(0, max_freq_cm)
    ax2.grid(True, alpha=0.3)
    
    # Mark key water bands
    water_bands = [
        (1595, 'H-O-H bend'),
        (3200, 'O-H stretch (symmetric)'),
        (3400, 'O-H stretch (asymmetric)')
    ]
    for freq, label in water_bands:
        if freq <= max_freq_cm:
            ax2.axvline(freq, color='gray', linestyle='--', alpha=0.5)
            ax2.text(freq, ax2.get_ylim()[1]*0.95, label, 
                    rotation=90, verticalalignment='top', fontsize=9, alpha=0.7)
    
    plt.tight_layout()
    plot_file = f"{output_prefix}_spectrum.png"
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"\nSpectrum saved: {plot_file}")
    print(f"Data saved: {spectrum_file}")
    
    return freqs_plot, intensity_plot

def monitor_realtime(dipole_file, output_prefix="realtime_ir", 
                     update_interval=30, max_freq_cm=4500):
    """
    Monitor the dipole file and update the IR spectrum periodically.
    
    Args:
        dipole_file: Path to real-time dipole data file
        output_prefix: Prefix for output files
        update_interval: Update interval in seconds
        max_freq_cm: Maximum frequency to plot (cm^-1)
    """
    print("=" * 80)
    print("Real-time IR Spectrum Monitor")
    print("=" * 80)
    print(f"Monitoring: {dipole_file}")
    print(f"Update interval: {update_interval} seconds")
    print(f"Press Ctrl+C to stop")
    print("=" * 80)
    
    last_size = 0
    iteration = 0
    
    while True:
        try:
            # Check if file has grown
            if Path(dipole_file).exists():
                current_size = Path(dipole_file).stat().st_size
                
                if current_size > last_size:
                    iteration += 1
                    print(f"\n[{time.strftime('%H:%M:%S')}] Update #{iteration}")
                    
                    result = compute_ir_spectrum_realtime(
                        dipole_file, 
                        output_prefix=f"{output_prefix}_iter{iteration:03d}",
                        max_freq_cm=max_freq_cm
                    )
                    
                    last_size = current_size
                else:
                    print(f"[{time.strftime('%H:%M:%S')}] Waiting for new data...", end='\r')
            else:
                print(f"[{time.strftime('%H:%M:%S')}] Waiting for file to be created...", end='\r')
            
            time.sleep(update_interval)
            
        except KeyboardInterrupt:
            print("\n\nMonitoring stopped by user.")
            break
        except Exception as e:
            print(f"\nError: {e}")
            time.sleep(update_interval)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Real-time IR spectrum calculator')
    parser.add_argument('dipole_file', type=str, help='Path to real-time dipole data file')
    parser.add_argument('--output', type=str, default='realtime_ir', help='Output prefix')
    parser.add_argument('--monitor', action='store_true', help='Continuously monitor and update')
    parser.add_argument('--interval', type=int, default=30, help='Update interval in seconds (for monitor mode)')
    parser.add_argument('--max-freq', type=float, default=4500, help='Maximum frequency (cm^-1)')
    parser.add_argument('--window', type=float, default=None, help='Use only last N ps of data')
    
    args = parser.parse_args()
    
    if args.monitor:
        monitor_realtime(args.dipole_file, args.output, args.interval, args.max_freq)
    else:
        compute_ir_spectrum_realtime(args.dipole_file, args.output, args.max_freq, args.window)
