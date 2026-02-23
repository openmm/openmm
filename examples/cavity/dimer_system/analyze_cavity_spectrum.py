#!/usr/bin/env python3
"""
Calculate Cavity Position (q) Spectrum using Maximum Entropy Method
Reference: https://maximum-entropy-spectrum.readthedocs.io/en/latest/

This script analyzes the cavity mode position q(t) to compute its power spectrum,
which shows the cavity resonance and any driving frequencies.
"""

import numpy as np
from scipy import signal
import memspectrum
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import sys
import os

print("=" * 60)
print("Cavity Position Spectrum Calculation using FFT of Autocorrelation")
print("=" * 60)

# Get filename from command line
if len(sys.argv) < 2:
    print("\nUsage: python analyze_cavity_spectrum.py <trajectory_file.npz>")
    print("\nSearching for cavity_driven_dipole_trajectory.npz files...")
    import glob
    npz_files = glob.glob("cavity_driven_dipole_trajectory.npz")
    if npz_files:
        # Sort by modification time, get most recent
        npz_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
        input_file = npz_files[0]
        print(f"  Found {len(npz_files)} file(s), using most recent: {input_file}")
    else:
        print("  No cavity_driven_dipole_trajectory.npz files found!")
        sys.exit(1)
else:
    input_file = sys.argv[1]

# Load the data
print(f"\n1. Loading cavity position trajectory data from: {input_file}")
data = np.load(input_file, allow_pickle=True)
time = data['time_ps']

# Check if cavity_q_nm exists
if 'cavity_q_nm' not in data:
    print("   ERROR: 'cavity_q_nm' not found in file!")
    print(f"   Available keys: {list(data.keys())}")
    sys.exit(1)

cavity_q = data['cavity_q_nm']  # Shape: (N, 3) in nm

print(f"   Total data points: {len(time)}")
print(f"   Time range: {time[0]:.2f} - {time[-1]:.2f} ps")
print(f"   Cavity position shape: {cavity_q.shape}")

# Cut off equilibration (first 20 ps or 10% of data, whichever is smaller)
print("\n2. Removing equilibration period...")
total_time = time[-1] - time[0]
cutoff_time = min(20.0, total_time * 0.1)  # 20 ps or 10% of total time
print(f"   Total simulation time: {total_time:.2f} ps")
print(f"   Cutoff time: {cutoff_time:.2f} ps")

if total_time > cutoff_time:
    idx_cutoff = np.where(time >= cutoff_time)[0][0]
    time_cut = time[idx_cutoff:]
    cavity_q_cut = cavity_q[idx_cutoff:]
    print(f"   Data points after cutoff: {len(time_cut)}")
    print(f"   Time range: {time_cut[0]:.2f} - {time_cut[-1]:.2f} ps")
else:
    print(f"   ⚠ Warning: Simulation too short ({total_time:.1f} ps), using all data")
    time_cut = time
    cavity_q_cut = cavity_q
    print(f"   Data points: {len(time_cut)}")
    print(f"   Time range: {time_cut[0]:.2f} - {time_cut[-1]:.2f} ps")

if len(time_cut) < 100:
    print(f"\n   ERROR: Not enough data points ({len(time_cut)}) for spectral analysis!")
    print(f"   Need at least 100 points. Please run simulation longer.")
    sys.exit(1)

# Compute cavity position with mean subtracted (δq = q - <q>)
print("\n3. Computing cavity position with mean subtracted (δq = q - <q>)...")
# Use x,y components (cavity mode is in xy plane)
cavity_q_xy = cavity_q_cut[:, :2]  # Only x,y components
cavity_q_mean = np.mean(cavity_q_xy, axis=0)
cavity_q_centered = cavity_q_xy - cavity_q_mean

print(f"   Mean cavity position (x,y): ({cavity_q_mean[0]:.6f}, {cavity_q_mean[1]:.6f}) nm")
print(f"   RMS cavity position fluctuation: {np.sqrt(np.mean(np.sum(cavity_q_centered**2, axis=1))):.6f} nm")

# Calculate autocorrelation for x,y components
print("\n4. Calculating cavity position autocorrelation function...")
print("   Computing <δq(0)·δq(t)> = <δqx(0)·δqx(t)> + <δqy(0)·δqy(t)> (x,y only)")

autocorr = np.zeros(len(cavity_q_centered))
for i in range(2):  # x, y components
    component_autocorr = signal.correlate(cavity_q_centered[:, i], cavity_q_centered[:, i], 
                                          mode='full', method='fft')
    autocorr += component_autocorr[len(component_autocorr)//2:]  # Keep only positive lags

autocorr = autocorr / autocorr[0]  # Normalize

dt_ps = time_cut[1] - time_cut[0]
time_autocorr = np.arange(len(autocorr)) * dt_ps

print(f"   Autocorrelation length: {len(autocorr)} points")
print(f"   Time step: {dt_ps:.4f} ps")

# Compute spectrum from autocorrelation using FFT (Wiener-Khinchin theorem)
print("\n5. Computing cavity position spectrum from autocorrelation (Wiener-Khinchin theorem)...")
print("   S(ω) = FT[<q(0)·q(t)>]")
print("   Using autocorrelation: <qx(0)qx(t)> + <qy(0)qy(t)>")

# Use FFT on the autocorrelation function (Wiener-Khinchin theorem)
dt_seconds = dt_ps * 1e-12  # Convert ps to seconds

# Compute FFT of autocorrelation
spectrum_fft = np.fft.rfft(autocorr)
freq_hz = np.fft.rfftfreq(len(autocorr), dt_seconds)
freq_cm = freq_hz / (3e10)  # Convert Hz to cm^-1

# Power spectrum from FFT
spectrum_power = np.real(spectrum_fft * np.conj(spectrum_fft))

# Multiply by ω² for power spectrum
# I(ω) ∝ ω² × S(ω) where S(ω) is the power spectrum
omega = 2 * np.pi * freq_hz  # Angular frequency in rad/s
spectrum = spectrum_power * (omega**2)

print(f"   Applied ω² factor for power spectrum intensity")
print(f"   Frequency range: 0 - {freq_cm[-1]:.1f} cm^-1")
print(f"   Nyquist frequency: {freq_cm[-1]:.1f} cm^-1")

# Extract cavity position trajectory for a representative window
print("\n6. Extracting cavity position trajectory for representative window...")
# Use last 20% of data or 20 ps window, whichever is smaller
total_time_available = time[-1] - time[0]
window_size = min(20.0, total_time_available * 0.2)
window_start = max(time[0], time[-1] - window_size)
window_end = time[-1]

print(f"   Total time available: {total_time_available:.1f} ps")
print(f"   Using window: {window_start:.1f} - {window_end:.1f} ps")

idx_window = np.where((time >= window_start) & (time <= window_end))[0]
time_window = time[idx_window]
cavity_q_window = cavity_q[idx_window]

print(f"   Window points: {len(time_window)}")
if len(time_window) > 0:
    print(f"   Window range: {time_window[0]:.2f} - {time_window[-1]:.2f} ps")
else:
    print(f"   ⚠ Warning: No data in window, using all available data")
    time_window = time
    cavity_q_window = cavity_q

# Create plot with 3 rows
print("\n7. Creating plot...")

fig, axes = plt.subplots(3, 1, figsize=(14, 14))

# Plot 1: Autocorrelation function
ax1 = axes[0]
max_time_plot = 10  # ps - focus on first 10 ps to see decay
idx_max = int(max_time_plot / dt_ps)
ax1.plot(time_autocorr[:idx_max], autocorr[:idx_max], 'b-', linewidth=2)
ax1.set_xlabel('Time Lag (ps)', fontsize=12)
ax1.set_ylabel('CACF C(t)', fontsize=12)
ax1.set_title('Cavity Position Autocorrelation Function', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, max_time_plot)
ax1.axhline(y=0, color='k', linestyle='--', alpha=0.3)

# Add decay info
decay_idx = np.where(autocorr[:idx_max] < np.exp(-1))[0]
if len(decay_idx) > 0:
    decay_time = time_autocorr[decay_idx[0]]
    ax1.axvline(x=decay_time, color='r', linestyle='--', alpha=0.5, label=f'1/e decay: {decay_time:.2f} ps')
    ax1.legend(fontsize=10)

# Plot 2: Cavity position spectrum from autocorrelation FFT
ax2 = axes[1]
max_freq = 4000  # cm^-1
idx_freq = np.where(freq_cm <= max_freq)[0]
ax2.plot(freq_cm[idx_freq], spectrum[idx_freq], 'r-', linewidth=2)
ax2.set_xlabel('Frequency (cm$^{-1}$)', fontsize=12)
ax2.set_ylabel('Power Spectral Density (arb. units)', fontsize=12)
ax2.set_title('Cavity Position Spectrum (FFT of Autocorrelation with ω² weighting)', fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0, max_freq)
ax2.set_ylim(bottom=0)

# Add vertical line for expected cavity frequency (if known)
# Default to 1560 cm⁻¹ (O-O stretch)
ax2.axvline(x=1560, color='blue', linestyle='--', alpha=0.5, linewidth=1.5, label='Expected cavity (~1560 cm⁻¹)')
ax2.legend(fontsize=9, loc='upper right')

# Add text box with simulation info
info_text = (f'Method: FFT of Autocorrelation\n'
             f'Data points analyzed: {len(time_cut):,}\n'
             f'Autocorrelation: <qx(0)qx(t)> + <qy(0)qy(t)>\n'
             f'Time step: {dt_ps*1000:.2f} fs')
ax2.text(0.98, 0.97, info_text, transform=ax2.transAxes,
         fontsize=10, verticalalignment='top', horizontalalignment='right',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# Plot 3: Cavity position trajectory (x,y components)
ax3 = axes[2]
ax3.plot(time_window, cavity_q_window[:, 0], 'r-', linewidth=1, alpha=0.8, label='$q_x$')
ax3.plot(time_window, cavity_q_window[:, 1], 'g-', linewidth=1, alpha=0.8, label='$q_y$')

# Also plot x,y magnitude (for reference, but note it has harmonics)
cavity_q_mag_window = np.sqrt(cavity_q_window[:, 0]**2 + cavity_q_window[:, 1]**2)
ax3.plot(time_window, cavity_q_mag_window, 'k--', linewidth=1, alpha=0.5, label='|q_xy| (has harmonics)')

ax3.set_xlabel('Time (ps)', fontsize=12)
ax3.set_ylabel('Cavity Position q (nm)', fontsize=12)
ax3.set_title('Cavity Position Trajectory (representative window)', fontsize=14, fontweight='bold')
ax3.grid(True, alpha=0.3)
ax3.legend(fontsize=10, loc='upper right')
ax3.axhline(y=0, color='k', linestyle='-', alpha=0.2, linewidth=0.5)

plt.tight_layout()
output_prefix = "cavity_position_spectrum"
plt.savefig(f'{output_prefix}_full.png', dpi=300, bbox_inches='tight')
print(f"   Saved: {output_prefix}_full.png")

# Also save a high-resolution version of just the spectrum
fig2, ax = plt.subplots(figsize=(14, 6))
ax.plot(freq_cm[idx_freq], spectrum[idx_freq], 'r-', linewidth=2)
ax.set_xlabel('Frequency (cm$^{-1}$)', fontsize=14)
ax.set_ylabel('Power Spectral Density (arb. units)', fontsize=14)
ax.set_title('Cavity Position Spectrum - FFT of Autocorrelation', 
             fontsize=16, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.set_xlim(0, max_freq)
ax.set_ylim(bottom=0)
ax.axvline(x=1560, color='blue', linestyle='--', alpha=0.5, linewidth=1.5, label='Expected cavity (~1560 cm⁻¹)')
ax.legend(fontsize=10)
plt.tight_layout()
plt.savefig(f'{output_prefix}_hires.png', dpi=300, bbox_inches='tight')
print(f"   Saved: {output_prefix}_hires.png")

# Save processed data
print("\n8. Saving processed data...")
np.savez(f'{output_prefix}_data.npz',
         frequencies_cm=freq_cm[idx_freq],
         spectrum=spectrum[idx_freq],
         autocorr_time=time_autocorr[:idx_max],
         autocorr=autocorr[:idx_max],
         cavity_q_window_time=time_window,
         cavity_q_window=cavity_q_window)
print(f"   Saved: {output_prefix}_data.npz")

print("\n" + "=" * 60)
print("Cavity Position Spectrum Calculation Complete!")
print("=" * 60)
print("\nOutput files:")
print(f"  1. {output_prefix}_full.png - 3-panel plot:")
print("     - Top: Autocorrelation function (0-10 ps)")
print("     - Middle: Cavity position spectrum via FFT (0-4000 cm⁻¹)")
print("     - Bottom: Cavity position trajectory (representative window)")
print(f"  2. {output_prefix}_hires.png - High-res spectrum")
print(f"  3. {output_prefix}_data.npz - Processed data")
print("\nMethod: FFT of Autocorrelation (Wiener-Khinchin theorem)")
print("Autocorrelation: <qx(0)qx(t)> + <qy(0)qy(t)>")
print("=" * 60)
