#!/usr/bin/env python3
"""
Calculate IR Spectrum from Dipole Trajectory using Maximum Entropy Method
Reference: https://maximum-entropy-spectrum.readthedocs.io/en/latest/
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
print("IR Spectrum Calculation using Maximum Entropy Method")
print("=" * 60)

# Get filename from command line
if len(sys.argv) < 2:
    print("\nUsage: python analyze_spectrum.py <dipole_file.npz>")
    print("\nSearching for cavity_diamer_lambda*.npz files...")
    import glob
    npz_files = glob.glob("cavity_diamer_lambda*.npz")
    if npz_files:
        # Sort by modification time, get most recent
        npz_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
        input_file = npz_files[0]
        print(f"  Found {len(npz_files)} file(s), using most recent: {input_file}")
    else:
        print("  No cavity_diamer_lambda*.npz files found!")
        sys.exit(1)
else:
    input_file = sys.argv[1]

# Load the data
print(f"\n1. Loading dipole trajectory data from: {input_file}")
data = np.load(input_file, allow_pickle=True)
time = data['time_ps']
dipole = data['dipole_nm']

# Extract lambda value for output filenames
try:
    metadata = data['metadata'].item()
    lambda_coupling = metadata.get('lambda_coupling', 0.001)
except Exception:
    # Fallback if metadata is not accessible (e.g., file still being written)
    lambda_coupling = 0.001
    print("   ⚠ Warning: Could not read metadata, using default lambda=0.001")
output_prefix = f"cavity_ir_spectrum_lambda{lambda_coupling:.4f}"

print(f"   Total data points: {len(time)}")
print(f"   Time range: {time[0]:.2f} - {time[-1]:.2f} ps")
print(f"   Dipole shape: {dipole.shape}")

# Cut off equilibration (first 20 ps or 10% of data, whichever is smaller)
print("\n2. Removing equilibration period...")
total_time = time[-1] - time[0]
cutoff_time = min(20.0, total_time * 0.1)  # 20 ps or 10% of total time
print(f"   Total simulation time: {total_time:.2f} ps")
print(f"   Cutoff time: {cutoff_time:.2f} ps")

if total_time > cutoff_time:
    idx_cutoff = np.where(time >= cutoff_time)[0][0]
    time_cut = time[idx_cutoff:]
    dipole_cut = dipole[idx_cutoff:]
    print(f"   Data points after cutoff: {len(time_cut)}")
    print(f"   Time range: {time_cut[0]:.2f} - {time_cut[-1]:.2f} ps")
else:
    print(f"   ⚠ Warning: Simulation too short ({total_time:.1f} ps), using all data")
    time_cut = time
    dipole_cut = dipole
    print(f"   Data points: {len(time_cut)}")
    print(f"   Time range: {time_cut[0]:.2f} - {time_cut[-1]:.2f} ps")

if len(time_cut) < 100:
    print(f"\n   ERROR: Not enough data points ({len(time_cut)}) for spectral analysis!")
    print(f"   Need at least 100 points. Please run simulation longer.")
    sys.exit(1)

# Compute dipole with mean subtracted (δM = M - <M>)
print("\n3. Computing dipole with mean subtracted (δM = M - <M>)...")
dipole_mean = np.mean(dipole_cut, axis=0)
dipole_centered = dipole_cut - dipole_mean

print(f"   Mean dipole: ({dipole_mean[0]:.4f}, {dipole_mean[1]:.4f}, {dipole_mean[2]:.4f}) e·nm")
print(f"   RMS dipole fluctuation: {np.sqrt(np.mean(np.sum(dipole_centered**2, axis=1))):.4f} e·nm")

# Calculate autocorrelation for each component
print("\n4. Calculating dipole autocorrelation function (DACF)...")
print("   Computing <δM(0)·δM(t)> = <δMx(0)·δMx(t)> + <δMy(0)·δMy(t)> (x,y only)")

autocorr = np.zeros(len(dipole_centered))
for i in range(2):  # x, y components
    component_autocorr = signal.correlate(dipole_centered[:, i], dipole_centered[:, i], 
                                          mode='full', method='fft')
    autocorr += component_autocorr[len(component_autocorr)//2:]  # Keep only positive lags

autocorr = autocorr / autocorr[0]  # Normalize

dt_ps = time_cut[1] - time_cut[0]
time_autocorr = np.arange(len(autocorr)) * dt_ps

print(f"   Autocorrelation length: {len(autocorr)} points")
print(f"   Time step: {dt_ps:.4f} ps")

# Compute spectrum using Maximum Entropy Method (MESA)
print("\n5. Computing IR spectrum using Maximum Entropy Method (MESA)...")
print("   Reference: https://maximum-entropy-spectrum.readthedocs.io/en/latest/")
print("   Using FULL 1 fs resolution - NO SUBSAMPLING!")

# Use x,y dipole magnitude for MESA (mean already subtracted)
dipole_magnitude = np.sqrt(dipole_centered[:, 0]**2 + dipole_centered[:, 1]**2)

print(f"   Time series points for MESA: {len(dipole_magnitude)}")
print(f"   Time step: {dt_ps:.6f} ps = {dt_ps*1000:.2f} fs")

# Use MESA on the FULL dipole time series
M = memspectrum.MESA()

# For large datasets, MESA can be slow. Use a reasonable order.
# But DON'T subsample the data!
print("   Computing MESA (this may take a moment for 980k points)...")
M.solve(dipole_magnitude, optimisation_method='CAT', early_stop=True)

print(f"   MESA order selected: {M.P}")

# Create frequency grid
dt_seconds = dt_ps * 1e-12  # Convert ps to seconds - FULL 1 fs resolution
# Frequency points
n_freq = 20000
freq_max_hz = 1.0 / (2 * dt_seconds)  # Nyquist with 1 fs
freq_hz = np.linspace(0, freq_max_hz, n_freq)
freq_cm = freq_hz / (3e10)  # Convert Hz to cm^-1

# Evaluate spectrum
spectrum_mem = M.spectrum(dt_seconds, freq_hz)

# Multiply by ω² for IR intensity
# I(ω) ∝ ω² × S(ω) where S(ω) is the power spectrum
omega = 2 * np.pi * freq_hz  # Angular frequency in rad/s
spectrum_mem = spectrum_mem * (omega**2)

print(f"   Applied ω² factor for IR absorption intensity")
print(f"   Frequency range: 0 - {freq_cm[-1]:.1f} cm^-1")
print(f"   Nyquist frequency: {freq_cm[-1]:.1f} cm^-1 (full 1 fs resolution!)")

# Extract dipole trajectory for a representative window
print("\n6. Extracting dipole trajectory for representative window...")
# Use last 20% of data or 20 ps window, whichever is smaller
total_time_available = time[-1] - time[0]
window_size = min(20.0, total_time_available * 0.2)
window_start = max(time[0], time[-1] - window_size)
window_end = time[-1]

print(f"   Total time available: {total_time_available:.1f} ps")
print(f"   Using window: {window_start:.1f} - {window_end:.1f} ps")

idx_window = np.where((time >= window_start) & (time <= window_end))[0]
time_window = time[idx_window]
dipole_window = dipole[idx_window]

print(f"   Window points: {len(time_window)}")
if len(time_window) > 0:
    print(f"   Window range: {time_window[0]:.2f} - {time_window[-1]:.2f} ps")
else:
    print(f"   ⚠ Warning: No data in window, using all available data")
    time_window = time
    dipole_window = dipole

# Create plot with 3 rows
print("\n7. Creating plot...")

fig, axes = plt.subplots(3, 1, figsize=(14, 14))

# Plot 1: Autocorrelation function
ax1 = axes[0]
max_time_plot = 10  # ps - focus on first 10 ps to see decay
idx_max = int(max_time_plot / dt_ps)
ax1.plot(time_autocorr[:idx_max], autocorr[:idx_max], 'b-', linewidth=2)
ax1.set_xlabel('Time Lag (ps)', fontsize=12)
ax1.set_ylabel('DACF C(t)', fontsize=12)
ax1.set_title('Dipole Moment Autocorrelation Function', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, max_time_plot)
ax1.axhline(y=0, color='k', linestyle='--', alpha=0.3)

# Add decay info
decay_idx = np.where(autocorr[:idx_max] < np.exp(-1))[0]
if len(decay_idx) > 0:
    decay_time = time_autocorr[decay_idx[0]]
    ax1.axvline(x=decay_time, color='r', linestyle='--', alpha=0.5, label=f'1/e decay: {decay_time:.2f} ps')
    ax1.legend(fontsize=10)

# Plot 2: IR spectrum using MESA
ax2 = axes[1]
max_freq = 4000  # cm^-1 (increased to see N-N peak at ~3290 cm⁻¹)
idx_freq = np.where(freq_cm <= max_freq)[0]
ax2.plot(freq_cm[idx_freq], spectrum_mem[idx_freq], 'r-', linewidth=2)
ax2.set_xlabel('Frequency (cm$^{-1}$)', fontsize=12)
ax2.set_ylabel('IR Absorption Intensity (arb. units)', fontsize=12)
ax2.set_title('IR Spectrum (Maximum Entropy Method with ω² weighting)', fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0, max_freq)
ax2.set_ylim(bottom=0)

# Add vertical lines for expected O-O and N-N peaks
ax2.axvline(x=1555, color='blue', linestyle='--', alpha=0.5, linewidth=1.5, label='Expected O-O (~1560 cm⁻¹)')
ax2.axvline(x=2330, color='green', linestyle='--', alpha=0.5, linewidth=1.5, label='Expected N-N (~2233 cm⁻¹)')
ax2.legend(fontsize=9, loc='upper right')

# Add text box with simulation info
info_text = (f'Method: MESA (Maximum Entropy)\n'
             f'Simulation: 1 ns @ 300 K\n'
             f'Equilibration: 100 ps (coupling OFF)\n'
             f'Production: 900 ps (coupling ON)\n'
             f'Data cutoff: 20 ps\n'
             f'Points analyzed: {len(time_cut):,}\n'
             f'MESA order: {M.P}')
ax2.text(0.98, 0.97, info_text, transform=ax2.transAxes,
         fontsize=10, verticalalignment='top', horizontalalignment='right',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# Plot 3: Dipole moment trajectory (800-820 ps, x/y only)
ax3 = axes[2]
ax3.plot(time_window, dipole_window[:, 0], 'r-', linewidth=1, alpha=0.8, label='$M_x$')
ax3.plot(time_window, dipole_window[:, 1], 'g-', linewidth=1, alpha=0.8, label='$M_y$')

# Also plot x,y magnitude
dipole_mag_window = np.sqrt(dipole_window[:, 0]**2 + dipole_window[:, 1]**2)
ax3.plot(time_window, dipole_mag_window, 'k-', linewidth=1.5, label='|M_xy|')

ax3.set_xlabel('Time (ps)', fontsize=12)
ax3.set_ylabel('Dipole Moment (e·nm)', fontsize=12)
ax3.set_title('Dipole Moment Trajectory (800-820 ps window)', fontsize=14, fontweight='bold')
ax3.grid(True, alpha=0.3)
ax3.set_xlim(800, 820)
ax3.legend(fontsize=10, loc='upper right')
ax3.axhline(y=0, color='k', linestyle='-', alpha=0.2, linewidth=0.5)

plt.tight_layout()
plt.savefig(f'{output_prefix}_full.png', dpi=300, bbox_inches='tight')
print(f"   Saved: {output_prefix}_full.png")

# Also save a high-resolution version of just the MESA spectrum
fig2, ax = plt.subplots(figsize=(6, 6))
ax.plot(freq_cm[idx_freq], spectrum_mem[idx_freq], 'r-', linewidth=2)
# Add vertical lines for expected O-O and N-N peaks
ax.axvline(x=1555, color='blue', linestyle='--', alpha=0.5, linewidth=1.5, label='A-A stretch (~1555 cm⁻¹)')
ax.axvline(x=2330, color='green', linestyle='--', alpha=0.5, linewidth=1.5, label='B-B stretch (~2330 cm⁻¹)')
ax.legend(fontsize=9, loc='upper right')
ax.set_xlim([0, 3000])
ax.set_ylim(bottom=0)
ax.set_xlabel('Frequency (cm$^{-1}$)', fontsize=14)
ax.set_ylabel('Power Spectral Density (arb. units)', fontsize=14)
ax.set_title('IR Spectrum with ExternalLaser Driving at $\omega_d$ = 1555 cm⁻¹\nand Cavity Frequency $\omega_c$ = 1555 cm⁻¹', 
             fontsize=16, fontweight='bold')
plt.savefig(f'{output_prefix}_hires.png', dpi=300, bbox_inches='tight')
print(f"   Saved: {output_prefix}_hires.png")