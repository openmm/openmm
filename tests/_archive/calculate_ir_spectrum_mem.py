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

print("=" * 60)
print("IR Spectrum Calculation using Maximum Entropy Method")
print("=" * 60)

# 1. Load the data
print("\n1. Loading dipole trajectory data...")
data = np.load('cavity_diamer_dipole.npz')
time = data['time_ps']
dipole = data['dipole_nm']

print(f"   Total data points: {len(time)}")
print(f"   Time range: {time[0]:.2f} - {time[-1]:.2f} ps")
print(f"   Dipole shape: {dipole.shape}")

# 2. Cut off first 20 ps as requested
print("\n2. Cutting off first 20 ps...")
idx_cutoff = np.where(time >= 20.0)[0][0]
time_cut = time[idx_cutoff:]
dipole_cut = dipole[idx_cutoff:]

print(f"   Data points after cutoff: {len(time_cut)}")
print(f"   Time range: {time_cut[0]:.2f} - {time_cut[-1]:.2f} ps")

# 3. Compute dipole with mean subtracted (δM = M - <M>)
print("\n3. Computing dipole with mean subtracted (δM = M - <M>)...")
dipole_mean = np.mean(dipole_cut, axis=0)
dipole_centered = dipole_cut - dipole_mean

print(f"   Mean dipole: ({dipole_mean[0]:.4f}, {dipole_mean[1]:.4f}, {dipole_mean[2]:.4f}) e·nm")
print(f"   RMS dipole fluctuation: {np.sqrt(np.mean(np.sum(dipole_centered**2, axis=1))):.4f} e·nm")

# 4. Calculate autocorrelation for each component
print("\n4. Calculating dipole autocorrelation function (DACF)...")
print("   Computing <δM(0)·δM(t)> = <δMx(0)·δMx(t)> + <δMy(0)·δMy(t)> + <δMz(0)·δMz(t)>")

autocorr = np.zeros(len(dipole_centered))
for i in range(3):  # x, y, z components
    component_autocorr = signal.correlate(dipole_centered[:, i], dipole_centered[:, i], 
                                          mode='full', method='fft')
    autocorr += component_autocorr[len(component_autocorr)//2:]  # Keep only positive lags

autocorr = autocorr / autocorr[0]  # Normalize

dt_ps = time_cut[1] - time_cut[0]
time_autocorr = np.arange(len(autocorr)) * dt_ps

print(f"   Autocorrelation length: {len(autocorr)} points")
print(f"   Time step: {dt_ps:.4f} ps")

# 5. Compute spectrum using Maximum Entropy Method (MESA)
print("\n5. Computing IR spectrum using Maximum Entropy Method (MESA)...")
print("   Reference: https://maximum-entropy-spectrum.readthedocs.io/en/latest/")
print("   Using FULL 1 fs resolution - NO SUBSAMPLING!")

# Use total dipole magnitude for MESA (mean already subtracted)
dipole_magnitude = np.sqrt(np.sum(dipole_centered**2, axis=1))

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

print(f"   Frequency range: 0 - {freq_cm[-1]:.1f} cm^-1")
print(f"   Nyquist frequency: {freq_cm[-1]:.1f} cm^-1 (full 1 fs resolution!)")

# 6. Extract dipole trajectory for 800-820 ps window
print("\n6. Extracting dipole trajectory for 800-820 ps window...")
idx_window = np.where((time >= 800.0) & (time <= 820.0))[0]
time_window = time[idx_window]
dipole_window = dipole[idx_window]

print(f"   Window points: {len(time_window)}")
print(f"   Window range: {time_window[0]:.2f} - {time_window[-1]:.2f} ps")

# 7. Create comprehensive plot with 3 rows
print("\n7. Creating comprehensive plot...")

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
max_freq = 3000  # cm^-1
idx_freq = np.where(freq_cm <= max_freq)[0]
ax2.plot(freq_cm[idx_freq], spectrum_mem[idx_freq], 'r-', linewidth=2)
ax2.set_xlabel('Frequency (cm$^{-1}$)', fontsize=12)
ax2.set_ylabel('Power Spectral Density (arb. units)', fontsize=12)
ax2.set_title('IR Spectrum (Maximum Entropy Method)', fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0, max_freq)
ax2.set_ylim(bottom=0)

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

# Plot 3: Dipole moment trajectory (800-820 ps)
ax3 = axes[2]
ax3.plot(time_window, dipole_window[:, 0], 'r-', linewidth=1, alpha=0.8, label='$M_x$')
ax3.plot(time_window, dipole_window[:, 1], 'g-', linewidth=1, alpha=0.8, label='$M_y$')
ax3.plot(time_window, dipole_window[:, 2], 'b-', linewidth=1, alpha=0.8, label='$M_z$')

# Also plot total magnitude
dipole_mag_window = np.linalg.norm(dipole_window, axis=1)
ax3.plot(time_window, dipole_mag_window, 'k-', linewidth=1.5, label='|M|')

ax3.set_xlabel('Time (ps)', fontsize=12)
ax3.set_ylabel('Dipole Moment (e·nm)', fontsize=12)
ax3.set_title('Dipole Moment Trajectory (800-820 ps window)', fontsize=14, fontweight='bold')
ax3.grid(True, alpha=0.3)
ax3.set_xlim(800, 820)
ax3.legend(fontsize=10, loc='upper right')
ax3.axhline(y=0, color='k', linestyle='-', alpha=0.2, linewidth=0.5)

plt.tight_layout()
plt.savefig('cavity_ir_spectrum_mem.png', dpi=300, bbox_inches='tight')
print("   ✓ Saved: cavity_ir_spectrum_mem.png")

# Also save a high-resolution version of just the MESA spectrum
fig2, ax = plt.subplots(figsize=(14, 6))
ax.plot(freq_cm[idx_freq], spectrum_mem[idx_freq], 'r-', linewidth=2)
ax.set_xlabel('Frequency (cm$^{-1}$)', fontsize=14)
ax.set_ylabel('Power Spectral Density (arb. units)', fontsize=14)
ax.set_title('IR Spectrum - Maximum Entropy Method (MESA)', 
             fontsize=16, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.set_xlim(0, max_freq)
ax.set_ylim(bottom=0)
plt.tight_layout()
plt.savefig('cavity_ir_spectrum_mem_hires.png', dpi=300, bbox_inches='tight')
print("   ✓ Saved: cavity_ir_spectrum_mem_hires.png")

# 8. Save processed data
print("\n8. Saving processed data...")
np.savez('ir_spectrum_mem_data.npz',
         frequencies_cm=freq_cm[idx_freq],
         spectrum=spectrum_mem[idx_freq],
         autocorr_time=time_autocorr[:idx_max],
         autocorr=autocorr[:idx_max],
         dipole_window_time=time_window,
         dipole_window=dipole_window,
         mesa_order=M.P)
print("   ✓ Saved: ir_spectrum_mem_data.npz")

print("\n" + "=" * 60)
print("IR Spectrum Calculation Complete!")
print("=" * 60)
print("\nOutput files:")
print("  1. cavity_ir_spectrum_mem.png - 3-panel plot:")
print("     - Top: Autocorrelation function (0-10 ps)")
print("     - Middle: IR spectrum via MESA (0-3000 cm⁻¹)")
print("     - Bottom: Dipole trajectory (800-820 ps)")
print("  2. cavity_ir_spectrum_mem_hires.png - High-res MESA spectrum")
print("  3. ir_spectrum_mem_data.npz - Processed data")
print("\nMethod: Maximum Entropy Spectral Analysis (MESA)")
print("Reference: https://maximum-entropy-spectrum.readthedocs.io/en/latest/")
print("\nTo view:")
print("  xdg-open cavity_ir_spectrum_mem.png")
print("=" * 60)
