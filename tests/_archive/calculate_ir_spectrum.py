#!/usr/bin/env python3
"""
Calculate IR Spectrum from Dipole Trajectory
"""

import numpy as np
from scipy import signal
from scipy.fft import rfft, rfftfreq
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt

print("=" * 60)
print("IR Spectrum Calculation from Dipole Trajectory")
print("=" * 60)

# 1. Load the data
print("\n1. Loading dipole trajectory data...")
data = np.load('rabi_250mol_40bohr.npz')
time = data['time_ps']
dipole = data['dipole_nm']

print(f"   Total data points: {len(time)}")
print(f"   Time range: {time[0]:.2f} - {time[-1]:.2f} ps")
print(f"   Dipole shape: {dipole.shape}")

# 2. Cut off first 20 ps if available, otherwise use all data
print("\n2. Cutting off equilibration period...")
if time[-1] >= 20.0:
    idx_cutoff = np.where(time >= 20.0)[0][0]
    time_cut = time[idx_cutoff:]
    dipole_cut = dipole[idx_cutoff:]
    print(f"   Cut off first 20 ps")
else:
    time_cut = time
    dipole_cut = dipole
    print(f"   Data < 20 ps ({time[-1]:.1f} ps), using all data")

print(f"   Data points after cutoff: {len(time_cut)}")
print(f"   Time range: {time_cut[0]:.2f} - {time_cut[-1]:.2f} ps")

# 3. Compute dipole components (subtract mean for proper DACF)
print("\n3. Computing dipole with mean subtracted (δM = M - <M>)...")
dipole_mean = np.mean(dipole_cut, axis=0)
dipole_centered = dipole_cut - dipole_mean

print(f"   Mean dipole: ({dipole_mean[0]:.4f}, {dipole_mean[1]:.4f}, {dipole_mean[2]:.4f}) e·nm")
print(f"   RMS dipole fluctuation: {np.sqrt(np.mean(np.sum(dipole_centered**2, axis=1))):.4f} e·nm")

# 4. Calculate autocorrelation for X and Y components only (cavity-coupled directions)
print("\n4. Calculating dipole autocorrelation function (DACF)...")
print("   Computing <δM(0)·δM(t)> = <δMx(0)·δMx(t)> + <δMy(0)·δMy(t)>  [XY only - cavity coupled]")

autocorr = np.zeros(len(dipole_centered))
for i in range(2):  # x, y components only (cavity-coupled)
    component_autocorr = signal.correlate(dipole_centered[:, i], dipole_centered[:, i], 
                                          mode='full', method='fft')
    autocorr += component_autocorr[len(component_autocorr)//2:]  # Keep only positive lags

autocorr = autocorr / autocorr[0]  # Normalize

dt_ps = time_cut[1] - time_cut[0]
time_autocorr = np.arange(len(autocorr)) * dt_ps

print(f"   Autocorrelation length: {len(autocorr)} points")
print(f"   Time step: {dt_ps:.4f} ps")

# 5. Apply window and FFT with ω² weighting for IR absorption
print("\n5. Computing IR spectrum via FFT with ω² weighting...")
print("   IR absorption: α(ω) ∝ ω² × FT[⟨M(0)·M(t)⟩]")
window = np.hanning(len(autocorr))
autocorr_windowed = autocorr * window

spectrum_raw = np.abs(rfft(autocorr_windowed))

# Frequency axis - correct calculation
dt_seconds = dt_ps * 1e-12  # Convert ps to seconds
frequencies_Hz = rfftfreq(len(autocorr_windowed), d=dt_seconds)
frequencies_cm = frequencies_Hz / (3e10)  # Convert Hz to cm^-1 (wavenumber = freq/c)

# Apply ω² weighting for proper IR absorption spectrum
# This suppresses low-frequency modes and enhances vibrational modes
omega = 2 * np.pi * frequencies_Hz  # Angular frequency in rad/s
spectrum = spectrum_raw * omega**2

# Normalize to max = 1 for plotting convenience
spectrum = spectrum / np.max(spectrum[frequencies_cm > 100])  # Normalize ignoring DC

print(f"   Spectrum points: {len(spectrum)}")
print(f"   Time step: {dt_ps:.6f} ps = {dt_ps*1000:.2f} fs")
print(f"   Nyquist frequency: {frequencies_cm[-1]:.1f} cm^-1")
print(f"   Frequency resolution: {frequencies_cm[1]:.4f} cm^-1")
print(f"   Total duration: {len(autocorr)*dt_ps:.1f} ps")

# 6. Plot and save
print("\n6. Creating plots...")

# Create figure with two subplots
fig, axes = plt.subplots(2, 1, figsize=(12, 10))

# Plot 1: Autocorrelation function
ax1 = axes[0]
max_time_plot = 50  # ps
idx_max = int(max_time_plot / dt_ps)
ax1.plot(time_autocorr[:idx_max], autocorr[:idx_max], 'b-', linewidth=1.5)
ax1.set_xlabel('Time Lag (ps)', fontsize=12)
ax1.set_ylabel('Autocorrelation C(t)', fontsize=12)
ax1.set_title('Dipole Moment Autocorrelation Function', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, max_time_plot)
ax1.axhline(y=0, color='k', linestyle='--', alpha=0.3)

# Plot 2: IR spectrum
ax2 = axes[1]
max_freq = 3000  # cm^-1
idx_freq = np.where(frequencies_cm <= max_freq)[0]
ax2.plot(frequencies_cm[idx_freq], spectrum[idx_freq], 'r-', linewidth=1.5)
ax2.set_xlabel('Frequency (cm$^{-1}$)', fontsize=12)
ax2.set_ylabel('Intensity (arb. units)', fontsize=12)
ax2.set_title('IR Absorption Spectrum: α(ω) ∝ ω² × FT[DACF]', 
              fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0, max_freq)

# Add text box with simulation info
info_text = (f'Simulation: 1 ns @ 300 K\n'
             f'Equilibration: 100 ps (coupling OFF)\n'
             f'Production: 900 ps (coupling ON)\n'
             f'Data cutoff: 20 ps\n'
             f'Points analyzed: {len(time_cut):,}')
ax2.text(0.98, 0.97, info_text, transform=ax2.transAxes,
         fontsize=10, verticalalignment='top', horizontalalignment='right',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig('cavity_ir_spectrum.png', dpi=300, bbox_inches='tight')
print("   ✓ Saved: cavity_ir_spectrum.png")

# Also save a high-resolution version of just the spectrum
fig2, ax = plt.subplots(figsize=(14, 6))
ax.plot(frequencies_cm[idx_freq], spectrum[idx_freq], 'r-', linewidth=2)
ax.set_xlabel('Frequency (cm$^{-1}$)', fontsize=14)
ax.set_ylabel('Intensity (arb. units)', fontsize=14)
ax.set_title('IR Spectrum - Cavity-Coupled Diamer System (300 K)', 
             fontsize=16, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.set_xlim(0, max_freq)
plt.tight_layout()
plt.savefig('cavity_ir_spectrum_hires.png', dpi=300, bbox_inches='tight')
print("   ✓ Saved: cavity_ir_spectrum_hires.png")

# 7. Save processed data
print("\n7. Saving processed data...")
np.savez('ir_spectrum_data.npz',
         frequencies_cm=frequencies_cm[idx_freq],
         spectrum=spectrum[idx_freq],
         autocorr_time=time_autocorr[:idx_max],
         autocorr=autocorr[:idx_max])
print("   ✓ Saved: ir_spectrum_data.npz")

print("\n" + "=" * 60)
print("IR Spectrum Calculation Complete!")
print("=" * 60)
print("\nOutput files:")
print("  1. cavity_ir_spectrum.png - Combined plot (autocorr + spectrum)")
print("  2. cavity_ir_spectrum_hires.png - High-res spectrum only")
print("  3. ir_spectrum_data.npz - Processed data for further analysis")
print("\nTo view the images:")
print("  xdg-open cavity_ir_spectrum.png")
print("=" * 60)
