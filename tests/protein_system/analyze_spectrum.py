#!/usr/bin/env python3
"""
Calculate IR Spectrum from Protein Dipole Trajectory using Maximum Entropy Method
Reference: https://maximum-entropy-spectrum.readthedocs.io/en/latest/
"""

import numpy as np
from scipy import signal
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import sys
import os
import json

# Use LaTeX-style fonts (Computer Modern)
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Computer Modern Roman', 'DejaVu Serif', 'Times New Roman'],
    'mathtext.fontset': 'cm',
    'axes.formatter.use_mathtext': True,
})

print("=" * 60)
print("Protein IR Spectrum Calculation using FFT")
print("=" * 60)

def _filter_valid_data(time, dipole):
    """Filter out zero/uninitialized entries from pre-allocated arrays.
    
    When streaming to disk, arrays are pre-allocated and may contain zeros
    for entries that were never written (e.g., if simulation stopped early).
    
    Parameters
    ----------
    time : np.ndarray
        Time array in ps.
    dipole : np.ndarray
        Dipole array of shape (N, 3) in e·nm.
        
    Returns
    -------
    time_valid : np.ndarray
        Filtered time array with only valid entries.
    dipole_valid : np.ndarray
        Filtered dipole array with only valid entries.
    n_filtered : int
        Number of invalid entries that were removed.
    """
    # Identify rows where dipole is all zeros (uninitialized)
    # A real dipole from a charged system should never be exactly [0, 0, 0]
    valid_mask = np.any(dipole != 0, axis=1)
    
    # Also check for time consistency (time should be monotonically increasing)
    if np.sum(valid_mask) > 1:
        valid_indices = np.where(valid_mask)[0]
        time_valid_temp = time[valid_mask]
        dt = np.diff(time_valid_temp)
        
        # If there are large negative jumps in time, data is corrupted
        if np.any(dt < -0.1):  # Allow small numerical noise
            print("   ⚠ Warning: Detected non-monotonic time, filtering by contiguous block")
            # Find the first discontinuity and truncate there
            bad_idx = np.where(dt < -0.1)[0]
            if len(bad_idx) > 0:
                cutoff = valid_indices[bad_idx[0]]
                valid_mask[cutoff:] = False
    
    time_valid = time[valid_mask]
    dipole_valid = dipole[valid_mask]
    n_filtered = len(time) - len(time_valid)
    
    return time_valid, dipole_valid, n_filtered


def _load_npz(input_file):
    data = np.load(input_file, allow_pickle=True)
    time = data["time_ps"]
    dipole = data["dipole_nm"]
    metadata = {}
    try:
        metadata = data["metadata"].item()
    except Exception:
        metadata = {}
    
    # Filter out invalid (zero) entries
    time, dipole, n_filtered = _filter_valid_data(time, dipole)
    if n_filtered > 0:
        print(f"   ⚠ Filtered {n_filtered} invalid (zero) entries from data")
    
    return time, dipole, metadata


def _load_stream_npy(input_file):
    stream = np.load(input_file, allow_pickle=False, mmap_mode="r")
    time = np.array(stream["time_ps"])  # Copy to allow filtering
    dipole = np.array(stream["dipole_nm"])
    metadata = {}
    meta_path = input_file.replace("_stream.npy", "_metadata.json")
    if os.path.exists(meta_path):
        with open(meta_path, "r") as handle:
            metadata = json.load(handle)
    
    # Filter out invalid (zero) entries - critical for stream files!
    time, dipole, n_filtered = _filter_valid_data(time, dipole)
    if n_filtered > 0:
        print(f"   ⚠ Filtered {n_filtered} invalid (zero) entries from stream")
        print(f"   (This is normal if simulation was interrupted)")
    
    return time, dipole, metadata


def _discover_input():
    import glob
    npz_files = glob.glob("protein_cavity_lambda*.npz")
    stream_files = glob.glob("protein_cavity_lambda*_stream.npy")
    
    # Prefer npz files over stream files (npz is typically more complete)
    # Only use stream if npz doesn't exist
    if npz_files:
        candidates = npz_files
        print("  Found .npz file(s), preferring over stream")
    elif stream_files:
        candidates = stream_files
        print("  Using stream file(s) (no .npz found)")
    else:
        return None
    
    candidates.sort(key=lambda x: os.path.getmtime(x), reverse=True)
    return candidates[0]


# Get filename from command line
if len(sys.argv) < 2:
    print("\nUsage: python analyze_spectrum.py <dipole_file.npz|stream.npy>")
    print("\nSearching for protein_cavity_lambda*.npz or *_stream.npy files...")
    input_file = _discover_input()
    if input_file:
        print(f"  Using most recent: {input_file}")
    else:
        print("  No protein_cavity_lambda*.npz or *_stream.npy files found!")
        sys.exit(1)
else:
    input_file = sys.argv[1]

# Load the data
print(f"\n1. Loading dipole trajectory data from: {input_file}")
if input_file.endswith("_stream.npy"):
    time, dipole, metadata = _load_stream_npy(input_file)
else:
    time, dipole, metadata = _load_npz(input_file)

lambda_coupling = metadata.get("lambda_coupling", 0.001)
cavity_freq_cm = metadata.get("cavity_freq_cm", 3663.0)
temperature_K = metadata.get("temperature_K", 300.0)
dipole_dt_ps_meta = metadata.get("dipole_dt_ps", None)
dipole_interval = metadata.get("dipole_interval", 1)
if not metadata:
    print("   ⚠ Warning: Metadata missing, using defaults")
if dipole_dt_ps_meta is not None:
    print(f"   Dipole sampling interval: {dipole_dt_ps_meta:.4f} ps (every {dipole_interval} MD steps)")

output_prefix = f"protein_ir_spectrum_lambda{lambda_coupling:.4f}"

print(f"   Total data points: {len(time)}")
print(f"   Time range: {time[0]:.2f} - {time[-1]:.2f} ps")
print(f"   Dipole shape: {dipole.shape}")

# Cut off equilibration (first 20 ps or 10% of data)
print("\n2. Removing equilibration period...")
total_time = time[-1] - time[0]
cutoff_time = min(20.0, total_time * 0.1)
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
    print("   Need at least 100 points. Please run simulation longer.")
    sys.exit(1)

# Compute dipole with mean subtracted
print("\n3. Computing dipole with mean subtracted (δM = M - <M>)...")
dipole_mean = np.mean(dipole_cut, axis=0)
dipole_centered = dipole_cut - dipole_mean

print(f"   Mean dipole: ({dipole_mean[0]:.4f}, {dipole_mean[1]:.4f}, {dipole_mean[2]:.4f}) e·nm")
print(f"   RMS dipole fluctuation: {np.sqrt(np.mean(np.sum(dipole_centered**2, axis=1))):.4f} e·nm")

# Autocorrelation (x,y only)
print("\n4. Calculating dipole autocorrelation function (DACF)...")
print("   Computing <δM(0)·δM(t)> using x,y components")

autocorr = np.zeros(len(dipole_centered))
for i in range(2):
    component_autocorr = signal.correlate(
        dipole_centered[:, i], dipole_centered[:, i], mode="full", method="fft"
    )
    autocorr += component_autocorr[len(component_autocorr) // 2:]

autocorr = autocorr / autocorr[0]

dt_ps = time_cut[1] - time_cut[0]
time_autocorr = np.arange(len(autocorr)) * dt_ps

print(f"   Autocorrelation length: {len(autocorr)} points")
print(f"   Time step: {dt_ps:.4f} ps")

# Compute IR spectrum using Welch's method on centered dipole
print("\n5. Computing IR spectrum using Welch's method...")
print("   Welch's method provides better spectral resolution than FFT of DACF")

dt_seconds = dt_ps * 1e-12
fs = 1.0 / dt_ps  # Sampling frequency in ps^-1

# Use Welch's method on each dipole component
# nperseg controls frequency resolution vs noise tradeoff
# Larger nperseg = better frequency resolution but more noise
nperseg = min(8192, len(dipole_centered) // 4)
print(f"   Welch nperseg: {nperseg} (freq resolution: {fs/nperseg:.1f} ps^-1)")

# Compute power spectrum for x and y components and sum
freqs_welch, psd_x = signal.welch(dipole_centered[:, 0], fs=fs, nperseg=nperseg, 
                                   window='hann', scaling='density')
_, psd_y = signal.welch(dipole_centered[:, 1], fs=fs, nperseg=nperseg,
                         window='hann', scaling='density')
psd_total = psd_x + psd_y

# Convert frequency from ps^-1 to cm^-1
# freq [ps^-1] * 1e12 [1/s per 1/ps] / 3e10 [cm/s] = freq [cm^-1]
freq_cm = freqs_welch * 1e12 / 3e10

# Apply ω² weighting for IR absorption intensity
# I(ω) ∝ ω² × S(ω) where S(ω) is the power spectrum
omega = 2 * np.pi * freqs_welch * 1e12  # Convert to rad/s
spectrum_welch = psd_total * (omega**2)

# Also compute FFT-based spectrum for comparison (normalized autocorrelation method)
n = len(autocorr)
freq_hz_fft = np.fft.rfftfreq(n, d=dt_seconds)
freq_cm_fft = freq_hz_fft / 3e10
autocorr_fft = np.fft.rfft(autocorr)
power_spectrum_fft = np.real(autocorr_fft)
omega_fft = 2 * np.pi * freq_hz_fft
spectrum_fft = power_spectrum_fft * (omega_fft**2)

print("   Applied ω² factor for IR absorption intensity")
print(f"   Frequency range: 0 - {freq_cm[-1]:.1f} cm^-1")
print(f"   Nyquist frequency: {freq_cm[-1]:.1f} cm^-1")

# Use full production time range for dipole plot
print("\n6. Using full dipole trajectory for fluctuation plot...")
time_window = time
dipole_window = dipole
print(f"   Window points: {len(time_window)}")
print(f"   Time range: {time_window[0]:.2f} - {time_window[-1]:.2f} ps")

# Plots
print("\n7. Creating plots...")
fig, axes = plt.subplots(3, 1, figsize=(6, 6))

ax1 = axes[0]
max_time_plot = 10
idx_max = int(max_time_plot / dt_ps)
ax1.plot(time_autocorr[:idx_max], autocorr[:idx_max], "b-", linewidth=1)
ax1.set_xlabel("Time Lag (ps)", fontsize=8)
ax1.set_ylabel("DACF C(t)", fontsize=8)
ax1.set_title("Dipole Moment Autocorrelation Function", fontsize=9, fontweight="bold")
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, max_time_plot)
ax1.axhline(y=0, color="k", linestyle="--", alpha=0.3)
ax1.tick_params(axis='both', labelsize=7)

decay_idx = np.where(autocorr[:idx_max] < np.exp(-1))[0]
if len(decay_idx) > 0:
    decay_time = time_autocorr[decay_idx[0]]
    ax1.axvline(
        x=decay_time, color="r", linestyle="--", alpha=0.5,
        label=f"1/e decay: {decay_time:.2f} ps"
    )
    ax1.legend(fontsize=7)

ax2 = axes[1]
max_freq = 6000
idx_freq = np.where(freq_cm <= max_freq)[0]
ax2.plot(freq_cm[idx_freq], spectrum_welch[idx_freq], "r-", linewidth=1)
ax2.set_xlabel("Frequency (cm$^{-1}$)", fontsize=8)
ax2.set_ylabel("IR Intensity (arb. units)", fontsize=8)
ax2.set_title("IR Spectrum (Welch method with ω² weighting)", fontsize=9, fontweight="bold")
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0, max_freq)
ax2.set_ylim(bottom=0)
ax2.axvline(x=cavity_freq_cm, color="blue", linestyle="--", alpha=0.5, linewidth=1,
            label=f"Cavity ({cavity_freq_cm:.0f} cm$^{{-1}}$)")
ax2.legend(fontsize=7, loc="upper right")
ax2.tick_params(axis='both', labelsize=7)

# info_text = (
#     f"Method: Welch PSD\n"
#     f"Temperature: {temperature_K:.0f} K\n"
#     f"Cavity: {cavity_freq_cm:.0f} cm$^{{-1}}$\n"
#     f"Data cutoff: {cutoff_time:.1f} ps\n"
#     f"Points analyzed: {len(time_cut):,}\n"
#     f"Welch nperseg: {nperseg:,}"
# )
# ax2.text(
#     0.98, 0.97, info_text, transform=ax2.transAxes,
#     fontsize=10, verticalalignment="top", horizontalalignment="right",
#     bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
# )

ax3 = axes[2]
# Subtract mean from window to show fluctuations more clearly
dipole_window_mean = np.mean(dipole_window, axis=0)
dipole_window_centered = dipole_window - dipole_window_mean
ax3.plot(time_window, dipole_window_centered[:, 0], "r-", linewidth=0.5, alpha=0.8, label="$δM_x$")
ax3.plot(time_window, dipole_window_centered[:, 1], "g-", linewidth=0.5, alpha=0.8, label="$δM_y$")
dipole_mag_centered = np.sqrt(dipole_window_centered[:, 0]**2 + dipole_window_centered[:, 1]**2)
ax3.plot(time_window, dipole_mag_centered, "k-", linewidth=0.5, alpha=0.6, label="|δM_xy|")
ax3.set_xlabel("Time (ps)", fontsize=8)
ax3.set_ylabel("Dipole Fluctuation (e·nm)", fontsize=8)
ax3.set_title("Dipole Moment Fluctuations (mean subtracted)", fontsize=9, fontweight="bold")
ax3.grid(True, alpha=0.3)
ax3.set_xlim(time_window[0], time_window[-1])
ax3.legend(fontsize=7, loc="upper right")
ax3.axhline(y=0, color="k", linestyle="-", alpha=0.3, linewidth=0.5)
ax3.tick_params(axis='both', labelsize=7)

plt.tight_layout()
plt.savefig(f"{output_prefix}_full.png", dpi=300, bbox_inches="tight")
print(f"   Saved: {output_prefix}_full.png")

fig2, ax = plt.subplots(figsize=(6, 3))
ax.plot(freq_cm[idx_freq], spectrum_welch[idx_freq], "r-", linewidth=1)
ax.set_xlabel("Frequency (cm$^{-1}$)", fontsize=9)
ax.set_ylabel("IR Intensity (arb. units)", fontsize=9)
ax.set_title("IR Spectrum - Welch Method with ω² weighting", fontsize=10, fontweight="bold")
ax.grid(True, alpha=0.3)
ax.set_xlim(0, max_freq)
ax.set_ylim(bottom=0)
ax.axvline(x=cavity_freq_cm, color="blue", linestyle="--", alpha=0.5, linewidth=1,
           label=f"Cavity ({cavity_freq_cm:.0f} cm$^{{-1}}$)")
ax.legend(fontsize=8)
ax.tick_params(axis='both', labelsize=8)
plt.tight_layout()
plt.savefig(f"{output_prefix}_hires.png", dpi=300, bbox_inches="tight")
print(f"   Saved: {output_prefix}_hires.png")

# Save processed data
print("\n8. Saving processed data...")
np.savez(
    f"{output_prefix}_data.npz",
    frequencies_cm=freq_cm[idx_freq],
    spectrum_welch=spectrum_welch[idx_freq],
    spectrum_fft=spectrum_fft[np.where(freq_cm_fft <= max_freq)[0]],
    frequencies_cm_fft=freq_cm_fft[np.where(freq_cm_fft <= max_freq)[0]],
    autocorr_time=time_autocorr[:idx_max],
    autocorr=autocorr[:idx_max],
    dipole_window_time=time_window,
    dipole_window=dipole_window,
    welch_nperseg=nperseg,
)
print(f"   Saved: {output_prefix}_data.npz")

print("\n" + "=" * 60)
print("Protein IR Spectrum Calculation Complete!")
print("=" * 60)
print("\nOutput files:")
print(f"  1. {output_prefix}_full.png - 3-panel plot")
print(f"  2. {output_prefix}_hires.png - High-res FFT spectrum")
print(f"  3. {output_prefix}_data.npz - Processed data")
print("=" * 60)
