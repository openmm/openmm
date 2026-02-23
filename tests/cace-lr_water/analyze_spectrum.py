#!/usr/bin/env python3
"""
IR Spectrum Analysis for UMA Water MD
======================================

Computes IR absorption spectra from dipole moment trajectories and identifies:
- OH stretch vibration (~3200-3700 cm⁻¹)
- HOH bending mode (~1645 cm⁻¹)
- Librational modes (400-900 cm⁻¹)

Based on dipole autocorrelation function:
  C_μμ(t) = ⟨μ(t)·μ(0)⟩
  α(ω) ∝ ω² × FT[C_μμ(t)]
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.signal import find_peaks
from scipy.fft import rfft, rfftfreq

# Physical constants
KB_KJMOL = 0.008314463  # kJ/(mol·K)
C_LIGHT = 299792458  # m/s


def load_trajectory(npz_file):
    """
    Load trajectory from NPZ file.
    
    Returns
    -------
    time_ps : ndarray
        Time points in ps
    dipole_nm : ndarray (N, 3)
        Dipole moments in e·nm
    metadata : dict
        Simulation metadata
    """
    data = np.load(npz_file, allow_pickle=True)
    time_ps = data['time_ps']
    dipole_nm = data['dipole_nm']
    metadata = data['metadata'].item() if 'metadata' in data else {}
    
    print(f"Loaded: {npz_file}")
    print(f"  Time range: {time_ps[0]:.3f} - {time_ps[-1]:.3f} ps")
    print(f"  Data points: {len(time_ps)}")
    print(f"  Model: {metadata.get('model', 'unknown')}")
    print(f"  Lambda: {metadata.get('lambda_coupling', 'N/A')}")
    print(f"  Status: {metadata.get('status', 'unknown')}")
    
    return time_ps, dipole_nm, metadata


def compute_dipole_autocorrelation(dipole_nm, subtract_mean=True):
    """
    Compute dipole autocorrelation function C_μμ(t).
    
    Parameters
    ----------
    dipole_nm : ndarray (N, 3)
        Dipole moment trajectory in e·nm
    subtract_mean : bool
        Whether to subtract mean (removes DC component)
        
    Returns
    -------
    C_mumu : ndarray (N,)
        Autocorrelation function
    """
    dipole = dipole_nm.copy()
    
    if subtract_mean:
        dipole = dipole - dipole.mean(axis=0)
    
    N = len(dipole)
    C_mumu = np.zeros(N)
    
    # Compute autocorrelation using FFT (fast method)
    for alpha in range(3):  # All components
        mu_alpha = dipole[:, alpha]
        
        # FFT-based autocorrelation
        nfft = 2 * N
        mu_fft = np.fft.fft(mu_alpha, n=nfft)
        acf = np.fft.ifft(mu_fft * np.conj(mu_fft))
        C_mumu += acf[:N].real
    
    # Normalize by number of samples
    norm = np.arange(N, 0, -1)
    C_mumu /= norm
    
    return C_mumu


def compute_ir_spectrum(time_ps, C_mumu, temperature_K=300.0):
    """
    Compute IR absorption spectrum from autocorrelation function.
    
    Parameters
    ----------
    time_ps : ndarray
        Time points in ps
    C_mumu : ndarray
        Dipole autocorrelation function
    temperature_K : float
        Temperature in K
        
    Returns
    -------
    freq_cm : ndarray
        Frequencies in cm⁻¹
    spectrum : ndarray
        IR absorption spectrum (arbitrary units)
    """
    # Get timestep
    dt_ps = time_ps[1] - time_ps[0]
    dt_s = dt_ps * 1e-12  # Convert to seconds
    
    # Apply windowing to reduce spectral leakage
    window = np.hanning(len(C_mumu))
    C_windowed = C_mumu * window
    
    # FFT (only positive frequencies)
    N = len(C_windowed)
    fft_result = rfft(C_windowed)
    freq_hz = rfftfreq(N, dt_s)
    
    # Convert to wavenumbers (cm⁻¹)
    freq_cm = freq_hz / (C_LIGHT * 100)  # Hz -> cm⁻¹
    
    # Compute spectrum with ω² prefactor
    omega_rad = 2 * np.pi * freq_hz
    spectrum = np.abs(fft_result) * omega_rad**2
    
    # Normalize to max = 1 for plotting
    if spectrum.max() > 0:
        spectrum /= spectrum.max()
    
    return freq_cm, spectrum


def find_water_peaks(freq_cm, spectrum):
    """
    Identify expected water vibrational peaks.
    
    Returns
    -------
    peaks : dict
        Dictionary with peak positions for water modes
    """
    peaks = {}
    
    # OH stretch region (3000-4000 cm⁻¹)
    mask_oh = (freq_cm > 3000) & (freq_cm < 4000)
    if mask_oh.sum() > 0:
        freq_oh = freq_cm[mask_oh]
        spec_oh = spectrum[mask_oh]
        peaks_idx, properties = find_peaks(spec_oh, height=0.05, prominence=0.03)
        if len(peaks_idx) > 0:
            max_idx = np.argmax(spec_oh[peaks_idx])
            peaks['oh_stretch'] = freq_oh[peaks_idx[max_idx]]
            peaks['oh_height'] = spec_oh[peaks_idx[max_idx]]
        else:
            peaks['oh_stretch'] = None
            peaks['oh_height'] = 0
    
    # HOH bend region (1500-1800 cm⁻¹)
    mask_bend = (freq_cm > 1500) & (freq_cm < 1800)
    if mask_bend.sum() > 0:
        freq_bend = freq_cm[mask_bend]
        spec_bend = spectrum[mask_bend]
        peaks_idx, properties = find_peaks(spec_bend, height=0.03, prominence=0.02)
        if len(peaks_idx) > 0:
            max_idx = np.argmax(spec_bend[peaks_idx])
            peaks['hoh_bend'] = freq_bend[peaks_idx[max_idx]]
            peaks['bend_height'] = spec_bend[peaks_idx[max_idx]]
        else:
            peaks['hoh_bend'] = None
            peaks['bend_height'] = 0
    
    # Librational modes (400-900 cm⁻¹)
    mask_lib = (freq_cm > 400) & (freq_cm < 900)
    if mask_lib.sum() > 0:
        freq_lib = freq_cm[mask_lib]
        spec_lib = spectrum[mask_lib]
        peaks_idx, properties = find_peaks(spec_lib, height=0.03, prominence=0.02)
        if len(peaks_idx) > 0:
            peaks['librational'] = freq_lib[peaks_idx]
        else:
            peaks['librational'] = None
    
    return peaks


def plot_spectrum(freq_cm, spectrum, metadata, water_peaks=None,
                 output_file=None, xlim=(0, 4500)):
    """
    Plot IR spectrum with peak identification.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot spectrum
    ax.plot(freq_cm, spectrum, 'b-', linewidth=1.5, label='IR Spectrum')
    
    # Mark expected water peaks
    ax.axvline(3500, color='green', linestyle='--', alpha=0.3, label='OH stretch (~3500 cm⁻¹)')
    ax.axvline(1645, color='orange', linestyle='--', alpha=0.3, label='HOH bend (~1645 cm⁻¹)')
    
    # Mark detected peaks
    if water_peaks:
        if water_peaks.get('oh_stretch'):
            ax.axvline(water_peaks['oh_stretch'], color='red', linestyle=':',
                      label=f'Found OH: {water_peaks["oh_stretch"]:.0f} cm⁻¹')
        if water_peaks.get('hoh_bend'):
            ax.axvline(water_peaks['hoh_bend'], color='purple', linestyle=':',
                      label=f'Found HOH: {water_peaks["hoh_bend"]:.0f} cm⁻¹')
    
    # Formatting
    ax.set_xlabel('Frequency (cm⁻¹)', fontsize=12)
    ax.set_ylabel('Normalized Intensity', fontsize=12)
    ax.set_xlim(xlim)
    ax.set_ylim(0, 1.1)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=9)
    
    # Title with metadata
    model = metadata.get('model', 'UMA')
    num_mol = metadata.get('num_molecules', 'N/A')
    temp = metadata.get('temperature_K', 'N/A')
    title = f'IR Spectrum: {model}, {num_mol} waters, T={temp}K'
    ax.set_title(title, fontsize=14)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"  Saved plot: {output_file}")
    else:
        plt.show()
    
    plt.close()


def analyze_trajectory(npz_file, plot_output=None, xlim=(0, 4500)):
    """
    Complete analysis pipeline for a single trajectory.
    """
    print(f"\n{'='*80}")
    print(f"Analyzing: {npz_file}")
    print(f"{'='*80}\n")
    
    # Load trajectory
    time_ps, dipole_nm, metadata = load_trajectory(npz_file)
    
    # Check if we have enough data
    if len(time_ps) < 100:
        print(f"\n⚠ Warning: Only {len(time_ps)} data points")
        print(f"  Spectral resolution will be poor")
        print(f"  Recommend at least 10000 points for good resolution")
    
    # Compute autocorrelation
    print(f"\nComputing dipole autocorrelation...")
    C_mumu = compute_dipole_autocorrelation(dipole_nm, subtract_mean=True)
    print(f"  C_μμ(0) = {C_mumu[0]:.2e} (e·nm)²")
    print(f"  C_μμ(t_max) = {C_mumu[-1]:.2e} (e·nm)²")
    
    # Compute spectrum
    print(f"\nComputing IR spectrum...")
    temperature_K = metadata.get('temperature_K', 300.0)
    freq_cm, spectrum = compute_ir_spectrum(time_ps, C_mumu, temperature_K)
    print(f"  Frequency range: {freq_cm[0]:.1f} - {freq_cm[-1]:.1f} cm⁻¹")
    print(f"  Spectral resolution: {freq_cm[1]-freq_cm[0]:.2f} cm⁻¹")
    
    # Find water peaks
    print(f"\nIdentifying spectral peaks...")
    water_peaks = find_water_peaks(freq_cm, spectrum)
    
    print(f"\n--- Water IR Peak Analysis ---")
    print(f"Expected peaks for liquid water:")
    print(f"  OH stretch (H-bonded): 3200-3400 cm⁻¹")
    print(f"  Free OH stretch:       ~3700 cm⁻¹")
    print(f"  HOH bending mode:      ~1645 cm⁻¹")
    print(f"  Librational modes:     400-900 cm⁻¹")
    
    print(f"\nFound peaks:")
    if water_peaks.get('oh_stretch'):
        oh_freq = water_peaks['oh_stretch']
        print(f"  OH stretch:     {oh_freq:.0f} cm⁻¹ (height={water_peaks['oh_height']:.2f})")
        if 3000 < oh_freq < 4000:
            print(f"    PASS: OH stretch in expected range")
        else:
            print(f"    FAIL: OH stretch outside expected range!")
    else:
        print(f"  OH stretch:     NOT FOUND")
        print(f"    → May need longer simulation for better resolution")
    
    if water_peaks.get('hoh_bend'):
        bend_freq = water_peaks['hoh_bend']
        print(f"  HOH bend:       {bend_freq:.0f} cm⁻¹ (height={water_peaks['bend_height']:.2f})")
        if 1500 < bend_freq < 1800:
            print(f"    PASS: HOH bend in expected range")
        else:
            print(f"    FAIL: HOH bend outside expected range!")
    else:
        print(f"  HOH bend:       NOT FOUND")
    
    if water_peaks.get('librational') is not None and len(water_peaks['librational']) > 0:
        print(f"  Librational:    {len(water_peaks['librational'])} peak(s)")
    else:
        print(f"  Librational:    NOT FOUND")
    
    # Plot
    if plot_output is None:
        plot_output = str(Path(npz_file).with_suffix('.png'))
    
    print(f"\nGenerating plot...")
    plot_spectrum(freq_cm, spectrum, metadata, water_peaks, plot_output, xlim)
    
    results = {
        'freq_cm': freq_cm,
        'spectrum': spectrum,
        'C_mumu': C_mumu,
        'water_peaks': water_peaks,
        'metadata': metadata
    }
    
    print(f"\nAnalysis complete!")
    return results


def main():
    """Main entry point with command-line argument parsing."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze UMA Water MD IR Spectra')
    parser.add_argument('files', nargs='+', help='NPZ file(s) to analyze')
    parser.add_argument('--xlim', nargs=2, type=float, default=[0, 4500],
                       help='X-axis limits for plotting (default: 0 4500)')
    parser.add_argument('--output', type=str, default=None,
                       help='Output filename for plot')
    
    args = parser.parse_args()
    
    for npz_file in args.files:
        analyze_trajectory(npz_file, args.output, tuple(args.xlim))


if __name__ == "__main__":
    main()
