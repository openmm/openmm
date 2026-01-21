#!/usr/bin/env python3
"""
IR Spectrum Analysis for Water Cavity MD
=========================================

Computes IR absorption spectra from dipole moment trajectories and identifies:
- Polariton splitting (upper/lower branches)
- Rabi splitting Ω_R vs coupling strength λ
- Vibrational modes (OH stretch, HOH bend)

Based on Equations 169-170 from the paper:
  C_μμ(t) = ⟨μ(t)·μ(0)⟩
  n(ω)α(ω) = (β ω²)/(4ε₀Vc) × FT[C_μμ(t)]
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.signal import find_peaks
from scipy.fft import rfft, rfftfreq

# Physical constants
KB_KJMOL = 0.008314463  # kJ/(mol·K)
ELEM_CHARGE = 1.60217663e-19  # C
EPSILON_0 = 8.8541878128e-12  # F/m
C_LIGHT = 299792458  # m/s
HARTREE_TO_CM = 219474.63  # cm⁻¹/Hartree


def load_trajectory(npz_file):
    """
    Load trajectory from NPZ file.
    
    Returns
    -------
    time_ps : ndarray
        Time points in ps
    dipole_nm : ndarray (N, 3)
        Dipole moments in e·nm
    cavity_nm : ndarray (N, 3) or None
        Cavity positions in nm (None for baseline/cavity-free simulations)
    metadata : dict
        Simulation metadata
    """
    data = np.load(npz_file, allow_pickle=True)
    time_ps = data['time_ps']
    dipole_nm = data['dipole_nm']
    cavity_nm = data.get('cavity_nm', None)  # May be None for baseline simulations
    metadata = data['metadata'].item() if 'metadata' in data else {}
    
    print(f"Loaded: {npz_file}")
    print(f"  Time range: {time_ps[0]:.1f} - {time_ps[-1]:.1f} ps")
    print(f"  Data points: {len(time_ps)}")
    print(f"  Lambda: {metadata.get('lambda_coupling', 'N/A')}")
    print(f"  Flexible water: {metadata.get('flexible_water', False)}")
    print(f"  Status: {metadata.get('status', 'unknown')}")
    if cavity_nm is None:
        print(f"  Type: Baseline (cavity-free) simulation")
    
    return time_ps, dipole_nm, cavity_nm, metadata


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
    # Process each component
    dipole = dipole_nm.copy()
    
    if subtract_mean:
        dipole = dipole - dipole.mean(axis=0)
    
    N = len(dipole)
    C_mumu = np.zeros(N)
    
    # Compute autocorrelation using FFT (fast method)
    for alpha in range(2):  # x, y components only (cavity-coupled)
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


def compute_ir_spectrum(time_ps, C_mumu, temperature_K, box_size_nm):
    """
    Compute IR absorption spectrum from autocorrelation function.
    
    From paper Eq. 169-170:
      n(ω)α(ω) = (β ω²)/(4ε₀Vc) × ∫ dt e^(-iωt) C_μμ(t)
    
    Parameters
    ----------
    time_ps : ndarray
        Time points in ps
    C_mumu : ndarray
        Dipole autocorrelation function
    temperature_K : float
        Temperature in K
    box_size_nm : float
        Box size in nm
        
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
    
    # Compute spectrum with proper normalization
    # n(ω)α(ω) ∝ β ω² × FT[C_μμ(t)]
    beta = 1.0 / (KB_KJMOL * temperature_K)  # mol/kJ
    volume_m3 = (box_size_nm * 1e-9)**3  # m³
    
    # Prefactor (convert units appropriately)
    # This is approximate - main goal is to see peak positions
    omega_rad = 2 * np.pi * freq_hz
    spectrum = np.abs(fft_result) * omega_rad**2
    
    # Normalize to max = 1 for plotting
    if spectrum.max() > 0:
        spectrum /= spectrum.max()
    
    return freq_cm, spectrum


def find_water_peaks(freq_cm, spectrum):
    """
    Identify expected water vibrational peaks for baseline validation.
    
    Parameters
    ----------
    freq_cm : ndarray
        Frequencies in cm⁻¹
    spectrum : ndarray
        IR spectrum
        
    Returns
    -------
    peaks : dict
        Dictionary with peak positions for water modes
    """
    from scipy.signal import find_peaks
    
    peaks = {}
    
    # OH stretch region (3000-4000 cm⁻¹)
    mask_oh = (freq_cm > 3000) & (freq_cm < 4000)
    if mask_oh.sum() > 0:
        freq_oh = freq_cm[mask_oh]
        spec_oh = spectrum[mask_oh]
        peaks_idx, properties = find_peaks(spec_oh, height=0.1, prominence=0.05)
        if len(peaks_idx) > 0:
            # Find highest peak in OH region
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
        peaks_idx, properties = find_peaks(spec_bend, height=0.05, prominence=0.03)
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
        peaks_idx, properties = find_peaks(spec_lib, height=0.05, prominence=0.03)
        if len(peaks_idx) > 0:
            peaks['librational'] = freq_lib[peaks_idx]
        else:
            peaks['librational'] = None
    
    return peaks


def find_polariton_peaks(freq_cm, spectrum, cavity_freq_cm, search_window=500):
    """
    Identify polariton peaks (upper and lower branches).
    
    Parameters
    ----------
    freq_cm : ndarray
        Frequencies in cm⁻¹
    spectrum : ndarray
        IR spectrum
    cavity_freq_cm : float
        Cavity frequency in cm⁻¹
    search_window : float
        Search window around cavity frequency in cm⁻¹
        
    Returns
    -------
    peaks : dict
        Dictionary with peak positions and Rabi splitting
    """
    # Find region around cavity frequency
    mask = (freq_cm > cavity_freq_cm - search_window) & \
           (freq_cm < cavity_freq_cm + search_window)
    
    freq_region = freq_cm[mask]
    spec_region = spectrum[mask]
    
    # Find peaks
    peaks_idx, properties = find_peaks(spec_region, height=0.1, prominence=0.05)
    
    if len(peaks_idx) == 0:
        return {'lower': None, 'upper': None, 'rabi_splitting': 0}
    
    peak_freqs = freq_region[peaks_idx]
    peak_heights = spec_region[peaks_idx]
    
    # Sort by frequency
    sorted_idx = np.argsort(peak_freqs)
    peak_freqs = peak_freqs[sorted_idx]
    peak_heights = peak_heights[sorted_idx]
    
    # Identify upper and lower polaritons
    cavity_idx = np.argmin(np.abs(peak_freqs - cavity_freq_cm))
    
    if len(peak_freqs) >= 2:
        # Find peaks on either side of cavity frequency
        lower_peaks = peak_freqs[peak_freqs < cavity_freq_cm]
        upper_peaks = peak_freqs[peak_freqs > cavity_freq_cm]
        
        lower_pol = lower_peaks[-1] if len(lower_peaks) > 0 else None
        upper_pol = upper_peaks[0] if len(upper_peaks) > 0 else None
        
        if lower_pol is not None and upper_pol is not None:
            rabi_splitting = upper_pol - lower_pol
        else:
            rabi_splitting = 0
    else:
        lower_pol = None
        upper_pol = None
        rabi_splitting = 0
    
    return {
        'lower': lower_pol,
        'upper': upper_pol,
        'rabi_splitting': rabi_splitting,
        'all_peaks': peak_freqs,
        'all_heights': peak_heights
    }


def plot_spectrum(freq_cm, spectrum, metadata, polariton_info=None, 
                 output_file=None, xlim=(0, 5000)):
    """
    Plot IR spectrum with polariton peak identification.
    
    Parameters
    ----------
    freq_cm : ndarray
        Frequencies in cm⁻¹
    spectrum : ndarray
        IR spectrum
    metadata : dict
        Simulation metadata
    polariton_info : dict, optional
        Polariton peak information
    output_file : str, optional
        Output filename for saving plot
    xlim : tuple
        X-axis limits (min_freq, max_freq) in cm⁻¹
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot spectrum
    ax.plot(freq_cm, spectrum, 'b-', linewidth=1.5, label='IR Spectrum')
    
    # Mark cavity frequency
    cavity_freq = metadata.get('cavity_freq_cm', 3656)
    ax.axvline(cavity_freq, color='gray', linestyle='--', alpha=0.5,
               label=f'Cavity: {cavity_freq} cm⁻¹')
    
    # Mark polariton peaks if available
    if polariton_info:
        if polariton_info['lower']:
            ax.axvline(polariton_info['lower'], color='red', linestyle=':',
                      label=f'Lower: {polariton_info["lower"]:.0f} cm⁻¹')
        if polariton_info['upper']:
            ax.axvline(polariton_info['upper'], color='orange', linestyle=':',
                      label=f'Upper: {polariton_info["upper"]:.0f} cm⁻¹')
        
        # Show Rabi splitting
        if polariton_info['rabi_splitting'] > 0:
            omega_R = polariton_info['rabi_splitting']
            ax.text(0.98, 0.98, f'Ω_R = {omega_R:.0f} cm⁻¹',
                   transform=ax.transAxes, ha='right', va='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Formatting
    ax.set_xlabel('Frequency (cm⁻¹)', fontsize=12)
    ax.set_ylabel('Normalized Intensity', fontsize=12)
    ax.set_xlim(xlim)
    ax.set_ylim(0, 1.1)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper left')
    
    # Title with metadata
    lambda_val = metadata.get('lambda_coupling', 'N/A')
    temp = metadata.get('temperature_K', 'N/A')
    title = f'IR Spectrum: λ={lambda_val}, T={temp}K'
    ax.set_title(title, fontsize=14)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"  Saved plot: {output_file}")
    else:
        plt.show()
    
    plt.close()


def analyze_trajectory(npz_file, plot_output=None, xlim=(2500, 4500)):
    """
    Complete analysis pipeline for a single trajectory.
    
    Parameters
    ----------
    npz_file : str
        Path to NPZ file
    plot_output : str, optional
        Output filename for plot
    xlim : tuple
        Frequency range for plotting
        
    Returns
    -------
    results : dict
        Analysis results including spectrum and polariton info
    """
    print(f"\n{'='*80}")
    print(f"Analyzing: {npz_file}")
    print(f"{'='*80}\n")
    
    # Load trajectory
    time_ps, dipole_nm, cavity_nm, metadata = load_trajectory(npz_file)
    
    # Compute autocorrelation
    print(f"\nComputing dipole autocorrelation...")
    C_mumu = compute_dipole_autocorrelation(dipole_nm, subtract_mean=True)
    print(f"  C_μμ(0) = {C_mumu[0]:.2e} (e·nm)²")
    print(f"  C_μμ(t_max) = {C_mumu[-1]:.2e} (e·nm)²")
    
    # Compute spectrum
    print(f"\nComputing IR spectrum...")
    temperature_K = metadata.get('temperature_K', 300.0)
    box_size_nm = metadata.get('box_size_nm', 3.0)
    freq_cm, spectrum = compute_ir_spectrum(time_ps, C_mumu, temperature_K, box_size_nm)
    print(f"  Frequency range: {freq_cm[0]:.1f} - {freq_cm[-1]:.1f} cm⁻¹")
    print(f"  Spectral resolution: {freq_cm[1]-freq_cm[0]:.2f} cm⁻¹")
    
    # Find polariton peaks (if cavity simulation)
    print(f"\nIdentifying spectral peaks...")
    cavity_freq_cm = metadata.get('cavity_freq_cm', 3656)
    lambda_coupling = metadata.get('lambda_coupling', 0.0)
    
    if lambda_coupling > 0 and cavity_nm is not None:
        # Cavity simulation - look for polaritons
        polariton_info = find_polariton_peaks(freq_cm, spectrum, cavity_freq_cm)
        
        if polariton_info['lower'] and polariton_info['upper']:
            print(f"  Lower polariton: {polariton_info['lower']:.1f} cm⁻¹")
            print(f"  Upper polariton: {polariton_info['upper']:.1f} cm⁻¹")
            print(f"  Rabi splitting Ω_R: {polariton_info['rabi_splitting']:.1f} cm⁻¹")
        else:
            print(f"  Could not identify distinct polariton peaks")
            if 'all_peaks' in polariton_info and len(polariton_info['all_peaks']) > 0:
                print(f"  Found peaks at: {polariton_info['all_peaks']}")
    else:
        # Baseline simulation - validate water peaks
        polariton_info = None
        water_peaks = find_water_peaks(freq_cm, spectrum)
        
        print(f"\n--- Water IR Peak Validation ---")
        print(f"Expected peaks for liquid water:")
        print(f"  OH stretch (H-bonded): 3200-3400 cm⁻¹")
        print(f"  Free OH stretch:       ~3700 cm⁻¹")
        print(f"  HOH bending mode:      ~1645 cm⁻¹")
        print(f"  Librational modes:     400-900 cm⁻¹")
        
        print(f"\nFound peaks:")
        if water_peaks.get('oh_stretch'):
            oh_freq = water_peaks['oh_stretch']
            print(f"  OH stretch:     {oh_freq:.0f} cm⁻¹ (height={water_peaks['oh_height']:.2f})")
            # Validate
            if 3000 < oh_freq < 4000:
                print(f"    ✓ PASS: OH stretch in expected range (3000-4000 cm⁻¹)")
            else:
                print(f"    ✗ FAIL: OH stretch outside expected range!")
        else:
            print(f"  OH stretch:     NOT FOUND")
            print(f"    ✗ FAIL: No OH stretch peak detected!")
            print(f"    → Check if water bonds are truly flexible (constraints=None)")
        
        if water_peaks.get('hoh_bend'):
            bend_freq = water_peaks['hoh_bend']
            print(f"  HOH bend:       {bend_freq:.0f} cm⁻¹ (height={water_peaks['bend_height']:.2f})")
            if 1500 < bend_freq < 1800:
                print(f"    ✓ PASS: HOH bend in expected range (1500-1800 cm⁻¹)")
            else:
                print(f"    ✗ FAIL: HOH bend outside expected range!")
        else:
            print(f"  HOH bend:       NOT FOUND")
        
        if water_peaks.get('librational') is not None and len(water_peaks['librational']) > 0:
            print(f"  Librational:    {len(water_peaks['librational'])} peak(s) at {water_peaks['librational'][:3]} cm⁻¹")
        else:
            print(f"  Librational:    NOT FOUND (this is optional)")
        
        # Overall validation
        print(f"\n--- Validation Summary ---")
        has_oh = water_peaks.get('oh_stretch') is not None
        oh_valid = has_oh and (3000 < water_peaks['oh_stretch'] < 4000)
        has_bend = water_peaks.get('hoh_bend') is not None
        
        if oh_valid:
            print(f"✓ BASELINE VALIDATION PASSED")
            print(f"  Flexible water bonds are working correctly")
            print(f"  IR spectrum shows expected OH stretch peak")
            print(f"  Ready to proceed with cavity simulations")
        else:
            print(f"✗ BASELINE VALIDATION FAILED")
            if not has_oh:
                print(f"  No OH stretch peak found → Bonds may still be constrained")
                print(f"  Check: rigidWater=False and constraints=None in simulation")
            else:
                print(f"  OH stretch peak at wrong frequency")
                print(f"  Force field parameters may need adjustment")
            print(f"  DO NOT proceed with cavity simulations until this passes")
    
    # Plot
    if plot_output is None:
        plot_output = str(Path(npz_file).with_suffix('.png'))
    
    print(f"\nGenerating plot...")
    plot_spectrum(freq_cm, spectrum, metadata, polariton_info, plot_output, xlim)
    
    results = {
        'freq_cm': freq_cm,
        'spectrum': spectrum,
        'C_mumu': C_mumu,
        'polariton_info': polariton_info,
        'metadata': metadata
    }
    
    print(f"\nAnalysis complete!")
    return results


def compare_coupling_strengths(npz_files, output_file='rabi_splitting_vs_lambda.png'):
    """
    Compare spectra at different coupling strengths and plot Rabi splitting vs λ.
    
    Parameters
    ----------
    npz_files : list of str
        List of NPZ files with different λ values
    output_file : str
        Output filename for comparison plot
    """
    print(f"\n{'='*80}")
    print(f"Comparing coupling strengths")
    print(f"{'='*80}\n")
    
    results = []
    for npz_file in npz_files:
        result = analyze_trajectory(npz_file, plot_output=None, xlim=(2500, 4500))
        results.append(result)
    
    # Extract lambda and Rabi splitting
    lambdas = []
    rabi_splittings = []
    
    for r in results:
        lambda_val = r['metadata'].get('lambda_coupling', 0)
        omega_R = r['polariton_info']['rabi_splitting']
        lambdas.append(lambda_val)
        rabi_splittings.append(omega_R)
    
    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Panel 1: Overlaid spectra
    for r in results:
        lambda_val = r['metadata'].get('lambda_coupling', 0)
        freq = r['freq_cm']
        spec = r['spectrum']
        mask = (freq > 2500) & (freq < 4500)
        ax1.plot(freq[mask], spec[mask], label=f'λ={lambda_val:.3f}', linewidth=1.5)
    
    ax1.axvline(3656, color='gray', linestyle='--', alpha=0.5, label='Cavity')
    ax1.set_xlabel('Frequency (cm⁻¹)', fontsize=12)
    ax1.set_ylabel('Normalized Intensity', fontsize=12)
    ax1.set_title('IR Spectra at Different Coupling Strengths', fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Rabi splitting vs lambda
    if len(lambdas) > 1:
        ax2.plot(lambdas, rabi_splittings, 'bo-', markersize=8, linewidth=2)
        
        # Linear fit
        if len(lambdas) >= 2:
            coeffs = np.polyfit(lambdas, rabi_splittings, 1)
            fit_lambdas = np.linspace(min(lambdas), max(lambdas), 100)
            fit_omega = np.polyval(coeffs, fit_lambdas)
            ax2.plot(fit_lambdas, fit_omega, 'r--', label=f'Linear fit: Ω_R = {coeffs[0]:.0f}λ + {coeffs[1]:.0f}')
        
        ax2.set_xlabel('Coupling strength λ', fontsize=12)
        ax2.set_ylabel('Rabi splitting Ω_R (cm⁻¹)', fontsize=12)
        ax2.set_title('Rabi Splitting vs Coupling Strength', fontsize=14)
        ax2.legend()
        ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nSaved comparison plot: {output_file}")
    plt.close()


def main():
    """Main entry point with command-line argument parsing."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze Water Cavity MD IR Spectra')
    parser.add_argument('files', nargs='+', help='NPZ file(s) to analyze')
    parser.add_argument('--compare', action='store_true',
                       help='Compare multiple files and plot Rabi splitting vs λ')
    parser.add_argument('--xlim', nargs=2, type=float, default=[2500, 4500],
                       help='X-axis limits for plotting (default: 2500 4500)')
    parser.add_argument('--output', type=str, default=None,
                       help='Output filename for plot')
    
    args = parser.parse_args()
    
    if args.compare and len(args.files) > 1:
        compare_coupling_strengths(args.files, args.output or 'rabi_comparison.png')
    else:
        for npz_file in args.files:
            analyze_trajectory(npz_file, args.output, tuple(args.xlim))


if __name__ == "__main__":
    main()
