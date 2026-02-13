#!/usr/bin/env python3
"""
Plot Fictive Temperatures from Energy Components

This script creates a multi-panel plot of fictive temperatures derived from different
energy components using U(T) relationships from temperature sweep data.

Each energy component gets its own panel for clear comparison.

Usage:
    python plot_fictive_temperatures.py <energy_tracker_file>
    
Example:
    python plot_fictive_temperatures.py cavity_coupling_1eneg03_switch_1.0ps/prod-0_energy_tracker.txt
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
from pathlib import Path
from scipy.interpolate import interp1d
import warnings

def load_energy_data(filename):
    """Load energy tracker data from file using numpy."""
    try:
        # First, read all lines to find the header
        with open(filename, 'r') as f:
            lines = f.readlines()
            
        # Find the header line with column names (contains 'time(ps)' but is not a comment)
        header_line = None
        data_start_line = 0
        
        for i, line in enumerate(lines):
            line = line.strip()
            if not line.startswith('#') and 'time(ps)' in line:
                header_line = line
                data_start_line = i + 1
                break
        
        if header_line is None:
            raise ValueError("Could not find column headers in file")
        
        column_names = header_line.split()
        
        # Create a dictionary mapping column names to indices
        columns = {name: i for i, name in enumerate(column_names)}
        
        # Read the data starting from the line after the header
        data = np.loadtxt(filename, skiprows=data_start_line)
        
        return data, columns
        
    except Exception as e:
        print(f"Error loading data from {filename}: {e}")
        sys.exit(1)

def load_temperature_energy_components(filename="potential_energy_components_vs_temperature.txt"):
    """Load the potential energy components vs temperature data."""
    try:
        # Read the header to get column names
        with open(filename, 'r') as f:
            lines = f.readlines()
            
        # Find the header line (first non-comment line)
        header_line = None
        data_start_line = 0
        
        for i, line in enumerate(lines):
            line = line.strip()
            if not line.startswith('#') and line:
                header_line = line
                data_start_line = i + 1
                break
        
        if header_line is None:
            raise ValueError("Could not find column headers in temperature components file")
        
        column_names = header_line.split()
        temp_columns = {name: i for i, name in enumerate(column_names)}
        
        # Read the data
        temp_data = np.loadtxt(filename, skiprows=data_start_line)
        
        return temp_data, temp_columns
        
    except Exception as e:
        print(f"Error: Could not load temperature components file {filename}: {e}")
        print("This file is required for fictive temperature calculation.")
        sys.exit(1)

def create_energy_temperature_interpolators(temp_data, temp_columns):
    """Create interpolating functions U(T) for energy components."""
    interpolators = {}
    temperatures = temp_data[:, temp_columns['temperature']]
    
    # Create interpolators for different energy components
    energy_components = {
        'total_PE': ('total_PE_hartree', 'Total Potential Energy'),
        'harmonic': ('harmonic_hartree', 'Harmonic Bond Energy'), 
        'lj': ('lj_hartree', 'Lennard-Jones Energy'),
        'coulombic': ('coulombic_hartree', 'Coulombic Energy')
    }
    
    for comp_name, (col_name, display_name) in energy_components.items():
        if col_name in temp_columns:
            energies = temp_data[:, temp_columns[col_name]]
            
            # Create interpolator (T -> U) and inverse interpolator (U -> T)
            # Sort by temperature for proper interpolation
            sort_idx = np.argsort(temperatures)
            T_sorted = temperatures[sort_idx]
            U_sorted = energies[sort_idx]
            
            # Remove any duplicates in temperature
            unique_mask = np.diff(T_sorted, prepend=T_sorted[0]-1) != 0
            T_unique = T_sorted[unique_mask]
            U_unique = U_sorted[unique_mask]
            
            if len(T_unique) > 3:  # Need sufficient points for interpolation
                # T -> U interpolator with extrapolation
                U_of_T = interp1d(T_unique, U_unique, kind='cubic', 
                                bounds_error=False, fill_value='extrapolate')
                
                # U -> T interpolator (inverse function)
                # Sort by energy for the inverse
                sort_idx_inv = np.argsort(U_unique)
                U_sorted_inv = U_unique[sort_idx_inv]
                T_sorted_inv = T_unique[sort_idx_inv]
                
                # Remove duplicates in energy
                unique_mask_inv = np.diff(U_sorted_inv, prepend=U_sorted_inv[0]-1) != 0
                U_unique_inv = U_sorted_inv[unique_mask_inv]
                T_unique_inv = T_sorted_inv[unique_mask_inv]
                
                if len(U_unique_inv) > 3:
                    # Create a robust interpolator with linear extrapolation fallback
                    T_of_U = interp1d(U_unique_inv, T_unique_inv, kind='cubic',
                                    bounds_error=False, fill_value='extrapolate')
                    
                    # Create linear extrapolation functions for extreme cases
                    if len(U_unique_inv) >= 3:
                        low_coeffs = np.polyfit(U_unique_inv[:3], T_unique_inv[:3], 1)
                        low_linear = np.poly1d(low_coeffs)
                        
                        high_coeffs = np.polyfit(U_unique_inv[-3:], T_unique_inv[-3:], 1)
                        high_linear = np.poly1d(high_coeffs)
                    else:
                        low_linear = None
                        high_linear = None
                    
                    interpolators[comp_name] = {
                        'U_of_T': U_of_T,
                        'T_of_U': T_of_U,
                        'T_range': (T_unique.min(), T_unique.max()),
                        'U_range': (U_unique.min(), U_unique.max()),
                        'low_linear': low_linear,
                        'high_linear': high_linear,
                        'display_name': display_name
                    }
                    
                    print(f"  Created interpolators for {comp_name} ({display_name})")
                    print(f"    T range: {T_unique.min():.1f} - {T_unique.max():.1f} K")
                    print(f"    U range: {U_unique.min():.6f} - {U_unique.max():.6f} Hartree")
    
    return interpolators

def calculate_energy_components(data, columns):
    """Calculate the energy components we want to analyze."""
    
    # Total molecular potential energy
    molecular_potential = (data[:, columns['harmonic_energy']] + 
                          data[:, columns['lj_energy']] + 
                          data[:, columns['ewald_short_energy']] + 
                          data[:, columns['ewald_long_energy']])
    
    # Individual components for fictive temperature calculation
    harmonic = data[:, columns['harmonic_energy']]
    lj = data[:, columns['lj_energy']]
    coulombic = (data[:, columns['ewald_short_energy']] + 
                data[:, columns['ewald_long_energy']])
    
    return {
        'total_PE': molecular_potential,
        'harmonic': harmonic,
        'lj': lj,
        'coulombic': coulombic
    }

def calculate_fictive_temperatures(energy_dict, interpolators):
    """Calculate fictive temperatures from energy components with improved extrapolation."""
    fictive_temps = {}
    
    for comp_name, interp_data in interpolators.items():
        if comp_name in energy_dict:
            energies = energy_dict[comp_name]
            T_of_U = interp_data['T_of_U']
            U_range = interp_data['U_range']
            low_linear = interp_data.get('low_linear')
            high_linear = interp_data.get('high_linear')
            
            # Calculate fictive temperatures with improved extrapolation
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fictive_T = T_of_U(energies)
            
            # Use linear extrapolation for extreme cases
            if low_linear is not None and high_linear is not None:
                # Identify extreme extrapolation regions (beyond 50% extension of original range)
                U_span = U_range[1] - U_range[0]
                low_extreme = energies < (U_range[0] - 0.5 * U_span)
                high_extreme = energies > (U_range[1] + 0.5 * U_span)
                
                # Apply linear extrapolation for extreme cases
                if np.any(low_extreme):
                    fictive_T[low_extreme] = low_linear(energies[low_extreme])
                
                if np.any(high_extreme):
                    fictive_T[high_extreme] = high_linear(energies[high_extreme])
            
            # Only mask physically unreasonable temperatures (negative or extremely high)
            reasonable_mask = (fictive_T > 0) & (fictive_T < 10000)  # 0-10000 K seems reasonable
            fictive_T[~reasonable_mask] = np.nan
            
            # Count extrapolated vs interpolated points
            interpolated_mask = (energies >= U_range[0]) & (energies <= U_range[1])
            n_interpolated = np.sum(interpolated_mask & reasonable_mask)
            n_extrapolated = np.sum(~interpolated_mask & reasonable_mask)
            n_total_valid = np.sum(reasonable_mask)
            
            fictive_temps[comp_name] = {
                'temperature': fictive_T,
                'interpolated_mask': interpolated_mask,
                'valid_mask': reasonable_mask,
                'n_interpolated': n_interpolated,
                'n_extrapolated': n_extrapolated,
                'n_total_valid': n_total_valid
            }
            
            print(f"  {comp_name}: {n_total_valid}/{len(fictive_T)} valid fictive temperatures")
            print(f"    Interpolated: {n_interpolated}, Extrapolated: {n_extrapolated}")
            
            if n_extrapolated > 0:
                pct_extrapolated = 100 * n_extrapolated / n_total_valid
                print(f"    {pct_extrapolated:.1f}% of valid points are extrapolated")
    
    return fictive_temps

def plot_fictive_temperatures_multi_panel(data, columns, interpolators, fictive_temps, 
                                        actual_temperature, output_file=None):
    """Create a multi-panel plot with each energy component in its own panel."""
    
    time_ps = data[:, columns['time(ps)']]
    
    # Determine grid layout based on number of components
    n_components = len(fictive_temps)
    if n_components == 1:
        rows, cols = 1, 1
    elif n_components == 2:
        rows, cols = 1, 2
    elif n_components <= 4:
        rows, cols = 2, 2
    else:
        rows = int(np.ceil(n_components / 3))
        cols = 3
    
    # Create figure with subplots
    fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 4*rows), sharex=True)
    if n_components == 1:
        axes = [axes]
    elif rows == 1 or cols == 1:
        axes = axes.flatten()
    else:
        axes = axes.flatten()
    
    # Colors for different components
    colors = {'total_PE': 'purple', 'harmonic': 'brown', 'lj': 'pink', 'coulombic': 'cyan'}
    
    # Plot each component in its own panel
    for i, (comp_name, temp_data) in enumerate(fictive_temps.items()):
        if i >= len(axes):
            break
            
        ax = axes[i]
        fictive_T = temp_data['temperature']
        interpolated_mask = temp_data['interpolated_mask']
        valid_mask = temp_data['valid_mask']
        
        # Get display name and color
        display_name = interpolators[comp_name]['display_name']
        color = colors.get(comp_name, 'blue')
        
        # Plot interpolated points as a continuous line
        interp_and_valid = interpolated_mask & valid_mask
        if np.any(interp_and_valid):
            ax.plot(time_ps[interp_and_valid], fictive_T[interp_and_valid], 
                   color=color, linewidth=2, alpha=0.8, label='Interpolated')
        
        # Plot extrapolated points as scattered points
        extrap_and_valid = (~interpolated_mask) & valid_mask
        if np.any(extrap_and_valid):
            ax.scatter(time_ps[extrap_and_valid], fictive_T[extrap_and_valid],
                      color=color, alpha=0.4, s=0.5, label='Extrapolated')
        
        # Add actual temperature reference line
        temp_mean = np.mean(actual_temperature)
        ax.axhline(y=temp_mean, color='red', linestyle='--', alpha=0.7, 
                  label=f'Actual Mean T: {temp_mean:.1f} K')
        
        # Formatting
        ax.set_ylabel('Fictive Temperature (K)', fontsize=11)
        ax.set_title(f'{display_name}\nFictive Temperature', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9)
        
        # Add statistics text box
        n_total_valid = temp_data['n_total_valid']
        n_interpolated = temp_data['n_interpolated']
        n_extrapolated = temp_data['n_extrapolated']
        
        if n_total_valid > 0:
            pct_extrapolated = 100 * n_extrapolated / n_total_valid
            stats_text = f'{n_total_valid}/{len(fictive_T)} valid\n{pct_extrapolated:.1f}% extrapolated'
        else:
            stats_text = 'No valid data'
            
        props = dict(boxstyle='round', facecolor='lightgray', alpha=0.8)
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=8,
               verticalalignment='top', bbox=props)
        
        # Set reasonable y-axis limits
        if np.any(valid_mask):
            T_valid = fictive_T[valid_mask]
            T_min, T_max = np.percentile(T_valid, [1, 99])  # Use 1st and 99th percentiles
            T_range = T_max - T_min
            ax.set_ylim(max(0, T_min - 0.1*T_range), T_max + 0.1*T_range)
    
    # Remove unused subplots
    for i in range(n_components, len(axes)):
        fig.delaxes(axes[i])
    
    # Set x-axis label for bottom row
    for i in range(max(0, len(axes) - cols), len(axes)):
        if i < n_components:
            axes[i].set_xlabel('Time (ps)', fontsize=11)
    
    # Overall title
    plt.suptitle('Fictive Temperatures from Energy Components\n(Derived using U(T) relationships)', 
                 fontsize=14, fontweight='bold')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save or show plot
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Fictive temperature plot saved to {output_file}")
    else:
        plt.show()
    
    return fig

def print_fictive_temperature_summary(fictive_temps, actual_temperature):
    """Print summary statistics for fictive temperatures."""
    
    print("\n" + "="*60)
    print("FICTIVE TEMPERATURE SUMMARY")
    print("="*60)
    
    temp_mean = np.mean(actual_temperature)
    temp_std = np.std(actual_temperature)
    print(f"Actual Temperature: {temp_mean:.2f} ± {temp_std:.2f} K")
    
    for comp_name, temp_data in fictive_temps.items():
        fictive_T = temp_data['temperature']
        valid_mask = temp_data['valid_mask']
        
        if np.any(valid_mask):
            T_valid = fictive_T[valid_mask]
            fictive_mean = np.mean(T_valid)
            fictive_std = np.std(T_valid)
            fictive_min = np.min(T_valid)
            fictive_max = np.max(T_valid)
            
            print(f"\n{comp_name.upper()} Fictive Temperature:")
            print(f"  Mean: {fictive_mean:.2f} ± {fictive_std:.2f} K")
            print(f"  Range: {fictive_min:.2f} - {fictive_max:.2f} K")
            print(f"  Difference from actual: {fictive_mean - temp_mean:.2f} K")
            print(f"  Valid data points: {temp_data['n_total_valid']}/{len(fictive_T)}")
        else:
            print(f"\n{comp_name.upper()}: No valid fictive temperature data")
    
    print("="*60)

def main():
    parser = argparse.ArgumentParser(
        description='Create multi-panel fictive temperature plots from energy components',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument('input_file', type=str, 
                       help='Energy tracker output file (e.g., prod-0_energy_tracker.txt)')
    parser.add_argument('--output', '-o', type=str, 
                       help='Output plot file (default: show plot)')
    parser.add_argument('--summary', action='store_true', 
                       help='Print fictive temperature summary statistics')
    parser.add_argument('--temp_components', type=str, 
                       default='potential_energy_components_vs_temperature.txt',
                       help='Temperature components file for fictive temperature calculation')
    
    args = parser.parse_args()
    
    # Check if input file exists
    input_path = Path(args.input_file)
    if not input_path.exists():
        print(f"Error: Input file '{args.input_file}' not found")
        sys.exit(1)
    
    print(f"Loading energy data from: {args.input_file}")
    
    # Load energy tracker data
    data, columns = load_energy_data(args.input_file)
    
    print(f"Loaded {len(data)} data points")
    print(f"Time range: {data[0, columns['time(ps)']]:.6f} - {data[-1, columns['time(ps)']]:.6f} ps")
    
    # Load temperature components data
    print(f"\nLoading temperature components from: {args.temp_components}")
    temp_data, temp_columns = load_temperature_energy_components(args.temp_components)
    
    # Create interpolators for fictive temperature calculation
    print("Creating energy-temperature interpolators...")
    interpolators = create_energy_temperature_interpolators(temp_data, temp_columns)
    
    if not interpolators:
        print("Error: No interpolators could be created")
        sys.exit(1)
    
    print(f"Successfully created {len(interpolators)} interpolators")
    
    # Calculate energy components
    print("\nCalculating energy components...")
    energy_dict = calculate_energy_components(data, columns)
    
    # Calculate fictive temperatures
    print("Calculating fictive temperatures...")
    fictive_temps = calculate_fictive_temperatures(energy_dict, interpolators)
    
    # Get actual temperature
    actual_temperature = data[:, columns['temperature']]
    
    # Print summary if requested
    if args.summary:
        print_fictive_temperature_summary(fictive_temps, actual_temperature)
    
    # Create multi-panel plot
    print("\nCreating multi-panel fictive temperature plot...")
    fig = plot_fictive_temperatures_multi_panel(data, columns, interpolators, 
                                              fictive_temps, actual_temperature, 
                                              args.output)
    
    if not args.output:
        print("Close the plot window to exit.")

if __name__ == '__main__':
    main()

