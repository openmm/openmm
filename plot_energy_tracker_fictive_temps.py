#!/usr/bin/env python3
"""
Plot Energy Tracker Output with Fictive Temperatures

This script reads the output from EnergyTracker and creates a three-panel plot:
1. Energy components vs time (molecular, cavity, reservoir, universe)
2. Temperature vs time
3. Fictive temperatures vs time (derived from energy components using U(T) relationships)

The fictive temperature calculation uses interpolating functions derived from 
potential_energy_components_vs_temperature.txt to invert energy->temperature.

Usage:
    python plot_energy_tracker_fictive_temps.py <energy_tracker_file>
    
Example:
    python plot_energy_tracker_fictive_temps.py cavity_coupling_1eneg03_switch_1.0ps/prod-0_energy_tracker.txt
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
        print(f"Warning: Could not load temperature components file {filename}: {e}")
        print("Fictive temperature calculation will be skipped.")
        return None, None

def create_energy_temperature_interpolators(temp_data, temp_columns):
    """Create interpolating functions U(T) for energy components."""
    if temp_data is None or temp_columns is None:
        return {}
    
    interpolators = {}
    temperatures = temp_data[:, temp_columns['temperature']]
    
    # Create interpolators for different energy components
    energy_components = {
        'total_PE': 'total_PE_hartree',
        'harmonic': 'harmonic_hartree', 
        'lj': 'lj_hartree',
        'coulombic': 'coulombic_hartree'
    }
    
    for comp_name, col_name in energy_components.items():
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
                # T -> U interpolator with better extrapolation
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
                    # Create a more robust interpolator with linear extrapolation fallback
                    T_of_U = interp1d(U_unique_inv, T_unique_inv, kind='cubic',
                                    bounds_error=False, fill_value='extrapolate')
                    
                    # Create linear extrapolation functions for extreme cases
                    # Linear fit to first 3 points for low energy extrapolation
                    if len(U_unique_inv) >= 3:
                        low_coeffs = np.polyfit(U_unique_inv[:3], T_unique_inv[:3], 1)
                        low_linear = np.poly1d(low_coeffs)
                        
                        # Linear fit to last 3 points for high energy extrapolation  
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
                        'high_linear': high_linear
                    }
                    
                    print(f"  Created interpolators for {comp_name}")
                    print(f"    T range: {T_unique.min():.1f} - {T_unique.max():.1f} K")
                    print(f"    U range: {U_unique.min():.6f} - {U_unique.max():.6f} Hartree")
    
    return interpolators

def calculate_energy_components(data, columns):
    """Calculate the energy components we want to plot."""
    
    # Total molecular energy = molecular kinetic + molecular potential
    # Molecular potential = harmonic + lj + ewald_short + ewald_long
    molecular_potential = (data[:, columns['harmonic_energy']] + 
                          data[:, columns['lj_energy']] + 
                          data[:, columns['ewald_short_energy']] + 
                          data[:, columns['ewald_long_energy']])
    molecular_total = data[:, columns['molecular_kinetic_energy']] + molecular_potential
    
    # Individual components for fictive temperature calculation
    harmonic = data[:, columns['harmonic_energy']]
    lj = data[:, columns['lj_energy']]
    coulombic = (data[:, columns['ewald_short_energy']] + 
                data[:, columns['ewald_long_energy']])
    
    # Total cavity energy = cavity kinetic + cavity potential
    cavity_total = (data[:, columns['cavity_kinetic_energy']] + 
                   data[:, columns['cavity_total_potential_energy']])
    
    # Total reservoir energy (already calculated in file)
    reservoir_total = data[:, columns['total_reservoir_energy']]
    
    # Universe energy (already calculated in file - should be conserved)
    universe_total = data[:, columns['universe_total_energy']]
    
    return {
        'molecular_total': molecular_total,
        'cavity_total': cavity_total, 
        'reservoir_total': reservoir_total,
        'universe_total': universe_total,
        'molecular_potential': molecular_potential,
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
            
            fictive_temps[comp_name] = fictive_T
            
            print(f"  {comp_name}: {n_total_valid}/{len(fictive_T)} valid fictive temperatures")
            print(f"    Interpolated: {n_interpolated}, Extrapolated: {n_extrapolated}")
            
            if n_extrapolated > 0:
                pct_extrapolated = 100 * n_extrapolated / n_total_valid
                print(f"    {pct_extrapolated:.1f}% of valid points are extrapolated")
    
    return fictive_temps

def plot_energy_temperature_and_fictive(data, columns, interpolators, output_file=None):
    """Create a three-panel plot of energy components, temperature, and fictive temperatures."""
    
    # Calculate energy components
    energy_dict = calculate_energy_components(data, columns)
    
    # Calculate fictive temperatures
    fictive_temps = calculate_fictive_temperatures(energy_dict, interpolators)
    
    # Create figure with three subplots
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 14), sharex=True)
    
    # Panel 1: Energy components
    time_ps = data[:, columns['time(ps)']]
    
    ax1.plot(time_ps, energy_dict['molecular_total'], label='Total Molecular Energy', 
             linewidth=2, color='blue')
    ax1.plot(time_ps, energy_dict['cavity_total'], label='Total Cavity Energy', 
             linewidth=2, color='red')
    ax1.plot(time_ps, energy_dict['reservoir_total'], label='Total Reservoir Energy', 
             linewidth=2, color='green')
    ax1.plot(time_ps, energy_dict['universe_total'], label='Universe Energy (Conserved)', 
             linewidth=2, color='black', linestyle='--')
    
    ax1.set_ylabel('Energy (Hartree)', fontsize=12)
    ax1.set_title('Energy Components vs Time', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # Add text box with energy conservation info
    energy_drift = energy_dict['universe_total'][-1] - energy_dict['universe_total'][0]
    energy_std = np.std(energy_dict['universe_total'])
    textstr = f'Universe Energy Drift: {energy_drift:.2e} Hartree\nStd Dev: {energy_std:.2e} Hartree'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax1.text(0.02, 0.98, textstr, transform=ax1.transAxes, fontsize=9,
             verticalalignment='top', bbox=props)
    
    # Panel 2: Temperature
    temperature = data[:, columns['temperature']]
    ax2.plot(time_ps, temperature, label='System Temperature', linewidth=2, color='orange')
    
    ax2.set_ylabel('Temperature (K)', fontsize=12)
    ax2.set_title('Temperature vs Time', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    
    # Add horizontal line at target temperature if it's roughly constant
    temp_mean = np.mean(temperature)
    temp_std = np.std(temperature)
    ax2.axhline(y=temp_mean, color='red', linestyle=':', alpha=0.7, 
                label=f'Mean: {temp_mean:.1f} ± {temp_std:.1f} K')
    ax2.legend(fontsize=10)
    
    # Panel 3: Fictive temperatures
    if fictive_temps:
        colors = {'total_PE': 'purple', 'harmonic': 'brown', 'lj': 'pink', 'coulombic': 'cyan'}
        labels = {'total_PE': 'Total PE', 'harmonic': 'Harmonic', 'lj': 'Lennard-Jones', 'coulombic': 'Coulombic'}
        
        for comp_name, fictive_T in fictive_temps.items():
            if not np.all(np.isnan(fictive_T)):
                # Plot all valid fictive temperatures
                valid_mask = ~np.isnan(fictive_T)
                ax3.plot(time_ps[valid_mask], fictive_T[valid_mask], 
                        label=f'{labels.get(comp_name, comp_name)} Fictive T',
                        linewidth=2, color=colors.get(comp_name, 'gray'), alpha=0.8)
                
                # Mark extrapolated regions if we have interpolator info
                if comp_name in interpolators:
                    U_range = interpolators[comp_name]['U_range']
                    energy_data = energy_dict[comp_name]
                    extrapolated_mask = ((energy_data < U_range[0]) | (energy_data > U_range[1])) & valid_mask
                    
                    if np.any(extrapolated_mask):
                        # Plot extrapolated points with different style
                        ax3.scatter(time_ps[extrapolated_mask], fictive_T[extrapolated_mask],
                                  color=colors.get(comp_name, 'gray'), alpha=0.3, s=1,
                                  label=f'{labels.get(comp_name, comp_name)} (Extrapolated)')
        
        ax3.set_ylabel('Fictive Temperature (K)', fontsize=12)
        ax3.set_title('Fictive Temperatures vs Time\n(Derived from Energy Components using U(T) relationships)', 
                     fontsize=14, fontweight='bold')
        ax3.legend(fontsize=9, ncol=2)
        ax3.grid(True, alpha=0.3)
        
        # Add reference line at actual mean temperature
        ax3.axhline(y=temp_mean, color='red', linestyle=':', alpha=0.7, 
                   label=f'Actual Mean T: {temp_mean:.1f} K')
        ax3.legend(fontsize=9, ncol=2)
        
        # Add text box with extrapolation info
        total_points = len(time_ps)
        total_extrapolated = 0
        for comp_name in fictive_temps:
            if comp_name in interpolators:
                U_range = interpolators[comp_name]['U_range']
                energy_data = energy_dict[comp_name]
                extrapolated_mask = (energy_data < U_range[0]) | (energy_data > U_range[1])
                total_extrapolated += np.sum(extrapolated_mask)
        
        if total_extrapolated > 0:
            avg_extrapolated_pct = 100 * total_extrapolated / (len(fictive_temps) * total_points)
            extrap_textstr = f'Avg. {avg_extrapolated_pct:.1f}% extrapolated\n(shown as scattered points)'
            extrap_props = dict(boxstyle='round', facecolor='lightyellow', alpha=0.8)
            ax3.text(0.98, 0.02, extrap_textstr, transform=ax3.transAxes, fontsize=8,
                    verticalalignment='bottom', horizontalalignment='right', bbox=extrap_props)
    else:
        ax3.text(0.5, 0.5, 'Fictive temperature calculation unavailable\n(missing temperature components file)', 
                ha='center', va='center', transform=ax3.transAxes, fontsize=12)
        ax3.set_ylabel('Fictive Temperature (K)', fontsize=12)
        ax3.set_title('Fictive Temperatures vs Time', fontsize=14, fontweight='bold')
    
    ax3.set_xlabel('Time (ps)', fontsize=12)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save or show plot
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_file}")
    else:
        plt.show()
    
    return fig

def print_energy_statistics(data, columns):
    """Print useful statistics about the energy data."""
    
    energy_dict = calculate_energy_components(data, columns)
    time_ps = data[:, columns['time(ps)']]
    temperature = data[:, columns['temperature']]
    
    print("\n" + "="*60)
    print("ENERGY STATISTICS")
    print("="*60)
    
    print(f"Simulation time: {time_ps[0]:.6f} - {time_ps[-1]:.6f} ps")
    print(f"Number of data points: {len(data)}")
    
    print(f"\nEnergy Components (Hartree):")
    for comp_name, energies in energy_dict.items():
        if comp_name.endswith('_total') or comp_name == 'molecular_potential':
            print(f"  {comp_name.replace('_', ' ').title()}:")
            print(f"    Initial: {energies[0]:.6f}")
            print(f"    Final:   {energies[-1]:.6f}")
            print(f"    Change:  {energies[-1] - energies[0]:.6f}")
    
    print(f"\nTemperature Statistics:")
    print(f"  Mean: {np.mean(temperature):.2f} K")
    print(f"  Std:  {np.std(temperature):.2f} K")
    print(f"  Min:  {np.min(temperature):.2f} K")
    print(f"  Max:  {np.max(temperature):.2f} K")
    
    print("="*60)

def main():
    parser = argparse.ArgumentParser(
        description='Plot energy tracker output with energy components, temperature, and fictive temperatures',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument('input_file', type=str, 
                       help='Energy tracker output file (e.g., prod-0_energy_tracker.txt)')
    parser.add_argument('--output', '-o', type=str, 
                       help='Output plot file (default: show plot)')
    parser.add_argument('--stats', action='store_true', 
                       help='Print detailed energy statistics')
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
    interpolators = {}
    if temp_data is not None:
        print("Creating energy-temperature interpolators...")
        interpolators = create_energy_temperature_interpolators(temp_data, temp_columns)
        if interpolators:
            print(f"Successfully created {len(interpolators)} interpolators")
        else:
            print("No interpolators could be created")
    
    # Print statistics if requested
    if args.stats:
        print_energy_statistics(data, columns)
    
    # Create plot
    print("\nCreating plot...")
    fig = plot_energy_temperature_and_fictive(data, columns, interpolators, args.output)
    
    if not args.output:
        print("Close the plot window to exit.")

if __name__ == '__main__':
    main()
