#!/usr/bin/env python3
"""
Modified initialization with equilibrium bond length sampling.

This version samples initial bond lengths from the equilibrium Boltzmann
distribution for harmonic bonds: r ~ Normal(r0, sqrt(k_B*T/k))
"""

import numpy as np
import math
import gsd.hoomd
import argparse
import os

def generate_equilibrium_particles_lattice(box_size, num_particles, 
                                          temperature_K=100.0,
                                          bond_params_OO={'k': 0.73204, 'r0': 2.281655158},
                                          bond_params_NN={'k': 1.4325, 'r0': 2.0743522177},
                                          fraction_O=0.8,
                                          seed=None):
    """
    Generate particles with bond lengths sampled from equilibrium distribution.
    
    For a harmonic bond V(r) = (1/2) * k * (r - r0)², the equilibrium distribution is:
    P(r) ∝ exp(-k*(r-r0)²/(2*k_B*T))
    
    This is a Gaussian: r ~ Normal(r0, σ) where σ = sqrt(k_B*T/k)
    
    Parameters
    ----------
    box_size : float
        Box size in Bohr
    num_particles : int
        Total number of particles (must be even)
    temperature_K : float
        Temperature for sampling bond lengths (K)
    bond_params_OO : dict
        O-O bond parameters: {'k': spring_constant, 'r0': equilibrium_length}
    bond_params_NN : dict
        N-N bond parameters: {'k': spring_constant, 'r0': equilibrium_length}
    fraction_O : float
        Fraction of O-O dimers (vs N-N)
    seed : int, optional
        Random seed for reproducibility
        
    Returns
    -------
    particles : list
        List of particle positions
    bond_types : list
        List of bond types (0=O-O, 1=N-N)
    """
    if seed is not None:
        np.random.seed(seed)
    
    if num_particles % 2 != 0:
        raise ValueError("Number of particles must be even")
    
    # Physical constants
    kB_hartree = 3.167e-6  # Hartree/K
    kBT = kB_hartree * temperature_K
    
    # Calculate standard deviations for bond length sampling
    sigma_OO = np.sqrt(kBT / bond_params_OO['k'])
    sigma_NN = np.sqrt(kBT / bond_params_NN['k'])
    
    print(f"Equilibrium bond length sampling at T = {temperature_K} K:")
    print(f"  O-O: r0 = {bond_params_OO['r0']:.6f} Bohr, σ = {sigma_OO:.6f} Bohr ({sigma_OO*0.529177:.6f} Å)")
    print(f"  N-N: r0 = {bond_params_NN['r0']:.6f} Bohr, σ = {sigma_NN:.6f} Bohr ({sigma_NN*0.529177:.6f} Å)")
    
    # Generate lattice of dimer centers
    side_length = int(math.ceil((num_particles // 2)**(1/3)))
    particles_per_layer = side_length**2
    
    dimer_centers = []
    spacing = box_size / side_length
    offset = -box_size / 2 + spacing / 2
    
    for i in range(num_particles // 2):
        layer = i // particles_per_layer
        row = (i % particles_per_layer) // side_length
        col = (i % particles_per_layer) % side_length
        
        x_center = col * spacing + offset
        y_center = row * spacing + offset
        z_center = layer * spacing + offset
        
        dimer_centers.append((x_center, y_center, z_center))
    
    particles = []
    bond_types = []
    
    # Determine which dimers are O-O vs N-N
    num_dimers = len(dimer_centers)
    num_OO = int(fraction_O * num_dimers)
    
    # Randomly assign dimer types
    dimer_type_indices = np.arange(num_dimers)
    np.random.shuffle(dimer_type_indices)
    OO_indices = set(dimer_type_indices[:num_OO])
    
    # Generate dimers with equilibrium bond lengths
    for idx, center in enumerate(dimer_centers):
        # Determine bond type
        is_OO = (idx in OO_indices)
        
        # Sample bond length from equilibrium distribution
        if is_OO:
            bond_length = np.random.normal(bond_params_OO['r0'], sigma_OO)
            bond_types.append(0)  # O-O
        else:
            bond_length = np.random.normal(bond_params_NN['r0'], sigma_NN)
            bond_types.append(1)  # N-N
        
        # Ensure bond length is positive (very rare to be negative with these parameters)
        bond_length = abs(bond_length)
        
        # Random orientation
        theta = np.random.rand() * 2 * np.pi
        phi = np.arccos(2 * np.random.rand() - 1)
        
        direction = np.array([
            np.sin(phi) * np.cos(theta),
            np.sin(phi) * np.sin(theta),
            np.cos(phi)
        ])
        
        # Place atoms at sampled separation
        r1 = np.array(center) - 0.5 * bond_length * direction
        r2 = np.array(center) + 0.5 * bond_length * direction
        
        particles.append(r1)
        particles.append(r2)
    
    return particles, bond_types


def main(job_dir, replica, incavity, couplstr, frequency, Nmol, temperature_K=100.0, seed=None):
    """
    Create initial configuration with equilibrium bond lengths.
    
    Parameters
    ----------
    job_dir : str
        Output directory
    replica : int
        Replica number
    incavity : bool
        Whether to include cavity photon
    couplstr : float
        Coupling strength
    frequency : float
        Cavity frequency (cm^-1)
    Nmol : int
        Number of molecules
    temperature_K : float
        Temperature for equilibrium sampling (K)
    seed : int, optional
        Random seed
    """
    cwd = os.getcwd()
    os.makedirs(job_dir, exist_ok=True)
    os.chdir(job_dir)
    
    # System parameters
    N_particles = 2 * Nmol
    
    # Use the SAME density as original initlattice.py
    # For 500 atoms (250 molecules), box size is 40 Bohr
    # Density = 500 / (40^3) = 0.0078125 particles/Bohr³
    density = 500.0 / (40.0**3)
    
    volume = N_particles / density
    box_size = volume**(1/3)
    
    print(f"System size: {N_particles} particles")
    print(f"Box size: {box_size:.2f} Bohr")
    print(f"Density: {density:.6f} particles/Bohr³ (same as original initlattice.py)")
    
    # Bond parameters from force field
    bond_params_OO = {'k': 2 * 0.36602, 'r0': 2.281655158}  # Hartree/Bohr², Bohr
    bond_params_NN = {'k': 2 * 0.71625, 'r0': 2.0743522177}
    
    frac = 0.8  # Fraction of O-O dimers
    qcharge = 0.3
    
    # Generate particles with equilibrium bond lengths
    print(f"\nGenerating {Nmol} molecules at T = {temperature_K} K")
    position, bond_types = generate_equilibrium_particles_lattice(
        box_size=box_size,
        num_particles=N_particles,
        temperature_K=temperature_K,
        bond_params_OO=bond_params_OO,
        bond_params_NN=bond_params_NN,
        fraction_O=frac,
        seed=seed
    )
    
    # Create GSD frame
    frame = gsd.hoomd.Frame()
    frame.particles.N = N_particles
    frame.particles.position = position
    
    # Initialize arrays
    frame.bonds.group = []
    frame.particles.typeid = []
    frame.bonds.typeid = []
    frame.particles.diameter = []
    frame.particles.charge = []
    frame.particles.mass = []
    
    # Assign particle properties based on bond types
    for i in range(Nmol):
        bond_type = bond_types[i]
        frame.bonds.group.append([2*i, 2*i+1])
        
        if bond_type == 0:  # O-O
            frame.particles.typeid.extend([0, 0])
            frame.particles.diameter.extend([1.0, 1.0])
            frame.bonds.typeid.append(0)
            frame.particles.mass.extend([29150, 29150])
        else:  # N-N
            frame.particles.typeid.extend([1, 1])
            frame.particles.diameter.extend([1.0, 1.0])
            frame.bonds.typeid.append(1)
            frame.particles.mass.extend([25527, 25527])
        
        frame.particles.charge.extend([-qcharge, qcharge])
    
    # Box and types
    frame.configuration.box = [box_size, box_size, box_size, 0, 0, 0]
    frame.configuration.dimensions = 3
    frame.particles.types = ['O', 'N']
    frame.bonds.N = Nmol
    frame.bonds.types = ['O-O', 'N-N']
    
    # Add cavity photon if needed
    if incavity:
        # Calculate dipole moment
        positions = np.array(frame.particles.position)
        charges = np.array(frame.particles.charge)
        dipmom = np.sum(charges[:, np.newaxis] * positions, axis=0)
        
        # Calculate photon position
        hartree_to_cm_minus1 = 219474.63
        omegac = frequency / hartree_to_cm_minus1
        newpos = -dipmom * couplstr / omegac**2
        newpos[-1] = 0.0
        
        # Wrap to box
        box = np.array(frame.configuration.box[:3])
        image_flags = np.floor((newpos + box/2) / box).astype(int)
        newpos = newpos - image_flags * box
        
        # Add photon particle
        frame.particles.typeid.append(2)
        frame.particles.diameter.append(1.0)
        frame.particles.charge.append(0.0)
        frame.particles.position.append(newpos.tolist())
        frame.particles.mass.append(1.0)
        frame.particles.types.append('L')
        
        # Update images
        images = [[0, 0, 0]] * N_particles
        images.append(image_flags.tolist())
        frame.particles.image = images
        frame.particles.N = N_particles + 1
    
    # Write GSD file
    output_file = f'molecular-{replica}.gsd'
    with gsd.hoomd.open(name=output_file, mode='w') as f:
        f.append(frame)
    
    print(f"\nWrote {output_file}")
    print(f"  Box size: {box_size:.2f} Bohr")
    print(f"  Molecules: {Nmol}")
    print(f"  O-O dimers: {bond_types.count(0)}")
    print(f"  N-N dimers: {bond_types.count(1)}")
    
    os.chdir(cwd)


def parse_args():
    parser = argparse.ArgumentParser(description='Generate molecular lattice with equilibrium bond lengths')
    parser.add_argument('--job-dir', type=str, default='.', help='Output directory')
    parser.add_argument('--replica', type=int, default=0, help='Replica number')
    parser.add_argument('--incavity', action='store_true', help='Include cavity photon')
    parser.add_argument('--coupling', type=float, default=0.0, help='Coupling strength')
    parser.add_argument('--frequency', type=float, default=1560.0, help='Cavity frequency (cm^-1)')
    parser.add_argument('--nmol', type=int, default=250, help='Number of molecules')
    parser.add_argument('--temperature', type=float, default=100.0, help='Temperature for equilibrium sampling (K)')
    parser.add_argument('--seed', type=int, default=None, help='Random seed')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(
        job_dir=args.job_dir,
        replica=args.replica,
        incavity=args.incavity,
        couplstr=args.coupling,
        frequency=args.frequency,
        Nmol=args.nmol,
        temperature_K=args.temperature,
        seed=args.seed
    )
