#!/usr/bin/env python3
"""
Rabi Splitting Simulation
- 250 molecules
- 40 Bohr periodic box
- Cavity at 1560 cm⁻¹
- λ = 0.07 (strong coupling)
"""

from openmm import openmm
from openmm import unit
import numpy as np
import time

# Unit conversions
BOHR_TO_NM = 0.0529177
HARTREE_TO_KJMOL = 2625.5
HARTREE_TO_CM = 219474.63

# Simulation parameters
NUM_MOLECULES = 250
BOX_SIZE_BOHR = 40.0
BOX_SIZE_NM = BOX_SIZE_BOHR * BOHR_TO_NM  # ~2.12 nm

# Force constant for O-O at 1560 cm⁻¹ (NO doubling!)
K_OO_AU = 0.73204  # Hartree/Bohr² - gives 1560 cm⁻¹
R0_OO_AU = 2.281655158  # Bohr
K_OO = K_OO_AU * HARTREE_TO_KJMOL / (BOHR_TO_NM**2)  # kJ/(mol·nm²)
R0_OO = R0_OO_AU * BOHR_TO_NM  # nm

# Cavity parameters
CAVITY_FREQ_CM = 1560  # cm⁻¹ - resonant with O-O stretch
OMEGAC_AU = CAVITY_FREQ_CM / HARTREE_TO_CM
LAMBDA_COUPLING = 0.07
# Photon mass: cav-hoomd uses 1.0 a.u. (electron masses)
# 1 a.u. of mass = 1/1822.888 amu
AMU_TO_AU = 1822.888
PHOTON_MASS = 1.0 / AMU_TO_AU  # a.u. -> amu (matches cav-hoomd's phmass=1.0 a.u.)

# Molecular parameters
MASS_O = 16.0  # amu
CHARGE_MAGNITUDE = 0.3  # e
TEMPERATURE_K = 100.0

# Timing
DT_PS = 0.001  # 1 fs timestep
FRICTION = 0.01  # ps⁻¹
EQUIL_STEPS = 50000   # 50 ps equilibration
PROD_STEPS = 300000   # 300 ps production


def create_system():
    """Create the molecular system with cavity coupling."""
    
    system = openmm.System()
    system.setDefaultPeriodicBoxVectors(
        openmm.Vec3(BOX_SIZE_NM, 0, 0) * unit.nanometer,
        openmm.Vec3(0, BOX_SIZE_NM, 0) * unit.nanometer,
        openmm.Vec3(0, 0, BOX_SIZE_NM) * unit.nanometer
    )
    
    positions = []
    charges = []
    bond_force = openmm.HarmonicBondForce()
    nonbonded = openmm.NonbondedForce()
    nonbonded.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
    nonbonded.setCutoffDistance(1.0 * unit.nanometer)
    
    # Lattice placement
    np.random.seed(42)
    side = int(np.ceil(NUM_MOLECULES**(1/3)))
    spacing = BOX_SIZE_NM / side
    
    print(f"Lattice: {side}x{side}x{side}, spacing = {spacing:.4f} nm")
    
    mol_idx = 0
    for i in range(side):
        for j in range(side):
            for k in range(side):
                if mol_idx >= NUM_MOLECULES:
                    break
                
                # Center of molecule
                x = (i + 0.5) * spacing
                y = (j + 0.5) * spacing
                z = (k + 0.5) * spacing
                
                # Random orientation
                theta = np.random.rand() * 2 * np.pi
                phi = np.random.rand() * np.pi
                direction = np.array([
                    np.sin(phi) * np.cos(theta),
                    np.sin(phi) * np.sin(theta),
                    np.cos(phi)
                ])
                
                # Two atoms of the dimer
                pos1 = np.array([x, y, z]) - 0.5 * R0_OO * direction
                pos2 = np.array([x, y, z]) + 0.5 * R0_OO * direction
                
                # Add particles
                system.addParticle(MASS_O)
                system.addParticle(MASS_O)
                
                positions.append(openmm.Vec3(*pos1) * unit.nanometer)
                positions.append(openmm.Vec3(*pos2) * unit.nanometer)
                
                # Add bond
                bond_force.addBond(mol_idx*2, mol_idx*2+1, R0_OO, K_OO)
                
                # Add charges (opposite signs for dipole)
                nonbonded.addParticle(-CHARGE_MAGNITUDE, 0.3, 0.5)
                nonbonded.addParticle(+CHARGE_MAGNITUDE, 0.3, 0.5)
                nonbonded.addException(mol_idx*2, mol_idx*2+1, 0.0, 1.0, 0.0)
                
                charges.extend([-CHARGE_MAGNITUDE, +CHARGE_MAGNITUDE])
                mol_idx += 1
                
            if mol_idx >= NUM_MOLECULES:
                break
        if mol_idx >= NUM_MOLECULES:
            break
    
    system.addForce(bond_force)
    system.addForce(nonbonded)
    num_molecular = len(positions)
    
    print(f"Created {mol_idx} O-O molecules ({num_molecular} atoms)")
    
    # Add cavity particle (field variable, not a real particle)
    cavity_index = system.addParticle(PHOTON_MASS)
    positions.append(openmm.Vec3(0, 0, 0) * unit.nanometer)
    nonbonded.addParticle(0.0, 0.1, 0.0)  # No charge, small sigma, no LJ
    
    # CavityForce - coupling is ON from the start (lambda = 0.07)
    # Equilibration runs under full cavity potential
    cavity_force = openmm.CavityForce(cavity_index, OMEGAC_AU, LAMBDA_COUPLING, PHOTON_MASS)
    system.addForce(cavity_force)
    
    # CavityParticleDisplacer for finite-q correction (active from step 0 with lambda)
    displacer = openmm.CavityParticleDisplacer(cavity_index, OMEGAC_AU, PHOTON_MASS)
    displacer.setSwitchOnStep(0)  # Active from the start
    displacer.setSwitchOnLambda(LAMBDA_COUPLING)
    system.addForce(displacer)
    
    return system, positions, charges, num_molecular, cavity_index


def compute_dipole(state, charges, n_atoms):
    """Compute total molecular dipole moment."""
    positions = state.getPositions(asNumpy=True)
    dipole = np.zeros(3)
    for i in range(n_atoms):
        pos = positions[i].value_in_unit(unit.nanometer)
        dipole += charges[i] * np.array(pos)
    return dipole


def main():
    print("=" * 70)
    print("Rabi Splitting Simulation")
    print("=" * 70)
    print(f"\nParameters:")
    print(f"  Molecules: {NUM_MOLECULES}")
    print(f"  Box size: {BOX_SIZE_BOHR} Bohr = {BOX_SIZE_NM:.3f} nm")
    print(f"  O-O force constant: k = {K_OO_AU} a.u. = {K_OO:.0f} kJ/(mol·nm²)")
    print(f"  Cavity: ω_c = {CAVITY_FREQ_CM} cm⁻¹, λ = {LAMBDA_COUPLING}")
    print(f"  Photon mass: {PHOTON_MASS:.6f} amu = {PHOTON_MASS * AMU_TO_AU:.1f} a.u. (matches cav-hoomd)")
    print(f"  Temperature: {TEMPERATURE_K} K")
    
    # Create system
    system, positions, charges, num_molecular, cavity_index = create_system()
    
    # Create integrator - Langevin throughout
    integrator = openmm.LangevinMiddleIntegrator(
        TEMPERATURE_K * unit.kelvin,
        FRICTION / unit.picosecond,
        DT_PS * unit.picoseconds
    )
    
    # Create context - use CUDA for GPU acceleration
    platform = openmm.Platform.getPlatformByName('CUDA')
    properties = {'Precision': 'mixed'}
    context = openmm.Context(system, integrator, platform, properties)
    context.setPositions(positions)
    context.setVelocitiesToTemperature(TEMPERATURE_K * unit.kelvin)
    
    # Minimize
    print("\nMinimizing energy...", flush=True)
    openmm.LocalEnergyMinimizer.minimize(context, maxIterations=1000)
    
    state = context.getState(getEnergy=True)
    pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    print(f"Initial PE after minimization: {pe:.1f} kJ/mol", flush=True)
    
    # Storage for trajectory
    dipole_times = []
    dipole_values = []
    cavity_values = []  # Also save cavity position
    
    total_steps = EQUIL_STEPS + PROD_STEPS
    print(f"\nRunning {total_steps} steps...", flush=True)
    print(f"  Equilibration: {EQUIL_STEPS} steps ({EQUIL_STEPS*DT_PS} ps) - coupling ON from start!", flush=True)
    print(f"  Production: {PROD_STEPS} steps ({PROD_STEPS*DT_PS} ps)", flush=True)
    
    start = time.time()
    
    for step in range(total_steps):
        integrator.step(1)
        
        # Record dipole and cavity during production only
        if step >= EQUIL_STEPS:
            state = context.getState(getPositions=True)
            pos = state.getPositions(asNumpy=True)
            dipole = compute_dipole(state, charges, num_molecular)
            q_cav = pos[cavity_index].value_in_unit(unit.nanometer)
            
            dipole_times.append((step - EQUIL_STEPS) * DT_PS)
            dipole_values.append(dipole.copy())
            cavity_values.append(q_cav.copy())
        
        # Progress report every 1 ps (1000 steps at 1 fs timestep)
        # Fix error: Only print and save every 100 or 1000 steps, not every step!
        # Check rate to avoid division by zero
        if (step + 1) % 1000 == 0 or (step + 1) == total_steps:  # Print every 100 steps or at the last step
            pct = 100 * (step + 1) / total_steps
            elapsed = time.time() - start
            rate = (step + 1) / max(elapsed, 1e-8)  # Avoid division by zero
            eta = (total_steps - step - 1) / max(rate, 1e-8)
            phase = "EQUIL" if step < EQUIL_STEPS else "PROD"

            state = context.getState(getEnergy=True, getPositions=True)
            pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            pos = state.getPositions(asNumpy=True)
            q_cav = pos[cavity_index].value_in_unit(unit.nanometer)
            dipole_now = compute_dipole(state, charges, num_molecular)

            if phase == "PROD":
                np.savez('rabi_250mol_40bohr.npz',
                        time_ps=np.array(dipole_times),
                        dipole_nm=np.array(dipole_values),
                        cavity_nm=np.array(cavity_values),
                        metadata={
                            'lambda_coupling': LAMBDA_COUPLING,
                            'cavity_freq_cm': CAVITY_FREQ_CM,
                            'molecular_freq_cm': 1560,
                            'dt_ps': DT_PS,
                            'num_molecules': NUM_MOLECULES,
                            'box_size_bohr': BOX_SIZE_BOHR,
                            'k_OO_au': K_OO_AU,
                            'temperature_K': TEMPERATURE_K,
                            'status': 'running'
                        })
            print(f"[{pct:5.1f}%] {phase} | PE: {pe:.1f} | "
                  f"d_xy: ({dipole_now[0]:.3f}, {dipole_now[1]:.3f}) | "
                  f"q_cav: ({q_cav[0]:.3f}, {q_cav[1]:.3f}) | "
                  f"ETA: {eta/60:.1f}m", flush=True)
                
    elapsed = time.time() - start
    print(f"\nCompleted in {elapsed/60:.1f} min ({elapsed/total_steps*1e6:.1f} μs/step)", flush=True)
    
    # Final save
    np.savez('rabi_250mol_40bohr.npz',
             time_ps=np.array(dipole_times),
             dipole_nm=np.array(dipole_values),
             cavity_nm=np.array(cavity_values),
             metadata={
                 'lambda_coupling': LAMBDA_COUPLING,
                 'cavity_freq_cm': CAVITY_FREQ_CM,
                 'molecular_freq_cm': 1560,
                 'dt_ps': DT_PS,
                 'num_molecules': NUM_MOLECULES,
                 'box_size_bohr': BOX_SIZE_BOHR,
                 'k_OO_au': K_OO_AU,
                 'temperature_K': TEMPERATURE_K,
                 'status': 'complete'
             })
    
    print(f"Saved: rabi_250mol_40bohr.npz", flush=True)


if __name__ == "__main__":
    main()
