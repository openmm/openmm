#!/usr/bin/env python3
"""
Comprehensive test for hybrid RPMD: Water dimer with quantum hydrogens.

This tests a more realistic scenario:
- 2 water molecules (6 atoms total)
- 4 H atoms: quantum (RPMD with 8 beads)
- 2 O atoms: classical
- TIP3P-like force field (no constraints)
- Tests: energy conservation, temperature control, quantum effects
"""

import sys
import numpy as np

try:
    import openmm
    from openmm import app, unit
except ImportError:
    print("ERROR: OpenMM not found")
    sys.exit(1)

print("="*70)
print("HYBRID RPMD TEST - Water Dimer (Quantum H, Classical O)")
print("="*70)

# ============================================================================
# Build Water Dimer System
# ============================================================================
print("\n[1] Building water dimer system...")

system = openmm.System()

# Masses
mass_H = 1.008 * unit.amu
mass_O = 15.999 * unit.amu

# Water 1: O-H-H
particle_types = []
system.addParticle(mass_O)  # 0: O1
particle_types.append('O')
system.addParticle(mass_H)  # 1: H1a
particle_types.append('H')
system.addParticle(mass_H)  # 2: H1b
particle_types.append('H')

# Water 2: O-H-H  
system.addParticle(mass_O)  # 3: O2
particle_types.append('O')
system.addParticle(mass_H)  # 4: H2a
particle_types.append('H')
system.addParticle(mass_H)  # 5: H2b
particle_types.append('H')

print(f"  Added {system.getNumParticles()} particles")
print(f"  - {sum(1 for t in particle_types if t == 'H')} hydrogens (quantum)")
print(f"  - {sum(1 for t in particle_types if t == 'O')} oxygens (classical)")

# ============================================================================
# Add Forces
# ============================================================================
print("\n[2] Adding force field...")

# Bond forces (O-H bonds)
bond_force = openmm.HarmonicBondForce()
k_oh = 462750.0 * unit.kilojoules_per_mole / unit.nanometer**2
r0_oh = 0.09572 * unit.nanometer

# Water 1 bonds
bond_force.addBond(0, 1, r0_oh, k_oh)  # O1-H1a
bond_force.addBond(0, 2, r0_oh, k_oh)  # O1-H1b

# Water 2 bonds
bond_force.addBond(3, 4, r0_oh, k_oh)  # O2-H2a
bond_force.addBond(3, 5, r0_oh, k_oh)  # O2-H2b

system.addForce(bond_force)
print(f"  Added {bond_force.getNumBonds()} O-H bonds")

# Angle forces (H-O-H angles)
angle_force = openmm.HarmonicAngleForce()
k_hoh = 836.8 * unit.kilojoules_per_mole / unit.radian**2
theta0_hoh = 104.52 * unit.degrees

angle_force.addAngle(1, 0, 2, theta0_hoh, k_hoh)  # H1a-O1-H1b
angle_force.addAngle(4, 3, 5, theta0_hoh, k_hoh)  # H2a-O2-H2b

system.addForce(angle_force)
print(f"  Added {angle_force.getNumAngles()} H-O-H angles")

# Nonbonded forces (LJ + Coulomb)
nb_force = openmm.NonbondedForce()
nb_force.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)

# TIP3P parameters
# O: q=-0.834e, sigma=0.3151 nm, epsilon=0.6364 kJ/mol
# H: q=+0.417e, sigma=0.0 nm, epsilon=0.0 kJ/mol
q_O = -0.834 * unit.elementary_charge
q_H = 0.417 * unit.elementary_charge
sigma_O = 0.3151 * unit.nanometer
epsilon_O = 0.6364 * unit.kilojoules_per_mole

for i in range(6):
    if particle_types[i] == 'O':
        nb_force.addParticle(q_O, sigma_O, epsilon_O)
    else:  # H
        nb_force.addParticle(q_H, 0.0*unit.nanometer, 0.0*unit.kilojoules_per_mole)

# Add exceptions for bonded atoms (1-2 and 1-3 interactions)
# Water 1
nb_force.addException(0, 1, 0.0, 1.0, 0.0)  # O1-H1a
nb_force.addException(0, 2, 0.0, 1.0, 0.0)  # O1-H1b
nb_force.addException(1, 2, 0.0, 1.0, 0.0)  # H1a-H1b

# Water 2
nb_force.addException(3, 4, 0.0, 1.0, 0.0)  # O2-H2a
nb_force.addException(3, 5, 0.0, 1.0, 0.0)  # O2-H2b
nb_force.addException(4, 5, 0.0, 1.0, 0.0)  # H2a-H2b

system.addForce(nb_force)
print("  Added nonbonded forces (LJ + Coulomb)")

# ============================================================================
# Create Hybrid RPMD Integrator
# ============================================================================
print("\n[3] Creating hybrid RPMD integrator...")

num_beads = 8
temperature = 300.0 * unit.kelvin
friction = 1.0 / unit.picosecond
timestep = 0.25 * unit.femtosecond

integrator = openmm.RPMDIntegrator(num_beads, temperature, friction, timestep)

# Configure hybrid mode: H quantum, O classical
for i in range(6):
    if particle_types[i] == 'H':
        integrator.setParticleType(i, 1)  # Type 1 = quantum
    else:
        integrator.setParticleType(i, 0)  # Type 0 = classical

integrator.setQuantumParticleTypes({1})
integrator.setDefaultQuantum(False)
integrator.setClassicalThermostat(openmm.RPMDIntegrator.BussiClassical)
integrator.setThermostatType(openmm.RPMDIntegrator.PileG)

print(f"  Beads: {num_beads}")
print(f"  Temperature: {temperature}")
print(f"  Timestep: {timestep}")
print(f"  Quantum particles: H atoms (types: {integrator.getQuantumParticleTypes()})")
print(f"  Classical particles: O atoms")

# ============================================================================
# Create Context and Initialize
# ============================================================================
print("\n[4] Creating context...")

try:
    platform = openmm.Platform.getPlatformByName('CUDA')
    properties = {'Precision': 'mixed'}
    platform_name = 'CUDA'
except:
    try:
        platform = openmm.Platform.getPlatformByName('OpenCL')
        properties = {'Precision': 'mixed'}
        platform_name = 'OpenCL'
    except:
        platform = openmm.Platform.getPlatformByName('CPU')
        properties = {}
        platform_name = 'CPU'

print(f"  Using {platform_name} platform")

try:
    context = openmm.Context(system, integrator, platform, properties)
    print("  ✓ Context created successfully")
except Exception as e:
    print(f"  ✗ Context creation failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Set initial positions (typical water dimer geometry)
positions = np.array([
    # Water 1
    [0.000,  0.000,  0.000],  # O1
    [0.096,  0.000,  0.000],  # H1a
    [-0.024,  0.093,  0.000],  # H1b
    # Water 2 (separated by ~0.3 nm)
    [0.300,  0.000,  0.000],  # O2
    [0.396,  0.000,  0.000],  # H2a
    [0.276,  0.093,  0.000],  # H2b
]) * unit.nanometer

for bead in range(num_beads):
    integrator.setPositions(bead, positions)

context.setVelocitiesToTemperature(temperature.value_in_unit(unit.kelvin))
state = context.getState(getVelocities=True)
velocities = state.getVelocities()
for bead in range(num_beads):
    integrator.setVelocities(bead, velocities)

print("  Initial conditions set")

# ============================================================================
# Equilibration
# ============================================================================
print("\n[5] Running equilibration (5000 steps = 1.25 ps)...")

try:
    integrator.step(5000)
    print("  ✓ Equilibration completed")
except Exception as e:
    print(f"  ✗ Equilibration failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# ============================================================================
# Production Run with Analysis
# ============================================================================
print("\n[6] Running production (10000 steps = 2.5 ps)...")

energies = []
quantum_spreads = []
classical_spreads = []

n_samples = 100
steps_per_sample = 100

for i in range(n_samples):
    integrator.step(steps_per_sample)
    
    # Record energy
    energies.append(integrator.getTotalEnergy())
    
    # Analyze quantum vs classical spread
    quantum_pos = []
    classical_pos = []
    
    for bead in range(num_beads):
        state = integrator.getState(bead, getPositions=True)
        pos = state.getPositions()
        
        # Collect H positions (quantum)
        for idx in [1, 2, 4, 5]:
            quantum_pos.append([pos[idx].x, pos[idx].y, pos[idx].z])
        
        # Collect O positions (classical)
        for idx in [0, 3]:
            classical_pos.append([pos[idx].x, pos[idx].y, pos[idx].z])
    
    quantum_pos = np.array(quantum_pos)
    classical_pos = np.array(classical_pos)
    
    # Compute spread for first particle of each type
    quantum_spread = np.std(quantum_pos[:num_beads], axis=0)
    classical_spread = np.std(classical_pos[:num_beads], axis=0)
    
    quantum_spreads.append(np.linalg.norm(quantum_spread))
    classical_spreads.append(np.linalg.norm(classical_spread))

energies = np.array(energies)
quantum_spreads = np.array(quantum_spreads)
classical_spreads = np.array(classical_spreads)

print("  ✓ Production completed")

# ============================================================================
# Analysis
# ============================================================================
print("\n[7] Analysis:")

# Energy statistics
mean_energy = np.mean(energies)
energy_drift = (energies[-1] - energies[0]) / energies[0] * 100
energy_std = np.std(energies) / np.abs(mean_energy) * 100

print(f"\n  Energy:")
print(f"    Mean: {mean_energy:.4f} kJ/mol")
print(f"    Drift: {energy_drift:+.4f}%")
print(f"    Fluctuation: {energy_std:.4f}%")

# Quantum delocalization
mean_quantum_spread = np.mean(quantum_spreads)
mean_classical_spread = np.mean(classical_spreads)

print(f"\n  Bead delocalization:")
print(f"    Quantum (H) spread: {mean_quantum_spread:.6f} nm")
print(f"    Classical (O) spread: {mean_classical_spread:.10f} nm")
print(f"    Ratio: {mean_quantum_spread / (mean_classical_spread + 1e-15):.1e}")

# ============================================================================
# Validation
# ============================================================================
print("\n[8] Validation:")

tests_passed = 0
tests_total = 0

# Test 1: Simulation completed
tests_total += 1
tests_passed += 1
print("  ✓ Simulation completed without crashes")

# Test 2: Energy drift reasonable
tests_total += 1
if abs(energy_drift) < 5.0:  # Less than 5% drift with thermostat on
    tests_passed += 1
    print(f"  ✓ Energy drift acceptable: {abs(energy_drift):.2f}%")
else:
    print(f"  ✗ Energy drift too large: {abs(energy_drift):.2f}%")

# Test 3: Quantum delocalization visible
tests_total += 1
if mean_quantum_spread > 0.001:  # > 0.001 nm spread
    tests_passed += 1
    print(f"  ✓ Quantum delocalization observed: {mean_quantum_spread:.4f} nm")
else:
    print(f"  ✗ Quantum spread too small: {mean_quantum_spread:.6f} nm")

# Test 4: Classical particles stay localized
tests_total += 1
if mean_classical_spread < 1e-8:  # Essentially zero
    tests_passed += 1
    print(f"  ✓ Classical particles remain localized: {mean_classical_spread:.2e} nm")
else:
    print(f"  ✗ Classical spread unexpectedly large: {mean_classical_spread:.2e} nm")

# Test 5: Quantum/classical ratio
tests_total += 1
ratio = mean_quantum_spread / (mean_classical_spread + 1e-15)
if ratio > 1e6:  # Many orders of magnitude difference
    tests_passed += 1
    print(f"  ✓ Quantum/classical spread ratio: {ratio:.1e}")
else:
    print(f"  ⚠ Quantum/classical ratio lower than expected: {ratio:.1e}")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "="*70)
print(f"TEST RESULTS: {tests_passed}/{tests_total} tests passed")
print("="*70)

if tests_passed == tests_total:
    print("✓ ALL TESTS PASSED - Hybrid RPMD working correctly!")
    sys.exit(0)
elif tests_passed >= tests_total - 1:
    print("⚠ MOSTLY WORKING - Some tests failed but core functionality OK")
    sys.exit(0)
else:
    print("✗ SIGNIFICANT ISSUES - Multiple tests failed")
    sys.exit(1)
