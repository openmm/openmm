#!/usr/bin/env python3
"""
Quick diagnostic test to check cavity force and dipole
"""
import sys
import numpy as np
from openmm import openmm, unit

BOHR_TO_NM = 0.0529177
HARTREE_TO_KJMOL = 2625.5
AU_TIME_TO_PS = 0.02418884254

# Create simple 2-particle system
system = openmm.System()
positions = []

# Two charged particles
system.addParticle(1.0)  # amu
system.addParticle(1.0)

positions.append(openmm.Vec3(0.0, 0.0, 0.0) * unit.nanometer)
positions.append(openmm.Vec3(0.1, 0.0, 0.0) * unit.nanometer)

# Add nonbonded force
nb = openmm.NonbondedForce()
nb.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
nb.addParticle(-0.3, 0.1, 0.0)  # charge=-0.3e, sigma=0.1nm, epsilon=0
nb.addParticle(+0.3, 0.1, 0.0)  # charge=+0.3e
system.addForce(nb)

# Add cavity particle
cavity_idx = system.addParticle(1.0)
positions.append(openmm.Vec3(0.0, 0.0, 0.0) * unit.nanometer)
nb.addParticle(0.0, 0.1, 0.0)

# Add cavity force
omegac_au = 0.00913
lambda_coupling = 0.001
cavity_force = openmm.CavityForce(cavity_idx, omegac_au, lambda_coupling, 1.0)
system.addForce(cavity_force)

# Periodic box
system.setDefaultPeriodicBoxVectors(
    openmm.Vec3(10, 0, 0),
    openmm.Vec3(0, 10, 0),
    openmm.Vec3(0, 0, 10)
)

# Create context
integrator = openmm.VerletIntegrator(0.001 * unit.picosecond)
platform = openmm.Platform.getPlatformByName('CUDA')
context = openmm.Context(system, integrator, platform)
context.setPositions(positions)

# Get initial state
state = context.getState(getForces=True, getEnergy=True, getPositions=True)
forces = state.getForces(asNumpy=True)
energy = state.getPotentialEnergy()
pos = state.getPositions(asNumpy=True)

print(f"Positions:")
for i, p in enumerate(pos):
    print(f"  Particle {i}: {p}")

print(f"\nForces:")
for i, f in enumerate(forces):
    print(f"  Particle {i}: {f}")

print(f"\nTotal energy: {energy}")

# Calculate expected dipole
dipole_x = -0.3 * 0.0 + 0.3 * 0.1  # e*nm
dipole_y = 0.0
print(f"\nExpected dipole: ({dipole_x:.6f}, {dipole_y:.6f}) e*nm")

# Expected force on cavity:
# F = -K * q - epsilon * d
K_au = 1.0 * omegac_au**2
K_openmm = K_au * HARTREE_TO_KJMOL / (BOHR_TO_NM**2)
epsilon_au = lambda_coupling * omegac_au
epsilon_openmm = epsilon_au * HARTREE_TO_KJMOL / BOHR_TO_NM

print(f"\nK = {K_openmm:.4f} kJ/(mol*nm^2)")
print(f"epsilon = {epsilon_openmm:.4f} kJ/(mol*nm*e)")

q_cavity = 0.0  # At origin
F_expected_x = -K_openmm * q_cavity - epsilon_openmm * dipole_x
print(f"\nExpected force on cavity: F_x = {F_expected_x:.6f} kJ/(mol*nm)")
print(f"Actual force on cavity: F_x = {forces[cavity_idx][0]:.6f} kJ/(mol*nm)")
