#!/usr/bin/env python3
"""
Minimal test to verify CavityForce can read charges and compute coupling energy.
"""

import numpy as np
try:
    from openmm import openmm
    from openmm import unit
except:
    import openmm
    from openmm import unit

print("="*60)
print("Testing CavityForce Charge Reading")
print("="*60)

# Create minimal system: 2 particles + cavity
system = openmm.System()
system.addParticle(16.0)  # Particle 0: O atom with -0.5e
system.addParticle(16.0)  # Particle 1: O atom with +0.5e
system.addParticle(1.0)   # Particle 2: Cavity photon

# Set box
system.setDefaultPeriodicBoxVectors(
    openmm.Vec3(2.0, 0, 0),
    openmm.Vec3(0, 2.0, 0),
    openmm.Vec3(0, 0, 2.0)
)

# Add NonbondedForce with charges
nonbonded = openmm.NonbondedForce()
nonbonded.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
nonbonded.setCutoffDistance(0.9)  # nm

# Add particles with charges
nonbonded.addParticle(-0.5, 0.3, 0.5)  # Particle 0: -0.5e
nonbonded.addParticle(+0.5, 0.3, 0.5)  # Particle 1: +0.5e
nonbonded.addParticle(0.0, 0.1, 0.0)   # Particle 2: cavity (no charge)

system.addForce(nonbonded)

print("\nSystem setup:")
print(f"  3 particles (2 charged + 1 cavity)")
print(f"  Particle 0: -0.5e at position to be set")
print(f"  Particle 1: +0.5e at position to be set")
print(f"  Particle 2: cavity photon")

# Add CavityForce
cavity_index = 2
omegac_au = 0.007  # ~1540 cm-1
lambda_coupling = 0.1
photon_mass = 1.0

cavity_force = openmm.CavityForce(cavity_index, omegac_au, lambda_coupling, photon_mass)
system.addForce(cavity_force)

print(f"\nCavityForce setup:")
print(f"  Cavity index: {cavity_index}")
print(f"  ω_c: {omegac_au} a.u.")
print(f"  λ: {lambda_coupling}")
print(f"  Photon mass: {photon_mass} amu")

# Create context
integrator = openmm.LangevinMiddleIntegrator(100*unit.kelvin, 0.01/unit.picosecond, 0.001*unit.picosecond)
platform = openmm.Platform.getPlatformByName('CUDA')
context = openmm.Context(system, integrator, platform)

# Test Case 1: Particles separated along x-axis → dipole along x
print("\n" + "="*60)
print("Test Case 1: Dipole along x-axis")
print("="*60)
positions = [
    openmm.Vec3(0.0, 0.0, 0.0),  # -0.5e at origin
    openmm.Vec3(0.2, 0.0, 0.0),  # +0.5e at x=0.2 nm
    openmm.Vec3(0.05, 0.0, 0.0), # Cavity at x=0.05 nm (displaced)
]
context.setPositions(positions)

# Calculate expected dipole manually
dipole_manual = np.array([
    -0.5 * 0.0 + 0.5 * 0.2,  # x-component
    0.0,                      # y-component
    0.0                       # z-component
])
print(f"\nManual dipole calculation:")
print(f"  μ = q₀*r₀ + q₁*r₁")
print(f"  μ = (-0.5)*(0.0) + (+0.5)*(0.2)")
print(f"  μ = ({dipole_manual[0]:.3f}, {dipole_manual[1]:.3f}, {dipole_manual[2]:.3f}) e·nm")
print(f"  |μ| = {np.linalg.norm(dipole_manual):.3f} e·nm")

# Get energies from CavityForce
harmonic_e = cavity_force.getHarmonicEnergy(context).value_in_unit(unit.kilojoule_per_mole)
coupling_e = cavity_force.getCouplingEnergy(context).value_in_unit(unit.kilojoule_per_mole)
dipole_e = cavity_force.getDipoleSelfEnergy(context).value_in_unit(unit.kilojoule_per_mole)

print(f"\nCavityForce energies:")
print(f"  Harmonic: {harmonic_e:.6f} kJ/mol")
print(f"  Coupling: {coupling_e:.6f} kJ/mol")
print(f"  Dipole self: {dipole_e:.6f} kJ/mol")

# Expected coupling energy:
# E_coupling = λ * ω_c * (μ · q_photon)
# where q_photon is the cavity displacement (x, y components matter for 2D cavity)
cavity_pos_x = 0.05  # nm
expected_coupling = lambda_coupling * omegac_au * dipole_manual[0] * cavity_pos_x
print(f"\nExpected coupling (rough estimate):")
print(f"  E ≈ λ * ω_c * μ_x * q_x")
print(f"  E ≈ {lambda_coupling} * {omegac_au} * {dipole_manual[0]:.3f} * {cavity_pos_x}")
print(f"  E ≈ {expected_coupling:.6f} (atomic units)")

print("\n" + "="*60)
if abs(coupling_e) < 0.01:
    print("❌ FAIL: Coupling energy is zero!")
    print("   CavityForce is NOT reading charges from NonbondedForce")
else:
    print("✓ PASS: Coupling energy is non-zero")
    print("   CavityForce IS reading charges from NonbondedForce")
print("="*60)

del context
