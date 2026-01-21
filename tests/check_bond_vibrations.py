#!/usr/bin/env python3
"""
Check if bonds are actually vibrating in the saved trajectory
"""
import numpy as np

from openmm import openmm
from openmm import unit

print("="*60)
print("Diagnosing: Are bonds actually present and vibrating?")
print("="*60)

# Minimal test: 2-particle O-O system with Bussi thermostat
system = openmm.System()
system.addParticle(16.0)  # O1
system.addParticle(16.0)  # O2

# Add bond
BOHR_TO_NM = 0.0529177
HARTREE_TO_KJMOL = 2625.5
k_OO_au = 0.732
r0_OO_au = 2.28
k_OO = k_OO_au * HARTREE_TO_KJMOL / (BOHR_TO_NM**2)
r0_OO = r0_OO_au * BOHR_TO_NM

bond_force = openmm.HarmonicBondForce()
bond_force.addBond(0, 1, r0_OO, k_OO)
system.addForce(bond_force)

# Add Bussi thermostat
tau = 0.1
bussi = openmm.BussiThermostat(100.0, tau)  # 100 K
system.addForce(bussi)

# Integrator
integrator = openmm.LangevinMiddleIntegrator(100.0 * unit.kelvin, 0.01 / unit.picosecond, 0.001 * unit.picoseconds)
platform = openmm.Platform.getPlatformByName('Reference')
context = openmm.Context(system, integrator, platform)

# Set positions
positions = [
    openmm.Vec3(0, 0, 0) * unit.nanometer,
    openmm.Vec3(r0_OO, 0, 0) * unit.nanometer
]
context.setPositions(positions)

# Give it thermal velocities
context.setVelocitiesToTemperature(100.0 * unit.kelvin)

# Run and measure bond lengths
print(f"\nRunning 1000 steps with Bussi thermostat at 100 K...")
bond_lengths = []
for i in range(1000):
    integrator.step(1)
    if i % 10 == 0:  # Every 10 fs
        state = context.getState(getPositions=True)
        pos = state.getPositions(asNumpy=True)
        r = np.linalg.norm(pos[1] - pos[0])
        bond_lengths.append(r)

bond_lengths = np.array(bond_lengths)

print(f"\nBond length statistics:")
print(f"  Mean: {np.mean(bond_lengths):.6f} nm")
print(f"  Std: {np.std(bond_lengths):.6f} nm")
print(f"  Min: {np.min(bond_lengths):.6f} nm")
print(f"  Max: {np.max(bond_lengths):.6f} nm")
print(f"  r0: {r0_OO:.6f} nm")
print(f"  Oscillation amplitude: {np.std(bond_lengths):.6f} nm")

# Check if it's vibrating
if np.std(bond_lengths) < 1e-5:
    print("\n❌ FATAL: Bond is NOT vibrating at all!")
    print("The bond is frozen at its equilibrium length.")
elif np.std(bond_lengths) < 1e-3:
    print("\n⚠️  WARNING: Very small bond oscillations!")
    print(f"Expected thermal amplitude at 100 K: ~0.001-0.002 nm")
else:
    print("\n✓ Bond is vibrating")

# FFT to check frequency
from scipy.fft import rfft, rfftfreq
bond_centered = bond_lengths - np.mean(bond_lengths)
fft_result = np.abs(rfft(bond_centered))
freqs_hz = rfftfreq(len(bond_centered), 10e-15)  # 10 fs sampling
freqs_cm = freqs_hz / 3e10

if len(fft_result) > 1:
    peak_idx = np.argmax(fft_result[1:]) + 1
    print(f"\nDominant frequency: {freqs_cm[peak_idx]:.1f} cm⁻¹")
    print(f"Expected: ~1555 cm⁻¹")
