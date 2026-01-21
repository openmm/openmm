#!/usr/bin/env python3
"""
Diagnose why bonds are not vibrating
"""
import numpy as np
import matplotlib.pyplot as plt

# Load trajectory
data = np.load('cavity_diamer_dipole.npz')
time = data['time_ps']

# We need to load positions to check bond lengths
# But we only saved dipole... Let's run a quick test simulation instead

from openmm import openmm
from openmm import unit

print("Testing bond vibrations in a minimal system...")
print("="*60)

# Create minimal 2-atom O-O system
system = openmm.System()
system.addParticle(16.0)  # O atom 1
system.addParticle(16.0)  # O atom 2

# Add harmonic bond
BOHR_TO_NM = 0.0529177
HARTREE_TO_KJMOL = 2625.5
k_OO_au = 0.732
r0_OO_au = 2.28
k_OO = k_OO_au * HARTREE_TO_KJMOL / (BOHR_TO_NM**2)
r0_OO = r0_OO_au * BOHR_TO_NM

bond_force = openmm.HarmonicBondForce()
bond_force.addBond(0, 1, r0_OO, k_OO)
system.addForce(bond_force)

print(f"Bond parameters:")
print(f"  r0 = {r0_OO:.4f} nm")
print(f"  k = {k_OO:.2f} kJ/(mol·nm²)")

# No periodic box for this test
integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
platform = openmm.Platform.getPlatformByName('Reference')
context = openmm.Context(system, integrator, platform)

# Set initial positions slightly displaced from equilibrium
positions = [
    openmm.Vec3(0, 0, 0) * unit.nanometer,
    openmm.Vec3(r0_OO * 1.05, 0, 0) * unit.nanometer  # 5% stretched
]
context.setPositions(positions)

# Set zero initial velocities
velocities = [
    openmm.Vec3(0, 0, 0) * unit.nanometer/unit.picosecond,
    openmm.Vec3(0, 0, 0) * unit.nanometer/unit.picosecond
]
context.setVelocities(velocities)

# Run for a few periods
n_steps = 100  # 100 fs
bond_lengths = []
times = []

for i in range(n_steps):
    integrator.step(1)
    state = context.getState(getPositions=True)
    pos = state.getPositions(asNumpy=True)
    r = np.linalg.norm(pos[1] - pos[0])
    bond_lengths.append(r)
    times.append(i * 0.001)  # ps

bond_lengths = np.array(bond_lengths)
times = np.array(times)

print(f"\nBond length oscillation:")
print(f"  Mean: {np.mean(bond_lengths):.4f} nm")
print(f"  Std: {np.std(bond_lengths):.4f} nm")
print(f"  Min: {np.min(bond_lengths):.4f} nm")
print(f"  Max: {np.max(bond_lengths):.4f} nm")
print(f"  r0: {r0_OO:.4f} nm")

# FFT to find frequency
bond_centered = bond_lengths - np.mean(bond_lengths)
fft = np.fft.rfft(bond_centered)
freqs_hz = np.fft.rfftfreq(len(bond_centered), 0.001e-12)
freqs_cm = freqs_hz / 3e10
power = np.abs(fft)**2

peak_idx = np.argmax(power[1:]) + 1  # Skip DC
print(f"\nMeasured vibrational frequency: {freqs_cm[peak_idx]:.1f} cm⁻¹")
print(f"Expected: 1554 cm⁻¹")

# Plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

ax1.plot(times, bond_lengths * 10, 'b-', linewidth=2)  # Convert to Angstroms
ax1.axhline(r0_OO * 10, color='r', linestyle='--', label=f'r₀ = {r0_OO*10:.3f} Å')
ax1.set_xlabel('Time (ps)')
ax1.set_ylabel('Bond Length (Å)')
ax1.set_title('O-O Bond Vibration (Isolated Dimer, No Thermostat)')
ax1.grid(True, alpha=0.3)
ax1.legend()

ax2.plot(freqs_cm[:len(freqs_cm)//2], power[:len(power)//2], 'r-', linewidth=2)
ax2.set_xlabel('Frequency (cm⁻¹)')
ax2.set_ylabel('Power')
ax2.set_title('Vibrational Spectrum')
ax2.set_xlim(0, 3000)
ax2.grid(True, alpha=0.3)
ax2.axvline(1554, color='g', linestyle='--', label='Expected: 1554 cm⁻¹')
ax2.legend()

plt.tight_layout()
plt.savefig('bond_vibration_test.png', dpi=150)
print(f"\nSaved diagnostic plot: bond_vibration_test.png")

if np.std(bond_lengths) < 1e-5:
    print("\n⚠️  WARNING: Bond is NOT vibrating!")
    print("This suggests the bond force is not working properly.")
elif abs(freqs_cm[peak_idx] - 1554) > 100:
    print(f"\n⚠️  WARNING: Frequency mismatch!")
    print(f"Measured {freqs_cm[peak_idx]:.1f} cm⁻¹ vs expected 1554 cm⁻¹")
else:
    print("\n✓ Bond vibration looks correct!")
    print("The problem must be in the full system simulation.")
