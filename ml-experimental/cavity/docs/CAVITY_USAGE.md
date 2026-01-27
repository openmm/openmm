# Cavity Particle Usage Guide

## Quick Start

### Installation

The cavity particle implementation is built into OpenMM as part of the core library. No additional installation is needed if you have compiled OpenMM from source in the ml-experimental branch.

### Verify Installation

```python
import openmm
print(hasattr(openmm, 'CavityForce'))  # Should print True
print(hasattr(openmm, 'CavityParticleDisplacer'))  # Should print True
```

## Basic Usage

### Example 1: Water with Cavity

```python
from openmm import app, unit
from ml_experimental.cavity.utils import (
    add_cavity_particle,
    setup_cavity_coupling,
    wavenumber_to_hartree
)

# Load your system
pdb = app.PDBFile('water.pdb')
forcefield = app.ForceField('tip4pew.xml')
system = forcefield.createSystem(pdb.topology)

# Add cavity particle
omegac_cm = 3663.0  # OH stretch frequency
omegac_au = wavenumber_to_hartree(omegac_cm)
positions = list(pdb.positions)

cavity_idx = add_cavity_particle(system, positions, omegac_au)

# Setup coupling
lambda_coupling = 0.01  # Coupling strength
cavity_force, displacer = setup_cavity_coupling(
    system, cavity_idx, omegac_au, lambda_coupling
)

# Create integrator and context
integrator = openmm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 
                                       0.5*unit.femtoseconds)
context = openmm.Context(system, integrator)
context.setPositions(positions)

# Run simulation
integrator.step(10000)
```

### Example 2: With RPMD

```python
from openmm import RPMDIntegrator

# Setup system with cavity (as above)
# ...

# Use RPMD integrator for quantum nuclear effects
num_beads = 8
temperature = 300  # K
friction = 1.0  # ps^-1
timestep = 0.5  # fs

integrator = RPMDIntegrator(num_beads, temperature*unit.kelvin, 
                            friction/unit.picosecond, 
                            timestep*unit.femtoseconds)

# RPMD treats cavity particle classically (single bead)
context = openmm.Context(system, integrator)
context.setPositions(positions)

# Run RPMD simulation
integrator.step(100000)
```

## Command-Line Examples

### Water System

```bash
cd ml-experimental/cavity/examples/water_system

# Baseline (no cavity)
python run_simulation.py --no-cavity --molecules 1000 --prod 100

# With cavity coupling
python run_simulation.py --cavity-freq 1600 --lambda 0.01 --molecules 1000 --prod 100

# Different coupling strengths
python run_simulation.py --cavity-freq 1600 --lambda 0.001 --prod 100
python run_simulation.py --cavity-freq 1600 --lambda 0.010 --prod 100
python run_simulation.py --cavity-freq 1600 --lambda 0.050 --prod 100

# Analyze spectrum
python analyze_spectrum.py water_*.npz --xlim 1400 1800
```

### Dimer System (Toy Model)

```bash
cd ml-experimental/cavity/examples/dimer_system

# Run simulation
python run_simulation.py --lambda 0.001 --temp 100 --prod 900

# Analyze spectrum
python analyze_spectrum.py cavity_diamer_lambda*.npz

# Benchmark performance
python benchmark_speed.py
```

## Parameters

### Cavity Frequency (`omegac`)

Choose based on the molecular vibration you want to couple:

```python
# Water modes
OH_STRETCH = 3663.0  # cm^-1
HOH_BEND = 1600.0    # cm^-1

# Convert to atomic units
omegac_au = wavenumber_to_hartree(OH_STRETCH)
```

### Coupling Strength (`lambda`)

Dimensionless parameter controlling light-matter interaction strength:

- **λ < 0.01**: Weak coupling (perturbative)
- **0.01 ≤ λ ≤ 0.1**: Strong coupling (Rabi splitting visible)
- **λ > 0.1**: Ultra-strong coupling (requires careful numerical treatment)

### Photon Mass

Default is 1 a.u. (very light):

```python
photon_mass_amu = 1.0 / 1822.888  # 1 atomic unit
```

You can specify custom values:

```python
cavity_idx = add_cavity_particle(system, positions, omegac_au, 
                                photon_mass_amu=0.001)
```

## Advanced Usage

### Monitoring Cavity Displacement

```python
from ml_experimental.cavity.utils import get_cavity_displacement

# During simulation
state = context.getState(getPositions=True)
disp, mag = get_cavity_displacement(state, cavity_idx)
print(f"Cavity displacement: {mag:.4f} nm")
```

### Computing Dipole Moment

```python
from ml_experimental.cavity.utils import compute_dipole_moment

# Extract charges from system
charges = []
for force in system.getForces():
    if isinstance(force, openmm.NonbondedForce):
        for i in range(num_particles):
            charge, sigma, epsilon = force.getParticleParameters(i)
            charges.append(charge)
        break

# Compute dipole
state = context.getState(getPositions=True)
dipole = compute_dipole_moment(state, charges, num_particles)
dipole_mag = np.linalg.norm(dipole)
print(f"Dipole moment: {dipole_mag:.4f} e·nm")
```

### Thermostating

The cavity particle should NOT be thermostated (it's a quantum field mode):

```python
from ml_experimental.cavity.utils import setup_bussi_thermostat

# Only thermostat molecular particles
bussi = setup_bussi_thermostat(system, temperature_K=300, 
                               num_molecular_particles=n_atoms,
                               tau_ps=1.0)
```

### Saving Trajectories

```python
from openmm.app import PDBReporter

# Save trajectory (including cavity particle position)
reporter = PDBReporter('trajectory.pdb', 1000)
simulation.reporters.append(reporter)
```

## Common Patterns

### Pattern 1: Scan Coupling Strengths

```python
lambda_values = [0.001, 0.005, 0.01, 0.05, 0.1]

for lam in lambda_values:
    # Create system
    system = create_system()
    cavity_idx = add_cavity_particle(system, positions, omegac_au)
    setup_cavity_coupling(system, cavity_idx, omegac_au, lam)
    
    # Run simulation
    run_simulation(system, positions, output_file=f"lambda_{lam:.4f}.npz")
```

### Pattern 2: Scan Cavity Frequencies

```python
frequencies_cm = [1500, 1600, 1700, 1800]  # Around HOH bend

for freq in frequencies_cm:
    omegac_au = wavenumber_to_hartree(freq)
    
    # Create system
    system = create_system()
    cavity_idx = add_cavity_particle(system, positions, omegac_au)
    setup_cavity_coupling(system, cavity_idx, omegac_au, lambda_coupling)
    
    # Run simulation
    run_simulation(system, positions, output_file=f"freq_{freq}.npz")
```

### Pattern 3: Temperature Scan

```python
temperatures = [100, 200, 300, 400]  # K

for T in temperatures:
    # Create system with cavity
    system = create_system_with_cavity()
    
    # Thermostat at temperature T
    setup_bussi_thermostat(system, T, num_molecular_particles)
    
    # Run simulation
    run_simulation(system, positions, temperature=T, 
                  output_file=f"T_{T}K.npz")
```

## Troubleshooting

### Issue: No Rabi Splitting Observed

**Causes**:
1. Coupling strength too weak (increase λ)
2. Production time too short (need >100 ps for good spectrum)
3. Cavity frequency far from molecular resonance
4. Too few molecules (Rabi splitting scales as √N)

**Solutions**:
- Increase λ from 0.01 to 0.05
- Run longer production (1-10 ns)
- Match cavity frequency to molecular peak
- Use more molecules (>500 for water)

### Issue: Simulation Unstable

**Causes**:
1. Timestep too large
2. Coupling strength too large
3. Cavity particle not excluded from barostat

**Solutions**:
- Reduce timestep (use ≤0.5 fs)
- Reduce λ below 0.1
- Exclude cavity from barostat/thermostat

### Issue: Cavity Particle Drifts

**Expected behavior**: The cavity particle position represents the field amplitude and can be non-zero.

If magnitude is very large (>10 nm):
- Check that coupling is set correctly
- Verify photon mass is reasonable
- Check for energy drift

## Performance Tips

1. **Use GPU**: CUDA/OpenCL platforms are 10-50× faster
2. **Exclude cavity from cutoffs**: No nonbonded interactions needed
3. **Mixed precision**: Use mixed precision for GPU platforms
4. **Appropriate timestep**: 0.5 fs for flexible water, 1 fs for rigid

## Example Output

```
--- Adding Cavity Particle ---
  Cavity particle index: 3000
  Omega_c: 0.016693 a.u. = 3663.0 cm⁻¹
  Photon mass: 0.000549 amu

--- Setting Up Cavity Coupling ---
  Lambda: 0.01
  Coupling is ON from t=0
  CavityForce added with lambda=0.01
  CavityParticleDisplacer added with lambda=0.01

--- Setting Up Bussi Thermostat ---
  Temperature: 300.0 K
  Time constant: 1.0 ps
  Applied to 3000 molecular particles
  Cavity particle excluded from thermostat
```

## Next Steps

- See `CAVITY_THEORY.md` for theoretical background
- Check `examples/` for complete working examples
- Read `CAVITY_IMPLEMENTATION.md` for implementation details
