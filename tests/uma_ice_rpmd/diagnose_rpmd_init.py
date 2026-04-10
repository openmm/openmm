#!/usr/bin/env python3
"""
diagnose_rpmd_init.py: Diagnose the RPMDIntegrator temperature explosion.

Step 1: Check particle masses and whether integrator.getState(b).getKineticEnergy()
        matches a manual computation from the Python-side set velocities.
        Uses a harmonic force (no UMA) to keep runtime under 10 seconds.

Step 2 (optional): NVE drift check with the harmonic force.
"""

import sys
import os

# Set plugin dir before openmm import
if os.getenv('OPENMM_PLUGIN_DIR') is None:
    for _base in [os.getenv('CONDA_PREFIX'), os.path.expanduser('~/miniconda3')]:
        if _base:
            _plugins = os.path.join(_base, 'lib', 'plugins')
            if os.path.isdir(_plugins):
                os.environ['OPENMM_PLUGIN_DIR'] = _plugins
                break

import numpy as np
from pathlib import Path
import openmm
from openmm import app, unit, Vec3, RPMDIntegrator, Context, Platform

# ─── Reproduce the exact ice system geometry without UMA ─────────────────────
# Use the same 64-molecule ice structure that the actual run uses
INPUT_FILE = Path(__file__).parent / "pipeline_out" / "init_openmm_rpmd_64.xyz"

NUM_BEADS   = 8
TEMPERATURE = 243.0  # K
DT_FS       = 0.1    # fs - same as the failing run
SEED        = 284759

kT = 8.314e-3 * TEMPERATURE  # kJ/mol


def load_xyz(path):
    """Read an XYZ file, return (symbols, positions_nm, comment)."""
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].strip())
    comment = lines[1]
    symbols, pos = [], []
    for line in lines[2:2 + n]:
        parts = line.split()
        symbols.append(parts[0])
        pos.append([float(x) * 0.1 for x in parts[1:4]])  # Å → nm
    return symbols, np.array(pos), comment


def build_system_no_ml(symbols, positions_nm, box_nm):
    """Build a minimal OpenMM system with harmonic bond forces (no ML). """
    from openmm.app import element as elem
    _sym_map = {"O": elem.oxygen, "H": elem.hydrogen}

    topology = app.Topology()
    chain = topology.addChain()
    n_mol = len(symbols) // 3
    for m in range(n_mol):
        res = topology.addResidue("HOH", chain)
        topology.addAtom("O", _sym_map["O"], res)
        topology.addAtom("H1", _sym_map["H"], res)
        topology.addAtom("H2", _sym_map["H"], res)

    # Orthorhombic box (approximate)
    a, b, c = box_nm
    topology.setPeriodicBoxVectors((
        Vec3(a, 0, 0) * unit.nanometer,
        Vec3(0, b, 0) * unit.nanometer,
        Vec3(0, 0, c) * unit.nanometer,
    ))

    system = openmm.System()
    system.setDefaultPeriodicBoxVectors(
        Vec3(a, 0, 0) * unit.nanometer,
        Vec3(0, b, 0) * unit.nanometer,
        Vec3(0, 0, c) * unit.nanometer,
    )

    # Add particles with correct masses
    masses_da = []
    for sym in symbols:
        m = _sym_map[sym].mass.value_in_unit(unit.dalton)
        system.addParticle(m)
        masses_da.append(m)

    # Simple harmonic force to keep atoms near initial positions (prevents drift)
    # This is just to allow NVE drift check - it doesn't affect mass/velocity tests
    harmonic = openmm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    harmonic.addGlobalParameter("k", 100.0)  # kJ/(mol nm^2)
    harmonic.addPerParticleParameter("x0")
    harmonic.addPerParticleParameter("y0")
    harmonic.addPerParticleParameter("z0")
    for i, pos in enumerate(positions_nm):
        harmonic.addParticle(i, [pos[0], pos[1], pos[2]])
    system.addForce(harmonic)

    # CMMotionRemover (same as MLPotential.createSystem default)
    system.addForce(openmm.CMMotionRemover())

    return system, topology, np.array(masses_da)


def manual_ke(velocities_nm_ps, masses_da):
    """KE = 0.5 * sum(m * v^2), result in kJ/mol."""
    return 0.5 * float(np.sum(masses_da[:, None] * velocities_nm_ps**2))


# ─── Main diagnostic ──────────────────────────────────────────────────────────
print("=" * 70)
print("RPMD Temperature Explosion Diagnostic")
print("=" * 70)

# Load structure
if not INPUT_FILE.exists():
    sys.exit(f"Input file not found: {INPUT_FILE}")

symbols, positions_nm, comment = load_xyz(INPUT_FILE)
n_atoms = len(symbols)
n_O = sum(1 for s in symbols if s == "O")
n_H = sum(1 for s in symbols if s == "H")
print(f"\nSystem: {n_atoms} atoms ({n_O} O, {n_H} H), {n_atoms//3} molecules")

# Parse box from XYZ comment line (GenIce format: "Lattice=...")
# Fallback: use approximate values from the simulation output
box_nm = np.array([0.899, 1.558, 1.464])  # from terminal output (Å→nm: 8.99/10, ...)
try:
    for token in comment.split():
        if token.startswith("Lattice="):
            vals = [float(x) for x in token[9:].replace('"', '').split()]
            # vals is 3x3 flattened: a_x, a_y, a_z, b_x, ...
            if len(vals) == 9:
                box_nm = np.array([vals[0], vals[4], vals[8]]) * 0.1  # Å→nm
    print(f"Box: {box_nm * 10} Å (approx)")
except Exception as e:
    print(f"Box parse failed ({e}), using fallback: {box_nm * 10} Å")

# Build system
system, topology, masses_da = build_system_no_ml(symbols, positions_nm, box_nm)
n_massive = int(np.sum(masses_da > 0))
ndof_context = 3 * n_massive - 3  # CMMotionRemover removes 3 DOF
expected_ke = 0.5 * ndof_context * kT

print(f"\n--- Particle masses (first 12 atoms) ---")
for i in range(min(12, n_atoms)):
    m = system.getParticleMass(i).value_in_unit(unit.dalton)
    print(f"  Atom {i:3d} ({symbols[i]}): mass = {m:.4f} Da")
print(f"  ...")
print(f"  Total massive atoms: {n_massive}, ndof (context) = {ndof_context}")
print(f"  Expected KE at {TEMPERATURE} K: {expected_ke:.2f} kJ/mol")
print(f"  kT = {kT:.4f} kJ/mol")

# Build RPMD integrator on CPU for reproducibility
platform_name = "CUDA"
try:
    platform = Platform.getPlatformByName(platform_name)
    properties = {"Precision": "mixed", "DeviceIndex": "0"}
    print(f"\nUsing {platform_name} platform (mixed precision)")
except Exception:
    platform_name = "CPU"
    platform = Platform.getPlatformByName("CPU")
    properties = {}
    print(f"\n{platform_name} platform")

integrator = RPMDIntegrator(
    NUM_BEADS,
    TEMPERATURE * unit.kelvin,
    1.0 / unit.picosecond,   # internal friction
    DT_FS * unit.femtoseconds,
)
integrator.setRandomNumberSeed(SEED)
integrator.setThermostatType(RPMDIntegrator.PileG)
integrator.setCentroidFriction(0.5 / unit.picosecond)

context = Context(system, integrator, platform, properties)

# Set positions and box
pos_with_units = [Vec3(*p) * unit.nanometer for p in positions_nm]
context.setPositions(pos_with_units)
context.setPeriodicBoxVectors(
    Vec3(box_nm[0], 0, 0) * unit.nanometer,
    Vec3(0, box_nm[1], 0) * unit.nanometer,
    Vec3(0, 0, box_nm[2]) * unit.nanometer,
)

for i in range(NUM_BEADS):
    integrator.setPositions(i, pos_with_units)

# ─── STEP 0: Energy minimization (mirror the actual simulation path) ─────────
print("\n--- Step 0: Energy Minimization (mirror actual sim) ---")
# Set bead 0 positions in context for minimization
context.setPositions(pos_with_units)
context.setPeriodicBoxVectors(
    Vec3(box_nm[0], 0, 0) * unit.nanometer,
    Vec3(0, box_nm[1], 0) * unit.nanometer,
    Vec3(0, 0, box_nm[2]) * unit.nanometer,
)
state_pre = context.getState(getEnergy=True)
pe_pre = state_pre.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
print(f"  PE before minimization: {pe_pre:.2f} kJ/mol")
# Check context KE BEFORE minimization
ctx_ke_before_min = context.getState(getEnergy=True).getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
print(f"  KE before minimization (no vel set yet): {ctx_ke_before_min:.2f} kJ/mol")

from openmm import LocalEnergyMinimizer
LocalEnergyMinimizer.minimize(context, tolerance=1.0, maxIterations=200)
minimized_state = context.getState(getPositions=True, getEnergy=True)
minimized_positions = minimized_state.getPositions()
pe_post = minimized_state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
print(f"  PE after  minimization: {pe_post:.2f} kJ/mol")

# Copy minimized positions to all beads
for i in range(NUM_BEADS):
    integrator.setPositions(i, minimized_positions)

# ─── STEP 1A: Check context KE directly after setVelocitiesToTemperature ─────
print("\n--- Step 1A: Context KE after setVelocitiesToTemperature ---")
np.random.seed(SEED)
context.setPositions(minimized_positions)
context.setVelocitiesToTemperature(TEMPERATURE * unit.kelvin)
ctx_state = context.getState(getVelocities=True, getEnergy=True)
context_ke_direct = ctx_state.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
base_velocities = ctx_state.getVelocities()
vel_array = np.array([[v.x, v.y, v.z] for v in base_velocities])  # nm/ps

manual_ke_base = manual_ke(vel_array, masses_da)
print(f"  Context KE (direct):    {context_ke_direct:.2f} kJ/mol  ← should match expected")
print(f"  Manual KE (Python vel): {manual_ke_base:.2f} kJ/mol")
print(f"  Expected KE:            {expected_ke:.2f} kJ/mol")

# ─── STEP 1B: Set velocities to integrator beads (same code as actual sim) ───
print("\n--- Step 1B: Set velocities to integrator (actual sim code path) ---")
bead_vels_set = []  # Track what we actually set
for i in range(NUM_BEADS):
    velocities = []
    for v in base_velocities:
        vx = v[0].value_in_unit(unit.nanometer / unit.picosecond) + 0.01 * np.random.randn()
        vy = v[1].value_in_unit(unit.nanometer / unit.picosecond) + 0.01 * np.random.randn()
        vz = v[2].value_in_unit(unit.nanometer / unit.picosecond) + 0.01 * np.random.randn()
        velocities.append(Vec3(vx, vy, vz) * unit.nanometer / unit.picosecond)
    integrator.setVelocities(i, velocities)
    vel_set = np.array([[v[0].value_in_unit(unit.nanometer / unit.picosecond),
                         v[1].value_in_unit(unit.nanometer / unit.picosecond),
                         v[2].value_in_unit(unit.nanometer / unit.picosecond)]
                        for v in velocities])
    bead_vels_set.append(vel_set)

# Now check what integrator reports
print(f"  {'Bead':>4s}  {'getState KE':>14s}  {'Manual KE (set)':>16s}  {'Diff':>8s}")
print(f"  {'-'*4}  {'-'*14}  {'-'*16}  {'-'*8}")
for i in range(NUM_BEADS):
    state = integrator.getState(i, getEnergy=True)
    ke_from_state = state.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
    ke_manual = manual_ke(bead_vels_set[i], masses_da)
    diff = ke_from_state - ke_manual
    flag = "  *** MISMATCH" if abs(diff) > 5.0 else ""
    print(f"  {i:4d}  {ke_from_state:14.2f}  {ke_manual:16.2f}  {diff:+8.2f}{flag}")

# ─── STEP 1C: Check velocities returned by getState vs what we set ───────────
print("\n--- Step 1C: Velocities returned by getState vs set ---")
bead = 0
state_vels = integrator.getState(bead, getVelocities=True).getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond)
vel_set_b0 = bead_vels_set[bead]
v_diff = state_vels - vel_set_b0
rms_diff = np.sqrt(np.mean(v_diff**2))
max_diff = np.abs(v_diff).max()
print(f"  Bead 0 velocity RMS diff (getState - set): {rms_diff:.6f} nm/ps")
print(f"  Bead 0 velocity max diff (getState - set): {max_diff:.6f} nm/ps")
if rms_diff > 0.01:
    print("  *** VELOCITY SCRAMBLING DETECTED: getState returns different velocities than set!")

# Show first few atoms
print(f"\n  Atom  sym   set_vx    set_vy    set_vz    got_vx    got_vy    got_vz")
for i in range(min(9, n_atoms)):
    print(f"  {i:4d}  {symbols[i]:3s}  "
          f"{vel_set_b0[i,0]:8.4f}  {vel_set_b0[i,1]:8.4f}  {vel_set_b0[i,2]:8.4f}  "
          f"{state_vels[i,0]:8.4f}  {state_vels[i,1]:8.4f}  {state_vels[i,2]:8.4f}")

# ─── STEP 1D: Verify masses seen by KE computation ───────────────────────────
print("\n--- Step 1D: Effective mass inferred from KE and velocity ---")
# KE = 0.5 * sum(m_eff[i] * v[i]^2). We can compute per-atom KE contribution
# by getting velocities from getState and comparing with manual masses.
# KE contribution per atom from getState vs from correct masses
per_atom_v2 = np.sum(state_vels**2, axis=1)  # v^2 per atom (nm/ps)^2
ke_total_fromstate = integrator.getState(0, getEnergy=True).getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
ke_manual_fromgot = manual_ke(state_vels, masses_da)
print(f"  KE from getState:           {ke_total_fromstate:.2f} kJ/mol")
print(f"  Manual KE (state vels, correct masses): {ke_manual_fromgot:.2f} kJ/mol")
print(f"  Expected:                   {expected_ke:.2f} kJ/mol")
print(f"  Diff (getState - manual):   {ke_total_fromstate - ke_manual_fromgot:+.2f} kJ/mol")

# ─── STEP 2: NVE energy conservation (no thermostat) ─────────────────────────
print("\n" + "=" * 70)
print("STEP 2: NVE Energy Conservation (100 steps, thermostat off)")
print("=" * 70)

# Re-setup with fresh context (NVE)
integrator_nve = RPMDIntegrator(
    NUM_BEADS,
    TEMPERATURE * unit.kelvin,
    1.0 / unit.picosecond,
    DT_FS * unit.femtoseconds,
)
integrator_nve.setRandomNumberSeed(SEED)
integrator_nve.setThermostatType(RPMDIntegrator.NoneThermo)

context_nve = Context(system, integrator_nve, platform, properties)
context_nve.setPositions(pos_with_units)
context_nve.setPeriodicBoxVectors(
    Vec3(box_nm[0], 0, 0) * unit.nanometer,
    Vec3(0, box_nm[1], 0) * unit.nanometer,
    Vec3(0, 0, box_nm[2]) * unit.nanometer,
)
for i in range(NUM_BEADS):
    integrator_nve.setPositions(i, pos_with_units)

np.random.seed(SEED)
context_nve.setVelocitiesToTemperature(TEMPERATURE * unit.kelvin)
base_vels = context_nve.getState(getVelocities=True).getVelocities()
for i in range(NUM_BEADS):
    velocities = []
    for v in base_vels:
        vx = v[0].value_in_unit(unit.nanometer / unit.picosecond) + 0.01 * np.random.randn()
        vy = v[1].value_in_unit(unit.nanometer / unit.picosecond) + 0.01 * np.random.randn()
        vz = v[2].value_in_unit(unit.nanometer / unit.picosecond) + 0.01 * np.random.randn()
        velocities.append(Vec3(vx, vy, vz) * unit.nanometer / unit.picosecond)
    integrator_nve.setVelocities(i, velocities)


def sum_bead_total_energy(intgr, n_beads):
    total = 0.0
    for b in range(n_beads):
        s = intgr.getState(b, getEnergy=True)
        total += s.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        total += s.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
    return total


E0 = sum_bead_total_energy(integrator_nve, NUM_BEADS)
print(f"  Initial sum_bead(PE+KE): {E0:.4f} kJ/mol")

N_NVE = 100
print(f"  Running {N_NVE} NVE steps (dt={DT_FS} fs, no thermostat)...")
integrator_nve.step(N_NVE)

E1 = sum_bead_total_energy(integrator_nve, NUM_BEADS)
dE = E1 - E0
dt_total_ps = N_NVE * DT_FS * 1e-3
n_mol = n_atoms // 3
drift_rate = abs(dE) / n_mol / dt_total_ps
print(f"  Final   sum_bead(PE+KE): {E1:.4f} kJ/mol")
print(f"  ΔE = {dE:+.4f} kJ/mol over {dt_total_ps:.4f} ps")
print(f"  |ΔE|/(N_mol·Δt) = {drift_rate:.1f} kJ/(mol·ps)")
if drift_rate > 100:
    print("  *** HIGH ENERGY DRIFT: likely non-conservative forces or bad integration")
else:
    print("  Energy drift is acceptable")

# Temperature after NVE steps
print(f"\n--- Bead energies after NVE steps ---")
bead_ke_nve = []
for b in range(NUM_BEADS):
    s = integrator_nve.getState(b, getEnergy=True)
    ke = s.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
    bead_ke_nve.append(ke)
mean_ke_nve = np.mean(bead_ke_nve)
T_eff = (2.0 * mean_ke_nve) / (ndof_context * 8.314e-3)
print(f"  Mean bead KE: {mean_ke_nve:.2f} kJ/mol  → T_eff = {T_eff:.1f} K (using ndof={ndof_context})")

print("\n" + "=" * 70)
print("DIAGNOSTIC COMPLETE")
print("=" * 70)
