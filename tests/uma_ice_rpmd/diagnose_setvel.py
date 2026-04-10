#!/usr/bin/env python3
"""
diagnose_setvel.py: Deep-dive into setVelocitiesToTemperature KE after energy minimization.

Tests:
 A) Distribution of KE from setVelocitiesToTemperature (statistical vs systematic)
 B) RPMDIntegrator vs LangevinMiddleIntegrator: does RPMD context compute KE differently?
 C) Minimal water dimer test to confirm mass-less behavior
"""
import os
import sys

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
from openmm import app, unit, Vec3, RPMDIntegrator, LangevinMiddleIntegrator, Context, Platform, LocalEnergyMinimizer

INPUT_FILE = Path(__file__).parent / "pipeline_out" / "init_openmm_rpmd_64.xyz"
NUM_BEADS  = 8
TEMPERATURE = 243.0
DT_FS = 0.1
kT = 8.314e-3 * TEMPERATURE


def load_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].strip())
    symbols = []
    pos = []
    for line in lines[2:2 + n]:
        parts = line.split()
        symbols.append(parts[0])
        pos.append([float(x) * 0.1 for x in parts[1:4]])
    return symbols, np.array(pos)


def build_system_harmonic(symbols, positions_nm, box_nm):
    from openmm.app import element as elem
    _sym = {"O": elem.oxygen, "H": elem.hydrogen}
    topo = app.Topology()
    chain = topo.addChain()
    n_mol = len(symbols) // 3
    for m in range(n_mol):
        res = topo.addResidue("HOH", chain)
        topo.addAtom("O",  _sym["O"], res)
        topo.addAtom("H1", _sym["H"], res)
        topo.addAtom("H2", _sym["H"], res)
    a, b, c = box_nm
    topo.setPeriodicBoxVectors((
        Vec3(a, 0, 0) * unit.nanometer,
        Vec3(0, b, 0) * unit.nanometer,
        Vec3(0, 0, c) * unit.nanometer,
    ))
    masses_da = np.array([_sym[s].mass.value_in_unit(unit.dalton) for s in symbols])
    sys = openmm.System()
    sys.setDefaultPeriodicBoxVectors(
        Vec3(a,0,0)*unit.nanometer, Vec3(0,b,0)*unit.nanometer, Vec3(0,0,c)*unit.nanometer)
    for s in symbols:
        sys.addParticle(_sym[s].mass)
    harm = openmm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    harm.addGlobalParameter("k", 100.0)
    harm.addPerParticleParameter("x0")
    harm.addPerParticleParameter("y0")
    harm.addPerParticleParameter("z0")
    for i, pos in enumerate(positions_nm):
        harm.addParticle(i, list(pos))
    sys.addForce(harm)
    sys.addForce(openmm.CMMotionRemover())
    return sys, topo, masses_da


symbols, positions_nm = load_xyz(INPUT_FILE)
box_nm = np.array([0.899, 1.558, 1.464])
n_atoms = len(symbols)
n_mol = n_atoms // 3
ndof = 3 * n_atoms - 3
expected_ke = 0.5 * ndof * kT
print(f"System: {n_atoms} atoms, ndof={ndof}, expected_KE={expected_ke:.2f} kJ/mol")

system, topology, masses_da = build_system_harmonic(symbols, positions_nm, box_nm)

try:
    platform = Platform.getPlatformByName("CUDA")
    props = {"Precision": "mixed", "DeviceIndex": "0"}
except:
    platform = Platform.getPlatformByName("CPU")
    props = {}

pos_with_units = [Vec3(*p) * unit.nanometer for p in positions_nm]

# ─── Test A: Distribution of KE with RPMD integrator (multiple calls) ────────
print("\n=== Test A: KE distribution from setVelocitiesToTemperature (RPMD, after minimize) ===")
ke_samples = []
N_SAMPLES = 20

for trial in range(N_SAMPLES):
    intgr = RPMDIntegrator(NUM_BEADS, TEMPERATURE*unit.kelvin, 1.0/unit.picosecond, DT_FS*unit.femtoseconds)
    intgr.setRandomNumberSeed(100 + trial)  # different seed each trial
    intgr.setThermostatType(RPMDIntegrator.PileG)
    intgr.setCentroidFriction(0.5/unit.picosecond)
    ctx = Context(system, intgr, platform, props)
    ctx.setPositions(pos_with_units)
    ctx.setPeriodicBoxVectors(Vec3(box_nm[0],0,0)*unit.nanometer, Vec3(0,box_nm[1],0)*unit.nanometer, Vec3(0,0,box_nm[2])*unit.nanometer)
    for i in range(NUM_BEADS):
        intgr.setPositions(i, pos_with_units)
    # Minimize (same as actual sim path)
    state0 = intgr.getState(0, getPositions=True)
    ctx.setPositions(state0.getPositions())
    LocalEnergyMinimizer.minimize(ctx, tolerance=1.0, maxIterations=200)
    minimized_state = ctx.getState(getPositions=True)
    minimized_positions = minimized_state.getPositions()
    for i in range(NUM_BEADS):
        intgr.setPositions(i, minimized_positions)
    # Set velocities
    ctx.setPositions(minimized_positions)
    ctx.setVelocitiesToTemperature(TEMPERATURE * unit.kelvin)
    ke_ctx = ctx.getState(getEnergy=True).getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
    ke_samples.append(ke_ctx)
    del ctx, intgr

ke_arr = np.array(ke_samples)
print(f"  N={N_SAMPLES} trials: mean={ke_arr.mean():.2f} +/- {ke_arr.std():.2f} kJ/mol")
print(f"  Expected: {expected_ke:.2f} kJ/mol")
print(f"  Min={ke_arr.min():.2f}  Max={ke_arr.max():.2f}  Median={np.median(ke_arr):.2f}")
print(f"  Deviation from expected: {(ke_arr.mean()-expected_ke):.2f} kJ/mol = {(ke_arr.mean()-expected_ke)/ke_arr.std():.1f} sigma")
print(f"  All samples: {ke_arr.round(1)}")

# ─── Test B: RPMDIntegrator vs LangevinMiddleIntegrator ─────────────────────
print("\n=== Test B: RPMDIntegrator vs LangevinMiddleIntegrator KE ===")

# Same system, no minimization
for intgr_name, intgr_obj in [
    ("RPMD no-minimize",       RPMDIntegrator(NUM_BEADS, TEMPERATURE*unit.kelvin, 1.0/unit.picosecond, DT_FS*unit.femtoseconds)),
    ("Langevin no-minimize",   LangevinMiddleIntegrator(TEMPERATURE*unit.kelvin, 1.0/unit.picosecond, DT_FS*unit.femtoseconds)),
]:
    if hasattr(intgr_obj, 'setThermostatType'):
        intgr_obj.setThermostatType(RPMDIntegrator.PileG)
        intgr_obj.setCentroidFriction(0.5/unit.picosecond)
    intgr_obj.setRandomNumberSeed(42)
    ctx = Context(system, intgr_obj, platform, props)
    ctx.setPositions(pos_with_units)
    ctx.setPeriodicBoxVectors(Vec3(box_nm[0],0,0)*unit.nanometer, Vec3(0,box_nm[1],0)*unit.nanometer, Vec3(0,0,box_nm[2])*unit.nanometer)
    ctx.setVelocitiesToTemperature(TEMPERATURE * unit.kelvin)
    ke = ctx.getState(getEnergy=True).getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"  {intgr_name:30s}: KE = {ke:.2f} kJ/mol")
    del ctx, intgr_obj

# With minimization
for intgr_name, intgr_obj in [
    ("RPMD with-minimize",     RPMDIntegrator(NUM_BEADS, TEMPERATURE*unit.kelvin, 1.0/unit.picosecond, DT_FS*unit.femtoseconds)),
    ("Langevin with-minimize", LangevinMiddleIntegrator(TEMPERATURE*unit.kelvin, 1.0/unit.picosecond, DT_FS*unit.femtoseconds)),
]:
    if hasattr(intgr_obj, 'setThermostatType'):
        intgr_obj.setThermostatType(RPMDIntegrator.PileG)
        intgr_obj.setCentroidFriction(0.5/unit.picosecond)
    intgr_obj.setRandomNumberSeed(42)
    ctx = Context(system, intgr_obj, platform, props)
    ctx.setPositions(pos_with_units)
    ctx.setPeriodicBoxVectors(Vec3(box_nm[0],0,0)*unit.nanometer, Vec3(0,box_nm[1],0)*unit.nanometer, Vec3(0,0,box_nm[2])*unit.nanometer)
    if hasattr(intgr_obj, 'setPositions'):
        for i in range(NUM_BEADS):
            intgr_obj.setPositions(i, pos_with_units)
    # Minimize
    if hasattr(intgr_obj, 'getState'):
        s0 = intgr_obj.getState(0, getPositions=True)
        ctx.setPositions(s0.getPositions())
    LocalEnergyMinimizer.minimize(ctx, tolerance=1.0, maxIterations=200)
    min_state = ctx.getState(getPositions=True)
    min_pos = min_state.getPositions()
    if hasattr(intgr_obj, 'setPositions'):
        for i in range(NUM_BEADS):
            intgr_obj.setPositions(i, min_pos)
    ctx.setPositions(min_pos)
    ctx.setVelocitiesToTemperature(TEMPERATURE * unit.kelvin)
    ke = ctx.getState(getEnergy=True).getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"  {intgr_name:30s}: KE = {ke:.2f} kJ/mol")
    del ctx, intgr_obj

# ─── Test C: RPMD with minimize - check KE from context (no copyToContext) ───
print("\n=== Test C: Manual velocity check vs getState (context directly) ===")
intgr2 = RPMDIntegrator(NUM_BEADS, TEMPERATURE*unit.kelvin, 1.0/unit.picosecond, DT_FS*unit.femtoseconds)
intgr2.setRandomNumberSeed(77)
intgr2.setThermostatType(RPMDIntegrator.PileG)
intgr2.setCentroidFriction(0.5/unit.picosecond)
ctx2 = Context(system, intgr2, platform, props)
ctx2.setPositions(pos_with_units)
ctx2.setPeriodicBoxVectors(Vec3(box_nm[0],0,0)*unit.nanometer, Vec3(0,box_nm[1],0)*unit.nanometer, Vec3(0,0,box_nm[2])*unit.nanometer)
for i in range(NUM_BEADS):
    intgr2.setPositions(i, pos_with_units)

# Minimize
s0 = intgr2.getState(0, getPositions=True)
ctx2.setPositions(s0.getPositions())
LocalEnergyMinimizer.minimize(ctx2, tolerance=1.0, maxIterations=200)
min_s = ctx2.getState(getPositions=True)
min_pos2 = min_s.getPositions()
for i in range(NUM_BEADS):
    intgr2.setPositions(i, min_pos2)
ctx2.setPositions(min_pos2)

# Get velocities from context BEFORE setVelocitiesToTemperature
ke_pre = ctx2.getState(getEnergy=True).getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
print(f"  KE before setVelocitiesToTemperature: {ke_pre:.4f} kJ/mol")

# Set velocities
ctx2.setVelocitiesToTemperature(TEMPERATURE * unit.kelvin)

# Check KE from context (pure context, no RPMD copyToContext)
ke_ctx = ctx2.getState(getEnergy=True).getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
vel_from_ctx = ctx2.getState(getVelocities=True).getVelocities(asNumpy=True).value_in_unit(unit.nanometer/unit.picosecond)
manual_ke_ctx = float(0.5 * np.sum(masses_da[:, None] * vel_from_ctx**2))

print(f"  KE after setVelocitiesToTemperature (context.getState): {ke_ctx:.2f} kJ/mol")
print(f"  Manual KE from context velocities:                      {manual_ke_ctx:.2f} kJ/mol")
print(f"  Diff:                                                   {ke_ctx - manual_ke_ctx:+.4f} kJ/mol")
print(f"  Expected KE:                                            {expected_ke:.2f} kJ/mol")

# Now check what integrator.getState(b) returns for the SAME velocities
ke_intgr_b = intgr2.getState(0, getEnergy=True).getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
vel_from_intgr = intgr2.getState(0, getVelocities=True).getVelocities(asNumpy=True).value_in_unit(unit.nanometer/unit.picosecond)
manual_ke_intgr = float(0.5 * np.sum(masses_da[:, None] * vel_from_intgr**2))
print(f"\n  After one integrator.getState(0) [calls copyToContext]:")
print(f"  KE from integrator.getState(0):             {ke_intgr_b:.2f} kJ/mol")
print(f"  Manual KE (intgr.getState vels):            {manual_ke_intgr:.2f} kJ/mol")
print(f"  Vel RMS diff (ctx vs intgr):                {np.sqrt(np.mean((vel_from_ctx - vel_from_intgr)**2)):.6f} nm/ps")

# Check whether context KE changed AFTER integrator.getState
ke_ctx_after = ctx2.getState(getEnergy=True).getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
vel_ctx_after = ctx2.getState(getVelocities=True).getVelocities(asNumpy=True).value_in_unit(unit.nanometer/unit.picosecond)
print(f"  Context KE AFTER integrator.getState(0):    {ke_ctx_after:.2f} kJ/mol  (changed? was {ke_ctx:.2f})")

print(f"\n  Context vel changed by copyToContext: {np.sqrt(np.mean((vel_from_ctx - vel_ctx_after)**2)):.6f} nm/ps RMS")

print("\nDone.")
