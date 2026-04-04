#!/usr/bin/env python3
"""
Benchmark OpenMM UMA against ASE + FAIRChem (reference implementation).

UMA has no native LAMMPS interface; ASE + FAIRChemCalculator is the canonical
reference. This script compares single-point energy/forces and short trajectory
between ASE and OpenMM using a single water molecule to avoid OOM.

Run: python tests/uma_ice_rpmd/benchmark_ase_reference.py
     python tests/uma_ice_rpmd/benchmark_ase_reference.py -o results/

Expects: OpenMM with RPMD plugin, openmm-ml, fairchem-core, ase
"""

import os
import sys
from pathlib import Path

if os.getenv('OPENMM_PLUGIN_DIR') is None:
    for _base in [os.getenv('CONDA_PREFIX'), os.path.expanduser('~/miniconda3')]:
        if _base:
            _plugins = os.path.join(_base, 'lib', 'plugins')
            if os.path.isdir(_plugins):
                os.environ['OPENMM_PLUGIN_DIR'] = _plugins
                break
_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
for _plugdir in (os.path.join(_root, 'build'), os.path.join(_root, 'build', 'lib', 'plugins')):
    if os.path.isdir(_plugdir) and any(f.startswith('libOpenMMRPMD') for f in os.listdir(_plugdir)):
        os.environ['OPENMM_PLUGIN_DIR'] = _plugdir
        break

import numpy as np
from ase import Atoms
from ase.build import molecule
from ase import units
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

# EV/angstrom to kJ/(mol*nm)
EV_ANG_TO_KJ_NM = 96.4853


def _single_water_in_box(cell_size=12.0):
    """Create one water molecule in a cubic box."""
    atoms = molecule('H2O')
    atoms.set_cell([cell_size, cell_size, cell_size])
    atoms.set_pbc(True)
    atoms.center()
    return atoms


def run_ase_reference(atoms, n_steps=50, dt_fs=0.5):
    """Run ASE MD with FAIRChem UMA and return energies, forces, final positions."""
    try:
        from fairchem.core import pretrained_mlip, FAIRChemCalculator
    except ImportError:
        from fairchem.core.calculate import pretrained_mlip
        from fairchem.core.calculate.ase_calculator import FAIRChemCalculator

    # Use CUDA when available for FAIRChem
    import torch
    _dev = 'cuda' if torch.cuda.is_available() else 'cpu'
    predictor = pretrained_mlip.get_predict_unit('uma-s-1', device=_dev)
    calc = FAIRChemCalculator(predictor, task_name='omol')
    atoms.calc = calc

    # Single-point at t=0
    e0 = atoms.get_potential_energy()
    f0 = atoms.get_forces()

    # MD
    MaxwellBoltzmannDistribution(atoms, temperature_K=300.0)
    dyn = VelocityVerlet(atoms, dt_fs * units.fs)
    for _ in range(n_steps):
        dyn.run(1)

    e_final = atoms.get_potential_energy()
    f_final = atoms.get_forces()
    pos_final = atoms.get_positions().copy()
    return {
        'energy_0': float(e0),
        'forces_0': np.array(f0),
        'energy_final': float(e_final),
        'forces_final': np.array(f_final),
        'positions_final': pos_final,
    }


def run_openmm_1bead(n_steps=50, dt_fs=0.5):
    """Run OpenMM NVT with 1 bead (classical) on same structure."""
    from openmm import app, unit, Vec3, Context, Platform
    from openmm import RPMDIntegrator
    from openmmml import MLPotential

    atoms_ase = _single_water_in_box()
    pos_ang = atoms_ase.get_positions()
    cell = atoms_ase.get_cell()

    topology = app.Topology()
    chain = topology.addChain()
    res = topology.addResidue('HOH', chain)
    topology.addAtom('O', app.Element.getBySymbol('O'), res)
    topology.addAtom('H', app.Element.getBySymbol('H'), res)
    topology.addAtom('H', app.Element.getBySymbol('H'), res)

    cell_nm = np.array(cell) * 0.1
    box_vectors = (
        Vec3(cell_nm[0, 0], cell_nm[0, 1], cell_nm[0, 2]) * unit.nanometer,
        Vec3(cell_nm[1, 0], cell_nm[1, 1], cell_nm[1, 2]) * unit.nanometer,
        Vec3(cell_nm[2, 0], cell_nm[2, 1], cell_nm[2, 2]) * unit.nanometer,
    )
    topology.setPeriodicBoxVectors(box_vectors)
    positions = (pos_ang * 0.1).tolist()
    pos_with_units = [Vec3(float(p[0]), float(p[1]), float(p[2])) * unit.nanometer for p in positions]

    potential = MLPotential('uma-s-1-pythonforce-batch')
    # Prefer CUDA platform for OpenMM
    try:
        _platform = Platform.getPlatformByName('CUDA')
        _ml_dev = 'cuda'
    except Exception:
        _platform = Platform.getPlatformByName('Reference')
        _ml_dev = 'cpu'
    system = potential.createSystem(
        topology,
        task_name='omol',
        charge=0,
        spin=1,
        device=_ml_dev,
        use_atom_wrap_for_lammps_parity=True,
    )

    integrator = RPMDIntegrator(
        1,
        300.0 * unit.kelvin,
        1.0 / unit.picosecond,
        dt_fs * unit.femtoseconds,
    )
    integrator.setThermostatType(RPMDIntegrator.Pile)
    context = Context(system, integrator, _platform)
    context.setPositions(pos_with_units)
    context.setPeriodicBoxVectors(*box_vectors)
    context.setVelocitiesToTemperature(300.0 * unit.kelvin)
    v = context.getState(getVelocities=True).getVelocities()
    integrator.setPositions(0, pos_with_units)
    integrator.setVelocities(0, v)

    state0 = integrator.getState(0, getEnergy=True, getForces=True)
    e0_kj = state0.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    f0_kj_nm = state0.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole / unit.nanometer)

    for _ in range(n_steps):
        integrator.step(1)

    state_final = integrator.getState(0, getEnergy=True, getForces=True, getPositions=True)
    e_final_kj = state_final.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    f_final_kj_nm = state_final.getForces(asNumpy=True).value_in_unit(
        unit.kilojoules_per_mole / unit.nanometer
    )
    pos_final_nm = state_final.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    pos_final_ang = pos_final_nm * 10.0

    return {
        'energy_0': e0_kj,
        'forces_0': f0_kj_nm,
        'energy_final': e_final_kj,
        'forces_final': f_final_kj_nm,
        'positions_final': pos_final_ang,
    }


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Benchmark OpenMM UMA vs ASE+FAIRChem')
    parser.add_argument('-o', '--output', type=str, default='.',
                        help='Output directory for metrics PNG (default: current directory)')
    args = parser.parse_args()
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    png_path = output_dir / 'benchmark_ase_reference_metrics.png'

    print("=" * 60)
    print("UMA OpenMM vs ASE+FAIRChem benchmark (single water molecule)")
    print("=" * 60)

    n_steps = 50
    dt_fs = 0.5

    print("\n[1] Running ASE + FAIRChem reference...")
    ase_result = run_ase_reference(_single_water_in_box(), n_steps=n_steps, dt_fs=dt_fs)

    print("[2] Running OpenMM 1-bead NVT...")
    omm_result = run_openmm_1bead(n_steps=n_steps, dt_fs=dt_fs)

    # Units: ASE uses eV and eV/angstrom
    e_ase_0 = ase_result['energy_0'] * 96.4853  # eV -> kJ/mol
    f_ase_0 = ase_result['forces_0'] * 96.4853 / 10.0  # eV/A -> kJ/(mol*nm)
    e_ase_final = ase_result['energy_final'] * 96.4853
    f_ase_final = ase_result['forces_final'] * 96.4853 / 10.0

    e_omm_0 = omm_result['energy_0']
    f_omm_0 = omm_result['forces_0']
    e_omm_final = omm_result['energy_final']
    f_omm_final = omm_result['forces_final']

    print("\n--- Single-point (t=0) ---")
    print(f"  ASE energy:   {e_ase_0:.4f} kJ/mol")
    print(f"  OpenMM energy: {e_omm_0:.4f} kJ/mol")
    print(f"  |ΔE| = {abs(e_ase_0 - e_omm_0):.4f} kJ/mol")

    print("\n--- After {} steps @ {} fs ---".format(n_steps, dt_fs))
    print(f"  ASE energy:   {e_ase_final:.4f} kJ/mol")
    print(f"  OpenMM energy: {e_omm_final:.4f} kJ/mol")

    max_f_diff_0 = np.max(np.abs(f_ase_0 - f_omm_0))
    print(f"\n  Max |ΔF| at t=0: {max_f_diff_0:.4f} kJ/(mol·nm)")

    # RMSD of final positions (wrap into same cell for comparison)
    from ase.geometry import wrap_positions
    cell = np.diag([12.0, 12.0, 12.0])
    p_ase = wrap_positions(ase_result['positions_final'], cell)
    p_omm = wrap_positions(omm_result['positions_final'], cell)
    rmsd = np.sqrt(np.mean(np.sum((p_ase - p_omm) ** 2, axis=1)))
    print(f"  Position RMSD (final): {rmsd:.4f} Å")

    # Stability check: O-H bonds
    o_ase, h1_ase, h2_ase = p_ase[0], p_ase[1], p_ase[2]
    oh1_ase = np.linalg.norm(o_ase - h1_ase)
    oh2_ase = np.linalg.norm(o_ase - h2_ase)
    o_omm, h1_omm, h2_omm = p_omm[0], p_omm[1], p_omm[2]
    oh1_omm = np.linalg.norm(o_omm - h1_omm)
    oh2_omm = np.linalg.norm(o_omm - h2_omm)
    print(f"\n  ASE O-H bonds: {oh1_ase:.3f}, {oh2_ase:.3f} Å")
    print(f"  OpenMM O-H bonds: {oh1_omm:.3f}, {oh2_omm:.3f} Å")

    ok = True
    if abs(e_ase_0 - e_omm_0) > 0.1:
        print("\n  WARNING: Energy at t=0 differs by > 0.1 kJ/mol")
        ok = False
    if max_f_diff_0 > 100.0:
        print("\n  WARNING: Max force diff at t=0 > 100 kJ/(mol·nm) - possible unit/sign error")
        ok = False
    if oh1_omm > 1.5 or oh2_omm > 1.5:
        print("\n  WARNING: OpenMM O-H bond too long (molecule unstable)")
        ok = False

    if ok:
        print("\n✓ Benchmark PASSED: ASE and OpenMM agree; single water stable.")
    else:
        print("\n✗ Benchmark issues detected.")
        sys.exit(1)

    # Save error metrics PNG
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    fig.suptitle('UMA OpenMM vs ASE+FAIRChem Error Metrics (single water)', fontsize=12)

    # 1. Energy comparison
    ax = axes[0, 0]
    labels = ['ASE t=0', 'OpenMM t=0', 'ASE final', 'OpenMM final']
    vals = [e_ase_0 / 1000, e_omm_0 / 1000, e_ase_final / 1000, e_omm_final / 1000]
    colors = ['#2ecc71', '#3498db', '#2ecc71', '#3498db']
    ax.bar(labels, vals, color=colors)
    ax.set_ylabel('Energy (×1000 kJ/mol)')
    ax.set_title('Energy comparison')
    ax.tick_params(axis='x', rotation=15)

    # 2. Max |ΔF| at t=0
    ax = axes[0, 1]
    ax.bar(['Max |ΔF|'], [max_f_diff_0], color='#e74c3c')
    ax.set_ylabel('kJ/(mol·nm)')
    ax.set_title('Force diff at t=0')

    # 3. Position RMSD (final)
    ax = axes[1, 0]
    ax.bar(['Position RMSD'], [rmsd], color='#9b59b6')
    ax.set_ylabel('Å')
    ax.set_title('Final position RMSD')

    # 4. O-H bond lengths
    ax = axes[1, 1]
    x = np.arange(4)
    oh_vals = [oh1_ase, oh2_ase, oh1_omm, oh2_omm]
    bars = ax.bar(['ASE O-H1', 'ASE O-H2', 'OpenMM O-H1', 'OpenMM O-H2'], oh_vals,
                  color=['#2ecc71', '#2ecc71', '#3498db', '#3498db'])
    ax.axhline(y=0.96, color='gray', linestyle='--', alpha=0.7, label='ref ~0.96 Å')
    ax.set_ylabel('Å')
    ax.set_title('O-H bond lengths')
    ax.tick_params(axis='x', rotation=15)

    plt.tight_layout()
    plt.savefig(png_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {png_path}")


if __name__ == '__main__':
    main()
