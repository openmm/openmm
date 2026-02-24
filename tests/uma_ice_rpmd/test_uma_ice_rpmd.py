#!/usr/bin/env python3
"""
RPMD simulation testing if ice stays frozen at 243K with UMA or eSEN potentials.
Uses proper ice Ih structure (GenIce2 or embedded CIF). PILE thermostat, NPT/NVT.
Run: python test_uma_ice_rpmd.py --molecules 32 --beads 8 --temperature 243 --dt 1.0 --equil 5 --prod 50
"""

import sys
import os
import subprocess
import tempfile

# Set plugin dir before openmm import. Try CONDA_PREFIX, then miniconda3, then build.
if os.getenv('OPENMM_PLUGIN_DIR') is None:
    for _base in [os.getenv('CONDA_PREFIX'), os.path.expanduser('~/miniconda3')]:
        if _base:
            _plugins = os.path.join(_base, 'lib', 'plugins')
            if os.path.isdir(_plugins):
                os.environ['OPENMM_PLUGIN_DIR'] = _plugins
                break

import numpy as np
import time
from pathlib import Path

from openmm import app, unit, Vec3
from openmm import RPMDIntegrator, Context, Platform, RPMDMonteCarloBarostat
from openmmml import MLPotential

# Embedded ice Ih CIF from Avogadro/COD (12 water molecules) - fallback when GenIce2 unavailable
# P 63 c m, a=b=7.82 A, c=7.36 A, gamma=120
_ICE_IH_CIF = """
data_1011023
_cell_length_a 7.82
_cell_length_b 7.82
_cell_length_c 7.36
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 120
_symmetry_space_group_name_H-M 'P 63 c m'
loop_
_symmetry_equiv_pos_as_xyz
 x,y,z
 -y,x-y,z
 y-x,-x,z
 y,x,z
 x-y,-y,z
 -x,y-x,z
 -x,-y,1/2+z
 y,y-x,1/2+z
 x-y,x,1/2+z
 -y,-x,1/2+z
 y-x,y,1/2+z
 x,x-y,1/2+z
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
 O1 O 0.3333 0. 0.0625
 O2 O 0.6667 0. 0.9375
 H1 H 0.3333 0. 0.174
 H2 H 0.438 0. 0.026
 H3 H 0.772 0.105 0.975
"""


def _create_ice_from_ase(atoms_ase, num_molecules):
    """Convert ASE Atoms (ice) to OpenMM topology and positions. Trim to num_molecules."""
    # Ensure we have O and H in correct order (O, H, H per molecule)
    symbols = atoms_ase.get_chemical_symbols()
    pos = atoms_ase.get_positions()
    cell = atoms_ase.get_cell()
    mol_list = []
    seen = set()
    for i, s in enumerate(symbols):
        if s == 'O' and i not in seen:
            oidx = i
            # Find closest 2 H to this O
            dists = []
            for j, (sj, pj) in enumerate(zip(symbols, pos)):
                if sj == 'H':
                    d = np.linalg.norm(pos[oidx] - pj)
                    dists.append((d, j))
            dists.sort(key=lambda x: x[0])
            h1, h2 = dists[0][1], dists[1][1]
            mol_list.append([oidx, h1, h2])
            seen.update([oidx, h1, h2])

    mol_count_full = len(mol_list)
    mol_list = mol_list[:num_molecules]

    # When we trim molecules, scale the box to preserve ice Ih density (~0.92 g/cm³)
    # Box volume should scale with number of molecules: V_new = V_old * (n_new / n_full)
    scale = (num_molecules / mol_count_full) ** (1.0 / 3.0) if num_molecules < mol_count_full else 1.0

    topology = app.Topology()
    positions = []
    for oidx, h1, h2 in mol_list:
        chain = topology.addChain()
        residue = topology.addResidue('HOH', chain)
        o_atom = topology.addAtom('O', app.Element.getBySymbol('O'), residue)
        h1_atom = topology.addAtom('H', app.Element.getBySymbol('H'), residue)
        h2_atom = topology.addAtom('H', app.Element.getBySymbol('H'), residue)
        topology.addBond(o_atom, h1_atom)
        topology.addBond(o_atom, h2_atom)
        for idx in [oidx, h1, h2]:
            positions.append((pos[idx] * scale).tolist())  # Scale positions with box

    # Box in nm (ASE cell: rows are lattice vectors); scale when trimmed
    v = np.array(cell) * scale * 0.1  # 0.1: Angstrom to nm
    box_vectors = (
        [v[0, 0], v[0, 1], v[0, 2]] * unit.nanometer,
        [v[1, 0], v[1, 1], v[1, 2]] * unit.nanometer,
        [v[2, 0], v[2, 1], v[2, 2]] * unit.nanometer,
    )

    topology.setPeriodicBoxVectors(box_vectors)
    return topology, positions, box_vectors


def _create_ice_genice(num_molecules):
    """Try GenIce2 subprocess. Returns (topology, positions, box_vectors) or None on failure."""
    try:
        # Replication to get >= num_molecules: 12 * r^3 >= n  =>  r = ceil((n/12)^(1/3))
        r = max(1, int(np.ceil((num_molecules / 12) ** (1 / 3))))
        with tempfile.NamedTemporaryFile(suffix='.y', delete=False) as f:
            outpath = f.name
        result = subprocess.run(
            ['genice2', '1h', '--rep', str(r), str(r), str(r), '-f', 'y', '-o', outpath, '--seed', '42'],
            capture_output=True,
            text=True,
            timeout=30,
        )
        if result.returncode != 0:
            return None
        from ase.io import read
        atoms = read(outpath, format='yaml')
        os.unlink(outpath)
        return _create_ice_from_ase(atoms, num_molecules)
    except Exception:
        return None


def _create_ice_embedded(num_molecules):
    """Create ice Ih from embedded CIF (12 molecules) + replication."""
    from ase.io import read
    from io import StringIO

    atoms = read(StringIO(_ICE_IH_CIF), format='cif')
    n_base = len(atoms) // 3
    if num_molecules <= n_base:
        return _create_ice_from_ase(atoms, num_molecules)
    r = max(1, int(np.ceil((num_molecules / n_base) ** (1 / 3))))
    from ase.build import make_supercell
    P = np.diag([r, r, r])
    supercell = make_supercell(atoms, P)
    return _create_ice_from_ase(supercell, num_molecules)


def create_ice_structure(num_molecules=32):
    """
    Create ice Ih structure with tetrahedral hydrogen-bonding network.

    Tries GenIce2 first (proton-disordered); falls back to embedded CIF (ice Ih from COD).
    Ice Ih: O-O ~2.75 A, density ~0.92 g/cm3.

    Returns
    -------
    topology : app.Topology
    positions : list of Vec3
    box_vectors : tuple of Vec3
    """
    print(f"\n--- Creating Ice Ih Structure ---")
    print(f"  Target molecules: {num_molecules}")

    result = _create_ice_genice(num_molecules)
    if result is not None:
        topology, positions, box_vectors = result
        print(f"  Source: GenIce2 (proton-disordered)")
    else:
        topology, positions, box_vectors = _create_ice_embedded(num_molecules)
        print(f"  Source: Embedded CIF (ice Ih from Avogadro/COD)")

    n_atoms = len(positions)
    mol_count = n_atoms // 3
    print(f"  Actual molecules: {mol_count}")
    print(f"  Total atoms: {n_atoms}")

    # Use determinant for correct volume (hexagonal/triclinic cells)
    box_matrix = np.array([[box_vectors[i][j].value_in_unit(unit.nanometer) for j in range(3)] for i in range(3)])
    vol_nm3 = np.abs(np.linalg.det(box_matrix))
    total_mass_g = (mol_count * 18.015) * 1.66054e-24
    density = total_mass_g / (vol_nm3 * 1e-21)  # 1 nm³ = 1e-21 cm³
    print(f"  Initial density: {density:.3f} g/cm³")
    print(f"  Expected ice Ih: ~0.92 g/cm³")

    return topology, positions, box_vectors


print("All required packages loaded successfully")


def run_simulation(num_molecules=32, num_beads=8, temperature_K=243.0,
                   pressure_bar=1.0, dt_fs=1.0, equilibration_ps=5.0,
                   production_ps=100.0, model_name='uma-s-1-pythonforce-batch',
                   output_dir='.', report_interval_ps=1.0, pdb_interval_ps=1.0,
                   platform_name='cuda', ml_device=None, precision_override=None):
    """
    Run RPMD simulation of ice.
    
    Parameters
    ----------
    num_molecules : int
        Number of water molecules
    num_beads : int
        Number of RPMD beads
    temperature_K : float
        Temperature in Kelvin
    pressure_bar : float
        Pressure in bar (0 = NVT, >0 = NPT)
    dt_fs : float
        Timestep in femtoseconds
    equilibration_ps : float
        Equilibration time in picoseconds
    production_ps : float
        Production time in picoseconds
    model_name : str
        UMA or eSEN model name
    output_dir : str
        Output directory for trajectories
    report_interval_ps : float
        Console reporting interval in picoseconds (default: 1.0)
    pdb_interval_ps : float
        PDB frame save interval in picoseconds (default: 1.0)
    platform_name : str
        OpenMM platform: 'cuda' (GPU) or 'cpu' (default: 'cuda')
    ml_device : str or None
        ML model device: 'cpu', 'cuda', or None. When None with CUDA platform,
        defaults to 'cpu' to avoid PyTorch/OpenMM CUDA context conflict.
        
    Returns
    -------
    success : bool
        True if simulation completed successfully
    """
    
    print("=" * 80)
    print("UMA/eSEN Ice RPMD Simulation")
    print("=" * 80)
    print(f"Model: {model_name}")
    print(f"Temperature: {temperature_K} K")
    print(f"Pressure: {pressure_bar} bar {'(NPT)' if pressure_bar > 0 else '(NVT)'}")
    print(f"RPMD beads: {num_beads}")
    print(f"Timestep: {dt_fs} fs")
    print(f"Molecules: {num_molecules}")
    print(f"Equilibration: {equilibration_ps} ps")
    print(f"Production: {production_ps} ps")
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Create structure
    topology, positions, box_vectors = create_ice_structure(num_molecules)
    n_atoms = len(positions)
    
    # Create potential
    print(f"\n--- Creating ML Potential ---")
    print("Creating OpenMM system...")
    potential = MLPotential(model_name)
    
    if ml_device is not None:
        _ml_device = ml_device if ml_device != 'auto' else None
    elif platform_name.lower() == 'cuda':
        _ml_device = 'cuda'  # Lazy-load allows single GPU; use --ml-device cpu to force CPU
    else:
        _ml_device = None

    print(f"  ML model device: {_ml_device or 'auto (library default)'}")

    system = potential.createSystem(
        topology,
        task_name='omol',
        charge=0,
        spin=1,
        device=_ml_device
    )
    
    print(f"  System uses PBC: {system.usesPeriodicBoundaryConditions()}")
    
    # Add barostat if NPT
    if pressure_bar > 0:
        barostat = RPMDMonteCarloBarostat(
            pressure_bar * unit.bar,
            25  # Attempt every 25 steps
        )
        system.addForce(barostat)
        print(f"  Added RPMDMonteCarloBarostat (NPT ensemble)")
        print(f"  Pressure: {pressure_bar} bar")
    else:
        print(f"  Running NVT (no barostat)")
    
    # Calculate density (use determinant for hexagonal/triclinic boxes)
    mass_per_molecule = 18.015
    total_mass_g = (num_molecules * mass_per_molecule) * 1.66054e-24
    box_matrix = np.array([[box_vectors[i][j].value_in_unit(unit.nanometer) for j in range(3)] for i in range(3)])
    vol_nm3 = np.abs(np.linalg.det(box_matrix))
    density = total_mass_g / (vol_nm3 * 1e-21)
    print(f"  Initial density: {density:.3f} g/cm³")
    print(f"  Expected ice density: ~0.92 g/cm³")
    
    # Create RPMD integrator with PILE_G thermostat
    print(f"\n--- Creating RPMD Integrator ---")
    print(f"  Beads: {num_beads}")
    print(f"  Temperature: {temperature_K} K")
    print(f"  Timestep: {dt_fs} fs")
    
    friction = 1.0  # ps^-1
    centroid_friction = 0.5  # ps^-1 for Bussi thermostat
    
    integrator = RPMDIntegrator(
        num_beads,
        temperature_K * unit.kelvin,
        friction / unit.picosecond,
        dt_fs * unit.femtoseconds
    )
    
    # Set PILE thermostat (standard Langevin for all modes)
    integrator.setThermostatType(RPMDIntegrator.Pile)
    
    print(f"  Thermostat: PILE (Langevin for all modes)")
    print(f"  Friction: {friction} ps^-1")
    
    # Create context
    print(f"\n--- Creating Simulation Context ---")
    # CRITICAL: Initialize PyTorch CUDA before OpenMM creates its Context.
    # torch.cuda.init() pops existing contexts from the stack; if OpenMM creates
    # first, PyTorch later corrupts the context (CUDA_ERROR_ILLEGAL_ADDRESS).
    # Init PyTorch first so both share the primary context (pytorch#75025).
    if platform_name.lower() == 'cuda' and _ml_device == 'cuda':
        import torch
        if torch.cuda.is_available():
            torch.cuda.init()
            print(f"  PyTorch CUDA initialized before OpenMM Context (context-order fix)")
    platform = Platform.getPlatformByName(platform_name.upper())
    if platform_name.lower() == 'cuda':
        # NPT: use double precision to avoid posqCorrection CUDA_ERROR_ILLEGAL_ADDRESS
        precision = precision_override or ('double' if pressure_bar > 0 else 'mixed')
        properties = {
            'Precision': precision,
            'DeviceIndex': '0',
            'DisablePmeStream': 'true',
        }
        print(f"  Using CUDA platform (GPU acceleration, DeviceIndex=0)")
        if precision == 'double' and pressure_bar > 0:
            print(f"  Precision: double (NPT workaround for posqCorrection crash)")
    else:
        properties = {}
        print(f"  Using CPU platform")
    context = Context(system, integrator, platform, properties)
    
    # Initialize beads
    print(f"\n--- Initializing RPMD Beads ---")
    positions_nm = np.array(positions)
    
    # Set context positions first (required before getting velocities)
    pos_with_units = [Vec3(float(p[0]), float(p[1]), float(p[2])) * unit.nanometer for p in positions_nm]
    context.setPositions(pos_with_units)
    
    for i in range(num_beads):
        perturbation = np.random.randn(n_atoms, 3) * 0.001
        perturbed_pos = positions_nm + perturbation
        integrator.setPositions(i, perturbed_pos * unit.nanometer)
        print(f"  Bead {i}: positions set")
    
    # Set velocities
    context.setVelocitiesToTemperature(temperature_K * unit.kelvin)
    state = context.getState(getVelocities=True)
    base_velocities = state.getVelocities()
    
    for i in range(num_beads):
        velocities = []
        for v in base_velocities:
            vx = v[0].value_in_unit(unit.nanometer/unit.picosecond) + 0.01 * np.random.randn()
            vy = v[1].value_in_unit(unit.nanometer/unit.picosecond) + 0.01 * np.random.randn()
            vz = v[2].value_in_unit(unit.nanometer/unit.picosecond) + 0.01 * np.random.randn()
            velocities.append(Vec3(vx, vy, vz) * unit.nanometer/unit.picosecond)
        integrator.setVelocities(i, velocities)
    
    print(f"  Velocities initialized to {temperature_K} K")
    
    # Energy minimization to relax the structure
    print(f"\n--- Energy Minimization ---")
    print("  Minimizing energy to relax initial structure...")
    initial_state = integrator.getState(0, getEnergy=True)
    initial_pe = initial_state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"  Initial PE: {initial_pe:.2f} kJ/mol")

    # Sync context with bead 0 positions before minimization (RPMD stores positions in integrator)
    state0 = integrator.getState(0, getPositions=True)
    bead0_positions = state0.getPositions()
    context.setPositions(bead0_positions)

    # Perform energy minimization on bead 0, then copy to all beads
    from openmm import LocalEnergyMinimizer
    LocalEnergyMinimizer.minimize(context, tolerance=1.0, maxIterations=2000)
    
    # Get minimized positions and copy to all beads
    minimized_state = context.getState(getPositions=True, getEnergy=True)
    minimized_positions = minimized_state.getPositions()
    minimized_pe = minimized_state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"  Minimized PE: {minimized_pe:.2f} kJ/mol")
    print(f"  Energy reduction: {initial_pe - minimized_pe:.2f} kJ/mol")
    
    for i in range(num_beads):
        integrator.setPositions(i, minimized_positions)
    print(f"  Minimized structure copied to all {num_beads} beads")

    # Store for NVT RMSD check
    minimized_positions_nm = np.array([[p[0].value_in_unit(unit.nanometer), p[1].value_in_unit(unit.nanometer), p[2].value_in_unit(unit.nanometer)] for p in minimized_positions])
    
    # Get initial energies
    print(f"\n--- Initial Bead Energies ---")
    print("  Getting initial bead energies...")
    
    bead_energies = []
    for i in range(num_beads):
        state = integrator.getState(i, getEnergy=True)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        ke = state.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
        total = pe + ke
        bead_energies.append((pe, ke, total))
        print(f"  Bead {i:2d}: PE = {pe:10.2f} kJ/mol, KE = {ke:8.2f} kJ/mol, Total = {total:10.2f} kJ/mol")
    
    mean_pe = np.mean([e[0] for e in bead_energies])
    mean_ke = np.mean([e[1] for e in bead_energies])
    mean_total = np.mean([e[2] for e in bead_energies])
    print(f"\n  Mean: PE = {mean_pe:10.2f} kJ/mol, KE = {mean_ke:8.2f} kJ/mol, Total = {mean_total:10.2f} kJ/mol")
    
    # Setup PDB reporter (centroid only)
    print(f"\n--- Setting up PDB Reporter ---")
    pdb_file = output_path / f"ice_rpmd_{model_name.split('-')[0]}_T{temperature_K:.0f}_b{num_beads}.pdb"
    print(f"  PDB output: {pdb_file}")
    
    # Calculate steps
    dt_ps = dt_fs / 1000.0
    equilibration_steps = int(equilibration_ps / dt_ps)
    production_steps = int(production_ps / dt_ps)
    total_steps = equilibration_steps + production_steps
    
    report_interval_steps = max(1, int(report_interval_ps / dt_ps))
    pdb_interval_steps = max(1, int(pdb_interval_ps / dt_ps))
    
    print(f"  Report interval: {report_interval_steps} steps ({report_interval_steps * dt_ps:.2f} ps)")
    print(f"  PDB save interval: {pdb_interval_steps} steps ({pdb_interval_steps * dt_ps:.2f} ps)")
    
    # Run simulation
    print("\n" + "=" * 80)
    print("RUNNING SIMULATION")
    print("=" * 80)
    
    start_time = time.time()
    last_report_time = start_time
    
    step_counter = 0
    phase = "Production" if equilibration_steps == 0 else "Equilibration"

    energies = []
    temperatures = []
    densities = []
    rmsds = []  # For NVT ice stability
    minimized_positions_nm = None  # Store for RMSD
    
    pdb_file_handle = open(pdb_file, 'w')
    
    # Write PDB header with proper unit handling
    box_size_nm = box_vectors[0][0].value_in_unit(unit.nanometer)
    box_size_angstrom = box_size_nm * 10.0  # Convert nm to Angstrom for PDB
    pdb_file_handle.write("REMARK   Ice RPMD simulation - centroid trajectory\n")
    pdb_file_handle.write("REMARK   Model: %s\n" % model_name)
    pdb_file_handle.write("REMARK   Temperature: %.1f K, Beads: %d\n" % (temperature_K, num_beads))
    pdb_file_handle.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n" % (
        box_size_angstrom, box_size_angstrom, box_size_angstrom, 90.0, 90.0, 90.0))
    
    # Cache masses for temperature calculation
    masses = [system.getParticleMass(j).value_in_unit(unit.dalton) for j in range(n_atoms)]
    
    # Pre-allocate velocity buffer for reporting
    all_bead_velocities_buffer = np.zeros((num_beads, n_atoms, 3), dtype=np.float64)
    
    step = 0
    while step < total_steps:
        # Calculate next checkpoint
        next_report = ((step // report_interval_steps) + 1) * report_interval_steps
        next_pdb = ((step // pdb_interval_steps) + 1) * pdb_interval_steps if step >= equilibration_steps else total_steps
        next_stop = min(next_report, next_pdb, total_steps)
        chunk = next_stop - step
        
        # Run chunked steps
        integrator.step(chunk)
        step = next_stop
        step_counter = step
        
        # Switch to production phase (use >= so we catch it even if a chunk skips over equilibration_steps)
        if step >= equilibration_steps and phase == "Equilibration":
            phase = "Production"
            print(f"\n{'=' * 80}")
            print(f"PRODUCTION PHASE STARTED")
            print(f"{'=' * 80}\n")
        
        # Report progress
        if step % report_interval_steps == 0:
            current_wall_time = time.time()
            elapsed_since_last = current_wall_time - last_report_time
            last_report_time = current_wall_time
            
            # Get state from bead 0 for PE and positions
            state0 = integrator.getState(0, getEnergy=True, getPositions=True)
            current_pe = state0.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            
            # Calculate physical temperature from centroid kinetic energy
            # Physical T = (2 * KE_centroid) / (3 * N_atoms * k_B)
            centroid_ke = 0.0
            
            # Get velocities for all beads (write into pre-allocated buffer)
            for b in range(num_beads):
                state_b = integrator.getState(b, getVelocities=True)
                all_bead_velocities_buffer[b] = state_b.getVelocities(asNumpy=True).value_in_unit(unit.nanometer/unit.picosecond)
            
            centroid_velocities = np.mean(all_bead_velocities_buffer, axis=0)
            
            # KE = 0.5 * sum(m * v^2)
            # v is in nm/ps, m is in dalton (g/mol)
            # 1 dalton * (nm/ps)^2 = 1 g/mol * (1e-9 m / 1e-12 s)^2 = 1e-3 kg/mol * 1e6 m^2/s^2 = 1000 J/mol = 1 kJ/mol
            # So the result is directly in kJ/mol
            for j in range(n_atoms):
                v2 = np.sum(centroid_velocities[j]**2)
                centroid_ke += 0.5 * masses[j] * v2
            
            current_ke = centroid_ke # This is the physical KE
            current_total = current_pe + current_ke
            
            # Calculate temperature from physical KE
            # k_B = 8.314e-3 kJ/(mol*K)
            current_temp = (2.0 * current_ke) / (3.0 * n_atoms * 8.314e-3)
            
            # Calculate density if NPT
            if pressure_bar > 0:
                box = state0.getPeriodicBoxVectors()
                volume_nm3 = (box[0][0] * box[1][1] * box[2][2]).value_in_unit(unit.nanometer**3)
                volume_cm3 = volume_nm3 * 1e-21
                current_density = total_mass_g / volume_cm3
                densities.append(current_density)
            else:
                current_density = density

            # NVT: compute centroid RMSD from minimized structure (proxy for melting)
            if pressure_bar == 0 and minimized_positions_nm is not None:
                centroid_positions = np.zeros((n_atoms, 3))
                for b in range(num_beads):
                    state_b = integrator.getState(b, getPositions=True)
                    centroid_positions += state_b.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                centroid_positions /= num_beads
                diff = centroid_positions - minimized_positions_nm
                rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
                rmsds.append(rmsd)

            energies.append(current_total)
            temperatures.append(current_temp)
            
            # Calculate speed
            sim_time_ps = step * dt_ps
            progress_pct = (step / total_steps) * 100
            steps_per_sec = report_interval_steps / elapsed_since_last
            ns_per_day = steps_per_sec * dt_ps * 86400.0 / 1000.0
            
            print(f"[{progress_pct:5.1f}%] {phase:13s} | "
                  f"Step {step:7d}/{total_steps} ({sim_time_ps:6.1f} ps) | "
                  f"{ns_per_day:6.2f} ns/day")
            print(f"         PE={current_pe:10.1f} KE={current_ke:8.1f} Total={current_total:10.1f} kJ/mol | "
                  f"T={current_temp:6.1f} K | ρ={current_density:.3f} g/cm³")
            print("")
        
        # Save PDB frame (production only, at specified interval)
        if step > equilibration_steps and (step - equilibration_steps) % pdb_interval_steps == 0:
            # Get centroid positions
            centroid_positions = np.zeros((n_atoms, 3))
            for bead in range(num_beads):
                state = integrator.getState(bead, getPositions=True)
                pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                centroid_positions += pos
            centroid_positions /= num_beads
            
            # Write PDB frame manually
            model_num = (step - equilibration_steps) // pdb_interval_steps
            pdb_file_handle.write("MODEL     %4d\n" % model_num)
            
            atom_idx = 1
            for residue in topology.residues():
                for atom in residue.atoms():
                    pos_angstrom = centroid_positions[atom.index] * 10.0  # nm to Angstrom
                    # PDB chain ID must be 1 char; OpenMM uses "1","2",..."10"... so use first char
                    chain_char = (str(residue.chain.id) or ' ')[:1]
                    pdb_file_handle.write("ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n" % (
                        atom_idx,
                        atom.name,
                        residue.name,
                        chain_char,
                        residue.index + 1,
                        pos_angstrom[0],
                        pos_angstrom[1],
                        pos_angstrom[2],
                        1.0,  # occupancy
                        0.0,  # temperature factor
                        atom.element.symbol
                    ))
                    atom_idx += 1
            
            pdb_file_handle.write("ENDMDL\n")
    
    # Close PDB file
    pdb_file_handle.write("END\n")
    pdb_file_handle.close()
    
    total_elapsed = time.time() - start_time
    
    print("=" * 80)
    print("SIMULATION COMPLETE")
    print("=" * 80)
    print(f"Total time: {total_elapsed/60:.1f} minutes ({total_elapsed/3600:.2f} hours)")
    print(f"Average speed: {(total_steps * dt_ps * 86400.0) / (total_elapsed * 1000.0):.2f} ns/day")
    
    # Analysis
    print(f"\n--- Simulation Analysis ---")
    energies = np.array(energies)
    temperatures = np.array(temperatures)
    
    # Production phase only
    prod_start_idx = len(energies) // 2
    prod_energies = energies[prod_start_idx:]
    prod_temps = temperatures[prod_start_idx:]
    
    print(f"\nProduction phase statistics:")
    print(f"  Energy:")
    print(f"    Mean: {np.mean(prod_energies):.2f} kJ/mol")
    print(f"    Std:  {np.std(prod_energies):.2f} kJ/mol")
    print(f"    Drift: {(prod_energies[-1] - prod_energies[0])/prod_energies[0]*100:+.2f}%")
    print(f"  Temperature:")
    print(f"    Mean: {np.mean(prod_temps):.2f} K (target: {temperature_K} K)")
    print(f"    Std:  {np.std(prod_temps):.2f} K")
    
    if pressure_bar > 0 and len(densities) > 0:
        densities = np.array(densities)
        prod_densities = densities[prod_start_idx:]
        print(f"  Density:")
        print(f"    Mean: {np.mean(prod_densities):.3f} g/cm³")
        print(f"    Std:  {np.std(prod_densities):.3f} g/cm³")
        print(f"    Ice Ih reference: ~0.92 g/cm³")
        
        is_frozen = np.mean(prod_densities) > 0.85
        print(f"\n  Ice status: {'FROZEN' if is_frozen else 'MELTED'} (NPT: density > 0.85 g/cm³)")
    elif len(rmsds) > 0:
        rmsds_arr = np.array(rmsds)
        prod_rmsds = rmsds_arr[prod_start_idx:]
        mean_rmsd = np.mean(prod_rmsds)
        print(f"  RMSD from minimized (NVT):")
        print(f"    Mean: {mean_rmsd:.4f} nm")
        print(f"    Std:  {np.std(prod_rmsds):.4f} nm")
        print(f"    Threshold for frozen: < 0.15 nm")
        
        is_frozen = mean_rmsd < 0.15
        print(f"\n  Ice status: {'FROZEN' if is_frozen else 'MELTED'} (NVT: RMSD < 0.15 nm)")
    else:
        is_frozen = True
        print(f"\n  Ice status: Cannot determine (no density or RMSD data)")
    
    print(f"\n  Output files:")
    print(f"    Trajectory: {pdb_file}")
    
    print("\n" + "=" * 80)
    
    return is_frozen


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='UMA/eSEN Ice RPMD Simulation')
    parser.add_argument('--molecules', type=int, default=32,
                       help='Number of water molecules (default: 32)')
    parser.add_argument('--beads', type=int, default=8,
                       help='Number of RPMD beads (default: 8)')
    parser.add_argument('--temperature', type=float, default=243.0,
                       help='Temperature in K (default: 243.0)')
    parser.add_argument('--pressure', type=float, default=0.0,
                       help='Pressure in bar (default: 1.0, use 0 for NVT)')
    parser.add_argument('--dt', type=float, default=1.0,
                       help='Timestep in fs (default: 1.0, use 0.5-1.0 for stability)')
    parser.add_argument('--equil', type=float, default=5.0,
                       help='Equilibration time in ps (default: 5.0)')
    parser.add_argument('--prod', type=float, default=100.0,
                       help='Production time in ps (default: 100.0)')
    parser.add_argument('--model', type=str, default='uma-s-1-pythonforce-batch',
                       help='Model name (default: uma-s-1-pythonforce-batch, smallest UMA model)')
    parser.add_argument('--output', type=str, default='.',
                       help='Output directory (default: current directory)')
    parser.add_argument('--report-interval', type=float, default=1.0,
                       help='Console report interval in ps (default: 1.0)')
    parser.add_argument('--pdb-interval', type=float, default=1.0,
                       help='PDB frame save interval in ps (default: 1.0)')
    parser.add_argument('--platform', type=str, default='cuda', choices=['cuda', 'cpu'],
                       help='OpenMM platform: cuda (GPU) or cpu (default: cuda)')
    parser.add_argument('--ml-device', type=str, default=None,
                       help='ML model device: cuda (default with CUDA platform), cpu, or auto. '
                            'Lazy-load enables single-GPU use.')
    parser.add_argument('--precision', type=str, default=None,
                       help='OpenMM precision: single, mixed, or double. Default: double for NPT+CUDA, mixed otherwise.')
    
    args = parser.parse_args()
    
    try:
        success = run_simulation(
            num_molecules=args.molecules,
            num_beads=args.beads,
            temperature_K=args.temperature,
            pressure_bar=args.pressure,
            dt_fs=args.dt,
            equilibration_ps=args.equil,
            production_ps=args.prod,
            model_name=args.model,
            output_dir=args.output,
            report_interval_ps=args.report_interval,
            pdb_interval_ps=args.pdb_interval,
            platform_name=args.platform,
            ml_device=args.ml_device,
            precision_override=args.precision
        )
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"\nSIMULATION FAILED")
        print(f"Error: {type(e).__name__}: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
