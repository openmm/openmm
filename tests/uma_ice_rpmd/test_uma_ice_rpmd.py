#!/usr/bin/env python3
"""
RPMD simulation testing if ice stays frozen at 243K with UMA or eSEN potentials.
Uses proper ice Ih structure (GenIce CIF, GenIce2, or embedded CIF) in a periodic box.
PILE thermostat, NPT/NVT.

Run:
  python test_uma_ice_rpmd.py --input ice.cif --beads 8 --temperature 243 --dt 1.0 --equil 5 --prod 50
  python test_uma_ice_rpmd.py --box 1.0 --beads 8 ...
  python test_uma_ice_rpmd.py --molecules 32 --beads 8 ...

Input from GenIce: genice 1h --rep 2 2 2 --format cif > ice.cif
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

# Set DEBUG_LOGS before openmmml import so UMA batch profiling is enabled
if '--debug-logs' in sys.argv:
    os.environ['OPENMMML_DEBUG_LOGS'] = '1'

import numpy as np
import time
from pathlib import Path

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    _HAS_MATPLOTLIB = True
except ImportError:
    _HAS_MATPLOTLIB = False

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


def _molecules_for_box(box_nm, density_g_cm3=0.92):
    """Number of water molecules for cubic box L×L×L at ice Ih density."""
    vol_cm3 = (box_nm * 1e-7) ** 3  # 1 nm = 1e-7 cm, so 1 nm³ = 1e-21 cm³
    return max(1, int(round(density_g_cm3 * vol_cm3 / 18.015 * 6.022e23)))


def _hexagonal_to_orthorhombic(atoms_ase):
    """Transform hexagonal ice Ih cell to orthorhombic. Preserves Cartesian positions."""
    cell = np.array(atoms_ase.get_cell())
    pos = atoms_ase.get_positions()
    # Orthorhombic: a' = a, b' = 2b - a, c' = c (IUCr convention)
    cell_ortho = np.array([
        cell[0],
        2.0 * cell[1] - cell[0],
        cell[2],
    ])
    # frac_new @ cell_ortho = frac_old @ cell  =>  frac_new = pos @ inv(cell_ortho)
    frac_new = pos @ np.linalg.inv(cell_ortho.T)
    frac_wrapped = frac_new % 1.0
    from ase import Atoms
    out = Atoms(
        symbols=atoms_ase.get_chemical_symbols(),
        positions=frac_wrapped @ cell_ortho,  # back to Cartesian (same as pos, wrapped)
        cell=cell_ortho,
        pbc=True,
    )
    return out


def _ice_to_orthorhombic_then_scale(atoms_ase, num_molecules, box_nm=None, ice_type='1h'):
    """
    Transform ice to orthorhombic (1h only), replicate for roughly cubic supercell, uniformly scale.
    Preserves O-O distances (~2.75 A) and tetrahedral angles.
    For ice_type='1c', cell is already cubic; skip orthorhombic transform.
    Returns (topology, positions_nm, box_vectors).
    """
    from ase.build import make_supercell

    if ice_type == '1c':
        atoms_ortho = atoms_ase.copy()
    else:
        atoms_ortho = _hexagonal_to_orthorhombic(atoms_ase)
    cell = np.array(atoms_ortho.get_cell())
    Lx, Ly, Lz = np.linalg.norm(cell[0]), np.linalg.norm(cell[1]), np.linalg.norm(cell[2])

    # Replicate to get roughly cubic supercell with at least num_molecules
    n_mol_per_cell = len(atoms_ortho) // 3
    r_min = max(1, int(np.ceil((num_molecules / n_mol_per_cell) ** (1.0 / 3.0))))
    # Choose na, nb, nc so na*Lx ≈ nb*Ly ≈ nc*Lz
    ratios = np.array([Lx, Ly, Lz])
    ratios = ratios / np.min(ratios)
    r_approx = max(r_min, 1)
    na = max(1, int(round(r_approx / ratios[0])))
    nb = max(1, int(round(r_approx / ratios[1])))
    nc = max(1, int(round(r_approx / ratios[2])))
    P = np.diag([na, nb, nc])
    supercell = make_supercell(atoms_ortho, P)

    # Trim to num_molecules (first N by index preserves local order)
    mol_list = []
    symbols = supercell.get_chemical_symbols()
    pos_full = supercell.get_positions()
    seen = set()
    for i, s in enumerate(symbols):
        if s == 'O' and i not in seen:
            oidx = i
            dists = [(np.linalg.norm(pos_full[oidx] - pos_full[j]), j)
                     for j, sj in enumerate(symbols) if sj == 'H']
            dists.sort(key=lambda x: x[0])
            h1, h2 = dists[0][1], dists[1][1]
            mol_list.append([oidx, h1, h2])
            seen.update([oidx, h1, h2])
    mol_list = mol_list[:num_molecules]

    pos_ang = np.array([pos_full[idx] for mol in mol_list for idx in mol])
    cell_super = np.array(supercell.get_cell())
    vol_angstrom3 = np.abs(np.linalg.det(cell_super))
    L_current = vol_angstrom3 ** (1.0 / 3.0)

    if box_nm is not None:
        vol_target = (box_nm * 10.0) ** 3
        scale = (vol_target / vol_angstrom3) ** (1.0 / 3.0)
    else:
        scale = 1.0

    pos_scaled = pos_ang * scale
    cell_scaled = cell_super * scale

    # OpenMM requires reduced-form box vectors: a=(x,0,0), b=(bx,by,0), c=(cx,cy,cz)
    raw_vectors = (
        Vec3(cell_scaled[0, 0] * 0.1, cell_scaled[0, 1] * 0.1, cell_scaled[0, 2] * 0.1),
        Vec3(cell_scaled[1, 0] * 0.1, cell_scaled[1, 1] * 0.1, cell_scaled[1, 2] * 0.1),
        Vec3(cell_scaled[2, 0] * 0.1, cell_scaled[2, 1] * 0.1, cell_scaled[2, 2] * 0.1),
    )
    from openmm.app.internal.unitcell import reducePeriodicBoxVectors
    reduced = reducePeriodicBoxVectors(raw_vectors)
    vals = reduced.value_in_unit(unit.nanometer)
    box_vectors = tuple(Vec3(v[0], v[1], v[2]) * unit.nanometer for v in vals)
    positions_nm = (pos_scaled * 0.1).tolist()

    topology = app.Topology()
    for _ in mol_list:
        chain = topology.addChain()
        residue = topology.addResidue('HOH', chain)
        o_atom = topology.addAtom('O', app.Element.getBySymbol('O'), residue)
        h1_atom = topology.addAtom('H', app.Element.getBySymbol('H'), residue)
        h2_atom = topology.addAtom('H', app.Element.getBySymbol('H'), residue)
        topology.addBond(o_atom, h1_atom)
        topology.addBond(o_atom, h2_atom)

    topology.setPeriodicBoxVectors(box_vectors)
    return topology, positions_nm, box_vectors


def _ase_water_to_openmm(atoms_ase):
    """
    Convert ASE Atoms (water only) to OpenMM (topology, positions_nm, box_vectors).
    Uses structure as-is; no replication, trimming, or scaling.
    Assigns H to O by distance (nearest 2 H per O).
    """
    symbols = np.array(atoms_ase.get_chemical_symbols())
    pos_ang = np.array(atoms_ase.get_positions())
    cell = np.array(atoms_ase.get_cell())

    mol_list = []
    seen = set()
    for i, s in enumerate(symbols):
        if s == 'O' and i not in seen:
            dists = [(np.linalg.norm(pos_ang[i] - pos_ang[j]), j)
                     for j, sj in enumerate(symbols) if sj == 'H' and j not in seen]
            dists.sort(key=lambda x: x[0])
            if len(dists) < 2:
                raise ValueError(f"O at index {i} has < 2 H neighbors")
            h1, h2 = dists[0][1], dists[1][1]
            mol_list.append([i, h1, h2])
            seen.update([i, h1, h2])

    pos_ordered = np.array([pos_ang[idx] for mol in mol_list for idx in mol])

    raw_vectors = (
        Vec3(cell[0, 0] * 0.1, cell[0, 1] * 0.1, cell[0, 2] * 0.1),
        Vec3(cell[1, 0] * 0.1, cell[1, 1] * 0.1, cell[1, 2] * 0.1),
        Vec3(cell[2, 0] * 0.1, cell[2, 1] * 0.1, cell[2, 2] * 0.1),
    )
    from openmm.app.internal.unitcell import reducePeriodicBoxVectors
    reduced = reducePeriodicBoxVectors(raw_vectors)
    vals = reduced.value_in_unit(unit.nanometer)
    box_vectors = tuple(Vec3(v[0], v[1], v[2]) * unit.nanometer for v in vals)

    topology = app.Topology()
    for _ in mol_list:
        chain = topology.addChain()
        residue = topology.addResidue('HOH', chain)
        o_atom = topology.addAtom('O', app.Element.getBySymbol('O'), residue)
        h1_atom = topology.addAtom('H', app.Element.getBySymbol('H'), residue)
        h2_atom = topology.addAtom('H', app.Element.getBySymbol('H'), residue)
        topology.addBond(o_atom, h1_atom)
        topology.addBond(o_atom, h2_atom)

    topology.setPeriodicBoxVectors(box_vectors)
    positions_nm = (pos_ordered * 0.1).tolist()
    return topology, positions_nm, box_vectors


def _load_structure_from_file(filepath, max_molecules=None):
    """
    Load ice/water structure from CIF or PDB file.
    Returns (topology, positions_nm, box_vectors).
    Supports GenIce CIF output (fractional coords) and standard PDB.
    max_molecules: if set, trim to first N molecules and scale box (preserves density). For CUDA OOM.
    """
    from ase.io import read

    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"Structure file not found: {filepath}")

    suffix = path.suffix.lower()
    if suffix in ('.cif', '.pdb', '.xyz'):
        atoms = read(str(path))
    else:
        raise ValueError(f"Unsupported structure format: {suffix}. Use .cif, .pdb, or .xyz")

    # XYZ and some formats lack cell; set cubic box for gas-phase structures.
    # ASE set_cell() uses Angstroms.
    # Single molecule (3 atoms): use 5 Å box to keep molecule stable (avoids "no edges" when H drifts in large box).
    # Multi-molecule: use 12 Å (1.2 nm) per molecule region.
    cell = np.array(atoms.get_cell())
    if cell.size < 9 or np.abs(np.linalg.det(cell)) < 1e-10:
        box_angstrom = 5.0 if len(atoms) == 3 else 12.0
        atoms.set_cell([box_angstrom, box_angstrom, box_angstrom])
        atoms.set_pbc(True)
        atoms.center()

    topo, positions, box = _ase_water_to_openmm(atoms)
    if max_molecules is not None:
        n_mol = len(positions) // 3
        if n_mol > max_molecules:
            n_keep = max_molecules
            n_atoms_keep = n_keep * 3
            positions = positions[:n_atoms_keep]
            box_matrix = np.array([[box[i][j].value_in_unit(unit.nanometer) for j in range(3)] for i in range(3)])
            scale = (n_keep / n_mol) ** (1.0 / 3.0)
            box_matrix *= scale
            positions_nm = np.array([[p[0], p[1], p[2]] for p in positions]) * scale
            box_vectors = (
                Vec3(box_matrix[0, 0], box_matrix[0, 1], box_matrix[0, 2]) * unit.nanometer,
                Vec3(box_matrix[1, 0], box_matrix[1, 1], box_matrix[1, 2]) * unit.nanometer,
                Vec3(box_matrix[2, 0], box_matrix[2, 1], box_matrix[2, 2]) * unit.nanometer,
            )
            topo = app.Topology()
            for _ in range(n_keep):
                chain = topo.addChain()
                res = topo.addResidue('HOH', chain)
                o_atom = topo.addAtom('O', app.Element.getBySymbol('O'), res)
                h1 = topo.addAtom('H', app.Element.getBySymbol('H'), res)
                h2 = topo.addAtom('H', app.Element.getBySymbol('H'), res)
                topo.addBond(o_atom, h1)
                topo.addBond(o_atom, h2)
            topo.setPeriodicBoxVectors(box_vectors)
            return topo, positions_nm.tolist(), box_vectors
    return topo, positions, box


def _create_ice_from_ase(atoms_ase, num_molecules, box_nm=None, ice_type='1h'):
    """Convert ASE Atoms (ice Ih or Ic) to OpenMM topology and positions.
    Uses orthorhombic transformation (1h only) + uniform scale to preserve structure.
    """
    return _ice_to_orthorhombic_then_scale(atoms_ase, num_molecules, box_nm=box_nm, ice_type=ice_type)


def _create_ice_genice(num_molecules, box_nm=None, ice_type='1h'):
    """Try GenIce (or GenIce2) subprocess. Returns (topology, positions, box_vectors) or None on failure.
    ice_type: '1h' (hexagonal) or '1c' (cubic ice, natural cubic box).
    Tries GenIce 1.x first (genice --format cif), then GenIce2, then embedded CIF.
    """
    from ase.io import read

    n_per_cell = 12 if ice_type == '1h' else 8  # ice Ic has 8 molecules per cubic cell (Fd-3m)
    r = max(1, int(np.ceil((num_molecules / n_per_cell) ** (1 / 3))))
    rep_args = [str(r), str(r), str(r)]

    # Try GenIce 1.x (genice): outputs CIF to stdout
    for cmd, out_format in [
        (['genice', ice_type, '--rep'] + rep_args + ['--format', 'cif'], 'cif'),
        (['genice2', ice_type, '--rep'] + rep_args + ['--format', 'cif'], 'cif'),
        (['genice2', ice_type, '--rep'] + rep_args + ['-f', 'y'], 'yaml'),
    ]:
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=30,
            )
            stdout = result.stdout.strip()
            if result.returncode != 0 or not stdout:
                continue
            with tempfile.NamedTemporaryFile(suffix='.cif' if out_format == 'cif' else '.y', delete=False, mode='w') as f:
                f.write(stdout)
                outpath = f.name
            try:
                atoms = read(outpath, format=out_format)
            finally:
                os.unlink(outpath)
            if len(atoms) >= num_molecules * 3:
                return _create_ice_from_ase(atoms, num_molecules, box_nm=box_nm, ice_type=ice_type)
        except Exception:
            continue
    return None


def _create_ice_embedded(num_molecules, box_nm=None):
    """Create ice Ih from embedded CIF (12 molecules) + replication."""
    from ase.io import read
    from io import StringIO

    atoms = read(StringIO(_ICE_IH_CIF), format='cif')
    n_base = len(atoms) // 3
    if num_molecules <= n_base:
        return _create_ice_from_ase(atoms, num_molecules, box_nm=box_nm, ice_type='1h')
    r = max(1, int(np.ceil((num_molecules / n_base) ** (1 / 3))))
    from ase.build import make_supercell
    P = np.diag([r, r, r])
    supercell = make_supercell(atoms, P)
    return _create_ice_from_ase(supercell, num_molecules, box_nm=box_nm, ice_type='1h')


def create_ice_structure(num_molecules=None, box_nm=None, ice_type='1h', input_file=None,
                         max_molecules=None):
    """
    Create ice structure with tetrahedral hydrogen-bonding network.

    If input_file is given, load structure directly from CIF/PDB/XYZ (e.g. GenIce output).
    max_molecules: when loading from file, trim to first N molecules and scale box (for CUDA OOM).
    Otherwise specify num_molecules OR box_nm.
    ice_type: '1h' (hexagonal Ih) or '1c' (cubic Ic, natural cubic box).

    Returns
    -------
    topology : app.Topology
    positions : list of Vec3
    box_vectors : tuple of Vec3
    """
    if input_file is not None:
        print(f"\n--- Loading Ice Structure from File ---")
        print(f"  Input: {input_file}")
        topology, positions, box_vectors = _load_structure_from_file(input_file, max_molecules=max_molecules)
        n_atoms = len(positions)
        if max_molecules is not None:
            print(f"  Trimmed to {max_molecules} molecules (--molecules, for GPU memory)")
    else:
        if box_nm is not None:
            num_molecules = _molecules_for_box(box_nm)
        elif num_molecules is None:
            num_molecules = 32

        print(f"\n--- Creating Ice Structure ---")
        print(f"  Ice type: {ice_type}")
        if box_nm is not None:
            print(f"  Target box volume: {box_nm**3:.4f} nm³")
        print(f"  Target molecules: {num_molecules}")

        result = _create_ice_genice(num_molecules, box_nm=box_nm, ice_type=ice_type)
        if result is not None:
            topology, positions, box_vectors = result
            print(f"  Source: GenIce (proton-disordered)")
        else:
            if ice_type == '1c':
                print(f"  GenIce2 failed for 1c; falling back to 1h")
                ice_type = '1h'
            topology, positions, box_vectors = _create_ice_embedded(num_molecules, box_nm=box_nm)
            print(f"  Source: Embedded CIF (ice Ih from Avogadro/COD)")

        n_atoms = len(positions)
    mol_count = n_atoms // 3
    print(f"  Actual molecules: {mol_count}")
    print(f"  Total atoms: {n_atoms}")

    box_matrix = np.array([[box_vectors[i][j].value_in_unit(unit.nanometer) for j in range(3)] for i in range(3)])
    vol_nm3 = np.abs(np.linalg.det(box_matrix))
    L_sides = np.linalg.norm(box_matrix, axis=1) * 10  # nm to Angstrom
    total_mass_g = (mol_count * 18.015) * 1.66054e-24
    density = total_mass_g / (vol_nm3 * 1e-21)  # 1 nm³ = 1e-21 cm³
    print(f"  Box: {L_sides[0]:.2f} x {L_sides[1]:.2f} x {L_sides[2]:.2f} Å")
    print(f"  Initial density: {density:.3f} g/cm³")
    print(f"  Expected ice Ih: ~0.92 g/cm³")

    return topology, positions, box_vectors


print("All required packages loaded successfully")


def run_simulation(num_molecules=None, box_nm=None, ice_type='1h', input_file=None,
                   max_molecules=None,
                   num_beads=8, temperature_K=243.0,
                   pressure_bar=1.0, dt_fs=1.0, equilibration_ps=5.0,
                   production_ps=100.0, model_name='uma-s-1-pythonforce-batch',
                   output_dir='.', report_interval_ps=1.0, pdb_interval_ps=0.2,
                   platform_name='cuda', ml_device=None, precision_override=None,
                   minimal_report=False, optimize_inference=False,
                   optimize_inference_tf32_only=False):
    """
    Run RPMD simulation of ice.
    
    Parameters
    ----------
    num_molecules : int, optional
        Number of water molecules (ignored if box_nm given)
    box_nm : float, optional
        Cubic box side in nm; overrides num_molecules
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
        PDB frame save interval in picoseconds (default: 0.2, 5 frames/ps)
    platform_name : str
        OpenMM platform: 'cuda' (GPU) or 'cpu' (default: 'cuda')
    ml_device : str or None
        ML model device: 'cpu', 'cuda', or None. When None with CUDA platform,
        defaults to 'cpu' to avoid PyTorch/OpenMM CUDA context conflict.
    minimal_report : bool
        If True, skip getState(Energy) in report block; only print progress and
        ns/day. Reduces extra UMA inference for faster runs (default: False).
    optimize_inference : bool
        If True, use fast inference settings (no activation checkpointing,
        TF32) to speed up UMA force evaluation. Increases GPU memory; reduce
        beads/molecules if OOM. Use optimize_inference_tf32_only on constrained GPUs.
    optimize_inference_tf32_only : bool
        If True, enable only TF32 (keeps activation checkpointing). Lower memory,
        some speedup (default: False).

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
    if input_file is not None:
        print(f"Input structure: {input_file}")
    else:
        if box_nm is not None:
            num_molecules = _molecules_for_box(box_nm)
        if num_molecules is None:
            num_molecules = 32
        print(f"Molecules: {num_molecules}" + (f" (from box {box_nm} nm)" if box_nm else ""))
    print(f"Equilibration: {equilibration_ps} ps")
    print(f"Production: {production_ps} ps")
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Create structure
    topology, positions, box_vectors = create_ice_structure(
        num_molecules=num_molecules, box_nm=box_nm, ice_type=ice_type, input_file=input_file,
        max_molecules=max_molecules
    )
    n_atoms = len(positions)
    num_molecules = n_atoms // 3
    
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

    create_kwargs = {'task_name': 'omol', 'charge': 0, 'spin': 1, 'device': _ml_device}
    if optimize_inference or optimize_inference_tf32_only:
        from dataclasses import replace
        from fairchem.core.units.mlip_unit.api.inference import inference_settings_default
        d = inference_settings_default()
        if optimize_inference_tf32_only and not optimize_inference:
            create_kwargs['inference_settings'] = replace(d, compile=False, tf32=True)
            print(f"  Inference optimizations: tf32=on only (lower memory)")
        else:
            create_kwargs['inference_settings'] = replace(
                d, activation_checkpointing=False, compile=False, tf32=True
            )
            print(f"  Inference optimizations: checkpoint=off, tf32=on (compile=off, UMA equivariant layers incompatible)")
    system = potential.createSystem(topology, **create_kwargs)
    
    print(f"  System uses PBC: {system.usesPeriodicBoundaryConditions()}")
    
    # Add barostat if NPT. Use on all platforms including CUDA.
    _use_barostat = pressure_bar > 0
    if _use_barostat:
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
        # Barostat disabled on CUDA (see above). Use mixed precision for NVT.
        precision = precision_override or 'mixed'
        properties = {
            'Precision': precision,
            'DeviceIndex': '0',
            'DisablePmeStream': 'true',
        }
        print(f"  Using CUDA platform (GPU acceleration, DeviceIndex=0)")
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
    max_min_iter = 200 if n_atoms <= 3 else 2000  # Faster for single-molecule smoke tests
    LocalEnergyMinimizer.minimize(context, tolerance=1.0, maxIterations=max_min_iter)
    
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
    if minimal_report and report_interval_steps < 100:
        report_interval_steps = 100
        print(f"  Minimal mode: report interval clamped to {report_interval_steps} steps for speed")
    pdb_interval_steps = max(1, int(pdb_interval_ps / dt_ps))
    
    print(f"  Report interval: {report_interval_steps} steps ({report_interval_steps * dt_ps:.2f} ps)")
    print(f"  PDB save interval: {pdb_interval_steps} steps ({pdb_interval_steps * dt_ps:.2f} ps)")
    
    # Run simulation
    print("\n" + "=" * 80)
    print("RUNNING SIMULATION")
    print("=" * 80)
    if minimal_report and platform_name.lower() == 'cuda' and os.environ.get('OPENMM_CUDA_FAST_EXTERNAL_CALL') != '1':
        print("  Tip: OPENMM_CUDA_FAST_EXTERNAL_CALL=1 may improve speed (reduces GPU sync before Python callback)")
    
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
    
    # Write PDB header with proper unit handling (support orthorhombic boxes from GenIce CIF)
    _bm = box_matrix  # nm, 3x3
    a_ang = np.linalg.norm(_bm[0]) * 10.0
    b_ang = np.linalg.norm(_bm[1]) * 10.0
    c_ang = np.linalg.norm(_bm[2]) * 10.0
    alpha = np.degrees(np.arccos(np.clip(np.dot(_bm[1], _bm[2]) / (np.linalg.norm(_bm[1]) * np.linalg.norm(_bm[2]) + 1e-10), -1, 1)))
    beta = np.degrees(np.arccos(np.clip(np.dot(_bm[0], _bm[2]) / (np.linalg.norm(_bm[0]) * np.linalg.norm(_bm[2]) + 1e-10), -1, 1)))
    gamma = np.degrees(np.arccos(np.clip(np.dot(_bm[0], _bm[1]) / (np.linalg.norm(_bm[0]) * np.linalg.norm(_bm[1]) + 1e-10), -1, 1)))
    pdb_file_handle.write("REMARK   Ice RPMD simulation - centroid trajectory\n")
    pdb_file_handle.write("REMARK   Model: %s\n" % model_name)
    pdb_file_handle.write("REMARK   Temperature: %.1f K, Beads: %d\n" % (temperature_K, num_beads))
    pdb_file_handle.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n" % (
        a_ang, b_ang, c_ang, alpha, beta, gamma))
    
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
            
            sim_time_ps = step * dt_ps
            progress_pct = (step / total_steps) * 100
            steps_per_sec = report_interval_steps / elapsed_since_last
            ns_per_day = steps_per_sec * dt_ps * 86400.0 / 1000.0
            
            if minimal_report:
                print(f"[{progress_pct:5.1f}%] {phase:13s} | "
                      f"Step {step:7d}/{total_steps} ({sim_time_ps:6.1f} ps) | "
                      f"{ns_per_day:6.2f} ns/day")
            else:
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
                
                # Calculate density if NPT (barostat active)
                if _use_barostat:
                    box = state0.getPeriodicBoxVectors()
                    volume_nm3 = (box[0][0] * box[1][1] * box[2][2]).value_in_unit(unit.nanometer**3)
                    volume_cm3 = volume_nm3 * 1e-21
                    current_density = total_mass_g / volume_cm3
                    densities.append(current_density)
                else:
                    current_density = density

                # NVT: compute centroid RMSD from minimized structure (proxy for melting)
                if not _use_barostat and minimized_positions_nm is not None:
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
                
                print(f"[{progress_pct:5.1f}%] {phase:13s} | "
                      f"Step {step:7d}/{total_steps} ({sim_time_ps:6.1f} ps) | "
                      f"{ns_per_day:6.2f} ns/day")
                print(f"         PE={current_pe:10.1f} KE={current_ke:8.1f} Total={current_total:10.1f} kJ/mol | "
                      f"T={current_temp:6.1f} K | ρ={current_density:.3f} g/cm³")
            print("")
        
        # Save PDB frame (production only, at specified interval)
        if step > equilibration_steps and (step - equilibration_steps) % pdb_interval_steps == 0:
            # Get centroid positions (OpenMM returns unwrapped; wrap into [0,1) frac for PBC)
            centroid_positions = np.zeros((n_atoms, 3))
            for bead in range(num_beads):
                state = integrator.getState(bead, getPositions=True)
                pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                centroid_positions += pos
            centroid_positions /= num_beads
            # Wrap into primary cell so PDB viewers show compact structure
            frac = centroid_positions @ np.linalg.inv(_bm)
            frac_wrapped = frac - np.floor(frac)
            centroid_positions = frac_wrapped @ _bm
            
            # Write PDB frame manually
            model_num = (step - equilibration_steps) // pdb_interval_steps
            pdb_file_handle.write("MODEL     %4d\n" % model_num)
            
            atom_idx = 1
            for residue in topology.residues():
                for atom in residue.atoms():
                    pos_angstrom = centroid_positions[atom.index] * 10.0  # nm to Angstrom
                    # PDB chain ID must be 1 char; use %s (not %c) to avoid TypeError with non-char chain.id
                    chain_str = (str(residue.chain.id) or ' ')[:1]
                    pdb_file_handle.write("ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n" % (
                        atom_idx,
                        atom.name,
                        residue.name,
                        chain_str,
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
    if minimal_report:
        print("  Minimal mode: no energy/temp analysis.")
        is_frozen = True
    else:
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
        
        if _use_barostat and len(densities) > 0:
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
    
    # NPZ and analysis PNG (skip in minimal mode)
    npz_path = None
    png_path = None
    if not minimal_report and len(energies) > 0:
        base_name = f"uma_ice_rpmd_T{temperature_K:.0f}K_b{num_beads}"
        npz_path = output_path / f"{base_name}.npz"
        
        times_ps = np.arange(len(energies)) * report_interval_steps * dt_ps
        rmsds_arr = np.array(rmsds) if rmsds else np.array([], dtype=np.float64)
        densities_arr = np.array(densities) if densities else np.array([], dtype=np.float64)
        
        # Final centroid positions
        final_centroid = np.zeros((n_atoms, 3))
        for bead in range(num_beads):
            state = integrator.getState(bead, getPositions=True)
            final_centroid += state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        final_centroid /= num_beads
        
        np.savez(
            npz_path,
            times=times_ps,
            rmsds=rmsds_arr,
            densities=densities_arr,
            energies=np.array(energies),
            temperatures=np.array(temperatures),
            initial_positions=minimized_positions_nm if minimized_positions_nm is not None else np.zeros((n_atoms, 3)),
            final_positions=final_centroid,
            is_frozen=np.array(is_frozen),
        )
        
        # Analysis PNG (2x2 subplots: RMSD, density, energy, temperature)
        if _HAS_MATPLOTLIB:
            png_path = output_path / f"{base_name}_analysis.png"
            fig, axes = plt.subplots(2, 2, figsize=(10, 8))
            # axes[0]: RMSD, [1]: Density, [2]: Energy, [3]: Temperature
            if len(rmsds_arr) > 0:
                axes[0, 0].plot(times_ps[:len(rmsds_arr)], rmsds_arr * 10, 'b-')  # nm -> Å
                axes[0, 0].set_xlabel("Time (ps)")
                axes[0, 0].set_ylabel("RMSD (Å)")
                axes[0, 0].set_title("RMSD vs time")
                axes[0, 0].grid(True, alpha=0.3)
            else:
                axes[0, 0].set_visible(False)
            if len(densities_arr) > 0:
                axes[0, 1].plot(times_ps[:len(densities_arr)], densities_arr, 'g-')
                axes[0, 1].set_xlabel("Time (ps)")
                axes[0, 1].set_ylabel("Density (g/cm³)")
                axes[0, 1].set_title("Density vs time")
                axes[0, 1].axhline(0.92, color='gray', linestyle='--', alpha=0.7, label="Ice Ih ref")
                axes[0, 1].grid(True, alpha=0.3)
            else:
                axes[0, 1].set_visible(False)
            axes[1, 0].plot(times_ps, energies, 'r-')
            axes[1, 0].set_xlabel("Time (ps)")
            axes[1, 0].set_ylabel("Total energy (kJ/mol)")
            axes[1, 0].set_title("Total energy vs time")
            axes[1, 0].grid(True, alpha=0.3)
            axes[1, 1].plot(times_ps, temperatures, 'm-')
            axes[1, 1].axhline(temperature_K, color='gray', linestyle='--', alpha=0.7)
            axes[1, 1].set_xlabel("Time (ps)")
            axes[1, 1].set_ylabel("Temperature (K)")
            axes[1, 1].set_title("Temperature vs time")
            axes[1, 1].grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(png_path, dpi=150, bbox_inches='tight')
            plt.close()
    
    print(f"\n  Output files:")
    print(f"    Trajectory: {pdb_file}")
    if npz_path is not None:
        print(f"    NPZ data: {npz_path}")
    if png_path is not None:
        print(f"    Analysis PNG: {png_path}")
    
    print("\n" + "=" * 80)
    
    return is_frozen


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='UMA/eSEN Ice RPMD Simulation')
    parser.add_argument('--input', '-i', type=str, default=None,
                       help='Input structure file (CIF, PDB, or XYZ). e.g. ice.cif from GenIce. Overrides --molecules/--box.')
    parser.add_argument('--molecules', type=int, default=None,
                       help='With --input: trim to N molecules (for CUDA OOM). Else: target count. Default: from --box or 32')
    parser.add_argument('--box', type=float, default=None,
                       help='Cubic box side length in nm (e.g. 0.5, 1.0). Overrides --molecules.')
    parser.add_argument('--ice-type', type=str, default='1h', choices=['1h', '1c'],
                       help='Ice type: 1h (hexagonal Ih) or 1c (cubic Ic, natural cubic box). Default: 1h')
    parser.add_argument('--beads', type=int, default=8,
                       help='RPMD beads (default: 8). Use 1 for classical MD / LAMMPS-like UMA cost per step.')
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
    parser.add_argument('--model', type=str, default='uma-s-1p1-pythonforce-batch',
                       help='Model name (default: uma-s-1p1; HF no longer has uma-s-1.pt)')
    parser.add_argument('--output', type=str, default='.',
                       help='Output directory (default: current directory)')
    parser.add_argument('--report-interval', type=float, default=1.0,
                       help='Console report interval in ps (default: 1.0)')
    parser.add_argument('--pdb-interval', type=float, default=0.2,
                       help='PDB frame save interval in ps (default: 0.2, 5 frames/ps)')
    parser.add_argument('--platform', type=str, default='cuda', choices=['cuda', 'cpu'],
                       help='OpenMM platform: cuda (GPU) or cpu (default: cuda)')
    parser.add_argument('--ml-device', type=str, default=None,
                       help='ML model device: cuda (default with CUDA platform), cpu, or auto. '
                            'Lazy-load enables single-GPU use.')
    parser.add_argument('--precision', type=str, default=None,
                       help='OpenMM precision: single, mixed, or double. Default: double for NPT+CUDA, mixed otherwise.')
    parser.add_argument('--minimal', action='store_true',
                       help='Minimal reporting: only progress and ns/day. Avoids getState(Energy) for faster runs. '
                            'With CUDA, set OPENMM_CUDA_FAST_EXTERNAL_CALL=1 for extra speed. Use --report-interval 1.0 or higher.')
    parser.add_argument('--debug-logs', action='store_true',
                       help='Enable OPENMMML_DEBUG_LOGS for UMA batch timing breakdown (prep/atomicdata/inference). Use for profiling.')
    parser.add_argument('--optimize-inference', action='store_true',
                       help='Fast inference: no activation checkpointing, TF32. Increases GPU memory; use --optimize-inference-tf32-only if OOM.')
    parser.add_argument('--optimize-inference-tf32-only', action='store_true',
                       help='Enable only TF32 (keeps checkpointing). Lower memory, some speedup. Use if --optimize-inference causes OOM.')
    
    args = parser.parse_args()

    if args.debug_logs:
        import openmmml.models.umapotential_pythonforce_batch as _umab
        _umab.DEBUG_LOGS = True

    if args.input is not None and args.box is not None:
        print("Warning: --input given; ignoring --box")
    # --molecules with --input: trim to N molecules (for CUDA OOM). Without --input: target count.
    num_mol = args.molecules
    max_mol = args.molecules if args.input is not None else None

    try:
        success = run_simulation(
            num_molecules=num_mol if args.input is None else None,
            box_nm=args.box if args.input is None else None,
            input_file=args.input,
            max_molecules=max_mol,
            ice_type=args.ice_type,
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
            minimal_report=args.minimal,
            platform_name=args.platform,
            ml_device=args.ml_device,
            precision_override=args.precision,
            optimize_inference=args.optimize_inference,
            optimize_inference_tf32_only=args.optimize_inference_tf32_only
        )
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"\nSIMULATION FAILED")
        print(f"Error: {type(e).__name__}: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
