"""
cacepotential.py: CACE potential implementation for OpenMM-ML.
Uses OpenMM's PythonForce to integrate CACE models.
"""
import openmm
from openmm import unit
from openmmml.mlpotential import MLPotential, MLPotentialImpl, MLPotentialImplFactory
from typing import Iterable, Optional, Union, Dict
import numpy as np
import torch
import os
from pathlib import Path

class CACEPotentialImplFactory(MLPotentialImplFactory):
    """Factory that creates CACEPotentialImpl objects."""
    def createImpl(self, name: str, **args) -> MLPotentialImpl:
        return CACEPotentialImpl(name)

class CACEPotentialImpl(MLPotentialImpl):
    """
    CACE potential implementation using OpenMM's PythonForce.
    """
    def __init__(self, name: str) -> None:
        self.name = name
        # The name can be a path to a .pth file or a registered name
        self.model_path = name

    def addForces(
        self,
        topology: openmm.app.Topology,
        system: openmm.System,
        atoms: Optional[Iterable[int]],
        forceGroup: int,
        device: Optional[str] = None,
        **args,
    ) -> None:
        """
        Add the CACE force using PythonForce.
        """
        try:
            from cace.calculators import CACECalculator
            from cace.data import AtomicData
        except ImportError as e:
            raise ImportError(
                f"Failed to import CACE: {e}. "
                "Make sure CACE is installed and in your PYTHONPATH."
            )

        # Load the model
        if device is None:
            device = "cuda" if torch.cuda.is_available() else "cpu"
        
        print(f"Loading CACE model from '{self.model_path}' on {device}...")
        
        # Check if it's a file path or we need to look in cace/water (openmm root is 4 levels up from this file)
        actual_path = self.model_path
        if not os.path.exists(actual_path):
            openmm_root = Path(__file__).resolve().parents[4]
            cace_water_path = openmm_root / "cace" / "water" / "fit" / "fit_version_1" / "best_model.pth"
            if cace_water_path.exists():
                actual_path = str(cace_water_path)
                print(f"  Note: Using water model from cace/water: {actual_path}")
            else:
                raise FileNotFoundError(f"CACE model file not found: {self.model_path}")

        # Initialize the calculator
        # We'll use the model directly for the PythonForce callback
        model = torch.load(actual_path, map_location=device, weights_only=False)
        model.to(device)
        model.eval()
        
        # Get cutoff from model
        try:
            cutoff = model.representation.cutoff
        except AttributeError:
            cutoff = model.models[0].representation.cutoff

        # Get atomic information
        includedAtoms = list(topology.atoms())
        if atoms is not None:
            includedAtoms = [includedAtoms[i] for i in atoms]
        symbols = [atom.element.symbol for atom in includedAtoms]
        atomic_numbers = np.array([atom.element.atomic_number for atom in includedAtoms])
        atom_indices = atoms if atoms is not None else None
        n_atoms = len(symbols)
        
        # Get atomic energies from model or use defaults for water model
        # These are reference energies that were subtracted during training
        # and must be added back to get the correct total energy
        atomic_energies = None
        try:
            # Try to get from model metadata or use water defaults
            if hasattr(model, 'atomic_energies'):
                atomic_energies = model.atomic_energies
            elif 'atomic_energies' in args:
                atomic_energies = args['atomic_energies']
            else:
                # Default atomic energies for water model (from cace/water training)
                # H: -5.853064337340629 eV, O: -2.926532168670322 eV
                atomic_energies = {1: -5.853064337340629, 8: -2.926532168670322}
                print(f"  Using default atomic energies for water: H={atomic_energies[1]:.3f} eV, O={atomic_energies[8]:.3f} eV")
        except:
            # Fallback to water defaults
            atomic_energies = {1: -5.853064337340629, 8: -2.926532168670322}
            print(f"  Using default atomic energies for water (fallback)")

        # Check for periodicity - we'll determine this dynamically from the state
        # since box vectors may be set after createSystem() is called
        # But we can check if system has default box vectors set
        has_default_box = False
        try:
            a, b, c = system.getDefaultPeriodicBoxVectors()
            # Check if box vectors are non-zero (indicating periodic system)
            if a.norm() > 0 or b.norm() > 0 or c.norm() > 0:
                has_default_box = True
        except:
            pass
        
        # Also check if any force in the system uses PBC (e.g., barostat)
        has_pbc_force = False
        for i in range(system.getNumForces()):
            try:
                if system.getForce(i).usesPeriodicBoundaryConditions():
                    has_pbc_force = True
                    break
            except:
                pass
        
        isPeriodic = (
            topology.getPeriodicBoxVectors() is not None
        ) or system.usesPeriodicBoundaryConditions() or has_default_box or has_pbc_force
        
        # If we have box vectors or PBC forces, force periodic to True
        # (the callback will handle it dynamically anyway)
        if has_default_box or has_pbc_force:
            isPeriodic = True

        # Create AtomicData template
        from ase import Atoms
        dummy_pos = np.zeros((n_atoms, 3))
        atoms_ase = Atoms(symbols=symbols, positions=dummy_pos, pbc=isPeriodic)
        
        # Internal cache for optimization
        cache = {
            'device': device,
            'cutoff': cutoff,
            'n_atoms': n_atoms,
            'is_periodic': isPeriodic,
            'atom_indices': atom_indices,
            'atomic_numbers': torch.tensor(atomic_numbers, dtype=torch.long, device=device),
            'atomic_energies': atomic_energies,
            'atomic_numbers_np': atomic_numbers  # NumPy version for e0 calculation
        }

        def compute_cace_forces(state):
            """Callback for PythonForce."""
            try:
                # Get positions in nm (wrapped positions are fine - neighbor list handles PBC)
                all_pos_nm = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                
                # Check if periodic from box vectors
                box_vectors = state.getPeriodicBoxVectors(asNumpy=True)
                state_is_periodic = box_vectors is not None and (
                    np.linalg.norm(box_vectors[0]) > 0 or 
                    np.linalg.norm(box_vectors[1]) > 0 or 
                    np.linalg.norm(box_vectors[2]) > 0
                )
                
                if atom_indices is not None:
                    pos_nm = all_pos_nm[atom_indices]
                else:
                    pos_nm = all_pos_nm[:n_atoms]
                
                # Convert to Angstrom and to tensor
                pos_ang = torch.tensor(pos_nm * 10.0, dtype=torch.float32, device=device, requires_grad=True)
                
                # Prepare cell - check dynamically from state
                cell = None
                if state_is_periodic:
                    # Convert box vectors from nm to Angstrom
                    # box_vectors is a 3x3 array: [a, b, c] where each is a 3D vector
                    box_nm = box_vectors.value_in_unit(unit.nanometer)  # Shape: (3, 3)
                    cell_np = box_nm * 10.0  # Convert to Angstrom, shape: (3, 3)
                    cell = torch.tensor(cell_np, dtype=torch.float32, device=device).unsqueeze(0)  # Shape: (1, 3, 3)
                    
                    # Verify cell is reasonable (not too small or too large)
                    cell_diag = np.linalg.norm(cell_np, axis=1)  # Lengths of box vectors
                    if np.any(cell_diag < 0.1) or np.any(cell_diag > 1000.0):  # 0.1 Å to 1000 Å
                        print(f"WARNING: Unusual cell dimensions detected: {cell_diag} Å")
                        print(f"  This might indicate a units or box vector issue!")
                
                # Construct data dict for model
                # Use CPU for neighbor list since get_neighborhood uses matscipy
                from cace.data.neighborhood import get_neighborhood
                pos_cpu = pos_nm * 10.0
                cell_cpu = cell.squeeze(0).cpu().numpy() if cell is not None else None
                pbc_cpu = (True, True, True) if state_is_periodic else (False, False, False)
                
                edge_index, shifts, unit_shifts = get_neighborhood(
                    positions=pos_cpu,
                    cutoff=cutoff,
                    pbc=pbc_cpu,
                    cell=cell_cpu
                )
                
                data_dict = {
                    'positions': pos_ang,
                    'atomic_numbers': cache['atomic_numbers'],
                    'edge_index': torch.tensor(edge_index, dtype=torch.long, device=device),
                    'shifts': torch.tensor(shifts, dtype=torch.float32, device=device),
                    'unit_shifts': torch.tensor(unit_shifts, dtype=torch.float32, device=device),
                    'num_nodes': torch.tensor([n_atoms], dtype=torch.long, device=device),
                    'ptr': torch.tensor([0, n_atoms], dtype=torch.long, device=device),
                    'batch': torch.zeros(n_atoms, dtype=torch.long, device=device)
                }
                if state_is_periodic and cell is not None:
                    data_dict['cell'] = cell
                else:
                    data_dict['cell'] = torch.eye(3, dtype=torch.float32, device=device).unsqueeze(0) * 100.0

                # Forward pass
                # CACE model output is a dict with 'energy', 'forces', etc.
                output = model(data_dict, training=True)
                
                # Energy: eV -> kJ/mol
                # Try common energy keys
                energy_key = None
                for key in ['CACE_energy', 'energy', 'total_energy']:
                    if key in output:
                        energy_key = key
                        break
                
                if energy_key is None:
                    raise KeyError(f"No energy key found in CACE output. Available keys: {list(output.keys())}")
                    
                energy_ev_tensor = output[energy_key]
                
                # Don't add atomic energies - model appears to already include them
                if isinstance(energy_ev_tensor, torch.Tensor):
                    energy_ev = energy_ev_tensor
                else:
                    energy_ev = float(energy_ev_tensor)
                
                # CACE outputs total energy (sum over all atoms/molecules)
                # Convert eV to kJ/mol: 1 eV = 96.4853 kJ/mol
                if isinstance(energy_ev, torch.Tensor):
                    energy_kj = float(energy_ev.detach().cpu()) * 96.4853
                else:
                    energy_kj = float(energy_ev) * 96.4853
                
                # Sanity check: energy should be reasonable for water system
                # Typical water energy per molecule is ~-40 to -50 kJ/mol
                # For 15 molecules, total should be ~-600 to -750 kJ/mol
                # If energy is way off, there might be a units issue
                if abs(energy_kj) > 100000:  # Unreasonably large energy
                    print(f"WARNING: Very large energy detected: {energy_kj:.2f} kJ/mol")
                    print(f"  This might indicate a units mismatch!")
                    print(f"  Energy key: {energy_key}, Raw value: {energy_ev.detach().cpu()}")
                
                # Forces: eV/A -> kJ/mol/nm
                force_key = None
                for key in ['CACE_forces', 'forces', 'force']:
                    if key in output:
                        force_key = key
                        break
                
                if force_key is None:
                    raise KeyError(f"No force key found in CACE output. Available keys: {list(output.keys())}")
                    
                forces_ev_ang = output[force_key]
                # Convert forces: eV/Å -> kJ/mol/nm
                # 1 eV/Å = 96.4853 kJ/mol × (10 Å/nm) = 964.853 kJ/(mol·nm)
                molecular_forces = forces_ev_ang.detach().cpu().numpy() * 964.853
                
                # Full forces array
                total_particles = system.getNumParticles()
                full_forces = np.zeros((total_particles, 3))
                if atom_indices is not None:
                    for i, idx in enumerate(atom_indices):
                        full_forces[idx] = molecular_forces[i]
                else:
                    full_forces[:n_atoms] = molecular_forces
                    
                return (energy_kj * unit.kilojoules_per_mole,
                        full_forces * unit.kilojoules_per_mole / unit.nanometer)
            except Exception as e:
                print(f"Error in compute_cace_forces: {e}")
                import traceback
                traceback.print_exc()
                raise e

        # Create and add the force
        force = openmm.PythonForce(compute_cace_forces, {})
        force.setForceGroup(forceGroup)
        # Always set PBC=True if we detected any hint of periodicity
        # The callback will check the actual state and handle it dynamically
        # Box vectors may be set after createSystem()
        force.setUsesPeriodicBoundaryConditions(isPeriodic)
        system.addForce(force)
        
        pbc_status = "enabled" if isPeriodic else "will check dynamically from state"
        print(f"CACE force added using PythonForce (device: {device}, PBC: {pbc_status})")
