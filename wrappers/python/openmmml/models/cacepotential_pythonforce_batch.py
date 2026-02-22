"""
cacepotential_pythonforce_batch.py: CACE-LR potential with batched RPMD support

Implements both single and batched force evaluation for CACE-LR models.
The batched version evaluates all RPMD beads in a single ML inference call.
"""
import openmm
from openmm import unit
from openmmml.mlpotential import MLPotential, MLPotentialImpl, MLPotentialImplFactory
from typing import Iterable, Optional
import numpy as np
import torch
from ase import Atoms

# CACE imports are done lazily inside addForces/compute_* so that openmmml
# can be imported without the cace package when using only UMA or other models.

# Debug controls
DEBUG_LOGS = False

class CACEPotentialBatchedImplFactory(MLPotentialImplFactory):
    """Factory that creates batched CACE potential implementations."""
    def createImpl(self, name: str, **args) -> MLPotentialImpl:
        return CACEPotentialBatchedImpl(name)

class CACEPotentialBatchedImpl(MLPotentialImpl):
    """
    CACE-LR potential with batched RPMD support.
    
    Features:
    - Single-copy evaluation for standard MD
    - Batched evaluation for RPMD (all beads in one inference call)
    - GPU-accelerated via PyTorch
    - Automatic periodic boundary condition handling
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
        atomic_energies: Optional[dict] = None,
        **args,
    ) -> None:
        """Add CACE force with batched RPMD support."""
        
        try:
            from cace.data.neighborhood import get_neighborhood
        except ImportError as e:
            raise ImportError(
                f"Failed to import CACE: {e}. "
                "Make sure CACE is installed and in your PYTHONPATH."
            )

        # Load the model
        if device is None:
            device = "cuda" if torch.cuda.is_available() else "cpu"
        
        print(f"Loading CACE-LR model from '{self.model_path}' on {device}...")
        
        # Check if it's a file path or we need to look in cace/water (openmm root is 4 levels up from this file)
        import os
        from pathlib import Path
        actual_path = self.model_path
        if not os.path.exists(actual_path):
            openmm_root = Path(__file__).resolve().parents[4]
            cace_water_path = openmm_root / "cace" / "water" / "fit" / "fit_version_1" / "best_model.pth"
            if cace_water_path.exists():
                actual_path = str(cace_water_path)
                print(f"  Note: Using water model from cace/water: {actual_path}")
            else:
                raise FileNotFoundError(f"CACE model file not found: {self.model_path}")

        # Load model
        model = torch.load(actual_path, map_location=device, weights_only=False)
        model.to(device)
        model.eval()
        
        # Enable PyTorch optimizations
        torch.backends.cudnn.benchmark = True
        if device == "cuda":
            torch.backends.cuda.matmul.allow_tf32 = True
            torch.backends.cudnn.allow_tf32 = True
        
        # Get cutoff from model
        try:
            cutoff = model.representation.cutoff
        except AttributeError:
            cutoff = model.models[0].representation.cutoff

        print(f"  Model cutoff: {cutoff} Å")
        
        # Get atomic information
        includedAtoms = list(topology.atoms())
        if atoms is not None:
            includedAtoms = [includedAtoms[i] for i in atoms]
        symbols = [atom.element.symbol for atom in includedAtoms]
        atomic_numbers = np.array([atom.element.atomic_number for atom in includedAtoms])
        atom_indices = atoms if atoms is not None else None
        n_atoms = len(symbols)
        
        # Get atomic energies
        if atomic_energies is None:
            # Default atomic energies for water model (from cace/water training)
            atomic_energies = {1: -5.853064337340629, 8: -2.926532168670322}
            print(f"  Using default atomic energies: H={atomic_energies[1]:.3f} eV, O={atomic_energies[8]:.3f} eV")

        # Check for periodicity
        has_default_box = False
        try:
            a, b, c = system.getDefaultPeriodicBoxVectors()
            if a.norm() > 0 or b.norm() > 0 or c.norm() > 0:
                has_default_box = True
        except:
            pass
        
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
        
        if has_default_box or has_pbc_force:
            isPeriodic = True

        # Cache for closures
        cache = {
            'device': device,
            'cutoff': cutoff,
            'n_atoms': n_atoms,
            'is_periodic': isPeriodic,
            'atom_indices': atom_indices,
            'atomic_numbers': torch.tensor(atomic_numbers, dtype=torch.long, device=device),
            'atomic_energies': atomic_energies,
            'atomic_numbers_np': atomic_numbers,
            'total_particles': system.getNumParticles(),
            'symbols': symbols,  # Needed for batched function
        }

        def compute_cace_forces_single(state):
            """Compute forces for a single copy (standard MD or single RPMD bead)."""
            try:
                if DEBUG_LOGS and not cache.get("warned_single", False):
                    print("⚠️ CACE single-copy function called", flush=True)
                    cache["warned_single"] = True
                
                # Initialize PyTorch CUDA if needed
                if torch.cuda.is_available() and not torch.cuda.is_initialized():
                    torch.cuda.init()
                
                # Get positions in nm
                all_pos_nm = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                
                # Check for PBC from box vectors
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
                
                # Prepare cell
                cell = None
                if state_is_periodic:
                    box_nm = box_vectors.value_in_unit(unit.nanometer)
                    cell_np = box_nm * 10.0  # Convert to Angstrom
                    cell = torch.tensor(cell_np, dtype=torch.float32, device=device).unsqueeze(0)
                
                # Build neighbor list on CPU
                pos_cpu = pos_nm * 10.0
                cell_cpu = cell.squeeze(0).cpu().numpy() if cell is not None else None
                pbc_cpu = (True, True, True) if state_is_periodic else (False, False, False)
                
                edge_index, shifts, unit_shifts = get_neighborhood(
                    positions=pos_cpu,
                    cutoff=cutoff,
                    pbc=pbc_cpu,
                    cell=cell_cpu
                )
                
                # Construct data dict
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

                # Forward pass - CACE needs training=True to compute forces via autograd
                # Do NOT use torch.no_grad() - CACE computes forces via autograd.grad()
                output = model(data_dict, training=True)
                
                # Extract energy (eV -> kJ/mol)
                energy_key = None
                for key in ['CACE_energy', 'energy', 'total_energy']:
                    if key in output:
                        energy_key = key
                        break
                
                if energy_key is None:
                    raise KeyError(f"No energy key found in CACE output. Available keys: {list(output.keys())}")
                    
                energy_ev_tensor = output[energy_key]
                if isinstance(energy_ev_tensor, torch.Tensor):
                    energy_ev = float(energy_ev_tensor.detach().cpu())
                else:
                    energy_ev = float(energy_ev_tensor)
                
                # Convert eV to kJ/mol
                energy_kj = energy_ev * 96.4853
                
                # Extract forces (eV/Å -> kJ/mol/nm)
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
                full_forces = np.zeros((cache['total_particles'], 3))
                if atom_indices is not None:
                    for i, idx in enumerate(atom_indices):
                        full_forces[idx] = molecular_forces[i]
                else:
                    full_forces[:n_atoms] = molecular_forces
                
                # Clean up tensors (minimal, avoid CUDA operations that can fail)
                del output, pos_ang, data_dict
                    
                return (energy_kj * unit.kilojoules_per_mole,
                        full_forces * unit.kilojoules_per_mole / unit.nanometer)
            
            except Exception as e:
                print(f"Error in compute_cace_forces_single: {e}")
                import traceback
                traceback.print_exc()
                raise e

        def compute_cace_forces_batched(all_states):
            """Compute forces for multiple RPMD beads using true tensor batching.
            
            Uses CACE's native AtomicData and Batch classes for proper graph batching.
            This enables single model forward pass for all beads simultaneously.
            """
            try:
                from cace.data import AtomicData
                from cace.tools.torch_geometric import Batch
                if DEBUG_LOGS and not cache.get("warned_batch", False):
                    print(f"✓ CACE batched function called with TRUE tensor batching (n_beads={len(all_states)})", flush=True)
                    cache["warned_batch"] = True
                
                # Initialize PyTorch CUDA if needed
                if torch.cuda.is_available() and not torch.cuda.is_initialized():
                    torch.cuda.init()
                
                n_beads = len(all_states)
                total_particles = cache['total_particles']
                symbols = cache['symbols']  # Get symbols from cache
                
                # Step 1: Create AtomicData for each bead using CACE's native infrastructure
                data_list = []
                for state in all_states:
                    # Get positions for this bead
                    all_pos_nm = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                    
                    if atom_indices is not None:
                        pos_nm = all_pos_nm[atom_indices]
                    else:
                        pos_nm = all_pos_nm[:n_atoms]
                    
                    # Convert to Angstrom
                    pos_ang = pos_nm * 10.0
                    
                    # Check for PBC from box vectors
                    box_vectors = state.getPeriodicBoxVectors(asNumpy=True)
                    state_is_periodic = box_vectors is not None and (
                        np.linalg.norm(box_vectors[0]) > 0 or 
                        np.linalg.norm(box_vectors[1]) > 0 or 
                        np.linalg.norm(box_vectors[2]) > 0
                    )
                    
                    # Create ASE Atoms object
                    if state_is_periodic:
                        box_nm = box_vectors.value_in_unit(unit.nanometer)
                        cell_ang = box_nm * 10.0  # Convert to Angstrom
                        atoms_ase = Atoms(
                            symbols=symbols,
                            positions=pos_ang,
                            pbc=True,
                            cell=cell_ang
                        )
                    else:
                        atoms_ase = Atoms(
                            symbols=symbols,
                            positions=pos_ang,
                            pbc=False
                        )
                    
                    # AtomicData.from_atoms handles neighbor list construction automatically
                    # CRITICAL: Pass atomic_energies to ensure proper energy shifts
                    data = AtomicData.from_atoms(
                        atoms_ase, 
                        cutoff=cutoff,
                        atomic_energies=cache['atomic_energies']
                    )
                    data_list.append(data)
                
                # Step 2: Batch all AtomicData objects using CACE's Batch class
                # This automatically:
                # - Concatenates positions, atomic_numbers, etc.
                # - Offsets edge_index properly for each graph
                # - Creates batch tensor [0,0,...,1,1,...,n-1,n-1,...]
                # - Creates ptr tensor [0, n0, n0+n1, ...]
                batch_data = Batch.from_data_list(data_list)
                
                # Step 3: Ensure positions require grad for force computation
                # CACE computes forces via autograd, so we need gradients
                batch_data.positions = batch_data.positions.clone().requires_grad_(True)
                
                # Step 4: Move to device (GPU)
                batch_data = batch_data.to(device)
                
                # Step 5: Convert to dict format expected by CACE model
                # Save batch indices before conversion
                batch_indices_cpu = batch_data.batch.cpu().numpy()
                
                # Convert Batch object to dict
                data_dict = {}
                for key in batch_data.keys:
                    data_dict[key] = batch_data[key]
                
                # Step 6: Single forward pass for all beads
                # CACE will use scatter_sum with batch tensor to aggregate per-graph energies
                output = model(data_dict, training=True)
                
                # Step 7: Extract energies (one per bead)
                energy_key = None
                for key in ['CACE_energy', 'energy', 'total_energy']:
                    if key in output:
                        energy_key = key
                        break
                
                if energy_key is None:
                    raise KeyError(f"No energy key found in CACE output. Available keys: {list(output.keys())}")
                
                if DEBUG_LOGS:
                    print(f"  Energy key: {energy_key}")
                    print(f"  Energy shape: {output[energy_key].shape if hasattr(output[energy_key], 'shape') else 'scalar'}")
                
                energies_ev_tensor = output[energy_key]
                if isinstance(energies_ev_tensor, torch.Tensor):
                    energies_ev = energies_ev_tensor.detach().cpu().numpy()
                else:
                    energies_ev = np.array([float(energies_ev_tensor)])
                
                if DEBUG_LOGS:
                    print(f"  Energies (eV): {energies_ev}")
                    print(f"  N beads: {n_beads}, Energies shape: {energies_ev.shape}")
                
                # Total energy is sum of all bead energies
                total_energy = np.sum(energies_ev) * 96.4853  # Convert eV to kJ/mol
                
                if DEBUG_LOGS:
                    print(f"  Total energy (kJ/mol): {total_energy}")
                
                # Step 8: Extract forces using batch indices
                force_key = None
                for key in ['CACE_forces', 'forces', 'force']:
                    if key in output:
                        force_key = key
                        break
                
                if force_key is None:
                    raise KeyError(f"No force key found in CACE output. Available keys: {list(output.keys())}")
                
                forces_ev_ang = output[force_key].detach().cpu().numpy()  # (n_beads*n_atoms, 3)
                
                # Step 9: Split forces by bead using batch indices
                all_forces = np.zeros((n_beads, total_particles, 3), dtype=np.float32)
                
                for bead_idx in range(n_beads):
                    # Select forces for this bead
                    mask = (batch_indices_cpu == bead_idx)
                    bead_forces_ev_ang = forces_ev_ang[mask]
                    
                    # Convert eV/Å -> kJ/mol/nm
                    bead_forces_kj_nm = bead_forces_ev_ang * 964.853
                    
                    # Place in full forces array
                    if atom_indices is not None:
                        for i, idx in enumerate(atom_indices):
                            all_forces[bead_idx, idx] = bead_forces_kj_nm[i]
                    else:
                        all_forces[bead_idx, :n_atoms] = bead_forces_kj_nm
                
                # Step 10: Clean up GPU resources
                # Delete tensors (minimal cleanup, avoid operations that can fail on corrupted GPU state)
                del output, batch_data, data_dict, energies_ev_tensor, data_list
                
                # Return format: (total_energy, forces_array)
                return (total_energy * unit.kilojoules_per_mole,
                        all_forces * unit.kilojoules_per_mole / unit.nanometer)
                
            except Exception as e:
                print(f"Error in compute_cace_forces_batched: {e}")
                import traceback
                traceback.print_exc()
                raise e

        # Create PythonForce with both single and batched callbacks
        # Constructor: PythonForce(single_callback, dict, batched_callback)
        force = openmm.PythonForce(compute_cace_forces_single, {}, compute_cace_forces_batched)
        force.setForceGroup(forceGroup)
        force.setUsesPeriodicBoundaryConditions(isPeriodic)
        system.addForce(force)
        
        pbc_status = "enabled" if isPeriodic else "dynamic"
        print(f"✓ CACE-LR force added with batched RPMD support (device: {device}, PBC: {pbc_status})")
        print(f"  Single-copy callback: compute_cace_forces_single")
        print(f"  Batched callback: compute_cace_forces_batched (for RPMD)")


# Register the implementation
MLPotential.registerImplFactory('cace-lr-batch', CACEPotentialBatchedImplFactory())
MLPotential.registerImplFactory('cace-pythonforce-batch', CACEPotentialBatchedImplFactory())
