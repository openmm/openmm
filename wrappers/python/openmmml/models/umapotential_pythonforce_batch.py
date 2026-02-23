"""
umapotential_pythonforce_batch.py: UMA potential with batched RPMD support

Implements both single and batched force evaluation for UMA models.
The batched version evaluates all RPMD beads in a single ML inference call.
"""
import openmm
from openmm import unit
from openmmml.mlpotential import MLPotential, MLPotentialImpl, MLPotentialImplFactory
from typing import Iterable, Optional, Union
import numpy as np
import torch

# Debug controls (leave code paths intact, just silence output)
DEBUG_LOGS = False
DEBUG_FILE_LOGS = False


class UMAPotentialPythonForceBatchedImplFactory(MLPotentialImplFactory):
    """Factory that creates UMAPotentialPythonForceBatchedImpl objects."""

    def createImpl(self, name: str, **args) -> MLPotentialImpl:
        return UMAPotentialPythonForceBatchedImpl(name)


class UMAPotentialPythonForceBatchedImpl(MLPotentialImpl):
    """UMA potential with RPMD batching support."""

    def __init__(self, name: str) -> None:
        self.name = name
        # Extract model name: uma-s-1p1-pythonforce-batch -> uma-s-1p1
        if name.endswith('-pythonforce-batch'):
            self.model_name = name.replace('-pythonforce-batch', '')
        elif name.endswith('-batch'):
            self.model_name = name.replace('-batch', '')
        else:
            self.model_name = name

    def addForces(
        self,
        topology: openmm.app.Topology,
        system: openmm.System,
        atoms: Optional[Iterable[int]],
        forceGroup: int,
        task_name: Optional[str] = None,
        inference_settings: Union[str, object] = "default",
        charge: int = 0,
        spin: int = 1,
        device: Optional[str] = None,
        **args,
    ) -> None:
        """Add UMA force with batched RPMD support."""
        import torch
        
        try:
            from fairchem.core.calculate import pretrained_mlip
            from fairchem.core.datasets.atomic_data import AtomicData
        except ImportError as e:
            raise ImportError(f"Failed to import fairchem: {e}")

        try:
            from ase import Atoms
        except ImportError as e:
            raise ImportError(f"Failed to import ase: {e}")

        # Load UMA model
        if not hasattr(self, '_predict_unit'):
            if device is None:
                device = "cuda" if torch.cuda.is_available() else "cpu"

            print(f"Loading UMA model '{self.model_name}' on {device}...")
            
            # Enable PyTorch optimizations for faster inference
            torch.backends.cudnn.benchmark = True  # Auto-tune kernel selection
            torch.backends.cuda.matmul.allow_tf32 = True  # Enable TF32 for matmul
            torch.backends.cudnn.allow_tf32 = True  # Enable TF32 for cudnn
            
            self._predict_unit = pretrained_mlip.get_predict_unit(
                self.model_name,
                inference_settings=inference_settings,
                device=device,
            )
            self._device = device
            
            # Note: torch.compile() optimization attempted but causes compilation issues  
            # with UMA's rotation layers. Keeping eager mode for stability.
            
            print(f"PyTorch optimizations enabled:")
            print(f"  - cudnn.benchmark: {torch.backends.cudnn.benchmark}")
            print(f"  - TF32 (matmul): {torch.backends.cuda.matmul.allow_tf32}")
            print(f"  - TF32 (cudnn): {torch.backends.cudnn.allow_tf32}")

        predict_unit = self._predict_unit

        # Validate task name
        if task_name is None:
            valid_datasets = list(predict_unit.dataset_to_tasks.keys())
            if len(valid_datasets) == 1:
                task_name = valid_datasets[0]
            else:
                raise ValueError(f"Multiple tasks available: {valid_datasets}. Specify task_name.")
        
        # Store the valid dataset name for the model
        valid_dataset_name = task_name  # Dataset ID model expects

        # Extract atomic symbols from topology
        symbols = []
        for atom in topology.atoms():
            symbols.append(atom.element.symbol)

        n_atoms = len(symbols)
        atom_indices = atoms if atoms is not None else None

        # Check for periodicity
        isPeriodic = (
            topology.getPeriodicBoxVectors() is not None
        ) or system.usesPeriodicBoundaryConditions()

        # Create initial AtomicData template
        n_atoms = len(symbols)
        initial_positions = np.zeros((n_atoms, 3), dtype=np.float32)
        
        atoms_ase = Atoms(
            symbols=symbols,
            positions=initial_positions,
            pbc=isPeriodic
        )
        atoms_ase.info['charge'] = charge
        atoms_ase.info['spin'] = spin
        
        atomic_data_template = AtomicData.from_ase(
            atoms_ase,
            task_name=task_name,
            r_edges=False,
            r_data_keys=['spin', 'charge']
        )
        atomic_data_template = atomic_data_template.to(device)
        
        # Fix dataset field
        if hasattr(atomic_data_template, 'dataset'):
            if isinstance(atomic_data_template.dataset, str):
                atomic_data_template.dataset = [atomic_data_template.dataset]
        
        # Pre-allocate buffers
        total_particles = system.getNumParticles()
        full_forces_buffer = np.zeros((total_particles, 3), dtype=np.float32) if total_particles > n_atoms or atom_indices is not None else None
        
        # Cache for closures
        cache = {
            'atomic_data': atomic_data_template,
            'full_forces_buffer': full_forces_buffer,
            'batch_pos_buffer': None,  # Pre-allocated on first use
            'batch_box_buffer': None,  # Pre-allocated on first use
            'max_beads': 0,  # Track max beads seen for buffer sizing
            'ase_atoms_templates': []  # Cache ASE Atoms objects (created on first use)
        }
        
        # Import modules once for closure scope
        from fairchem.core.datasets.atomic_data import atomicdata_list_to_batch

        # Single-copy computation function
        def compute_uma_forces_single(state):
            """Compute forces for a single copy - use shared CUDA context with OpenMM."""
            try:
                if DEBUG_LOGS and not cache.get("warned_single", False):
                    print(
                        "Single-copy function called (likely from getState/reporters). "
                        "Batch path is still used for RPMD force evaluation.",
                        flush=True,
                    )
                    cache["warned_single"] = True
                
                # Get positions from OpenMM State and ensure it's contiguous numpy with no shared memory
                all_pos_nm = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                if atom_indices is not None:
                    pos_nm = all_pos_nm[atom_indices].copy()  # Explicit copy
                else:
                    pos_nm = all_pos_nm[:n_atoms].copy()  # Explicit copy

                pos_angstrom = np.ascontiguousarray(pos_nm * 10.0)  # Ensure contiguous
                
                # Create ASE atoms (this is CPU-side)
                atoms_ase = Atoms(symbols=symbols, positions=pos_angstrom, pbc=isPeriodic)
                atoms_ase.info['charge'] = charge
                atoms_ase.info['spin'] = spin

                # Create AtomicData (initially CPU-side)
                data = AtomicData.from_ase(atoms_ase, task_name=task_name, r_edges=False, r_data_keys=['spin', 'charge'])

                
                # Set cell before moving to GPU (matches batch path)
                # Now using shared CUDA context, so tensors should work correctly
                with torch.no_grad():
                    if isPeriodic:
                        try:
                            box = state.getPeriodicBoxVectors(asNumpy=True)
                            if box is not None:
                                # Explicit copy to break any OpenMM memory connection
                                cell_np = np.ascontiguousarray(box.value_in_unit(unit.nanometer) * 10.0)
                                
                                                
                                # Allocate cell in PyTorch CUDA context (not OpenMM's)
                                # Create empty tensor on GPU first (PyTorch allocates in its own context)
                                cell_gpu = torch.empty((1, 3, 3), dtype=torch.float32, device=device)
                                # Copy data from numpy (CPU) to GPU tensor
                                cell_gpu[0] = torch.from_numpy(cell_np).float()
                                
                                                
                                # Set cell BEFORE .to(device) - matches batch path exactly (line 378)
                                data.cell = cell_gpu
                                
                                                
                                del cell_gpu
                        except Exception as box_err:
                                        print(f"Warning: Box vector error in single-copy: {box_err}")
                    
                    # Now move entire data to GPU (matches batch path line 382)
                    data_device = data.to(device)
                    
                        
                    # dataset must be list (FAIRChem iterates over it)
                    data_device.dataset = [valid_dataset_name]
                    
                        
                    try:
                        pred = predict_unit.predict(data_device)
                    except Exception as pred_err:
                                raise
                    
                    # Extract results and move to CPU immediately
                    energy_ev = float(pred['energy'].cpu().item())
                    forces_ev_ang = pred['forces'].cpu().numpy()

                # Clean up GPU tensors immediately
                del data_device, pred, data
                torch.cuda.empty_cache()

                # Convert units (all CPU now)
                energy_kj = energy_ev * 96.4853
                molecular_forces = forces_ev_ang * (96.4853 / 10.0)

                if cache['full_forces_buffer'] is not None:
                    forces_kj_nm = cache['full_forces_buffer'].copy()
                    forces_kj_nm.fill(0.0)
                    if atom_indices is not None:
                        for i, atom_idx in enumerate(atom_indices):
                            forces_kj_nm[atom_idx] = molecular_forces[i]
                    else:
                        forces_kj_nm[:n_atoms] = molecular_forces
                else:
                    forces_kj_nm = molecular_forces

                return (energy_kj * unit.kilojoules_per_mole,
                        forces_kj_nm * unit.kilojoules_per_mole / unit.nanometer)
            except Exception as e:
                import traceback
                error_msg = f"UMA force computation failed: {type(e).__name__}: {str(e)}\n{traceback.format_exc()}"
                raise openmm.OpenMMException(error_msg)

        # Batched computation function for RPMD - TRUE tensor batching
        def compute_uma_forces_batch(states, _from_single=False):
            """Compute forces for all RPMD beads using TRUE tensor batching."""
            try:
                import time
                batch_start = time.time()
                if DEBUG_LOGS:
                    print(f"Batch function called: num_beads={len(states)} _from_single={_from_single}", flush=True)
                
                num_copies = len(states)
                
                # Pre-allocate or reuse position and box buffers
                if cache['batch_pos_buffer'] is None or cache['max_beads'] < num_copies:
                    cache['batch_pos_buffer'] = np.zeros((num_copies, n_atoms, 3), dtype=np.float64)
                    if isPeriodic:
                        cache['batch_box_buffer'] = np.zeros((num_copies, 3, 3), dtype=np.float64)
                    cache['max_beads'] = num_copies
                
                all_pos_nm = cache['batch_pos_buffer'][:num_copies]
                all_boxes = cache['batch_box_buffer'][:num_copies] if isPeriodic else None
                
                # =================================================================
                # TRUE TENSOR BATCHING: Build a batched AtomicData object
                # =================================================================
                
                # Collect all positions and boxes (write directly into pre-allocated buffers)
                for bead_idx, state in enumerate(states):
                    # Positions are already wrapped by C++ when State is built in executeBatch
                    pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                    if atom_indices is not None:
                        all_pos_nm[bead_idx] = pos[atom_indices]
                    else:
                        all_pos_nm[bead_idx] = pos[:n_atoms]
                    
                    # Validate positions: check for NaN, Inf, or extreme values
                    if np.any(np.isnan(all_pos_nm[bead_idx])) or np.any(np.isinf(all_pos_nm[bead_idx])):
                        raise ValueError(f"Invalid positions (NaN/Inf) detected in bead {bead_idx}")
                    if np.any(np.abs(all_pos_nm[bead_idx]) > 1000.0):
                        raise ValueError(f"Extreme positions detected in bead {bead_idx}: max={np.max(np.abs(all_pos_nm[bead_idx])):.2f} nm")
                    
                    if isPeriodic:
                        try:
                            box = state.getPeriodicBoxVectors(asNumpy=True)
                            if box is not None:
                                all_boxes[bead_idx] = box.value_in_unit(unit.nanometer)
                                # Validate box vectors
                                if np.any(np.isnan(all_boxes[bead_idx])) or np.any(np.isinf(all_boxes[bead_idx])):
                                    raise ValueError(f"Invalid box vectors (NaN/Inf) detected in bead {bead_idx}")
                                # Ensure box is reasonable (not zero or negative)
                                box_diag = np.array([all_boxes[bead_idx][0,0], all_boxes[bead_idx][1,1], all_boxes[bead_idx][2,2]])
                                if np.any(box_diag <= 0) or np.any(box_diag > 1000.0):
                                    raise ValueError(f"Invalid box size in bead {bead_idx}: {box_diag}")
                        except Exception as e:
                            raise ValueError(f"Error getting box vectors for bead {bead_idx}: {e}")
                
                # Convert to Angstroms (views into buffer, no copy)
                all_pos_angstrom = all_pos_nm * 10.0
                
                # PROFILING
                t_prep = time.time() - batch_start
                
                # =================================================================
                # BATCHED EXECUTION: Single batched model call
                # =================================================================
                t_inference_start = time.time()
                
                # Create or reuse ASE Atoms templates
                if len(cache['ase_atoms_templates']) < num_copies:
                    # First time or increased number of beads - create new templates
                    cache['ase_atoms_templates'] = []
                    for i in range(num_copies):
                        if isPeriodic:
                            atoms_ase = Atoms(symbols=symbols, pbc=True)
                        else:
                            atoms_ase = Atoms(symbols=symbols, pbc=False)
                        atoms_ase.info['charge'] = charge
                        atoms_ase.info['spin'] = spin
                        cache['ase_atoms_templates'].append(atoms_ase)
                
                # Update positions (and cell if periodic) in cached templates
                data_list = []
                for i in range(num_copies):
                    atoms_ase = cache['ase_atoms_templates'][i]
                    atoms_ase.set_positions(all_pos_angstrom[i])
                    if isPeriodic and all_boxes is not None:
                        atoms_ase.set_cell(all_boxes[i] * 10.0)
                    
                    data = AtomicData.from_ase(
                        atoms_ase,
                        task_name=task_name,
                        r_edges=False,
                        r_data_keys=['spin', 'charge'],
                    )
                    data.sid = [f"bead-{i}"]
                    data_list.append(data)
                
                batch_data = atomicdata_list_to_batch(data_list)
                batch_indices = batch_data.batch.numpy()
                batch_data = batch_data.to(device)
                
                with torch.no_grad():
                    try:
                        pred = predict_unit.predict(batch_data)
                    except RuntimeError as e:
                        # If CUDA error occurs, provide detailed diagnostics
                        if "CUDA" in str(e) or "illegal" in str(e).lower():
                            error_info = {
                                "model": self.model_name,
                                "num_beads": num_copies,
                                "num_atoms": n_atoms,
                                "batch_size": batch_data.num_graphs if hasattr(batch_data, 'num_graphs') else 'unknown',
                                "has_cell": hasattr(batch_data, 'cell'),
                                "has_edges": hasattr(batch_data, 'edge_index'),
                                "device": str(device),
                            }
                            print(f"CUDA error in batched inference: {error_info}")
                        raise
                
                energies_ev = pred['energy'].detach().cpu().numpy()
                forces_ev_ang = pred['forces'].detach().cpu().numpy()
                
                t_inference = time.time() - t_inference_start
                if DEBUG_LOGS:
                    print(f"  ⏱️ Batched GPU inference: {t_inference*1000:.1f} ms ({num_copies} systems)", flush=True)
                
                # =================================================================
                # Collect results and reshape
                # =================================================================
                energies_kj = energies_ev * 96.4853
                total_energy = float(np.sum(energies_kj))
                
                forces_ev_ang_list = [forces_ev_ang[batch_indices == i] for i in range(num_copies)]
                
                forces_kj_nm_list = []
                for i, forces_ev_ang_i in enumerate(forces_ev_ang_list):
                    forces_kj_nm = forces_ev_ang_i * (96.4853 / 10.0)
                    
                    if cache['full_forces_buffer'] is not None:
                        forces_full = cache['full_forces_buffer'].copy()
                        forces_full.fill(0.0)
                        if atom_indices is not None:
                            for j, atom_idx in enumerate(atom_indices):
                                forces_full[atom_idx] = forces_kj_nm[j]
                        else:
                            forces_full[:n_atoms] = forces_kj_nm
                        forces_kj_nm_list.append(forces_full)
                    else:
                        forces_kj_nm_list.append(forces_kj_nm)
                
                forces_cpu = np.array(forces_kj_nm_list, dtype=np.float32)
                
                batch_time = time.time() - batch_start
                if DEBUG_LOGS:
                    print(f"BATCH COMPLETED in {batch_time:.3f}s for {num_copies} beads ({batch_time/num_copies*1000:.1f} ms/bead)", flush=True)
                    print(f"   Breakdown: prep={t_prep*1000:.1f}ms, inference={t_inference*1000:.1f}ms", flush=True)
                    print(f"   Speedup vs sequential: {num_copies * 41.4 / (t_inference*1000):.2f}x (ideal: {num_copies}x)", flush=True)
                
                return (total_energy * unit.kilojoules_per_mole,
                        forces_cpu * unit.kilojoules_per_mole / unit.nanometer)
            except Exception as e:
                if _from_single:
                    raise openmm.OpenMMException(f"UMA batch computation failed (from single-copy): {e}")
                raise openmm.OpenMMException(f"UMA batch computation failed: {e}")

        # Create PythonForce with both single and batch callbacks
        force = openmm.PythonForce(compute_uma_forces_single, {}, compute_uma_forces_batch)
        force.setForceGroup(forceGroup)
        force.setUsesPeriodicBoundaryConditions(isPeriodic)
        system.addForce(force)

        print(f"UMA force added with batched RPMD support (model: {self.model_name}, task: {task_name}, device: {device})")


# No need to register at module level - MLPotential will auto-register from entry points
