"""
umapotential_pythonforce_batch.py: UMA potential with batched RPMD support

Implements both single and batched force evaluation for UMA models.
The batched version evaluates all RPMD beads in a single ML inference call.

The model is lazy-loaded on first force computation (when OpenMM has pushed its
CUDA context), allowing PyTorch and OpenMM to share a single GPU without
CUDA_ERROR_ILLEGAL_ADDRESS. Pass device='cuda' for GPU; device='cpu' or None
for CPU.

Profiling: set OPENMMML_DEBUG_LOGS=1 to enable per-callback timing breakdown
(prep, atomicdata, inference) for bottleneck analysis.

Memory: for many RPMD beads, set OPENMMML_UMA_RPMD_CHUNK to a positive integer
(e.g. 4 or 8) to evaluate UMA on bead subsets and cap peak GPU memory.
"""
import os
import openmm
from openmm import unit
from openmmml.mlpotential import MLPotential, MLPotentialImpl, MLPotentialImplFactory
from typing import Iterable, Optional, Union
import numpy as np
import torch

# Debug controls (leave code paths intact, just silence output)
# Set OPENMMML_DEBUG_LOGS=1 to enable profiling timings for bottleneck analysis
DEBUG_LOGS = os.environ.get('OPENMMML_DEBUG_LOGS', '0').lower() in ('1', 'true', 'yes')
DEBUG_FILE_LOGS = False

# Unit conversion: eV/Å -> kJ/(mol·nm) for forces
# 1 eV/Å = 96.4853 kJ/mol / 0.1 nm = 964.853 kJ/(mol·nm)  (1 Å = 0.1 nm)
EV_ANG_TO_KJ_NM = 96.4853 * 10.0


def _wrap_water_molecules_per_bead(pos_angstrom: np.ndarray, cell_3x3: np.ndarray) -> np.ndarray:
    """
    Molecular PBC wrap for one full atomic configuration (one RPMD bead or classical step).

    **RPMD:** Call with **that bead’s** Cartesian positions and **that bead’s** cell only.
    The batched path does ``for bead_idx: wrap(all_pos[bead_idx], all_boxes[bead_idx])`` —
    no mixing across beads. Algorithm per H2O: **O** into primary cell (fractional floor);
    each **H = O_new + minimum-image(Δ)** from the **same bead’s** O/H (not centroid across
    beads; not per-atom independent wrap).

    Per-atom ``wrap_positions()`` breaks O–H when OpenMM unwraps coordinates.
    """
    n = int(pos_angstrom.shape[0])
    if n % 3 != 0:
        return np.asarray(pos_angstrom, dtype=np.float64)
    cell = np.asarray(cell_3x3, dtype=np.float64)
    inv = np.linalg.inv(cell)
    out = np.array(pos_angstrom, dtype=np.float64, copy=True)
    for m in range(n // 3):
        o = out[3 * m]
        fo = o @ inv.T
        fo = fo - np.floor(fo)
        o_n = fo[0] * cell[0] + fo[1] * cell[1] + fo[2] * cell[2]
        out[3 * m] = o_n
        for k in (1, 2):
            d = out[3 * m + k] - o
            df = d @ inv.T
            df = df - np.round(df)
            out[3 * m + k] = o_n + df[0] * cell[0] + df[1] * cell[1] + df[2] * cell[2]
    return out


def _is_water_ohm_topology(symbols: list) -> bool:
    if len(symbols) < 3 or len(symbols) % 3 != 0:
        return False
    for m in range(len(symbols) // 3):
        if symbols[3 * m] != "O" or symbols[3 * m + 1] != "H" or symbols[3 * m + 2] != "H":
            return False
    return True


def _fast_inference_settings():
    """InferenceSettings with checkpoint disabled and tf32 enabled for speed.

    torch.compile is not enabled because UMA equivariant layers have known
    incompatibilities (shape mismatches in MoLE layers).
    """
    from dataclasses import replace
    from fairchem.core.units.mlip_unit.api.inference import inference_settings_default
    d = inference_settings_default()
    return replace(d, activation_checkpointing=False, compile=False, tf32=True)


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
        use_atom_wrap_for_lammps_parity: bool = False,
        conservative_forces: bool = False,
        **args,
    ) -> None:
        """Add UMA force with batched RPMD support."""
        import torch

        if inference_settings == "optimized":
            inference_settings = _fast_inference_settings()
        
        try:
            from fairchem.core.calculate import pretrained_mlip
            from fairchem.core.datasets.atomic_data import AtomicData
        except ImportError as e:
            raise ImportError(f"Failed to import fairchem: {e}")

        try:
            from ase import Atoms
            from ase.geometry import wrap_positions
        except ImportError as e:
            raise ImportError(f"Failed to import ase: {e}")

        # Device: default to CPU if not specified. Lazy-load on first compute allows
        # device='cuda' to work (model loads when OpenMM's context is pushed).
        if device is None:
            device = "cpu"

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

        # Pre-allocate buffers (no device dependency)
        total_particles = system.getNumParticles()
        full_forces_buffer = np.zeros((total_particles, 3), dtype=np.float32) if total_particles > n_atoms or atom_indices is not None else None
        
        # Cache for closures. Model and valid_dataset_name loaded lazily on first compute.
        cache = {
            'predict_unit': None,
            'valid_dataset_name': None,
            'full_forces_buffer': full_forces_buffer,
            'batch_pos_buffer': None,
            'batch_box_buffer': None,
            'max_beads': 0,
            'ase_atoms_templates': [],
            'atomic_data_templates': [],  # Reused across steps; only pos/cell updated
        }

        # Config for lazy load (stored in closure scope)
        _use_atom_wrap = use_atom_wrap_for_lammps_parity
        _config = {
            'model_name': self.model_name,
            'inference_settings': inference_settings,
            'task_name': task_name,
            'device': device,
            'symbols': symbols,
            'charge': charge,
            'spin': spin,
            'isPeriodic': isPeriodic,
            'conservative_forces': conservative_forces,
        }

        def _ensure_model_loaded():
            """Load UMA model on first compute, when OpenMM's CUDA context is current."""
            if cache['predict_unit'] is not None:
                return
            dev = _config['device']
            if dev == 'cuda':
                torch.backends.cudnn.benchmark = True
                torch.backends.cuda.matmul.allow_tf32 = True
                torch.backends.cudnn.allow_tf32 = True
                # Workaround for gradual simulation slowdown (NVFuser fallback). Set
                # PYTORCH_NVFUSER_DISABLE=fallback or torch._C._jit_set_nvfuser_enabled(False).
                try:
                    torch._C._jit_set_nvfuser_enabled(False)
                except Exception:
                    pass
            predict_unit = pretrained_mlip.get_predict_unit(
                _config['model_name'],
                inference_settings=_config['inference_settings'],
                device=_config['device'],
            )
            task = _config['task_name']
            if task is None:
                valid = list(predict_unit.dataset_to_tasks.keys())
                if len(valid) == 1:
                    task = valid[0]
                else:
                    raise ValueError(f"Multiple datasets: {valid}. Specify task_name.")
            elif task not in predict_unit.dataset_to_tasks:
                raise ValueError(f"Invalid task_name '{task}'. Valid: {list(predict_unit.dataset_to_tasks.keys())}")
            if _config.get('conservative_forces', False):
                _backbone = getattr(predict_unit.model.module, 'backbone', predict_unit.model.module)
                if hasattr(_backbone, 'regress_config'):
                    _was_direct = _backbone.regress_config.direct_forces
                    _backbone.regress_config.direct_forces = False
                    print(
                        f"  OpenMM UMA: conservative_forces=True — "
                        f"backbone.regress_config.direct_forces {_was_direct} → False",
                        flush=True,
                    )
                else:
                    print(
                        "  OpenMM UMA: WARNING conservative_forces set but regress_config missing",
                        flush=True,
                    )
            cache['predict_unit'] = predict_unit
            cache['valid_dataset_name'] = task
            try:
                df = bool(predict_unit.direct_forces)
            except Exception as exc:
                print(
                    f"Lazy-loaded UMA model '{_config['model_name']}' on {_config['device']} "
                    f"(OpenMM context); could not read predict_unit.direct_forces: {exc}",
                    flush=True,
                )
            else:
                if df:
                    mode_desc = (
                        "learned direct force head; FairChem predict() uses torch.no_grad()"
                    )
                else:
                    mode_desc = "autograd forces from model energy (grad-enabled inference path)"
                print(
                    f"Lazy-loaded UMA model '{_config['model_name']}' on {_config['device']} "
                    f"(OpenMM context).\n"
                    f"  predict_unit.direct_forces={df} — {mode_desc}.",
                    flush=True,
                )
        
        # Import modules once for closure scope
        from fairchem.core.datasets.atomic_data import atomicdata_list_to_batch

        # Single-copy computation function
        def compute_uma_forces_single(state):
            """Compute forces for a single copy - use shared CUDA context with OpenMM."""
            try:
                _ensure_model_loaded()
                predict_unit = cache['predict_unit']
                valid_dataset_name = cache['valid_dataset_name']
                device = _config['device']

                all_pos_nm = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                if atom_indices is not None:
                    pos_nm = all_pos_nm[atom_indices].copy()
                else:
                    pos_nm = all_pos_nm[:n_atoms].copy()

                pos_angstrom = np.ascontiguousarray(pos_nm * 10.0, dtype=np.float64)

                # Detect periodicity from the runtime state, not from the
                # closure captured at addForces time, because box vectors may
                # be set on the Context after system creation.
                box_vectors = state.getPeriodicBoxVectors(asNumpy=True)
                has_periodic_box = box_vectors is not None and np.any(
                    box_vectors.value_in_unit(unit.nanometer)
                )

                cell_angstrom = None
                if has_periodic_box:
                    cell_angstrom = np.ascontiguousarray(
                        box_vectors.value_in_unit(unit.nanometer) * 10.0,
                        dtype=np.float64,
                    )
                    if _use_atom_wrap or not _is_water_ohm_topology(symbols):
                        pos_angstrom = wrap_positions(
                            pos_angstrom, cell=cell_angstrom, pbc=True, eps=0
                        )
                    else:
                        pos_angstrom = _wrap_water_molecules_per_bead(
                            pos_angstrom, cell_angstrom
                        )

                atoms_ase = Atoms(
                    symbols=symbols,
                    positions=pos_angstrom,
                    pbc=has_periodic_box,
                )
                atoms_ase.info['charge'] = charge
                atoms_ase.info['spin'] = spin
                if cell_angstrom is not None:
                    atoms_ase.set_cell(cell_angstrom)

                data = AtomicData.from_ase(
                    atoms_ase,
                    task_name=valid_dataset_name,
                    r_edges=False,
                    r_data_keys=['spin', 'charge'],
                )

                pred = predict_unit.predict(data)

                energy_ev = float(pred['energy'].cpu().item())
                forces_ev_ang = pred['forces'].cpu().numpy()

                del pred, data

                energy_kj = energy_ev * 96.4853
                molecular_forces = forces_ev_ang * EV_ANG_TO_KJ_NM

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
                _ensure_model_loaded()
                predict_unit = cache['predict_unit']
                valid_dataset_name = cache['valid_dataset_name']
                device = _config['device']

                import time
                batch_start = time.time()
                
                num_copies = len(states)
                
                if cache['batch_pos_buffer'] is None or cache['max_beads'] < num_copies:
                    cache['batch_pos_buffer'] = np.zeros((num_copies, n_atoms, 3), dtype=np.float64)
                    cache['batch_box_buffer'] = np.zeros((num_copies, 3, 3), dtype=np.float64)
                    cache['max_beads'] = num_copies
                
                all_pos_nm = cache['batch_pos_buffer'][:num_copies]
                all_boxes = cache['batch_box_buffer'][:num_copies]
                has_box = False
                
                for bead_idx, bead_state in enumerate(states):
                    pos = bead_state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                    if atom_indices is not None:
                        all_pos_nm[bead_idx] = pos[atom_indices]
                    else:
                        all_pos_nm[bead_idx] = pos[:n_atoms]
                    
                    if np.any(np.isnan(all_pos_nm[bead_idx])) or np.any(np.isinf(all_pos_nm[bead_idx])):
                        raise ValueError(f"Invalid positions (NaN/Inf) detected in bead {bead_idx}")
                    
                    box_vectors = bead_state.getPeriodicBoxVectors(asNumpy=True)
                    if box_vectors is not None:
                        box_nm = box_vectors.value_in_unit(unit.nanometer)
                        if np.any(box_nm):
                            all_boxes[bead_idx] = box_nm
                            has_box = True
                
                all_pos_angstrom = all_pos_nm * 10.0

                if has_box:
                    for i in range(num_copies):
                        cell_angstrom = np.ascontiguousarray(all_boxes[i] * 10.0, dtype=np.float64)
                        if _use_atom_wrap or not _is_water_ohm_topology(symbols):
                            all_pos_angstrom[i] = wrap_positions(
                                all_pos_angstrom[i],
                                cell=cell_angstrom,
                                pbc=True,
                                eps=0,
                            )
                        else:
                            all_pos_angstrom[i] = _wrap_water_molecules_per_bead(
                                all_pos_angstrom[i], cell_angstrom
                            )
                
                t_prep = time.time() - batch_start
                t_atomic_start = time.time()

                if len(cache['ase_atoms_templates']) < num_copies:
                    cache['ase_atoms_templates'] = []
                    for i in range(num_copies):
                        atoms_ase = Atoms(symbols=symbols, pbc=has_box)
                        atoms_ase.info['charge'] = charge
                        atoms_ase.info['spin'] = spin
                        cache['ase_atoms_templates'].append(atoms_ase)

                atomic_templates = cache['atomic_data_templates']
                if len(atomic_templates) < num_copies:
                    for i in range(len(atomic_templates), num_copies):
                        atoms_ase = cache['ase_atoms_templates'][i]
                        atoms_ase.set_positions(all_pos_angstrom[i])
                        if has_box:
                            atoms_ase.set_cell(all_boxes[i] * 10.0)
                        data = AtomicData.from_ase(
                            atoms_ase,
                            task_name=valid_dataset_name,
                            r_edges=False,
                            r_data_keys=['spin', 'charge'],
                        )
                        data.sid = [f"bead-{i}"]
                        atomic_templates.append(data)

                for i in range(num_copies):
                    atomic_templates[i].pos = torch.from_numpy(
                        np.ascontiguousarray(all_pos_angstrom[i], dtype=np.float32)
                    ).float()
                    if has_box:
                        cell_np = np.ascontiguousarray(all_boxes[i] * 10.0, dtype=np.float32)
                        atomic_templates[i].cell = torch.from_numpy(cell_np).float().unsqueeze(0)
                    co = getattr(atomic_templates[i], "cell_offsets", None)
                    if co is not None and co.numel() > 0:
                        atomic_templates[i].cell_offsets = co.float()

                data_list = atomic_templates[:num_copies]
                t_atomicdata = time.time() - t_atomic_start
                t_inference_start = time.time()

                # Optional: split bead batch to cap peak GPU memory (32+ beads on consumer GPUs).
                # Example: OPENMMML_UMA_RPMD_CHUNK=8
                _chunk_s = os.environ.get("OPENMMML_UMA_RPMD_CHUNK", "").strip()
                chunk_max = int(_chunk_s) if _chunk_s.isdigit() and int(_chunk_s) > 0 else None
                if chunk_max is not None and num_copies > chunk_max:
                    n_chunks = (num_copies + chunk_max - 1) // chunk_max
                else:
                    n_chunks = 1
                    chunk_max = num_copies

                total_energy = 0.0
                forces_kj_nm_list: list = []

                for chunk_idx in range(n_chunks):
                    start_b = chunk_idx * chunk_max
                    end_b = min(start_b + chunk_max, num_copies)
                    sub_list = data_list[start_b:end_b]
                    n_sub = end_b - start_b

                    batch_data = atomicdata_list_to_batch(sub_list)
                    batch_indices = batch_data.batch.numpy()

                    pred = predict_unit.predict(batch_data)

                    energies_t = pred['energy'].detach().to('cpu', non_blocking=True)
                    forces_t = pred['forces'].detach().to('cpu', non_blocking=True)
                    energies_ev = energies_t.numpy()
                    forces_ev_ang = forces_t.numpy()

                    energies_kj = energies_ev * 96.4853
                    total_energy += float(np.sum(energies_kj))

                    forces_ev_ang_list = [forces_ev_ang[batch_indices == i] for i in range(n_sub)]

                    for i, forces_ev_ang_i in enumerate(forces_ev_ang_list):
                        forces_kj_nm = forces_ev_ang_i * EV_ANG_TO_KJ_NM

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

                t_inference = time.time() - t_inference_start

                forces_cpu = np.array(forces_kj_nm_list, dtype=np.float32)

                return (total_energy * unit.kilojoules_per_mole,
                        forces_cpu * unit.kilojoules_per_mole / unit.nanometer)
            except Exception as e:
                if _from_single:
                    raise openmm.OpenMMException(f"UMA batch computation failed (from single-copy): {e}")
                raise openmm.OpenMMException(f"UMA batch computation failed: {e}")

        # Create PythonForce with both single and batch callbacks
        force = openmm.PythonForce(compute_uma_forces_single, {}, compute_uma_forces_batch)
        force.setForceGroup(forceGroup)
        force.setUsesPeriodicBoundaryConditions(True)
        system.addForce(force)

        wrap_mode = "atom (LAMMPS parity)" if _use_atom_wrap else "molecular (water O-H preserving)"
        print(
            f"UMA force added with batched RPMD support (model: {self.model_name}, "
            f"task: {task_name or 'auto'}, device: {device}, wrap: {wrap_mode}, "
            f"conservative_forces CLI={conservative_forces}). "
            f"Effective predict_unit.direct_forces is printed on first force evaluation.",
            flush=True,
        )


# No need to register at module level - MLPotential will auto-register from entry points
