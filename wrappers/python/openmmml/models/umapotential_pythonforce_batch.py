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
import ctypes
import ctypes.util

# Debug controls (leave code paths intact, just silence output)
DEBUG_LOGS = False
DEBUG_FILE_LOGS = False

# CUDA Driver API types and functions for context sharing
_cuda = None
_cuda_context_handle = None

def _init_cuda_driver_api():
    """Initialize CUDA Driver API bindings for context management."""
    global _cuda
    if _cuda is not None:
        return _cuda
    
    try:
        # Try to load libcuda.so (Linux) or cuda.dll (Windows)
        lib_path = ctypes.util.find_library('cuda')
        if lib_path is None:
            # Try common paths
            import platform
            if platform.system() == 'Linux':
                lib_path = 'libcuda.so.1'
            elif platform.system() == 'Darwin':
                lib_path = '/usr/local/cuda/lib/libcuda.dylib'
            else:
                lib_path = 'cuda.dll'
        
        _cuda = ctypes.CDLL(lib_path)
        
        # Define CUDA types - CUcontext is a void pointer
        CUcontext = ctypes.c_void_p
        CUresult = ctypes.c_int
        
        # cuCtxGetCurrent - takes pointer to CUcontext
        _cuda.cuCtxGetCurrent.argtypes = [ctypes.POINTER(CUcontext)]
        _cuda.cuCtxGetCurrent.restype = CUresult
        
        # cuCtxSetCurrent - takes CUcontext (void pointer)
        _cuda.cuCtxSetCurrent.argtypes = [CUcontext]
        _cuda.cuCtxSetCurrent.restype = CUresult
        
        # cuCtxPushCurrent
        _cuda.cuCtxPushCurrent.argtypes = [CUcontext]
        _cuda.cuCtxPushCurrent.restype = CUresult
        
        # cuCtxPopCurrent
        _cuda.cuCtxPopCurrent.argtypes = [ctypes.POINTER(CUcontext)]
        _cuda.cuCtxPopCurrent.restype = CUresult
        
        # cuInit
        _cuda.cuInit.argtypes = [ctypes.c_uint]
        _cuda.cuInit.restype = CUresult
        
        # Initialize CUDA Driver API
        result = _cuda.cuInit(0)
        if result != 0:
            raise RuntimeError(f"cuInit failed with error {result}")
        
    except Exception as e:
        # If we can't load CUDA Driver API, we'll fall back to default behavior
        try:
            if DEBUG_FILE_LOGS:
                import json
                with open('/media/extradrive/Trajectories/openmm/.cursor/debug.log', 'a') as f:
                    f.write(json.dumps({"sessionId":"debug-session","runId":"run1","hypothesisId":"CONTEXT","location":"_init_cuda_driver_api","message":"Failed to init CUDA Driver API","data":{"error":str(e)},"timestamp":int(__import__('time').time()*1000)})+'\n')
        except:
            pass
        _cuda = None
    
    return _cuda

def _ensure_openmm_cuda_context():
    """Ensure OpenMM's CUDA context is active for PyTorch operations.
    
    When called from OpenMM's PythonForce, OpenMM's CUDA context is already current.
    This function ensures PyTorch will use that same context by:
    1. Getting OpenMM's current CUDA context
    2. Initializing PyTorch CUDA (which may create its own context)
    3. Switching back to OpenMM's context and making it current
    4. PyTorch operations will then use OpenMM's shared context
    """
    global _cuda_context_handle
    
    if not torch.cuda.is_available():
        return
    
    cuda = _init_cuda_driver_api()
    if cuda is None:
        # Fall back to default PyTorch behavior
        try:
            if DEBUG_FILE_LOGS:
                import json
                with open('/media/extradrive/Trajectories/openmm/.cursor/debug.log', 'a') as f:
                    f.write(json.dumps({"sessionId":"debug-session","runId":"run1","hypothesisId":"CONTEXT","location":"_ensure_openmm_cuda_context","message":"CUDA Driver API not available, using default","data":{},"timestamp":int(__import__('time').time()*1000)})+'\n')
        except:
            pass
        torch.cuda.set_device(0)
        return
    
    try:
        # Get the current CUDA context (should be OpenMM's when called from PythonForce)
        CUcontext = ctypes.c_void_p
        openmm_ctx = CUcontext()
        result = cuda.cuCtxGetCurrent(ctypes.byref(openmm_ctx))
        
        try:
            if DEBUG_FILE_LOGS:
                import json
                with open('/media/extradrive/Trajectories/openmm/.cursor/debug.log', 'a') as f:
                    f.write(json.dumps({"sessionId":"debug-session","runId":"run1","hypothesisId":"CONTEXT","location":"_ensure_openmm_cuda_context","message":"Got CUDA context","data":{"result":result,"ctx_ptr":hex(openmm_ctx.value) if openmm_ctx.value else None},"timestamp":int(__import__('time').time()*1000)})+'\n')
        except:
            pass
        
        if result == 0 and openmm_ctx.value is not None:
            # Store OpenMM's context
            _cuda_context_handle = openmm_ctx
            
            # Push OpenMM's context onto the stack BEFORE initializing PyTorch
            # This ensures PyTorch will use OpenMM's context when it initializes
            push_result = cuda.cuCtxPushCurrent(_cuda_context_handle)
            
            try:
                if DEBUG_FILE_LOGS:
                    import json
                    with open('/media/extradrive/Trajectories/openmm/.cursor/debug.log', 'a') as f:
                        f.write(json.dumps({"sessionId":"debug-session","runId":"run1","hypothesisId":"CONTEXT","location":"_ensure_openmm_cuda_context","message":"Pushed OpenMM context","data":{"push_result":push_result,"pytorch_already_init":torch.cuda.is_initialized()},"timestamp":int(__import__('time').time()*1000)})+'\n')
            except:
                pass
            
            # Initialize PyTorch CUDA if not already initialized - it should now use OpenMM's context
            if not torch.cuda.is_initialized():
                torch.cuda.init()
            else:
                # PyTorch already initialized - ensure we're using OpenMM's context
                # The context should already be current after cuCtxPushCurrent
                pass
            
            # Set the device to match OpenMM's device
            torch.cuda.set_device(0)
            
            try:
                if DEBUG_FILE_LOGS:
                    import json
                    with open('/media/extradrive/Trajectories/openmm/.cursor/debug.log', 'a') as f:
                        f.write(json.dumps({"sessionId":"debug-session","runId":"run1","hypothesisId":"CONTEXT","location":"_ensure_openmm_cuda_context","message":"Context sharing complete","data":{"is_initialized":torch.cuda.is_initialized(),"current_device":torch.cuda.current_device()},"timestamp":int(__import__('time').time()*1000)})+'\n')
            except:
                pass
        else:
            try:
                if DEBUG_FILE_LOGS:
                    import json
                    with open('/media/extradrive/Trajectories/openmm/.cursor/debug.log', 'a') as f:
                        f.write(json.dumps({"sessionId":"debug-session","runId":"run1","hypothesisId":"CONTEXT","location":"_ensure_openmm_cuda_context","message":"Failed to get OpenMM context","data":{"result":result},"timestamp":int(__import__('time').time()*1000)})+'\n')
            except:
                pass
            torch.cuda.set_device(0)
            
    except Exception as e:
        # If context management fails, fall back to default behavior
        try:
            if DEBUG_FILE_LOGS:
                import json
                with open('/media/extradrive/Trajectories/openmm/.cursor/debug.log', 'a') as f:
                    f.write(json.dumps({"sessionId":"debug-session","runId":"run1","hypothesisId":"CONTEXT","location":"_ensure_openmm_cuda_context","message":"Exception in context management","data":{"error":str(e),"error_type":type(e).__name__},"timestamp":int(__import__('time').time()*1000)})+'\n')
        except:
            pass
        try:
            torch.cuda.set_device(0)
        except:
            pass


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
            
            print(f"✓ PyTorch optimizations enabled:")
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
        valid_dataset_name = task_name  # This is the correct dataset ID the model expects

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
            'full_forces_buffer': full_forces_buffer
        }

        # Single-copy computation function
        def compute_uma_forces_single(state):
            """Compute forces for a single copy - use shared CUDA context with OpenMM."""
            try:
                # CRITICAL DEBUG: Log single-copy calls (rate-limited)
                if DEBUG_LOGS and not cache.get("warned_single", False):
                    print(
                        "⚠️ SINGLE-COPY FUNCTION CALLED (likely from getState/reporters). "
                        "Batch path is still used for RPMD force evaluation.",
                        flush=True,
                    )
                    cache["warned_single"] = True
                
                # #region agent log
                try:
                    import json
                    with open('/media/extradrive/Trajectories/openmm/.cursor/debug.log', 'a') as f:
                        f.write(json.dumps({"sessionId":"debug-session","runId":"run1","hypothesisId":"ALL","location":"umapotential_pythonforce_batch.py:137","message":"Function entry","data":{"cuda_available":torch.cuda.is_available()},"timestamp":int(__import__('time').time()*1000)})+'\n')
                except: pass
                # #endregion
                
                # CRITICAL: With C++ fix, OpenMM's CUDA context is now active before Python is called
                # We need to ensure PyTorch can work with OpenMM's context.
                # Initialize PyTorch CUDA if needed - it should detect and use the current (OpenMM's) context
                if torch.cuda.is_available():
                    if not torch.cuda.is_initialized():
                        torch.cuda.init()
                    torch.cuda.set_device(0)
                    torch.cuda.synchronize()
                
                from fairchem.core.datasets.atomic_data import AtomicData
                from ase import Atoms
                
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

                # #region agent log
                try:
                    import json
                    with open('/media/extradrive/Trajectories/openmm/.cursor/debug.log', 'a') as f:
                        f.write(json.dumps({"sessionId":"debug-session","runId":"run1","hypothesisId":"B","location":"umapotential_pythonforce_batch.py:175","message":"Before setting cell (shared CUDA context)","data":{"cuda_available":torch.cuda.is_available(),"current_device":torch.cuda.current_device() if torch.cuda.is_available() else None},"timestamp":int(__import__('time').time()*1000)})+'\n')
                except: pass
                # #endregion
                
                # CRITICAL: Set cell BEFORE moving to GPU (matches batch path behavior)
                # Now using shared CUDA context, so tensors should work correctly
                with torch.no_grad():
                    if isPeriodic:
                        try:
                            box = state.getPeriodicBoxVectors(asNumpy=True)
                            if box is not None:
                                # Explicit copy to break any OpenMM memory connection
                                cell_np = np.ascontiguousarray(box.value_in_unit(unit.nanometer) * 10.0)
                                
                                # #region agent log
                                try:
                                    with open('/media/extradrive/Trajectories/openmm/.cursor/debug.log', 'a') as f:
                                        f.write(json.dumps({"sessionId":"debug-session","runId":"run1","hypothesisId":"C","location":"umapotential_pythonforce_batch.py:190","message":"After creating cell_np","data":{"cell_shape":list(cell_np.shape),"cell_dtype":str(cell_np.dtype)},"timestamp":int(__import__('time').time()*1000)})+'\n')
                                except: pass
                                # #endregion
                                
                                # CRITICAL: Create cell tensor directly on GPU using PyTorch's allocation
                                # This ensures it's allocated in PyTorch's CUDA context, not OpenMM's
                                # Create empty tensor on GPU first (PyTorch allocates in its own context)
                                cell_gpu = torch.empty((1, 3, 3), dtype=torch.float32, device=device)
                                # Copy data from numpy (CPU) to GPU tensor
                                cell_gpu[0] = torch.from_numpy(cell_np).float()
                                
                                # #region agent log
                                try:
                                    with open('/media/extradrive/Trajectories/openmm/.cursor/debug.log', 'a') as f:
                                        f.write(json.dumps({"sessionId":"debug-session","runId":"run1","hypothesisId":"C","location":"umapotential_pythonforce_batch.py:200","message":"After creating cell_gpu directly on device","data":{"cell_gpu_shape":list(cell_gpu.shape),"cell_gpu_device":str(cell_gpu.device),"cell_gpu_ptr":hex(cell_gpu.data_ptr()) if hasattr(cell_gpu,'data_ptr') else None},"timestamp":int(__import__('time').time()*1000)})+'\n')
                                except: pass
                                # #endregion
                                
                                # Set cell BEFORE .to(device) - matches batch path exactly (line 378)
                                data.cell = cell_gpu
                                
                                # #region agent log
                                try:
                                    with open('/media/extradrive/Trajectories/openmm/.cursor/debug.log', 'a') as f:
                                        f.write(json.dumps({"sessionId":"debug-session","runId":"run1","hypothesisId":"B","location":"umapotential_pythonforce_batch.py:220","message":"After setting data.cell (before .to(device))","data":{"has_cell":hasattr(data,'cell'),"cell_shape":list(data.cell.shape) if hasattr(data,'cell') and hasattr(data.cell,'shape') else None,"cell_device":str(data.cell.device) if hasattr(data,'cell') and hasattr(data.cell,'device') else None,"cell_ptr":hex(data.cell.data_ptr()) if hasattr(data,'cell') and hasattr(data.cell,'data_ptr') else None},"timestamp":int(__import__('time').time()*1000)})+'\n')
                                except: pass
                                # #endregion
                                
                                del cell_gpu
                        except Exception as box_err:
                            # #region agent log
                            try:
                                with open('/media/extradrive/Trajectories/openmm/.cursor/debug.log', 'a') as f:
                                    f.write(json.dumps({"sessionId":"debug-session","runId":"run1","hypothesisId":"ALL","location":"umapotential_pythonforce_batch.py:225","message":"Box vector error","data":{"error":str(box_err),"error_type":type(box_err).__name__},"timestamp":int(__import__('time').time()*1000)})+'\n')
                            except: pass
                            # #endregion
                            print(f"Warning: Box vector error in single-copy: {box_err}")
                    
                    # Now move entire data to GPU (matches batch path line 382)
                    data_device = data.to(device)
                    
                    # #region agent log
                    try:
                        with open('/media/extradrive/Trajectories/openmm/.cursor/debug.log', 'a') as f:
                            f.write(json.dumps({"sessionId":"debug-session","runId":"run1","hypothesisId":"B","location":"umapotential_pythonforce_batch.py:239","message":"After data.to(device)","data":{"has_cell":hasattr(data_device,'cell'),"cell_shape":list(data_device.cell.shape) if hasattr(data_device,'cell') and hasattr(data_device.cell,'shape') else None,"cell_type":str(type(data_device.cell)) if hasattr(data_device,'cell') else None},"timestamp":int(__import__('time').time()*1000)})+'\n')
                    except: pass
                    # #endregion
                    
                    # CRITICAL: dataset must be a LIST (FAIRChem iterates over it)
                    data_device.dataset = [valid_dataset_name]
                    
                    # #region agent log
                    try:
                        with open('/media/extradrive/Trajectories/openmm/.cursor/debug.log', 'a') as f:
                            f.write(json.dumps({"sessionId":"debug-session","runId":"run1","hypothesisId":"ALL","location":"umapotential_pythonforce_batch.py:250","message":"Before predict() call","data":{"has_cell":hasattr(data_device,'cell'),"cell_shape":list(data_device.cell.shape) if hasattr(data_device,'cell') and hasattr(data_device.cell,'shape') else None,"cell_device":str(data_device.cell.device) if hasattr(data_device,'cell') and hasattr(data_device.cell,'device') else None},"timestamp":int(__import__('time').time()*1000)})+'\n')
                    except: pass
                    # #endregion
                    
                    try:
                        pred = predict_unit.predict(data_device)
                    except Exception as pred_err:
                        # #region agent log
                        try:
                            with open('/media/extradrive/Trajectories/openmm/.cursor/debug.log', 'a') as f:
                                f.write(json.dumps({"sessionId":"debug-session","runId":"run1","hypothesisId":"ALL","location":"umapotential_pythonforce_batch.py:258","message":"predict() failed","data":{"error":str(pred_err),"error_type":type(pred_err).__name__,"has_cell":hasattr(data_device,'cell'),"cell_info":str(data_device.cell) if hasattr(data_device,'cell') else None},"timestamp":int(__import__('time').time()*1000)})+'\n')
                        except: pass
                        # #endregion
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
                # CRITICAL DEBUG: Log that batch function is being called
                import time
                batch_start = time.time()
                if DEBUG_LOGS:
                    print(f"🚀 BATCH FUNCTION CALLED! num_beads={len(states)} _from_single={_from_single}", flush=True)
                    print(f"   Using TRUE tensor batching (independent systems)", flush=True)
                
                # Ensure we're using the correct CUDA device
                if torch.cuda.is_available():
                    torch.cuda.set_device(0)  # Use device 0 (matches OpenMM's default)
                
                # Ensure CUDA is synchronized before starting
                torch.cuda.synchronize()
                
                num_copies = len(states)
                
                # =================================================================
                # TRUE TENSOR BATCHING: Build a batched AtomicData object
                # =================================================================
                
                # Collect all positions and boxes
                all_pos_nm = []
                all_boxes = []
                for bead_idx, state in enumerate(states):
                    # Use enforcePeriodicBox to get wrapped positions
                    pos = state.getPositions(asNumpy=True, enforcePeriodicBox=True).value_in_unit(unit.nanometer)
                    if atom_indices is not None:
                        pos = pos[atom_indices].copy()
                    else:
                        pos = pos[:n_atoms].copy()
                    
                    # Validate positions: check for NaN, Inf, or extreme values
                    if np.any(np.isnan(pos)) or np.any(np.isinf(pos)):
                        raise ValueError(f"Invalid positions (NaN/Inf) detected in bead {bead_idx}")
                    if np.any(np.abs(pos) > 1000.0):  # Sanity check: positions > 1000 nm are likely wrong
                        raise ValueError(f"Extreme positions detected in bead {bead_idx}: max={np.max(np.abs(pos)):.2f} nm")
                    
                    all_pos_nm.append(pos)
                    
                    if isPeriodic:
                        try:
                            box = state.getPeriodicBoxVectors(asNumpy=True)
                            if box is not None:
                                box_nm = box.value_in_unit(unit.nanometer).copy()
                                # Validate box vectors
                                if np.any(np.isnan(box_nm)) or np.any(np.isinf(box_nm)):
                                    raise ValueError(f"Invalid box vectors (NaN/Inf) detected in bead {bead_idx}")
                                # Ensure box is reasonable (not zero or negative)
                                box_diag = np.array([box_nm[0][0], box_nm[1][1], box_nm[2][2]])
                                if np.any(box_diag <= 0) or np.any(box_diag > 1000.0):
                                    raise ValueError(f"Invalid box size in bead {bead_idx}: {box_diag}")
                                all_boxes.append(box_nm)
                            else:
                                all_boxes.append(None)
                        except Exception as e:
                            raise ValueError(f"Error getting box vectors for bead {bead_idx}: {e}")
                
                # Convert to Angstroms
                all_pos_angstrom = [pos * 10.0 for pos in all_pos_nm]
                
                # PROFILING
                t_prep = time.time() - batch_start
                
                # =================================================================
                # BATCHED EXECUTION: Single batched model call
                # =================================================================
                t_inference_start = time.time()
                
                from fairchem.core.datasets.atomic_data import AtomicData, atomicdata_list_to_batch
                from ase import Atoms
                
                data_list = []
                for i in range(num_copies):
                    if isPeriodic and i < len(all_boxes) and all_boxes[i] is not None:
                        cell_angstrom = all_boxes[i] * 10.0
                        atoms_ase = Atoms(
                            symbols=symbols,
                            positions=all_pos_angstrom[i],
                            pbc=True,
                            cell=cell_angstrom,
                        )
                    else:
                        atoms_ase = Atoms(
                            symbols=symbols,
                            positions=all_pos_angstrom[i],
                            pbc=isPeriodic,
                        )
                    atoms_ase.info['charge'] = charge
                    atoms_ase.info['spin'] = spin
                    
                    data = AtomicData.from_ase(
                        atoms_ase,
                        task_name=task_name,
                        r_edges=False,
                        r_data_keys=['spin', 'charge'],
                    )
                    data.sid = [f"bead-{i}"]
                    # NOTE: Keep dataset as string - atomicdata_list_to_batch will combine into list
                    data_list.append(data)
                
                batch_data = atomicdata_list_to_batch(data_list)
                batch_indices = batch_data.batch.numpy()
                batch_data = batch_data.to(device)
                
                # Synchronize CUDA before inference to catch async errors
                if torch.cuda.is_available():
                    torch.cuda.synchronize()
                
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
                    print(f"✅ BATCH COMPLETED in {batch_time:.3f}s for {num_copies} beads ({batch_time/num_copies*1000:.1f} ms/bead)", flush=True)
                    print(f"   Breakdown: prep={t_prep*1000:.1f}ms, inference={t_inference*1000:.1f}ms", flush=True)
                    print(f"   Speedup vs sequential: {num_copies * 41.4 / (t_inference*1000):.2f}x (ideal: {num_copies}x)", flush=True)
                
                return (total_energy * unit.kilojoules_per_mole,
                        forces_cpu * unit.kilojoules_per_mole / unit.nanometer)
            except Exception as e:
                # Fallback to sequential if batching fails - but NOT if called from single-copy (avoid recursion!)
                if _from_single:
                    # We're already in fallback mode, just raise the error
                    raise openmm.OpenMMException(f"UMA batch computation failed (from single-copy): {e}")
                
                print(f"Warning: Batch computation failed ({e}), falling back to sequential")
                torch.cuda.empty_cache()  # Clear GPU memory before sequential calls
                
                total_energy = 0.0
                all_forces = []
                
                # Call the ACTUAL single-copy implementation (not the wrapper)
                for state in states:
                    # Direct implementation without delegation to avoid recursion
                    from fairchem.core.datasets.atomic_data import AtomicData
                    from ase import Atoms
                    
                    all_pos_nm = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                    if atom_indices is not None:
                        pos_nm = all_pos_nm[atom_indices]
                    else:
                        pos_nm = all_pos_nm[:n_atoms]
                    
                    pos_angstrom = pos_nm * 10.0
                    atoms_ase = Atoms(symbols=symbols, positions=pos_angstrom, pbc=isPeriodic)
                    atoms_ase.info['charge'] = charge
                    atoms_ase.info['spin'] = spin
                    
                    data = AtomicData.from_ase(atoms_ase, task_name=task_name, r_edges=False, r_data_keys=['spin', 'charge']).to(device)
                    if hasattr(data, 'dataset') and isinstance(data.dataset, str):
                        data.dataset = [data.dataset]
                    
                    if isPeriodic:
                        try:
                            box = state.getPeriodicBoxVectors(asNumpy=True)
                            if box is not None:
                                cell_np = box.value_in_unit(unit.nanometer) * 10.0
                                data.cell = torch.from_numpy(cell_np.copy()).float().to(device).unsqueeze(0)
                        except:
                            pass
                    
                    with torch.no_grad():
                        pred = predict_unit.predict(data)
                    
                    energy_kj = float(pred['energy'].item()) * 96.4853
                    molecular_forces = (pred['forces'] * (96.4853 / 10.0)).cpu().numpy()
                    
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
                    
                    total_energy += energy_kj
                    all_forces.append(forces_kj_nm)
                
                return (total_energy * unit.kilojoules_per_mole,
                        np.array(all_forces) * unit.kilojoules_per_mole / unit.nanometer)

        # Create PythonForce with both single and batch callbacks
        force = openmm.PythonForce(compute_uma_forces_single, {}, compute_uma_forces_batch)
        force.setForceGroup(forceGroup)
        force.setUsesPeriodicBoundaryConditions(isPeriodic)
        system.addForce(force)

        print(f"✓ UMA force added with batched RPMD support (model: {self.model_name}, task: {task_name}, device: {device})")


# No need to register at module level - MLPotential will auto-register from entry points
