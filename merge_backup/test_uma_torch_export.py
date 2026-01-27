"""
Test script for torch.export feasibility with FAIRChem UMA models.

This script tests whether torch.export can capture the UMA model backbone,
and if so, whether AOTInductor can compile it for C++ deployment.

Based on the plan: uma_torch.export_feasibility_afcd7e58.plan.md
"""

import os
import sys
import traceback
import copy

import torch
from torch.export import export
from torch.utils._pytree import register_pytree_node


def register_atomic_data_pytree():
    """
    Register AtomicData as a pytree node so torch.export can handle it.
    
    This enables torch.export to properly trace through AtomicData inputs
    by defining how to flatten (extract tensors) and unflatten (reconstruct).
    """
    from fairchem.core.datasets.atomic_data import AtomicData, _REQUIRED_KEYS, _OPTIONAL_KEYS
    
    # Define which keys contain tensors vs non-tensor data
    TENSOR_KEYS = [
        'pos', 'atomic_numbers', 'cell', 'pbc', 'natoms', 'edge_index',
        'cell_offsets', 'nedges', 'charge', 'spin', 'fixed', 'tags',
        'batch', 'energy', 'forces', 'stress'
    ]
    NON_TENSOR_KEYS = ['sid', 'dataset']
    
    def atomic_data_flatten(data: AtomicData):
        """
        Flatten AtomicData into (leaves, context).
        
        leaves: list of tensors (the pytree leaves that torch.export will trace)
        context: tuple of (keys_present, non_tensor_values) for reconstruction
        """
        # Get all keys present in this AtomicData
        keys_present = list(data.keys())
        
        # Separate tensors from non-tensors
        tensor_leaves = []
        tensor_keys = []
        non_tensor_values = {}
        
        for key in keys_present:
            value = getattr(data, key, None)
            if value is None:
                continue
            if torch.is_tensor(value):
                tensor_leaves.append(value)
                tensor_keys.append(key)
            else:
                # Non-tensor data goes in context
                non_tensor_values[key] = copy.deepcopy(value)
        
        # Context contains everything needed to reconstruct
        context = (tuple(tensor_keys), non_tensor_values)
        
        return tensor_leaves, context
    
    def atomic_data_unflatten(leaves, context):
        """
        Unflatten (leaves, context) back into AtomicData.
        
        IMPORTANT: We construct AtomicData without validation to avoid
        data-dependent operations (like batch.max()) that torch.export cannot trace.
        """
        tensor_keys, non_tensor_values = context
        
        # Build full dict of all values
        all_values = {}
        
        # Add tensor values
        for key, tensor in zip(tensor_keys, leaves):
            all_values[key] = tensor
        
        # Add non-tensor values
        for key, value in non_tensor_values.items():
            all_values[key] = value
        
        # Create AtomicData without calling __init__ to skip validation
        # This is needed because validate() uses data-dependent ops
        data = object.__new__(AtomicData)
        
        # Initialize __keys__ set
        data.__keys__ = set()
        
        # Set all attributes directly
        for key, value in all_values.items():
            object.__setattr__(data, key, value)
            data.__keys__.add(key)
        
        # Provide required defaults if missing
        required_defaults = {
            'pos': torch.empty(0, 3),
            'atomic_numbers': torch.empty(0, dtype=torch.long),
            'cell': torch.eye(3).unsqueeze(0),
            'pbc': torch.tensor([[False, False, False]]),
            'natoms': torch.tensor([0], dtype=torch.long),
            'edge_index': torch.empty(2, 0, dtype=torch.long),
            'cell_offsets': torch.empty(0, 3),
            'nedges': torch.tensor([0], dtype=torch.long),
            'charge': torch.tensor([0], dtype=torch.long),
            'spin': torch.tensor([0], dtype=torch.long),
            'fixed': torch.empty(0, dtype=torch.long),
            'tags': torch.empty(0, dtype=torch.long),
        }
        
        for key, default in required_defaults.items():
            if not hasattr(data, key):
                object.__setattr__(data, key, default)
                data.__keys__.add(key)
        
        # Initialize internal state needed by AtomicData
        data.__slices__ = None
        data.__cumsum__ = None
        data.__cat_dims__ = None
        data.__natoms_list__ = None
        
        return data
    
    def atomic_data_flatten_with_keys(data: AtomicData):
        """
        Flatten with keys for better error messages (PyTorch 2.x).
        """
        from torch.utils._pytree import MappingKey
        
        keys_present = list(data.keys())
        tensor_leaves = []
        tensor_keys = []
        non_tensor_values = {}
        
        for key in keys_present:
            value = getattr(data, key, None)
            if value is None:
                continue
            if torch.is_tensor(value):
                tensor_leaves.append((MappingKey(key), value))
                tensor_keys.append(key)
            else:
                non_tensor_values[key] = copy.deepcopy(value)
        
        context = (tuple(tensor_keys), non_tensor_values)
        return tensor_leaves, context
    
    # Register with PyTorch's pytree system
    try:
        # Try the newer API with flatten_with_keys
        register_pytree_node(
            AtomicData,
            flatten_fn=atomic_data_flatten,
            unflatten_fn=atomic_data_unflatten,
            flatten_with_keys_fn=atomic_data_flatten_with_keys,
            serialized_type_name="fairchem.core.datasets.atomic_data.AtomicData"
        )
        print("  Registered AtomicData with flatten_with_keys")
    except TypeError:
        # Fall back to older API without flatten_with_keys
        register_pytree_node(
            AtomicData,
            flatten_fn=atomic_data_flatten,
            unflatten_fn=atomic_data_unflatten,
        )
        print("  Registered AtomicData (basic registration)")
    
    return True


def test_uma_export(model_name: str = "uma-s-1p1", device: str = "cuda"):
    """
    Test whether torch.export can capture the UMA model backbone.
    
    Parameters
    ----------
    model_name : str
        Name of the UMA model to test (default: "uma-s-1p1")
    device : str
        Device to run on (default: "cuda")
        
    Returns
    -------
    tuple
        (success: bool, result: str or path to .pt2 file)
    """
    print(f"=" * 60)
    print(f"Testing torch.export with UMA model: {model_name}")
    print(f"Device: {device}")
    print(f"PyTorch version: {torch.__version__}")
    print(f"=" * 60)
    
    # Step 1: Load the FULL UMA model (including MOLE)
    print("\n[Step 1] Loading FULL UMA model (including MOLE)...")
    try:
        from fairchem.core.calculate import pretrained_mlip
        from fairchem.core.datasets.atomic_data import AtomicData
        from ase import Atoms
        
        # Use default settings - keeps MOLE (Mixture of Linear Experts)
        # We need to export the FULL model, not a merged version
        predict_unit = pretrained_mlip.get_predict_unit(
            model_name,
            inference_settings="default",  # Keep full MOLE model
            device=device
        )
        print(f"  Model loaded successfully")
        print(f"  Model type: {type(predict_unit.model)}")
        print(f"  Inference settings: default (full MOLE model)")
    except Exception as e:
        print(f"  FAILED to load model: {e}")
        traceback.print_exc()
        return False, f"Model loading failed: {e}"
    
    # Step 2: Create sample AtomicData FIRST (needed for lazy init)
    print("\n[Step 2] Creating sample AtomicData...")
    try:
        # Create a simple water molecule
        atoms = Atoms(
            symbols=['O', 'H', 'H'],
            positions=[
                [0.0, 0.0, 0.0],
                [0.957, 0.0, 0.0],
                [-0.24, 0.93, 0.0]
            ]
        )
        atoms.info['charge'] = 0
        atoms.info['spin'] = 1
        
        # Convert to AtomicData
        data = AtomicData.from_ase(
            atoms,
            task_name="omol",
            r_edges=False,
            r_data_keys=['spin', 'charge']
        )
        
        print(f"  AtomicData created successfully")
        print(f"  AtomicData type: {type(data)}")
        print(f"  AtomicData fields: {list(data.keys()) if hasattr(data, 'keys') else 'N/A'}")
        
        # Print tensor shapes
        if hasattr(data, 'pos'):
            print(f"  pos shape: {data.pos.shape}")
        if hasattr(data, 'atomic_numbers'):
            print(f"  atomic_numbers shape: {data.atomic_numbers.shape}")
    except Exception as e:
        print(f"  FAILED to create AtomicData: {e}")
        traceback.print_exc()
        return False, f"AtomicData creation failed: {e}"
    
    # Step 3: Trigger lazy initialization via predict() call
    print("\n[Step 3] Initializing model via predict_unit.predict()...")
    try:
        # This triggers lazy initialization: model is moved to device, 
        # optionally compiled, and MOLE model is merged
        result = predict_unit.predict(data)
        print(f"  Initialization successful!")
        print(f"  Result keys: {list(result.keys())}")
        for key, val in result.items():
            if isinstance(val, torch.Tensor):
                print(f"    {key}: shape={val.shape}, device={val.device}")
    except Exception as e:
        print(f"  FAILED initialization: {e}")
        traceback.print_exc()
        return False, f"Model initialization failed: {e}"
    
    # Step 4: Extract the backbone model (now properly initialized on device)
    print("\n[Step 4] Extracting backbone model...")
    try:
        # The model is typically wrapped in AveragedModel or DataParallel
        if hasattr(predict_unit.model, 'module'):
            backbone = predict_unit.model.module
        else:
            backbone = predict_unit.model
        
        print(f"  Backbone type: {type(backbone)}")
        print(f"  Backbone class: {backbone.__class__.__name__}")
        
        # Check if it's a HydraModel (typical for UMA)
        if hasattr(backbone, 'backbone'):
            print(f"  Inner backbone type: {type(backbone.backbone)}")
            
        # Verify device placement
        sample_param = next(backbone.parameters())
        print(f"  Model device: {sample_param.device}")
    except Exception as e:
        print(f"  FAILED to extract backbone: {e}")
        traceback.print_exc()
        return False, f"Backbone extraction failed: {e}"
    
    # Step 5: Move data to device and test direct forward pass
    print("\n[Step 5] Testing direct forward pass on backbone...")
    try:
        data_device = data.to(device)
        backbone.eval()
        with torch.no_grad():
            output = backbone(data_device)
        print(f"  Forward pass successful!")
        print(f"  Output type: {type(output)}")
        if isinstance(output, dict):
            print(f"  Output keys: {list(output.keys())}")
            for key, val in output.items():
                if isinstance(val, torch.Tensor):
                    print(f"    {key}: shape={val.shape}, dtype={val.dtype}")
                elif isinstance(val, dict):
                    print(f"    {key}: dict with keys {list(val.keys())}")
    except Exception as e:
        print(f"  FAILED forward pass: {e}")
        traceback.print_exc()
        return False, f"Forward pass failed: {e}"
    
    # Step 6: Register AtomicData as pytree node
    print("\n[Step 6] Registering AtomicData as pytree node...")
    try:
        register_atomic_data_pytree()
        print("  Registration successful!")
    except Exception as e:
        print(f"  WARNING: Registration failed: {e}")
        # Continue anyway - might already be registered
    
    # Step 7: Try torch.export with strict=False
    print("\n[Step 7] Attempting torch.export.export(strict=False)...")
    try:
        # Ensure model is in eval mode
        backbone.eval()
        
        # First try with strict=False (more lenient for complex models)
        print("  Trying strict=False mode...")
        with torch.no_grad():
            ep = export(backbone, (data_device,), strict=False)
        
        print(f"  SUCCESS: Model exported!")
        print(f"  ExportedProgram type: {type(ep)}")
        print(f"  Graph signature: {ep.graph_signature}")
        
    except Exception as e1:
        print(f"  strict=False FAILED: {type(e1).__name__}: {e1}")
        
        # Try draft_export for more debugging info
        print("\n  Trying draft_export() for diagnostics...")
        try:
            from torch.export import draft_export
            with torch.no_grad():
                draft_ep = draft_export(backbone, (data_device,))
            print(f"  draft_export succeeded!")
            print(f"  Draft failures: {draft_ep.failures if hasattr(draft_ep, 'failures') else 'N/A'}")
            # Use the draft result
            ep = draft_ep.exported_program if hasattr(draft_ep, 'exported_program') else draft_ep
            print(f"  Proceeding with draft export result...")
        except Exception as e2:
            print(f"  draft_export also FAILED: {type(e2).__name__}: {e2}")
            traceback.print_exc()
            
            # Provide specific error analysis
            error_msg = str(e1)
            if "graph break" in error_msg.lower():
                print("\n  Analysis: Graph breaks detected. May need to use non-strict mode or refactor.")
            elif "dynamic" in error_msg.lower() or "symbol" in error_msg.lower():
                print("\n  Analysis: Dynamic shape/symbol issues. Model has data-dependent sizes.")
            elif "unsupported" in error_msg.lower():
                print("\n  Analysis: Unsupported operations. May need custom decompositions.")
            elif "pytree" in error_msg.lower() or "flatten" in error_msg.lower():
                print("\n  Analysis: Input data structure not supported. Check AtomicData registration.")
            elif "control flow" in error_msg.lower() or "cond" in error_msg.lower():
                print("\n  Analysis: Data-dependent control flow. May need torch.cond().")
            elif "could not be traced" in error_msg.lower() or "guard" in error_msg.lower():
                print("\n  Analysis: Tracing issues. Model may have data-dependent behavior.")
            elif "keyerror" in error_msg.lower():
                print("\n  Analysis: KeyError - likely symbolic shape resolution failed.")
            
            return False, f"torch.export failed: {e1}"
    
    # Step 8: Try AOTInductor compilation
    print("\n[Step 8] Attempting AOTInductor compilation...")
    try:
        output_path = os.path.join(os.getcwd(), f"{model_name.replace('-', '_')}_aot.pt2")
        
        with torch.no_grad():
            pt2_path = torch._inductor.aoti_compile_and_package(
                ep,
                package_path=output_path
            )
        
        print(f"  SUCCESS: Compiled to {pt2_path}")
        
        # Verify the file exists
        if os.path.exists(pt2_path):
            file_size = os.path.getsize(pt2_path) / (1024 * 1024)  # MB
            print(f"  File size: {file_size:.2f} MB")
        
        return True, pt2_path
        
    except Exception as e:
        print(f"  FAILED AOTInductor: {type(e).__name__}: {e}")
        traceback.print_exc()
        return False, f"AOTInductor failed (but torch.export worked): {e}"


def test_torch_compile_fallback(model_name: str = "uma-s-1p1", device: str = "cuda"):
    """
    Test torch.compile as a fallback - this is what FAIRChem already supports.
    
    This confirms the model CAN be compiled, just not exported for C++ deployment.
    """
    print(f"\n{'=' * 60}")
    print(f"Testing torch.compile (fallback) with: {model_name}")
    print(f"{'=' * 60}")
    
    try:
        from fairchem.core.calculate import pretrained_mlip
        from fairchem.core.datasets.atomic_data import AtomicData
        from ase import Atoms
        
        # Load FULL model with default settings (keeps MOLE)
        predict_unit = pretrained_mlip.get_predict_unit(
            model_name,
            inference_settings="default",  # Keep full MOLE model
            device=device
        )
        
        # Create sample data
        atoms = Atoms(
            symbols=['O', 'H', 'H'],
            positions=[
                [0.0, 0.0, 0.0],
                [0.957, 0.0, 0.0],
                [-0.24, 0.93, 0.0]
            ]
        )
        atoms.info['charge'] = 0
        atoms.info['spin'] = 1
        
        data = AtomicData.from_ase(
            atoms,
            task_name="omol",
            r_edges=False,
            r_data_keys=['spin', 'charge']
        )
        
        # IMPORTANT: First call predict() to trigger lazy initialization
        print("  Initializing model via predict()...")
        predict_unit.predict(data)
        
        # Now get backbone (properly initialized on device)
        if hasattr(predict_unit.model, 'module'):
            backbone = predict_unit.model.module
        else:
            backbone = predict_unit.model
        
        data_device = data.to(device)
        
        # Compile with dynamic=True (same as FAIRChem does)
        print("  Compiling model with torch.compile(dynamic=True)...")
        torch.compiler.reset()
        compiled_backbone = torch.compile(backbone, dynamic=True)
        
        # Run compiled model
        print("  Running compiled model...")
        compiled_backbone.eval()
        with torch.no_grad():
            output = compiled_backbone(data_device)
        
        print(f"  SUCCESS: torch.compile works!")
        print(f"  Output keys: {list(output.keys()) if isinstance(output, dict) else type(output)}")
        
        return True, "torch.compile works as fallback"
        
    except Exception as e:
        print(f"  FAILED: {e}")
        traceback.print_exc()
        return False, str(e)


def main():
    """Main entry point for the feasibility test."""
    print("\n" + "=" * 70)
    print("UMA torch.export Feasibility Test")
    print("=" * 70)
    
    # Check CUDA availability
    if torch.cuda.is_available():
        device = "cuda"
        print(f"\nCUDA available: {torch.cuda.get_device_name(0)}")
    else:
        device = "cpu"
        print("\nCUDA not available, using CPU")
    
    # Test torch.export
    success, result = test_uma_export(model_name="uma-s-1p1", device=device)
    
    if success:
        print("\n" + "=" * 70)
        print("RESULT: torch.export + AOTInductor WORKS!")
        print("=" * 70)
        print(f"\nProceeding with C++ plugin is feasible.")
        print(f"Compiled model: {result}")
        print("\nNext steps:")
        print("  1. Create export script for all UMA models")
        print("  2. Implement C++ plugin using AOTIModelPackageLoader")
        print("  3. Build AtomicData construction in C++")
    else:
        print("\n" + "=" * 70)
        print("RESULT: torch.export FAILED")
        print("=" * 70)
        print(f"\nError: {result}")
        
        # Test fallback
        print("\nTesting torch.compile as fallback...")
        fb_success, fb_result = test_torch_compile_fallback(device=device)
        
        if fb_success:
            print("\n" + "=" * 70)
            print("FALLBACK: torch.compile works!")
            print("=" * 70)
            print("\nRecommendation: Use Optimized PythonForce approach")
            print("  - torch.compile the model for faster inference")
            print("  - Keep using PythonForce for OpenMM integration")
            print("  - Optimize data conversion paths")
        else:
            print("\n" + "=" * 70)
            print("CRITICAL: Both torch.export and torch.compile failed!")
            print("=" * 70)
            print("\nThis indicates fundamental issues with the model.")
    
    return success


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
