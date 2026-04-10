#!/usr/bin/env python3
"""
Diagnostic script to compare UMA vs eSEN model behavior with batched inference.
"""

import torch
import numpy as np
from ase import Atoms
from fairchem.core.calculate import pretrained_mlip
from fairchem.core.datasets.atomic_data import AtomicData, atomicdata_list_to_batch

def test_model_batch(model_name, num_systems=2):
    """Test batched inference for a given model."""
    print(f"\n{'='*80}")
    print(f"Testing: {model_name}")
    print(f"{'='*80}")
    
    try:
        # Load model
        print("Loading model...")
        predictor = pretrained_mlip.get_predict_unit(model_name, device="cuda")
        
        # Check model properties
        model = predictor.model
        if hasattr(model, 'module'):
            inner = model.module
            print(f"Model type: {type(inner).__name__}")
            if hasattr(inner, 'otf_graph'):
                print(f"  otf_graph: {inner.otf_graph}")
            if hasattr(inner, 'use_pbc'):
                print(f"  use_pbc: {inner.use_pbc}")
            if hasattr(inner, 'cutoff'):
                print(f"  cutoff: {inner.cutoff}")
        
        # Create test water molecules
        print(f"\nCreating {num_systems} water molecules...")
        data_list = []
        for i in range(num_systems):
            # Water molecule with slight perturbation
            pos = np.array([
                [0.0 + i*0.01, 0.0, 0.0],      # O
                [0.9572 + i*0.01, 0.0, 0.0],   # H1
                [-0.24, 0.93, 0.0]              # H2
            ])
            
            atoms = Atoms(
                symbols=['O', 'H', 'H'],
                positions=pos,
                pbc=True,
                cell=[10.0, 10.0, 10.0]
            )
            atoms.info['charge'] = 0
            atoms.info['spin'] = 1
            
            # Create AtomicData
            print(f"  System {i}: Creating AtomicData...", end=" ")
            data = AtomicData.from_ase(
                atoms,
                task_name='omol',
                r_edges=False,  # Let model handle edges
                r_data_keys=['spin', 'charge']
            )
            data.sid = [f"mol-{i}"]
            # NOTE: Keep dataset as string for individual data
            # atomicdata_list_to_batch will combine into ['omol', 'omol']
            print(f"(atoms: {data.pos.shape[0]}, dataset: {repr(data.dataset)})")
            data_list.append(data)
        
        # Test 1: Sequential prediction
        print("\n--- Test 1: Sequential Prediction ---")
        sequential_energies = []
        sequential_forces = []
        
        for i, data in enumerate(data_list):
            print(f"  Predicting system {i}...", end=" ")
            # Make a copy to avoid modifying original
            import copy
            data_copy = copy.deepcopy(data)
            data_gpu = data_copy.to("cuda")
            # For single prediction, dataset must be a list
            if isinstance(data_gpu.dataset, str):
                data_gpu.dataset = [data_gpu.dataset]
            
            with torch.no_grad():
                pred = predictor.predict(data_gpu)
            
            energy = pred['energy'].cpu().item()
            forces = pred['forces'].cpu().numpy()
            sequential_energies.append(energy)
            sequential_forces.append(forces)
            print(f"E={energy:.6f} eV")
            
            del data_gpu, pred
            torch.cuda.empty_cache()
        
        # Test 2: Batched prediction
        print("\n--- Test 2: Batched Prediction ---")
        print("Creating batch...", end=" ")
        batch_data = atomicdata_list_to_batch(data_list)
        print(f"(num_graphs: {batch_data.num_graphs}, total_atoms: {batch_data.pos.shape[0]})")
        
        print("Batch properties:")
        print(f"  batch indices: {batch_data.batch}")
        print(f"  natoms: {batch_data.natoms}")
        print(f"  dataset: {batch_data.dataset}")
        print(f"  has edge_index: {hasattr(batch_data, 'edge_index')}")
        if hasattr(batch_data, 'edge_index'):
            print(f"  edge_index shape: {batch_data.edge_index.shape if batch_data.edge_index is not None else None}")
        print(f"  has cell: {hasattr(batch_data, 'cell')}")
        if hasattr(batch_data, 'cell'):
            print(f"  cell shape: {batch_data.cell.shape}")
        
        print("\nMoving batch to GPU...", end=" ")
        batch_data_gpu = batch_data.to("cuda")
        print("")
        
        print("Running batched prediction...", end=" ", flush=True)
        torch.cuda.synchronize()
        
        try:
            with torch.no_grad():
                batch_pred = predictor.predict(batch_data_gpu)
            
            torch.cuda.synchronize()
            print("")
            
            batch_energies = batch_pred['energy'].cpu().numpy()
            batch_forces = batch_pred['forces'].cpu().numpy()
            batch_indices = batch_data.batch.cpu().numpy()
            
            print(f"\nBatch results:")
            print(f"  energies shape: {batch_energies.shape}")
            print(f"  forces shape: {batch_forces.shape}")
            print(f"  batch_indices: {set(batch_indices)}")
            
            # Compare sequential vs batched
            print("\n--- Comparison ---")
            for i in range(num_systems):
                seq_e = sequential_energies[i]
                batch_e = batch_energies[i]
                seq_f = sequential_forces[i]
                batch_f = batch_forces[batch_indices == i]
                
                e_diff = abs(seq_e - batch_e)
                f_diff = np.max(np.abs(seq_f - batch_f))
                
                match = "" if e_diff < 1e-4 and f_diff < 1e-3 else ""
                print(f"  System {i}: {match} ΔE={e_diff:.2e} eV, ΔF_max={f_diff:.2e} eV/Å")
            
            print(f"\nSUCCESS: {model_name} works with batched inference!")
            return True
            
        except Exception as e:
            print(f"FAILED")
            print(f"\nError during batched prediction:")
            print(f"  Type: {type(e).__name__}")
            print(f"  Message: {str(e)}")
            
            # Try to get more info
            if torch.cuda.is_available():
                print(f"  CUDA available: {torch.cuda.is_available()}")
                print(f"  CUDA device: {torch.cuda.current_device()}")
                try:
                    torch.cuda.synchronize()
                except Exception as cuda_e:
                    print(f"  CUDA sync error: {cuda_e}")
            
            return False
            
    except Exception as e:
        print(f"\nFAILED to load or test model")
        print(f"  Error: {type(e).__name__}: {str(e)}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        # Cleanup
        try:
            del predictor, model
            torch.cuda.empty_cache()
        except:
            pass

if __name__ == '__main__':
    models_to_test = [
        'uma-s-1p1',
        'esen-sm-conserving-all-omol',
    ]
    
    results = {}
    for model_name in models_to_test:
        results[model_name] = test_model_batch(model_name, num_systems=2)
    
    print(f"\n{'='*80}")
    print("Summary")
    print(f"{'='*80}")
    for model_name, success in results.items():
        status = "PASS" if success else "FAIL"
        print(f"  {model_name:<40} {status}")
