#!/usr/bin/env python3
"""
Test eSEN with larger system (similar to actual simulation).
"""

import torch
import numpy as np
from ase import Atoms
from fairchem.core.calculate import pretrained_mlip
from fairchem.core.datasets.atomic_data import AtomicData, atomicdata_list_to_batch

def test_larger_system(model_name, num_molecules=50, num_beads=16):
    """Test eSEN with a larger water system similar to actual simulation."""
    print(f"\n{'='*80}")
    print(f"Testing {model_name} with {num_molecules} molecules, {num_beads} beads")
    print(f"{'='*80}")
    
    # Load model
    print("Loading model...")
    predictor = pretrained_mlip.get_predict_unit(model_name, device="cuda")
    
    # Create a larger water system
    print(f"Creating {num_molecules} water molecules...")
    num_atoms = num_molecules * 3
    
    # Simple grid arrangement
    molecules_per_side = int(np.ceil(num_molecules ** (1/3)))
    spacing = 0.3  # nm
    
    positions = []
    symbols = []
    mol_count = 0
    
    for i in range(molecules_per_side):
        for j in range(molecules_per_side):
            for k in range(molecules_per_side):
                if mol_count >= num_molecules:
                    break
                
                # Water molecule at this grid point
                base = np.array([i * spacing, j * spacing, k * spacing])
                positions.extend([
                    base,
                    base + np.array([0.09572, 0, 0]),
                    base + np.array([-0.024, 0.093, 0])
                ])
                symbols.extend(['O', 'H', 'H'])
                mol_count += 1
            if mol_count >= num_molecules:
                break
        if mol_count >= num_molecules:
            break
    
    positions = np.array(positions[:num_atoms])
    symbols = symbols[:num_atoms]
    
    box_size = molecules_per_side * spacing + 0.5
    print(f"  Total atoms: {len(positions)}")
    print(f"  Box size: {box_size:.3f} nm")
    
    # Test batched prediction with multiple beads
    print(f"\nCreating {num_beads} RPMD beads...")
    data_list = []
    
    for bead_idx in range(num_beads):
        # Add small perturbation to each bead
        perturbed_pos = positions + np.random.randn(*positions.shape) * 0.001
        
        atoms = Atoms(
            symbols=symbols,
            positions=perturbed_pos * 10,  # Convert to Angstroms
            pbc=True,
            cell=[box_size*10, box_size*10, box_size*10]
        )
        atoms.info['charge'] = 0
        atoms.info['spin'] = 1
        
        data = AtomicData.from_ase(
            atoms,
            task_name='omol',
            r_edges=False,
            r_data_keys=['spin', 'charge']
        )
        data.sid = [f"bead-{bead_idx}"]
        data_list.append(data)
    
    print(f"  Created {len(data_list)} beads")
    
    # Create batch
    print("Creating batched data...")
    batch_data = atomicdata_list_to_batch(data_list)
    print(f"  Batch: {batch_data.num_graphs} graphs, {batch_data.pos.shape[0]} total atoms")
    print(f"  dataset type: {type(batch_data.dataset)}")
    print(f"  dataset value: {batch_data.dataset[:3]}...")
    
    # Move to GPU and predict
    print("Moving to GPU and running prediction...")
    batch_data_gpu = batch_data.to("cuda")
    torch.cuda.synchronize()
    
    try:
        with torch.no_grad():
            pred = predictor.predict(batch_data_gpu)
        
        torch.cuda.synchronize()
        
        energies = pred['energy'].cpu().numpy()
        forces = pred['forces'].cpu().numpy()
        
        print(f"\n✓ SUCCESS!")
        print(f"  Energy shape: {energies.shape}")
        print(f"  Forces shape: {forces.shape}")
        print(f"  Mean energy: {np.mean(energies):.2f} eV")
        
        return True
        
    except Exception as e:
        print(f"\n✗ FAILED: {type(e).__name__}: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == '__main__':
    # Test with uma first (should work)
    test_larger_system('uma-s-1p1', num_molecules=50, num_beads=16)
    
    # Then test esen
    test_larger_system('esen-sm-conserving-all-omol', num_molecules=50, num_beads=16)
