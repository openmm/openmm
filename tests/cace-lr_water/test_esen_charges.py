#!/usr/bin/env python3
"""Test eSEN charge prediction on water molecule."""

import numpy as np
import torch
from ase import Atoms
from ase.build import molecule
from fairchem.core import FAIRChemCalculator
from fairchem.core.calculate import pretrained_mlip

def test_esen_charges():
    """Test charge extraction from eSEN model."""
    
    # Load model (try different variants)
    model_names = [
        "uma-s-1p1",
        "uma-s-1",
        "uma-m-1p1",
        "esen-md-direct-all-omol",
        "esen-sm-conserving-all-omol",
        "esen-sm-direct-all-omol"
    ]
    
    # Create single water molecule
    water = molecule('H2O')
    water.info['charge'] = 0
    water.info['spin'] = 1  # Singlet
    
    print("Testing charge prediction on water molecule")
    print("=" * 60)
    print(f"Water molecule:")
    print(f"  Atoms: {water.get_chemical_symbols()}")
    print(f"  Positions (Å):")
    for i, (sym, pos) in enumerate(zip(water.get_chemical_symbols(), water.get_positions())):
        print(f"    {i}: {sym:2s} [{pos[0]:7.3f}, {pos[1]:7.3f}, {pos[2]:7.3f}]")
    
    for model_name in model_names:
        print(f"\n{'='*60}")
        print(f"Trying model: {model_name}")
        print('='*60)
        try:
            # Load predictor
            predictor = pretrained_mlip.get_predict_unit(model_name, device='cuda')
            print(f"  Model loaded")
            
            # Create calculator for omol task
            if "omol" in predictor.dataset_to_tasks:
                calc = FAIRChemCalculator(predictor, task_name="omol")
            else:
                print(f"  Model does not support 'omol' task")
                continue
            
            # Attach calculator
            water.calc = calc
            
            # Get energy and forces
            energy = water.get_potential_energy()
            forces = water.get_forces()
            
            print(f"  Energy: {energy:.4f} eV")
            print(f"  Forces shape: {forces.shape}")
            
            # Check if model has charge prediction
            head = predictor.model.heads[0] if hasattr(predictor.model, 'heads') else predictor.model
            
            if hasattr(head, 'get_charges'):
                print(f"  Model has get_charges() method")
                
                # We need to call the model directly to get charges
                # Convert atoms to AtomicData
                from fairchem.core.datasets.atomic_data import AtomicData
                
                # Create minimal data dict
                data = AtomicData(
                    pos=torch.tensor(water.get_positions(), dtype=torch.float32),
                    atomic_numbers=torch.tensor(water.get_atomic_numbers(), dtype=torch.long),
                    natoms=torch.tensor([len(water)], dtype=torch.long),
                    charge=torch.tensor([water.info['charge']], dtype=torch.long),
                    spin=torch.tensor([water.info['spin']], dtype=torch.long),
                    cell=torch.eye(3).unsqueeze(0) * 100.0,  # Large box for non-periodic
                    pbc=torch.tensor([False, False, False]),
                ).to('cuda')
                
                # Run forward pass
                with torch.no_grad():
                    emb = predictor.model.backbone(data)
                    
                    # Try to get charges
                    try:
                        charge_dict = head.get_charges(
                            emb["node_embedding"].narrow(1, 0, 1).squeeze(),
                            data
                        )
                        
                        charges = charge_dict["charges"].cpu().numpy()
                        print(f"  Extracted charges (shape: {charges.shape}):")
                        
                        # Flatten charges
                        if len(charges.shape) == 3:
                            charges_flat = charges[:, 0, 0]
                        elif len(charges.shape) == 2:
                            charges_flat = charges[:, 0]
                        else:
                            charges_flat = charges
                        
                        for i, (sym, q) in enumerate(zip(water.get_chemical_symbols(), charges_flat)):
                            print(f"    {i}: {sym:2s} = {q:+.4f} e")
                        print(f"    Sum: {charges_flat.sum():.6f} e (should ≈ 0)")
                        
                        # Compute dipole
                        positions_nm = water.get_positions() / 10.0  # Å → nm
                        dipole = np.sum(charges_flat[:, None] * positions_nm, axis=0)
                        dipole_magnitude = np.linalg.norm(dipole)
                        
                        print(f"  Dipole moment:")
                        print(f"    μ = [{dipole[0]:+.4f}, {dipole[1]:+.4f}, {dipole[2]:+.4f}] e·nm")
                        print(f"    |μ| = {dipole_magnitude:.4f} e·nm")
                        print(f"    |μ| = {dipole_magnitude * 4.803:.4f} Debye")
                        print(f"    (Expected for water: ~1.85 Debye)")
                        
                        return True, model_name, charges_flat
                        
                    except Exception as e:
                        print(f"  Error calling get_charges(): {e}")
                        import traceback
                        traceback.print_exc()
                        
            else:
                print(f"  Model does not have get_charges() method")
                print(f"    Head type: {type(head)}")
                print(f"    Head methods: {[m for m in dir(head) if not m.startswith('_')][:10]}")
                
        except Exception as e:
            print(f"  Error: {e}")
            continue
    
    print("\n" + "=" * 60)
    print("No suitable model found with charge prediction!")
    return False, None, None

if __name__ == "__main__":
    success, model_name, charges = test_esen_charges()
    if success:
        print(f"\n{'='*60}")
        print(f"SUCCESS! Use model: {model_name}")
        print('='*60)
    else:
        print(f"\n{'='*60}")
        print("FAILED - No models with charge prediction found")
        print('='*60)
