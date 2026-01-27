#!/usr/bin/env python3
"""
Validation test for CACE-LR water model.
Checks charges and dipole of a single molecule.
"""

import numpy as np
import torch
from ase.build import molecule
from cace_dipole_calculator import CACEDipoleCalculator
import os

def validate_single_water():
    print("=" * 60)
    print("Validating CACE-LR Single Water Molecule")
    print("=" * 60)
    
    model_path = "/media/extradrive/Trajectories/openmm/LES-BEC/water/fit/fit_version_1/best_model.pth"
    if not os.path.exists(model_path):
        print(f"Error: Model not found at {model_path}")
        return
        
    calc = CACEDipoleCalculator(model_path, device='cpu')
    
    # Create water molecule
    water = molecule('H2O')
    print(f"\nWater molecule positions (A):")
    for i, (sym, pos) in enumerate(zip(water.get_chemical_symbols(), water.get_positions())):
        print(f"  {i}: {sym:2s} [{pos[0]:7.3f}, {pos[1]:7.3f}, {pos[2]:7.3f}]")
        
    # Predict dipole and charges
    dipole, charges = calc.compute_dipole(water)
    
    print(f"\nPredicted Charges (e):")
    for i, (sym, q) in enumerate(zip(water.get_chemical_symbols(), charges)):
        print(f"  {i}: {sym:2s} = {q:+.4f} e")
    print(f"  Sum = {charges.sum():.6f} e (should ≈ 0)")
    
    dipole_mag_eA = np.linalg.norm(dipole)
    dipole_mag_debye = dipole_mag_eA * 4.803
    
    print(f"\nDipole Moment:")
    print(f"  μ = [{dipole[0]:+.4f}, {dipole[1]:+.4f}, {dipole[2]:+.4f}] e·A")
    print(f"  |μ| = {dipole_mag_eA:.4f} e·A")
    print(f"  |μ| = {dipole_mag_debye:.4f} Debye")
    print(f"  (Expected for TIP3P: ~2.3 Debye, Gas phase: ~1.85 Debye)")
    
    if np.abs(charges.sum()) < 0.05 and dipole_mag_debye > 1.0:
        print("\n✓ SUCCESS: Charges and dipole look reasonable")
    else:
        print("\n✗ FAILED: Charges or dipole outside expected range")

if __name__ == "__main__":
    validate_single_water()
