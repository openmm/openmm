#!/usr/bin/env python3
"""
Validate UMA+LES model on OMol25 test set.

This script loads a trained UMA+LES checkpoint and evaluates it on the OMol25
test set, computing energy and force MAE metrics.

Usage:
    python validate_uma_les.py --checkpoint /path/to/checkpoint.pt --test-data /path/to/test.aselmdb
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import torch
from ase.db import connect
from tqdm import tqdm

try:
    from fairchem.core.units.mlip_unit import load_predict_unit
    from fairchem.core import FAIRChemCalculator
except ImportError:
    print("Error: fairchem-core not installed. Install with: pip install fairchem-core")
    sys.exit(1)


def validate_model(checkpoint_path: str, test_data_path: str, limit: int = 1000, task_name: str = "omol"):
    """
    Validate trained model on test set.
    
    Args:
        checkpoint_path: Path to trained model checkpoint
        test_data_path: Path to test data (ASE DB)
        limit: Number of samples to evaluate
        task_name: Task name for FAIRChemCalculator
    """
    print("=" * 60)
    print("UMA+LES Model Validation")
    print("=" * 60)
    print(f"Checkpoint: {checkpoint_path}")
    print(f"Test data: {test_data_path}")
    print(f"Task: {task_name}")
    print(f"Limit: {limit}")
    print("=" * 60)
    
    # Load model
    print("\nLoading model...")
    try:
        predictor = load_predict_unit(checkpoint_path)
        calc = FAIRChemCalculator(predictor, task_name=task_name)
        print("✓ Model loaded successfully")
    except Exception as e:
        print(f"✗ Failed to load model: {e}")
        return None
    
    # Load test data
    print(f"\nLoading test data from {test_data_path}...")
    try:
        db = connect(test_data_path)
        print(f"✓ Test database loaded")
    except Exception as e:
        print(f"✗ Failed to load test data: {e}")
        return None
    
    # Validate
    print(f"\nValidating on {limit} samples...")
    errors = {
        "energy": [],
        "forces": [],
        "energy_per_atom": []
    }
    
    failed = 0
    
    for i, row in enumerate(tqdm(db.select(limit=limit), total=limit)):
        try:
            atoms = row.toatoms()
            atoms.calc = calc
            
            # Predict
            pred_energy = atoms.get_potential_energy()
            pred_forces = atoms.get_forces()
            
            # Get ground truth
            true_energy = row.energy
            true_forces = np.array(row.data.forces) if hasattr(row.data, 'forces') else atoms.get_forces()
            
            # Compute errors
            n_atoms = len(atoms)
            energy_error = abs(pred_energy - true_energy)
            force_error = np.mean(np.abs(pred_forces - true_forces))
            
            errors["energy"].append(energy_error)
            errors["forces"].append(force_error)
            errors["energy_per_atom"].append(energy_error / n_atoms)
            
        except Exception as e:
            failed += 1
            if failed < 10:  # Only print first 10 failures
                print(f"\n✗ Failed on sample {i}: {e}")
            continue
    
    # Compute statistics
    print("\n" + "=" * 60)
    print("Validation Results")
    print("=" * 60)
    print(f"Samples evaluated: {len(errors['energy'])}")
    print(f"Failures: {failed}")
    print()
    print(f"Energy MAE: {np.mean(errors['energy']):.4f} eV")
    print(f"Energy per atom MAE: {np.mean(errors['energy_per_atom']):.4f} eV/atom")
    print(f"Forces MAE: {np.mean(errors['forces']):.4f} eV/Å")
    print()
    print(f"Energy Median AE: {np.median(errors['energy']):.4f} eV")
    print(f"Forces Median AE: {np.median(errors['forces']):.4f} eV/Å")
    print("=" * 60)
    
    return errors


def main():
    parser = argparse.ArgumentParser(
        description="Validate UMA+LES model on OMol25 test set",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--checkpoint", "-c",
        type=str,
        required=True,
        help="Path to trained model checkpoint (.pt file)"
    )
    
    parser.add_argument(
        "--test-data", "-t",
        type=str,
        required=True,
        help="Path to test data (ASE DB)"
    )
    
    parser.add_argument(
        "--limit", "-n",
        type=int,
        default=1000,
        help="Number of samples to evaluate"
    )
    
    parser.add_argument(
        "--task-name",
        type=str,
        default="omol",
        choices=["omol", "omat", "oc20", "odac", "omc"],
        help="Task name for FAIRChemCalculator"
    )
    
    parser.add_argument(
        "--output", "-o",
        type=str,
        default=None,
        help="Output file for detailed results (optional)"
    )
    
    args = parser.parse_args()
    
    # Validate inputs
    if not Path(args.checkpoint).exists():
        print(f"Error: Checkpoint file not found: {args.checkpoint}")
        sys.exit(1)
    
    if not Path(args.test_data).exists():
        print(f"Error: Test data not found: {args.test_data}")
        sys.exit(1)
    
    # Run validation
    errors = validate_model(
        checkpoint_path=args.checkpoint,
        test_data_path=args.test_data,
        limit=args.limit,
        task_name=args.task_name
    )
    
    # Save results if requested
    if args.output and errors:
        print(f"\nSaving results to {args.output}...")
        np.savez(
            args.output,
            energy_errors=errors["energy"],
            force_errors=errors["forces"],
            energy_per_atom_errors=errors["energy_per_atom"]
        )
        print("✓ Results saved")


if __name__ == "__main__":
    main()
