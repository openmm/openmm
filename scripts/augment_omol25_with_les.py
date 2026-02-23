#!/usr/bin/env python3
"""
Augment OMol25 dataset with LES (Latent Ewald Summation) features.

This script adds long-range electrostatic corrections to the OMol25 dataset
using the LES method implemented in CACE. The augmented dataset can then be
used to fine-tune UMA models for improved long-range electrostatics.

Usage:
    python augment_omol25_with_les.py \
        --input ./omol25_data/train_4M \
        --output ./omol25_les_augmented \
        --num-samples 10000 \
        --sigma 1.0 \
        --decomposition-level 2

Author: Scientific ML Pipeline
Date: 2026-01-23
"""

import argparse
import os
import sys
from pathlib import Path
from typing import Optional, Dict, Any
import logging

import numpy as np
import torch
from tqdm import tqdm
from ase import Atoms
from ase.db import connect

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def setup_paths():
    """Ensure CACE and FAIRChem are available."""
    try:
        from cace.modules import EwaldPotential
        from fairchem.core.datasets import AseDBDataset
        logger.info("Successfully imported CACE and FAIRChem modules")
        return True
    except ImportError as e:
        logger.error(f"Import error: {e}")
        logger.error("Please install: pip install fairchem-core && cd cace && pip install -e .")
        return False


def estimate_charges_from_composition(atoms: Atoms) -> np.ndarray:
    """
    Estimate atomic charges based on elemental composition.
    
    This is a simple heuristic. For production, you would want to:
    1. Use UMA's charge prediction if available
    2. Use a pre-trained charge model (e.g., from CACE)
    3. Use electronegativity-based methods
    
    Args:
        atoms: ASE Atoms object
        
    Returns:
        Estimated charges as numpy array
    """
    symbols = atoms.get_chemical_symbols()
    charges = np.zeros(len(atoms))
    
    # Simple electronegativity-based estimate
    # Pauling electronegativities (simplified)
    electronegativity = {
        'H': 2.20, 'C': 2.55, 'N': 3.04, 'O': 3.44, 'F': 3.98,
        'S': 2.58, 'Cl': 3.16, 'P': 2.19, 'Br': 2.96, 'I': 2.66,
        'Na': 0.93, 'K': 0.82, 'Ca': 1.00, 'Mg': 1.31,
        'Fe': 1.83, 'Cu': 1.90, 'Zn': 1.65, 'Ni': 1.91,
        'Si': 1.90, 'Al': 1.61,
    }
    
    # Calculate mean electronegativity
    electroneg_values = [electronegativity.get(s, 2.0) for s in symbols]
    mean_electroneg = np.mean(electroneg_values)
    
    # Assign partial charges based on deviation from mean
    for i, symbol in enumerate(symbols):
        en = electronegativity.get(symbol, 2.0)
        charges[i] = (mean_electroneg - en) * 0.1  # Scale factor
    
    # Ensure total charge matches system charge
    total_charge = atoms.info.get('charge', 0)
    current_total = charges.sum()
    if len(charges) > 0:
        charges += (total_charge - current_total) / len(charges)
    
    return charges


def compute_les_features(
    atoms: Atoms,
    ewald_module: 'EwaldPotential',
    device: torch.device
) -> Dict[str, Any]:
    """
    Compute LES features for a molecular structure.
    
    Args:
        atoms: ASE Atoms object
        ewald_module: Initialized EwaldPotential module
        device: PyTorch device
        
    Returns:
        Dictionary with LES features (ewald_energy, ewald_field, etc.)
    """
    # Get atomic properties
    positions = torch.tensor(atoms.get_positions(), dtype=torch.float32).to(device)
    charges = estimate_charges_from_composition(atoms)
    charges = torch.tensor(charges, dtype=torch.float32).unsqueeze(1).to(device)
    
    # Get cell (use large box if non-periodic)
    if atoms.pbc.any():
        cell = torch.tensor(atoms.get_cell(), dtype=torch.float32).to(device)
    else:
        # For non-periodic systems, use large box
        max_extent = positions.abs().max().item()
        box_size = max(max_extent * 3, 50.0)  # At least 50 Å
        cell = torch.eye(3, dtype=torch.float32).to(device) * box_size
    
    # Compute Ewald potential and field
    try:
        ewald_energy, ewald_field = ewald_module.compute_potential_optimized(
            positions, charges, cell, compute_field=True
        )
        
        # Also compute real-space component for comparison
        ewald_energy_real, ewald_field_real = ewald_module.compute_potential_realspace(
            positions, charges, compute_field=True
        )
        
        return {
            'ewald_energy': ewald_energy.cpu().item(),
            'ewald_field': ewald_field.cpu().numpy(),
            'ewald_energy_real': ewald_energy_real.cpu().item(),
            'ewald_field_real': ewald_field_real.cpu().numpy(),
            'charges': charges.squeeze().cpu().numpy(),
        }
    except Exception as e:
        logger.warning(f"Failed to compute LES features: {e}")
        return None


def augment_dataset(
    input_path: str,
    output_path: str,
    num_samples: Optional[int] = None,
    sigma: float = 1.0,
    decomposition_level: int = 2,
    device: str = 'cuda'
) -> None:
    """
    Augment OMol25 dataset with LES features.
    
    Args:
        input_path: Path to input LMDB or ASE DB
        output_path: Path for output augmented dataset
        num_samples: Number of samples to process (None = all)
        sigma: Ewald parameter (width of Gaussian)
        decomposition_level: Decomposition level for LES
        device: Device to use ('cuda' or 'cpu')
    """
    from cace.modules import EwaldPotential
    from fairchem.core.datasets import AseDBDataset
    
    logger.info(f"Loading dataset from {input_path}")
    
    # Load dataset
    try:
        dataset = AseDBDataset({"src": input_path})
        logger.info(f"Loaded dataset with {len(dataset)} structures")
    except Exception as e:
        logger.error(f"Failed to load dataset: {e}")
        return
    
    # Initialize Ewald module
    device_torch = torch.device(device if torch.cuda.is_available() else 'cpu')
    logger.info(f"Using device: {device_torch}")
    
    ewald = EwaldPotential(
        dl=decomposition_level,
        sigma=sigma,
        exponent=1,  # 1/r potential
        feature_key='q',
        aggregation_mode='sum',
        compute_field=True
    ).to(device_torch)
    
    # Create output directory
    os.makedirs(output_path, exist_ok=True)
    output_db = Path(output_path) / 'augmented.db'
    
    # Connect to output database
    db = connect(str(output_db))
    
    # Process structures
    num_to_process = min(num_samples, len(dataset)) if num_samples else len(dataset)
    logger.info(f"Processing {num_to_process} structures...")
    
    success_count = 0
    failure_count = 0
    
    for idx in tqdm(range(num_to_process)):
        try:
            # Get atoms
            atoms = dataset.get_atoms(idx)
            
            # Compute LES features
            les_features = compute_les_features(atoms, ewald, device_torch)
            
            if les_features is None:
                failure_count += 1
                continue
            
            # Add LES features to atoms.info
            atoms.info['ewald_energy'] = les_features['ewald_energy']
            atoms.info['ewald_energy_real'] = les_features['ewald_energy_real']
            atoms.info['ewald_contribution'] = (
                les_features['ewald_energy'] - les_features['ewald_energy_real']
            )
            
            # Store as arrays (ASE DB can handle this)
            atoms.arrays['ewald_field'] = les_features['ewald_field']
            atoms.arrays['estimated_charges'] = les_features['charges']
            
            # Write to database
            db.write(
                atoms,
                data={
                    'original_idx': idx,
                    'source': atoms.info.get('source', 'omol25'),
                    'charge': atoms.info.get('charge', 0),
                    'spin': atoms.info.get('spin', 1),
                }
            )
            
            success_count += 1
            
        except Exception as e:
            logger.debug(f"Error processing structure {idx}: {e}")
            failure_count += 1
            continue
    
    # Summary
    logger.info(f"\n{'='*60}")
    logger.info(f"Augmentation complete!")
    logger.info(f"Successfully processed: {success_count}/{num_to_process}")
    logger.info(f"Failures: {failure_count}/{num_to_process}")
    logger.info(f"Output saved to: {output_db}")
    logger.info(f"{'='*60}\n")
    
    # Write metadata
    metadata_file = Path(output_path) / 'metadata.txt'
    with open(metadata_file, 'w') as f:
        f.write(f"Input dataset: {input_path}\n")
        f.write(f"Number of structures: {success_count}\n")
        f.write(f"Ewald sigma: {sigma}\n")
        f.write(f"Decomposition level: {decomposition_level}\n")
        f.write(f"Device: {device_torch}\n")
        f.write(f"Failures: {failure_count}\n")
    
    logger.info(f"Metadata written to: {metadata_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Augment OMol25 dataset with LES features",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '--input', '-i',
        type=str,
        required=True,
        help='Path to input OMol25 dataset (LMDB or ASE DB)'
    )
    
    parser.add_argument(
        '--output', '-o',
        type=str,
        required=True,
        help='Path for output augmented dataset'
    )
    
    parser.add_argument(
        '--num-samples', '-n',
        type=int,
        default=None,
        help='Number of samples to process (default: all)'
    )
    
    parser.add_argument(
        '--sigma',
        type=float,
        default=1.0,
        help='Ewald parameter (width of Gaussian)'
    )
    
    parser.add_argument(
        '--decomposition-level', '-dl',
        type=int,
        default=2,
        help='Decomposition level for LES'
    )
    
    parser.add_argument(
        '--device',
        type=str,
        default='cuda',
        choices=['cuda', 'cpu'],
        help='Device to use for computation'
    )
    
    parser.add_argument(
        '--check-imports',
        action='store_true',
        help='Only check if required packages are installed'
    )
    
    args = parser.parse_args()
    
    # Check imports
    if not setup_paths():
        sys.exit(1)
    
    if args.check_imports:
        logger.info("All required packages are available")
        sys.exit(0)
    
    # Run augmentation
    augment_dataset(
        input_path=args.input,
        output_path=args.output,
        num_samples=args.num_samples,
        sigma=args.sigma,
        decomposition_level=args.decomposition_level,
        device=args.device
    )


if __name__ == '__main__':
    main()
