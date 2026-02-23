#!/usr/bin/env python3
"""Test CACE model outputs for a simple water configuration."""

import sys
import torch
import numpy as np
from pathlib import Path

# Add CACE to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / 'cace'))

from cace.models.atomistic import NeuralNetworkPotential
from cace.data.neighborhood import get_neighborhood

# Load model (openmm root is 2 levels up from tests/cace-lr_water/)
openmm_root = Path(__file__).resolve().parents[2]
model_path = str(openmm_root / "cace" / "water" / "fit" / "fit_version_1" / "best_model.pth")
print(f"Loading CACE model from {model_path}...")
device = torch.device('cpu')
model = torch.load(model_path, map_location=device, weights_only=False)
model.eval()
print(f"Model loaded")
print(f"Model type: {type(model)}")
print(f"Cutoff: {model.representation.cutoff if hasattr(model, 'representation') else 'unknown'}")

# Create a simple water molecule at equilibrium geometry
# O at origin, H at typical positions
oh_bond = 0.9572  # Angstrom
hoh_angle = 104.52 * np.pi / 180.0

positions = np.array([
    [0.0, 0.0, 0.0],  # O
    [oh_bond * np.sin(hoh_angle/2), 0.0, oh_bond * np.cos(hoh_angle/2)],  # H1
    [-oh_bond * np.sin(hoh_angle/2), 0.0, oh_bond * np.cos(hoh_angle/2)],  # H2
])

atomic_numbers = np.array([8, 1, 1])

# Put in a box
box_size = 10.0  # Angstrom
cell = np.eye(3) * box_size

print(f"\nTest configuration:")
print(f"  Single water molecule")
print(f"  O-H bond: {oh_bond:.4f} Å")
print(f"  H-O-H angle: {hoh_angle * 180 / np.pi:.2f}°")
print(f"  Box size: {box_size:.1f} Å")

# Get neighborhood
print(f"\nComputing neighborhood with cutoff {model.representation.cutoff} Å...")
edge_index, shifts, unit_shifts = get_neighborhood(
    positions=positions,
    cutoff=model.representation.cutoff,
    pbc=np.array([True, True, True]),
    cell=cell
)
print(f"  Edges: {len(edge_index[0])}")

# Prepare data
data_dict = {
    'positions': torch.tensor(positions, dtype=torch.float32, device=device, requires_grad=True),
    'atomic_numbers': torch.tensor(atomic_numbers, dtype=torch.long, device=device),
    'edge_index': torch.tensor(edge_index, dtype=torch.long, device=device),
    'shifts': torch.tensor(shifts, dtype=torch.float32, device=device),
    'unit_shifts': torch.tensor(unit_shifts, dtype=torch.float32, device=device),
    'num_nodes': torch.tensor([3], dtype=torch.long, device=device),
    'ptr': torch.tensor([0, 3], dtype=torch.long, device=device),
    'batch': torch.zeros(3, dtype=torch.long, device=device),
    'cell': torch.tensor(cell, dtype=torch.float32, device=device).unsqueeze(0)
}

# Run model
print(f"\nRunning model (inference mode for forces)...")
output = model(data_dict, training=True)

print(f"\nModel outputs:")
for key, value in output.items():
    if isinstance(value, torch.Tensor):
        print(f"  {key}: shape={value.shape}, dtype={value.dtype}")
        if value.numel() <= 10:
            print(f"    value={value.detach().cpu().numpy()}")
        else:
            print(f"    min={value.min():.6f}, max={value.max():.6f}, mean={value.mean():.6f}")

# Check energy
if 'CACE_energy' in output:
    energy_ev = float(output['CACE_energy'].detach().cpu())
    print(f"\nEnergy: {energy_ev:.6f} eV")
    print(f"  = {energy_ev * 96.4853:.2f} kJ/mol")
    print(f"  Expected: ~0 eV for equilibrium geometry (without atomic energies)")
    
# Check forces
if 'CACE_forces' in output:
    forces = output['CACE_forces'].detach().cpu().numpy()
    print(f"\nForces (eV/Å):")
    for i, (Z, pos, force) in enumerate(zip(atomic_numbers, positions, forces)):
        element = 'O' if Z == 8 else 'H'
        force_mag = np.linalg.norm(force)
        print(f"  {element}{i}: {force} (magnitude: {force_mag:.6f} eV/Å)")
    
    max_force = np.max(np.abs(forces))
    print(f"\n  Max force magnitude: {max_force:.6f} eV/Å")
    print(f"  = {max_force * 964.853:.2f} kJ/mol/nm")
    print(f"  Expected: < 0.01 eV/Å for equilibrium geometry")
    
    if max_force > 0.1:
        print(f"\n⚠ WARNING: Forces are too large for equilibrium geometry!")
        print(f"  This suggests the model may have issues.")
