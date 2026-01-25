import numpy as np
import torch
from cace.calculators import CACECalculator
from cace.data import AtomicData
from cace.data.neighborhood import get_neighborhood
from ase import Atoms

class CACEDipoleCalculator:
    """
    Calculator for CACE charges and dipole moments.
    Used during MD simulation to record dipole trajectory.
    """
    def __init__(self, model_path, device='cpu', charge_unit=None):
        """
        Parameters
        ----------
        model_path : str
            Path to pretrained CACE-LR model (.pth)
        device : str
            Device to run on ('cpu' or 'cuda')
        charge_unit : float, optional
            Scaling factor for charges. If None, uses CACE default 1/sqrt(90.0474).
        """
        if charge_unit is None:
            # Standard normalization factor in CACE convention
            self.charge_unit = 1.0 / np.sqrt(90.0474)
        else:
            self.charge_unit = charge_unit
            
        print(f"Loading CACE dipole model from '{model_path}' on {device}...")
        self.device = torch.device(device)
        self.model = torch.load(model_path, map_location=self.device, weights_only=False)
        self.model.eval()
        
        # Get cutoff from model
        try:
            self.cutoff = self.model.representation.cutoff
        except AttributeError:
            self.cutoff = self.model.models[0].representation.cutoff
            
        print(f"  Model cutoff: {self.cutoff} A")

    def compute_dipole(self, atoms):
        """
        Compute dipole moment from CACE predicted charges.
        
        Parameters
        ----------
        atoms : ase.Atoms
            Current atomic configuration
            
        Returns
        -------
        dipole : np.ndarray
            Dipole moment vector [dx, dy, dz] in e·A
        charges : np.ndarray
            Predicted atomic charges in e
        """
        # Get neighborhood
        positions = atoms.get_positions()
        pbc = tuple(atoms.get_pbc())
        cell = np.array(atoms.get_cell())
        
        edge_index, shifts, unit_shifts = get_neighborhood(
            positions=positions,
            cutoff=self.cutoff,
            pbc=pbc,
            cell=cell
        )
        
        # Construct data dict
        n_atoms = len(atoms)
        data_dict = {
            'positions': torch.tensor(positions, dtype=torch.float32, device=self.device, requires_grad=True),
            'atomic_numbers': torch.tensor(atoms.get_atomic_numbers(), dtype=torch.long, device=self.device),
            'edge_index': torch.tensor(edge_index, dtype=torch.long, device=self.device),
            'shifts': torch.tensor(shifts, dtype=torch.float32, device=self.device),
            'unit_shifts': torch.tensor(unit_shifts, dtype=torch.float32, device=self.device),
            'num_nodes': torch.tensor([n_atoms], dtype=torch.long, device=self.device),
            'ptr': torch.tensor([0, n_atoms], dtype=torch.long, device=self.device),
            'batch': torch.zeros(n_atoms, dtype=torch.long, device=self.device)
        }
        if any(pbc):
            data_dict['cell'] = torch.tensor(cell, dtype=torch.float32, device=self.device).unsqueeze(0)
        else:
            # Provide a dummy large cell for non-periodic systems
            data_dict['cell'] = torch.eye(3, dtype=torch.float32, device=self.device).unsqueeze(0) * 100.0

        # Forward pass
        # We need gradients for forces even if we don't use them, 
        # because the model has a Forces module that calls autograd.grad
        output = self.model(data_dict, training=True)
            
        # Extract charges
        # CACE-LR typically stores charges in the output dict
        # We need to find the right key. Default is 'charge' or 'charges'
        charge_key = None
        for key in ['q', 'charge', 'charges', 'CACE_charge', 'CACE_charges']:
            if key in output:
                charge_key = key
                break
        
        if charge_key is None:
            # Fallback: check if atomwise module has output_key
            # This is model structure dependent
            print("  Warning: No charge key found in output. Available keys:", list(output.keys()))
            return np.zeros(3), np.zeros(n_atoms)
            
        charges_raw = output[charge_key].detach().cpu().numpy()
        
        # Charges are usually (N, 1) or (N,)
        if len(charges_raw.shape) > 1:
            charges_raw = charges_raw.flatten()
            
        # Scale charges to elementary charge e
        charges_e = charges_raw * self.charge_unit
        
        # Compute dipole moment μ = Σ q_i * r_i
        # Units: e * A
        dipole = np.sum(charges_e[:, None] * positions, axis=0)
        
        return dipole, charges_e
