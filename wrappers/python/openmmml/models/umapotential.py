"""
umapotential.py: Implements the UMA (Universal Molecular Atomistic) potential function.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2024 Stanford University and the Authors.
Authors: Peter Eastman
Contributors:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
import warnings
warnings.warn(
    "umapotential.py (TorchScript backend) is deprecated. "
    "Use the '-pythonforce-batch' variants instead (e.g. 'uma-s-1p1-pythonforce-batch').",
    DeprecationWarning,
    stacklevel=2,
)
import openmm
from openmmml.mlpotential import MLPotential, MLPotentialImpl, MLPotentialImplFactory
from typing import Iterable, Optional, Union
import torch
import numpy as np
from ase import Atoms


class UMAForce(torch.nn.Module):
    """
    A PyTorch module that wraps FAIRChem's UMA model for use with OpenMM-Torch.
    
    This module converts OpenMM positions/box vectors to ASE Atoms objects,
    calls the FAIRChem predict unit, and returns energies in OpenMM units.
    """

    def __init__(
        self,
        predict_unit,
        atomicNumbers: np.ndarray,
        symbols: list,
        atoms: Optional[Iterable[int]],
        task_name: str,
        periodic: bool,
        charge: int,
        spin: int,
    ):
        """
        Parameters
        ----------
        predict_unit : MLIPPredictUnit
            The FAIRChem prediction unit.
        atomicNumbers : np.ndarray
            The atomic numbers.
        symbols : list
            The element symbols.
        atoms : iterable of int
            The indices of the atoms. If ``None``, all atoms are included.
        task_name : str
            The task name.
        periodic : bool
            Whether the system is periodic.
        charge : int
            Total charge.
        spin : int
            Spin multiplicity.
        """
        super(UMAForce, self).__init__()

        self.predict_unit = predict_unit
        self.task_name = task_name
        self.periodic = periodic
        
        # Energy conversion factor: eV to kJ/mol
        self.energyScale = 96.4853
        # Length conversion factor: nm to Angstrom
        self.lengthScale = 10.0

        if atoms is None:
            self.indices = None
        else:
            self.indices = torch.tensor(sorted(atoms), dtype=torch.int64)

        # Store atomic information
        self.register_buffer("atomicNumbers", torch.tensor(atomicNumbers, dtype=torch.long, requires_grad=False))
        self.symbols = symbols
        
        # Store charge and spin as buffers
        self.register_buffer("charge", torch.tensor([charge], dtype=torch.long, requires_grad=False))
        self.register_buffer("spin", torch.tensor([spin], dtype=torch.long, requires_grad=False))

    def forward(
        self, positions: torch.Tensor, boxvectors: Optional[torch.Tensor] = None
    ) -> torch.Tensor:
        """
        Forward pass of the model.

        Parameters
        ----------
        positions : torch.Tensor
            The positions of the atoms in nm.
        boxvectors : torch.Tensor, optional
            The box vectors in nm.

        Returns
        -------
        energy : torch.Tensor
            The predicted energy in kJ/mol.
        """
        # NOTE: This forward method is NOT TorchScript compatible due to ASE/FAIRChem dependencies
        # However, it works with OpenMM-Torch through eager execution mode
        
        # Extract subset of atoms if needed
        if self.indices is not None:
            positions = positions[self.indices]

        # Convert positions from nm to Angstrom
        positions_angstrom = positions * self.lengthScale
        
        # Convert box vectors if periodic
        cell = None
        pbc = [False, False, False]
        if boxvectors is not None:
            cell = boxvectors * self.lengthScale
            pbc = [True, True, True] if self.periodic else [False, False, False]

        # Create ASE Atoms object
        from ase import Atoms
        atoms_ase = Atoms(
            symbols=self.symbols,
            positions=positions_angstrom.detach().cpu().numpy(),
            cell=cell.detach().cpu().numpy() if cell is not None else None,
            pbc=pbc,
        )

        # Set charge and spin in atoms.info for omol tasks
        atoms_ase.info["charge"] = int(self.charge.item())
        atoms_ase.info["spin"] = int(self.spin.item())

        # Convert ASE Atoms to AtomicData
        from fairchem.core.datasets.atomic_data import AtomicData
        atomic_data = AtomicData.from_ase(
            atoms_ase,
            task_name=self.task_name,
            r_edges=False,  # Let the model handle edge generation
            r_data_keys=["spin", "charge"],
        )

        # Predict energy using the predict unit
        pred = self.predict_unit.predict(atomic_data)

        # Extract energy and convert from eV to kJ/mol
        energy_ev = pred["energy"]
        energy_kj_mol = energy_ev * self.energyScale

        return energy_kj_mol


class UMAPotentialImplFactory(MLPotentialImplFactory):
    """This is the factory that creates UMAPotentialImpl objects."""

    def createImpl(self, name: str, **args) -> MLPotentialImpl:
        return UMAPotentialImpl(name)


class UMAPotentialImpl(MLPotentialImpl):
    """This is the MLPotentialImpl implementing the UMA (Universal Molecular Atomistic) potential.

    The UMA potential uses FAIRChem's pre-trained models from the facebook/UMA HuggingFace repository
    and integrates them into the OpenMM System using a TorchForce.

    To use a UMA model, specify the model name. For example:

    >>> potential = MLPotential('uma-s-1p1')

    Available UMA models include:
    - 'uma-s-1': UMA small model version 1
    - 'uma-s-1p1': UMA small model version 1.1
    - 'uma-m-1p1': UMA medium model version 1.1
    - 'esen-md-direct-all-omol': eSEN model for molecules (OMol25)
    - 'esen-sm-conserving-all-omol': eSEN small model for molecules (energy conserving)
    - 'esen-sm-direct-all-omol': eSEN small model for molecules
    - 'esen-sm-conserving-all-oc25': eSEN small model for catalysis (OC25)
    - 'esen-md-direct-all-oc25': eSEN medium model for catalysis (OC25)

    During system creation, you must specify the task name using the ``task_name`` keyword argument.
    Supported tasks are 'omol', 'omat', 'oc20', 'odac', and 'omc'. For example:

    >>> system = potential.createSystem(topology, task_name='omol')

    For 'omol' tasks (molecular systems), you can optionally specify charge and spin multiplicity:

    >>> system = potential.createSystem(topology, task_name='omol', charge=0, spin=1)

    You can also control inference settings for performance tuning:

    >>> system = potential.createSystem(topology, task_name='omol', inference_settings='turbo')

    Supported inference_settings are:
    - 'default': General purpose (slower but works with changing composition)
    - 'turbo': Optimized for speed (requires fixed atomic composition)

    Attributes
    ----------
    name : str
        The name of the UMA model.
    """

    def __init__(self, name: str) -> None:
        """
        Initialize the UMAPotentialImpl.

        Parameters
        ----------
        name : str
            The name of the UMA model.
        """
        self.name = name

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
        """
        Add the UMAForce to the OpenMM System.

        Parameters
        ----------
        topology : openmm.app.Topology
            The topology of the system.
        system : openmm.System
            The system to which the force will be added.
        atoms : iterable of int
            The indices of the atoms to include in the model. If ``None``, all atoms are included.
        forceGroup : int
            The force group to which the force should be assigned.
        task_name : str, optional
            The task name for the UMA model. Options are 'omol', 'omat', 'oc20', 'odac', 'omc'.
            If not provided, will be auto-detected from the model.
        inference_settings : str or InferenceSettings, optional
            Inference settings. Can be 'default' or 'turbo'. Default is 'default'.
        charge : int, optional
            Total system charge (for omol tasks). Default is 0.
        spin : int, optional
            Spin multiplicity (for omol tasks). Default is 1.
        device : str, optional
            Device to use for inference ('cuda' or 'cpu'). If None, uses default.
        """
        import torch
        import openmmtorch
        import numpy as np

        try:
            from fairchem.core.calculate import pretrained_mlip
            from fairchem.core.datasets.atomic_data import AtomicData
        except ImportError as e:
            raise ImportError(
                f"Failed to import fairchem with error: {e}. "
                "Install fairchem-core with 'pip install fairchem-core'."
            )

        try:
            from ase import Atoms
        except ImportError as e:
            raise ImportError(
                f"Failed to import ase with error: {e}. "
                "Install ase with 'pip install ase'."
            )

        # Load the UMA model using FAIRChem's pretrained model loader
        # Note: device="cpu" here because OpenMM-Torch handles device placement
        predict_unit = pretrained_mlip.get_predict_unit(
            self.name,
            inference_settings=inference_settings,
            device="cpu" if device is None else device,
        )

        # Determine task name
        if task_name is None:
            # Auto-detect from model's available tasks
            valid_datasets = list(predict_unit.dataset_to_tasks.keys())
            if len(valid_datasets) == 1:
                task_name = valid_datasets[0]
            else:
                raise ValueError(
                    f"Multiple tasks available for model {self.name}: {valid_datasets}. "
                    "Please specify task_name explicitly."
                )

        # Validate task name
        valid_datasets = list(predict_unit.dataset_to_tasks.keys())
        if task_name not in valid_datasets:
            raise ValueError(
                f"Invalid task_name '{task_name}' for model {self.name}. "
                f"Valid options are: {valid_datasets}"
            )

        # Get the atomic numbers of the ML region
        includedAtoms = list(topology.atoms())
        if atoms is not None:
            includedAtoms = [includedAtoms[i] for i in atoms]
        atomicNumbers = np.array([atom.element.atomic_number for atom in includedAtoms])
        symbols = [atom.element.symbol for atom in includedAtoms]

        isPeriodic = (
            topology.getPeriodicBoxVectors() is not None
        ) or system.usesPeriodicBoundaryConditions()

        umaForce = UMAForce(
            predict_unit,
            atomicNumbers,
            symbols,
            atoms,
            task_name,
            isPeriodic,
            charge,
            spin,
        )

        # NOTE: UMAForce uses ASE/FAIRChem which aren't TorchScript compatible
        # We need to trace it with an example input to create a TorchScript module
        import tempfile
        import os
        
        # Create example inputs for tracing
        n_atoms = len(atomicNumbers)
        example_positions = torch.zeros(n_atoms, 3, dtype=torch.float64)
        example_boxvectors = torch.eye(3, dtype=torch.float64) * 10.0 if isPeriodic else None
        
        # Trace the model with example inputs
        if isPeriodic:
            traced_model = torch.jit.trace(umaForce, (example_positions, example_boxvectors))
        else:
            traced_model = torch.jit.trace(umaForce, (example_positions,))
        
        # Save the traced model to a temporary file
        with tempfile.NamedTemporaryFile(mode='wb', suffix='.pt', delete=False) as f:
            model_file = f.name
            torch.jit.save(traced_model, f)
        
        # Create the TorchForce from the saved model file
        force = openmmtorch.TorchForce(model_file)
        
        # Clean up the temporary file
        os.unlink(model_file)
        force.setForceGroup(forceGroup)
        force.setUsesPeriodicBoundaryConditions(isPeriodic)
        system.addForce(force)
