"""
macepotential.py: Implements the MACE potential function.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2021 Stanford University and the Authors.
Authors: Peter Eastman
Contributors: Stephen Farr, Joao Morado

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
import openmm
from openmmml.mlpotential import MLPotential, MLPotentialImpl, MLPotentialImplFactory
from typing import Iterable, Optional, Tuple


class MACEPotentialImplFactory(MLPotentialImplFactory):
    """This is the factory that creates MACEPotentialImpl objects."""

    def createImpl(
        self, name: str, modelPath: Optional[str] = None, **args
    ) -> MLPotentialImpl:
        return MACEPotentialImpl(name, modelPath)


class MACEPotentialImpl(MLPotentialImpl):
    """This is the MLPotentialImpl implementing the MACE potential.

    The MACE potential is constructed using MACE to build a PyTorch model,
    and then integrated into the OpenMM System using a TorchForce.
    This implementation supports both foundation models and locally trained MACE models.

    To use one of the pre-trained MACE foundation models, specify the model name. For example:

    >>> potential = MLPotential('mace-off23-small')

    Other available models include 'mace-off23-medium', 'mace-off23-large', and 'mace-off24-medium'.

    To use a locally trained MACE model, provide the path to the model file. For example:

    >>> potential = MLPotential('mace', modelPath='MACE.model')

    During system creation, you can optionally specify the precision of the model using the
    ``precision`` keyword argument. Supported options are 'single' and 'double'. For example:

    >>> system = potential.createSystem(topology, precision='single')

    By default, the implementation uses the precision of the loaded MACE model.
    According to the MACE documentation, 'single' precision is recommended for MD (faster but
    less accurate), while 'double' precision is recommended for geometry optimization.

    Additionally, you can request computation of the full atomic energy, including the atom
    self-energy, instead of the default interaction energy, by setting ``returnEnergyType`` to
    'energy'. For example:
    
    >>> system = potential.createSystem(topology, returnEnergyType='energy')

    The default is to compute the interaction energy, which can be made explicit by setting
    ``returnEnergyType='interaction_energy'``.

    Attributes
    ----------
    name : str
        The name of the MACE model.
    modelPath : str
        The path to the locally trained MACE model if ``name`` is 'mace'.
    """

    def __init__(self, name: str, modelPath) -> None:
        """
        Initialize the MACEPotentialImpl.

        Parameters
        ----------
        name : str
            The name of the MACE model.
            Options include 'mace-off23-small', 'mace-off23-medium', 'mace-off23-large',
            'mace-off24-medium', and 'mace'.
        modelPath : str, optional
            The path to the locally trained MACE model if ``name`` is 'mace'.
        """
        self.name = name
        self.modelPath = modelPath

    def addForces(
        self,
        topology: openmm.app.Topology,
        system: openmm.System,
        atoms: Optional[Iterable[int]],
        forceGroup: int,
        precision: Optional[str] = None,
        returnEnergyType: str = "interaction_energy",
        **args,
    ) -> None:
        """
        Add the MACEForce to the OpenMM System.

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
        precision : str, optional
            The precision of the model. Supported options are 'single' and 'double'.
            If ``None``, the default precision of the model is used.
        returnEnergyType : str, optional
            Whether to return the interaction energy or the energy including the self-energy.
            Default is 'interaction_energy'. Supported options are 'interaction_energy' and 'energy'.
        """
        import torch
        import openmmtorch

        try:
            from mace.tools import utils, to_one_hot, atomic_numbers_to_indices
            from mace.calculators.foundations_models import mace_off, mace_mp, mace_omol
        except ImportError as e:
            raise ImportError(
                f"Failed to import mace with error: {e}. "
                "Install mace with 'pip install mace-torch'."
            )
        try:
            from e3nn.util import jit
        except ImportError as e:
            raise ImportError(
                f"Failed to import e3nn with error: {e}. "
                "Install e3nn with 'pip install e3nn'."
            )
        try:
            from NNPOps.neighbors import getNeighborPairs
        except ImportError as e:
            raise ImportError(
                f"Failed to import NNPOps with error: {e}. "
                "Install NNPOps with 'conda install -c conda-forge nnpops'."
            )

        assert returnEnergyType in [
            "interaction_energy",
            "energy",
        ], f"Unsupported returnEnergyType: '{returnEnergyType}'. Supported options are 'interaction_energy' or 'energy'."

        models = {
            'mace-off23-small': (mace_off, 'small', True),
            'mace-off23-medium': (mace_off, 'medium', True),
            'mace-off23-large': (mace_off, 'large', True),
            'mace-off24-medium': (mace_off, 'https://github.com/ACEsuit/mace-off/blob/main/mace_off24/MACE-OFF24_medium.model?raw=true', True),
            'mace-mpa-0-medium': (mace_mp, 'medium-mpa-0', False),
            'mace-omat-0-small': (mace_mp, 'small-omat-0', True),
            'mace-omat-0-medium': (mace_mp, 'medium-omat-0', True),
            'mace-omol-0-extra-large': (mace_omol, 'extra_large', True)
        }

        # Load the model to the CPU (OpenMM-Torch takes care of loading to the right devices)
        if self.name in models:
            fn, name, warn = models[self.name]
            model = fn(model=name, device="cpu", return_raw_model=True)
            if warn:
                import logging
                logging.warning(f'The model {self.name} is distributed under the restrictive ASL license.  Commercial use is not permitted.')
        elif self.name == "mace":
            if self.modelPath is not None:
                model = torch.load(self.modelPath, map_location="cpu")
            else:
                raise ValueError("No modelPath provided for local MACE model.")
        else:
            raise ValueError(f"Unsupported MACE model: {self.name}")

        # Compile the model.
        model = jit.compile(model)  

        # Get the atomic numbers of the ML region.
        includedAtoms = list(topology.atoms())
        if atoms is not None:
            includedAtoms = [includedAtoms[i] for i in atoms]
        atomicNumbers = [atom.element.atomic_number for atom in includedAtoms]

        # Set the precision that the model will be used with.
        modelDefaultDtype = next(model.parameters()).dtype
        if precision is None:
            dtype = modelDefaultDtype
        elif precision == "single":
            dtype = torch.float32
        elif precision == "double":
            dtype = torch.float64
        else:
            raise ValueError(
                f"Unsupported precision {precision} for the model. "
                "Supported values are 'single' and 'double'."
            )
        if dtype != modelDefaultDtype:
            print(
                f"Model dtype is {modelDefaultDtype} "
                f"and requested dtype is {dtype}. "
                "The model will be converted to the requested dtype."
            )

        # One hot encoding of atomic numbers
        zTable = utils.AtomicNumberTable([int(z) for z in model.atomic_numbers])
        nodeAttrs = to_one_hot(
            torch.tensor(
                atomic_numbers_to_indices(atomicNumbers, z_table=zTable),
                dtype=torch.long,
            ).unsqueeze(-1),
            num_classes=len(zTable),
        )

        class MACEForce(torch.nn.Module):
            """
            MACEForce class to be used with TorchForce.

            Parameters
            ----------
            model : torch.jit._script.RecursiveScriptModule
                The compiled MACE model.
            dtype : torch.dtype
                The precision with which the model will be used.
            energyScale : float
                Conversion factor for the energy, viz. eV to kJ/mol.
            lengthScale : float
                Conversion factor for the length, viz. nm to Angstrom.
            indices : torch.Tensor
                The indices of the atoms to calculate the energy for.
            charge : float
                Total charge of the system
            multiplicity : float
                Spin multiplicity of the system
            returnEnergyType : str
                Whether to return the interaction energy or the energy including the self-energy.
            inputDict : dict
                The input dictionary passed to the model.
            """

            def __init__(
                self,
                model: torch.jit._script.RecursiveScriptModule,
                nodeAttrs: torch.Tensor,
                atoms: Optional[Iterable[int]],
                charge: float,
                multiplicity: float,
                periodic: bool,
                dtype: torch.dtype,
                returnEnergyType: str,
            ) -> None:
                """
                Initialize the MACEForce.

                Parameters
                ----------
                model : torch.jit._script.RecursiveScriptModule
                    The MACE model.
                nodeAttrs : torch.Tensor
                    The one-hot encoded atomic numbers.
                atoms : iterable of int
                    The indices of the atoms. If ``None``, all atoms are included.
                periodic : bool
                    Whether the system is periodic.
                dtype : torch.dtype
                    The precision of the model.
                returnEnergyType : str
                    Whether to return the interaction energy or the energy including the self-energy.
                """
                super(MACEForce, self).__init__()

                self.dtype = dtype
                self.model = model.to(self.dtype)
                self.energyScale = 96.4853
                self.lengthScale = 10.0
                self.returnEnergyType = returnEnergyType

                if atoms is None:
                    self.indices = None
                else:
                    self.indices = torch.tensor(sorted(atoms), dtype=torch.int64)
                
                # Create the default input dict.
                self.register_buffer("ptr", torch.tensor([0, nodeAttrs.shape[0]], dtype=torch.long, requires_grad=False))
                self.register_buffer("node_attrs", nodeAttrs.to(self.dtype))
                self.register_buffer("batch", torch.zeros(nodeAttrs.shape[0], dtype=torch.long, requires_grad=False))
                self.register_buffer("pbc", torch.tensor([periodic, periodic, periodic], dtype=torch.bool, requires_grad=False))
                self.register_buffer("charge", torch.tensor([charge], dtype=dtype, requires_grad=False))
                self.register_buffer("multiplicity", torch.tensor([multiplicity], dtype=dtype, requires_grad=False))

            def _getNeighborPairs(
                self, positions: torch.Tensor, cell: Optional[torch.Tensor]
            ) -> Tuple[torch.Tensor, torch.Tensor]:
                """
                Get the shifts and edge indices.

                Notes
                -----
                This method calculates the shifts and edge indices by determining neighbor pairs (``neighbors``)
                and respective wrapped distances (``wrappedDeltas``) using ``NNPOps.neighbors.getNeighborPairs``.
                After obtaining the ``neighbors`` and ``wrappedDeltas``, the pairs with negative indices (r>cutoff)
                are filtered out, and the edge indices and shifts are finally calculated.

                Parameters
                ----------
                positions : torch.Tensor
                    The positions of the atoms.
                cell : torch.Tensor
                    The cell vectors.

                Returns
                -------
                edgeIndex : torch.Tensor
                    The edge indices.
                shifts : torch.Tensor
                    The shifts.
                """
                # Get the neighbor pairs, shifts and edge indices.
                neighbors, wrappedDeltas, _, _ = getNeighborPairs(
                    positions, self.model.r_max, -1, cell
                )
                mask = neighbors >= 0
                neighbors = neighbors[mask].view(2, -1)
                wrappedDeltas = wrappedDeltas[mask[0], :]

                edgeIndex = torch.hstack((neighbors, neighbors.flip(0))).to(torch.int64)
                if cell is not None:
                    deltas = positions[edgeIndex[0]] - positions[edgeIndex[1]]
                    wrappedDeltas = torch.vstack((wrappedDeltas, -wrappedDeltas))
                    shiftsIdx = torch.mm(deltas - wrappedDeltas, torch.linalg.inv(cell))
                    shifts = torch.mm(shiftsIdx, cell)
                else:
                    shifts = torch.zeros((edgeIndex.shape[1], 3), dtype=self.dtype, device=positions.device)

                return edgeIndex, shifts

            def forward(
                self, positions: torch.Tensor, boxvectors: Optional[torch.Tensor] = None
            ) -> torch.Tensor:
                """
                Forward pass of the model.

                Parameters
                ----------
                positions : torch.Tensor
                    The positions of the atoms.
                box_vectors : torch.Tensor
                    The box vectors.

                Returns
                -------
                energy : torch.Tensor
                    The predicted energy in kJ/mol.
                """
                # Setup positions and cell.
                if self.indices is not None:
                    positions = positions[self.indices]

                positions = positions.to(self.dtype) * self.lengthScale

                if boxvectors is not None:
                    cell = boxvectors.to(self.dtype) * self.lengthScale
                else:
                    cell = None

                # Get the shifts and edge indices.
                edgeIndex, shifts = self._getNeighborPairs(positions, cell)

                # Update input dictionary.
                inputDict = {
                    "ptr": self.ptr,
                    "node_attrs": self.node_attrs,
                    "batch": self.batch,
                    "pbc": self.pbc,
                    "positions": positions,
                    "edge_index": edgeIndex,
                    "shifts": shifts,
                    "cell": cell if cell is not None else torch.zeros(3, 3, dtype=self.dtype),
                    "total_charge": self.charge,
                    "total_spin": self.multiplicity
                }

                # Predict the energy.
                energy = self.model(inputDict, compute_force=False)[
                    self.returnEnergyType
                ]

                assert (
                    energy is not None
                ), "The model did not return any energy. Please check the input."

                return energy * self.energyScale

        isPeriodic = (
            topology.getPeriodicBoxVectors() is not None
        ) or system.usesPeriodicBoundaryConditions()

        maceForce = MACEForce(
            model,
            nodeAttrs,
            atoms,
            float(args.get('charge', 0)),
            float(args.get('multiplicity', 1)),
            isPeriodic,
            dtype,
            returnEnergyType,
        )

        # Convert it to TorchScript.
        module = torch.jit.script(maceForce)

        # Create the TorchForce and add it to the System.
        force = openmmtorch.TorchForce(module)
        force.setForceGroup(forceGroup)
        force.setUsesPeriodicBoundaryConditions(isPeriodic)
        system.addForce(force)
