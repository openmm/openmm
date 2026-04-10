"""
anipotential.py: Implements the ANI potential function using TorchANI.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2021-2023 Stanford University and the Authors.
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

from openmmml.mlpotential import MLPotential, MLPotentialImpl, MLPotentialImplFactory
import openmm
from typing import Iterable, Optional

class ANIPotentialImplFactory(MLPotentialImplFactory):
    """This is the factory that creates ANIPotentialImpl objects."""

    def createImpl(self, name: str, **args) -> MLPotentialImpl:
        return ANIPotentialImpl(name)


class ANIPotentialImpl(MLPotentialImpl):
    """This is the MLPotentialImpl implementing the ANI potential.

    The potential is implemented using TorchANI to build a PyTorch model.  A
    TorchForce is used to add it to the OpenMM System.  The ANI1ccx and ANI2x
    versions are currently supported.

    Both ANI1ccx and ANI2x are ensembles of eight models.  Averaging over all eight
    models leads to slightly more accurate results than any one individually.  You
    can optionally use only a single model by specifying the modelIndex argument to
    select which one to use.  This leads to a large improvement in speed, at the
    cost of a small decrease in accuracy.
    """

    def __init__(self, name):
        self.name = name


    def addForces(self,
                  topology: openmm.app.Topology,
                  system: openmm.System,
                  atoms: Optional[Iterable[int]],
                  forceGroup: int,
                  implementation: str = 'torchani',
                  modelIndex: Optional[int] = None,
                  **args):
        # Create the TorchANI model.

        try:
            import torchani
        except ImportError as e:
            raise ImportError(f"Failed to import torchani with error: {e}. Make sure torchani is installed.")
        import torch
        import openmmtorch

        # `nnpops` throws error if `periodic_table_index`=False if one passes `species` as `species_to_tensor` from `element`
        _kwarg_dict = {'periodic_table_index': True}
        if self.name == 'ani1ccx':
            model = torchani.models.ANI1ccx(**_kwarg_dict)
        elif self.name == 'ani2x':
            model = torchani.models.ANI2x(**_kwarg_dict)
        else:
            raise ValueError('Unsupported ANI model: '+self.name)
        if modelIndex is not None:
            model = model[modelIndex]

        # Create the PyTorch model that will be invoked by OpenMM.
        includedAtoms = list(topology.atoms())
        if atoms is not None:
            includedAtoms = [includedAtoms[i] for i in atoms]
        species = torch.tensor([[atom.element.atomic_number for atom in includedAtoms]])
        if implementation == 'nnpops':
            try:
                from NNPOps import OptimizedTorchANI
                model = OptimizedTorchANI(model, species)
            except Exception as e:
                print(f"failed to equip `nnpops` with error: {e}")
        elif implementation == "torchani":
            pass # do nothing
        else:
            raise NotImplementedError(f"implementation {implementation} is not supported")

        class ANIForce(torch.nn.Module):

            def __init__(self, model, species, atoms, periodic, implementation):
                super(ANIForce, self).__init__()
                self.model = model
                self.species = torch.nn.Parameter(species, requires_grad=False)
                self.energyScale = torchani.units.hartree2kjoulemol(1)
                self.implementation = implementation
                if atoms is None:
                    self.indices = None
                else:
                    self.indices = torch.tensor(sorted(atoms), dtype=torch.int64)
                if periodic:
                    self.pbc = torch.nn.Parameter(torch.tensor([True, True, True], dtype=torch.bool), requires_grad=False)
                else:
                    self.pbc = None

            def forward(self, positions, boxvectors: Optional[torch.Tensor] = None):
                positions = positions.to(torch.float32)
                if self.indices is not None:
                    positions = positions[self.indices]
                if boxvectors is None:
                    _, energy = self.model((self.species, 10.0*positions.unsqueeze(0)))
                else:
                    boxvectors = boxvectors.to(torch.float32)
                    if self.implementation == "torchani":
                        positions = positions - torch.outer(torch.floor(positions[:,2]/boxvectors[2,2]), boxvectors[2])
                        positions = positions - torch.outer(torch.floor(positions[:,1]/boxvectors[1,1]), boxvectors[1])
                        positions = positions - torch.outer(torch.floor(positions[:,0]/boxvectors[0,0]), boxvectors[0])
                    _, energy = self.model((self.species, 10.0*positions.unsqueeze(0)), cell=10.0*boxvectors, pbc=self.pbc)

                return self.energyScale*energy

        # is_periodic...
        is_periodic = (topology.getPeriodicBoxVectors() is not None) or system.usesPeriodicBoundaryConditions()
        aniForce = ANIForce(model, species, atoms, is_periodic, implementation)

        # Convert it to TorchScript.

        module = torch.jit.script(aniForce)

        # Create the TorchForce and add it to the System.

        force = openmmtorch.TorchForce(module)
        force.setForceGroup(forceGroup)
        force.setUsesPeriodicBoundaryConditions(is_periodic)
        system.addForce(force)
