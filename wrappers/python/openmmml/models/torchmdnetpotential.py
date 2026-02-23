"""
torchmdnetpotential.py: Implements the TorchMD-Net potential function.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2021-2025 Stanford University and the Authors.
Authors: Peter Eastman
Contributors: Stephen Farr

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
from openmmml.mlpotential import MLPotentialImpl, MLPotentialImplFactory
from typing import Iterable, Optional

class TorchMDNetPotentialImplFactory(MLPotentialImplFactory):
    """This is the factory that creates TorchMDNetPotentialImpl objects."""

    def createImpl(
        self, 
        name: str, 
        modelPath: Optional[str] = None,
        lengthScale: float = 0.1, # angstrom -> nm
        energyScale: float = 96.4916,  # eV -> kJ/mol
    ) -> MLPotentialImpl:
        return TorchMDNetPotentialImpl(name, modelPath, lengthScale, energyScale)

class TorchMDNetPotentialImpl(MLPotentialImpl):
    """This is the MLPotentialImpl implementing the TorchMDNet potential.

    The TorchMDNet potential is constructed using `torchmdnet` to build a PyTorch model,
    and then integrated into the OpenMM System using a TorchForce.  To use it, specify the model by name
    and provide the path to a model.

    >>> potential = MLPotential('torchmdnet', modelPath=<model_file_path>)

    The default energy and length scales assume a model is trained with positions in angstroms and energies in eV.
    If this is not the case you can specify the length and energy scales by passing the factors that convert the model
    distance to nm and the energy to kJ/mol, for example:

    >>> potential = MLPotential('torchmdnet', modelPath=<model_file_path>, 
                                lengthScale=0.1 # angstrom to nm, 
                                energyScale=4.184 # kcal/mol to kJ/mol)

    During system creation you can enable CUDA graphs for a speed-up for small molecules:

    >>>  system = potential.createSystem(pdb.topology, cudaGraphs=True)

    The default is to enable this for TensorNet models.

    You can also specify the molecule's total charge:

    >>>  system = potential.createSystem(pdb.topology, charge=0)

    Pretained AceFF models can be used directly:

    >>> potential = MLPotential('aceff-2.0')

    >>> potential = MLPotential('aceff-1.1')

    >>> potential = MLPotential('aceff-1.0')

    """

    def __init__(self, 
                 name: str,
                 modelPath: str,
                 lengthScale: float,
                 energyScale: float
    ) -> None:
        """
        Initialize the TorchMDNetPotentialImpl.

        Parameters
        ----------
        name : str
            The name of the model.
            'torchmdnet' for a local model file, or pretrained models are available: 'aceff-1.0' or 'aceff-1.1'. 
        modelPath : str, optional
            The path to the locally trained torchmdnet model if ``name`` is 'torchmdnet'.
        lengthScale : float
            The length conversion factor from the model units to nanometers. 
            If not specified the default is 0.1 which corresponds to a model in angstrom
        energyScale : float
            The energy conversion factor from the model units to kJ/mol.
            If not specified the default is 96.4916 which corresponds to a model in eV.
   
        """
        self.name = name
        self.modelPath = modelPath
        self.lengthScale = lengthScale
        self.energyScale = energyScale

    def addForces(self,
                  topology: openmm.app.Topology,
                  system: openmm.System,
                  atoms: Optional[Iterable[int]],
                  forceGroup: int,
                  **args):
        # Load the TorchMDNet model.
        try:
            import torchmdnet
            from torchmdnet.models.model import load_model
        except ImportError as e:
            raise ImportError(f"Failed to import torchmdnet please install from https://torchmd-net.readthedocs.io/en/latest/installation.html")
        import torch
        import openmmtorch

        includedAtoms = list(topology.atoms())
        if atoms is not None:
            includedAtoms = [includedAtoms[i] for i in atoms]
        numbers = torch.tensor([atom.element.atomic_number for atom in includedAtoms])
        charge = torch.tensor([args.get('charge', 0)], dtype=torch.float32)

        if self.name == 'torchmdnet':
            # a local path to a torchmdnet checkpoint must be provided 
            model_file_path = self.modelPath
        else:
            try:
                from huggingface_hub import hf_hub_download
            except ImportError as e:
                raise ImportError(f"Failed to import huggingface_hub please install from https://huggingface.co/docs/huggingface_hub/en/installation")
            
            if self.name == 'aceff-1.0':
                repo_id="Acellera/AceFF-1.0"
                filename="aceff_v1.0.ckpt"
            elif self.name == 'aceff-1.1':
                repo_id="Acellera/AceFF-1.1"
                filename="aceff_v1.1.ckpt"
            elif self.name == 'aceff-2.0':
                repo_id="Acellera/AceFF-2.0"
                filename="aceff_v2.0.ckpt"
            else:
                raise ValueError(f'Model name {self.name} does not exist.')

            model_file_path = hf_hub_download(
                repo_id=repo_id,
                filename=filename,
            )

        model = load_model(
            model_file_path,
            derivative=False,
            remove_ref_energy = args.get('remove_ref_energy', True),
            max_num_neighbors = min(args.get('max_num_neighbors', 64), numbers.shape[0]),
            static_shapes = True,
            check_errors = False
        )
        for parameter in model.parameters():
            parameter.requires_grad = False

        batch = args.get('batch', None)
        if batch is None:
            batch = torch.zeros_like(numbers)
        else:
            batch = torch.tensor(batch, dtype=torch.long)

        # TensorNet models can use CUDA graphs and the default is to use them.
        use_cudagraphs = args.get('cudaGraphs', 
                                True if (isinstance(model.representation_model, torchmdnet.models.tensornet.TensorNet)
                                or isinstance(model.representation_model, torchmdnet.models.tensornet2.TensorNet2))
                                else False
                            )

        class TorchMDNetForce(torch.nn.Module):
            def __init__(self, model, numbers, charge, atoms, batch, lengthScale, energyScale):
                super(TorchMDNetForce, self).__init__()
                self.model = model
                self.numbers = torch.nn.Parameter(numbers, requires_grad=False)
                self.charge = torch.nn.Parameter(charge, requires_grad=False)
                self.batch = torch.nn.Parameter(batch, requires_grad=False)
                self.lengthScale = lengthScale
                self.energyScale = energyScale
                if atoms is None:
                    self.subset = False
                    self.indices = torch.empty(1) # for torchscript
                else:
                    self.subset = True
                    self.indices = torch.nn.Parameter(torch.tensor(sorted(atoms), dtype=torch.int64), requires_grad=False)
                    
            def forward(self, positions: torch.Tensor, boxvectors: Optional[torch.Tensor] = None):
                positions = positions.to(torch.float32).to(self.numbers.device)
                if self.subset:
                    positions = positions[self.indices]

                energy = self.model(z=self.numbers, pos=positions/self.lengthScale, batch=self.batch, q=self.charge)[0]
                return energy*self.energyScale

        # Create the TorchForce and add it to the System.
        module = torch.jit.script(TorchMDNetForce(model, numbers, charge, atoms, batch, self.lengthScale, self.energyScale)).to(torch.device('cpu'))
        force = openmmtorch.TorchForce(module)
        if use_cudagraphs:
            force.setProperty("useCUDAGraphs", "true")
        force.setForceGroup(forceGroup)
        system.addForce(force)
