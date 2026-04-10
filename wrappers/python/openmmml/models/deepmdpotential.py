"""
deepmdpotential.py: Implements the DeePMD potential (https://github.com/deepmodeling/deepmd-kit) function using OpenMMDeepmdPlugin (https://github.com/JingHuangLab/openmm_deepmd_plugin). 

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

class DeepmdPotentialImplFactory(MLPotentialImplFactory):
    """This is the factory that creates DeepmdPotentialImpl objects."""

    def createImpl(self, name: str, **args) -> MLPotentialImpl:
        return DeepmdPotentialImpl(name)


class DeepmdPotentialImpl(MLPotentialImpl):
    """Implementation of the DeePMD (Deep Potential Molecular Dynamics) potential for OpenMM.
    
    This class provides an interface to use DeePMD neural network potentials in OpenMM simulations through the OpenMMDeepmdPlugin. 
    DeePMD potentials are trained using the DeepMD-kit framework and can provide highly accurate quantum mechanical-level forces and energies for molecular dynamics simulations.
    
    The DeePMD potential supports various advanced features including:
    - Alchemical free energy calculations with lambda parameters  
    - NNP/MM-style region selection
    - Custom unit conversions between different simulation packages
    
    Prerequisites:
        - OpenMMDeepmdPlugin must be installed:
          `conda install -c conda-forge ye-ding::openmm_deepmd_plugin`
        - A trained DeePMD model file (.pb format from DeepMD-kit)
    
    Usage:
        The DeepmdPotentialImpl is typically used through the MLPotential interface:
        
        .. code-block:: python

            from openmmml import MLPotential
            import openmm.app as app

            # Load your system
            pdb = app.PDBFile('system.pdb')

            # Create DeePMD potential
            potential = MLPotential('deepmd', model='model.pb',
                                    coordinatesCoefficient=10.0,
                                    forceCoefficient=964.8792534459,
                                    energyCoefficient=96.48792534459)

            # Create system with DeePMD forces
            system = potential.createSystem(
                topology=pdb.topology,
                atoms=None,  # Use all atoms, or specify subset
            )

        For alchemical simulations:

        .. code-block:: python

            # Create system with lambda parameter for free energy calculations
            potential = MLPotential('deepmd', model='model.pb')
            system = potential.createSystem(
                topology=pdb.topology,
                lambdaName='lambda_alchemical',
                lambdaValue=0.5,  # Lambda value between 0 and 1
            )

            # During simulation, you can also update lambda value through the context:
            context.setParameter('lambda_alchemical', new_lambda_value)

        For NNP/MM-style simulations with atom subsets:
        
        .. code-block:: python

            forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
            mm_system = forcefield.createSystem(topology)

            # Define NNP region atoms (0-indexed)
            nnp_atoms = [0, 1, 2, 3, 4]  # First 5 atoms

            potential = MLPotential('deepmd', model='model.pb')
            nnp_mm_system = potential.createMixedSystem(
                topology=pdb.topology,
                system=mm_system,
                atoms=nnp_atoms,  # Only these atoms use DeePMD
            )

    Note:
        The default unit conversion coefficients are set up for converting between
        DeePMD-kit's native units (Angstrom, eV) and OpenMM's units (nm, kJ/mol).
        Modify these coefficients if your model uses different units.
    """

    def __init__(self, 
                name, 
                model: str = None, 
                coordinatesCoefficient: float = 10.0,
                forceCoefficient: float = 964.8792534459,
                energyCoefficient: float = 96.48792534459,) -> None:
        try:
            from OpenMMDeepmdPlugin import DeepPotentialModel
            from OpenMMDeepmdPlugin import DeepmdForce
        except ImportError:
            raise ImportError(
                "OpenMMDeepmdPlugin is not installed."
                "Please install it with `conda install -c conda-forge ye-ding::openmm_deepmd_plugin`."
            )
        self.name = name
        self.modelPath = model
        if self.modelPath is None:
            raise ValueError("model must be specified for DeepmdPotentialImpl.")
        
        # Create the DeepPotentialModel object.
        self.dp_model = DeepPotentialModel(self.modelPath)
        self.dp_model.setUnitTransformCoefficients(
            coordinatesCoefficient, forceCoefficient, energyCoefficient
        )

    def addForces(self,
                  topology: openmm.app.Topology,
                  system: openmm.System,
                  atoms: Optional[Iterable[int]],
                  forceGroup: int = 0,
                  lambdaName: Optional[str] = None,
                  lambdaValue: Optional[float] = 1.0,
                  **args) -> None:    
        
        if atoms is not None:
            dp_force = self.dp_model.addParticlesToDPRegion(atoms, topology)
        else:
            atoms_all = [atom.index for atom in topology.atoms()]
            dp_force = self.dp_model.addParticlesToDPRegion(atoms_all, topology)
        
        if lambdaName is not None:
            dp_force.addLambdaParameter(lambdaName, lambdaValue)
        
        # Create the TorchForce and add it to the System.
        dp_force.setForceGroup(forceGroup)
        system.addForce(dp_force)
