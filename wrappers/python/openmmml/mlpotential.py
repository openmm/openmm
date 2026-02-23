"""
mlpotential.py: Provides a common API for creating OpenMM Systems with ML potentials.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2021-2024 Stanford University and the Authors.
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

import openmm
import openmm.app
import openmm.unit as unit
from copy import deepcopy
from typing import Dict, Iterable, Optional
import sys
if sys.version_info < (3, 10):
    from importlib_metadata import entry_points
else:
    from importlib.metadata import entry_points


class MLPotentialImplFactory(object):
    """Abstract interface for classes that create MLPotentialImpl objects.

    If you are defining a new potential function, you need to create subclasses
    of MLPotentialImpl and MLPotentialImplFactory, and register an instance of
    the factory by calling MLPotential.registerImplFactory().  Alternatively,
    if a Python package creates an entry point in the group "openmmml.potentials",
    the potential will be registered automatically.  The entry point name is the
    name of the potential function, and the value should be the name of the
    MLPotentialImplFactory subclass.
    """
    
    def createImpl(self, name: str, **args) -> "MLPotentialImpl":
        """Create a MLPotentialImpl that will be used to implement a MLPotential.

        When a MLPotential is created, it invokes this method to create an object
        implementing the requested potential.  Subclasses must implement this method
        to return an instance of the correct MLPotentialImpl subclass.

        Parameters
        ----------
        name: str
            the name of the potential that was specified to the MLPotential constructor
        args:
            any additional keyword arguments that were provided to the MLPotential
            constructor are passed to this method.  This allows subclasses to customize
            their behavior based on extra arguments.

        Returns
        -------
        a MLPotentialImpl that implements the potential
        """
        raise NotImplementedError('Subclasses must implement createImpl()')


class MLPotentialImpl(object):
    """Abstract interface for classes that implement potential functions.

    If you are defining a new potential function, you need to create subclasses
    of MLPotentialImpl and MLPotentialImplFactory.  When a user creates a
    MLPotential and specifies a name for the potential to use, it looks up the
    factory that has been registered for that name and uses it to create a
    MLPotentialImpl of the appropriate subclass.
    """
    
    def addForces(self,
                  topology: openmm.app.Topology,
                  system: openmm.System,
                  atoms: Optional[Iterable[int]],
                  forceGroup: int,
                  **args):
        """Add Force objects to a System to implement the potential function.

        This is invoked by MLPotential.createSystem().  Subclasses must implement
        it to create the requested potential function.

        Parameters
        ----------
        topology: Topology
            the Topology from which the System is being created
        system: System
            the System that is being created
        atoms: Optional[Iterable[int]]
            the indices of atoms the potential should be applied to, or None if
            it should be applied to the entire System
        forceGroup: int
            the force group that any newly added Forces should be in
        args:
            any additional keyword arguments that were provided to createSystem()
            are passed to this method.  This allows subclasses to customize their
            behavior based on extra arguments.
        """
        raise NotImplementedError('Subclasses must implement addForces()')


class MLPotential(object):
    """A potential function that can be used in simulations.

    To use this class, create a MLPotential, specifying the name of the potential
    function to use.  You can then call createSystem() to create a System object
    for a simulation.  For example,

    >>> potential = MLPotential('ani2x')
    >>> system = potential.createSystem(topology)

    Alternatively, you can use createMixedSystem() to create a System where part is
    modeled with this potential and the rest is modeled with a conventional force
    field.  As an example, suppose the Topology contains three chains.  Chain 0 is
    a protein, chain 1 is a ligand, and chain 2 is solvent.  The following code
    creates a System in which the internal energy of the ligand is computed with
    ANI2x, while everything else (including interactions between the ligand and the
    rest of the System) is computed with Amber14.

    >>> forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    >>> mm_system = forcefield.createSystem(topology)
    >>> chains = list(topology.chains())
    >>> ml_atoms = [atom.index for atom in chains[1].atoms()]
    >>> potential = MLPotential('ani2x')
    >>> ml_system = potential.createMixedSystem(topology, mm_system, ml_atoms)
    """

    _implFactories: Dict[str, MLPotentialImplFactory] = {}
    
    def __init__(self, name: str, **args):
        """Create a MLPotential.

        Parameters
        ----------
        name: str
            the name of the potential function to use.  Built in support is currently
            provided for the following: 'ani1ccx', 'ani2x'.  Others may be added by
            calling MLPotential.registerImplFactory().
        args:
            particular potential functions may define additional arguments that can
            be used to customize them.  See the documentation on the specific
            potential functions for more information.
        """
        self._impl = MLPotential._implFactories[name].createImpl(name, **args)

    def createSystem(self, topology: openmm.app.Topology, removeCMMotion: bool = True, **args) -> openmm.System:
        """Create a System for running a simulation with this potential function.

        Parameters
        ----------
        topology: Topology
            the Topology for which to create a System
        removeCMMotion: bool
            if true, a CMMotionRemover will be added to the System. 
        args:
            particular potential functions may define additional arguments that can
            be used to customize them.  See the documentation on the specific
            potential functions for more information.

        Returns
        -------
        a newly created System object that uses this potential function to model the Topology
        """
        system = openmm.System()
        if topology.getPeriodicBoxVectors() is not None:
            system.setDefaultPeriodicBoxVectors(*topology.getPeriodicBoxVectors())
        for atom in topology.atoms():
            if atom.element is None:
                system.addParticle(0)
            else:
                system.addParticle(atom.element.mass)
        self._impl.addForces(topology, system, None, 0, **args)
        if removeCMMotion:
            system.addForce(openmm.CMMotionRemover())
        return system

    def createMixedSystem(self,
                          topology: openmm.app.Topology,
                          system: openmm.System,
                          atoms: Iterable[int],
                          removeConstraints: bool = True,
                          forceGroup: int = 0,
                          interpolate: bool = False,
                          **args) -> openmm.System:
        """Create a System that is partly modeled with this potential and partly
        with a conventional force field.

        To use this method, first create a System that is entirely modeled with the
        conventional force field.  Pass it to this method, along with the indices of the
        atoms to model with this potential (the "ML subset").  It returns a new System
        that is identical to the original one except for the following changes.

        1. Removing all bonds, angles, and torsions for which *all* atoms are in the
           ML subset.
        2. For every NonbondedForce and CustomNonbondedForce, adding exceptions/exclusions
           to prevent atoms in the ML subset from interacting with each other.
        3. (Optional) Removing constraints between atoms that are both in the ML subset.
        4. Adding Forces as necessary to compute the internal energy of the ML subset
           with this potential.

        Alternatively, the System can include Forces to compute the energy both with the
        conventional force field and with this potential, and to smoothly interpolate
        between them.  In that case, it creates a CustomCVForce containing the following.

        1. The Forces to compute this potential.
        2. Forces to compute the bonds, angles, and torsions that were removed above.
        3. For every NonbondedForce, a corresponding CustomBondForce to compute the
           nonbonded interactions within the ML subset.

        The CustomCVForce defines a global parameter called "lambda_interpolate" that interpolates
        between the two potentials.  When lambda_interpolate=0, the energy is computed entirely with
        the conventional force field.  When lambda_interpolate=1, the energy is computed entirely with
        the ML potential.  You can set its value by calling setParameter() on the Context.

        Parameters
        ----------
        topology: Topology
            the Topology for which to create a System
        system: System
            a System that models the Topology with a conventional force field
        atoms: Iterable[int]
            the indices of all atoms whose interactions should be computed with
            this potential
        removeConstraints: bool
            if True, remove constraints between pairs of atoms whose interaction
            will be computed with this potential
        forceGroup: int
            the force group the ML potential's Forces should be placed in
        interpolate: bool
            if True, create a System that can smoothly interpolate between the conventional
            and ML potentials
        args:
            particular potential functions may define additional arguments that can
            be used to customize them.  See the documentation on the specific
            potential functions for more information.

        Returns
        -------
        a newly created System object that uses this potential function to model the Topology
        """
        # Create the new System, removing bonded interactions within the ML subset.

        newSystem = self._removeBonds(system, atoms, True, removeConstraints)

        # Add nonbonded exceptions and exclusions.

        atomList = list(atoms)
        for force in newSystem.getForces():
            if isinstance(force, openmm.NonbondedForce):
                for i in range(len(atomList)):
                    for j in range(i):
                        force.addException(atomList[i], atomList[j], 0, 1, 0, True)
            elif isinstance(force, openmm.CustomNonbondedForce):
                existing = set(tuple(force.getExclusionParticles(i)) for i in range(force.getNumExclusions()))
                for i in range(len(atomList)):
                    a1 = atomList[i]
                    for j in range(i):
                        a2 = atomList[j]
                        if (a1, a2) not in existing and (a2, a1) not in existing:
                            force.addExclusion(a1, a2)

        # Add the ML potential.

        if not interpolate:
            self._impl.addForces(topology, newSystem, atomList, forceGroup, **args)
        else:
            # Create a CustomCVForce and put the ML forces inside it.

            cv = openmm.CustomCVForce('')
            cv.addGlobalParameter('lambda_interpolate', 1)
            tempSystem = openmm.System()
            self._impl.addForces(topology, tempSystem, atomList, forceGroup, **args)
            mlVarNames = []
            for i, force in enumerate(tempSystem.getForces()):
                name = f'mlForce{i+1}'
                cv.addCollectiveVariable(name, deepcopy(force))
                mlVarNames.append(name)

            # Create Forces for all the bonded interactions within the ML subset and add them to the CustomCVForce.

            bondedSystem = self._removeBonds(system, atoms, False, removeConstraints)
            bondedForces = []
            for force in bondedSystem.getForces():
                if hasattr(force, 'addBond') or hasattr(force, 'addAngle') or hasattr(force, 'addTorsion'):
                    bondedForces.append(force)
            mmVarNames = []
            for i, force in enumerate(bondedForces):
                name = f'mmForce{i+1}'
                cv.addCollectiveVariable(name, deepcopy(force))
                mmVarNames.append(name)

            # Create a CustomBondForce that computes all nonbonded interactions within the ML subset.

            for force in system.getForces():
                if isinstance(force, openmm.NonbondedForce):
                    internalNonbonded = openmm.CustomBondForce('138.935456*chargeProd/r + 4*epsilon*((sigma/r)^12-(sigma/r)^6)')
                    internalNonbonded.addPerBondParameter('chargeProd')
                    internalNonbonded.addPerBondParameter('sigma')
                    internalNonbonded.addPerBondParameter('epsilon')
                    numParticles = system.getNumParticles()
                    atomCharge = [0]*numParticles
                    atomSigma = [0]*numParticles
                    atomEpsilon = [0]*numParticles
                    for i in range(numParticles):
                        charge, sigma, epsilon = force.getParticleParameters(i)
                        atomCharge[i] = charge
                        atomSigma[i] = sigma
                        atomEpsilon[i] = epsilon
                    exceptions = {}
                    for i in range(force.getNumExceptions()):
                        p1, p2, chargeProd, sigma, epsilon = force.getExceptionParameters(i)
                        exceptions[(p1, p2)] = (chargeProd, sigma, epsilon)
                    for p1 in atomList:
                        for p2 in atomList:
                            if p1 == p2:
                                break
                            if (p1, p2) in exceptions:
                                chargeProd, sigma, epsilon = exceptions[(p1, p2)]
                            elif (p2, p1) in exceptions:
                                chargeProd, sigma, epsilon = exceptions[(p2, p1)]
                            else:
                                chargeProd = atomCharge[p1]*atomCharge[p2]
                                sigma = 0.5*(atomSigma[p1]+atomSigma[p2])
                                epsilon = unit.sqrt(atomEpsilon[p1]*atomEpsilon[p2])
                            if chargeProd._value != 0 or epsilon._value != 0:
                                internalNonbonded.addBond(p1, p2, [chargeProd, sigma, epsilon])
                    if internalNonbonded.getNumBonds() > 0:
                        name = f'mmForce{len(mmVarNames)+1}'
                        cv.addCollectiveVariable(name, internalNonbonded)
                        mmVarNames.append(name)

            # Configure the CustomCVForce so lambda_interpolate interpolates between the conventional and ML potentials.

            mlSum = '+'.join(mlVarNames) if len(mlVarNames) > 0 else '0'
            mmSum = '+'.join(mmVarNames) if len(mmVarNames) > 0 else '0'
            cv.setEnergyFunction(f'lambda_interpolate*({mlSum}) + (1-lambda_interpolate)*({mmSum})')
            newSystem.addForce(cv)
        return newSystem

    def _removeBonds(self, system: openmm.System, atoms: Iterable[int], removeInSet: bool, removeConstraints: bool) -> openmm.System:
        """Copy a System, removing all bonded interactions between atoms in (or not in) a particular set.

        Parameters
        ----------
        system: System
            the System to copy
        atoms: Iterable[int]
            a set of atom indices
        removeInSet: bool
            if True, any bonded term connecting atoms in the specified set is removed.  If False,
            any term that does *not* connect atoms in the specified set is removed
        removeConstraints: bool
            if True, remove constraints between pairs of atoms in the set

        Returns
        -------
        a newly created System object in which the specified bonded interactions have been removed
        """
        atomSet = set(atoms)

        # Create an XML representation of the System.

        import xml.etree.ElementTree as ET
        xml = openmm.XmlSerializer.serialize(system)
        root = ET.fromstring(xml)

        # This function decides whether a bonded interaction should be removed.

        def shouldRemove(termAtoms):
            return all(a in atomSet for a in termAtoms) == removeInSet

        # Remove bonds, angles, and torsions.

        for bonds in root.findall('./Forces/Force/Bonds'):
            for bond in bonds.findall('Bond'):
                bondAtoms = [int(bond.attrib[p]) for p in ('p1', 'p2')]
                if shouldRemove(bondAtoms):
                    bonds.remove(bond)
        for angles in root.findall('./Forces/Force/Angles'):
            for angle in angles.findall('Angle'):
                angleAtoms = [int(angle.attrib[p]) for p in ('p1', 'p2', 'p3')]
                if shouldRemove(angleAtoms):
                    angles.remove(angle)
        for torsions in root.findall('./Forces/Force/Torsions'):
            for torsion in torsions.findall('Torsion'):
                torsionLabels =  ('p1', 'p2', 'p3', 'p4') if 'p1' in torsion.attrib else ('a1', 'a2', 'a3', 'a4', 'b1', 'b2', 'b3', 'b4')
                torsionAtoms = [int(torsion.attrib[p]) for p in torsionLabels]
                if shouldRemove(torsionAtoms):
                    torsions.remove(torsion)

        # Optionally remove constraints.

        if removeConstraints:
            for constraints in root.findall('./Constraints'):
                for constraint in constraints.findall('Constraint'):
                    constraintAtoms = [int(constraint.attrib[p]) for p in ('p1', 'p2')]
                    if shouldRemove(constraintAtoms):
                        constraints.remove(constraint)

        # Create a new System from it.

        return openmm.XmlSerializer.deserialize(ET.tostring(root, encoding='unicode'))

    @staticmethod
    def registerImplFactory(name: str, factory: MLPotentialImplFactory):
        """Register a new potential function that can be used with MLPotential.

        Parameters
        ----------
        name: str
            the name of the potential function that will be passed to the MLPotential constructor
        factory: MLPotentialImplFactory
            a factory object that will be used to create MLPotentialImpl objects
        """
        MLPotential._implFactories[name] = factory


# Register any potential functions defined by entry points.

for potential in entry_points(group='openmmml.potentials'):
    MLPotential.registerImplFactory(potential.name, potential.load()())
