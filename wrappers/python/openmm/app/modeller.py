"""
modeller.py: Provides tools for editing molecular models

This is part of the OpenMM molecular simulation toolkit.
See https://openmm.org/development.

Portions copyright (c) 2012-2026 Stanford University and the Authors.
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
from __future__ import division
from __future__ import absolute_import

__author__ = "Peter Eastman"
__version__ = "1.0"

from openmm.app import Topology, PDBFile, ForceField
from openmm.app.forcefield import AllBonds, CutoffNonPeriodic, CutoffPeriodic, DrudeGenerator, _getDataDirectories
from openmm.app.internal import compiled
from openmm.vec3 import Vec3
from openmm import System, Context, NonbondedForce, AmoebaVdwForce, AmoebaMultipoleForce, CustomNonbondedForce, HarmonicBondForce, CustomBondForce, HarmonicAngleForce, \
    VerletIntegrator, LangevinIntegrator, CustomIntegrator, MonteCarloBarostat, LocalEnergyMinimizer, TwoParticleAverageSite, OutOfPlaneSite
from openmm.unit import nanometer, molar, elementary_charge, degree, acos, is_quantity, dot, norm, kilojoules_per_mole
import openmm.unit as unit
from . import element as elem
import gc
import os
import random
import sys
import xml.etree.ElementTree as etree
from copy import deepcopy
from math import ceil, floor, sqrt
from collections import defaultdict, namedtuple
import numpy as np

class Modeller(object):
    """Modeller provides tools for editing molecular models, such as adding water or missing hydrogens.

    To use it, create a Modeller object, specifying the initial Topology and atom positions.  You can
    then call various methods to change the model in different ways.  Each time you do, a new Topology
    and list of coordinates is created to represent the changed model.  Finally, call getTopology()
    and getPositions() to get the results.
    """

    _residueHydrogens = {}
    _hasLoadedStandardHydrogens = False

    # From Alvarez, Dalton Trans. 42, 8617 (2023).  All values are in nm.  Elements for which no data are available are
    # given a default radius of 0.2 nm.  When no element is specified for an atom, index 0 will be accessed; this is
    # given a radius of 0.05 nm so virtual sites don't take up much volume.
    _PACKING_RADII = np.array([
        0.05 , 0.12 , 0.143, 0.212, 0.198, 0.191, 0.177, 0.166, 0.15 , 0.146,
        0.158, 0.25 , 0.251, 0.225, 0.219, 0.19 , 0.189, 0.182, 0.183, 0.273,
        0.262, 0.258, 0.246, 0.242, 0.245, 0.245, 0.244, 0.24 , 0.24 , 0.238,
        0.239, 0.232, 0.229, 0.188, 0.182, 0.186, 0.225, 0.321, 0.284, 0.275,
        0.252, 0.256, 0.245, 0.244, 0.246, 0.244, 0.215, 0.253, 0.249, 0.243,
        0.242, 0.247, 0.199, 0.204, 0.206, 0.348, 0.303, 0.298, 0.288, 0.292,
        0.295, 0.2  , 0.29 , 0.287, 0.283, 0.279, 0.287, 0.281, 0.283, 0.279,
        0.28 , 0.274, 0.263, 0.253, 0.257, 0.249, 0.248, 0.241, 0.229, 0.232,
        0.245, 0.247, 0.26 , 0.254, 0.2  , 0.2  , 0.2  , 0.2  , 0.2  , 0.28 ,
        0.293, 0.288, 0.271, 0.282, 0.281, 0.283, 0.305, 0.34 , 0.305, 0.27 ,
        0.2  , 0.2  , 0.2  , 0.2  , 0.2  , 0.2  , 0.2  , 0.2  , 0.2  , 0.2  ,
        0.2  , 0.2  , 0.2  , 0.2  , 0.2  , 0.2  , 0.2  , 0.2  ,
    ])

    def __init__(self, topology, positions):
        """Create a new Modeller object

        Parameters
        ----------
        topology : Topology
            the initial Topology of the model
        positions : list
            the initial atomic positions
        """
        ## The Topology describing the structure of the system
        self.topology = topology
        if not is_quantity(positions):
            positions = positions*nanometer
        ## The list of atom positions
        self.positions = positions

    def getTopology(self):
        """Get the Topology of the model."""
        return self.topology

    def getPositions(self):
        """Get the atomic positions."""
        return self.positions

    def add(self, addTopology, addPositions):
        """Add chains, residues, atoms, and bonds to the model.

        Specify what to add by providing a new Topology object and the
        corresponding atomic positions. All chains, residues, atoms, and bonds
        contained in the Topology are added to the model.

        Parameters
        ----------
        addTopology : Topology
            a Topology whose contents should be added to the model
        addPositions : list
            the positions of the atoms to add
        """
        # Copy over the existing model.

        newTopology = Topology()
        newTopology.setPeriodicBoxVectors(self.topology.getPeriodicBoxVectors())
        newAtoms = {}
        newPositions = []*nanometer
        for chain in self.topology.chains():
            newChain = newTopology.addChain(chain.id)
            for residue in chain.residues():
                newResidue = newTopology.addResidue(residue.name, newChain, residue.id, residue.insertionCode)
                for atom in residue.atoms():
                    newAtom = newTopology.addAtom(atom.name, atom.element, newResidue, atom.id, atom.formalCharge)
                    newAtoms[atom] = newAtom
                    newPositions.append(deepcopy(self.positions[atom.index]))
        for bond in self.topology.bonds():
            newTopology.addBond(newAtoms[bond[0]], newAtoms[bond[1]], bond.type, bond.order)

        # Add the new model

        newAtoms = {}
        for chain in addTopology.chains():
            newChain = newTopology.addChain(chain.id)
            for residue in chain.residues():
                newResidue = newTopology.addResidue(residue.name, newChain, residue.id, residue.insertionCode)
                for atom in residue.atoms():
                    newAtom = newTopology.addAtom(atom.name, atom.element, newResidue, atom.id, atom.formalCharge)
                    newAtoms[atom] = newAtom
                    newPositions.append(deepcopy(addPositions[atom.index]))
        for bond in addTopology.bonds():
            newTopology.addBond(newAtoms[bond[0]], newAtoms[bond[1]], bond.type, bond.order)
        self.topology = newTopology
        self.positions = newPositions

    def delete(self, toDelete):
        """Delete chains, residues, atoms, and bonds from the model.

        You can specify objects to delete at any granularity: atoms, residues, or chains.  Passing
        in an Atom object causes that Atom to be deleted.  Passing in a Residue object causes that
        Residue and all Atoms it contains to be deleted.  Passing in a Chain object causes that
        Chain and all Residues and Atoms it contains to be deleted.

        In all cases, when an Atom is deleted, any bonds it participates in are also deleted.
        You also can specify a bond (as a tuple of Atom objects) to delete just that bond without
        deleting the Atoms it connects.

        Parameters
        ----------
        toDelete : list
            a list of Atoms, Residues, Chains, and bonds (specified as tuples of
            Atoms) to delete
        """
        newTopology = Topology()
        newTopology.setPeriodicBoxVectors(self.topology.getPeriodicBoxVectors())
        newAtoms = {}
        newPositions = []*nanometer
        deleteSet = set(toDelete)
        for chain in self.topology.chains():
            if chain not in deleteSet:
                needNewChain = True;
                for residue in chain.residues():
                    if residue not in deleteSet:
                        needNewResidue = True
                        for atom in residue.atoms():
                            if atom not in deleteSet:
                                if needNewChain:
                                    newChain = newTopology.addChain(chain.id)
                                    needNewChain = False;
                                if needNewResidue:
                                    newResidue = newTopology.addResidue(residue.name, newChain, residue.id, residue.insertionCode)
                                    needNewResidue = False;
                                newAtom = newTopology.addAtom(atom.name, atom.element, newResidue, atom.id, atom.formalCharge)
                                newAtoms[atom] = newAtom
                                newPositions.append(deepcopy(self.positions[atom.index]))
        for bond in self.topology.bonds():
            if bond[0] in newAtoms and bond[1] in newAtoms:
                if bond not in deleteSet and (bond[1], bond[0]) not in deleteSet:
                    newTopology.addBond(newAtoms[bond[0]], newAtoms[bond[1]], bond.type, bond.order)
        self.topology = newTopology
        self.positions = newPositions

    def deleteWater(self):
        """Delete all water molecules from the model."""
        self.delete(res for res in self.topology.residues() if res.name == "HOH")

    def convertWater(self, model='tip3p'):
        """Convert all water molecules to a different water model.

        @deprecated Use addExtraParticles() instead.  It performs the same function but in a more general way.

        Parameters
        ----------
        model : string='tip3p'
            the water model to convert to.  Supported values are 'tip3p',
            'spce', 'tip4pew', and 'tip5p'.
        """
        if model in ('tip3p', 'spce'):
            sites = 3
        elif model == 'tip4pew':
            sites = 4
        elif model == 'tip5p':
            sites = 5
        else:
            raise ValueError('Unknown water model: %s' % model)
        newTopology = Topology()
        newTopology.setPeriodicBoxVectors(self.topology.getPeriodicBoxVectors())
        newAtoms = {}
        newPositions = []*nanometer
        for chain in self.topology.chains():
            newChain = newTopology.addChain(chain.id)
            for residue in chain.residues():
                newResidue = newTopology.addResidue(residue.name, newChain, residue.id, residue.insertionCode)
                if residue.name == "HOH":
                    # Copy the oxygen and hydrogens
                    oatom = [atom for atom in residue.atoms() if atom.element == elem.oxygen]
                    hatoms = [atom for atom in residue.atoms() if atom.element == elem.hydrogen]
                    if len(oatom) != 1 or len(hatoms) != 2:
                        raise ValueError('Illegal water molecule (residue %d): contains %d oxygen(s) and %d hydrogen(s)' % (residue.index, len(oatom), len(hatoms)))
                    o = newTopology.addAtom(oatom[0].name, oatom[0].element, newResidue)
                    h1 = newTopology.addAtom(hatoms[0].name, hatoms[0].element, newResidue)
                    h2 = newTopology.addAtom(hatoms[1].name, hatoms[1].element, newResidue)
                    newAtoms[oatom[0]] = o
                    newAtoms[hatoms[0]] = h1
                    newAtoms[hatoms[1]] = h2
                    po = deepcopy(self.positions[oatom[0].index])
                    ph1 = deepcopy(self.positions[hatoms[0].index])
                    ph2 = deepcopy(self.positions[hatoms[1].index])
                    newPositions.append(po)
                    newPositions.append(ph1)
                    newPositions.append(ph2)

                    # Add virtual sites.

                    if sites == 4:
                        newTopology.addAtom('M', None, newResidue)
                        newPositions.append(0.786646558*po + 0.106676721*ph1 + 0.106676721*ph2)
                    elif sites == 5:
                        newTopology.addAtom('M1', None, newResidue)
                        newTopology.addAtom('M2', None, newResidue)
                        v1 = (ph1-po).value_in_unit(nanometer)
                        v2 = (ph2-po).value_in_unit(nanometer)
                        cross = Vec3(v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0])
                        newPositions.append(po - (0.34490826*v1 - 0.34490826*v2 - 6.4437903*cross)*nanometer)
                        newPositions.append(po - (0.34490826*v1 - 0.34490826*v2 + 6.4437903*cross)*nanometer)
                else:
                    # Just copy the residue over.
                    for atom in residue.atoms():
                        newAtom = newTopology.addAtom(atom.name, atom.element, newResidue, atom.id, atom.formalCharge)
                        newAtoms[atom] = newAtom
                        newPositions.append(deepcopy(self.positions[atom.index]))
        for bond in self.topology.bonds():
            if bond[0] in newAtoms and bond[1] in newAtoms:
                newTopology.addBond(newAtoms[bond[0]], newAtoms[bond[1]], bond.type, bond.order)
        self.topology = newTopology
        self.positions = newPositions

    def _addIons(self, forcefield, numWaters, replaceableMols, ionCutoff=0.5*nanometer, positiveIon='Na+', negativeIon='Cl-', ionicStrength=0*molar, neutralize=True, residueTemplates=dict()):
        """Adds ions to the system by replacing certain molecules.

        Parameters
        ----------
        forcefield : ForceField
            the ForceField to use to determine the total charge of the system.
        numWaters : int
            the total number of water molecules in the simulation box, used to
            calculate the number of ions / concentration to add.
        replaceableMols : dict
            the molecules to replace by ions, as a dictionary of residue:positions
        ionCutoff: distance=0.5*nanometer
        positiveIon : string='Na+'
            the type of positive ion to add.  Allowed values are 'Cs+', 'K+', 'Li+', 'Na+', and 'Rb+'
        negativeIon : string='Cl-'
            the type of negative ion to add.  Allowed values are 'Cl-', 'Br-', 'F-', and 'I-'. Be aware
            that not all force fields support all ion types.
        ionicStrength : concentration=0*molar
            the total concentration of ions (both positive and negative) to add.  This
            does not include ions that are added to neutralize the system.
            Note that only monovalent ions are currently supported.
        neutralize : bool=True
            whether to add ions to neutralize the system
        residueTemplates : dict=dict()
            specifies which template the ForceField should use for particular residues.  The keys
            should be Residue objects from the Topology, and the values should be the names of the
            templates to use for them.  This is useful when a ForceField contains multiple templates
            that can match the same residue (e.g Fe2+ and Fe3+ templates in the ForceField for a
            monoatomic iron ion in the Topology).
        """

        if not neutralize and ionicStrength == 0 *molar:
            return

        posIonElements = {'Cs+': elem.cesium, 'K+': elem.potassium,
                          'Li+': elem.lithium, 'Na+': elem.sodium,
                          'Rb+': elem.rubidium}
        negIonElements = {'Cl-': elem.chlorine, 'Br-': elem.bromine,
                          'F-': elem.fluorine, 'I-': elem.iodine}

        ionPositions = []

        numReplaceableMols = len(replaceableMols)

        # Fetch ion elements from user input
        if positiveIon not in posIonElements:
            raise ValueError('Illegal value for positive ion: {}'.format(positiveIon))
        if negativeIon not in negIonElements:
            raise ValueError('Illegal value for negative ion: {}'.format(negativeIon))
        positiveElement = posIonElements[positiveIon]
        negativeElement = negIonElements[negativeIon]

        # Determine the total charge of the system
        system = forcefield.createSystem(self.topology, residueTemplates=residueTemplates)
        for i in range(system.getNumForces()):
            if isinstance(system.getForce(i), (NonbondedForce, AmoebaMultipoleForce)):
                nonbonded = system.getForce(i)
                break
        else:
            raise ValueError('The ForceField does not specify a NonbondedForce')

        totalCharge = 0.0
        if isinstance(nonbonded, AmoebaMultipoleForce):
            for i in range(nonbonded.getNumMultipoles()):
                totalCharge += nonbonded.getMultipoleParameters(i)[0].value_in_unit(elementary_charge)
        else:
            for i in range(nonbonded.getNumParticles()):
                totalCharge += nonbonded.getParticleParameters(i)[0].value_in_unit(elementary_charge)

        # Round to nearest integer
        totalCharge = int(floor(0.5 + totalCharge))

        # Figure out how many ions to add based on requested params/concentration
        numPositive, numNegative = 0, 0
        if neutralize:
            if abs(totalCharge) > numReplaceableMols:
                raise Exception('Cannot neutralize the system because the charge is greater than the number of available positions for ions')
            if totalCharge > 0:
                numNegative += totalCharge
            else:
                numPositive -= totalCharge

        if ionicStrength > 0 * molar:
            numIons = (numWaters - numPositive - numNegative) * ionicStrength / (55.4 * molar)  # Pure water is about 55.4 molar (depending on temperature)
            numPairs = int(floor(numIons + 0.5))
            numPositive += numPairs
            numNegative += numPairs
        totalIons = numPositive + numNegative

        if totalIons > 0:
            # Randomly select a set of waters
            # while ensuring ions are not placed too close to each other.
            modeller = Modeller(self.topology, self.positions)

            replaceableList = list(replaceableMols.keys())
            numAddedIons = 0
            numTrials = 10  # Attempts to add ions N times before quitting
            toReplace = []  # list of molecules to be replaced
            while numAddedIons < totalIons:
                pickedMol = random.choice(replaceableList)
                replaceableList.remove(pickedMol)
                # Check distance to other ions
                for pos in ionPositions:
                    distance = norm(pos - replaceableMols[pickedMol])
                    if distance <= ionCutoff:
                        numTrials -= 1
                        break
                else:
                    toReplace.append(pickedMol)
                    ionPositions.append(replaceableMols[pickedMol])
                    numAddedIons += 1

                    n_trials = 10

                if n_trials == 0:
                    raise ValueError('Could not add more than {} ions to the system'.format(numAddedIons))

            # Replace waters/ions in the topology
            modeller.delete(toReplace)
            ionChain = modeller.topology.addChain()
            for i, water in enumerate(toReplace):
                element = (positiveElement if i < numPositive else negativeElement)
                newResidue = modeller.topology.addResidue(element.symbol.upper(), ionChain)
                modeller.topology.addAtom(element.symbol, element, newResidue)
                modeller.positions.append(replaceableMols[water])

            # Update topology/positions
            self.topology = modeller.topology
            self.positions = modeller.positions

    def addSolvent(self, forcefield, model='tip3p', boxSize=None, boxVectors=None, padding=None, numAdded=None, boxShape='cube', positiveIon='Na+', negativeIon='Cl-', ionicStrength=0*molar, neutralize=True, residueTemplates=dict()):
        """Add solvent (both water and ions) to the model to fill a periodic box.

        The algorithm works as follows:

        1. Water molecules are added to fill the box.
        2. Water molecules are removed if their distance to any solute atom is less than the sum of their van der Waals radii.
        3. If the solute is charged and neutralize=True, enough positive or negative ions are added to neutralize it.  Each ion is added by
           randomly selecting a water molecule and replacing it with the ion.
        4. Ion pairs are added to give the requested total ionic strength.  Note that only monovalent ions are currently supported.

        The box size can be specified in any of several ways:

        1. You can explicitly give the vectors defining the periodic box to use.
        2. Alternatively, for a rectangular box you can simply give the dimensions of the unit cell.
        3. You can give a padding distance.  A bounding sphere containing the solute is determined, and the box size is
           set to (sphere diameter)+(padding).  This guarantees no atom in the solute will come closer than the padding
           distance to any atom of another periodic copy.  If the sphere diameter is less than the padding distance,
           the box size is set to 2*(padding) to ensure no atom is closer than the padding distance to two periodic
           copies of any other atom.
        4. You can specify the total number of molecules (both waters and ions) to add.  A box is then created whose size is
           just large enough hold the specified amount of solvent.
        5. Finally, if none of the above options is specified, the existing Topology's box vectors are used.

        When specifying either a padding distance or a number of molecules, you can specify a shape for the periodic box:
        cubic, rhombic dodecahedron, or truncated octahedron.  Using a non-rectangular box allows the same distance
        between periodic copies to be achieved with a smaller box.  The most compact option is a rhombic dodecahedron,
        for which the box volume is 70.7% the volume of a cubic box with the same amount of padding.

        There exist many different water models, many of which are very similar to each other. This method creates
        preequilibrated water boxes for a limited set of water models. In most cases, they work equally well for other models
        that involve the same number of particles. For example, to simulate a box of TIP3P-FB water, use this method to
        create a box of TIP3P water, construct a System using TIP3P-FB parameters, and perform a local energy minimization
        to correct for the small differences between the models. Likewise, a box of TIP4P-Ew water can be used for most
        four site water models.

        Parameters
        ----------
        forcefield : ForceField
            the ForceField to use for determining van der Waals radii and atomic charges
        model : str='tip3p'
            the water model to use.  Supported values are 'tip3p', 'spce', 'tip4pew', 'tip5p', and 'swm4ndp' (polarizable).
        boxSize : Vec3=None
            the size of the box to fill with water
        boxVectors : tuple of Vec3=None
            the vectors defining the periodic box to fill with water
        padding : distance=None
            the padding distance to use
        numAdded : int=None
            the total number of molecules (waters and ions) to add
        boxShape: str='cube'
            the box shape to use.  Allowed values are 'cube', 'dodecahedron', and 'octahedron'.  If padding and numAdded
            are both None, this is ignored.
        positiveIon : string='Na+'
            the type of positive ion to add.  Allowed values are 'Cs+', 'K+', 'Li+', 'Na+', and 'Rb+'
        negativeIon : string='Cl-'
            the type of negative ion to add.  Allowed values are 'Cl-', 'Br-', 'F-', and 'I-'. Be aware
            that not all force fields support all ion types.
        ionicStrength : concentration=0*molar
            the total concentration of ions (both positive and negative) to add.  This
            does not include ions that are added to neutralize the system.
            Note that only monovalent ions are currently supported.
        neutralize : bool=True
            whether to add ions to neutralize the system
        residueTemplates : dict=dict()
            specifies which template the ForceField should use for particular residues.  The keys
            should be Residue objects from the Topology, and the values should be the names of the
            templates to use for them.  This is useful when a ForceField contains multiple templates
            that can match the same residue (e.g Fe2+ and Fe3+ templates in the ForceField for a
            monoatomic iron ion in the Topology).
        """
        if len([x for x in (boxSize, boxVectors, padding, numAdded) if x is not None]) > 1:
            raise ValueError('At most one of the following arguments may be specified: boxSize, boxVectors, padding, numAdded')

        # Load the pre-equilibrated water box.

        vdwRadiusPerSigma = 0.5612310241546864907
        if model == 'tip3p':
            waterRadius = 0.31507524065751241*vdwRadiusPerSigma
        elif model == 'spce':
            waterRadius = 0.31657195050398818*vdwRadiusPerSigma
        elif model == 'tip4pew':
            waterRadius = 0.315365*vdwRadiusPerSigma
        elif model == 'tip5p':
            waterRadius = 0.312*vdwRadiusPerSigma
        elif model == 'swm4ndp':
            waterRadius = 0.318395*vdwRadiusPerSigma
        else:
            raise ValueError('Unknown water model: %s' % model)
        pdb = PDBFile(os.path.join(os.path.dirname(__file__), 'data', model+'.pdb'))
        pdbTopology = pdb.getTopology()
        pdbPositions = pdb.getPositions().value_in_unit(nanometer)
        pdbResidues = list(pdbTopology.residues())
        pdbBoxSize = pdbTopology.getUnitCellDimensions().value_in_unit(nanometer)

        # Pick a unit cell size.

        if numAdded is not None:
            # Select a padding distance which is guaranteed to give more than the specified number of molecules.

            padding = 2.2*(numAdded/((len(pdbResidues)/pdbBoxSize[0]**3)*8))**(1.0/3.0)
            if padding < 0.5:
                padding = 0.5 # Ensure we have enough when adding very small numbers of molecules
        if boxVectors is not None:
            if is_quantity(boxVectors[0]):
                boxVectors = (boxVectors[0].value_in_unit(nanometer), boxVectors[1].value_in_unit(nanometer), boxVectors[2].value_in_unit(nanometer))
            box = Vec3(boxVectors[0][0], boxVectors[1][1], boxVectors[2][2])
            vectors = boxVectors
        elif boxSize is not None:
            if is_quantity(boxSize):
                boxSize = boxSize.value_in_unit(nanometer)
            box = Vec3(boxSize[0], boxSize[1], boxSize[2])
            vectors = (Vec3(boxSize[0], 0, 0), Vec3(0, boxSize[1], 0), Vec3(0, 0, boxSize[2]))
        elif padding is not None:
            width = self._computeBoxWidth(padding)
            vectors = self._computeBoxVectors(width, boxShape)
            box = Vec3(vectors[0][0], vectors[1][1], vectors[2][2])
        else:
            box = self.topology.getUnitCellDimensions().value_in_unit(nanometer)
            vectors = self.topology.getPeriodicBoxVectors().value_in_unit(nanometer)
            if box is None:
                raise ValueError('Neither the box size, box vectors, nor padding was specified, and the Topology does not define unit cell dimensions')

        # Have the ForceField build a System for the solute from which we can determine van der Waals radii.

        system = forcefield.createSystem(self.topology, residueTemplates=residueTemplates)
        nonbonded = None
        for i in range(system.getNumForces()):
            if isinstance(system.getForce(i), (NonbondedForce, AmoebaVdwForce)):
                nonbonded = system.getForce(i)
      
        if nonbonded is None:
            raise ValueError('The ForceField does not specify a NonbondedForce')
        cutoff = [waterRadius]*system.getNumParticles()

        for i in range(system.getNumParticles()):
            params = nonbonded.getParticleParameters(i)
            if params[2] != 0*kilojoules_per_mole:
                cutoff[i] += params[1].value_in_unit(nanometer)*vdwRadiusPerSigma

        waterCutoff = waterRadius
        if len(cutoff) == 0:
            maxCutoff = waterCutoff
        else:
            maxCutoff = max(waterCutoff, max(cutoff))

        # Copy the solute over.

        newTopology = Topology()
        newTopology.setPeriodicBoxVectors(vectors*nanometer)
        newAtoms = {}
        newPositions = []*nanometer
        newResidueTemplates=dict()
        for chain in self.topology.chains():
            newChain = newTopology.addChain(chain.id)
            for residue in chain.residues():
                newResidue = newTopology.addResidue(residue.name, newChain, residue.id, residue.insertionCode)
                if residue in residueTemplates:
                    newResidueTemplates[newResidue] = residueTemplates[residue]
                for atom in residue.atoms():
                    newAtom = newTopology.addAtom(atom.name, atom.element, newResidue, atom.id, atom.formalCharge)
                    newAtoms[atom] = newAtom
                    newPositions.append(deepcopy(self.positions[atom.index]))
        for bond in self.topology.bonds():
            newTopology.addBond(newAtoms[bond[0]], newAtoms[bond[1]], bond.type, bond.order)

        # Sort the solute atoms into cells for fast lookup.

        if len(self.positions) == 0:
            positions = []
        else:
            positions = deepcopy(self.positions.value_in_unit(nanometer))
        cells = compiled.CellList(positions, maxCutoff, vectors, True)

        # Create a function to compute the distance between two points, taking periodic boundary conditions into account.

        periodicDistance = compiled.periodicDistance(vectors)

        # Find the list of water molecules to add.

        newChain = newTopology.addChain()
        if len(positions) == 0:
            center = Vec3(0, 0, 0)
        else:
            center = [(max((pos[i] for pos in positions))+min((pos[i] for pos in positions)))/2 for i in range(3)]
            center = Vec3(center[0], center[1], center[2])
        numBoxes = [int(ceil(box[i]/pdbBoxSize[i])) for i in range(3)]
        addedWaters = []
        for residue in pdbResidues:
            oxygen = [atom for atom in residue.atoms() if atom.element == elem.oxygen][0]
            for boxx in range(numBoxes[0]):
                for boxy in range(numBoxes[1]):
                    for boxz in range(numBoxes[2]):
                        offset = Vec3(boxx*pdbBoxSize[0], boxy*pdbBoxSize[1], boxz*pdbBoxSize[2])
                        atomPos = pdbPositions[oxygen.index]+offset
                        if not any([atomPos[i] > box[i] for i in range(3)]):
                            # This molecule is inside the box, so see how close to it is to the solute.

                            atomPos += center-box/2
                            for i in cells.neighbors(atomPos):
                                if periodicDistance(atomPos, positions[i]) < cutoff[i]:
                                    break
                            else:
                                # Record this water molecule as one to add.

                                addedWaters.append((residue.index, atomPos))

        if numAdded is not None:
            # We added many more waters than we actually want.  Sort them based on distance to the nearest box edge and
            # only keep the ones in the middle.

            lowerBound = center-box/2
            upperBound = center+box/2
            distToEdge = (min(min(pos-lowerBound), min(upperBound-pos)) for index, pos in addedWaters)
            sortedIndex = [i[0] for i in sorted(enumerate(distToEdge), key=lambda x: -x[1])]
            addedWaters = [addedWaters[i] for i in sortedIndex[:numAdded]]

            # Compute a new periodic box size.

            maxSize = max(max((pos[i] for index, pos in addedWaters))-min((pos[i] for index, pos in addedWaters)) for i in range(3))
            maxSize += 0.1  # Add padding to reduce clashes at the edge.
            newTopology.setPeriodicBoxVectors(self._computeBoxVectors(maxSize, boxShape))
        else:
            # There could be clashes between water molecules at the box edges.  Find ones to remove.

            upperCutoff = center+box/2-Vec3(waterCutoff, waterCutoff, waterCutoff)
            lowerCutoff = center-box/2+Vec3(waterCutoff, waterCutoff, waterCutoff)
            lowerSkinPositions = [pos for index, pos in addedWaters if pos[0] < lowerCutoff[0] or pos[1] < lowerCutoff[1] or pos[2] < lowerCutoff[2]]
            filteredWaters = []
            cells.cells = {}
            for i in range(len(lowerSkinPositions)):
                cell = cells.cellForPosition(*lowerSkinPositions[i])
                if cell in cells.cells:
                    cells.cells[cell].append(i)
                else:
                    cells.cells[cell] = [i]
            for entry in addedWaters:
                pos = entry[1]
                if pos[0] < upperCutoff[0] and pos[1] < upperCutoff[1] and pos[2] < upperCutoff[2]:
                    filteredWaters.append(entry)
                else:
                    if not any((periodicDistance(lowerSkinPositions[i], pos) < waterCutoff and norm(lowerSkinPositions[i]-pos) > waterCutoff for i in cells.neighbors(pos))):
                        filteredWaters.append(entry)
            addedWaters = filteredWaters

        # Add the water molecules.
        waterPos = {}
        for index, pos in addedWaters:
            newResidue = newTopology.addResidue(residue.name, newChain)
            residue = pdbResidues[index]
            oxygen = [atom for atom in residue.atoms() if atom.element == elem.oxygen][0]
            oPos = pdbPositions[oxygen.index]
            molAtoms = []
            for atom in residue.atoms():
                molAtoms.append(newTopology.addAtom(atom.name, atom.element, newResidue))
                newPositions.append((pos+pdbPositions[atom.index]-oPos)*nanometer)
                if atom.element == elem.oxygen:
                    waterPos[newResidue] = newPositions[-1]
            for atom1 in molAtoms:
                if atom1.element == elem.oxygen:
                    for atom2 in molAtoms:
                        if atom2.element == elem.hydrogen:
                            newTopology.addBond(atom1, atom2)

        self.topology = newTopology
        self.positions = newPositions

        # Total number of waters in the box
        numTotalWaters = len(waterPos)

        # Add ions to neutralize the system.
        self._addIons(forcefield, numTotalWaters, waterPos, positiveIon=positiveIon, negativeIon=negativeIon, ionicStrength=ionicStrength, neutralize=neutralize, residueTemplates=newResidueTemplates)

    class _ComponentData:
        """Inner class used to encapsulate data about a mixture component."""
        def __init__(self, topology, positions, count, radii):
            self.topology = topology
            self.positions = positions
            self.count = count
            self.radii = radii

    def addMixture(self, components, boxSize=None, boxVectors=None, padding=None, boxShape='cube', platform=None):
        """Adds copies of one or more components to the model to fill a periodic box.  The requested number of copies of
        each component will be added to the box.  Next, positions and orientations of each copy (treated as a rigid
        body) will be optimized to reduce overlaps between all of the existing and newly added atoms.  The method used
        is similar to that described in Martínez and Martínez, J. Comput. Chem. 24, 819 (2003), and Martínez et al.,
        J. Comput. Chem. 30, 2157 (2009), along with some heuristics for choosing a box size if none is specified.

        The box size can be specified in any of several ways:

        1. You can explicitly give the vectors defining the periodic box to use.
        2. Alternatively, for a rectangular box you can simply give the dimensions of the unit cell.
        3. You can give a padding distance.  A bounding sphere containing the solute is determined, and the box size is
           set to (sphere diameter)+(padding).  This guarantees no atom in the solute will come closer than the padding
           distance to any atom of another periodic copy.  If the sphere diameter is less than the padding distance,
           the box size is set to 2*(padding) to ensure no atom is closer than the padding distance to two periodic
           copies of any other atom.
        4. If none of the above options are specified, and the starting Topology has box vectors, they will be used.
        5. Otherwise, Modeller will try to guess a reasonable initial density and calculate box vectors based on it.

        When specifying a padding distance, you can specify a shape for the periodic box:
        cubic, rhombic dodecahedron, or truncated octahedron.  Using a non-rectangular box allows the same distance
        between periodic copies to be achieved with a smaller box.  The most compact option is a rhombic dodecahedron,
        for which the box volume is 70.7% the volume of a cubic box with the same amount of padding.

        This method can also be used to pack boxes of mixtures and pure components without a solute.  Simply start from
        a Modeller with an empty Topology and no positions.

        Parameters
        ----------
        components : list
            the components to add.  Each component must be provided as a dictionary, with the required keys 'topology'
            holding a Topology for the component, 'positions' holding positions to use for instances of the component,
            and 'count' specifying the number of copies of the component to add.
        boxSize : Vec3=None
            the size of the box to fill
        boxVectors : tuple of Vec3=None
            the vectors defining the periodic box to fill
        padding : distance=None
            the padding distance to use
        boxShape : str='cube'
            the box shape to use.  Allowed values are 'cube', 'dodecahedron', and 'octahedron'.  If boxSize or
            boxVectors is specified, this is ignored.
        platform : Platform=None
            the Platform to use for solving the packing problem.  If this is None, the default Platform will be used.
        """

        # Helper functions for addMixture().

        def addCopyToTopology(dstTopology, srcTopology):
            """Add a copy of the given source topology to a given destination topology."""
            dstAtoms = []
            for srcChain in srcTopology.chains():
                dstChain = dstTopology.addChain()
                for srcResidue in srcChain.residues():
                    dstResidue = dstTopology.addResidue(srcResidue.name, dstChain)
                    for srcAtom in srcResidue.atoms():
                        dstAtoms.append(dstTopology.addAtom(srcAtom.name, srcAtom.element, dstResidue, srcAtom.id, srcAtom.formalCharge))
            for srcBond in srcTopology.bonds():
                dstTopology.addBond(dstAtoms[srcBond.atom1.index], dstAtoms[srcBond.atom2.index], srcBond.type, srcBond.order)

        def addCopyToChain(dstChain, srcTopology):
            """Add a copy of all residues in the given source topology to a given destination chain."""
            dstTopology = dstChain.topology
            dstAtoms = []
            for srcResidue in srcTopology.residues():
                dstResidue = dstTopology.addResidue(srcResidue.name, dstChain)
                for srcAtom in srcResidue.atoms():
                    dstAtoms.append(dstTopology.addAtom(srcAtom.name, srcAtom.element, dstResidue, srcAtom.id, srcAtom.formalCharge))
            for srcBond in srcTopology.bonds():
                dstTopology.addBond(dstAtoms[srcBond.atom1.index], dstAtoms[srcBond.atom2.index], srcBond.type, srcBond.order)

        def processPositions(topology, positions):
            """Finds suitable reference atoms in a set of positions to constrain."""
            numAtoms = topology.getNumAtoms()

            # Center the component.

            positions = positions - np.mean(positions, axis=0)

            if numAtoms < 3:
                # There are no choices for the reference atoms.

                return positions, np.arange(numAtoms)

            # Choose the first reference atom as the one furthest from the center.

            index1 = np.argmax(np.sum(positions * positions, axis=1))

            # Choose the second reference atom as the one furthest from the first.

            delta1x = positions - positions[index1]
            objective2 = np.linalg.norm(delta1x, axis=1)
            index2 = np.argmax(objective2)

            # Choose the third reference atom as the one furthest from the line passing through the first two.

            dir12 = positions[index2] - positions[index1]
            dir12 /= np.linalg.norm(dir12)
            objective3 = np.linalg.norm(delta1x - dir12 * np.sum(dir12 * delta1x, axis=1, keepdims=True), axis=1)
            objective3[[index1, index2]] = np.nan
            index3 = np.nanargmax(objective3)

            # If the molecule is linear or very close to linear, don't return the third reference atom since it will create
            # a degenerate triangle.  This tolerance is fairly loose to keep the constraint solver from failing, so slightly
            # non-linear molecules will be inaccurately treated as linear, but this should be good enough for packing.

            if objective3[index3] < 0.01 * objective2[index2]:
                indices = np.array([index1, index2])
            else:
                indices = np.array([index1, index2, index3])
            return positions, indices

        def placeRandomly(positions, vectors, baseRadii, radii, cells, periodicDistance):
            """Place positions representing a copy of a component randomly in a periodic box."""

            # Search for a placement such that none of the atoms of the new copy overlap with the Modeller's existing atoms.
            # We also do not attempt to avoid overlaps between new copies since they will move during the optimization.

            for _ in range(100):
                # Try to find a reasonable position for the center of the new copy (away from existing atoms).

                centerPos = np.random.uniform(-0.5, 0.5, 3) @ vectors
                if cells is not None and cells.cellForPosition(*centerPos) in cells.cells:
                    continue

                # Generate a random rotation matrix: see scipy.stats.ortho_group.rv() and scipy.stats.special_ortho_group.rv().

                q, r = np.linalg.qr(np.random.normal(size=(3, 3)))
                d = np.diag(r)
                q *= d / np.abs(d)
                q[0] /= np.linalg.det(q)

                # Rotate the copy and place it in the box, then check its distances to existing atoms.

                newPositions = centerPos + positions @ q
                if cells is not None and any(
                    periodicDistance(newPosition, cells.positions[neighbor]) < baseRadii[neighbor] + radius
                    for newPosition, radius in zip(newPositions, radii)
                    for neighbor in cells.neighbors(newPosition)
                ):
                    continue

                return newPositions

            raise ValueError("No space found for component; use a larger box size or padding distance")

        def getPackingRadii(topology):
            """Return van der Waals radii for each atom in a Topology."""
            return Modeller._PACKING_RADII[[0 if atom.element is None else atom.element.atomic_number for atom in topology.atoms()]]

        # Check input arguments.

        if len([x for x in (boxSize, boxVectors, padding) if x is not None]) > 1:
            raise ValueError('At most one of the following arguments may be specified: boxSize, boxVectors, padding')

        # Process components.

        componentList = []
        for component in components:
            component = dict(component)
            componentTopology = component.pop('topology')
            componentPositions = component.pop('positions')
            componentCount = component.pop('count')
            if component:
                raise ValueError(f'Unsupported keys in component dictionary: {list(component.keys())}')

            if componentTopology.getNumAtoms() != len(componentPositions):
                raise ValueError('The number of atoms in the Topology of a component must match the number of position coordinates')
            if componentTopology.getNumAtoms() == 0 or componentCount == 0:
                continue
            if is_quantity(componentPositions[0]):
                componentPositions = [position.value_in_unit(nanometer) for position in componentPositions]
            componentList.append(Modeller._ComponentData(componentTopology, np.array(componentPositions), componentCount, getPackingRadii(componentTopology)))

        if self.topology.getNumAtoms():
            baseRadii = getPackingRadii(self.topology)
        else:
            baseRadii = None

        # Determine what box vectors to use for packing.

        vectorsGuessed = False
        if boxVectors is not None:
            if is_quantity(boxVectors[0]):
                boxVectors = (boxVectors[0].value_in_unit(nanometer), boxVectors[1].value_in_unit(nanometer), boxVectors[2].value_in_unit(nanometer))
            vectors = np.array(boxVectors)
        elif boxSize is not None:
            if is_quantity(boxSize):
                boxSize = boxSize.value_in_unit(nanometer)
            vectors = np.diag(boxSize)
        elif padding is not None:
            width = self._computeBoxWidth(padding)
            vectors = np.array(self._computeBoxVectors(width, boxShape))
        elif self.topology.getPeriodicBoxVectors() is not None:
            vectors = np.array(self.topology.getPeriodicBoxVectors().value_in_unit(nanometer))
        else:
            # Make a rough estimate of the initial density to be refined later.

            guessVol = sum(component.count * np.sum(component.radii ** 3) for component in componentList)
            if baseRadii is not None:
                guessVol += np.sum(baseRadii ** 3)
            vectors = np.array(self._computeBoxVectors(1.0, boxShape))

            # This constant is chosen to make the density of pure water come to 1 g/cm^3 exactly.

            vectors *= np.cbrt(4.37932210654291 * guessVol / np.linalg.det(vectors))
            vectorsGuessed = True

        # Set up a Topology, positions, and System for solving the packing problem.

        topology = Topology()
        positions = []
        system = System()

        force = CustomNonbondedForce("p^2/(1-lambda*p)^6; p=max(1-r/(radius1+radius2),0)")
        force.addGlobalParameter("lambda", 0)
        force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
        force.addPerParticleParameter("radius")
        system.addForce(force)

        topology.setPeriodicBoxVectors(vectors)
        system.setDefaultPeriodicBoxVectors(*vectors)

        # Add the existing model.  Prepare a cell list to make sure new components are not added inside it.

        if self.topology.getNumAtoms():
            addCopyToTopology(topology, self.topology)

            basePositions = np.array(self.positions.value_in_unit(nanometer))
            positions.append(basePositions)

            for radius in baseRadii:
                system.addParticle(0)
                force.addParticle([radius])

            cells = compiled.CellList(basePositions, np.amax(baseRadii), vectors, True)
            periodicDistance = compiled.periodicDistance(vectors)
        else:
            cells = None
            periodicDistance = None

        # Add new components.

        for component in componentList:
            dstChain = topology.addChain()
            for _ in range(component.count):
                addCopyToChain(dstChain, component.topology)

            # Get the component positions, and indices of the particles that should be real (not virtual sites).

            refPositions, refIndices = processPositions(component.topology, component.positions)

            # Place copies of the component randomly in the box.

            for _ in range(component.count):
                positions.append(placeRandomly(refPositions, vectors, baseRadii, component.radii, cells, periodicDistance))

            if refIndices.size == 1:
                for _ in range(component.count):
                    system.addParticle(1)

                    force.addParticle([component.radii[0]])

            elif refIndices.size == 2:
                # The component is linear; constrain two atoms and use virtual sites for all other atoms.

                refPos1, refPos2 = refPositions[refIndices]
                refVec = refPos2 - refPos1
                refDist = np.linalg.norm(refVec)

                refWeights = (refPositions - refPos1) @ refVec / (refDist * refDist)

                for _ in range(component.count):
                    offset = system.getNumParticles()
                    index1, index2 = offset + refIndices

                    for index, weight in enumerate(refWeights):
                        if index in refIndices:
                            system.addParticle(1)
                        else:
                            system.setVirtualSite(system.addParticle(0), TwoParticleAverageSite(index1, index2, 1 - weight, weight))

                        force.addParticle([component.radii[index]])

                    system.addConstraint(index1, index2, refDist)

            else:
                # The component has three or more atoms that are not all collinear; constrain a triangle and use virtual
                # sites for all other atoms.

                refPos1, refPos2, refPos3 = refPositions[refIndices]
                refVec1 = refPos2 - refPos1
                refVec2 = refPos3 - refPos1

                refWeights = np.linalg.solve(np.stack([refVec1, refVec2, np.cross(refVec1, refVec2)]).T, (refPositions - refPos1).T).T

                dist12 = np.linalg.norm(refVec1)
                dist13 = np.linalg.norm(refVec2)
                dist23 = np.linalg.norm(refPos3 - refPos2)

                # OpenMM can use SETTLE instead of CCMA for constraint triangles if two distances are exactly equal,
                # which speeds up context setup.  Try to make a pair of distances exactly equal if they are close.

                if abs(dist12 - dist13) <= 0.001:
                    dist12 = dist13 = (dist12 + dist13) / 2
                elif abs(dist12 - dist23) <= 0.001:
                    dist12 = dist23 = (dist12 + dist23) / 2
                elif abs(dist13 - dist23) <= 0.001:
                    dist13 = dist23 = (dist13 + dist23) / 2

                for _ in range(component.count):
                    offset = system.getNumParticles()
                    index1, index2, index3 = offset + refIndices

                    for index, (weight12, weight13, weightCross) in enumerate(refWeights):
                        if index in refIndices:
                            system.addParticle(1)
                        else:
                            system.setVirtualSite(system.addParticle(0), OutOfPlaneSite(index1, index2, index3, weight12, weight13, weightCross))

                        force.addParticle([component.radii[index]])

                    system.addConstraint(index1, index2, dist12)
                    system.addConstraint(index1, index3, dist13)
                    system.addConstraint(index2, index3, dist23)

        if not system.getNumParticles():
            return

        if vectorsGuessed:
            # If the initial box vectors were guesses, we will try to find better ones using a barostat to rescale them.
            # The barostat pressure is not physically meaningful but was chosen to yield a density of pure water close to 1 g/cm^3.

            barostat = MonteCarloBarostat(100, 0, 1)
            system.addForce(barostat)
            integrator = CustomIntegrator(0)
            integrator.addUpdateContextState()
            integrator.addConstrainPositions()

            # Add a dummy force to make the barostat treat any fixed atoms as a molecule so they won't be rescaled.

            if self.topology.getNumAtoms() > 1:
                bondForce = CustomBondForce("0")
                for index in range(1, self.topology.getNumAtoms()):
                    bondForce.addBond(0, index)
                system.addForce(bondForce)

        else:
            integrator = VerletIntegrator(0)

        # Solve the packing problem and update the Modeller with the output.

        force.setCutoffDistance(2 * max(force.getParticleParameters(index)[0] for index in range(force.getNumParticles())))
        if platform is None:
            context = Context(system, integrator)
        else:
            context = Context(system, integrator, platform)
        context.setPositions(np.concatenate(positions))

        minimizeTol = 0.01
        LocalEnergyMinimizer.minimize(context, tolerance=minimizeTol)
        context.setParameter("lambda", 1)
        LocalEnergyMinimizer.minimize(context, tolerance=minimizeTol)

        if vectorsGuessed:
            currentVol = context.getState().getPeriodicBoxVolume().value_in_unit(nanometer**3)
            for _ in range(100):
                integrator.step(100)
                LocalEnergyMinimizer.minimize(context, tolerance=minimizeTol)
                newVol = context.getState().getPeriodicBoxVolume().value_in_unit(nanometer**3)
                if abs(newVol - currentVol) <= 0.001 * currentVol:
                    break
                currentVol = newVol

            topology.setPeriodicBoxVectors(context.getState().getPeriodicBoxVectors())

        self.topology = topology
        self.positions = context.getState(positions=True).getPositions()

    def _computeBoxWidth(self, padding):
        """Compute the necessary periodic box width for a given padding distance."""
        if is_quantity(padding):
            padding = padding.value_in_unit(nanometer)
        if len(self.positions) == 0:
            radius = 0
        else:
            positions = self.positions.value_in_unit(nanometer)
            minRange = Vec3(*(min((pos[i] for pos in positions)) for i in range(3)))
            maxRange = Vec3(*(max((pos[i] for pos in positions)) for i in range(3)))
            center = 0.5*(minRange+maxRange)
            radius = max(unit.norm(center-pos) for pos in positions)
        return max(2*radius+padding, 2*padding)

    def _computeBoxVectors(self, width, boxShape):
        """Compute the periodic box vectors given a box width and shape."""
        if boxShape == 'cube':
            return (Vec3(width, 0, 0), Vec3(0, width, 0), Vec3(0, 0, width))
        elif boxShape == 'dodecahedron':
            return (Vec3(width, 0, 0), Vec3(0, width, 0), Vec3(0.5, 0.5, 0.5*sqrt(2))*width)
        elif boxShape == 'octahedron':
            return (Vec3(width, 0, 0), Vec3(1/3, 2*sqrt(2)/3, 0)*width, Vec3(-1/3, sqrt(2)/3, sqrt(6)/3)*width)
        else:
            raise ValueError(f'Illegal box shape: {boxShape}')

    class _ResidueData:
        """Inner class used to encapsulate data about the hydrogens for a residue."""
        def __init__(self, name):
            self.name = name
            self.variants = []
            self.hydrogens = []

    class _Hydrogen:
        """Inner class used to encapsulate data about a hydrogen atom."""
        def __init__(self, name, parent, maxph, variants, terminal):
            self.name = name
            self.parent = parent
            self.maxph = maxph
            self.variants = variants
            self.terminal = terminal

    @staticmethod
    def loadHydrogenDefinitions(file):
        """Load an XML file containing definitions of hydrogens that should be added by addHydrogens().

        The built in hydrogens.xml file containing definitions for standard amino acids and nucleotides is loaded automatically.
        This method can be used to load additional definitions for other residue types.  They will then be used in subsequent
        calls to addHydrogens().

        Parameters
        ----------
        file : string or file
            An XML file containing hydrogen definitions.  It may be either an
            absolute file path, a path relative to the current working
            directory, a path relative to this module's data subdirectory (for
            built in sets of definitions), or an open file-like object with a read()
            method from which the data can be loaded.
        """
        tree = None
        try:
            # this handles either filenames or open file-like objects
            tree = etree.parse(file)
        except IOError:
            for dataDir in _getDataDirectories():
                f = os.path.join(dataDir, file)
                if os.path.isfile(f):
                    tree = etree.parse(f)
                    break
        if tree is None:
            raise ValueError('Could not locate file')
        infinity = float('Inf')
        for residue in tree.getroot().findall('Residue'):
            resName = residue.attrib['name']
            data = Modeller._ResidueData(resName)
            Modeller._residueHydrogens[resName] = data
            for variant in residue.findall('Variant'):
                data.variants.append(variant.attrib['name'])
            for hydrogen in residue.findall('H'):
                maxph = infinity
                if 'maxph' in hydrogen.attrib:
                    maxph = float(hydrogen.attrib['maxph'])
                atomVariants = None
                if 'variant' in hydrogen.attrib:
                    atomVariants = hydrogen.attrib['variant'].split(',')
                terminal = None
                if 'terminal' in hydrogen.attrib:
                    terminal = hydrogen.attrib['terminal']
                data.hydrogens.append(Modeller._Hydrogen(hydrogen.attrib['name'], hydrogen.attrib['parent'], maxph, atomVariants, terminal))

    @staticmethod
    def _loadStandardHydrogenDefinitions():
        """Load the definitions of hydrogens for standard residues.  Normally there is no need to call this directly.
        It is automatically called by addHydrogens().  If the definitions have already been loaded, this returns without
        doing anything."""
        if not Modeller._hasLoadedStandardHydrogens:
            Modeller.loadHydrogenDefinitions(os.path.join(os.path.dirname(__file__), 'data', 'hydrogens.xml'))
            Modeller._hasLoadedStandardHydrogens = True

    def addHydrogens(self, forcefield=None, pH=7.0, variants=None, platform=None, residueTemplates=dict()):
        """Add missing hydrogens to the model.

        Some residues can exist in multiple forms depending on the pH and properties of the local environment.  These
        variants differ in the presence or absence of particular hydrogens.  In particular, the following variants
        are supported:

        Aspartic acid:
            ASH: Neutral form with a hydrogen on one of the delta oxygens
            ASP: Negatively charged form without a hydrogen on either delta oxygen

        Cysteine:
            CYS: Neutral form with a hydrogen on the sulfur
            CYX: No hydrogen on the sulfur (either negatively charged, or part of a disulfide bond)

        Glutamic acid:
            GLH: Neutral form with a hydrogen on one of the epsilon oxygens
            GLU: Negatively charged form without a hydrogen on either epsilon oxygen

        Histidine:
            HID: Neutral form with a hydrogen on the ND1 atom
            HIE: Neutral form with a hydrogen on the NE2 atom
            HIP: Positively charged form with hydrogens on both ND1 and NE2
            HIN: Negatively charged form without a hydrogen on either ND1 or NE2

        Lysine:
            LYN: Neutral form with two hydrogens on the zeta nitrogen
            LYS: Positively charged form with three hydrogens on the zeta nitrogen

        The variant to use for each residue is determined by the following rules:

        1. The most common variant at the specified pH is selected.
        2. Any Cysteine that participates in a disulfide bond uses the CYX variant regardless of pH.
        3. For a neutral Histidine residue, the HID or HIE variant is selected based on which one forms a better hydrogen bond.

        You can override these rules by explicitly specifying a variant for any residue.  To do that, provide a list for the
        'variants' parameter, and set the corresponding element to the name of the variant to use.

        A special case is when the model already contains a hydrogen that should not be present in the desired variant.
        If you explicitly specify a variant using the 'variants' parameter, the residue will be modified to match the
        desired variant, removing hydrogens if necessary.  On the other hand, for residues whose variant is selected
        automatically, this function will only add hydrogens.  It will never remove ones that are already present in the
        model, regardless of the specified pH.

        In all cases, the positions of existing atoms (including existing hydrogens) are not modified.

        Definitions for standard amino acids and nucleotides are built in.  You can call loadHydrogenDefinitions() to load
        additional definitions for other residue types.

        Parameters
        ----------
        forcefield : ForceField=None
            the ForceField to use for determining the positions of hydrogens.
            If this is None, positions will be picked which are generally
            reasonable but not optimized for any particular ForceField.
        pH : float=7.0
            the pH based on which to select variants
        variants : list=None
            an optional list of variants to use.  If this is specified, its
            length must equal the number of residues in the model.  variants[i]
            is the name of the variant to use for residue i (indexed starting at
            0). If an element is None, the standard rules will be followed to
            select a variant for that residue.  Alternatively, an element may specify
            exactly which hydrogens to add.  In that case, variants[i] should be
            a list of tuples [(name1, parent1), (name2, parent2), ...].  Each
            tuple specifies the name of a hydrogen and the name of the parent atom
            it should be bonded to.
        platform : Platform=None
            the Platform to use when computing the hydrogen atom positions.  If
            this is None, the default Platform will be used.
        residueTemplates : dict=dict()
            specifies which template the ForceField should use for particular residues.  The keys
            should be Residue objects from the Topology, and the values should be the names of the
            templates to use for them.  This is useful when a ForceField contains multiple templates
            that can match the same residue (e.g Fe2+ and Fe3+ templates in the ForceField for a
            monoatomic iron ion in the Topology).

        Returns
        -------
        list
             a list of what variant was actually selected for each residue,
             in the same format as the variants parameter
        """
        # Check the list of variants.

        residues = list(self.topology.residues())
        if variants is not None:
            if len(variants) != len(residues):
                raise ValueError("The length of the variants list must equal the number of residues")
        else:
            variants = [None]*len(residues)
        actualVariants = [None]*len(residues)

        # Load the residue specifications.

        Modeller._loadStandardHydrogenDefinitions()

        # Make a list of atoms bonded to each atom.

        bonded = {}
        for atom in self.topology.atoms():
            bonded[atom] = []
        for atom1, atom2 in self.topology.bonds():
            bonded[atom1].append(atom2)
            bonded[atom2].append(atom1)

        # Define a function that decides whether a set of atoms form a hydrogen bond, using fairly tolerant criteria.

        def isHbond(d, h, a):
            if norm(d-a) > 0.35*nanometer:
                return False
            deltaDH = h-d
            deltaHA = a-h
            deltaDH /= norm(deltaDH)
            deltaHA /= norm(deltaHA)
            return acos(dot(deltaDH, deltaHA)) < 50*degree

        # Loop over residues.

        newTopology = Topology()
        newTopology.setPeriodicBoxVectors(self.topology.getPeriodicBoxVectors())
        newAtoms = {}
        newPositions = []*nanometer
        newResidueTemplates = {}
        newIndices = []
        acceptors = [atom for atom in self.topology.atoms() if atom.element in (elem.oxygen, elem.nitrogen)]
        positions = self.positions.value_in_unit(nanometer)
        acceptorPositions = [positions[a.index] for a in acceptors]
        cells = compiled.CellList(acceptorPositions, 0.35, None, False)
        for chain in self.topology.chains():
            newChain = newTopology.addChain(chain.id)
            for residue in chain.residues():
                newResidue = newTopology.addResidue(residue.name, newChain, residue.id, residue.insertionCode)
                if residue in residueTemplates:
                    newResidueTemplates[newResidue] = residueTemplates[residue]
                isNTerminal = (residue == chain._residues[0])
                isCTerminal = (residue == chain._residues[-1])
                if residue.name in Modeller._residueHydrogens or isinstance(variants[residue.index], list):
                    # Add hydrogens.  First select which variant to use.

                    variant = variants[residue.index]
                    if variant is None:
                        if residue.name == 'CYS':
                            # If this is part of a disulfide, use CYX.

                            sulfur = [atom for atom in residue.atoms() if atom.element == elem.sulfur]
                            if len(sulfur) == 1 and any((atom.residue != residue for atom in bonded[sulfur[0]])):
                                variant = 'CYX'
                        if residue.name == 'HIS' and pH > 6.5:
                            # See if either nitrogen already has a hydrogen attached.

                            nd1 = [atom for atom in residue.atoms() if atom.name == 'ND1']
                            ne2 = [atom for atom in residue.atoms() if atom.name == 'NE2']
                            if len(nd1) != 1 or len(ne2) != 1:
                                raise ValueError('HIS residue (%d) has the wrong set of atoms' % residue.index)
                            nd1 = nd1[0]
                            ne2 = ne2[0]
                            nd1HasHydrogen = any((atom.element == elem.hydrogen for atom in bonded[nd1]))
                            ne2HasHydrogen = any((atom.element == elem.hydrogen for atom in bonded[ne2]))
                            if nd1HasHydrogen and ne2HasHydrogen:
                                variant = 'HIP'
                            elif nd1HasHydrogen:
                                variant = 'HID'
                            elif ne2HasHydrogen:
                                variant = 'HIE'
                            else:
                                # Estimate the hydrogen positions.

                                nd1Pos = self.positions[nd1.index]
                                ne2Pos = self.positions[ne2.index]
                                hd1Delta = Vec3(0, 0, 0)*nanometer
                                for other in bonded[nd1]:
                                    hd1Delta += nd1Pos-self.positions[other.index]
                                hd1Delta *= 0.1*nanometer/norm(hd1Delta)
                                hd1Pos = nd1Pos+hd1Delta
                                he2Delta = Vec3(0, 0, 0)*nanometer
                                for other in bonded[ne2]:
                                    he2Delta += ne2Pos-self.positions[other.index]
                                he2Delta *= 0.1*nanometer/norm(he2Delta)
                                he2Pos = ne2Pos+he2Delta

                                # See whether either hydrogen would form a hydrogen bond.

                                nd1IsBonded = False
                                ne2IsBonded = False
                                for acceptorIndex in cells.neighbors(nd1Pos.value_in_unit(nanometer)):
                                    acceptor = acceptors[acceptorIndex]
                                    if acceptor.residue != residue:
                                        acceptorPos = self.positions[acceptor.index]
                                        if isHbond(nd1Pos, hd1Pos, acceptorPos):
                                            nd1IsBonded = True
                                            break
                                if not nd1IsBonded:
                                    for acceptorIndex in cells.neighbors(ne2Pos.value_in_unit(nanometer)):
                                        acceptor = acceptors[acceptorIndex]
                                        if acceptor.residue != residue:
                                            acceptorPos = self.positions[acceptor.index]
                                            if isHbond(ne2Pos, he2Pos, acceptorPos):
                                                ne2IsBonded = True
                                                break
                                if ne2IsBonded and not nd1IsBonded:
                                    variant = 'HIE'
                                else:
                                    variant = 'HID'
                        elif residue.name == 'HIS':
                            variant = 'HIP'
                    if isinstance(variant, list):
                        spec = Modeller._ResidueData(residue.name)
                        infinity = float('Inf')
                        spec.hydrogens = [Modeller._Hydrogen(name, parent, infinity, None, None) for name, parent in variant]
                    else:
                        spec = Modeller._residueHydrogens[residue.name]
                        if variant is not None and variant not in spec.variants:
                            raise ValueError('Illegal variant for %s residue: %s' % (residue.name, variant))
                    actualVariants[residue.index] = variant
                    removeExtraHydrogens = (variants[residue.index] is not None)

                    # Make a list of hydrogens that should be present in the residue.

                    parents = [atom for atom in residue.atoms() if atom.element != elem.hydrogen]
                    parentNames = [atom.name for atom in parents]
                    hydrogens = [h for h in spec.hydrogens if (variant is None and pH <= h.maxph) or (h.variants is None and pH <= h.maxph) or (h.variants is not None and variant in h.variants)]
                    hydrogens = [h for h in hydrogens if h.terminal is None or (isNTerminal and 'N' in h.terminal) or (isCTerminal and 'C' in h.terminal) or
                                 ((not isNTerminal) and (not isCTerminal) and '-' in h.terminal)]
                    hydrogens = [h for h in hydrogens if h.parent in parentNames]

                    # Loop over atoms in the residue, adding them to the new topology along with required hydrogens.

                    for parent in residue.atoms():
                        # Check whether this is a hydrogen that should be removed.

                        if removeExtraHydrogens and parent.element == elem.hydrogen and not any(parent.name == h.name for h in hydrogens):
                            continue

                        # Add the atom.

                        newAtom = newTopology.addAtom(parent.name, parent.element, newResidue)
                        newAtoms[parent] = newAtom
                        newPositions.append(deepcopy(self.positions[parent.index]))
                        if parent in parents:
                            # Match expected hydrogens with existing ones and find which ones need to be added.

                            existing = [atom for atom in bonded[parent] if atom.element == elem.hydrogen]
                            expected = [h for h in hydrogens if h.parent == parent.name]
                            if len(existing) < len(expected):
                                # Try to match up existing hydrogens to expected ones.

                                matches = []
                                for e in existing:
                                    match = [h for h in expected if h.name == e.name]
                                    if len(match) > 0:
                                        matches.append(match[0])
                                        expected.remove(match[0])
                                    else:
                                        matches.append(None)

                                # If any hydrogens couldn't be matched by name, just match them arbitrarily.

                                for i in range(len(matches)):
                                    if matches[i] is None:
                                        matches[i] = expected[-1]
                                        expected.remove(expected[-1])

                                # Add the missing hydrogens.

                                for h in expected:
                                    newH = newTopology.addAtom(h.name, elem.hydrogen, newResidue)
                                    newIndices.append(newH.index)
                                    delta = Vec3(0, 0, 0)*nanometer
                                    if len(bonded[parent]) > 0:
                                        for other in bonded[parent]:
                                            delta += self.positions[parent.index]-self.positions[other.index]
                                    else:
                                        delta = Vec3(random.random(), random.random(), random.random())*nanometer
                                    delta *= 0.1*nanometer/norm(delta)
                                    delta += 0.05*Vec3(random.random(), random.random(), random.random())*nanometer
                                    delta *= 0.1*nanometer/norm(delta)
                                    newPositions.append(self.positions[parent.index]+delta)
                                    newTopology.addBond(newAtom, newH)
                else:
                    # Just copy over the residue.

                    for atom in residue.atoms():
                        newAtom = newTopology.addAtom(atom.name, atom.element, newResidue)
                        newAtoms[atom] = newAtom
                        newPositions.append(deepcopy(self.positions[atom.index]))
        for bond in self.topology.bonds():
            if bond[0] in newAtoms and bond[1] in newAtoms:
                newTopology.addBond(newAtoms[bond[0]], newAtoms[bond[1]], bond.type, bond.order)

        # The hydrogens were added at random positions.  Now perform an energy minimization to fix them up.

        addedH = set(newIndices)  # keep track of Hs added

        if forcefield is not None:
            # Use the ForceField the user specified.

            system = forcefield.createSystem(newTopology, rigidWater=False, nonbondedMethod=CutoffNonPeriodic, residueTemplates=newResidueTemplates)
            for i in range(system.getNumParticles()):
                if i not in addedH:
                    # Existing atom, make it immobile.
                    system.setParticleMass(i, 0)
        else:
            # Create a System that restrains the distance of each hydrogen from its parent atom
            # and causes hydrogens to spread out evenly.

            system = System()
            nonbonded = CustomNonbondedForce('100/(r/0.1)^4')
            nonbonded.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic);
            nonbonded.setCutoffDistance(1*nanometer)
            bonds = HarmonicBondForce()
            angles = HarmonicAngleForce()
            system.addForce(nonbonded)
            system.addForce(bonds)
            system.addForce(angles)
            bondedTo = []
            for atom in newTopology.atoms():
                nonbonded.addParticle([])
                if atom.index not in addedH:  # make immobile
                    system.addParticle(0.0)
                else:
                    system.addParticle(1.0)
                bondedTo.append([])
            for atom1, atom2 in newTopology.bonds():
                if atom1.element == elem.hydrogen or atom2.element == elem.hydrogen:
                    bonds.addBond(atom1.index, atom2.index, 0.1, 100000.0)
                bondedTo[atom1.index].append(atom2)
                bondedTo[atom2.index].append(atom1)
            for residue in newTopology.residues():
                if residue.name == 'HOH':
                    # Add an angle term to make the water geometry correct.

                    atoms = list(residue.atoms())
                    oindex = [i for i in range(len(atoms)) if atoms[i].element == elem.oxygen]
                    if len(atoms) == 3 and len(oindex) == 1:
                        hindex = list(set([0,1,2])-set(oindex))
                        angles.addAngle(atoms[hindex[0]].index, atoms[oindex[0]].index, atoms[hindex[1]].index, 1.824, 836.8)
                else:
                    # Add angle terms for any hydroxyls.

                    for atom in residue.atoms():
                        index = atom.index
                        if atom.element == elem.oxygen and len(bondedTo[index]) == 2 and elem.hydrogen in (a.element for a in bondedTo[index]):
                            angles.addAngle(bondedTo[index][0].index, index, bondedTo[index][1].index, 1.894, 460.24)

        if platform is None:
            context = Context(system, VerletIntegrator(0.0))
        else:
            context = Context(system, VerletIntegrator(0.0), platform)
        context.setPositions(newPositions)
        LocalEnergyMinimizer.minimize(context, 1.0, 50)
        self.topology = newTopology
        self.positions = context.getState(positions=True).getPositions()
        del context
        return actualVariants

    def addExtraParticles(self, forcefield, ignoreExternalBonds=False, residueTemplates=dict()):
        """Add missing extra particles to the model that are required by a force
        field.

        Some force fields use "extra particles" that do not represent
        actual atoms, but still need to be included in the System.  Examples
        include lone pairs, Drude particles, and the virtual sites used in some
        water models to adjust the charge distribution.  Extra particles can be
        recognized by the fact that their element is None.

        This method is primarily used to add extra particles, but it can also
        remove them.  It tries to match every residue in the Topology to a
        template in the force field.  If there is no match, it will both add
        and remove extra particles as necessary to make it match.

        Parameters
        ----------
        forcefield : ForceField
            the ForceField defining what extra particles should be present
        ignoreExternalBonds : boolean=False
            If true, ignore external bonds when matching residues to templates.
            This is useful when the Topology represents one piece of a larger
            molecule, so chains are not terminated properly.
        residueTemplates : dict=dict()
            specifies which template the ForceField should use for particular residues.  The keys
            should be Residue objects from the Topology, and the values should be the names of the
            templates to use for them.  This is useful when a ForceField contains multiple templates
            that can match the same residue (e.g Fe2+ and Fe3+ templates in the ForceField for a
            monoatomic iron ion in the Topology).
        """
        # Record which atoms are bonded to each other atom.

        bondedToAtom = [set() for _ in self.topology.atoms()]
        for atom1, atom2 in self.topology.bonds():
            bondedToAtom[atom1.index].add(atom2.index)
            bondedToAtom[atom2.index].add(atom1.index)

        # If the force field has a DrudeForce, record the types of Drude particles and their parents since we'll
        # need them for picking particle positions.

        drudeTypeMap = {}
        for force in forcefield._forces:
            if isinstance(force, DrudeGenerator):
                for type in force.typeMap:
                    drudeTypeMap[type] = force.typeMap[type][0]

        # Identify the template to use for each residue.

        templates = forcefield._matchAllResiduesToTemplates(ForceField._SystemData(self.topology), self.topology, residueTemplates, False, True, False)

        # Identify extra points to add.  Record information about each one we will need when adding it.

        ExtraPoint = namedtuple('ExtraPoint', ['matchedResidue', 'template', 'templateAtom', 'matchingAtoms', 'templateAtomPositions'])
        extraPoints = defaultdict(list)
        from openmm.app.topology import MergedResidue
        for matchedResidue, template in templates.items():
            if len(template.atoms) != len(list(matchedResidue.atoms())):
                matches = compiled.matchResidueToTemplate(matchedResidue, template, bondedToAtom, ignoreExternalBonds, True)
                atomsNoEP = [a for a in matchedResidue.atoms() if a.element is not None]
                templateAtomsNoEP = [a for a in template.atoms if a.element is not None]
                matchingAtoms = {}
                for atom, match in zip(atomsNoEP, matches):
                    templateAtomName = templateAtomsNoEP[match].name
                    for templateAtom in template.atoms:
                        if templateAtom.name == templateAtomName:
                            matchingAtoms[templateAtom] = atom
                templateAtomPositions = len(template.atoms)*[None]
                for index, atom in enumerate(template.atoms):
                    if atom in matchingAtoms:
                        templateAtomPositions[index] = self.positions[matchingAtoms[atom].index].value_in_unit(nanometer)
                for index, atom in enumerate(template.atoms):
                    if atom.element is None:
                        # This is an extra particle.  Figure out what residue it belongs to.
                        residue = None
                        if not isinstance(matchedResidue, MergedResidue):
                            residue = matchedResidue
                        else:
                            if atom.type in drudeTypeMap:
                                # It's a Drude particle.
                                for atom2 in template.atoms:
                                    if atom2.type in drudeTypeMap[atom.type]:
                                        residue = matchingAtoms[atom2].residue
                                        break
                            if residue is None:
                                for site in template.virtualSites:
                                    if site.index == index:
                                        # It's a virtual site.
                                        residue = matchingAtoms[template.atoms[site.atoms[0]]].residue
                                        break
                            if residue is None:
                                # Neither a virtual site nor a Drude particle???  Just guess at the residue.
                                residue = matchedResidue.residues[0]
                        extraPoints[residue].append(ExtraPoint(matchedResidue, template, atom, matchingAtoms, templateAtomPositions))

        # Create the new Topology.

        newTopology = Topology()
        newTopology.setPeriodicBoxVectors(self.topology.getPeriodicBoxVectors())
        newAtoms = {}
        newPositions = []*nanometer
        newResidueTemplates = {}
        missingPositions = set()
        for chain in self.topology.chains():
            newChain = newTopology.addChain(chain.id)
            for residue in chain.residues():
                newResidue = newTopology.addResidue(residue.name, newChain, residue.id, residue.insertionCode)
                if residue in residueTemplates:
                    newResidueTemplates[newResidue] = residueTemplates[residue]
                if residue not in extraPoints:
                    # Just copy the residue over.

                    for atom in residue.atoms():
                        newAtom = newTopology.addAtom(atom.name, atom.element, newResidue)
                        newAtoms[atom] = newAtom
                        newPositions.append(deepcopy(self.positions[atom.index]))
                else:
                    # Add the regular atoms.

                    for atom in residue.atoms():
                        if atom.element is not None:
                            newAtoms[atom] = newTopology.addAtom(atom.name, atom.element, newResidue)
                            newPositions.append(deepcopy(self.positions[atom.index]))

                    # Add the extra points.

                    newExtraPoints = {}
                    for extraPoint in extraPoints[residue]:
                        atom = extraPoint.templateAtom
                        template = extraPoint.template
                        templateAtomPositions = extraPoint.templateAtomPositions
                        newExtraPoints[atom] = newTopology.addAtom(atom.name, None, newResidue)
                        position = None
                        for site in extraPoint.template.virtualSites:
                            if template.atoms[site.index] == atom:
                                # This is a virtual site.  Compute its position by the correct rule.

                                try:
                                    if site.type == 'average2':
                                        position = site.weights[0]*templateAtomPositions[site.atoms[0]] + site.weights[1]*templateAtomPositions[site.atoms[1]]
                                    elif site.type == 'average3':
                                        position = site.weights[0]*templateAtomPositions[site.atoms[0]] + site.weights[1]*templateAtomPositions[site.atoms[1]] + site.weights[2]*templateAtomPositions[site.atoms[2]]
                                    elif site.type == 'outOfPlane':
                                        v1 = templateAtomPositions[site.atoms[1]] - templateAtomPositions[site.atoms[0]]
                                        v2 = templateAtomPositions[site.atoms[2]] - templateAtomPositions[site.atoms[0]]
                                        cross = Vec3(v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0])
                                        position = templateAtomPositions[site.atoms[0]] + site.weights[0]*v1 + site.weights[1]*v2 + site.weights[2]*cross
                                    elif site.type == 'localCoords':
                                        origin = unit.sum([templateAtomPositions[atom]*weight for atom, weight in zip(site.atoms, site.originWeights)])
                                        xdir = unit.sum([templateAtomPositions[atom]*weight for atom, weight in zip(site.atoms, site.xWeights)])
                                        ydir = unit.sum([templateAtomPositions[atom]*weight for atom, weight in zip(site.atoms, site.yWeights)])
                                        zdir = Vec3(xdir[1]*ydir[2]-xdir[2]*ydir[1], xdir[2]*ydir[0]-xdir[0]*ydir[2], xdir[0]*ydir[1]-xdir[1]*ydir[0])
                                        xdir /= norm(xdir)
                                        zdir /= norm(zdir)
                                        ydir = Vec3(zdir[1]*xdir[2]-zdir[2]*xdir[1], zdir[2]*xdir[0]-zdir[0]*xdir[2], zdir[0]*xdir[1]-zdir[1]*xdir[0])
                                        position = origin + xdir*site.localPos[0] + ydir*site.localPos[1] + zdir*site.localPos[2]
                                except:
                                    # This can happen if the virtual site depends on another virtual site whose position
                                    # hasn't been set yet.  Ignore the error.  We'll put it at a random position (see below),
                                    # which will get replaced with the correct position at the start of the simulation.

                                    pass
                        if position is None and atom.type in drudeTypeMap:
                            # This is a Drude particle.  Put it on top of its parent atom.

                            for atom2, pos in zip(template.atoms, templateAtomPositions):
                                if atom2.type in drudeTypeMap[atom.type]:
                                    position = deepcopy(pos)
                        if position is None:
                            # We couldn't figure out the correct position.  Put it at a random position near the center of the residue,
                            # and we'll try to fix it later based on bonds.

                            knownPositions = [x for x in templateAtomPositions if x is not None]
                            position = Vec3(random.gauss(0, 1), random.gauss(0, 1), random.gauss(0, 1))+(unit.sum(knownPositions)/len(knownPositions))
                            missingPositions.add(len(newPositions))
                        newPositions.append(position*nanometer)

                    # Add bonds involving the extra points.

                    template = extraPoints[residue][0].template
                    matchingAtoms = extraPoints[residue][0].matchingAtoms
                    for atom1, atom2 in template.bonds:
                        atom1 = template.atoms[atom1]
                        atom2 = template.atoms[atom2]
                        if atom1 in newExtraPoints or atom2 in newExtraPoints:
                            if atom1 in newExtraPoints:
                                a1 = newExtraPoints[atom1]
                            else:
                                a1 = newAtoms[matchingAtoms[atom1]]
                            if atom2 in newExtraPoints:
                                a2 = newExtraPoints[atom2]
                            else:
                                a2 = newAtoms[matchingAtoms[atom2]]
                            newTopology.addBond(a1, a2)

        for bond in self.topology.bonds():
            if bond[0] in newAtoms and bond[1] in newAtoms:
                newTopology.addBond(newAtoms[bond[0]], newAtoms[bond[1]], bond.type, bond.order)

        if len(missingPositions) > 0:
            # There were particles whose position we couldn't identify before, since they were neither virtual sites nor Drude particles.
            # Try to figure them out based on bonds.  First, use the ForceField to create a list of every bond involving one of them.

            system = forcefield.createSystem(newTopology, constraints=AllBonds, residueTemplates=newResidueTemplates)
            bonds = []
            for i in range(system.getNumConstraints()):
                bond = system.getConstraintParameters(i)
                if bond[0] in missingPositions or bond[1] in missingPositions:
                    bonds.append(bond)

            # Now run a few iterations of SHAKE to try to select reasonable positions.

            for iteration in range(15):
                for atom1, atom2, distance in bonds:
                    if atom1 in missingPositions:
                        if atom2 in missingPositions:
                            weights = (0.5, 0.5)
                        else:
                            weights = (1.0, 0.0)
                    else:
                        weights = (0.0, 1.0)
                    delta = newPositions[atom2]-newPositions[atom1]
                    length = norm(delta)
                    delta *= (distance-length)/length
                    newPositions[atom1] -= weights[0]*delta
                    newPositions[atom2] += weights[1]*delta

        self.topology = newTopology
        self.positions = newPositions


    def addMembrane(self, forcefield, lipidType='POPC', membraneCenterZ=0*nanometer, minimumPadding=1*nanometer, positiveIon='Na+', negativeIon='Cl-',
                    ionicStrength=0*molar, neutralize=True, residueTemplates=dict(), platform=None):
        """Add a lipid membrane to the model.

        This method actually adds both a membrane and a water box.  It is best to build them together,
        both to avoid adding waters inside the membrane and to ensure that lipid head groups are properly
        solvated.  For that reason, this method includes many of the same arguments as addSolvent().

        The membrane is added in the XY plane, and the existing protein is assumed to already be oriented
        and positioned correctly.  When possible, it is recommended to start with a model from the
        Orientations of Proteins in Membranes (OPM) database at http://opm.phar.umich.edu.  Otherwise, it
        is up to you to select the protein position yourself.

        The algorithm is based on the one described in Wolf et al., J. Comp. Chem. 31, pp. 2169-2174 (2010).
        It begins by tiling copies of a pre-equilibrated membrane patch to create a membrane of the desired
        size.  Next it scales down the protein by 50% along the X and Y axes.  Any lipid within a cutoff
        distance of the scaled protein is removed.  It also ensures that equal numbers of lipids are removed
        from each leaf of the membrane.  Finally, it runs molecular dynamics to let the membrane relax while
        gradually scaling the protein back up to its original size.

        The size of the membrane and water box are determined by the minimumPadding argument.  All
        pre-existing atoms are guaranteed to be at least this far from any edge of the periodic box.  It
        is also possible for the periodic box to have more padding than requested.  In particular, it only
        adds whole copies of the pre-equilibrated membrane patch, so the box dimensions will always be
        integer multiples of the patch size.  That may lead to a larger membrane than what you requested.

        This method has built in support for POPC, POPE, DLPC, DLPE, DMPC, DOPC and DPPC lipids.
        You can also build other types of membranes by providing a pre-equilibrated, solvated membrane patch
        that can be tiled in the XY plane to form the membrane.

        Parameters
        ----------
        forcefield : ForceField
            the ForceField to use for determining atomic charges and for relaxing the membrane
        lipidType : string or object
            the type of lipid to use.  Supported string values are 'POPC', 'POPE', 'DLPC', 'DLPE', 'DMPC',
            'DOPC', and 'DPPC'.  For other types of lipids, provide a PDBFile or PDBxFile object (or any
            other object with "topology" and "positions" fields) containing a membrane patch.
        membraneCenterZ: distance=0*nanometer
            the position along the Z axis of the center of the membrane
        minimumPadding : distance=1*nanometer
            the padding distance to use
        positiveIon : string='Na+'
            the type of positive ion to add.  Allowed values are 'Cs+', 'K+', 'Li+', 'Na+', and 'Rb+'
        negativeIon : string='Cl-'
            the type of negative ion to add.  Allowed values are 'Cl-', 'Br-', 'F-', and 'I-'. Be aware
            that not all force fields support all ion types.
        ionicStrength : concentration=0*molar
            the total concentration of ions (both positive and negative) to add.  This
            does not include ions that are added to neutralize the system.
            Note that only monovalent ions are currently supported.
        neutralize : bool=True
            whether to add ions to neutralize the system
        residueTemplates : dict=dict()
            specifies which template the ForceField should use for particular residues.  The keys
            should be Residue objects from the Topology, and the values should be the names of the
            templates to use for them.  This is useful when a ForceField contains multiple templates
            that can match the same residue (e.g Fe2+ and Fe3+ templates in the ForceField for a
            monoatomic iron ion in the Topology).
        platform : Platform=None
            the Platform to use when computing the hydrogen atom positions.  If
            this is None, the default Platform will be used.
        """
        if 'topology' in dir(lipidType) and 'positions' in dir(lipidType):
            patch = lipidType
        elif lipidType.upper() in ('POPC', 'POPE', 'DLPC', 'DLPE', 'DMPC', 'DOPC', 'DPPC'):
            patch = PDBFile(os.path.join(os.path.dirname(__file__), 'data', lipidType.upper()+'.pdb'))
        else:
            raise ValueError('Unsupported lipid type: '+lipidType)
        if is_quantity(membraneCenterZ):
            membraneCenterZ = membraneCenterZ.value_in_unit(nanometer)
        if is_quantity(minimumPadding):
            minimumPadding = minimumPadding.value_in_unit(nanometer)

        # Figure out how many copies of the membrane patch we need in each direction.

        proteinPos = self.positions.value_in_unit(nanometer)
        proteinMinPos = Vec3(*[min((p[i] for p in proteinPos)) for i in range(3)])
        proteinMaxPos = Vec3(*[max((p[i] for p in proteinPos)) for i in range(3)])
        proteinSize = proteinMaxPos-proteinMinPos
        proteinCenterPos = (proteinMinPos+proteinMaxPos)/2
        proteinCenterPos = Vec3(proteinCenterPos[0], proteinCenterPos[1], membraneCenterZ)
        patchPos = patch.positions.value_in_unit(nanometer)
        patchSize = patch.topology.getUnitCellDimensions().value_in_unit(nanometer)
        patchMinPos = Vec3(*[min((p[i] for p in patchPos)) for i in range(3)])
        patchMaxPos = Vec3(*[max((p[i] for p in patchPos)) for i in range(3)])
        patchCenterPos = (patchMinPos+patchMaxPos)/2
        nx = int(ceil((proteinSize[0]+2*minimumPadding)/patchSize[0]))
        ny = int(ceil((proteinSize[1]+2*minimumPadding)/patchSize[1]))

        # Record the bonds for each residue.

        resBonds = defaultdict(list)
        for bond in patch.topology.bonds():
            resBonds[bond[0].residue].append(bond)

        # Identify which leaf of the membrane each lipid is in.

        numLipidAtoms = 0
        resMeanZ = {}
        membraneMeanZ = 0.0
        for res in patch.topology.residues():
            if res.name != 'HOH':
                numResAtoms = 0
                sumZ = 0.0
                for atom in res.atoms():
                    numResAtoms += 1
                    sumZ += patchPos[atom.index][2]
                numLipidAtoms += numResAtoms
                membraneMeanZ += sumZ
                resMeanZ[res] = sumZ/numResAtoms
        membraneMeanZ /= numLipidAtoms
        lipidLeaf = dict((res, 0 if resMeanZ[res] < membraneMeanZ else 1) for res in resMeanZ)

        # Compute scaled positions for the protein.

        scaledProteinPos = [None]*len(proteinPos)
        for i, p in enumerate(proteinPos):
            p = p-proteinCenterPos
            p = Vec3(0.5*p[0], 0.5*p[1], p[2])
            scaledProteinPos[i] = p+proteinCenterPos

        # Create a new Topology for the membrane.

        membraneTopology = Topology()
        membranePos = []
        boxSizeZ = patchSize[2]
        if self.topology.getUnitCellDimensions() is not None:
            boxSizeZ = max(boxSizeZ, self.topology.getUnitCellDimensions()[2].value_in_unit(nanometer)+2*minimumPadding)
        else:
            boxSizeZ = max(boxSizeZ, proteinSize[2]+2*minimumPadding)
        membraneTopology.setUnitCellDimensions((nx*patchSize[0], ny*patchSize[1], boxSizeZ))

        # Add membrane patches.  We exclude any water that is within a cutoff distance of either the actual or scaled
        # protein, and any lipid that is within a cutoff distance of the scaled protein.  We also keep track of how
        # many lipids have been excluded from each leaf of the membrane, so we can make sure exactly the same
        # number get removed from each leaf.

        overlapCutoff = 0.22
        addedWater = []
        addedLipids = []
        removedFromLeaf = [0, 0]
        vectors = membraneTopology.getPeriodicBoxVectors().value_in_unit(nanometer)
        proteinCells = compiled.CellList(proteinPos, overlapCutoff, vectors, False)
        scaledProteinCells = compiled.CellList(scaledProteinPos, overlapCutoff, vectors, False)
        for x in range(nx):
            for y in range(ny):
                offset = proteinCenterPos - patchCenterPos + Vec3((x-0.5*(nx-1))*patchSize[0], (y-0.5*(ny-1))*patchSize[1], 0)
                for res in patch.topology.residues():
                    resPos = [patchPos[atom.index]+offset for atom in res.atoms()]
                    if res.name == 'HOH':
                        # Remove waters that are too close to either the original OR scaled protein positions.
                        referencePosLists = [proteinPos, scaledProteinPos]
                        cellLists = [proteinCells, scaledProteinCells]
                    else:
                        # Remove lipids that are too close to the scaled protein positions.
                        referencePosLists = [scaledProteinPos]
                        cellLists = [scaledProteinCells]
                    overlap = False
                    nearest = nx*patchSize[0]
                    for cells, referencePos in zip(cellLists, referencePosLists):
                        if overlap:
                            break
                        for index, atom in enumerate(res.atoms()):
                            pos = resPos[index]
                            for atom in cells.neighbors(pos):
                                distance = norm(pos-referencePos[atom])
                                if distance < overlapCutoff:
                                    overlap = True
                                    break
                                nearest = min(nearest, distance)
                            if overlap:
                                break
                    if res.name == 'HOH':
                        if not overlap:
                            addedWater.append((res, resPos))
                    else:
                        if overlap:
                            removedFromLeaf[lipidLeaf[res]] += 1
                        else:
                            addedLipids.append((nearest, res, resPos))
        skipFromLeaf = [max(removedFromLeaf)-removedFromLeaf[i] for i in (0,1)]
        del cellLists
        del cells
        del proteinCells

        # Add the lipids.

        newAtoms = {}
        lipidChain = membraneTopology.addChain()
        lipidResNum = 1  # renumber lipid residues to handle large patches
        for (nearest, residue, pos) in addedLipids:
            if skipFromLeaf[lipidLeaf[residue]] > 0:
                # Remove the same number of residues from each leaf.
                skipFromLeaf[lipidLeaf[residue]] -= 1
            else:
                newResidue = membraneTopology.addResidue(residue.name, lipidChain, str(lipidResNum), residue.insertionCode)
                lipidResNum += 1

                for atom in residue.atoms():
                    newAtom = membraneTopology.addAtom(atom.name, atom.element, newResidue, atom.id, atom.formalCharge)
                    newAtoms[atom] = newAtom
                membranePos += pos
                for bond in resBonds[residue]:
                    membraneTopology.addBond(newAtoms[bond[0]], newAtoms[bond[1]], bond.type, bond.order)

        del lipidLeaf
        del addedLipids

        # Add the solvent.

        solventChain = membraneTopology.addChain()
        for (residue, pos) in addedWater:
            newResidue = membraneTopology.addResidue(residue.name, solventChain, residue.id, residue.insertionCode)
            for atom in residue.atoms():
                newAtom = membraneTopology.addAtom(atom.name, atom.element, newResidue, atom.id, atom.formalCharge)
                newAtoms[atom] = newAtom
            membranePos += pos
            for bond in resBonds[residue]:
                membraneTopology.addBond(newAtoms[bond[0]], newAtoms[bond[1]], bond.type, bond.order)
        del newAtoms
        del addedWater
        del resBonds
        gc.collect()

        # Create a System for the lipids, then add in the protein as stationary particles.

        system = forcefield.createSystem(membraneTopology, nonbondedMethod=CutoffPeriodic, residueTemplates=residueTemplates)
        proteinSystem = forcefield.createSystem(self.topology, nonbondedMethod=CutoffNonPeriodic, residueTemplates=residueTemplates)
        numMembraneParticles = system.getNumParticles()
        numProteinParticles = proteinSystem.getNumParticles()
        for i in range(numProteinParticles):
            system.addParticle(0.0)
        nonbonded = None
        for f1, f2 in zip(system.getForces(), proteinSystem.getForces()):
            if isinstance(f1, NonbondedForce):
                nonbonded = f2
                for i in range(numProteinParticles):
                    f1.addParticle(*f2.getParticleParameters(i))
                    for j in scaledProteinCells.neighbors(scaledProteinPos[i]):
                        if j < i:
                            f1.addException(i+numMembraneParticles, j+numMembraneParticles, 0.0, 1.0, 0.0)
            elif isinstance(f1, CustomNonbondedForce):
                for i in range(numProteinParticles):
                    f1.addParticle(f2.getParticleParameters(i))
                    for j in scaledProteinCells.neighbors(scaledProteinPos[i]):
                        if j < i:
                            f1.addExclusion(i + numMembraneParticles, j + numMembraneParticles)
        if nonbonded is None:
            raise ValueError('The ForceField does not specify a NonbondedForce')
        mergedPositions = membranePos+scaledProteinPos
        del membranePos
        del scaledProteinCells
        gc.collect()

        # Run a simulation while slowly scaling up the protein so the membrane can relax.
        # Select the number of steps to ensure no atom will ever move more than 0.25 A
        # in one step.

        steps = int(max(proteinSize.x, proteinSize.y)*10) + 1
        integrator = LangevinIntegrator(10.0, 50.0, 0.001)
        if platform is None:
            context = Context(system, integrator)
        else:
            context = Context(system, integrator, platform)
        context.setPositions(mergedPositions)
        LocalEnergyMinimizer.minimize(context, 10.0, 30)
        proteinPosArray = np.array(proteinPos)
        scaledProteinPosArray = np.array(scaledProteinPos)
        for i in range(steps):
            weight1 = i/(steps-1)
            weight2 = 1.0-weight1
            mergedPositions = context.getState(positions=True).getPositions(asNumpy=True).value_in_unit(nanometer)
            mergedPositions[numMembraneParticles:] = weight1*proteinPosArray + weight2*scaledProteinPosArray
            context.setPositions(mergedPositions)
            integrator.step(20)

        # Add the membrane to the protein.

        modeller = Modeller(self.topology, self.positions)
        modeller.add(membraneTopology, context.getState(positions=True).getPositions()[:numMembraneParticles])
        modeller.topology.setPeriodicBoxVectors(membraneTopology.getPeriodicBoxVectors())
        del context
        del system
        del integrator

        # Depending on the box size, we may need to add more water beyond what was included with the membrane patch.

        needExtraWater = (boxSizeZ > patchSize[2])
        if needExtraWater:
            newResidueTemplates = {}
            for r1, r2 in zip(self.topology.residues(), modeller.topology.residues()):
                if r1 in residueTemplates:
                    newResidueTemplates[r2] = residueTemplates[r1]
            modeller.addSolvent(forcefield, neutralize=False, residueTemplates=newResidueTemplates)

        # Record the positions of all waters that have been added.

        waterPos = {}
        for chain in list(modeller.topology.chains())[-2:]:
            for residue in chain.residues():
                if residue.name == 'HOH':
                    for atom in residue.atoms():
                        if atom.element == elem.oxygen:
                            waterPos[residue] = modeller.positions[atom.index].value_in_unit(nanometer)

        # We may have added extra water molecules inside the membrane.  We really only wanted to extend the box
        # without adding more water in the existing box, so remove the unwanted ones.
        if needExtraWater:
            toDelete = []
            addedChain = list(modeller.topology.chains())[-1]
            patchMinZ = patchMinPos[2]-patchCenterPos[2]+membraneCenterZ
            patchMaxZ = patchMaxPos[2]-patchCenterPos[2]+membraneCenterZ
            for residue in addedChain.residues():
                z = waterPos[residue][2]
                if z > patchMinZ and z < patchMaxZ:
                    toDelete.append(residue)
                    del waterPos[residue]
            if len(toDelete) > 0:
                modeller.delete(toDelete)

        newResidueTemplates = {}
        for r1, r2 in zip(self.topology.residues(), modeller.topology.residues()):
            if r1 in residueTemplates:
                newResidueTemplates[r2] = residueTemplates[r1]
        self.topology = modeller.topology
        self.positions = modeller.positions

        # Select a subset of water molecules to replace with ions, ignoring
        # those within a certain distance from either leaflet of the membrane.

        waterPos = {}  # redo because modeller.delete changes chain indexes
        for chain in list(modeller.topology.chains())[-2:]:
            for residue in chain.residues():
                if residue.name == 'HOH':
                    for atom in residue.atoms():
                        if atom.element == elem.oxygen:
                            waterPos[residue] = modeller.positions[atom.index]

        # Total number of water molecules
        # Use this number to avoid underestimating the concentration of ions
        # in _addIons after we exclude waters close to lipids.
        numTotalWaters = len(waterPos)

        # Calculate lipid Z boundaries
        lipidNames = {res.name for res in patch.topology.residues() if res.name != 'HOH'}
        lipidZMax = sys.float_info.min
        lipidZMin = sys.float_info.max
        for res in modeller.topology.residues():
            if res.name in lipidNames:
                for atom in res.atoms():
                    atomZ = modeller.positions[atom.index][2].value_in_unit(nanometer)
                    lipidZMax = max(lipidZMax, atomZ)
                    lipidZMin = min(lipidZMin, atomZ)

        # Ignore waters that are within a certain distance of the membrane
        lipidOffset = 0.25
        upperZBoundary = (lipidZMax + lipidOffset)
        lowerZBoundary = (lipidZMin - lipidOffset)
        waterResidues = list(waterPos)
        for wRes in waterResidues:
            waterZ = waterPos[wRes][2]
            if lowerZBoundary < waterZ.value_in_unit(nanometer) < upperZBoundary:
                del waterPos[wRes]

        self._addIons(forcefield, numTotalWaters, waterPos, positiveIon=positiveIon, negativeIon=negativeIon, ionicStrength=ionicStrength, neutralize=neutralize, residueTemplates=newResidueTemplates)
