"""
modeller.py: Provides tools for editing molecular models

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012 Stanford University and the Authors.
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

__author__ = "Peter Eastman"
__version__ = "1.0"

from simtk.openmm.app import Topology, PDBFile
from simtk.openmm.app.forcefield import HAngles
from simtk.openmm.vec3 import Vec3
from simtk.openmm import System, Context, NonbondedForce, VerletIntegrator, LocalEnergyMinimizer
from simtk.unit import nanometer, molar, elementary_charge, amu, gram, liter, degree, sqrt, acos, is_quantity, dot, norm
import element as elem
import os
import random
import xml.etree.ElementTree as etree
from copy import deepcopy
from math import ceil, floor

class Modeller(object):
    """Modeller provides tools for editing molecular models, such as adding water or missing hydrogens.
    
    To use it, create a Modeller object, specifying the initial Topology and atom positions.  You can
    then call various methods to change the model in different ways.  Each time you do, a new Topology
    and list of coordinates is created to represent the changed model.  Finally, call getTopology()
    and getPositions() to get the results.
    """
        
    _residueHydrogens = None

    def __init__(self, topology, positions):
        """Create a new Modeller object
        
        Parameters:
         - topology (Topology) the initial Topology of the model
         - positions (list) the initial atomic positions
        """
        ## The Topology describing the structure of the system
        self.topology = topology
        if not is_quantity(positions):
            positions = positions*nanometers
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
        
        Specify what to add by providing a new Topology object and the corresponding atomic positions.
        All chains, residues, atoms, and bonds contained in the Topology are added to the model.
        
        Parameters:
         - addTopoology (Topology) a Topology whose contents should be added to the model
         - addPositions (list) the positions of the atoms to add
        """
        # Copy over the existing model.
        
        newTopology = Topology()
        newTopology.setUnitCellDimensions(deepcopy(self.topology.getUnitCellDimensions()))
        newAtoms = {}
        newPositions = []*nanometer
        for chain in self.topology.chains():
            newChain = newTopology.addChain()
            for residue in chain.residues():
                newResidue = newTopology.addResidue(residue.name, newChain)
                for atom in residue.atoms():
                    newAtom = newTopology.addAtom(atom.name, atom.element, newResidue)
                    newAtoms[atom] = newAtom
                    newPositions.append(deepcopy(self.positions[atom.index]))
        for bond in self.topology.bonds():
            newTopology.addBond(newAtoms[bond[0]], newAtoms[bond[1]])
        
        # Add the new model

        newAtoms = {}
        for chain in addTopology.chains():
            newChain = newTopology.addChain()
            for residue in chain.residues():
                newResidue = newTopology.addResidue(residue.name, newChain)
                for atom in residue.atoms():
                    newAtom = newTopology.addAtom(atom.name, atom.element, newResidue)
                    newAtoms[atom] = newAtom
                    newPositions.append(deepcopy(addPositions[atom.index]))
        for bond in addTopology.bonds():
            newTopology.addBond(newAtoms[bond[0]], newAtoms[bond[1]])
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
        
        Parameters:
         - toDelete (list) a list of Atoms, Residues, Chains, and bonds (specified as tuples of Atoms) to delete
        """
        newTopology = Topology()
        newTopology.setUnitCellDimensions(deepcopy(self.topology.getUnitCellDimensions()))
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
                                    newChain = newTopology.addChain()
                                    needNewChain = False;
                                if needNewResidue:
                                    newResidue = newTopology.addResidue(residue.name, newChain)
                                    needNewResidue = False;
                                newAtom = newTopology.addAtom(atom.name, atom.element, newResidue)
                                newAtoms[atom] = newAtom
                                newPositions.append(deepcopy(self.positions[atom.index]))
        for bond in self.topology.bonds():
            if bond[0] in newAtoms and bond[1] in newAtoms:
                if bond not in deleteSet and (bond[1], bond[0]) not in deleteSet:
                    newTopology.addBond(newAtoms[bond[0]], newAtoms[bond[1]])
        self.topology = newTopology
        self.positions = newPositions
    
    def deleteWater(self):
        """Delete all water molecules from the model."""
        self.delete(res for res in self.topology.residues() if res.name == "HOH")
    
    def convertWater(self, model='tip3p'):
        """Convert all water molecules to a different water model.
        
        Parameters:
         - model (string='tip3p') the water model to convert to.  Supported values are 'tip3p', 'spce', 'tip4pew', and 'tip5p'.
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
        newTopology.setUnitCellDimensions(deepcopy(self.topology.getUnitCellDimensions()))
        newAtoms = {}
        newPositions = []*nanometer
        for chain in self.topology.chains():
            newChain = newTopology.addChain()
            for residue in chain.residues():
                newResidue = newTopology.addResidue(residue.name, newChain)
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
                        newAtom = newTopology.addAtom(atom.name, atom.element, newResidue)
                        newAtoms[atom] = newAtom
                        newPositions.append(deepcopy(self.positions[atom.index]))
        for bond in self.topology.bonds():
            if bond[0] in newAtoms and bond[1] in newAtoms:
                newTopology.addBond(newAtoms[bond[0]], newAtoms[bond[1]])
        self.topology = newTopology
        self.positions = newPositions

    def addSolvent(self, forcefield, model='tip3p', boxSize=None, padding=None, positiveIon='Na+', negativeIon='Cl-', ionicStrength=0*molar):
        """Add solvent (both water and ions) to the model to fill a rectangular box.
        
        The algorithm works as follows:
        1. Water molecules are added to fill the box.
        2. Water molecules are removed if their distance to any solute atom is less than the sum of their van der Waals radii.
        3. If the solute is charged, enough positive or negative ions are added to neutralize it.  Each ion is added by
           randomly selecting a water molecule and replacing it with the ion.
        4. Ion pairs are added to give the requested total ionic strength.
        
        The box size can be specified in three ways.  First, you can explicitly give a box size to use.  Alternatively, you can
        give a padding distance.  The largest dimension of the solute (along the x, y, or z axis) is determined, and a cubic
        box of size (largest dimension)+2*padding is used.  Finally, if neither a box size nor a padding distance is specified,
        the existing Topology's unit cell dimensions are used.
        
        Parameters:
         - forcefield (ForceField) the ForceField to use for determining van der Waals radii and atomic charges
         - model (string='tip3p') the water model to use.  Supported values are 'tip3p', 'spce', 'tip4pew', and 'tip5p'.
         - boxSize (Vec3=None) the size of the box to fill with water
         - padding (distance=None) the padding distance to use
         - positiveIon (string='Na+') the type of positive ion to add.  Allowed values are 'Cs+', 'K+', 'Li+', 'Na+', and 'Rb+'
         - negativeIon (string='Cl-') the type of negative ion to add.  Allowed values are 'Cl-', 'Br-', 'F-', and 'I-'. Be aware
           that not all force fields support all ion types.
         - ionicString (concentration=0*molar) the total concentration of ions (both positive and negative) to add.  This
           does not include ions that are added to neutralize the system.
        """
        # Pick a unit cell size.
        
        if boxSize is not None:
            if is_quantity(boxSize):
                boxSize = boxSize.value_in_unit(nanometer)
            box = Vec3(boxSize[0], boxSize[1], boxSize[2])*nanometer
        elif padding is not None:
            maxSize = max(max((pos[i] for pos in self.positions))-min((pos[i] for pos in self.positions)) for i in range(3))
            box = (maxSize+2*padding)*Vec3(1, 1, 1)
        else:
            box = self.topology.getUnitCellDimensions()
            if box is None:
                raise ValueError('Neither the box size nor padding was specified, and the Topology does not define unit cell dimensions')
        box = box.value_in_unit(nanometer)
        invBox = Vec3(1.0/box[0], 1.0/box[1], 1.0/box[2])
        
        # Identify the ion types.
        
        posIonElements = {'Cs+':elem.cesium, 'K+':elem.potassium, 'Li+':elem.lithium, 'Na+':elem.sodium, 'Rb+':elem.rubidium}
        negIonElements = {'Cl-':elem.chlorine, 'Br-':elem.bromine, 'F-':elem.fluorine, 'I-':elem.iodine}
        if positiveIon not in posIonElements:
            raise ValueError('Illegal value for positive ion: %s' % positiveIon)
        if negativeIon not in negIonElements:
            raise ValueError('Illegal value for negative ion: %s' % negativeIon)
        positiveElement = posIonElements[positiveIon]
        negativeElement = negIonElements[negativeIon]
        
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
        else:
            raise ValueError('Unknown water model: %s' % model)
        pdb = PDBFile(os.path.join(os.path.dirname(__file__), 'data', model+'.pdb'))
        pdbTopology = pdb.getTopology()
        pdbPositions = pdb.getPositions().value_in_unit(nanometer)
        pdbResidues = list(pdbTopology.residues())
        pdbBoxSize = pdbTopology.getUnitCellDimensions().value_in_unit(nanometer)
        
        # Have the ForceField build a System for the solute from which we can determine van der Waals radii.
        
        system = forcefield.createSystem(self.topology)
        nonbonded = None
        for i in range(system.getNumForces()):
            if isinstance(system.getForce(i), NonbondedForce):
                nonbonded = system.getForce(i)
        if nonbonded is None:
            raise ValueError('The ForceField does not specify a NonbondedForce')
        cutoff = [nonbonded.getParticleParameters(i)[1].value_in_unit(nanometer)*vdwRadiusPerSigma+waterRadius for i in range(system.getNumParticles())]
        waterCutoff = waterRadius
        if len(cutoff) == 0:
            maxCutoff = waterCutoff
        else:
            maxCutoff = max(waterCutoff, max(cutoff))
        
        # Copy the solute over.

        newTopology = Topology()
        newTopology.setUnitCellDimensions(box)
        newAtoms = {}
        newPositions = []*nanometer
        for chain in self.topology.chains():
            newChain = newTopology.addChain()
            for residue in chain.residues():
                newResidue = newTopology.addResidue(residue.name, newChain)
                for atom in residue.atoms():
                    newAtom = newTopology.addAtom(atom.name, atom.element, newResidue)
                    newAtoms[atom] = newAtom
                    newPositions.append(deepcopy(self.positions[atom.index]))
        for bond in self.topology.bonds():
            newTopology.addBond(newAtoms[bond[0]], newAtoms[bond[1]])
        
        # Sort the solute atoms into cells for fast lookup.
        
        if len(self.positions) == 0:
            positions = []
        else:
            positions = self.positions.value_in_unit(nanometer)
        cells = {}
        numCells = tuple((max(1, int(floor(box[i]/maxCutoff))) for i in range(3)))
        cellSize = tuple((box[i]/numCells[i] for i in range(3)))
        for i in range(len(positions)):
            cell = tuple((int(floor(positions[i][j]/cellSize[j]))%numCells[j] for j in range(3)))
            if cell in cells:
                cells[cell].append(i)
            else:
                cells[cell] = [i]
        
        # Create a generator that loops over atoms close to a position.
        
        def neighbors(pos):
            centralCell = tuple((int(floor(pos[i]/cellSize[i])) for i in range(3)))
            offsets = (-1, 0, 1)
            for i in offsets:
                for j in offsets:
                    for k in offsets:
                        cell = ((centralCell[0]+i+numCells[0])%numCells[0], (centralCell[1]+j+numCells[1])%numCells[1], (centralCell[2]+k+numCells[2])%numCells[2])
                        if cell in cells:
                            for atom in cells[cell]:
                                yield atom
        
        # Define a function to compute the distance between two points, taking periodic boundary conditions into account.
        
        def periodicDistance(pos1, pos2):
            delta = pos1-pos2
            delta = [delta[i]-floor(delta[i]*invBox[i]+0.5)*box[i] for i in range(3)]
            return norm(delta)
        
        # Find the list of water molecules to add.
        
        newChain = newTopology.addChain()
        if len(positions) == 0:
            center = Vec3(0, 0, 0)
        else:
            center = [(max((pos[i] for pos in positions))+min((pos[i] for pos in positions)))/2 for i in range(3)]
            center = Vec3(center[0], center[1], center[2])
        numBoxes = [int(ceil(box[i]/pdbBoxSize[i])) for i in range(3)]
        addedWaters = []
        for boxx in range(numBoxes[0]):
            for boxy in range(numBoxes[1]):
                for boxz in range(numBoxes[2]):
                    offset = Vec3(boxx*pdbBoxSize[0], boxy*pdbBoxSize[1], boxz*pdbBoxSize[2])
                    for residue in pdbResidues:
                        oxygen = [atom for atom in residue.atoms() if atom.element == elem.oxygen][0]
                        atomPos = pdbPositions[oxygen.index]+offset
                        if not any((atomPos[i] > box[i] for i in range(3))):
                            # This molecule is inside the box, so see how close to it is to the solute.
                            
                            atomPos += center-box/2
                            for i in neighbors(atomPos):
                                if periodicDistance(atomPos, positions[i]) < cutoff[i]:
                                    break
                            else:
                                # Record this water molecule as one to add.
                            
                                addedWaters.append((residue.index, atomPos))
        
        # There could be clashes between water molecules at the box edges.  Find ones to remove.
        
        upperCutoff = center+box/2-Vec3(waterCutoff, waterCutoff, waterCutoff)
        lowerCutoff = center-box/2+Vec3(waterCutoff, waterCutoff, waterCutoff)
        lowerSkinPositions = [pos for index, pos in addedWaters if pos[0] < lowerCutoff[0] or pos[1] < lowerCutoff[1] or pos[2] < lowerCutoff[2]]
        filteredWaters = []
        cells = {}
        for i in range(len(lowerSkinPositions)):
            cell = tuple((int(floor(lowerSkinPositions[i][j]/cellSize[j]))%numCells[j] for j in range(3)))
            if cell in cells:
                cells[cell].append(i)
            else:
                cells[cell] = [i]
        for entry in addedWaters:
            pos = entry[1]
            if pos[0] < upperCutoff[0] and pos[1] < upperCutoff[1] and pos[2] < upperCutoff[2]:
                filteredWaters.append(entry)
            else:
                if not any((periodicDistance(lowerSkinPositions[i], pos) < waterCutoff and norm(lowerSkinPositions[i]-pos) > waterCutoff for i in neighbors(pos))):
                    filteredWaters.append(entry)
        addedWaters = filteredWaters
        
        # Add ions to neutralize the system.
        
        totalCharge = int(floor(0.5+sum((nonbonded.getParticleParameters(i)[0].value_in_unit(elementary_charge) for i in range(system.getNumParticles())))))
        if abs(totalCharge) > len(addedWaters):
            raise Exception('Cannot neutralize the system because the charge is greater than the number of available positions for ions')
        def addIon(element):
            # Replace a water by an ion.
            index = random.randint(0, len(addedWaters)-1)
            newResidue = newTopology.addResidue(element.symbol.upper(), newChain)
            newTopology.addAtom(element.symbol, element, newResidue)
            newPositions.append(addedWaters[index][1]*nanometer)
            del addedWaters[index]
        for i in range(abs(totalCharge)):
            addIon(positiveElement if totalCharge < 0 else negativeElement)
        
        # Add ions based on the desired ionic strength.
        
        numIons = len(addedWaters)*ionicStrength/(55.4*molar) # Pure water is about 55.4 molar (depending on temperature)
        numPairs = int(floor(numIons/2+0.5))
        for i in range(numPairs):
            addIon(positiveElement)
        for i in range(numPairs):
            addIon(negativeElement)
        
        # Add the water molecules.
        
        for index, pos in addedWaters:
            newResidue = newTopology.addResidue(residue.name, newChain)
            residue = pdbResidues[index]
            oxygen = [atom for atom in residue.atoms() if atom.element == elem.oxygen][0]
            oPos = pdbPositions[oxygen.index]
            molAtoms = []
            for atom in residue.atoms():
                molAtoms.append(newTopology.addAtom(atom.name, atom.element, newResidue))
                newPositions.append((pos+pdbPositions[atom.index]-oPos)*nanometer)
            for atom1 in molAtoms:
                if atom1.element == elem.oxygen:
                    for atom2 in molAtoms:
                        if atom2.element == elem.hydrogen:
                            newTopology.addBond(atom1, atom2)
        newTopology.setUnitCellDimensions(deepcopy(box)*nanometer)
        self.topology = newTopology
        self.positions = newPositions
    
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
    
    def addHydrogens(self, forcefield, pH=7.0, variants=None):
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
            
        Lysine:
            LYN: Neutral form with two hydrogens on the zeta nitrogen
            LYS: Positively charged form with three hydrogens on the zeta nitrogen
        
        The variant to use for each residue is determined by the following rules:
        
        1. The most common variant at the specified pH is selected.
        2. Any Cysteine that participates in a disulfide bond uses the CYX variant regardless of pH.
        3. For a neutral Histidine residue, the HID or HIE variant is selected based on which one forms a better hydrogen bond.
        
        You can override these rules by explicitly specifying a variant for any residue.  Also keep in mind that this
        function will only add hydrogens.  It will never remove ones that are already present in the model, regardless
        of the specified pH.
        
        Parameters:
         - forcefield (ForceField) the ForceField to use for determining the positions of hydrogens
         - pH (float=7.0) the pH based on which to select variants
         - variants (list=None) an optional list of variants to use.  If this is specified, its length must equal the number
           of residues in the model.  variants[i] is the name of the variant to use for residue i (indexed starting at 0).
           If an element is None, the standard rules will be followed to select a variant for that residue.
        """
        # Check the list of variants.
        
        residues = list(self.topology.residues())
        if variants is not None:
            if len(variants) != len(residues):
                raise ValueError("The length of the variants list must equal the number of residues")
        else:
            variants = [None]*len(residues)
        
        # Load the residue specifications.
        
        if Modeller._residueHydrogens is None:
            Modeller._residueHydrogens = {}
            tree = etree.parse(os.path.join(os.path.dirname(__file__), 'data', 'hydrogens.xml'))
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
        newTopology.setUnitCellDimensions(deepcopy(self.topology.getUnitCellDimensions()))
        newAtoms = {}
        newPositions = []*nanometer
        newIndices = []
        acceptors = [atom for atom in self.topology.atoms() if atom.element in (elem.oxygen, elem.nitrogen)]
        for chain in self.topology.chains():
            newChain = newTopology.addChain()
            for residue in chain.residues():
                newResidue = newTopology.addResidue(residue.name, newChain)
                isNTerminal = (residue == chain._residues[0])
                isCTerminal = (residue == chain._residues[-1])
                if residue.name in Modeller._residueHydrogens:
                    # Add hydrogens.  First select which variant to use.
                    
                    spec = Modeller._residueHydrogens[residue.name]
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
                                for acceptor in acceptors:
                                    if acceptor.residue != residue:
                                        acceptorPos = self.positions[acceptor.index]
                                        if isHbond(nd1Pos, hd1Pos, acceptorPos):
                                            nd1IsBonded = True
                                            break
                                        if isHbond(ne2Pos, he2Pos, acceptorPos):
                                            ne2IsBonded = True
                                if ne2IsBonded and not nd1IsBonded:
                                    variant = 'HIE'
                                else:
                                    variant = 'HID'
                        elif residue.name == 'HIS':
                            variant = 'HIP'
                    if variant is not None and variant not in spec.variants:
                        raise ValueError('Illegal variant for %s residue: %s' % (residue.name, variant))
                    
                    # Make a list of hydrogens that should be present in the residue.
                    
                    parents = [atom for atom in residue.atoms() if atom.element != elem.hydrogen]
                    parentNames = [atom.name for atom in parents]
                    hydrogens = [h for h in spec.hydrogens if (variant is None and pH <= h.maxph) or (h.variants is None and pH <= h.maxph) or (h.variants is not None and variant in h.variants)]
                    hydrogens = [h for h in hydrogens if h.terminal is None or (isNTerminal and h.terminal == 'N') or (isCTerminal and h.terminal == 'C')]
                    hydrogens = [h for h in hydrogens if h.parent in parentNames]
                    
                    # Loop over atoms in the residue, adding them to the new topology along with required hydrogens.
                    
                    for parent in residue.atoms():
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
                                    if len(expected) > 1:
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
                newTopology.addBond(newAtoms[bond[0]], newAtoms[bond[1]])
        
        # The hydrogens were added at random positions.  Now use the ForceField to fix them up.
        
        system = forcefield.createSystem(newTopology, rigidWater=False)
        atoms = list(newTopology.atoms())
        for i in range(system.getNumParticles()):
            if atoms[i].element != elem.hydrogen:
                # This is a heavy atom, so make it immobile.
                system.setParticleMass(i, 0)
        context = Context(system, VerletIntegrator(0.0))
        context.setPositions(newPositions)
        LocalEnergyMinimizer.minimize(context)
        self.topology = newTopology
        self.positions = context.getState(getPositions=True).getPositions()
