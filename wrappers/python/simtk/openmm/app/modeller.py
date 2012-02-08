"""
modeller.py: Provides tools for editing molecular models
"""
__author__ = "Peter Eastman"
__version__ = "1.0"

from simtk.openmm.app import Topology, PDBFile
from simtk.openmm.vec3 import Vec3
from simtk.openmm import NonbondedForce
from simtk.unit import nanometer, molar, elementary_charge, amu, gram, liter, sqrt, is_quantity
import element as elem
import os
import random
from copy import deepcopy
from math import ceil, floor

class Modeller(object):
    """Modeller provides tools for editing molecular models, such as adding water or missing hydrogens.
    
    To use it, create a Modeller object, specifying the initial Topology and atom positions.  You can
    then call various methods to change the model in different ways.  Each time you do, a new Topology
    and list of coordinates is created to represent the changed model.  Finally, call getTopology()
    and getPositions() to get the results.
    """
        
    def __init__(self, topology, positions):
        """Create a new Modeller object
        
        Parameters:
         - topology (Topology) the initial Topology of the model
         - positions (list) the initial atomic positions
        """
        self.topology = topology
        if not is_quantity(positions):
            positions = positions*nanometers
        self.positions = positions
        
    def getTopology(self):
        """Get the Topology of the model."""
        return self.topology
        
    def getPositions(self):
        """Get the atomic positions."""
        return self.positions

    def deleteWater(self):
        """Delete all water molecules from the model."""
        newTopology = Topology()
        newTopology.setUnitCellDimensions(deepcopy(self.topology.getUnitCellDimensions()))
        newAtoms = {}
        newPositions = []*nanometer
        for chain in self.topology.chains():
            newChain = newTopology.addChain()
            for residue in chain.residues():
                if residue.name != "HOH":
                    newResidue = newTopology.addResidue(residue.name, newChain)
                    for atom in residue.atoms():
                        newAtom = newTopology.addAtom(atom.name, atom.element, newResidue)
                        newAtoms[atom] = newAtom
                        newPositions.append(deepcopy(self.positions[atom.index]))
        for bond in self.topology.bonds():
            if bond[0] in newAtoms and bond[1] in newAtoms:
                newTopology.addBond(newAtoms[bond[0]], newAtoms[bond[1]])
        self.topology = newTopology
        self.positions = newPositions
    
    def convertWater(self, model='tip3p'):
        """Convert all water molecules to a different water model.
        
        Parameters:
         - model (string='tip3p') the water model to convert to.  Supported values are 'tip3p', 'tip4pew', and 'tip5p'.
        """
        if model == 'tip3p':
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
         - model (string='tip3p') the water model to use.  Supported values are 'tip3p', 'tip4pew', and 'tip5p'.
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
            box = boxSize
        elif padding is not None:
            maxSize = max(max((pos[i] for pos in self.positions))-min((pos[i] for pos in self.positions)) for i in range(3))
            box = (maxSize+2*padding)*Vec3(1, 1, 1)
        else:
            box = topology.getUnitCellDimensions()
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
        waterCutoff = 2*waterRadius
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
            if bond[0] in newAtoms and bond[1] in newAtoms:
                newTopology.addBond(newAtoms[bond[0]], newAtoms[bond[1]])
        
        # Sort the solute atoms into cells for fast lookup.
        
        positions = self.positions.value_in_unit(nanometer)
        cells = {}
        numCells = tuple((int(floor(box[i]/maxCutoff)) for i in range(3)))
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
            return sqrt(delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2])
        
        # Find the list of water molecules to add.
        
        newChain = newTopology.addChain()
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
                if not any((periodicDistance(lowerSkinPositions[i], pos) < waterCutoff for i in neighbors(pos))):
                    filteredWaters.append(entry)
        addedWaters = filteredWaters
        
        # Add ions to neutralize the system.
        
        totalCharge = int(sum((nonbonded.getParticleParameters(i)[0].value_in_unit(elementary_charge) for i in range(system.getNumParticles()))))
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
