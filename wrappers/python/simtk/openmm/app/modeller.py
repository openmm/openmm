"""
modeller.py: Provides tools for editing molecular models
"""
__author__ = "Peter Eastman"
__version__ = "1.0"

from simtk.openmm.app import Topology
from simtk.openmm.vec3 import Vec3
import simtk.unit as unit
import element as elem
import copy

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
        if not unit.is_quantity(positions):
            positions = positions*unit.nanometers
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
        newTopology.setUnitCellDimensions(copy.deepcopy(self.topology.getUnitCellDimensions()))
        newAtoms = {}
        newPositions = []*unit.nanometer
        for chain in self.topology.chains():
            newChain = newTopology.addChain()
            for residue in chain.residues():
                if residue.name != "HOH":
                    newResidue = newTopology.addResidue(residue.name, newChain)
                    for atom in residue.atoms():
                        newAtom = newTopology.addAtom(atom.name, atom.element, newResidue)
                        newAtoms[atom] = newAtom
                        newPositions.append(copy.deepcopy(self.positions[atom.index]))
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
        newTopology.setUnitCellDimensions(copy.deepcopy(self.topology.getUnitCellDimensions()))
        newAtoms = {}
        newPositions = []*unit.nanometer
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
                    po = copy.deepcopy(self.positions[oatom[0].index])
                    ph1 = copy.deepcopy(self.positions[hatoms[0].index])
                    ph2 = copy.deepcopy(self.positions[hatoms[1].index])
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
                        v1 = (ph1-po).value_in_unit(unit.nanometer)
                        v2 = (ph2-po).value_in_unit(unit.nanometer)
                        cross = Vec3(v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0])
                        newPositions.append(po - (0.34490826*v1 - 0.34490826*v2 - 6.4437903*cross)*unit.nanometer)
                        newPositions.append(po - (0.34490826*v1 - 0.34490826*v2 + 6.4437903*cross)*unit.nanometer)
                else:
                    # Just copy the residue over.
                    for atom in residue.atoms():
                        newAtom = newTopology.addAtom(atom.name, atom.element, newResidue)
                        newAtoms[atom] = newAtom
                        newPositions.append(copy.deepcopy(self.positions[atom.index]))
        for bond in self.topology.bonds():
            if bond[0] in newAtoms and bond[1] in newAtoms:
                newTopology.addBond(newAtoms[bond[0]], newAtoms[bond[1]])
        self.topology = newTopology
        self.positions = newPositions
