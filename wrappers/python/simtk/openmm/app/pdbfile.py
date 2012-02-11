"""
pdbfile.py: Used for loading PDB files.
"""
__author__ = "Peter Eastman"
__version__ = "1.0"

import os
import sys
import xml.etree.ElementTree as etree
from copy import copy
from simtk.openmm import Vec3
from simtk.openmm.app.internal.pdbstructure import PdbStructure
from simtk.openmm.app import Topology
from simtk.unit import nanometers, angstroms, is_quantity
import element as elem
try:
    import numpy
except:
    pass

class PDBFile(object):
    """PDBFile parses a Protein Data Bank (PDB) file and constructs a Topology and a set of atom positions from it.
    
    This class also provides methods for creating PDB files.  To write a file containing a single model, call
    writeFile().  You also can create files that contain multiple models.  To do this, first call writeHeader(),
    then writeModel() once for each model in the file, and finally writeFooter() to complete the file."""
    
    _residueNameReplacements = {}
    _atomNameReplacements = {}
    
    def __init__(self, file):
        """Load a PDB file.
        
        The atom positions and Topology can be retrieved by calling getPositions() and getTopology().
        
        Parameters:
         - file (string) the name of the file to load
        """
        top = Topology()
        coords = [];
        self.topology = top
        
        # Load the PDB file
        
        pdb = PdbStructure(open(file))
        PDBFile._loadNameReplacementTables()

        # Build the topology

        for chain in pdb.iter_chains():
            c = top.addChain()
            for residue in chain.iter_residues():
                resName = residue.get_name()
                if resName in PDBFile._residueNameReplacements:
                    resName = PDBFile._residueNameReplacements[resName]
                r = top.addResidue(resName, c)
                if resName in PDBFile._atomNameReplacements:
                    atomReplacements = PDBFile._atomNameReplacements[resName]
                else:
                    atomReplacements = {}
                for atom in residue.atoms:
                    atomName = atom.get_name()
                    if atomName in atomReplacements:
                        atomName = atomReplacements[atomName]
                    atomName = atomName.strip()
                    element = None

                    # Try to guess the element.
                    
                    upper = atomName.upper()
                    if upper.startswith('CL'):
                        element = elem.chlorine
                    elif upper.startswith('NA'):
                        element = elem.sodium
                    elif upper.startswith('MG'):
                        element = elem.magnesium
                    elif upper.startswith('BE'):
                        element = elem.beryllium
                    elif upper.startswith('LI'):
                        element = elem.lithium
                    elif upper.startswith('K'):
                        element = elem.potassium
                    elif( len( residue ) == 1 and upper.startswith('CA') ):
                        element = elem.calcium
                    else:
                        try:
                            element = elem.get_by_symbol(atomName[0])
                        except KeyError:
                            pass
                    top.addAtom(atomName, element, r)
                    pos = atom.get_position().value_in_unit(nanometers)
                    coords.append(Vec3(pos[0], pos[1], pos[2]))
        self.positions = coords*nanometers
        self.topology.setUnitCellDimensions(pdb.get_unit_cell_dimensions())
        self.topology.createStandardBonds()
        self.topology.createDisulfideBonds(self.positions)
        self._numpyPositions = None

    def getTopology(self):
        """Get the Topology of the model."""
        return self.topology
        
    def getPositions(self, asNumpy=False):
        """Get the atomic positions.
        
        Parameters:
         - asNumpy (boolean=False) if true, the values are returned as a numpy array instead of a list of Vec3s
         """
        if asNumpy:
            if self._numpyPositions is None:
                self._numpyPositions = numpy.array(self.positions.value_in_unit(nanometers))*nanometers
            return self._numpyPositions
        return self.positions

    @staticmethod
    def _loadNameReplacementTables():
        """Load the list of atom and residue name replacements."""
        if len(PDBFile._residueNameReplacements) == 0:
            tree = etree.parse(os.path.join(os.path.dirname(__file__), 'data', 'pdbNames.xml'))
            allResidues = {}
            proteinResidues = {}
            nucleicAcidResidues = {}
            for residue in tree.getroot().findall('Residue'):
                name = residue.attrib['name']
                if name == 'All':
                    PDBFile._parseResidueAtoms(residue, allResidues)
                elif name == 'Protein':
                    PDBFile._parseResidueAtoms(residue, proteinResidues)
                elif name == 'Nucleic':
                    PDBFile._parseResidueAtoms(residue, nucleicAcidResidues)
            for atom in allResidues:
                proteinResidues[atom] = allResidues[atom]
                nucleicAcidResidues[atom] = allResidues[atom]
            for residue in tree.getroot().findall('Residue'):
                name = residue.attrib['name']
                for id in residue.attrib:
                    if id == 'name' or id.startswith('alt'):
                        PDBFile._residueNameReplacements[residue.attrib[id]] = name
                if 'type' not in residue.attrib:
                    atoms = copy(allResidues)
                elif residue.attrib['type'] == 'Protein':
                    atoms = copy(proteinResidues)
                elif residue.attrib['type'] == 'Nucleic':
                    atoms = copy(nucleicAcidResidues)
                else:
                    atoms = copy(allResidues)
                PDBFile._parseResidueAtoms(residue, atoms)
                PDBFile._atomNameReplacements[name] = atoms

    @staticmethod
    def _parseResidueAtoms(residue, map):
        for atom in residue.findall('Atom'):
            name = atom.attrib['name']
            for id in atom.attrib:
                map[atom.attrib[id]] = name
    
    @staticmethod
    def writeFile(topology, positions, file=sys.stdout, modelIndex=None):
        """Write a PDB file containing a single model.
                
        Parameters:
         - topology (Topology) The Topology defining the model to write
         - positions (list) The list of atomic positions to write
         - file (file=stdout) A file to write to
        """
        PDBFile.writeHeader(topology, file)
        PDBFile.writeModel(topology, positions, file)
        PDBFile.writeFooter(topology, file)

    @staticmethod
    def writeHeader(topology, file=sys.stdout):
        """Write out the header for a PDB file.
        
        Parameters:
         - topology (Topology) The Topology defining the molecular system being written
         - file (file=stdout) A file to write the file to
        """
        boxSize = topology.getUnitCellDimensions()
        if boxSize is not None:
            size = boxSize.value_in_unit(angstroms)
            print >>file, "CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1 " % size
    
    @staticmethod
    def writeModel(topology, positions, file=sys.stdout, modelIndex=None):
        """Write out a model to a PDB file.
                
        Parameters:
         - topology (Topology) The Topology defining the model to write
         - positions (list) The list of atomic positions to write
         - file (file=stdout) A file to write the model to
         - modelIndex (int=None) If not None, the model will be surrounded by MODEL/ENDMDL records with this index
        """
        if len(list(topology.atoms())) != len(positions):
            raise ValueError('The number of positions must match the number of atoms') 
        if is_quantity(positions):
            positions = positions.value_in_unit(angstroms)
        atomIndex = 1
        posIndex = 0
        if modelIndex is not None:
            print >>file, "MODEL     %4d" % modelIndex
        for (chainIndex, chain) in enumerate(topology.chains()):
            chainName = chr(ord('A')+chainIndex)
            residues = list(chain.residues())
            for (resIndex, res) in enumerate(residues):
                if len(res.name) > 3:
                    resName = res.name[:3]
                else:
                    resName = res.name
                for atom in res.atoms():
                    if len(atom.name) < 4:
                        atomName = ' '+atom.name
                    elif len(atom.name) > 4:
                        atomName = atom.name[:4]
                    else:
                        atomName = atom.name
                    coords = positions[posIndex]
                    print >>file, "ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00" % (atomIndex, atomName, resName, chainName, resIndex+1, coords[0], coords[1], coords[2])
                    posIndex += 1
                    atomIndex += 1
                if resIndex == len(residues)-1:
                    print >>file, "TER   %5d      %3s %s%4d" % (atomIndex, resName, chainName, resIndex+1)
                    atomIndex += 1
        if modelIndex is not None:
            print >>file, "ENDMDL"

    @staticmethod
    def writeFooter(topology, file=sys.stdout):
        """Write out the footer for a PDB file.
        
        Parameters:
         - topology (Topology) The Topology defining the molecular system being written
         - file (file=stdout) A file to write the file to
        """
        print >>file, "END"

