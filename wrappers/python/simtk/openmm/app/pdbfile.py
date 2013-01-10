"""
pdbfile.py: Used for loading PDB files.

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
__author__ = "Peter Eastman"
__version__ = "1.0"

import os
import sys
import math
import xml.etree.ElementTree as etree
from copy import copy
from simtk.openmm import Vec3
from simtk.openmm.app.internal.pdbstructure import PdbStructure
from simtk.openmm.app import Topology
from simtk.unit import nanometers, angstroms, is_quantity, norm, Quantity
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
        ## The Topology read from the PDB file
        self.topology = top
        
        # Load the PDB file
        
        pdb = PdbStructure(open(file), load_all_models=True)
        PDBFile._loadNameReplacementTables()

        # Build the topology

        atomByNumber = {}
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
                    element = atom.element
                    if element is None:
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
                    newAtom = top.addAtom(atomName, element, r)
                    atomByNumber[atom.serial_number] = newAtom
        self._positions = []
        for model in pdb.iter_models(True):
            coords = []
            for atom in model.iter_atoms():
                pos = atom.get_position().value_in_unit(nanometers)
                coords.append(Vec3(pos[0], pos[1], pos[2]))
            self._positions.append(coords*nanometers)
        ## The atom positions read from the PDB file.  If the file contains multiple frames, these are the positions in the first frame.
        self.positions = self._positions[0]
        self.topology.setUnitCellDimensions(pdb.get_unit_cell_dimensions())
        self.topology.createStandardBonds()
        self.topology.createDisulfideBonds(self.positions)
        self._numpyPositions = None
        
        # Add bonds based on CONECT records.
        
        connectBonds = []
        for connect in pdb.models[0].connects:
            i = connect[0]
            for j in connect[1:]:
                connectBonds.append((atomByNumber[i], atomByNumber[j]))
        if len(connectBonds) > 0:
            # Only add bonds that don't already exist.
            existingBonds = set(top.bonds())
            for bond in connectBonds:
                if bond not in existingBonds and (bond[1], bond[0]) not in existingBonds:
                    top.addBond(bond[0], bond[1])
                    existingBonds.add(bond)

    def getTopology(self):
        """Get the Topology of the model."""
        return self.topology
        
    def getNumFrames(self):
        """Get the number of frames stored in the file."""
        return len(self._positions)
    
    def getPositions(self, asNumpy=False, frame=0):
        """Get the atomic positions.
        
        Parameters:
         - asNumpy (boolean=False) if true, the values are returned as a numpy array instead of a list of Vec3s
         - frame (int=0) the index of the frame for which to get positions
         """
        if asNumpy:
            if self._numpyPositions is None:
                self._numpyPositions = [None]*len(self._positions)
            if self._numpyPositions[frame] is None:
                self._numpyPositions[frame] = Quantity(numpy.array(self._positions[frame].value_in_unit(nanometers)), nanometers)
            return self._numpyPositions[frame]
        return self._positions[frame]

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
        if any(math.isnan(norm(pos)) for pos in positions):
            raise ValueError('Particle position is NaN')
        if any(math.isinf(norm(pos)) for pos in positions):
            raise ValueError('Particle position is infinite')
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
                    if len(atom.name) < 4 and atom.name[:1].isalpha() and (atom.element is None or len(atom.element.symbol) < 2):
                        atomName = ' '+atom.name
                    elif len(atom.name) > 4:
                        atomName = atom.name[:4]
                    else:
                        atomName = atom.name
                    coords = positions[posIndex]
                    print >>file, "ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00" % (atomIndex%100000, atomName, resName, chainName, (resIndex+1)%10000, coords[0], coords[1], coords[2])
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

