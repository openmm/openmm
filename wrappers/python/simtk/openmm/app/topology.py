"""
topology.py: Used for storing topological information about a system.

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
import xml.etree.ElementTree as etree
from simtk.unit import nanometers, sqrt

class Topology(object):
    """Topology stores the topological information about a system.
    
    The structure of a Topology object is similar to that of a PDB file.  It consists of a set of Chains
    (often but not always corresponding to polymer chains).  Each Chain contains a set of Residues,
    and each Residue contains a set of Atoms.  In addition, the Topology stores a list of which atom
    pairs are bonded to each other, and the dimensions of the crystallographic unit cell.
    
    Atom and residue names should follow the PDB 3.0 nomenclature for all molecules for which one exists.
    """
    
    _standardBonds = {}
    
    def __init__(self):
        """Create a new Topology object"""
        self._chains = []
        self._numResidues = 0
        self._numAtoms = 0
        self._bonds = []
        self._unitCellDimensions = None
        
    def addChain(self):
        """Create a new Chain and add it to the Topology.
        
        Returns: the newly created Chain
        """
        chain = Chain(len(self._chains), self)
        self._chains.append(chain)
        return chain
    
    def addResidue(self, name, chain):
        """Create a new Residue and add it to the Topology.
        
        Parameters:
         - name (string) The name of the residue to add
         - chain (Chain) The Chain to add it to
        Returns: the newly created Residue
        """
        residue = Residue(name, self._numResidues, chain)
        self._numResidues += 1
        chain._residues.append(residue)
        return residue
    
    def addAtom(self, name, element, residue):
        """Create a new Atom and add it to the Topology.
        
        Parameters:
         - name (string) The name of the atom to add
         - element (Element) The element of the atom to add
         - residue (Residue) The Residue to add it to
        Returns: the newly created Atom
        """
        atom = Atom(name, element, self._numAtoms, residue)
        self._numAtoms += 1
        residue._atoms.append(atom)
        return atom
    
    def addBond(self, atom1, atom2):
        """Create a new bond and add it to the Topology.
        
        Parameters:
         - atom1 (Atom) The first Atom connected by the bond
         - atom2 (Atom) The second Atom connected by the bond
        """
        self._bonds.append((atom1, atom2))
    
    def chains(self):
        """Iterate over all Chains in the Topology."""
        return iter(self._chains)
    
    def residues(self):
        """Iterate over all Residues in the Topology."""
        for chain in self._chains:
            for residue in chain._residues:
                yield residue
    
    def atoms(self):
        """Iterate over all Atoms in the Topology."""
        for chain in self._chains:
            for residue in chain._residues:
                for atom in residue._atoms:
                    yield atom
    
    def bonds(self):
        """Iterate over all bonds (each represented as a tuple of two Atoms) in the Topology."""
        return iter(self._bonds)
        
    def getUnitCellDimensions(self):
        """Get the dimensions of the crystallographic unit cell.
        
        The return value may be None if this Topology does not represent a periodic structure.
        """
        return self._unitCellDimensions
        
    def setUnitCellDimensions(self, dimensions):
        """Set the dimensions of the crystallographic unit cell."""
        self._unitCellDimensions = dimensions
    
    def createStandardBonds(self):
        """Create bonds based on the atom and residue names for all standard residue types."""
        if len(Topology._standardBonds) == 0:
            # Load the standard bond defitions.
            
            tree = etree.parse(os.path.join(os.path.dirname(__file__), 'data', 'residues.xml'))
            for residue in tree.getroot().findall('Residue'):
                bonds = []
                Topology._standardBonds[residue.attrib['name']] = bonds
                for bond in residue.findall('Bond'):
                    bonds.append((bond.attrib['from'], bond.attrib['to']))
        for chain in self._chains:
            # First build a map of atom names to atoms.
            
            atomMaps = []
            for residue in chain._residues:
                atomMap = {}
                atomMaps.append(atomMap)
                for atom in residue._atoms:
                    atomMap[atom.name] = atom
            
            # Loop over residues and construct bonds.
            
            for i in range(len(chain._residues)):
                name = chain._residues[i].name
                if name in Topology._standardBonds:
                    for bond in Topology._standardBonds[name]:
                        if bond[0].startswith('-') and i > 0:
                            fromResidue = i-1
                            fromAtom = bond[0][1:]
                        elif bond[0].startswith('+') and i <len(chain._residues):
                            fromResidue = i+1
                            fromAtom = bond[0][1:]
                        else:
                            fromResidue = i
                            fromAtom = bond[0]
                        if bond[1].startswith('-') and i > 0:
                            toResidue = i-1
                            toAtom = bond[1][1:]
                        elif bond[1].startswith('+') and i <len(chain._residues):
                            toResidue = i+1
                            toAtom = bond[1][1:]
                        else:
                            toResidue = i
                            toAtom = bond[1]
                        if fromAtom in atomMaps[fromResidue] and toAtom in atomMaps[toResidue]:
                            self.addBond(atomMaps[fromResidue][fromAtom], atomMaps[toResidue][toAtom])
                            
    def createDisulfideBonds(self, positions):
        """Identify disulfide bonds based on proximity and add them to the Topology.
        
        Parameters:
         - positions (list) The list of atomic positions based on which to identify bonded atoms
        """
        def isCyx(res):
            names = [atom.name for atom in res._atoms]
            return 'SG' in names and 'HG' not in names
        
        cyx = [res for res in self.residues() if res.name == 'CYS' and isCyx(res)]
        atomNames = [[atom.name for atom in res._atoms] for res in cyx]
        for i in range(len(cyx)):
            sg1 = cyx[i]._atoms[atomNames[i].index('SG')]
            pos1 = positions[sg1.index]
            for j in range(i):
                sg2 = cyx[j]._atoms[atomNames[j].index('SG')]
                pos2 = positions[sg2.index]
                delta = [x-y for (x,y) in zip(pos1, pos2)]
                distance = sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2])
                if distance < 0.3*nanometers:
                    self.addBond(sg1, sg2)

class Chain(object):
    """A Chain object represents a chain within a Topology."""
    def __init__(self, index, topology):
        """Construct a new Chain.  You should call addChain() on the Topology instead of calling this directly."""
        ## The index of the Chain within its Topology
        self.index = index
        ## The Topology this Chain belongs to
        self.topology = topology
        self._residues = []
    
    def residues(self):
        """Iterate over all Residues in the Chain."""
        return iter(self._residues)

    def atoms(self):
        """Iterate over all Atoms in the Chain."""
        for residue in self._residues:
            for atom in residue._atoms:
                yield atom

class Residue(object):
    """A Residue object represents a residue within a Topology."""
    def __init__(self, name, index, chain):
        """Construct a new Residue.  You should call addResidue() on the Topology instead of calling this directly."""
        ## The name of the Residue
        self.name = name
        ## The index of the Residue within its Topology
        self.index = index
        ## The Chain this Residue belongs to
        self.chain = chain
        self._atoms = []
    
    def atoms(self):
        """Iterate over all Atoms in the Residue."""
        return iter(self._atoms)

class Atom(object):
    """An Atom object represents a residue within a Topology."""
    
    def __init__(self, name, element, index, residue):
        """Construct a new Atom.  You should call addAtom() on the Topology instead of calling this directly."""
        ## The name of the Atom
        self.name = name
        ## That Atom's element
        self.element = element
        ## The index of the Atom within its Topology
        self.index = index
        ## The Residue this Atom belongs to
        self.residue = residue

