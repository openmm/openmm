"""
topology.py: Used for storing topological information about a system.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012-2018 Stanford University and the Authors.
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
from __future__ import absolute_import
__author__ = "Peter Eastman"
__version__ = "1.0"

from collections import namedtuple
import os
import xml.etree.ElementTree as etree
from simtk.openmm.vec3 import Vec3
from simtk.openmm.app.internal.singleton import Singleton
from simtk.unit import nanometers, sqrt, is_quantity
from copy import deepcopy

# Enumerated values for bond type

class Single(Singleton):
    def __repr__(self):
        return 'Single'
Single = Single()

class Double(Singleton):
    def __repr__(self):
        return 'Double'
Double = Double()

class Triple(Singleton):
    def __repr__(self):
        return 'Triple'
Triple = Triple()

class Aromatic(Singleton):
    def __repr__(self):
        return 'Aromatic'
Aromatic = Aromatic()

class Amide(Singleton):
    def __repr__(self):
        return 'Amide'
Amide = Amide()

class Topology(object):
    """Topology stores the topological information about a system.

    The structure of a Topology object is similar to that of a PDB file.  It consists of a set of Chains
    (often but not always corresponding to polymer chains).  Each Chain contains a set of Residues,
    and each Residue contains a set of Atoms.  In addition, the Topology stores a list of which atom
    pairs are bonded to each other, and the dimensions of the crystallographic unit cell.

    Atom and residue names should follow the PDB 3.0 nomenclature for all molecules for which one exists.
    """

    _standardBonds = {}
    _hasLoadedStandardBonds = False

    def __init__(self):
        """Create a new Topology object"""
        self._chains = []
        self._numResidues = 0
        self._numAtoms = 0
        self._bonds = []
        self._periodicBoxVectors = None

    def __repr__(self):
        nchains = len(self._chains)
        nres = self._numResidues
        natom = self._numAtoms
        nbond = len(self._bonds)
        return '<%s; %d chains, %d residues, %d atoms, %d bonds>' % (
                type(self).__name__, nchains, nres, natom, nbond)

    def getNumAtoms(self):
        """Return the number of atoms in the Topology.
        """
        return self._numAtoms

    def getNumResidues(self):
        """Return the number of residues in the Topology.
        """
        return self._numResidues

    def getNumChains(self):
        """Return the number of chains in the Topology.
        """
        return len(self._chains)

    def getNumBonds(self):
        """Return the number of bonds in the Topology.
        """
        return len(self._bonds)

    def addChain(self, id=None):
        """Create a new Chain and add it to the Topology.

        Parameters
        ----------
        id : string=None
            An optional identifier for the chain.  If this is omitted, an id is
            generated based on the chain index.

        Returns
        -------
        Chain
             the newly created Chain
        """
        if id is None:
            id = str(len(self._chains)+1)
        chain = Chain(len(self._chains), self, id)
        self._chains.append(chain)
        return chain

    def addResidue(self, name, chain, id=None, insertionCode=''):
        """Create a new Residue and add it to the Topology.

        Parameters
        ----------
        name : string
            The name of the residue to add
        chain : Chain
            The Chain to add it to
        id : string=None
            An optional identifier for the residue.  If this is omitted, an id
            is generated based on the residue index.
        insertionCode: string=''
            An optional insertion code for the residue.

        Returns
        -------
        Residue
             the newly created Residue
        """
        if len(chain._residues) > 0 and self._numResidues != chain._residues[-1].index+1:
            raise ValueError('All residues within a chain must be contiguous')
        if id is None:
            id = str(self._numResidues+1)
        residue = Residue(name, self._numResidues, chain, id, insertionCode)
        self._numResidues += 1
        chain._residues.append(residue)
        return residue

    def addAtom(self, name, element, residue, id=None):
        """Create a new Atom and add it to the Topology.

        Parameters
        ----------
        name : string
            The name of the atom to add
        element : Element
            The element of the atom to add
        residue : Residue
            The Residue to add it to
        id : string=None
            An optional identifier for the atom.  If this is omitted, an id is
            generated based on the atom index.

        Returns
        -------
        Atom
             the newly created Atom
        """
        if len(residue._atoms) > 0 and self._numAtoms != residue._atoms[-1].index+1:
            raise ValueError('All atoms within a residue must be contiguous')
        if id is None:
            id = str(self._numAtoms+1)
        atom = Atom(name, element, self._numAtoms, residue, id)
        self._numAtoms += 1
        residue._atoms.append(atom)
        return atom

    def addBond(self, atom1, atom2, type=None, order=None):
        """Create a new bond and add it to the Topology.

        Parameters
        ----------
        atom1 : Atom
            The first Atom connected by the bond
        atom2 : Atom
            The second Atom connected by the bond
        type : object=None
            The type of bond to add.  Allowed values are None, Single, Double, Triple,
            Aromatic, or Amide.
        order : int=None
            The bond order, or None if it is not specified
        """
        self._bonds.append(Bond(atom1, atom2, type, order))

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

    def getPeriodicBoxVectors(self):
        """Get the vectors defining the periodic box.

        The return value may be None if this Topology does not represent a periodic structure."""
        return self._periodicBoxVectors

    def setPeriodicBoxVectors(self, vectors):
        """Set the vectors defining the periodic box."""
        if vectors is not None:
            if not is_quantity(vectors[0][0]):
                vectors = vectors*nanometers
            if vectors[0][1] != 0*nanometers or vectors[0][2] != 0*nanometers:
                raise ValueError("First periodic box vector must be parallel to x.");
            if vectors[1][2] != 0*nanometers:
                raise ValueError("Second periodic box vector must be in the x-y plane.");
            if vectors[0][0] <= 0*nanometers or vectors[1][1] <= 0*nanometers or vectors[2][2] <= 0*nanometers or vectors[0][0] < 2*abs(vectors[1][0]) or vectors[0][0] < 2*abs(vectors[2][0]) or vectors[1][1] < 2*abs(vectors[2][1]):
                raise ValueError("Periodic box vectors must be in reduced form.");
        self._periodicBoxVectors = deepcopy(vectors)

    def getUnitCellDimensions(self):
        """Get the dimensions of the crystallographic unit cell.

        The return value may be None if this Topology does not represent a periodic structure.
        """
        if self._periodicBoxVectors is None:
            return None
        xsize = self._periodicBoxVectors[0][0].value_in_unit(nanometers)
        ysize = self._periodicBoxVectors[1][1].value_in_unit(nanometers)
        zsize = self._periodicBoxVectors[2][2].value_in_unit(nanometers)
        return Vec3(xsize, ysize, zsize)*nanometers

    def setUnitCellDimensions(self, dimensions):
        """Set the dimensions of the crystallographic unit cell.

        This method is an alternative to setPeriodicBoxVectors() for the case of a rectangular box.  It sets
        the box vectors to be orthogonal to each other and to have the specified lengths."""
        if dimensions is None:
            self._periodicBoxVectors = None
        else:
            if is_quantity(dimensions):
                dimensions = dimensions.value_in_unit(nanometers)
            self._periodicBoxVectors = (Vec3(dimensions[0], 0, 0), Vec3(0, dimensions[1], 0), Vec3(0, 0, dimensions[2]))*nanometers

    @staticmethod
    def loadBondDefinitions(file):
        """Load an XML file containing definitions of bonds that should be used by createStandardBonds().

        The built in residues.xml file containing definitions for standard amino acids and nucleotides is loaded automatically.
        This method can be used to load additional definitions for other residue types.  They will then be used in subsequent
        calls to createStandardBonds().  This is a static method, so it affects subsequent calls on all Topology objects.
        Also note that PDBFile calls createStandardBonds() automatically when a file is loaded, so the newly loaded definitions
        will be used for any PDB file loaded after this is called.
        """
        tree = etree.parse(file)
        for residue in tree.getroot().findall('Residue'):
            bonds = []
            Topology._standardBonds[residue.attrib['name']] = bonds
            for bond in residue.findall('Bond'):
                bonds.append((bond.attrib['from'], bond.attrib['to']))

    def createStandardBonds(self):
        """Create bonds based on the atom and residue names for all standard residue types.

        Definitions for standard amino acids and nucleotides are built in.  You can call loadBondDefinitions() to load
        additional definitions for other residue types.
        """
        if not Topology._hasLoadedStandardBonds:
            # Load the standard bond definitions.

            Topology.loadBondDefinitions(os.path.join(os.path.dirname(__file__), 'data', 'residues.xml'))
            Topology._hasLoadedStandardBonds = True
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
        """Identify disulfide bonds based on proximity and add them to the
        Topology.

        Parameters
        ----------
        positions : list
            The list of atomic positions based on which to identify bonded atoms
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
    def __init__(self, index, topology, id):
        """Construct a new Chain.  You should call addChain() on the Topology instead of calling this directly."""
        ## The index of the Chain within its Topology
        self.index = index
        ## The Topology this Chain belongs to
        self.topology = topology
        ## A user defined identifier for this Chain
        self.id = id
        self._residues = []

    def residues(self):
        """Iterate over all Residues in the Chain."""
        return iter(self._residues)

    def atoms(self):
        """Iterate over all Atoms in the Chain."""
        for residue in self._residues:
            for atom in residue._atoms:
                yield atom

    def __len__(self):
        return len(self._residues)

    def __repr__(self):
        return "<Chain %d>" % self.index

class Residue(object):
    """A Residue object represents a residue within a Topology."""
    def __init__(self, name, index, chain, id, insertionCode):
        """Construct a new Residue.  You should call addResidue() on the Topology instead of calling this directly."""
        ## The name of the Residue
        self.name = name
        ## The index of the Residue within its Topology
        self.index = index
        ## The Chain this Residue belongs to
        self.chain = chain
        ## A user defined identifier for this Residue
        self.id = id
        ## A user defined insertion code for this Residue
        self.insertionCode = insertionCode
        self._atoms = []

    def atoms(self):
        """Iterate over all Atoms in the Residue."""
        return iter(self._atoms)

    def bonds(self):
        """Iterate over all Bonds involving any atom in this residue."""
        return ( bond for bond in self.chain.topology.bonds() if ((bond[0] in self._atoms) or (bond[1] in self._atoms)) )

    def internal_bonds(self):
        """Iterate over all internal Bonds."""
        return ( bond for bond in self.chain.topology.bonds() if ((bond[0] in self._atoms) and (bond[1] in self._atoms)) )

    def external_bonds(self):
        """Iterate over all Bonds to external atoms."""
        return ( bond for bond in self.chain.topology.bonds() if ((bond[0] in self._atoms) != (bond[1] in self._atoms)) )

    def __len__(self):
        return len(self._atoms)

    def __repr__(self):
        return "<Residue %d (%s) of chain %d>" % (self.index, self.name, self.chain.index)

class Atom(object):
    """An Atom object represents an atom within a Topology."""

    def __init__(self, name, element, index, residue, id):
        """Construct a new Atom.  You should call addAtom() on the Topology instead of calling this directly."""
        ## The name of the Atom
        self.name = name
        ## That Atom's element
        self.element = element
        ## The index of the Atom within its Topology
        self.index = index
        ## The Residue this Atom belongs to
        self.residue = residue
        ## A user defined identifier for this Atom
        self.id = id

    def __repr__(self):
        return "<Atom %d (%s) of chain %d residue %d (%s)>" % (self.index, self.name, self.residue.chain.index, self.residue.index, self.residue.name)

class Bond(namedtuple('Bond', ['atom1', 'atom2'])):
    """A Bond object represents a bond between two Atoms within a Topology.

    This class extends tuple, and may be interpreted as a 2 element tuple of Atom objects.
    It also has fields that can optionally be used to describe the bond order and type of bond."""

    def __new__(cls, atom1, atom2, type=None, order=None):
        """Create a new Bond.  You should call addBond() on the Topology instead of calling this directly."""
        bond = super(Bond, cls).__new__(cls, atom1, atom2)
        bond.type = type
        bond.order = order
        return bond

    def __getnewargs__(self):
        "Support for pickle protocol 2: http://docs.python.org/2/library/pickle.html#pickling-and-unpickling-normal-class-instances"
        return self[0], self[1], self.type, self.order

    def __getstate__(self):
        """
        Additional support for pickle since parent class implements its own __getstate__
        so pickle does not store or restore the type and order, python 2 problem only
        https://www.python.org/dev/peps/pep-0307/#case-3-pickling-new-style-class-instances-using-protocol-2
        """
        return self.__dict__

    def __deepcopy__(self, memo):
        return Bond(self[0], self[1], self.type, self.order)

    def __repr__(self):
        s = "Bond(%s, %s" % (self[0], self[1])
        if self.type is not None:
            s = "%s, type=%s" % (s, self.type)
        if self.order is not None:
            s = "%s, order=%d" % (s, self.order)
        s += ")"
        return s
