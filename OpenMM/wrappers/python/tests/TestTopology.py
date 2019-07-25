import pickle
import sys
import unittest
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import simtk.openmm.app.element as elem
if sys.version_info >= (3, 0):
    from io import StringIO
else:
    from cStringIO import StringIO


class TestTopology(unittest.TestCase):
    """Test the Topology object"""


    def check_pdbfile(self, pdbfilename, natoms, nres, nchains):
        """Check that a PDB file has the specified number of atoms, residues, and chains."""
        pdb = PDBFile(pdbfilename)
        top = pdb.topology
        self.assertEqual(pdb.topology.getNumAtoms(), natoms)
        self.assertEqual(pdb.topology.getNumResidues(), nres)
        self.assertEqual(pdb.topology.getNumChains(), nchains)

    def test_getters(self):
        """Test getters for number of atoms, residues, chains."""
        self.check_pdbfile('systems/1T2Y.pdb', 271, 25, 1)

    def test_bondtype_singleton(self):
        """ Tests that the bond types are really singletons """
        self.assertIs(Single, pickle.loads(pickle.dumps(Single)))
        self.assertIs(Double, pickle.loads(pickle.dumps(Double)))
        self.assertIs(Triple, pickle.loads(pickle.dumps(Triple)))
        self.assertIs(Aromatic, pickle.loads(pickle.dumps(Aromatic)))
        self.assertIs(Amide, pickle.loads(pickle.dumps(Amide)))

    def test_residue_bonds(self):
        """Test retrieving bonds for a residue produces expected results."""
        # Create a test topology
        # atom connectivity = A1-|-B1-B2-|-C1
        topology = Topology()
        chain = topology.addChain(id='A')
        residue1 = topology.addResidue('AAA', chain)
        residue2 = topology.addResidue('BBB', chain)
        residue3 = topology.addResidue('CCC', chain)
        atom_A1 = topology.addAtom('A1', element.carbon, residue1)
        atom_B1 = topology.addAtom('B1', element.carbon, residue2)
        atom_B2 = topology.addAtom('B2', element.carbon, residue2)
        atom_C1 = topology.addAtom('C1', element.carbon, residue3)
        topology.addBond(atom_A1, atom_B1)
        topology.addBond(atom_B1, atom_B2)
        topology.addBond(atom_B2, atom_C1)
        # Check bonds
        all_bonds = [ bond for bond in residue2.bonds() ]
        internal_bonds = [ bond for bond in residue2.internal_bonds() ]
        external_bonds = [ bond for bond in residue2.external_bonds() ]
        self.assertEqual(all_bonds, [ (atom_A1, atom_B1), (atom_B1, atom_B2), (atom_B2, atom_C1) ])
        self.assertEqual(internal_bonds, [ (atom_B1, atom_B2) ])
        self.assertEqual(external_bonds, [ (atom_A1, atom_B1), (atom_B2, atom_C1) ])

if __name__ == '__main__':
    unittest.main()
