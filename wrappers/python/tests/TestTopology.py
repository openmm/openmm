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

if __name__ == '__main__':
    unittest.main()
