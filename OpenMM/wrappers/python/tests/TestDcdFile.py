import unittest
import tempfile
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from random import random
import os

class TestDCDFile(unittest.TestCase):
    def test_dcd(self):
        """ Test the DCD file """
        fname = tempfile.mktemp(suffix='.dcd')
        pdbfile = app.PDBFile('systems/alanine-dipeptide-implicit.pdb')
        natom = len(list(pdbfile.topology.atoms()))
        with open(fname, 'wb') as f:
            dcd = app.DCDFile(f, pdbfile.topology, 0.001)
            for i in range(5):
                dcd.writeModel([mm.Vec3(random(), random(), random()) for j in range(natom)]*unit.angstroms)
        os.remove(fname)
    
    def testLongTrajectory(self):
        """Test writing a trajectory that has more than 2^31 steps."""
        fname = tempfile.mktemp(suffix='.dcd')
        pdbfile = app.PDBFile('systems/alanine-dipeptide-implicit.pdb')
        natom = len(list(pdbfile.topology.atoms()))
        with open(fname, 'wb') as f:
            dcd = app.DCDFile(f, pdbfile.topology, 0.001, interval=1000000000)
            for i in range(5):
                dcd.writeModel([mm.Vec3(random(), random(), random()) for j in range(natom)]*unit.angstroms)
        os.remove(fname)

if __name__ == '__main__':
    unittest.main()
