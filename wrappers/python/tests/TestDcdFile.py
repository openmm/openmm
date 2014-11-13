import unittest
import tempfile
import numpy as np
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from random import random
import os

class TestDCDFile(unittest.TestCase):
    def test_dcd(self):
        """ Test the DCD file """
        pdbfile = app.PDBFile('systems/alanine-dipeptide-implicit.pdb')
        natom = len(list(pdbfile.topology.atoms()))
        dcd = app.DCDFile(open('test.dcd', 'w'), pdbfile.topology, 0.001)
        for i in range(5):
            dcd.writeModel([mm.Vec3(random(), random(), random()) for j in range(natom)]*unit.angstroms)
        dcd._file.close()
        os.remove('test.dcd')

if __name__ == '__main__':
    unittest.main()
