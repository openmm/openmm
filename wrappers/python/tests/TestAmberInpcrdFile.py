import unittest
from validateConstraints import *
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import simtk.openmm.app.element as elem

def compareByElement(array1, array2, cmp):
    for x, y in zip(array1, array2):
        cmp(x, y)

class TestAmberInpcrdFile(unittest.TestCase):
    """Test the Amber inpcrd file parser"""
 
    def test_CrdVelBox(self):
        """ Test parsing ASCII restarts with crds, vels, and box """
        cmp = self.assertAlmostEqual
        inpcrd = AmberInpcrdFile('systems/crds_vels_box.rst7')
        self.assertEqual(len(inpcrd.positions), 2101)
        compareByElement(inpcrd.positions[-1].value_in_unit(angstroms),
                         [3.5958082, 8.4176792, -8.2954064], cmp)
        compareByElement(inpcrd.velocities[-1].value_in_unit(angstroms/picoseconds),
                         [0.3091332*20.455, 0.7355925*20.455, -0.1206961*20.455], cmp)
        compareByElement(inpcrd.boxVectors[0].value_in_unit(angstroms),
                         [30.2642725, 0.0, 0.0], cmp)

    def test_NetCDF(self):
        """ Test NetCDF restart file parsing """
        cmp = self.assertAlmostEqual
        try:
            from scipy.io import netcdf
        except ImportError:
            print('Not testing NetCDF file parser... scipy cannot be found')
        else:
            inpcrd = AmberInpcrdFile('systems/amber.ncrst')
            self.assertEqual(len(inpcrd.positions), 2101)
            compareByElement(inpcrd.positions[0].value_in_unit(angstroms),
                             [6.82122492718229, 6.6276250662042, -8.51668999892245],
                             cmp)
            compareByElement(inpcrd.velocities[-1].value_in_unit(angstroms/picosecond),
                             [0.349702202733541*20.455, 0.391525333168534*20.455,
                              0.417941679767662*20.455], cmp)
            self.assertAlmostEqual(inpcrd.boxVectors[0][0].value_in_unit(angstroms),
                                   30.2642725, places=6)

    def test_CrdBox(self):
        """ Test parsing ASCII restarts with only crds and box """
        inpcrd = AmberInpcrdFile('systems/crds_box.rst7')
        self.assertEqual(len(inpcrd.positions), 18660)
        self.assertTrue(inpcrd.velocities is None)
        self.assertTrue(inpcrd.boxVectors is not None)

    def test_CrdVel(self):
        inpcrd = AmberInpcrdFile('systems/crds_vels.rst7')
        self.assertTrue(inpcrd.boxVectors is None)
        self.assertTrue(inpcrd.velocities is not None)

    def test_CrdOnly(self):
        inpcrd = AmberInpcrdFile('systems/crdsonly.rst7')
        self.assertTrue(inpcrd.boxVectors is None)
        self.assertTrue(inpcrd.velocities is None)

if __name__ == '__main__':
    unittest.main()
