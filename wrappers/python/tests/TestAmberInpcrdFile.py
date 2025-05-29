import unittest
from validateConstraints import *
from openmm.app import *
from openmm import *
from openmm.unit import *
import openmm.app.element as elem
try:
    from scipy.io import netcdf_file
    SCIPY_IMPORT_FAILED = False
except ImportError:
    SCIPY_IMPORT_FAILED = True

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

    @unittest.skipIf(SCIPY_IMPORT_FAILED, "Scipy is not installed")
    def test_NetCDF(self):
        """ Test NetCDF restart file parsing """
        cmp = self.assertAlmostEqual

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

    def test_CrdBoxTruncoct(self):
        # Check that the box vectors come out correct.
        inpcrd = AmberInpcrdFile('systems/tz2.truncoct.rst7')
        ac = Vec3(42.4388485, 0.0, 0.0) * angstroms
        bc = Vec3(-14.146281691908937, 40.011730483685835, 0.0) * angstroms
        cc = Vec3(-14.146281691908937, -20.0058628205162, 34.651176446201672) * angstroms
        a, b, c = inpcrd.getBoxVectors()
        diffa = ac - a
        diffb = bc - b
        diffc = cc - c
        self.assertAlmostEqual(norm(diffa)/angstroms, 0)
        self.assertAlmostEqual(norm(diffb)/angstroms, 0)
        self.assertAlmostEqual(norm(diffc)/angstroms, 0)
        # Make sure angles and lengths come out about right
        la = norm(a).in_units_of(angstroms)
        lb = norm(b).in_units_of(angstroms)
        lc = norm(c).in_units_of(angstroms)
        self.assertAlmostEqual(la/angstroms, 42.4388485, 6)
        self.assertAlmostEqual(lb/angstroms, 42.4388485, 6)
        self.assertAlmostEqual(lc/angstroms, 42.4388485, 6)
        self.assertAlmostEqual(dot(a,b)/la/lb, cos(109.4712190*degrees), 6)
        self.assertAlmostEqual(dot(a,c)/la/lc, cos(109.4712190*degrees), 6)
        self.assertAlmostEqual(dot(b,c)/lc/lb, cos(109.4712190*degrees), 6)

if __name__ == '__main__':
    unittest.main()
