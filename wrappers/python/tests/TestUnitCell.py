import unittest
import math
from simtk.openmm import Vec3
from simtk.unit import *
from simtk.openmm.app.internal.unitcell import computePeriodicBoxVectors
from simtk.openmm.app.internal.unitcell import computeLengthsAndAngles
from simtk.openmm.app.internal.unitcell import reducePeriodicBoxVectors

def strip_units(x):
    if is_quantity(x): return x.value_in_unit_system(md_unit_system)
    return x

class TestUnitCell(unittest.TestCase):

    """ Test the unitcell.py module """

    def testReducePBCVectors(self):
        """ Checks that reducePeriodicBoxVectors properly reduces vectors """
        a = Vec3(4.24388485, 0.0, 0.0)
        b = Vec3(-1.4146281691908937, 4.001173048368583, 0.0)
        c = Vec3(-1.4146281691908937, -2.0005862820516203, 3.4651176446201674)
        vecs = reducePeriodicBoxVectors((a, b, c)*nanometers)
        vecs2 = computePeriodicBoxVectors(4.24388485, 4.24388485, 4.24388485,
                109.4712190*degrees, 109.4712190*degrees, 109.4712190*degrees)

        # Check that the vectors are the same
        a1, a2, a3 = vecs
        b1, b2, b3 = vecs2
        for x, y in zip(a1, b1):
            self.assertAlmostEqual(strip_units(x), strip_units(y))
        for x, y in zip(a2, b2):
            self.assertAlmostEqual(strip_units(x), strip_units(y))
        for x, y in zip(a3, b3):
            self.assertAlmostEqual(strip_units(x), strip_units(y))

    def testComputePBCVectors(self):
        """ Tests computing periodic box vectors """
        deg90 = 90 * degrees
        vecs = computePeriodicBoxVectors(1, 2, 3, deg90, deg90, deg90)
        a, b, c = vecs
        self.assertAlmostEqual(a[0]/nanometers, 1)
        self.assertAlmostEqual(a[1]/nanometers, 0)
        self.assertAlmostEqual(a[2]/nanometers, 0)
        self.assertAlmostEqual(b[0]/nanometers, 0)
        self.assertAlmostEqual(b[1]/nanometers, 2)
        self.assertAlmostEqual(b[2]/nanometers, 0)
        self.assertAlmostEqual(c[0]/nanometers, 0)
        self.assertAlmostEqual(c[1]/nanometers, 0)
        self.assertAlmostEqual(c[2]/nanometers, 3)

        # Make sure round-trip works
        la, lb, lc, al, be, ga = computeLengthsAndAngles(vecs)
        self.assertAlmostEqual(la, 1)
        self.assertAlmostEqual(lb, 2)
        self.assertAlmostEqual(lc, 3)
        self.assertAlmostEqual(al, math.pi / 2)
        self.assertAlmostEqual(be, math.pi / 2)
        self.assertAlmostEqual(ga, math.pi / 2)

        # Now test a truncated octahedron. Can't do a simple round-trip though,
        # due to the reduced form. So test the *second* round-trip, which should
        # yield the same measurements
        vecs = computePeriodicBoxVectors(4.24388485, 4.24388485, 4.24388485,
                109.4712190*degrees, 109.4712190*degrees, 109.4712190*degrees)
        la, lb, lc, al, be, ga = computeLengthsAndAngles(vecs)
        vecs2 = computePeriodicBoxVectors(la, lb, lc, al, be, ga)
        la2, lb2, lc2, al2, be2, ga2 = computeLengthsAndAngles(vecs2)

        # Now make sure that the round-trip worked
        self.assertAlmostEqual(strip_units(la), strip_units(la2))
        self.assertAlmostEqual(strip_units(lb), strip_units(lb2))
        self.assertAlmostEqual(strip_units(lc), strip_units(lc2))
        self.assertAlmostEqual(strip_units(al), strip_units(al2))
        self.assertAlmostEqual(strip_units(be), strip_units(be2))
        self.assertAlmostEqual(strip_units(ga), strip_units(ga2))

        # Check that the vectors are the same
        a1, a2, a3 = vecs
        b1, b2, b3 = vecs2
        for x, y in zip(a1, b1):
            self.assertAlmostEqual(strip_units(x), strip_units(y))
        for x, y in zip(a2, b2):
            self.assertAlmostEqual(strip_units(x), strip_units(y))
        for x, y in zip(a3, b3):
            self.assertAlmostEqual(strip_units(x), strip_units(y))
