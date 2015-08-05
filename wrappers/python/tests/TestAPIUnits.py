import unittest
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
import math

class TestAPIUnits(unittest.TestCase):
    """Test the Simulation class"""

    def assertAlmostEqualUnit(self, x1, x2):
        self.assertAlmostEqual(x1._value, x2.value_in_unit(x1.unit))

    def testHarmonicBondForce(self):
        """ Tests HarmonicBondForce API features """
        force = HarmonicBondForce()
        force.addBond(0, 1, 1.0, 1.0)
        force.addBond(2, 3, 1.0*angstroms,
                      1.0*kilocalories_per_mole/angstroms**2)
        i, j, length, K = force.getBondParameters(0)
        self.assertEqual(force.getNumBonds(), 2)
        self.assertEqual(i, 0)
        self.assertEqual(j, 1)
        self.assertEqual(length, 1.0*nanometers)
        self.assertEqual(K, 1*kilojoules_per_mole/nanometers**2)
        self.assertIs(length.unit, nanometers)
        self.assertIs(K.unit, kilojoules_per_mole/nanometers**2)

        i, j, length, K = force.getBondParameters(1)
        self.assertEqual(i, 2)
        self.assertEqual(j, 3)
        self.assertEqual(length, 1.0*angstroms)
        self.assertAlmostEqualUnit(K, 1*kilocalories_per_mole/angstroms**2)
        self.assertIs(length.unit, nanometers)
        self.assertIs(K.unit, kilojoules_per_mole/nanometers**2)

    def testHarmonicAngleForce(self):
        """ Tests HarmonicAngleForce API features """
        force = HarmonicAngleForce()
        force.addAngle(0, 1, 2, 180*degrees,
                       1.0*kilocalories_per_mole/radians**2)
        force.addAngle(1, 2, 3, math.pi/2, 1.0)
        self.assertEqual(force.getNumAngles(), 2)
        i, j, k, angle, K = force.getAngleParameters(0)
        self.assertEqual(i, 0)
        self.assertEqual(j, 1)
        self.assertEqual(k, 2)
        self.assertEqual(angle, 180*degrees)
        self.assertEqual(K, 1.0*kilocalories_per_mole/radians**2)
        self.assertIs(angle.unit, radians)
        self.assertIs(K.unit, kilojoules_per_mole/radians**2)

        i, j, k, angle, K = force.getAngleParameters(1)
        self.assertEqual(i, 1)
        self.assertEqual(j, 2)
        self.assertEqual(k, 3)
        self.assertEqual(angle, math.pi/2*radians)
        self.assertEqual(K, 1.0*kilojoules_per_mole/radians**2)
        self.assertIs(angle.unit, radians)
        self.assertIs(K.unit, kilojoules_per_mole/radians**2)

    def testPeriodicTorsionForce(self):
        """ Tests PeriodicTorsionForce API features """
        force = PeriodicTorsionForce()
        force.addTorsion(0, 1, 2, 3, 1, math.pi, 1)
        force.addTorsion(1, 2, 3, 4, 2, 180*degrees, 1*kilocalories_per_mole)

        self.assertEqual(force.getNumTorsions(), 2)
        i, j, k, l, per, phase, K = force.getTorsionParameters(0)
        self.assertEqual(i, 0)
        self.assertEqual(j, 1)
        self.assertEqual(k, 2)
        self.assertEqual(l, 3)
        self.assertEqual(per, 1)
        self.assertFalse(is_quantity(per))
        self.assertEqual(phase, math.pi*radians)
        self.assertEqual(K, 1*kilojoules_per_mole)
        self.assertIs(phase.unit, radians)
        self.assertIs(K.unit, kilojoules_per_mole)

        i, j, k, l, per, phase, K = force.getTorsionParameters(1)
        self.assertEqual(i, 1)
        self.assertEqual(j, 2)
        self.assertEqual(k, 3)
        self.assertEqual(l, 4)
        self.assertEqual(per, 2)
        self.assertFalse(is_quantity(per))
        self.assertEqual(phase, 180*degrees)
        self.assertEqual(K, 1*kilocalories_per_mole)
        self.assertIs(phase.unit, radians)
        self.assertIs(K.unit, kilojoules_per_mole)

    def testRBTorsionForce(self):
        """ Tests the RBTorsionForce API features """
        force = RBTorsionForce()
        force.addTorsion(0, 1, 2, 3, 1, 2, 3, 4, 5, 6)
        force.addTorsion(1, 2, 3, 4, 1*kilocalories_per_mole,
                2*kilocalories_per_mole, 3*kilocalories_per_mole,
                4*kilocalories_per_mole, 5*kilocalories_per_mole,
                6*kilocalories_per_mole)

        self.assertEqual(force.getNumTorsions(), 2)
        i, j, k, l, c0, c1, c2, c3, c4, c5 = force.getTorsionParameters(0)
        self.assertEqual(i, 0)
        self.assertEqual(j, 1)
        self.assertEqual(k, 2)
        self.assertEqual(l, 3)
        self.assertEqual(c0, 1*kilojoules_per_mole)
        self.assertEqual(c1, 2*kilojoules_per_mole)
        self.assertEqual(c2, 3*kilojoules_per_mole)
        self.assertEqual(c3, 4*kilojoules_per_mole)
        self.assertEqual(c4, 5*kilojoules_per_mole)
        self.assertEqual(c5, 6*kilojoules_per_mole)
        self.assertIs(c0.unit, kilojoules_per_mole)
        self.assertIs(c1.unit, kilojoules_per_mole)
        self.assertIs(c2.unit, kilojoules_per_mole)
        self.assertIs(c3.unit, kilojoules_per_mole)
        self.assertIs(c4.unit, kilojoules_per_mole)
        self.assertIs(c5.unit, kilojoules_per_mole)

        i, j, k, l, c0, c1, c2, c3, c4, c5 = force.getTorsionParameters(1)
        self.assertEqual(i, 1)
        self.assertEqual(j, 2)
        self.assertEqual(k, 3)
        self.assertEqual(l, 4)
        self.assertAlmostEqualUnit(c0, 1*kilocalories_per_mole)
        self.assertAlmostEqualUnit(c1, 2*kilocalories_per_mole)
        self.assertAlmostEqualUnit(c2, 3*kilocalories_per_mole)
        self.assertAlmostEqualUnit(c3, 4*kilocalories_per_mole)
        self.assertAlmostEqualUnit(c4, 5*kilocalories_per_mole)
        self.assertAlmostEqualUnit(c5, 6*kilocalories_per_mole)
        self.assertIs(c0.unit, kilojoules_per_mole)
        self.assertIs(c1.unit, kilojoules_per_mole)
        self.assertIs(c2.unit, kilojoules_per_mole)
        self.assertIs(c3.unit, kilojoules_per_mole)
        self.assertIs(c4.unit, kilojoules_per_mole)
        self.assertIs(c5.unit, kilojoules_per_mole)

    def testNonbondedForce(self):
        """ Tests the NonbondedForce API features """
        force = NonbondedForce()
        force.addParticle(1.0, 1.0, 1.0)
        force.addParticle(1.0*coulombs, 1.0*angstroms,
                          1.0*kilocalories_per_mole)

        self.assertEqual(force.getNumParticles(), 2)
        charge, sigma, epsilon = force.getParticleParameters(0)
        self.assertEqual(charge, 1.0*elementary_charge)
        self.assertEqual(sigma, 1.0*nanometers)
        self.assertEqual(epsilon, 1.0*kilojoules_per_mole)
        self.assertIs(charge.unit, elementary_charge)
        self.assertIs(sigma.unit, nanometers)
        self.assertIs(epsilon.unit, kilojoules_per_mole)

        charge, sigma, epsilon = force.getParticleParameters(1)
        self.assertEqual(charge, 1.0*coulombs)
        self.assertEqual(sigma, 1.0*angstroms)
        self.assertEqual(epsilon, 1.0*kilocalories_per_mole)
        self.assertIs(charge.unit, elementary_charge)
        self.assertIs(sigma.unit, nanometers)
        self.assertIs(epsilon.unit, kilojoules_per_mole)

if __name__ == '__main__':
    unittest.main()
