from __future__ import division

import unittest
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
import random
import math

class TestAPIUnits(unittest.TestCase):
    """Test the Simulation class"""

    def assertAlmostEqualUnit(self, x1, x2):
        self.assertAlmostEqual(x1._value, x2.value_in_unit(x1.unit))

    def assertAlmostEqualUnitArray(self, a1, a2):
        for x1, x2 in zip(a1._value, a2.value_in_unit(a1.unit)):
            self.assertAlmostEqual(x1, x2)

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

        force.setCutoffDistance(10*angstroms)
        self.assertEqual(force.getCutoffDistance(), 1*nanometers)
        force.setCutoffDistance(1)
        self.assertEqual(force.getCutoffDistance(), 1*nanometers)

        force.setSwitchingDistance(8*angstroms)
        self.assertEqual(force.getSwitchingDistance(), 0.8*nanometers)
        self.assertIs(force.getSwitchingDistance().unit, nanometer)

        self.assertFalse(force.usesPeriodicBoundaryConditions())
        force.setNonbondedMethod(NonbondedForce.NoCutoff)
        self.assertFalse(force.usesPeriodicBoundaryConditions())
        force.setNonbondedMethod(NonbondedForce.CutoffNonPeriodic)
        self.assertFalse(force.usesPeriodicBoundaryConditions())
        force.setNonbondedMethod(NonbondedForce.CutoffPeriodic)
        self.assertTrue(force.usesPeriodicBoundaryConditions())
        force.setNonbondedMethod(NonbondedForce.Ewald)
        self.assertTrue(force.usesPeriodicBoundaryConditions())
        force.setNonbondedMethod(NonbondedForce.PME)
        self.assertTrue(force.usesPeriodicBoundaryConditions())
        force.setNonbondedMethod(NonbondedForce.LJPME)
        self.assertTrue(force.usesPeriodicBoundaryConditions())

    def testCmapForce(self):
        """ Tests the CMAPTorsionForce API features """
        map1 = [random.random() for i in range(24*24)]
        map2 = [random.random() for i in range(12*12)] * kilocalories_per_mole
        force = CMAPTorsionForce()
        force.addMap(24, map1)
        force.addMap(12, map2)
        force.addTorsion(0, 0, 1, 2, 3, 1, 2, 3, 4)
        force.addTorsion(1, 5, 6, 7, 8, 6, 7, 8, 9)
        force.addTorsion(0, 10, 11, 12, 13, 11, 12, 13, 14)
        force.addTorsion(1, 15, 16, 17, 18, 16, 17, 18, 19)

        self.assertEqual(force.getNumTorsions(), 4)
        self.assertEqual(force.getNumMaps(), 2)
        self.assertEqual(force.getMapParameters(0)[0], 24)
        self.assertEqual(force.getMapParameters(1)[0], 12)
        self.assertIs(force.getMapParameters(0)[1].unit, kilojoules_per_mole)
        self.assertIs(force.getMapParameters(1)[1].unit, kilojoules_per_mole)

        for x, y in zip(force.getMapParameters(0)[1], map1):
            self.assertAlmostEqual(x.value_in_unit(kilojoules_per_mole), y)

        for x, y in zip(force.getMapParameters(1)[1], map2):
            self.assertAlmostEqualUnit(x, y)

    def testCustomBondForce(self):
        """ Tests the CustomBondForce API features """
        force = CustomBondForce('1/2*k*(r-r0)^2')
        force.addPerBondParameter('r0')
        force.addBond(0, 1, [0.1])
        force.addBond(1, 2, [1.0*angstroms])

        self.assertEqual(force.getNumBonds(), 2)
        i, j, req = force.getBondParameters(0)
        self.assertEqual(i, 0)
        self.assertEqual(j, 1)
        self.assertEqual(req[0], 0.1)

        i, j, req = force.getBondParameters(1)
        self.assertEqual(i, 1)
        self.assertEqual(j, 2)
        self.assertEqual(req[0], 0.1)

    def testCustomAngleForce(self):
        """ Tests the CustomAngleForce API features """
        force = CustomAngleForce('1/2*k*(theta-theta0)^2')
        force.addPerAngleParameter('theta0')
        force.addAngle(0, 1, 2, [math.pi / 2])
        force.addAngle(3, 4, 5, [90*degrees])

        self.assertEqual(force.getNumAngles(), 2)
        i, j, k, theta = force.getAngleParameters(0)
        self.assertEqual(i, 0)
        self.assertEqual(j, 1)
        self.assertEqual(k, 2)
        self.assertEqual(theta[0], math.pi / 2)

        i, j, k, theta = force.getAngleParameters(1)
        self.assertEqual(i, 3)
        self.assertEqual(j, 4)
        self.assertEqual(k, 5)
        self.assertEqual(theta[0], math.pi / 2)

    def testCustomTorsionForce(self):
        """ Tests the CustomTorsionForce API features """
        force = CustomTorsionForce('1/2*k*(theta-theta0)^2')
        force.addTorsion(0, 1, 2, 3, [math.pi])
        force.addTorsion(4, 5, 6, 7, [180*degrees])

        self.assertEqual(force.getNumTorsions(), 2)
        i, j, k, l, theta = force.getTorsionParameters(0)
        self.assertEqual(i, 0)
        self.assertEqual(j, 1)
        self.assertEqual(k, 2)
        self.assertEqual(l, 3)
        self.assertEqual(theta[0], math.pi)

        i, j, k, l, theta = force.getTorsionParameters(1)
        self.assertEqual(i, 4)
        self.assertEqual(j, 5)
        self.assertEqual(k, 6)
        self.assertEqual(l, 7)
        self.assertEqual(theta[0], math.pi)

    def testCustomCompoundBondForce(self):
        """ Tests the CustomCompoundBondForce API features """
        force = CustomCompoundBondForce(4, 'kb*distance(p1, p2)*distance(p2, p3)*distance(p3, p4)+'
                                           'ka*angle(p1, p2, p3)*angle(p2, p3, p4)+'
                                           'kd*dihedral(p1, p2, p3, p4)')
        force.addPerBondParameter('kb')
        force.addPerBondParameter('ka')
        force.addPerBondParameter('kd')
        force.addBond([0, 1, 2, 3], [1.0, 2.0, 3.0])
        force.addBond([4, 5, 6, 7], [1.0*kilocalories_per_mole/angstroms,
                                     2.0*kilocalories_per_mole/radians,
                                     3.0*kilocalories_per_mole]
        )

        self.assertEqual(force.getNumBonds(), 2)
        (i, j, k, l), (kb, ka, kd) = force.getBondParameters(0)
        self.assertEqual(i, 0)
        self.assertEqual(j, 1)
        self.assertEqual(k, 2)
        self.assertEqual(l, 3)
        self.assertEqual(kb, 1)
        self.assertEqual(ka, 2)
        self.assertEqual(kd, 3)

        (i, j, k, l), (kb, ka, kd) = force.getBondParameters(1)
        self.assertEqual(i, 4)
        self.assertEqual(j, 5)
        self.assertEqual(k, 6)
        self.assertEqual(l, 7)
        self.assertEqual(kb, 1*10*4.184)
        self.assertEqual(ka, 2*4.184)
        self.assertEqual(kd, 3*4.184)

    def testCustomExternalForce(self):
        """ Tests the CustomExternalForce API features """
        force = CustomExternalForce('1/2*k*k2*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
        force.addGlobalParameter('k', 10*kilocalories_per_mole/angstroms**2)
        force.addGlobalParameter('k2', 20)
        force.addPerParticleParameter('x0')
        force.addPerParticleParameter('y0')
        force.addPerParticleParameter('z0')
        force.addParticle(0, [1.0, 2.0, 3.0])
        force.addParticle(1, [1.0*angstroms, 2.0*angstroms, 3.0*angstroms])

        self.assertEqual(force.getNumParticles(), 2)
        self.assertEqual(force.getNumGlobalParameters(), 2)
        self.assertEqual(force.getGlobalParameterName(0), 'k')
        self.assertEqual(force.getGlobalParameterName(1), 'k2')
        self.assertEqual(force.getGlobalParameterDefaultValue(0), 1000*4.184)
        self.assertEqual(force.getGlobalParameterDefaultValue(1), 20)

        i, (x0, y0, z0) = force.getParticleParameters(0)
        self.assertEqual(i, 0)
        self.assertEqual(x0, 1)
        self.assertEqual(y0, 2)
        self.assertEqual(z0, 3)

        i, (x0, y0, z0) = force.getParticleParameters(1)
        self.assertEqual(i, 1)
        self.assertAlmostEqual(x0, 1/10)
        self.assertAlmostEqual(y0, 2/10)
        self.assertAlmostEqual(z0, 3/10)

    def testCustomGBForce(self):
        """ Tests the CustomGBForce API features """
        force = CustomGBForce()
        force.addPerParticleParameter('q')
        force.addPerParticleParameter('radius')
        force.addPerParticleParameter('scale')
        force.addGlobalParameter('extDiel', 78.5)
        force.addGlobalParameter('intDiel', 1.0)
        force.addComputedValue('I',
                'step(r+sr2-or1)*0.5*(1/L-1/U+0.25*(1/U^2-1/L^2)*'
                '(r-sr2^2/r)+0.5*log(L/U)/r+C);'
                'U=r+sr2; C=2*(1/or1-1/L)*step(sr2-r-or1);'
                'L=max(or1, D); D=abs(r-sr2); sr2=scale2*or2;'
                'or1=radius1-0.009; or2=radius2-0.009',
                CustomGBForce.ParticlePairNoExclusions
        )
        force.addComputedValue('B',
                '1/(1/or-tanh(1*psi-0.8*psi^2+4.85*psi^3)/radius);'
                'psi=I*or; or=radius-0.009',
                CustomGBForce.SingleParticle
        )
        force.addEnergyTerm('-138.935456*(1/intDiel-1/extDiel)*q1*q2/f;'
                            'f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))',
                            CustomGBForce.ParticlePair
        )

        force.setNonbondedMethod(CustomGBForce.CutoffPeriodic)

        force.addParticle([1.0, 0.1, 0.5])
        force.addParticle([-1.0*coulombs, 1.0*angstroms, 0.5])

        self.assertEqual(force.getNumParticles(), 2)
        self.assertEqual(force.getNumComputedValues(), 2)
        self.assertEqual(force.getNumEnergyTerms(), 1)
        self.assertEqual(force.getNumPerParticleParameters(), 3)
        self.assertEqual(force.getPerParticleParameterName(0), 'q')
        self.assertEqual(force.getPerParticleParameterName(1), 'radius')
        self.assertEqual(force.getPerParticleParameterName(2), 'scale')
        self.assertTrue(force.usesPeriodicBoundaryConditions())

        force.setNonbondedMethod(CustomGBForce.NoCutoff)
        self.assertFalse(force.usesPeriodicBoundaryConditions())
        force.setNonbondedMethod(CustomGBForce.CutoffNonPeriodic)
        self.assertFalse(force.usesPeriodicBoundaryConditions())

        q, rad, scale = force.getParticleParameters(0)
        self.assertEqual(q, 1.0)
        self.assertEqual(rad, 0.1)
        self.assertEqual(scale, 0.5)

        q, rad, scale = force.getParticleParameters(1)
        self.assertEqual(q, -6.24150962915265e+18) # very electronegative
        self.assertEqual(rad, 0.1)
        self.assertEqual(scale, 0.5)

        force.setCutoffDistance(12*angstroms)
        self.assertAlmostEqualUnit(force.getCutoffDistance(), 1.2*nanometers)
        force.setCutoffDistance(1)
        self.assertEqual(force.getCutoffDistance(), 1*nanometers)

    def testCustomHBondForce(self):
        """ Tests the CustomHbondForce API features """
        force = CustomHbondForce('kd*(distance(a1,d1)-r0)^2 + '
                                 'ka*(angle(a1,d1,d2)-theta0)^2')
        force.addPerAcceptorParameter('r0')
        force.addPerAcceptorParameter('ka')
        force.addPerDonorParameter('theta0')
        force.addPerDonorParameter('kd')
        force.setCutoffDistance(10*angstroms)

        self.assertEqual(force.getNumPerAcceptorParameters(), 2)
        self.assertEqual(force.getNumPerDonorParameters(), 2)
        self.assertEqual(force.getCutoffDistance(), 1*nanometers)
        self.assertIs(force.getCutoffDistance().unit, nanometer)

        force.addAcceptor(0, 1, 2, [0.2, 10.0])
        force.addAcceptor(3, -1, -1, [4*angstroms,
                                      20.0*kilocalories_per_mole/angstroms**2])
        force.addDonor(4, 5, 6, [math.pi, 30])
        force.addDonor(7, 8, -1, [180*degrees,
                                  40*kilocalories_per_mole/radians**2])

        self.assertEqual(force.getNumAcceptors(), 2)
        self.assertEqual(force.getNumDonors(), 2)

        i, j, k, (r0, ka) = force.getAcceptorParameters(0)
        self.assertEqual(i, 0)
        self.assertEqual(j, 1)
        self.assertEqual(k, 2)
        self.assertEqual(r0, 0.2)
        self.assertEqual(ka, 10)

        i, j, k, (r0, ka) = force.getAcceptorParameters(1)
        self.assertEqual(i, 3)
        self.assertEqual(j, -1)
        self.assertEqual(k, -1)
        self.assertEqual(r0, 0.4)
        self.assertEqual(ka, 20*4.184*100)

        i, j, k, (theta0, kd) = force.getDonorParameters(0)
        self.assertEqual(i, 4)
        self.assertEqual(j, 5)
        self.assertEqual(k, 6)
        self.assertEqual(theta0, math.pi)
        self.assertEqual(kd, 30)

        i, j, k, (theta0, kd) = force.getDonorParameters(1)
        self.assertEqual(i, 7)
        self.assertEqual(j, 8)
        self.assertEqual(k, -1)
        self.assertEqual(theta0, math.pi)
        self.assertEqual(kd, 40*4.184)

        self.assertFalse(force.usesPeriodicBoundaryConditions())
        force.setNonbondedMethod(CustomHbondForce.NoCutoff)
        self.assertFalse(force.usesPeriodicBoundaryConditions())
        force.setNonbondedMethod(CustomHbondForce.CutoffNonPeriodic)
        self.assertFalse(force.usesPeriodicBoundaryConditions())
        force.setNonbondedMethod(CustomHbondForce.CutoffPeriodic)
        self.assertTrue(force.usesPeriodicBoundaryConditions())

    def testCustomNonbondedForce(self):
        """ Tests the CustomNonbondedForce API features """
        force = CustomNonbondedForce('4*diel*q1*q2/r+sqrt(p1*p2)/r^2-sqrt(m1*m2)/r^3')
        force.addGlobalParameter('diel', 1.0)
        force.addPerParticleParameter('q')
        force.addPerParticleParameter('p')
        force.addPerParticleParameter('m')
        force.addParticle([1, 2, 3])
        force.addParticle([1*coulombs, 2*kilocalories_per_mole*angstroms**2, 3*kilocalories_per_mole*angstroms**3])

        self.assertEqual(force.getNumParticles(), 2)
        charge, sigma, epsilon = force.getParticleParameters(0)
        self.assertEqual(charge, 1)
        self.assertEqual(sigma, 2)
        self.assertEqual(epsilon, 3)

        charge, sigma, epsilon = force.getParticleParameters(1)
        self.assertEqual(charge, 6.24150962915265e+18) # very electronegative
        self.assertEqual(sigma, 2*4.184/100)
        self.assertAlmostEqual(epsilon, 3*4.184/1000)

        force.setCutoffDistance(10*angstroms)
        self.assertEqual(force.getCutoffDistance(), 1*nanometers)
        force.setCutoffDistance(1)
        self.assertEqual(force.getCutoffDistance(), 1*nanometers)

        force.setSwitchingDistance(8*angstroms)
        self.assertEqual(force.getSwitchingDistance(), 0.8*nanometers)
        self.assertIs(force.getSwitchingDistance().unit, nanometer)

        self.assertFalse(force.usesPeriodicBoundaryConditions())
        force.setNonbondedMethod(CustomNonbondedForce.NoCutoff)
        self.assertFalse(force.usesPeriodicBoundaryConditions())
        force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
        self.assertFalse(force.usesPeriodicBoundaryConditions())
        force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
        self.assertTrue(force.usesPeriodicBoundaryConditions())

    def testCustomManyParticleForce(self):
        """ Tests the CustomManyParticleForce API features """
        force = CustomManyParticleForce(3,
                "C*(1+3*cos(theta1)*cos(theta2)*cos(theta3))/(r12*r13*r23)^3;"
                "theta1=k1*angle(p1,p2,p3); theta2=k2*angle(p2,p3,p1); theta3=k3*angle(p3,p1,p2);"
                "r12=distance(p1,p2); r13=distance(p1,p3); r23=distance(p2,p3)")
        force.setPermutationMode(CustomManyParticleForce.SinglePermutation)
        force.setTypeFilter(0, [0])
        force.setTypeFilter(1, [1])
        force.setTypeFilter(2, [2])

        force.addGlobalParameter('C', 1.0*kilocalories_per_mole)
        force.addPerParticleParameter('k')

        self.assertEqual(force.getNumGlobalParameters(), 1)
        self.assertEqual(force.getGlobalParameterName(0), 'C')
        self.assertEqual(force.getGlobalParameterDefaultValue(0), 4.184)
        self.assertEqual(force.getNumPerParticleParameters(), 1)

        force.addParticle([10], 0)
        force.addParticle([20], 1)
        force.addParticle([30*kilocalories_per_mole], 2)

        self.assertEqual(force.getNumParticles(), 3)
        self.assertEqual(force.getParticleParameters(0)[0][0], 10)
        self.assertEqual(force.getParticleParameters(1)[0][0], 20)
        self.assertEqual(force.getParticleParameters(2)[0][0], 30*4.184)

    def testDrudeForce(self):
        """ Tests the DrudeForce API features """
        force = DrudeForce()
        self.assertFalse(force.usesPeriodicBoundaryConditions())
        force.addParticle(0, 1, -1, -1, -1, 1, 1, 0, 0)
        force.addParticle(1, 2, 3, -1, -1, 1*elementary_charge, 1*angstrom**3, 0.5, 0)
        force.addParticle(2, 3, 4, 5, 6, 1*elementary_charge, 10*angstrom**3, 0.5, 0.5)
        force.addScreenedPair(0, 1, 0.5)
        force.addScreenedPair(1, 2, 0.25)
        force.addScreenedPair(0, 2, 0.125)

        self.assertEqual(force.getNumParticles(), 3)
        self.assertEqual(force.getNumScreenedPairs(), 3)

        i, j, k, l, m, q, a, an12, an34 = force.getParticleParameters(0)

        self.assertEqual(i, 0)
        self.assertEqual(j, 1)
        self.assertEqual(k, -1)
        self.assertEqual(l, -1)
        self.assertEqual(m, -1)
        self.assertEqual(q, 1*elementary_charge)
        self.assertEqual(a, 1*nanometer**3)
        self.assertEqual(an12, 0)
        self.assertEqual(an34, 0)

        i, j, k, l, m, q, a, an12, an34 = force.getParticleParameters(1)

        self.assertEqual(i, 1)
        self.assertEqual(j, 2)
        self.assertEqual(k, 3)
        self.assertEqual(l, -1)
        self.assertEqual(m, -1)
        self.assertEqual(q, 1*elementary_charge)
        self.assertAlmostEqualUnit(a, 1*angstrom**3)
        self.assertEqual(an12, 0.5)
        self.assertEqual(an34, 0)

        i, j, k, l, m, q, a, an12, an34 = force.getParticleParameters(2)

        self.assertEqual(i, 2)
        self.assertEqual(j, 3)
        self.assertEqual(k, 4)
        self.assertEqual(l, 5)
        self.assertEqual(m, 6)
        self.assertEqual(q, 1*elementary_charge)
        self.assertAlmostEqualUnit(a, 10*angstrom**3)
        self.assertEqual(an12, 0.5)
        self.assertEqual(an34, 0.5)

        i, j, thole = force.getScreenedPairParameters(0)
        self.assertEqual(i, 0)
        self.assertEqual(j, 1)
        self.assertEqual(thole, 0.5)

        i, j, thole = force.getScreenedPairParameters(1)
        self.assertEqual(i, 1)
        self.assertEqual(j, 2)
        self.assertEqual(thole, 0.25)

        i, j, thole = force.getScreenedPairParameters(2)
        self.assertEqual(i, 0)
        self.assertEqual(j, 2)
        self.assertEqual(thole, 0.125)

    def testAmoebaBondForce(self):
        """ Tests the AmoebaBondForce API features """
        force1 = AmoebaBondForce()
        force2 = AmoebaBondForce()
        force1.setAmoebaGlobalBondCubic(1.5)
        force2.setAmoebaGlobalBondCubic(1.5/angstrom)
        force1.setAmoebaGlobalBondQuartic(1.5)
        force2.setAmoebaGlobalBondQuartic(1.5/angstrom**2)

        self.assertEqual(force1.getAmoebaGlobalBondCubic(), 1.5/nanometer)
        self.assertEqual(force2.getAmoebaGlobalBondCubic(), 1.5/angstrom)
        self.assertEqual(force1.getAmoebaGlobalBondQuartic(), 1.5/nanometer**2)
        self.assertAlmostEqualUnit(force2.getAmoebaGlobalBondQuartic(), 1.5/angstrom**2)

        force1.addBond(0, 1, 0.15, 10.0)
        force1.addBond(1, 2, 1.5*angstroms, 10.0*kilocalories_per_mole/angstroms**2)

        self.assertEqual(force1.getNumBonds(), 2)
        self.assertEqual(force2.getNumBonds(), 0)

        i, j, req, k = force1.getBondParameters(0)
        self.assertEqual(req, 0.15*nanometers)
        self.assertEqual(k, 10.0*kilojoules_per_mole/nanometers**2)

        i, j, req, k = force1.getBondParameters(1)
        self.assertAlmostEqualUnit(req, 1.5*angstroms)
        self.assertAlmostEqualUnit(k, 10.0*kilocalories_per_mole/angstroms**2)

    def testAmoebaAngleForce(self):
        """ Tests the AmoebaAngleForce API features """
        force = AmoebaAngleForce()
        force.setAmoebaGlobalAngleCubic(1.0)
        force.setAmoebaGlobalAngleQuartic(2.0/radians**2)
        force.setAmoebaGlobalAnglePentic(3.0/degrees**3)
        force.setAmoebaGlobalAngleSextic(4.0)

        self.assertEqual(force.getAmoebaGlobalAngleCubic(), 1.0/radians)
        self.assertEqual(force.getAmoebaGlobalAngleQuartic(), 2.0/radians**2)
        self.assertAlmostEqualUnit(force.getAmoebaGlobalAnglePentic(), 3.0/degrees**3)
        self.assertEqual(force.getAmoebaGlobalAngleSextic(), 4.0/radians**4)

        force.addAngle(0, 1, 2, math.pi*radians, 1.5*kilocalories_per_mole/radians**2)
        force.addAngle(1, 2, 3, 180*degrees, 1.5*kilocalories_per_mole/radians**2)
        force.addAngle(2, 3, 4, 109.4, 1.5)

        self.assertEqual(force.getNumAngles(), 3)

        i, j, k, t, tk = force.getAngleParameters(0)
        self.assertEqual(i, 0)
        self.assertEqual(j, 1)
        self.assertEqual(k, 2)
        self.assertAlmostEqualUnit(t, math.pi*radians)
        self.assertIs(t.unit, degree)
        self.assertAlmostEqualUnit(tk, 1.5*kilocalories_per_mole/radians**2)
        self.assertIs(tk.unit, kilojoules_per_mole/radians**2)

        i, j, k, t, tk = force.getAngleParameters(1)
        self.assertEqual(i, 1)
        self.assertEqual(j, 2)
        self.assertEqual(k, 3)
        self.assertAlmostEqualUnit(t, 180*degrees)
        self.assertIs(t.unit, degree)
        self.assertAlmostEqualUnit(tk, 1.5*kilocalories_per_mole/radians**2)
        self.assertIs(tk.unit, kilojoules_per_mole/radians**2)

        i, j, k, t, tk = force.getAngleParameters(2)
        self.assertEqual(i, 2)
        self.assertEqual(j, 3)
        self.assertEqual(k, 4)
        self.assertAlmostEqualUnit(t, 109.4*degrees)
        self.assertIs(t.unit, degree)
        self.assertAlmostEqualUnit(tk, 1.5*kilojoules_per_mole/radians**2)
        self.assertIs(tk.unit, kilojoules_per_mole/radians**2)

    def testGeneralizedKirkwood(self):
        """ Tests the AmoebaGeneralizedKirkwoodForce API features """
        force = AmoebaGeneralizedKirkwoodForce()

        self.assertEqual(force.getProbeRadius(), 0.14*nanometer) # default
        force.setProbeRadius(0.16)
        self.assertEqual(force.getProbeRadius(), 0.16*nanometer)
        self.assertIs(force.getProbeRadius().unit, nanometer)
        force.setProbeRadius(1.4*angstrom)
        self.assertEqual(force.getProbeRadius(), 1.4*angstrom)
        self.assertIs(force.getProbeRadius().unit, nanometer)

        self.assertEqual(force.getSoluteDielectric(), 1.0) # default
        force.setSoluteDielectric(2.0)
        self.assertEqual(force.getSoluteDielectric(), 2.0)

        self.assertEqual(force.getSolventDielectric(), 78.3) # default
        force.setSolventDielectric(80)
        self.assertEqual(force.getSolventDielectric(), 80)

        self.assertEqual(force.getSurfaceAreaFactor(),
                -170.35173066268223*kilojoule_per_mole/nanometer**2) # default
        force.setSurfaceAreaFactor(-1.0*kilocalorie_per_mole/angstrom**2)
        self.assertAlmostEqualUnit(force.getSurfaceAreaFactor(),
                -1.0*kilocalorie_per_mole/angstrom**2) # default

        force.addParticle(1.0*coulomb, 1.0*angstroms, 0.5)
        force.addParticle(1.0, 1.0, 0.4)

        self.assertEqual(force.getNumParticles(), 2)

        q, r, s = force.getParticleParameters(0)
        self.assertAlmostEqualUnit(q, 1.0*coulomb)
        self.assertIs(q.unit, elementary_charge)
        self.assertEqual(r, 1.0*angstroms)
        self.assertIs(r.unit, nanometer)
        self.assertEqual(s, 0.5)

        q, r, s = force.getParticleParameters(1)
        self.assertAlmostEqualUnit(q, 1.0*elementary_charge)
        self.assertIs(q.unit, elementary_charge)
        self.assertEqual(r, 1.0*nanometer)
        self.assertIs(r.unit, nanometer)
        self.assertEqual(s, 0.4)

    def testAmoebaInPlaneAngleForce(self):
        """ Tests the AmoebaInPlaneAngleForce API features """
        force = AmoebaInPlaneAngleForce()

        force.setAmoebaGlobalInPlaneAngleCubic(1.0)
        self.assertEqual(force.getAmoebaGlobalInPlaneAngleCubic(), 1/radian)
        self.assertEqual(str(force.getAmoebaGlobalInPlaneAngleCubic().unit), '/radian')

        force.setAmoebaGlobalInPlaneAngleQuartic(1.0/degrees**2)
        self.assertAlmostEqualUnit(force.getAmoebaGlobalInPlaneAngleQuartic(), 1/degrees**2)
        self.assertEqual(str(force.getAmoebaGlobalInPlaneAngleQuartic().unit), '/(radian**2)')

        force.setAmoebaGlobalInPlaneAnglePentic(1.0/radians**3)
        self.assertEqual(force.getAmoebaGlobalInPlaneAnglePentic(), 1/radian**3)
        self.assertEqual(str(force.getAmoebaGlobalInPlaneAnglePentic().unit), '/(radian**3)')

        force.setAmoebaGlobalInPlaneAngleSextic(1.0/radians**4)
        self.assertEqual(force.getAmoebaGlobalInPlaneAngleSextic(), 1/radian**4)
        self.assertEqual(str(force.getAmoebaGlobalInPlaneAngleSextic().unit), '/(radian**4)')

        force.addAngle(0, 1, 2, 3, math.pi, 1.0)
        force.addAngle(1, 2, 3, 4, 180*degrees, 1.0*kilocalories_per_mole/radians**2)

        self.assertEqual(force.getNumAngles(), 2)

        i, j, k, l, t, tk = force.getAngleParameters(0)
        self.assertEqual(i, 0)
        self.assertEqual(j, 1)
        self.assertEqual(k, 2)
        self.assertEqual(l, 3)
        self.assertEqual(t, math.pi*radians)
        self.assertIs(t.unit, radians)
        self.assertEqual(tk, 1.0*kilojoules_per_mole/radians**2)
        self.assertIs(tk.unit, kilojoules_per_mole/radians**2)

        i, j, k, l, t, tk = force.getAngleParameters(1)
        self.assertEqual(i, 1)
        self.assertEqual(j, 2)
        self.assertEqual(k, 3)
        self.assertEqual(l, 4)
        self.assertEqual(t, 180*degrees)
        self.assertIs(t.unit, radians)
        self.assertEqual(tk, 1.0*kilocalorie_per_mole/radians**2)
        self.assertIs(tk.unit, kilojoules_per_mole/radians**2)

    def testAmoebaOutOfPlaneBendForce(self):
        """ Tests the AmoebaOutOfPlaneBendForce API features """
        force = AmoebaOutOfPlaneBendForce()

        force.setAmoebaGlobalOutOfPlaneBendCubic(1.0)
        self.assertEqual(force.getAmoebaGlobalOutOfPlaneBendCubic(), 1/radian)
        self.assertEqual(str(force.getAmoebaGlobalOutOfPlaneBendCubic().unit), '/radian')

        force.setAmoebaGlobalOutOfPlaneBendQuartic(1.0/degrees**2)
        self.assertAlmostEqualUnit(force.getAmoebaGlobalOutOfPlaneBendQuartic(), 1/degrees**2)
        self.assertEqual(str(force.getAmoebaGlobalOutOfPlaneBendQuartic().unit), '/(radian**2)')

        force.setAmoebaGlobalOutOfPlaneBendPentic(1.0/radians**3)
        self.assertEqual(force.getAmoebaGlobalOutOfPlaneBendPentic(), 1/radian**3)
        self.assertEqual(str(force.getAmoebaGlobalOutOfPlaneBendPentic().unit), '/(radian**3)')

        force.setAmoebaGlobalOutOfPlaneBendSextic(1.0/radians**4)
        self.assertEqual(force.getAmoebaGlobalOutOfPlaneBendSextic(), 1/radian**4)
        self.assertEqual(str(force.getAmoebaGlobalOutOfPlaneBendSextic().unit), '/(radian**4)')

        force.addOutOfPlaneBend(0, 1, 2, 3, 1.0)
        force.addOutOfPlaneBend(1, 2, 3, 4, 1.0*kilocalories_per_mole/radians**2)

        self.assertEqual(force.getNumOutOfPlaneBends(), 2)

        i, j, k, l, tk = force.getOutOfPlaneBendParameters(0)
        self.assertEqual(i, 0)
        self.assertEqual(j, 1)
        self.assertEqual(k, 2)
        self.assertEqual(l, 3)
        self.assertEqual(tk, 1.0*kilojoules_per_mole/radians**2)
        self.assertIs(tk.unit, kilojoules_per_mole/radians**2)

        i, j, k, l, tk = force.getOutOfPlaneBendParameters(1)
        self.assertEqual(i, 1)
        self.assertEqual(j, 2)
        self.assertEqual(k, 3)
        self.assertEqual(l, 4)
        self.assertEqual(tk, 1.0*kilocalorie_per_mole/radians**2)
        self.assertIs(tk.unit, kilojoules_per_mole/radians**2)

    def testAmoebaPiTorsionForce(self):
        """ Tests the AmoebaPiTorsionForce API features """
        force = AmoebaPiTorsionForce()

        force.addPiTorsion(0, 1, 2, 3, 4, 5, 1.0)
        force.addPiTorsion(1, 2, 3, 4, 5, 6, 1.0*kilocalories_per_mole)

        self.assertEqual(force.getNumPiTorsions(), 2)

        i, j, k, l, m, n, tk = force.getPiTorsionParameters(0)
        self.assertEqual(i, 0)
        self.assertEqual(j, 1)
        self.assertEqual(k, 2)
        self.assertEqual(l, 3)
        self.assertEqual(m, 4)
        self.assertEqual(n, 5)
        self.assertEqual(tk, 1.0*kilojoules_per_mole)
        self.assertIs(tk.unit, kilojoule_per_mole)

        i, j, k, l, m, n, tk = force.getPiTorsionParameters(1)
        self.assertEqual(i, 1)
        self.assertEqual(j, 2)
        self.assertEqual(k, 3)
        self.assertEqual(l, 4)
        self.assertEqual(m, 5)
        self.assertEqual(n, 6)
        self.assertEqual(tk, 1.0*kilocalories_per_mole)
        self.assertIs(tk.unit, kilojoule_per_mole)

    def testAmoebaStretchBendForce(self):
        """ Tests the AmoebaStretchBendForce API features """
        force = AmoebaStretchBendForce()

        force.addStretchBend(0, 1, 2, 0.10, 0.12, math.pi/2, 10.0, 12.0)
        force.addStretchBend(1, 2, 3, 1.0*angstroms, 1.2*angstroms, 60*degrees,
                10.0*kilocalories_per_mole/angstroms/radians,
                12.0*kilocalories_per_mole/angstroms/radians)

        self.assertEqual(force.getNumStretchBends(), 2)

        i, j, k, r1, r2, t, k1, k2 = force.getStretchBendParameters(0)
        self.assertEqual(i, 0)
        self.assertEqual(j, 1)
        self.assertEqual(k, 2)
        self.assertEqual(r1, 0.1*nanometers)
        self.assertIs(r1.unit, nanometers)
        self.assertEqual(r2, 0.12*nanometers)
        self.assertIs(r2.unit, nanometers)
        self.assertEqual(t, math.pi/2*radians)
        self.assertIs(t.unit, radians)
        self.assertEqual(k1, 10*kilojoules_per_mole/nanometers/radians)
        self.assertIs(k1.unit, kilojoules_per_mole/nanometers/radians)
        self.assertEqual(k2, 12*kilojoules_per_mole/nanometers/radians)
        self.assertIs(k2.unit, kilojoules_per_mole/nanometers/radians)

        i, j, k, r1, r2, t, k1, k2 = force.getStretchBendParameters(1)
        self.assertEqual(i, 1)
        self.assertEqual(j, 2)
        self.assertEqual(k, 3)
        self.assertEqual(r1, 1.0*angstroms)
        self.assertIs(r1.unit, nanometers)
        self.assertEqual(r2, 1.2*angstroms)
        self.assertIs(r2.unit, nanometers)
        self.assertAlmostEqualUnit(t, 60*degrees)
        self.assertIs(t.unit, radians)
        self.assertEqual(k1, 10*kilocalories_per_mole/angstroms/radians)
        self.assertIs(k1.unit, kilojoules_per_mole/nanometers/radians)
        self.assertEqual(k2, 12*kilocalories_per_mole/angstroms/radians)
        self.assertIs(k2.unit, kilojoules_per_mole/nanometers/radians)

    def testAmoebaTorsionTorsionForce(self):
        """ Tests the AmoebaTorsionTorsionForce API features """
        force = AmoebaTorsionTorsionForce()

        grid1 = [[[i, j, random.random(), random.random(), random.random(), random.random()]
                    for j in range(-180, 180, 24)] for i in range(-180, 180, 24)]
        kcal = kilocalories_per_mole
        grid2 = [[[i*math.pi/180*radians, j*math.pi/180*radians, random.random()*kcal,
                   random.random()*kcal/radians, random.random()*kcal/radians,
                   random.random()*kcal/radians**2]
            for j in range(-180, 0, 10)] for i in range(-180, 0, 10)]
        grid3 = [[[i, j, random.random()*kcal]
            for j in range(-180, 180, 24)] for i in range(-180, 180, 24)]
        force.setTorsionTorsionGrid(0, grid1)
        force.setTorsionTorsionGrid(1, grid2)
        force.setTorsionTorsionGrid(2, grid3)

        self.assertEqual(force.getNumTorsionTorsionGrids(), 3)

        # Check the grids
        g1 = force.getTorsionTorsionGrid(0)
        for row1, row2 in zip(g1, grid1):
            for column1, column2 in zip(row1, row2):
                self.assertEqual(len(column1), len(column2))
                for x1, x2 in zip(column1, column2):
                    self.assertEqual(x1, x2)
        g2 = force.getTorsionTorsionGrid(1)
        for row1, row2 in zip(g2, grid2):
            for column1, column2 in zip(row1, row2):
                self.assertEqual(column1[0], column2[0].value_in_unit(degree))
                self.assertEqual(column1[1], column2[1].value_in_unit(degree))
                self.assertEqual(column1[2], column2[2].value_in_unit(kilojoules_per_mole))
                self.assertEqual(column1[3], column2[3].value_in_unit(kilojoules_per_mole/radian))
                self.assertEqual(column1[4], column2[4].value_in_unit(kilojoules_per_mole/radian))
                self.assertEqual(column1[5], column2[5].value_in_unit(kilojoules_per_mole/radian**2))
        g3 = force.getTorsionTorsionGrid(2)
        for row1, row2 in zip(g3, grid3):
            for column1, column2 in zip(row1, row2):
                self.assertEqual(len(column1), 6)
                self.assertEqual(len(column2), 3)
                self.assertEqual(column1[0], column2[0])
                self.assertEqual(column1[1], column2[1])
                self.assertEqual(column1[2], column2[2].value_in_unit(kilojoules_per_mole))

        force.addTorsionTorsion(0, 1, 2, 3, 4, 5, 0)
        force.addTorsionTorsion(1, 2, 3, 4, 5, 6, 1)
        force.addTorsionTorsion(2, 3, 4, 5, 6, 7, 2)

        self.assertEqual(force.getNumTorsionTorsions(), 3)

        i, j, k, l, m, ch, g = force.getTorsionTorsionParameters(0)
        self.assertEqual(i, 0)
        self.assertEqual(j, 1)
        self.assertEqual(k, 2)
        self.assertEqual(l, 3)
        self.assertEqual(m, 4)
        self.assertEqual(ch, 5)
        self.assertEqual(g, 0)

        i, j, k, l, m, ch, g = force.getTorsionTorsionParameters(1)
        self.assertEqual(i, 1)
        self.assertEqual(j, 2)
        self.assertEqual(k, 3)
        self.assertEqual(l, 4)
        self.assertEqual(m, 5)
        self.assertEqual(ch, 6)
        self.assertEqual(g, 1)

    def testAmoebaVdwForce(self):
        """ Tests the AmoebaVdwForce API features """
        force = AmoebaVdwForce()

        self.assertEqual(force.getSigmaCombiningRule(), 'CUBIC-MEAN')
        force.setSigmaCombiningRule('ARITHMETIC')
        self.assertEqual(force.getSigmaCombiningRule(), 'ARITHMETIC')
        force.setSigmaCombiningRule('GEOMETRIC')
        self.assertEqual(force.getSigmaCombiningRule(), 'GEOMETRIC')

        self.assertEqual(force.getEpsilonCombiningRule(), 'HHG')
        force.setEpsilonCombiningRule('HARMONIC')
        self.assertEqual(force.getEpsilonCombiningRule(), 'HARMONIC')
        force.setEpsilonCombiningRule('GEOMETRIC')
        self.assertEqual(force.getEpsilonCombiningRule(), 'GEOMETRIC')
        force.setEpsilonCombiningRule('ARITHMETIC')
        self.assertEqual(force.getEpsilonCombiningRule(), 'ARITHMETIC')

        self.assertTrue(force.getUseDispersionCorrection())
        force.setUseDispersionCorrection(False)
        self.assertFalse(force.getUseDispersionCorrection())

        self.assertIs(force.getNonbondedMethod(), AmoebaVdwForce.NoCutoff)
        self.assertFalse(force.usesPeriodicBoundaryConditions())
        force.setNonbondedMethod(AmoebaVdwForce.CutoffPeriodic)
        self.assertTrue(force.usesPeriodicBoundaryConditions())
        self.assertIs(force.getNonbondedMethod(), AmoebaVdwForce.CutoffPeriodic)

        force.setCutoff(10.0*angstroms)
        self.assertEqual(force.getCutoff(), 10.0*angstroms)
        self.assertIs(force.getCutoff().unit, nanometers)

        force.addParticle(0, 0.1, 1.0, 1.0)
        force.addParticle(1, 1.0*angstroms, 1.0*kilocalories_per_mole, 0.5)
        force.addParticle(1, 0.8*angstroms, 2.0*kilocalories_per_mole, 0.25)

        self.assertEqual(force.getNumParticles(), 3)

        p, sig, eps, scale = force.getParticleParameters(0)
        self.assertEqual(p, 0)
        self.assertEqual(sig, 0.1*nanometers)
        self.assertIs(sig.unit, nanometers)
        self.assertEqual(eps, 1.0*kilojoules_per_mole)
        self.assertIs(eps.unit, kilojoules_per_mole)
        self.assertEqual(scale, 1.0)

        p, sig, eps, scale = force.getParticleParameters(1)
        self.assertEqual(p, 1)
        self.assertEqual(sig, 1.0*angstroms)
        self.assertIs(sig.unit, nanometers)
        self.assertEqual(eps, 1.0*kilocalories_per_mole)
        self.assertIs(eps.unit, kilojoules_per_mole)
        self.assertEqual(scale, 0.5)

        p, sig, eps, scale = force.getParticleParameters(2)
        self.assertEqual(p, 1)
        self.assertAlmostEqualUnit(sig, 0.8*angstroms)
        self.assertIs(sig.unit, nanometers)
        self.assertEqual(eps, 2.0*kilocalories_per_mole)
        self.assertIs(eps.unit, kilojoules_per_mole)
        self.assertEqual(scale, 0.25)

    def testAmoebaWcaDispersionForce(self):
        """ Tests the AmoebaWcaDispersionForce API features """
        force = AmoebaWcaDispersionForce()

        self.assertEqual(force.getDispoff(), 0.26*nanometer)
        self.assertEqual(force.getAwater(), 0.033428*nanometer**-3)
        self.assertEqual(force.getEpsh(), 0.0135*kilojoule_per_mole)
        self.assertEqual(force.getEpso(), 0.11*kilojoule_per_mole)
        self.assertEqual(force.getRminh(), 1.3275*nanometer)
        self.assertEqual(force.getRmino(), 1.7025*nanometer)
        self.assertEqual(force.getShctd(), 0.81)
        self.assertEqual(force.getSlevy(), 1.0)

        force.setDispoff(3*angstroms)
        self.assertAlmostEqualUnit(force.getDispoff(), 3*angstroms)
        self.assertIs(force.getDispoff().unit, nanometer)

        force.setAwater(3*angstroms**-3)
        self.assertAlmostEqualUnit(force.getAwater(), 3*angstroms**-3)
        self.assertEqual(1*force.getAwater().unit, 1*nanometer**-3)

        force.setEpsh(1*kilocalorie_per_mole)
        self.assertEqual(force.getEpsh(), 1*kilocalorie_per_mole)
        self.assertIs(force.getEpsh().unit, kilojoule_per_mole)

        force.setEpso(1*kilocalorie_per_mole)
        self.assertEqual(force.getEpso(), 1*kilocalorie_per_mole)
        self.assertIs(force.getEpso().unit, kilojoule_per_mole)

        force.setRminh(20*angstroms)
        self.assertEqual(force.getRminh(), 20*angstroms)
        self.assertIs(force.getRminh().unit, nanometer)

        force.setRmino(30*angstroms)
        self.assertEqual(force.getRmino(), 30*angstroms)
        self.assertIs(force.getRmino().unit, nanometer)

        force.setShctd(1)
        self.assertEqual(force.getShctd(), 1)

        force.setSlevy(2)
        self.assertEqual(force.getSlevy(), 2)

        force.addParticle(0.5, 1)
        force.addParticle(3*angstroms, 1*kilocalorie_per_mole)

        self.assertEqual(force.getNumParticles(), 2)

        sig, eps = force.getParticleParameters(0)
        self.assertEqual(sig, 0.5*nanometer)
        self.assertIs(sig.unit, nanometer)
        self.assertEqual(eps, 1*kilojoule_per_mole)
        self.assertIs(eps.unit, kilojoule_per_mole)

        sig, eps = force.getParticleParameters(1)
        self.assertAlmostEqualUnit(sig, 3*angstrom)
        self.assertIs(sig.unit, nanometer)
        self.assertEqual(eps, 1*kilocalorie_per_mole)
        self.assertIs(eps.unit, kilojoule_per_mole)

    def testAmoebaMultipoleForce(self):
        """ Tests the AmoebaMultipoleForce API features """
        force = AmoebaMultipoleForce()
        self.assertIs(force.getNonbondedMethod(), AmoebaMultipoleForce.NoCutoff)
        self.assertFalse(force.usesPeriodicBoundaryConditions())
        force.setNonbondedMethod(AmoebaMultipoleForce.PME)
        self.assertIs(force.getNonbondedMethod(), AmoebaMultipoleForce.PME)
        self.assertTrue(force.usesPeriodicBoundaryConditions())

        self.assertEqual(force.getCutoffDistance(), 1*nanometer)
        self.assertIs(force.getCutoffDistance().unit, nanometer)
        force.setCutoffDistance(8*angstroms)
        self.assertEqual(force.getCutoffDistance(), 8*angstrom)
        self.assertIs(force.getCutoffDistance().unit, nanometer)

        self.assertEqual(force.getPolarizationType(), AmoebaMultipoleForce.Mutual)
        force.setPolarizationType(AmoebaMultipoleForce.Direct)
        self.assertEqual(force.getPolarizationType(), AmoebaMultipoleForce.Direct)
        force.setPolarizationType(AmoebaMultipoleForce.Mutual)
        self.assertEqual(force.getPolarizationType(), AmoebaMultipoleForce.Mutual)

        force.addMultipole(1.0, [0.5, 0, -0.5], list(range(9)),
                AmoebaMultipoleForce.ZThenX, 1, 2, 3, 0.5, 0.5, 1.0)
        force.addMultipole(1.0*elementary_charge, [0.5, 0, -0.5]*elementary_charge*angstrom,
                list(range(9))*elementary_charge*angstrom**2,
                AmoebaMultipoleForce.Bisector, 2, 3, 4, 0.5, 0.5, 1.0*angstrom**3)

        self.assertEqual(force.getNumMultipoles(), 2)

        q, mu, quad, ax, i, j, k, thole, damp, polarity = force.getMultipoleParameters(0)

        self.assertEqual(q, 1.0*elementary_charge)
        self.assertEqual(mu, (0.5, 0, -0.5)*elementary_charge*nanometer)
        self.assertEqual(quad, tuple(range(9))*elementary_charge*nanometer**2)
        self.assertEqual(ax, AmoebaMultipoleForce.ZThenX)
        self.assertEqual(i, 1)
        self.assertEqual(j, 2)
        self.assertEqual(k, 3)
        self.assertEqual(thole, 0.5)
        self.assertEqual(damp, 0.5)
        self.assertEqual(polarity, 1*nanometer**3)

        q, mu, quad, ax, i, j, k, thole, damp, polarity = force.getMultipoleParameters(1)

        self.assertEqual(q, 1.0*elementary_charge)
        self.assertEqual(mu, (0.5, 0, -0.5)*elementary_charge*angstrom)
        self.assertAlmostEqualUnitArray(quad, tuple(range(9))*elementary_charge*angstrom**2)
        self.assertEqual(ax, AmoebaMultipoleForce.Bisector)
        self.assertEqual(i, 2)
        self.assertEqual(j, 3)
        self.assertEqual(k, 4)
        self.assertEqual(thole, 0.5)
        self.assertEqual(damp, 0.5)
        self.assertAlmostEqualUnit(polarity, 1*angstrom**3)

    def testVerletIntegrator(self):
        """ Tests the VerletIntegrator API features """
        integrator = VerletIntegrator(1.0)
        self.assertEqual(integrator.getStepSize(), 1.0*picosecond)
        self.assertEqual(integrator.getConstraintTolerance(), 1e-5)
        integrator.setConstraintTolerance(1e-6)
        self.assertEqual(integrator.getConstraintTolerance(), 1e-6)

        integrator = VerletIntegrator(1.0*femtoseconds)
        self.assertEqual(integrator.getStepSize(), 1.0*femtoseconds)

    def testVariableVerletIntegrator(self):
        """ Tests the VariableVerletIntegrator API features """
        integrator = VariableVerletIntegrator(0.1)
        self.assertEqual(integrator.getErrorTolerance(), 0.1)
        integrator.setErrorTolerance(0.01)
        self.assertEqual(integrator.getErrorTolerance(), 0.01)
        self.assertEqual(integrator.getConstraintTolerance(), 1e-5)
        integrator.setConstraintTolerance(1e-6)
        self.assertEqual(integrator.getConstraintTolerance(), 1e-6)

    def testLangevinIntegrator(self):
        """ Tests the LangevinIntegrator API features """
        integrator = LangevinIntegrator(300, 0.1, 1)
        self.assertEqual(integrator.getTemperature(), 300*kelvin)
        self.assertEqual(integrator.getFriction(), 0.1/picosecond)
        self.assertEqual(integrator.getStepSize(), 1*picosecond)

        integrator = LangevinIntegrator(300*kelvin, 0.1/microsecond, 1*femtosecond)
        self.assertEqual(integrator.getTemperature(), 300*kelvin)
        self.assertAlmostEqualUnit(integrator.getFriction(), 0.1/microsecond)
        self.assertEqual(integrator.getStepSize(), 1*femtosecond)

    def testVariableLangevinIntegrator(self):
        """ Tests the VariableLangevinIntegrator API features """
        integrator = VariableLangevinIntegrator(300, 0.1, 0.1)
        self.assertEqual(integrator.getTemperature(), 300*kelvin)
        self.assertEqual(integrator.getFriction(), 0.1/picosecond)
        self.assertEqual(integrator.getErrorTolerance(), 0.1)

        integrator = VariableLangevinIntegrator(300*kelvin, 0.1/microsecond, 0.01)
        self.assertEqual(integrator.getTemperature(), 300*kelvin)
        self.assertAlmostEqualUnit(integrator.getFriction(), 0.1/microsecond)
        self.assertEqual(integrator.getErrorTolerance(), 0.01)

    def testBrownianIntegrator(self):
        """ Tests the BrownianIntegrator API features """
        integrator = BrownianIntegrator(300, 0.1, 1)
        self.assertEqual(integrator.getTemperature(), 300*kelvin)
        self.assertEqual(integrator.getFriction(), 0.1/picosecond)
        self.assertEqual(integrator.getStepSize(), 1*picosecond)

        integrator = BrownianIntegrator(300*kelvin, 0.1/microsecond, 1*femtosecond)
        self.assertEqual(integrator.getTemperature(), 300*kelvin)
        self.assertAlmostEqualUnit(integrator.getFriction(), 0.1/microsecond)
        self.assertEqual(integrator.getStepSize(), 1*femtosecond)

    def testAndersenThermostat(self):
        """ Tests the AndersenThermostat API features """
        force = AndersenThermostat(300*kelvin, 1/microsecond)
        self.assertEqual(force.getDefaultTemperature(), 300*kelvin)
        self.assertAlmostEqualUnit(force.getDefaultCollisionFrequency(), 1/microsecond)

        force.setDefaultTemperature(298.15)
        force.setDefaultCollisionFrequency(1)
        self.assertEqual(force.getDefaultTemperature(), 298.15*kelvin)
        self.assertAlmostEqualUnit(force.getDefaultCollisionFrequency(), 1/picosecond)

    def testMonteCarloMembraneBarostat(self):
        """ Tests the MonteCarloMembraneBarostat API features """
        force = MonteCarloMembraneBarostat(1.0, 1.5, 300, MonteCarloMembraneBarostat.XYAnisotropic, MonteCarloMembraneBarostat.ZFixed, 25)
        self.assertEqual(force.getDefaultPressure(), 1.0*bar)
        self.assertEqual(force.getDefaultSurfaceTension(), 1.5*bar*nanometer)
        self.assertEqual(force.getDefaultTemperature(), 300*kelvin)
        self.assertEqual(force.getXYMode(), MonteCarloMembraneBarostat.XYAnisotropic)
        self.assertEqual(force.getZMode(), MonteCarloMembraneBarostat.ZFixed)
        self.assertEqual(force.getFrequency(), 25)

        force = MonteCarloMembraneBarostat(1.1*bar, 2.0*bar*nanometer, 350*kelvin, MonteCarloMembraneBarostat.XYAnisotropic, MonteCarloMembraneBarostat.ZFixed, 25)
        self.assertEqual(force.getDefaultPressure(), 1.1*bar)
        self.assertEqual(force.getDefaultSurfaceTension(), 2.0*bar*nanometer)
        self.assertEqual(force.getDefaultTemperature(), 350*kelvin)

        force.setDefaultPressure(1.2*bar)
        force.setDefaultSurfaceTension(2.5*bar*nanometer)
        force.setDefaultTemperature(298.15)
        self.assertEqual(force.getDefaultPressure(), 1.2*bar)
        self.assertEqual(force.getDefaultSurfaceTension(), 2.5*bar*nanometer)
        self.assertEqual(force.getDefaultTemperature(), 298.15*kelvin)

    def testDrudeSCFIntegrator(self):
        """ Tests the DrudeSCFIntegrator API features """
        integrator = DrudeSCFIntegrator(0.002)
        self.assertEqual(integrator.getStepSize(), 0.002*picoseconds)
        self.assertAlmostEqualUnit(integrator.getMinimizationErrorTolerance(),
                                   0.1*kilojoule_per_mole/nanometer)
        integrator.setMinimizationErrorTolerance(1*kilocalorie_per_mole/angstrom)
        self.assertAlmostEqualUnit(integrator.getMinimizationErrorTolerance(),
                                   1*kilocalorie_per_mole/angstrom)

    def testDrudeLangevinIntegrator(self):
        """ Tests the DrudeLangevinIntegrator API features """
        integrator = DrudeLangevinIntegrator(300*kelvin, 1.0/microsecond,
                                             100*kelvin, 1/picosecond, 0.1*femtosecond)
        self.assertEqual(integrator.getTemperature(), 300*kelvin)
        integrator.setTemperature(298.15)
        self.assertEqual(integrator.getTemperature(), 298.15*kelvin)
        self.assertAlmostEqualUnit(integrator.getFriction(), 1/microsecond)
        integrator.setFriction(1)
        self.assertEqual(integrator.getFriction(), 1/picosecond)
        self.assertEqual(integrator.getDrudeTemperature(), 100*kelvin)
        integrator.setDrudeTemperature(50)
        self.assertEqual(integrator.getDrudeTemperature(), 50*kelvin)
        self.assertEqual(integrator.getStepSize(), 0.1*femtosecond)
        integrator.setStepSize(0.0005)
        self.assertEqual(integrator.getStepSize(), 0.0005*picosecond)
        self.assertEqual(integrator.getMaxDrudeDistance(), 0*nanometer)
        integrator.setMaxDrudeDistance(0.05)
        self.assertEqual(integrator.getMaxDrudeDistance(), 0.05*nanometer)

    def testCustomIntegrator(self):
        """ Tests the CustomIntegrator API features """
        integrator = CustomIntegrator(1.0)
        self.assertEqual(integrator.getStepSize(), 1.0*picosecond)
        self.assertEqual(integrator.getConstraintTolerance(), 1e-5)
        integrator.setConstraintTolerance(1e-6)
        self.assertEqual(integrator.getConstraintTolerance(), 1e-6)

        integrator = CustomIntegrator(1.0*femtoseconds)
        self.assertEqual(integrator.getStepSize(), 1.0*femtoseconds)

if __name__ == '__main__':
    unittest.main()
