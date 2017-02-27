import unittest
from validateConstraints import *
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import simtk.openmm.app.element as elem
import simtk.openmm.app.forcefield as forcefield
import math
try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO
import os
import warnings

class TestGenerators(unittest.TestCase):
    """Test the various generators found in forcefield.py."""

    def setUp(self):
        # alanine dipeptide with implicit water
        self.pdb1 = PDBFile('systems/alanine-dipeptide-implicit.pdb')


    def test_CustomHbondGenerator(self):
        """Test the generator for CustomHbondForce."""

        for bondCutoff in range(4):
            xml = """
<ForceField>
 <CustomHbondForce energy="a*b*distance(a1,d1)" particlesPerDonor="3" particlesPerAcceptor="2" bondCutoff="%d">
  <PerDonorParameter name="a"/>
  <PerAcceptorParameter name="b"/>
  <Donor class1="C" class2="N" class3="H" a="3"/>
  <Acceptor class1="C" class2="O" b="2"/>
  <Function name="test" min="1" max="2" type="Continuous1D">
   0 1 2 3 4 5
  </Function>
 </CustomHbondForce>
</ForceField>""" % bondCutoff
            ff = ForceField('amber99sb.xml', StringIO(xml))
            system = ff.createSystem(self.pdb1.topology)
            hbond = [f for f in system.getForces() if isinstance(f, CustomHbondForce)][0]
            self.assertEqual(1, hbond.getNumPerDonorParameters())
            self.assertEqual(1, hbond.getNumPerAcceptorParameters())
            self.assertEqual('a', hbond.getPerDonorParameterName(0))
            self.assertEqual('b', hbond.getPerAcceptorParameterName(0))
            self.assertEqual(1, hbond.getNumTabulatedFunctions())
            expectedDonors = [(4,6,7), (14,16,17)]
            expectedAcceptors = [(4,5,-1), (14,15,-1)]
            self.assertEqual(len(expectedDonors), hbond.getNumDonors())
            self.assertEqual(len(expectedAcceptors), hbond.getNumAcceptors())
            for i in range(hbond.getNumDonors()):
                atom1, atom2, atom3, params = hbond.getDonorParameters(i)
                self.assertTrue((atom1, atom2, atom3) in expectedDonors)
                self.assertEqual((3.0,), params)
            for i in range(hbond.getNumAcceptors()):
                atom1, atom2, atom3, params = hbond.getAcceptorParameters(i)
                self.assertTrue((atom1, atom2, atom3) in expectedAcceptors)
                self.assertEqual((2.0,), params)
            expectedExclusions = [(0,0), (1,1)]
            if bondCutoff >= 2:
                expectedExclusions.append((0,1))
            if bondCutoff >= 3:
                expectedExclusions.append((1,0))
            self.assertEqual(len(expectedExclusions), hbond.getNumExclusions())
            for i in range(hbond.getNumExclusions()):
                self.assertTrue(tuple(hbond.getExclusionParticles(i)) in expectedExclusions)


    def test_CustomHbondGenerator2(self):
        """Test the generator for CustomHbondForce with different parameters."""

        xml = """
<ForceField>
 <CustomHbondForce energy="a*b*distance(a1,d1)" particlesPerDonor="2" particlesPerAcceptor="1" bondCutoff="0">
  <PerDonorParameter name="a"/>
  <PerAcceptorParameter name="b"/>
  <Donor class1="N" class2="H" a="3"/>
  <Acceptor class1="O" b="2"/>
 </CustomHbondForce>
</ForceField>"""
        ff = ForceField('amber99sb.xml', StringIO(xml))
        system = ff.createSystem(self.pdb1.topology)
        hbond = [f for f in system.getForces() if isinstance(f, CustomHbondForce)][0]
        self.assertEqual(1, hbond.getNumPerDonorParameters())
        self.assertEqual(1, hbond.getNumPerAcceptorParameters())
        self.assertEqual('a', hbond.getPerDonorParameterName(0))
        self.assertEqual('b', hbond.getPerAcceptorParameterName(0))
        expectedDonors = [(6,7,-1), (16,17,-1)]
        expectedAcceptors = [(5,-1,-1), (15,-1,-1)]
        self.assertEqual(len(expectedDonors), hbond.getNumDonors())
        self.assertEqual(len(expectedAcceptors), hbond.getNumAcceptors())
        for i in range(hbond.getNumDonors()):
            atom1, atom2, atom3, params = hbond.getDonorParameters(i)
            self.assertTrue((atom1, atom2, atom3) in expectedDonors)
            self.assertEqual((3.0,), params)
        for i in range(hbond.getNumAcceptors()):
            atom1, atom2, atom3, params = hbond.getAcceptorParameters(i)
            self.assertTrue((atom1, atom2, atom3) in expectedAcceptors)
            self.assertEqual((2.0,), params)


if __name__ == '__main__':
    unittest.main()
