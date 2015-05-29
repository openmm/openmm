import unittest
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import simtk.openmm.app.element as elem

class TestGromacsGroFile(unittest.TestCase):
    """Test the Gromacs GRO file parser"""
 
    def test_Triclinic(self):
        """Test parsing a file that describes a triclinic box."""
        gro = GromacsGroFile('systems/triclinic.gro')
        self.assertEqual(len(gro.positions), 8)
        expectedPositions = [
            Vec3(1.744, 2.788, 3.162),
            Vec3(1.048, 0.762, 2.340),
            Vec3(2.489, 1.570, 2.817),
            Vec3(1.027, 1.893, 3.271),
            Vec3(0.937, 0.825, 0.009),
            Vec3(2.290, 1.887, 3.352),
            Vec3(1.266, 1.111, 2.894),
            Vec3(0.933, 1.862, 3.490)]*nanometers
        for (p1, p2) in zip(expectedPositions, gro.positions):
            self.assertEqual(p1, p2)
        expectedVectors = [
            Vec3(2.5, 0, 0),
            Vec3(0.5, 3.0, 0),
            Vec3(0.7, 0.9, 3.5)]*nanometers
        for (v1, v2) in zip(expectedVectors, gro.getPeriodicBoxVectors()):
            self.assertEqual(v1, v2)
        self.assertEqual(Vec3(2.5, 3.0, 3.5)*nanometers, gro.getUnitCellDimensions())
        for i in range(4):
            self.assertEqual(elem.chlorine, gro.elements[i])
            self.assertEqual('Cl', gro.atomNames[i])
            self.assertEqual('Cl', gro.residueNames[i])
        for i in range(4, 8):
            self.assertEqual(elem.sodium, gro.elements[i])
            self.assertEqual('Na', gro.atomNames[i])
            self.assertEqual('Na', gro.residueNames[i])

if __name__ == '__main__':
    unittest.main()
