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


class TestPdbFile(unittest.TestCase):
    """Test the PDB file parser"""

    def build_bond_dict(topology):
        bonds = dict()
        for (atom1,atom2) in pdb.topology.bonds():
            if atom1.id in bonds:
                bonds[atom1.id].add(atom2.id)
            else:
                bonds[atom1.id] = set([atom2.id])

            if atom2.id in bonds:
                bonds[atom2.id].add(atom1.id)
            else:
                bonds[atom2.id] = set([atom1.id])

    def test_ReadFile_standardBonds(self):
        """Test parsing a file to see if standard bonds are correctly applied."""

        # phosphorylated p38 kinase
        pdb = PDBFile('systems/3py3.pdb')
        bonds = build_bond_dict(pdb.topology)
        # Check peptide bonds
        assert('3' in bonds['5'])
        assert('7' in bonds['11'])
        assert('1473' in bonds['1479']) # MET179-TPO180
        assert('1488' in bonds['1490']) # TPO180-GLY181
        # Check bonds in TPO180 residue
        assert('1479' in bonds['1480']) # N-CA
        assert('1480' in bonds['1488']) # CA-C
        assert('1484' in bonds['1485']) # P-O1P
        assert('1484' in bonds['1486']) # P-O2P
        assert('1484' in bonds['1487']) # P-O3P

        # metal-binding protein
        pdb = PDBFile('systems/1T2Y.pdb')
        bonds = build_bond_dict(pdb.topology)
        # Check N-terminus
        assert('1' in bonds['5']) # N-H1
        assert('1' in bonds['6']) # N-H2
        assert('1' in bonds['7']) # N-H3
        # Check C-terminus
        assert('251' in bonds['252']) # C-O
        assert('251' in bonds['258']) # C-OXT

    def test_Triclinic(self):
        """Test parsing a file that describes a triclinic box."""
        pdb = PDBFile('systems/triclinic.pdb')
        self.assertEqual(len(pdb.positions), 8)
        expectedPositions = [
            Vec3(1.744, 2.788, 3.162),
            Vec3(1.048, 0.762, 2.340),
            Vec3(2.489, 1.570, 2.817),
            Vec3(1.027, 1.893, 3.271),
            Vec3(0.937, 0.825, 0.009),
            Vec3(2.290, 1.887, 3.352),
            Vec3(1.266, 1.111, 2.894),
            Vec3(0.933, 1.862, 3.490)]*nanometers
        for (p1, p2) in zip(expectedPositions, pdb.positions):
            self.assertVecAlmostEqual(p1, p2)
        expectedVectors = [
            Vec3(2.5, 0, 0),
            Vec3(0.5, 3.0, 0),
            Vec3(0.7, 0.9, 3.5)]*nanometers
        for (v1, v2) in zip(expectedVectors, pdb.topology.getPeriodicBoxVectors()):
            self.assertVecAlmostEqual(v1, v2, 1e-4)
        self.assertVecAlmostEqual(Vec3(2.5, 3.0, 3.5)*nanometers, pdb.topology.getUnitCellDimensions(), 1e-4)
        for atom in pdb.topology.atoms():
            if atom.index < 4:
                self.assertEqual(elem.chlorine, atom.element)
                self.assertEqual('Cl', atom.name)
                self.assertEqual('Cl', atom.residue.name)
            else:
                self.assertEqual(elem.sodium, atom.element)
                self.assertEqual('Na', atom.name)
                self.assertEqual('Na', atom.residue.name)

    def test_WriteFile(self):
        """Write a file, read it back, and make sure it matches the original."""
        pdb1 = PDBFile('systems/triclinic.pdb')
        output = StringIO()
        PDBFile.writeFile(pdb1.topology, pdb1.positions, output)
        input = StringIO(output.getvalue())
        pdb2 = PDBFile(input)
        output.close();
        input.close();
        self.assertEqual(len(pdb2.positions), 8)
        for (p1, p2) in zip(pdb1.positions, pdb2.positions):
            self.assertVecAlmostEqual(p1, p2)
        for (v1, v2) in zip(pdb1.topology.getPeriodicBoxVectors(), pdb2.topology.getPeriodicBoxVectors()):
            self.assertVecAlmostEqual(v1, v2, 1e-4)
        self.assertVecAlmostEqual(pdb1.topology.getUnitCellDimensions(), pdb2.topology.getUnitCellDimensions(), 1e-4)
        for atom1, atom2 in zip(pdb1.topology.atoms(), pdb2.topology.atoms()):
            self.assertEqual(atom1.element, atom2.element)
            self.assertEqual(atom1.name, atom2.name)
            self.assertEqual(atom1.residue.name, atom2.residue.name)

    def test_BinaryStream(self):
        """Test reading a stream that was opened in binary mode."""
        pdb = PDBFile(open('systems/triclinic.pdb', 'rb'))
        self.assertEqual(len(pdb.positions), 8)

    def test_ExtraParticles(self):
        """Test reading, and writing and re-reading of a file containing extra particle atoms."""
        pdb = PDBFile('systems/tip5p.pdb')
        for atom in pdb.topology.atoms():
            if atom.index > 2:
                self.assertEqual(None, atom.element)
        output = StringIO()
        PDBFile.writeFile(pdb.topology, pdb.positions, output)
        input = StringIO(output.getvalue())
        pdb = PDBFile(input, extraParticleIdentifier = '')
        for atom in pdb.topology.atoms():
            if atom.index > 2:
                self.assertEqual(None, atom.element)


    def assertVecAlmostEqual(self, p1, p2, tol=1e-7):
        unit = p1.unit
        p1 = p1.value_in_unit(unit)
        p2 = p2.value_in_unit(unit)
        scale = max(1.0, norm(p1),)
        for i in range(3):
            diff = abs(p1[i]-p2[i])/scale
            self.assertTrue(diff < tol)


if __name__ == '__main__':
    unittest.main()
