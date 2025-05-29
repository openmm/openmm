import tempfile
import unittest
from openmm.app import *
from openmm import *
from openmm.unit import *
import openmm.app.element as elem
from io import StringIO


class TestPdbFile(unittest.TestCase):
    """Test the PDB file parser"""

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
        def compareFiles(pdb1, pdb2):
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

        pdb1 = PDBFile('systems/triclinic.pdb')

        # First try writing to an open file object.

        output = StringIO()
        PDBFile.writeFile(pdb1.topology, pdb1.positions, output)
        input = StringIO(output.getvalue())
        pdb2 = PDBFile(input)
        output.close()
        input.close()
        compareFiles(pdb1, pdb2)

        # Now try a filename.

        with tempfile.TemporaryDirectory() as tempdir:
            filename = os.path.join(tempdir, 'temp.pdb')
            PDBFile.writeFile(pdb1.topology, pdb1.positions, filename)
            pdb2 = PDBFile(filename)
            compareFiles(pdb1, pdb2)

    def test_BinaryStream(self):
        """Test reading a stream that was opened in binary mode."""
        with open('systems/triclinic.pdb', 'rb') as infile:
            pdb = PDBFile(infile)
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
        pdb = PDBFile(input)
        for atom in pdb.topology.atoms():
            if atom.index > 2:
                self.assertEqual(None, atom.element)

    def test_AltLocs(self):
        """Test reading a file that includes AltLocs"""
        for filename in ['altlocs.pdb', 'altlocs2.pdb']:
            pdb = PDBFile(f'systems/{filename}')
            self.assertEqual(1, pdb.topology.getNumResidues())
            self.assertEqual(19, pdb.topology.getNumAtoms())
            self.assertEqual(19, len(pdb.positions))
            self.assertEqual('ILE', list(pdb.topology.residues())[0].name)

    def test_FormalCharges(self):
        """Test reading, and writing and re-reading of a file containing formal charges."""
        pdb = PDBFile('systems/formal-charges.pdb')
        for atom in pdb.topology.atoms():
            if atom.index == 8:
                self.assertEqual(+1, atom.formalCharge)
            elif atom.index == 29:
                self.assertEqual(-1, atom.formalCharge)
            else:
                self.assertEqual(None, atom.formalCharge)
        output = StringIO()
        PDBFile.writeFile(pdb.topology, pdb.positions, output)
        input = StringIO(output.getvalue())
        pdb = PDBFile(input)
        for atom in pdb.topology.atoms():
            if atom.index == 8:
                self.assertEqual(+1, atom.formalCharge)
            elif atom.index == 29:
                self.assertEqual(-1, atom.formalCharge)
            else:
                self.assertEqual(None, atom.formalCharge)

    def test_LargeFile(self):
        """Write and read a file with more than 100,000 atoms"""
        topology = Topology()
        chain = topology.addChain('A')
        for i in range(20000):
            res = topology.addResidue('MOL', chain)
            atoms = []
            for j in range(6):
                atoms.append(topology.addAtom(f'AT{j}', elem.carbon, res))
            for j in range(5):
                topology.addBond(atoms[j], atoms[j+1])
        positions = [Vec3(0, 0, 0)]*topology.getNumAtoms()

        # The model has 20,000 residues and 120,000 atoms.

        output = StringIO()
        PDBFile.writeFile(topology, positions, output)
        input = StringIO(output.getvalue())
        pdb = PDBFile(input)
        output.close()
        input.close()
        self.assertEqual(len(positions), len(pdb.positions))
        self.assertEqual(topology.getNumAtoms(), pdb.topology.getNumAtoms())
        self.assertEqual(topology.getNumResidues(), pdb.topology.getNumResidues())
        self.assertEqual(topology.getNumChains(), pdb.topology.getNumChains())
        self.assertEqual(topology.getNumBonds(), pdb.topology.getNumBonds())
        for atom in pdb.topology.atoms():
            self.assertEqual(str(atom.index+1), atom.id)
        for res in pdb.topology.residues():
            self.assertEqual(str(res.index+1), res.id)

        # Make sure the CONECT records were interpreted correctly.

        bonds = set()
        for atom1, atom2 in topology.bonds():
            bonds.add(tuple(sorted((atom1.index, atom2.index))))
        for atom1, atom2 in pdb.topology.bonds():
            assert tuple(sorted((atom1.index, atom2.index))) in bonds

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
