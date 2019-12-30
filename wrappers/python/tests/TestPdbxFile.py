import unittest
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import simtk.openmm.app.element as elem
import os
if sys.version_info >= (3, 0):
    from io import StringIO
else:
    from cStringIO import StringIO

class TestPdbxFile(unittest.TestCase):
    """Test the PDBx/mmCIF file parser"""

    def test_FormatConversion(self):
        """Test conversion from PDB to PDBx"""

        mol = PDBFile('systems/ala_ala_ala.pdb')

        # Write to 'file'
        output = StringIO()
        PDBxFile.writeFile(mol.topology, mol.positions, output,
                           keepIds=True)

        # Read from 'file'
        input = StringIO(output.getvalue())
        try:
            pdbx = PDBxFile(input)
        except Exception:
            self.fail('Parser failed to read PDBx/mmCIF file')

        # Close file handles
        output.close()
        input.close()


    def test_Triclinic(self):
        """Test parsing a file that describes a triclinic box."""
        pdb = PDBxFile('systems/triclinic.pdbx')
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
            self.assertVecAlmostEqual(v1, v2, 1e-6)
        self.assertVecAlmostEqual(Vec3(2.5, 3.0, 3.5)*nanometers, pdb.topology.getUnitCellDimensions(), 1e-6)
        for atom in pdb.topology.atoms():
            if atom.index < 4:
                self.assertEqual(elem.chlorine, atom.element)
                self.assertEqual('Cl', atom.name)
                self.assertEqual('Cl', atom.residue.name)
            else:
                self.assertEqual(elem.sodium, atom.element)
                self.assertEqual('Na', atom.name)
                self.assertEqual('Na', atom.residue.name)

    def assertVecAlmostEqual(self, p1, p2, tol=1e-7):
        unit = p1.unit
        p1 = p1.value_in_unit(unit)
        p2 = p2.value_in_unit(unit)
        scale = max(1.0, norm(p1),)
        for i in range(3):
            diff = abs(p1[i]-p2[i])/scale
            self.assertTrue(diff < tol)

    def testReporterImplicit(self):
        """ Tests the PDBxReporter without PBC """
        parm = AmberPrmtopFile('systems/alanine-dipeptide-implicit.prmtop')
        system = parm.createSystem()
        sim = Simulation(parm.topology, system, VerletIntegrator(1*femtoseconds),
                         Platform.getPlatformByName('Reference'))
        sim.context.setPositions(PDBFile('systems/alanine-dipeptide-implicit.pdb').getPositions())
        sim.reporters.append(PDBxReporter('test.cif', 1))
        sim.step(10)
        pdb = PDBxFile('test.cif')
        self.assertEqual(len(list(pdb.topology.atoms())), len(list(parm.topology.atoms())))
        self.assertEqual(len(list(pdb.topology.residues())), len(list(parm.topology.residues())))
        for res1, res2 in zip(pdb.topology.residues(), parm.topology.residues()):
            self.assertEqual(res1.name, res2.name)
            for atom1, atom2 in zip(res1.atoms(), res2.atoms()):
                self.assertEqual(atom1.name, atom2.name)
        positions = pdb.getPositions(frame=9)
        self.assertFalse(all(x1 == x2 for x1, x2 in zip(positions, pdb.getPositions(frame=0))))
        # There should only be 10 frames (0 through 9)
        self.assertRaises(IndexError, lambda: pdb.getPositions(frame=10))
        self.assertIs(pdb.topology.getPeriodicBoxVectors(), None)
        del sim
        os.unlink('test.cif')

    def assertAlmostEqualVec(self, vec1, vec2, *args, **kwargs):
        if is_quantity(vec1):
            vec1 = vec1.value_in_unit_system(md_unit_system)
        if is_quantity(vec2):
            vec2 = vec2.value_in_unit_system(md_unit_system)
        for x, y in zip(vec1, vec2):
            self.assertAlmostEqual(x, y, *args, **kwargs)

    def testReporterExplicit(self):
        """ Tests the PDBxReporter with PBC """
        parm = AmberPrmtopFile('systems/alanine-dipeptide-explicit.prmtop')
        system = parm.createSystem(nonbondedCutoff=1.0, nonbondedMethod=PME)
        sim = Simulation(parm.topology, system, VerletIntegrator(1*femtoseconds),
                         Platform.getPlatformByName('Reference'))
        orig_pdb = PDBFile('systems/alanine-dipeptide-explicit.pdb')
        sim.context.setPositions(orig_pdb.getPositions())
        sim.context.setPeriodicBoxVectors(*parm.topology.getPeriodicBoxVectors())
        sim.reporters.append(PDBxReporter('test.cif', 1))
        sim.step(10)
        pdb = PDBxFile('test.cif')
        self.assertEqual(len(list(pdb.topology.atoms())), len(list(parm.topology.atoms())))
        self.assertEqual(len(list(pdb.topology.residues())), len(list(parm.topology.residues())))
        for res1, res2 in zip(pdb.topology.residues(), parm.topology.residues()):
            self.assertEqual(res1.name, res2.name)
            for atom1, atom2 in zip(res1.atoms(), res2.atoms()):
                self.assertEqual(atom1.name, atom2.name)
        positions = pdb.getPositions(frame=9)
        self.assertFalse(all(x1 == x2 for x1, x2 in zip(positions, pdb.getPositions(frame=0))))
        # There should only be 10 frames (0 through 9)
        self.assertRaises(IndexError, lambda: pdb.getPositions(frame=10))
        self.assertAlmostEqualVec(parm.topology.getPeriodicBoxVectors()[0],
                                  pdb.topology.getPeriodicBoxVectors()[0],
                                  places=5)
        self.assertAlmostEqualVec(parm.topology.getPeriodicBoxVectors()[1],
                                  pdb.topology.getPeriodicBoxVectors()[1],
                                  places=5)
        self.assertAlmostEqualVec(parm.topology.getPeriodicBoxVectors()[2],
                                  pdb.topology.getPeriodicBoxVectors()[2],
                                  places=5)
        del sim
        os.unlink('test.cif')

    def testBonds(self):
        """Test reading and writing a file that includes bonds."""
        pdb = PDBFile('systems/methanol_ions.pdb')
        output = StringIO()
        PDBxFile.writeFile(pdb.topology, pdb.positions, output)
        input = StringIO(output.getvalue())
        pdbx = PDBxFile(input)
        output.close()
        input.close()
        self.assertEqual(pdb.topology.getNumBonds(), pdbx.topology.getNumBonds())
        for bond1, bond2 in zip(pdb.topology.bonds(), pdbx.topology.bonds()):
            self.assertEqual(bond1[0].name, bond2[0].name)
            self.assertEqual(bond1[1].name, bond2[1].name)

    def testMultiChain(self):
        """Test reading and writing a file that includes multiple chains"""
        cif_ori = PDBxFile('systems/multichain.pdbx')

        output = StringIO()
        PDBxFile.writeFile(cif_ori.topology, cif_ori.positions, output, keepIds=True)
        input = StringIO(output.getvalue())
        cif_new = PDBxFile(input)
        output.close()
        input.close()

        self.assertEqual(cif_ori.topology.getNumChains(), cif_new.topology.getNumChains())

        for chain1, chain2 in zip(cif_ori.topology.chains(), cif_new.topology.chains()):
            self.assertEqual(chain1.id, chain2.id)

    def testInsertionCodes(self):
        """Test reading a file that uses insertion codes."""
        pdbx = PDBxFile('systems/insertions.pdbx')
        residues = list(pdbx.topology.residues())
        self.assertEqual(7, len(residues))
        names = ['PHE', 'ASP', 'LYS', 'ILE', 'LYS', 'ASN', 'TRP']
        ids = ['59', '60', '60', '60', '60', '60', '61']
        codes = ['', '', 'A', 'B', 'C', 'D', '']
        for res, name, id, code in zip(residues, names, ids, codes):
            self.assertEqual(name, res.name)
            self.assertEqual(id, res.id)
            self.assertEqual(code, res.insertionCode)

if __name__ == '__main__':
    unittest.main()
