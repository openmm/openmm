import unittest
import tempfile
from openmm import app
import openmm as mm
from openmm import unit
from openmm.unit.unit_math import norm
import os
import gc

def assertVecAlmostEqual(p1, p2, tol=1e-7):
    unit = p1.unit
    p1 = p1.value_in_unit(unit)
    p2 = p2.value_in_unit(unit)
    scale = max(1.0, norm(p1),)
    for i in range(3):
        diff = abs(p1[i]-p2[i])/scale
        assert(diff < tol)


class TestPDBReporter(unittest.TestCase):
    def setUp(self):
        self.pdb = app.PDBFile('systems/alanine-dipeptide-explicit.pdb')
        self.forcefield = app.ForceField('amber99sbildn.xml','tip3p.xml')
        self.system = self.forcefield.createSystem(self.pdb.topology, nonbondedMethod=app.CutoffNonPeriodic, constraints=app.HBonds)

    def testSubset(self):
        """Test writing out a subset of atoms"""
    
        with tempfile.TemporaryDirectory() as tempdir:
            filename = os.path.join(tempdir, 'temptraj.pdb')
            simulation = app.Simulation(self.pdb.topology, self.system, mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds))
            simulation.context.setPositions(self.pdb.positions)

            # just the alanine-dipeptide atoms
            subset = list(range(0,22))

            simulation.reporters.append(app.PDBReporter(filename, 1, atomSubset=subset))
            simulation.step(1)

            # clear reporters to ensure PDBReporter calls writeFooter and file.close
            simulation.reporters.clear()
            
            # check if the output out trajectory contains the expected number of atoms
            checkpdb = app.PDBFile(filename)
            self.assertEqual(len(checkpdb.positions), len(subset))

            # check the positions are correct
            validPositions  = [simulation.context.getState(getPositions=True).getPositions()[i] for i in subset]
            # round to 4 decimal places before comparing
            for i in range(len(validPositions)):
                validPositions[i] = [round(validPositions[i][j]._value, 4) for j in (0, 1, 2)]*unit.nanometer
                
            for (p1, p2) in zip(checkpdb.positions, validPositions):
                assertVecAlmostEqual(p1, p2)

            # check elements and residue names are correct
            validAtoms = [list(self.pdb.topology.atoms())[i] for i in subset]
            for atom1, atom2 in zip(checkpdb.topology.atoms(), validAtoms):
                self.assertEqual(atom1.element, atom2.element)
                self.assertEqual(atom1.name, atom2.name)
                self.assertEqual(atom1.residue.name, atom2.residue.name)


    def testWriteAll(self):
        """Test all atoms are written when atomSubset is not specified"""
        
        with tempfile.TemporaryDirectory() as tempdir:
            filename = os.path.join(tempdir, 'temptraj.pdb')
            simulation = app.Simulation(self.pdb.topology, self.system, mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds))
            simulation.context.setPositions(self.pdb.positions)

            simulation.reporters.append(app.PDBReporter(filename, 1))
            simulation.step(1)

            # clear reporters to ensure PDBReporter calls writeFooter and file.close
            simulation.reporters.clear()

            # check if the output out trajectory contains the expected number of atoms
            checkpdb = app.PDBFile(filename)
            self.assertEqual(len(checkpdb.positions), simulation.topology.getNumAtoms())
            
            # check the positions are correct
            validPositions  = simulation.context.getState(getPositions=True).getPositions()
            # round to 4 decimal places before comparing
            for i in range(len(validPositions)):
                validPositions[i] = [round(validPositions[i][j]._value, 4) for j in (0, 1, 2)]*unit.nanometer

            for (p1, p2) in zip(checkpdb.positions, validPositions):
                assertVecAlmostEqual(p1, p2)

            # check elements and residue names are correct
            validAtoms = list(self.pdb.topology.atoms())
            for atom1, atom2 in zip(checkpdb.topology.atoms(), validAtoms):
                self.assertEqual(atom1.element, atom2.element)
                self.assertEqual(atom1.name, atom2.name)
                self.assertEqual(atom1.residue.name, atom2.residue.name)

    
    def testInvalidSubsets(self):
        """Test that an exception is raised when the indices in atomSubset are invalid"""

        with tempfile.TemporaryDirectory() as tempdir:
            for i, subset in enumerate([[-1,10], [0,99999], [0,0,0,1], [0.1,0.2], [5,10,0,9], ["C", "H"],[]]):
                
                filename = os.path.join(tempdir, 'temptraj'+str(i)+'.pdb')

                simulation = app.Simulation(self.pdb.topology, self.system, mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds))
                simulation.context.setPositions(self.pdb.positions)

                simulation.reporters.append(app.PDBReporter(filename, 1, atomSubset=subset))
                
                self.assertRaises(ValueError, lambda: simulation.step(1))

    
    def testBondSubset(self):
        """ Test that CONECT records are output correctly when using atomSubset"""

        # use a testcase that has CONECT records in the input pdb file
        ff = app.ForceField('amber14/protein.ff14SB.xml', 'amber14/GLYCAM_06j-1.xml','amber14/tip3pfb.xml')
        pdb = app.PDBFile('systems/glycopeptide.pdb')

        # add in water molecules
        modeller = app.Modeller(pdb.topology, pdb.positions)
        modeller.addSolvent(ff, padding=1.0*unit.nanometer)

        system = ff.createSystem(modeller.topology)
        
        with tempfile.TemporaryDirectory() as tempdir:
            filename = os.path.join(tempdir, 'temptraj.pdb')
            simulation = app.Simulation(modeller.topology, system, mm.LangevinMiddleIntegrator(1.0*unit.kelvin, 1.0/unit.picosecond, 0.0000001*unit.picoseconds))
            simulation.context.setPositions(modeller.positions)

            # output just the glycopeptide atoms
            atomSubset = list(range(pdb.topology.getNumAtoms()))
            simulation.reporters.append(app.PDBReporter(filename, 1, atomSubset=atomSubset))

            simulation.step(1)

            # clear reporters to ensure PDBReporter calls writeFooter and file.close
            simulation.reporters.clear()

            # for PyPy the above line is not enough, we need to force garbage collection to ensure the 
            # PDBReporter.__del__ method has been called before we open the file for reading 
            gc.collect()

            validpdb = pdb
            testpdb = app.PDBFile(filename)

            validBonds = list(validpdb.topology.bonds())
            testBonds = list(testpdb.topology.bonds())
            
            self.assertEqual(len(validBonds), len(testBonds))

            for validBond, testBond in zip(validBonds, testBonds):
                self.assertEqual(str(validBond), str(testBond))


class TestPDBxReporter(unittest.TestCase):
    def setUp(self):
        self.pdb = app.PDBFile('systems/alanine-dipeptide-explicit.pdb')
        self.forcefield = app.ForceField('amber99sbildn.xml','tip3p.xml')
        self.system = self.forcefield.createSystem(self.pdb.topology, nonbondedMethod=app.CutoffNonPeriodic, constraints=app.HBonds)

    def testSubset(self):
        """Test writing out a subset of atoms"""
    
        with tempfile.TemporaryDirectory() as tempdir:
            filename = os.path.join(tempdir, 'temptraj.pdbx')
            simulation = app.Simulation(self.pdb.topology, self.system, mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds))
            simulation.context.setPositions(self.pdb.positions)

            # just the alanine-dipeptide atoms
            subset = list(range(0,22))

            simulation.reporters.append(app.PDBxReporter(filename, 1, atomSubset=subset))
            simulation.step(1)

            # clear reporters to ensure PDBxReporter calls file.close
            simulation.reporters.clear()

            # check if the output out trajectory contains the expected number of atoms
            checkpdb = app.PDBxFile(filename)
            self.assertEqual(len(checkpdb.positions), len(subset))

            # check the positions are correct
            validPositions  = [simulation.context.getState(getPositions=True).getPositions()[i] for i in subset]
            # round to 5 decimal places before comparing
            for i in range(len(validPositions)):
                validPositions[i] = [round(validPositions[i][j]._value, 5) for j in (0, 1, 2)]*unit.nanometer
                
            for (p1, p2) in zip(checkpdb.positions, validPositions):
                assertVecAlmostEqual(p1, p2)

            # check elements and residue names are correct
            validAtoms = [list(self.pdb.topology.atoms())[i] for i in subset]
            for atom1, atom2 in zip(checkpdb.topology.atoms(), validAtoms):
                self.assertEqual(atom1.element, atom2.element)
                self.assertEqual(atom1.name, atom2.name)
                self.assertEqual(atom1.residue.name, atom2.residue.name)


    def testWriteAll(self):
        """Test all atoms are written when atomSubset is not specified"""
        
        with tempfile.TemporaryDirectory() as tempdir:
            filename = os.path.join(tempdir, 'temptraj.pdbx')
            simulation = app.Simulation(self.pdb.topology, self.system, mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds))
            simulation.context.setPositions(self.pdb.positions)

            simulation.reporters.append(app.PDBxReporter(filename, 1))
            simulation.step(1)

            # clear reporters to ensure PDBxReporter calls file.close
            simulation.reporters.clear()

            # check if the output out trajectory contains the expected number of atoms
            checkpdb = app.PDBxFile(filename)
            self.assertEqual(len(checkpdb.positions), simulation.topology.getNumAtoms())

            # check the positions are correct
            validPositions  = simulation.context.getState(getPositions=True).getPositions()
            # round to 5 decimal places before comparing
            for i in range(len(validPositions)):
                validPositions[i] = [round(validPositions[i][j]._value, 5) for j in (0, 1, 2)]*unit.nanometer

            for (p1, p2) in zip(checkpdb.positions, validPositions):
                assertVecAlmostEqual(p1, p2)

            # check elements and residue names are correct
            validAtoms = list(self.pdb.topology.atoms())
            for atom1, atom2 in zip(checkpdb.topology.atoms(), validAtoms):
                self.assertEqual(atom1.element, atom2.element)
                self.assertEqual(atom1.name, atom2.name)
                self.assertEqual(atom1.residue.name, atom2.residue.name)
    
    
    def testInvalidSubsets(self):
        """Test that an exception is raised when the indices in atomSubset are invalid"""
    
        with tempfile.TemporaryDirectory() as tempdir:
            for i,subset in enumerate([[-1,10], [0,99999], [0,0,0,1], [0.1,0.2], [5,10,0,9], ["C", "H"],[]]):
                
                filename = os.path.join(tempdir, 'temptraj'+str(i)+'.pdbx')

                simulation = app.Simulation(self.pdb.topology, self.system, mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds))
                simulation.context.setPositions(self.pdb.positions)

                simulation.reporters.append(app.PDBxReporter(filename, 1, atomSubset=subset))
                
                self.assertRaises(ValueError, lambda: simulation.step(1))


    def testBondSubset(self):
        """ Test that struct_conn records are output correctly when using atomSubset"""

        # use a testcase that has CONECT records in the input pdb file
        ff = app.ForceField('amber14/protein.ff14SB.xml', 'amber14/GLYCAM_06j-1.xml','amber14/tip3pfb.xml')
        pdb = app.PDBFile('systems/glycopeptide.pdb')

        # add in water molecules
        modeller = app.Modeller(pdb.topology, pdb.positions)
        modeller.addSolvent(ff, padding=1.0*unit.nanometer)

        system = ff.createSystem(modeller.topology)
        
        with tempfile.TemporaryDirectory() as tempdir:
            filename = os.path.join(tempdir, 'temptraj.pdbx')
            simulation = app.Simulation(modeller.topology, system, mm.LangevinMiddleIntegrator(1.0*unit.kelvin, 1.0/unit.picosecond, 0.0000001*unit.picoseconds))
            simulation.context.setPositions(modeller.positions)

            # output just the glycopeptide atoms
            atomSubset = list(range(pdb.topology.getNumAtoms()))
            simulation.reporters.append(app.PDBxReporter(filename, 1, atomSubset=atomSubset))

            simulation.step(1)

            # clear reporters to ensure PDBxReporter calls file.close
            simulation.reporters.clear()

            validpdb = pdb
            testpdb = app.PDBxFile(filename)

            validBonds = set(tuple(sorted((bond[0].index, bond[1].index))) for bond in validpdb.topology.bonds())
            testBonds = set(tuple(sorted((bond[0].index, bond[1].index))) for bond in testpdb.topology.bonds())
            
            self.assertEqual(len(validBonds), len(testBonds))

            for bond in validBonds:
                self.assertTrue(bond in testBonds)




    