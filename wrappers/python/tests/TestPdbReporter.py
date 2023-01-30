import unittest
import tempfile
from openmm import app
import openmm as mm
from openmm import unit
import os

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

            simulation.reporters.append(app.PDBReporter(filename, 10, atomSubset=subset))
            simulation.step(10)

            # check if the output out trajectory contains the expected number of atoms
            checkpdb = app.PDBFile(filename)
            self.assertEqual(len(checkpdb.positions), len(subset))

    def testWriteAll(self):
        """Test all atoms are written when atomSubset is not specified"""
        
        with tempfile.TemporaryDirectory() as tempdir:
            filename = os.path.join(tempdir, 'temptraj.pdb')
            simulation = app.Simulation(self.pdb.topology, self.system, mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds))
            simulation.context.setPositions(self.pdb.positions)


            simulation.reporters.append(app.PDBReporter(filename, 10))
            simulation.step(10)

            # check if the output out trajectory contains the expected number of atoms
            checkpdb = app.PDBFile(filename)
            self.assertEqual(len(checkpdb.positions), simulation.topology.getNumAtoms())
    
    def testInvalidSubsets(self):
        """Test that an exception is raised when the indices in atomSubset are invalid"""
    
        with tempfile.TemporaryDirectory() as tempdir:
            filename = os.path.join(tempdir, 'temptraj.pdb')

            for subset in [[-1,10], [0,99999], [0,0,0,1], [0.1,0.2], [5,10,0,9], ["C", "H"],[]]:

                simulation = app.Simulation(self.pdb.topology, self.system, mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds))
                simulation.context.setPositions(self.pdb.positions)

                simulation.reporters.append(app.PDBReporter(filename, 10, atomSubset=subset))
                
                self.assertRaises(ValueError, lambda: simulation.step(10))
    
    def testBondSubset(self):
        """ Test that CONNECT records are output correctly when using atomSubset"""

        # use a testcase that has CONNECT records in the input pdb file
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

            simulation.minimizeEnergy(maxIterations=100)

            # output just the glycopeptide atoms
            atomSubset=list(range(pdb.topology.getNumAtoms()))
            simulation.reporters.append(app.PDBReporter(filename,1, atomSubset=atomSubset))

            simulation.step(1)

            # clear reporters to ensure PDB output is closed and footer written.
            simulation.reporters.clear()

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

            simulation.reporters.append(app.PDBxReporter(filename, 10, atomSubset=subset))
            simulation.step(10)

            # check if the output out trajectory contains the expected number of atoms
            checkpdb = app.PDBxFile(filename)
            self.assertEqual(len(checkpdb.positions), len(subset))

    def testWriteAll(self):
        """Test all atoms are written when atomSubset is not specified"""
        
        with tempfile.TemporaryDirectory() as tempdir:
            filename = os.path.join(tempdir, 'temptraj.pdbx')
            simulation = app.Simulation(self.pdb.topology, self.system, mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds))
            simulation.context.setPositions(self.pdb.positions)


            simulation.reporters.append(app.PDBxReporter(filename, 10))
            simulation.step(10)

            # check if the output out trajectory contains the expected number of atoms
            checkpdb = app.PDBxFile(filename)
            self.assertEqual(len(checkpdb.positions), simulation.topology.getNumAtoms())
    
    def testInvalidSubsets(self):
        """Test that an exception is raised when the indices in atomSubset are invalid"""
    
        with tempfile.TemporaryDirectory() as tempdir:
            filename = os.path.join(tempdir, 'temptraj.pdbx')

            for subset in [[-1,10], [0,99999], [0,0,0,1], [0.1,0.2], [5,10,0,9], ["C", "H"],[]]:

                simulation = app.Simulation(self.pdb.topology, self.system, mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds))
                simulation.context.setPositions(self.pdb.positions)

                simulation.reporters.append(app.PDBxReporter(filename, 10, atomSubset=subset))
                
                self.assertRaises(ValueError, lambda: simulation.step(10))


    def testBondSubset(self):
        """ Test that CONNECT records are output correctly when using atomSubset"""

        # use a testcase that has CONNECT records in the input pdb file
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

            simulation.minimizeEnergy(maxIterations=100)

            # output just the glycopeptide atoms
            atomSubset=list(range(pdb.topology.getNumAtoms()))
            simulation.reporters.append(app.PDBxReporter(filename,1, atomSubset=atomSubset))

            simulation.step(1)

            # clear reporters to ensure PDB output is closed and footer written.
            simulation.reporters.clear()

            validpdb = pdb
            testpdb = app.PDBxFile(filename)

            validBonds = list(validpdb.topology.bonds())
            testBonds = list(testpdb.topology.bonds())
            
            self.assertEqual(len(validBonds), len(testBonds))

            for validBond, testBond in zip(validBonds, testBonds):
                self.assertEqual(str(validBond), str(testBond))


    