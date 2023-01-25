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
    
    def testInvalidSubsets(self):
        """Test that an exception is raised when the indices in atomSubset are invalid"""
    
        with tempfile.TemporaryDirectory() as tempdir:
            filename = os.path.join(tempdir, 'temptraj.pdb')

            for subset in [[-1,10], [0,99999], [0,0,0,1], [0.1,0.2]]:

                simulation = app.Simulation(self.pdb.topology, self.system, mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds))
                simulation.context.setPositions(self.pdb.positions)

                simulation.reporters.append(app.PDBReporter(filename, 10, atomSubset=subset))
                
                self.assertRaises(ValueError, lambda: simulation.step(10))
    