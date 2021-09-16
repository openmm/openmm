import unittest
import tempfile
from openmm import app
import openmm as mm
from openmm import unit
import os


class TestStateDataReporter(unittest.TestCase):
    def setUp(self):
        self.pdb = app.PDBFile('systems/alanine-dipeptide-implicit.pdb')
        self.forcefield = app.ForceField('amber99sbildn.xml')
        self.system = self.forcefield.createSystem(self.pdb.topology, nonbondedMethod=app.CutoffNonPeriodic, constraints=app.HBonds)

    def testAppend(self):
        """Test appending to the StateDataReporter output."""
        with tempfile.TemporaryDirectory() as tempdir:
            # Create a Simulation and produce 10 steps of output.
            filename = os.path.join(tempdir, 'templog.txt')
            simulation = app.Simulation(self.pdb.topology, self.system, mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds))
            simulation.context.setPositions(self.pdb.positions)
            simulation.reporters.append(app.StateDataReporter(filename, 1, step=True))
            simulation.step(10)

            # Create a new Simulation and append 5 more steps of output.
            simulation = app.Simulation(self.pdb.topology, self.system, mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds))
            simulation.context.setPositions(self.pdb.positions)
            simulation.reporters.append(app.StateDataReporter(filename, 1, step=True, append=True))
            simulation.step(5)
            
            # See if the log contents are correct.
            del simulation
            lines = open(filename).read().split('\n')
            self.assertTrue(lines[0].startswith('#'))
            for i in range(10):
                self.assertEqual(lines[i+1], f'{i+1}')
            for i in range(5):
                self.assertEqual(lines[i+11], f'{i+1}')


if __name__ == '__main__':
    unittest.main()
