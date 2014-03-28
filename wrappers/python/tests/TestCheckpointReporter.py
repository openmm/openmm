import unittest
import tempfile
import numpy as np
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit


class TestCheckpointReporter(unittest.TestCase):
    def setUp(self):
        with open('systems/alanine-dipeptide-implicit.pdb') as f:
            pdb = app.PDBFile(f)
        forcefield = app.ForceField('amber99sbildn.xml')
        system = forcefield.createSystem(pdb.topology,
            nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=1.0*unit.nanometers,
            constraints=app.HBonds)
        self.simulation = app.Simulation(pdb.topology, system, mm.VerletIntegrator(0.002*unit.picoseconds))
        self.simulation.context.setPositions(pdb.positions)

    def test_1(self):
        file = tempfile.NamedTemporaryFile()
        self.simulation.reporters.append(app.CheckpointReporter(file, 1))
        self.simulation.step(1)

        # get the current positions
        positions = self.simulation.context.getState(getPositions=True).getPositions(asNumpy=True)._value
        # now set the positions into junk...
        self.simulation.context.setPositions(np.random.random(positions.shape))
        # then reload the right positions from the checkpoint
        with open(file.name, 'rb') as f:
            self.simulation.context.loadCheckpoint(f.read())
        file.close()

        newPositions = self.simulation.context.getState(getPositions=True).getPositions(asNumpy=True)._value
        np.testing.assert_array_equal(positions, newPositions)

if __name__ == '__main__':
    unittest.main()
