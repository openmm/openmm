import os
import unittest
import tempfile
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
        file = tempfile.NamedTemporaryFile(delete=False)
        self.simulation.reporters.append(app.CheckpointReporter(file, 1))
        self.simulation.step(1)

        # get the current positions
        positions = self.simulation.context.getState(getPositions=True).getPositions()

        # now set the positions into junk...
        self.simulation.context.setPositions([mm.Vec3(0, 0, 0)] * len(positions))

        # then reload the right positions from the checkpoint
        file.close()
        with open(file.name, 'rb') as f:
            self.simulation.context.loadCheckpoint(f.read())
        os.unlink(file.name)

        newPositions = self.simulation.context.getState(getPositions=True).getPositions()
        self.assertSequenceEqual(positions, newPositions)

if __name__ == '__main__':
    unittest.main()
