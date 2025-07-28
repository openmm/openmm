import os
import unittest
import tempfile
from io import BytesIO, StringIO
from openmm import app
import openmm as mm
from openmm import unit


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
        """Test CheckpointReporter."""
        for writeState in [True, False]:
            with tempfile.TemporaryDirectory() as tempdir:
                filename = os.path.join(tempdir, 'checkpoint')
                self.simulation.reporters.clear()
                self.simulation.reporters.append(app.CheckpointReporter(filename, 1, writeState=writeState))
                self.simulation.step(1)
        
                # get the current positions
                positions = self.simulation.context.getState(getPositions=True).getPositions()
        
                # now set the positions into junk...
                self.simulation.context.setPositions([mm.Vec3(0, 0, 0)] * len(positions))
        
                # then reload the right positions from the checkpoint
                if writeState:
                    self.simulation.loadState(filename)
                else:
                    self.simulation.loadCheckpoint(filename)
        
                newPositions = self.simulation.context.getState(getPositions=True).getPositions()
                self.assertSequenceEqual(positions, newPositions)

    def testFileObj(self):
        """Test writing to a file object.  This should truncate so that only the most recent frame is present in the output."""

        # Test checkpoint saving.

        checkpointBuffer = BytesIO()
        self.simulation.reporters.clear()
        self.simulation.reporters.append(app.CheckpointReporter(checkpointBuffer, 1, writeState=False))
        self.simulation.step(5)
        checkpointData = checkpointBuffer.getvalue()

        checkpointBuffer = BytesIO()
        self.simulation.saveCheckpoint(checkpointBuffer)
        self.assertSequenceEqual(checkpointData, checkpointBuffer.getvalue())

        # Test state saving.

        stateBuffer = StringIO()
        self.simulation.reporters.clear()
        self.simulation.reporters.append(app.CheckpointReporter(stateBuffer, 1, writeState=True))
        self.simulation.step(5)
        stateData = stateBuffer.getvalue()

        stateBuffer = StringIO()
        self.simulation.saveState(stateBuffer)
        self.assertSequenceEqual(stateData, stateBuffer.getvalue())

if __name__ == '__main__':
    unittest.main()
