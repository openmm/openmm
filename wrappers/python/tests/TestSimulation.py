import unittest
import tempfile
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *

class TestSimulation(unittest.TestCase):
    """Test the Simulation class"""

    def testCheckpointing(self):
        """Test that checkpointing works correctly."""
        pdb = PDBFile('systems/alanine-dipeptide-implicit.pdb')
        ff = ForceField('amber99sb.xml', 'tip3p.xml')
        system = ff.createSystem(pdb.topology)
        integrator = VerletIntegrator(0.001*picoseconds)

        # Create a Simulation.

        simulation = Simulation(pdb.topology, system, integrator, Platform.getPlatformByName('Reference'))
        simulation.context.setPositions(pdb.positions)
        simulation.context.setVelocitiesToTemperature(300*kelvin)
        initialState = simulation.context.getState(getPositions=True, getVelocities=True)

        # Create a checkpoint.

        filename = tempfile.mktemp()
        simulation.saveCheckpoint(filename)

        # Take a few steps so the positions and velocities will be different.

        simulation.step(2)
        state = simulation.context.getState(getPositions=True, getVelocities=True)
        self.assertNotEqual(initialState.getPositions(), state.getPositions())
        self.assertNotEqual(initialState.getVelocities(), state.getVelocities())

        # Reload the checkpoint and see if it resets them correctly.

        simulation.loadCheckpoint(filename)
        state = simulation.context.getState(getPositions=True, getVelocities=True)
        self.assertEqual(initialState.getPositions(), state.getPositions())
        self.assertEqual(initialState.getVelocities(), state.getVelocities())


    def testSaveState(self):
        """Test that saving States works correctly."""
        pdb = PDBFile('systems/alanine-dipeptide-implicit.pdb')
        ff = ForceField('amber99sb.xml', 'tip3p.xml')
        system = ff.createSystem(pdb.topology)
        integrator = VerletIntegrator(0.001*picoseconds)

        # Create a Simulation.

        simulation = Simulation(pdb.topology, system, integrator, Platform.getPlatformByName('Reference'))
        simulation.context.setPositions(pdb.positions)
        simulation.context.setVelocitiesToTemperature(300*kelvin)
        initialState = simulation.context.getState(getPositions=True, getVelocities=True)

        # Create a state.

        filename = tempfile.mktemp()
        simulation.saveState(filename)

        # Take a few steps so the positions and velocities will be different.

        simulation.step(2)
        state = simulation.context.getState(getPositions=True, getVelocities=True)
        self.assertNotEqual(initialState.getPositions(), state.getPositions())
        self.assertNotEqual(initialState.getVelocities(), state.getVelocities())

        # Reload the state and see if it resets them correctly.

        simulation.loadState(filename)
        state = simulation.context.getState(getPositions=True, getVelocities=True)
        self.assertEqual(initialState.getPositions(), state.getPositions())
        self.assertEqual(initialState.getVelocities(), state.getVelocities())


if __name__ == '__main__':
    unittest.main()
