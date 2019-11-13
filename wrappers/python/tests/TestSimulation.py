import unittest
import tempfile
from datetime import datetime, timedelta
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


    def testLoadFromXML(self):
        """ Test creating a Simulation from XML files """
        pdb = PDBFile('systems/alanine-dipeptide-implicit.pdb')
        ff = ForceField('amber99sb.xml', 'tip3p.xml')
        system = ff.createSystem(pdb.topology)
        integrator = VerletIntegrator(0.001*picoseconds)
        context = Context(system, integrator)
        context.setPositions(pdb.positions)
        state = context.getState(getPositions=True, getForces=True,
                                 getVelocities=True, getEnergy=True)
        systemfn = tempfile.mktemp()
        integratorfn = tempfile.mktemp()
        statefn = tempfile.mktemp()
        with open(systemfn, 'w') as f:
            f.write(XmlSerializer.serialize(system))
        with open(integratorfn, 'w') as f:
            f.write(XmlSerializer.serialize(integrator))
        with open(statefn, 'w') as f:
            f.write(XmlSerializer.serialize(state))

        # Now create a Simulation
        sim = Simulation(pdb.topology, systemfn, integratorfn, state=statefn)

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

    def testStep(self):
        """Test the step() method."""
        pdb = PDBFile('systems/alanine-dipeptide-implicit.pdb')
        ff = ForceField('amber99sb.xml', 'tip3p.xml')
        system = ff.createSystem(pdb.topology)
        integrator = VerletIntegrator(0.001*picoseconds)

        # Create a Simulation.

        simulation = Simulation(pdb.topology, system, integrator, Platform.getPlatformByName('Reference'))
        simulation.context.setPositions(pdb.positions)
        simulation.context.setVelocitiesToTemperature(300*kelvin)
        self.assertEqual(0, simulation.currentStep)
        self.assertEqual(0*picoseconds, simulation.context.getState().getTime())

        # Take some steps and verify the simulation has advanced by the correct amount.

        simulation.step(23)
        self.assertEqual(23, simulation.currentStep)
        self.assertAlmostEqual(0.023, simulation.context.getState().getTime().value_in_unit(picoseconds))

    def testRunForClockTime(self):
        """Test the runForClockTime() method."""
        pdb = PDBFile('systems/alanine-dipeptide-implicit.pdb')
        ff = ForceField('amber99sb.xml', 'tip3p.xml')
        system = ff.createSystem(pdb.topology)
        integrator = VerletIntegrator(0.001*picoseconds)

        # Create a Simulation.

        simulation = Simulation(pdb.topology, system, integrator, Platform.getPlatformByName('Reference'))
        simulation.context.setPositions(pdb.positions)
        simulation.context.setVelocitiesToTemperature(300*kelvin)
        self.assertEqual(0, simulation.currentStep)
        self.assertEqual(0*picoseconds, simulation.context.getState().getTime())

        # Run for five seconds, the save both a checkpoint and a state.

        checkpointFile = tempfile.mktemp()
        stateFile = tempfile.mktemp()
        startTime = datetime.now()
        simulation.runForClockTime(5*seconds, checkpointFile=checkpointFile, stateFile=stateFile)
        endTime = datetime.now()

        # Make sure at least five seconds have elapsed, but no more than ten.

        self.assertTrue(endTime >= startTime+timedelta(seconds=5))
        self.assertTrue(endTime < startTime+timedelta(seconds=10))
        
        # Check that the time and step count are consistent.
        
        time = simulation.context.getState().getTime().value_in_unit(picoseconds)
        expectedTime = simulation.currentStep*integrator.getStepSize().value_in_unit(picoseconds)
        self.assertAlmostEqual(expectedTime, time)

        # Load the checkpoint and state and make sure they are both correct.

        velocities = simulation.context.getState(getVelocities=True).getVelocities()
        simulation.context.setVelocitiesToTemperature(300*kelvin)
        simulation.loadCheckpoint(checkpointFile)
        self.assertEqual(velocities, simulation.context.getState(getVelocities=True).getVelocities())
        simulation.context.setVelocitiesToTemperature(300*kelvin)
        simulation.loadState(stateFile)
        self.assertEqual(velocities, simulation.context.getState(getVelocities=True).getVelocities())

    def testWrappedCoordinates(self):
        """Test generating reports with and without wrapped coordinates."""
        pdb = PDBFile('systems/alanine-dipeptide-explicit.pdb')
        ff = ForceField('amber99sb.xml', 'tip3p.xml')
        system = ff.createSystem(pdb.topology, nonbondedMethod=CutoffPeriodic, constraints=HBonds)
        integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

        class CompareCoordinatesReporter(object):
            def __init__(self, periodic):
                self.periodic = periodic
                self.interval = 100
                
            def describeNextReport(self, simulation):
                steps = self.interval - simulation.currentStep%self.interval
                return (steps, True, False, False, False, self.periodic)
        
            def report(self, simulation, state):
                state2 = simulation.context.getState(getPositions=True, enforcePeriodicBox=self.periodic)
                assert state.getPositions() == state2.getPositions()

        # Create a Simulation.

        simulation = Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.context.setVelocitiesToTemperature(300*kelvin)
        simulation.reporters.append(CompareCoordinatesReporter(False))
        simulation.reporters.append(CompareCoordinatesReporter(True))
        
        # Run for a little while and make sure the reporters don't find any problems.
        
        simulation.step(500)


if __name__ == '__main__':
    unittest.main()
