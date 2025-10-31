import unittest
from openmm import *
from openmm.unit import *
import numpy as np
import copy

def compute(state):
    """This is a computation function used by the test cases."""
    pos = state.getPositions(asNumpy=True).value_in_unit(nanometer)
    k = state.getParameters()['k']
    energy = k*np.sum(pos*pos)
    force = -0.5*k*pos
    return energy*kilojoules_per_mole, force*kilojoules_per_mole/nanometer

class TestPythonForce(unittest.TestCase):
    """Test the PythonForce class"""

    def testComputeForce(self):
        """Test using PythonForce to compute forces."""
        system = System()
        for i in range(5):
            system.addParticle(1.0)
        force = PythonForce(compute, {'k':2.5})
        system.addForce(force)
        positions = np.random.rand(5, 3)
        integrator = VerletIntegrator(0.001)
        context = Context(system, integrator, Platform.getPlatform('Reference'))
        context.setPositions(positions)
        state = context.getState(energy=True, forces=True)
        self.assertAlmostEqual(2.5*np.sum(positions*positions), state.getPotentialEnergy().value_in_unit(kilojoules_per_mole), places=5)
        self.assertTrue(np.allclose(-1.25*positions, state.getForces(asNumpy=True).value_in_unit(kilojoules_per_mole/nanometer)))

    def testExceptions(self):
        """Test that PythonForce handles exceptions correctly."""
        def compute2(state):
            raise ValueError('This should fail')

        system = System()
        system.addParticle(1.0)
        force = PythonForce(compute2)
        system.addForce(force)
        positions = np.random.rand(1, 3)
        integrator = VerletIntegrator(0.001)
        context = Context(system, integrator)
        context.setPositions(positions)
        with self.assertRaises(OpenMMException) as cm:
            context.getState(energy=True)
        self.assertEqual('This should fail', str(cm.exception))

    def testSerialize(self):
        """Test that PythonForce can be serialized."""
        force1 = PythonForce(compute, {'k':2.5})
        force1.setUsesPeriodicBoundaryConditions(True)

        # Make a copy by serializing and the deserializing it.

        force2 = copy.copy(force1)

        # They should be identical.

        self.assertEqual(XmlSerializer.serialize(force1), XmlSerializer.serialize(force2))
        self.assertEqual(dict(force2.getGlobalParameters()), {'k':2.5})
        self.assertTrue(force2.usesPeriodicBoundaryConditions())

        # A locally defined function cannot be pickled.  We should not be able to serialize a force
        # that uses it.

        def compute2(state):
            return 1.0, np.zeros(len(state.getPositions()), 3)

        force3 = PythonForce(compute2)
        with self.assertRaises(OpenMMException):
            XmlSerializer.serialize(force3)

    def testMinimization(self):
        """Test that PythonForce works correctly with the minimizer."""
        system = System()
        for i in range(5):
            system.addParticle(1.0)
        force = PythonForce(compute, {'k':2.5})
        system.addForce(force)
        positions = np.random.rand(5, 3)
        integrator = VerletIntegrator(0.001)
        context = Context(system, integrator, Platform.getPlatform('Reference'))
        context.setPositions(positions)

        # The PythonForce and the MinimizationReporter both involve calling back into Python code,
        # possibly from different threads.  Make sure it doesn't cause any problems.

        class Reporter(MinimizationReporter):
            count = 0
            def report(self, iteration, x, grad, args):
                self.count += 1
                return False

        reporter = Reporter()
        LocalEnergyMinimizer.minimize(context, tolerance=1e-3, reporter=reporter)
        self.assertTrue(reporter.count > 0)
        state = context.getState(energy=True, positions=True)
        self.assertAlmostEqual(0.0, state.getPotentialEnergy().value_in_unit(kilojoules_per_mole))

    def testMemory(self):
        """Test for memory leaks in the Python/C++ interface."""
        try:
            import resource
        except:
            # The resource module is not available on Windows.
            return
        system = System()
        for i in range(1000):
            system.addParticle(1.0)
        force = PythonForce(compute, {'k':2.5})
        system.addForce(force)
        positions = np.random.rand(1000, 3)
        integrator = VerletIntegrator(0.001)
        context = Context(system, integrator, Platform.getPlatform('Reference'))
        context.setPositions(positions)
        integrator.step(5000)
        memory1 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        integrator.step(5000)
        memory2 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        self.assertTrue(memory2 < 1.05*memory1)

    def testDtypes(self):
        """Test returning forces with different types."""
        for dtype in [np.float32, np.float64, int]:
            def compute2(state):
                return 0, np.array([[1,2,3],[4,5,6]], dtype=dtype)

            system = System()
            system.addParticle(1.0)
            system.addParticle(1.0)
            force = PythonForce(compute2)
            system.addForce(force)
            positions = np.random.rand(2, 3)
            integrator = VerletIntegrator(0.001)
            context = Context(system, integrator)
            context.setPositions(positions)
            forces = context.getState(forces=True).getForces().value_in_unit(kilojoules_per_mole/nanometer)
            self.assertEqual(Vec3(1,2,3), forces[0])
            self.assertEqual(Vec3(4,5,6), forces[1])

if __name__ == '__main__':
    unittest.main()
