from openmm import *
from openmm.app import *
from openmm.unit import *
import numpy as np
import os
import tempfile
import unittest

class TestExpandedEnsembleSampler(unittest.TestCase):
    def testTemperature(self):
        """Test a set of states that differ in temperature."""
        system = System()
        system.addParticle(1.0)
        force = CustomExternalForce('x*x+y*y+z*z')
        force.addParticle(0)
        system.addForce(force)
        states = [{'temperature':t*kelvin} for t in np.geomspace(300.0, 600.0, 5)]
        for reinitialize in [False, True]:
            integrator = LangevinIntegrator(300*kelvin, 10/picosecond, 0.01*picosecond)
            simulation = Simulation(Topology(), system, integrator, Platform.getPlatform('Reference'))
            simulation.context.setPositions([Vec3(0, 0, 0)])
            sampler = ExpandedEnsembleSampler(states, simulation, 10, reinitialize)

            # Run for a little while to let the weights stabilize.

            sampler.step(10000)

            # Run for a while and record the states and energies.

            energies = [[] for _ in range(len(states))]
            iterations = 20000
            for i in range(iterations):
                sampler.step(10)
                energies[sampler.currentStateIndex].append(simulation.context.getState(energy=True).getPotentialEnergy())

            # Check that it spent roughly equal time in each state, and that the energies are correct.

            for energy, state in zip(energies, states):
                n = len(energy)
                assert iterations/10 < n < iterations/2
                average = sum(energy)/n
                expected = 1.5*(state['temperature']*MOLAR_GAS_CONSTANT_R)
                self.assertTrue(0.7 < average/expected < 1.3)

    def testParameter(self):
        """Test a set of states that differ in a force parameter."""
        system = System()
        system.addParticle(1.0)
        force = CustomExternalForce('0.5*k*x*x')
        force.addGlobalParameter('k', 1.0)
        force.addParticle(0)
        system.addForce(force)
        states = [{'k':k*kilojoules_per_mole/(nanometer**2)} for k in np.geomspace(5.0, 100.0, 5)]
        for reinitialize in [False, True]:
            integrator = LangevinIntegrator(300*kelvin, 10/picosecond, 0.01*picosecond)
            simulation = Simulation(Topology(), system, integrator, Platform.getPlatform('Reference'))
            simulation.context.setPositions([Vec3(0, 0, 0)])
            sampler = ExpandedEnsembleSampler(states, simulation, 10, reinitialize)

            # Run for a little while to let the weights stabilize.

            sampler.step(10000)

            # Run for a while and record the states and displacements.

            r2 = [[] for _ in range(len(states))]
            iterations = 20000
            for i in range(iterations):
                sampler.step(10)
                x = simulation.context.getState(positions=True).getPositions()[0][0]
                r2[sampler.currentStateIndex].append(x*x)

            # Check that it spent roughly equal time in each state, and that the energies are correct.

            expected = 0.5*integrator.getTemperature()*MOLAR_GAS_CONSTANT_R
            for i in range(len(r2)):
                n = len(r2[i])
                assert iterations/10 < n < iterations/2
                average = 0.5*states[i]['k']*sum(r2[i])/n
                self.assertTrue(0.7 < average/expected < 1.3)

    def testReporter(self):
        """Test reporting output from an expanded ensemble simulation."""
        pdb = PDBFile('systems/alanine-dipeptide-implicit.pdb')
        ff = ForceField('amber19-all.xml')
        system = ff.createSystem(pdb.topology)
        integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.001*picosecond)
        simulation = Simulation(pdb.topology, system, integrator, Platform.getPlatform('Reference'))
        simulation.context.setPositions(pdb.positions)
        states = [{'temperature':t*kelvin} for t in [300.0, 320.0, 340.0]]
        with tempfile.TemporaryDirectory() as directory:
            logFile = os.path.join(directory, 'log.csv')
            sampler = ExpandedEnsembleSampler(states, simulation, 5, reportInterval=5, logFile=logFile)

            # Run a simulation.

            step = []
            iteration = []
            state = []
            weights = []
            for i in range(4):
                simulation.step(5)
                step.append(simulation.currentStep)
                iteration.append(sampler.currentIteration)
                state.append(sampler.currentStateIndex)
                weights.append(sampler.weights)

            # Check the log file.

            lines = open(logFile).readlines()[1:]
            for i, line in enumerate(lines):
                fields = line.split(',')
                self.assertEqual(int(fields[0]), step[i])
                self.assertEqual(int(fields[1]), iteration[i])
                self.assertEqual(int(fields[2]), state[i])
                self.assertTrue(np.allclose([float(x) for x in fields[3:]], weights[i]))
