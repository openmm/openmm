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
        states = [{'k':k*kilojoules_per_mole/(nanometer**2)} for k in np.geomspace(10.0, 100.0, 5)]
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
        system = System()
        force = CustomExternalForce('0.5*k*(x*x+y*y+z*z)')
        force.addGlobalParameter('k', 1.0)
        system.addForce(force)
        for i in range(3):
            system.addParticle(1.0)
            force.addParticle(0)
        states = [{'k':k} for k in (200.0, 300.0, 400.0)]
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as logFile:
            with tempfile.NamedTemporaryFile(mode='w', delete=False) as energyFile:
                with tempfile.NamedTemporaryFile(mode='w', delete=False) as checkpointFile:
                    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.001*picosecond)
                    simulation = Simulation(Topology(), system, integrator, Platform.getPlatform('Reference'))
                    simulation.context.setPositions([Vec3(0, 0, 0)]*3)
                    sampler = ExpandedEnsembleSampler(states, simulation, 5, reportInterval=5, logFile=logFile.name,
                                                      energyFile=energyFile.name, checkpointFile=checkpointFile.name)

                    # Run a simulation.

                    step = []
                    iteration = []
                    stateIndex = []
                    weights = []
                    energies = []

                    def runIteration():
                        simulation.step(5)
                        step.append(simulation.currentStep)
                        iteration.append(sampler.currentIteration)
                        stateIndex.append(sampler.currentStateIndex)
                        weights.append(sampler.weights)
                        kT = MOLAR_GAS_CONSTANT_R*simulation.integrator.getTemperature()
                        energies.append(sampler._sampler.computeAllEnergies()/kT)
                        sampler._sampler.applyState(sampler.currentStateIndex)

                    try:
                        for _ in range(4):
                            runIteration()
                    except PermissionError:
                        # tempfile is kind of broken on Windows.  Just skip the test.
                        return
                    state1 = simulation.context.getState(positions=True, velocities=True, parameters=True)

                    # Delete all objects from the simulation and create a new one, telling it to resume from the files.

                    del sampler
                    del simulation
                    del integrator
                    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.001*picosecond)
                    simulation = Simulation(Topology(), system, integrator, Platform.getPlatform('Reference'))
                    sampler = ExpandedEnsembleSampler(states, simulation, 5, reportInterval=5, logFile=logFile.name,
                                                      energyFile=energyFile.name, checkpointFile=checkpointFile.name,
                                                      resume=True)

                    # Make sure everything was loaded correctly.

                    state2 = simulation.context.getState(positions=True, velocities=True, parameters=True)
                    self.assertEqual(XmlSerializer.serialize(state1), XmlSerializer.serialize(state2))
                    self.assertEqual(step[-1], simulation.currentStep)
                    self.assertEqual(iteration[-1], sampler.currentIteration)
                    self.assertEqual(stateIndex[-1], sampler.currentStateIndex)
                    self.assertEqual(weights[-1], sampler.weights)

                    # Generate some more output.

                    for _ in range(4):
                        runIteration()

                    # Check the log file.

                    logFile.close()
                    with open(logFile.name) as input:
                        lines = input.readlines()[1:]
                    os.remove(logFile.name)
                    self.assertEqual(8, len(lines))
                    for i, line in enumerate(lines):
                        fields = line.split(',')
                        self.assertEqual(int(fields[0]), step[i])
                        self.assertEqual(int(fields[1]), iteration[i])
                        self.assertEqual(int(fields[2]), stateIndex[i])
                        self.assertTrue(np.allclose([float(x) for x in fields[3:]], weights[i]))

                    # Check the energy file.

                    energyFile.close()
                    with open(energyFile.name) as input:
                        lines = input.readlines()[1:]
                    os.remove(energyFile.name)
                    self.assertEqual(8, len(lines))
                    for i, line in enumerate(lines):
                        fields = line.split(',')
                        self.assertEqual(int(fields[0]), step[i])
                        self.assertTrue(np.allclose([float(x) for x in fields[1:]], energies[i]))
