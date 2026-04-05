from openmm import *
from openmm.app import *
from openmm.unit import *
import numpy as np
import unittest

class TestReplicaExchangeSampler(unittest.TestCase):
    def testTemperature(self):
        """Test a set of replicas that differ in temperature."""
        system = System()
        system.addParticle(1.0)
        force = CustomExternalForce('x*x+y*y+z*z')
        force.addParticle(0)
        system.addForce(force)
        integrator = LangevinIntegrator(300*kelvin, 10/picosecond, 0.01*picosecond)
        simulation = Simulation(Topology(), system, integrator, Platform.getPlatform('Reference'))
        states = [{'temperature':t*kelvin} for t in np.geomspace(300.0, 600.0, 5)]
        repex = ReplicaExchangeSampler(states, simulation, 20)
        energies = [0.0*kilojoules_per_mole]*len(states)
        exchanged = False

        def recordEnergies(repex):
            if repex.replicaStateIndex != list(range(len(states))):
                nonlocal exchanged
                exchanged = True
            for i in range(len(states)):
                simulation.context.setState(repex.replicaConformation[i])
                energies[repex.replicaStateIndex[i]] += simulation.context.getState(energy=True).getPotentialEnergy()

        repex.reporters.append(recordEnergies)
        for i in range(len(states)):
            repex.simulateReplica(i, 100)
        steps = 1000
        repex.simulate(steps)
        self.assertTrue(exchanged)
        for i, e in enumerate(energies):
            average = e/steps
            expected = 1.5*(states[i]['temperature']*MOLAR_GAS_CONSTANT_R)
            self.assertTrue(0.7 < average/expected < 1.3)

    def testParameter(self):
        """Test a set of replicas that differ in a force parameter."""
        system = System()
        system.addParticle(1.0)
        force = CustomExternalForce('0.5*k*x*x')
        force.addGlobalParameter('k', 1.0)
        force.addParticle(0)
        system.addForce(force)
        integrator = LangevinIntegrator(300*kelvin, 10/picosecond, 0.01*picosecond)
        simulation = Simulation(Topology(), system, integrator, Platform.getPlatform('Reference'))
        states = [{'k':k*kilojoules_per_mole/(nanometer**2)} for k in np.geomspace(5.0, 100.0, 5)]
        repex = ReplicaExchangeSampler(states, simulation, 20)
        r2 = [0.0*nanometer**2]*len(states)
        exchanged = False

        def recordDisplacements(repex):
            if repex.replicaStateIndex != list(range(len(states))):
                nonlocal exchanged
                exchanged = True
            for i in range(len(states)):
                x = repex.replicaConformation[i].getPositions()[0][0]
                r2[repex.replicaStateIndex[i]] += x*x

        repex.reporters.append(recordDisplacements)
        for i in range(len(states)):
            repex.simulateReplica(i, 100)
        steps = 2000
        repex.simulate(steps)
        self.assertTrue(exchanged)
        expected = 0.5*integrator.getTemperature()*MOLAR_GAS_CONSTANT_R
        for i in range(len(r2)):
            average = 0.5*states[i]['k']*r2[i]/steps
            self.assertTrue(0.7 < average/expected < 1.3)
