from openmm import *
from openmm.app import *
from openmm.unit import *
from openmm.app.internal.xtc_utils import read_xtc
import numpy as np
import os
import tempfile
import unittest

class TestReplicaExchangeSampler(unittest.TestCase):
    def testTemperature(self):
        """Test a set of replicas that differ in temperature."""
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
            repex = ReplicaExchangeSampler(states, simulation, 20, reinitialize)
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
                self.assertEqual(steps*20+100, repex.replicaConformation[i].getStepCount())

    def testParameter(self):
        """Test a set of replicas that differ in a force parameter."""
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
            repex = ReplicaExchangeSampler(states, simulation, 20, reinitialize)
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

    def testReporter(self):
        """Test reporting output from a replica exchange simulation."""
        # Set up a replica exchange simulation.

        pdb = PDBFile('systems/alanine-dipeptide-implicit.pdb')
        ff = ForceField('amber19-all.xml')
        system = ff.createSystem(pdb.topology)
        integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.001*picosecond)
        simulation = Simulation(pdb.topology, system, integrator, Platform.getPlatform('Reference'))
        simulation.context.setPositions(pdb.positions)
        states = [{'temperature':t*kelvin} for t in [300.0, 320.0, 340.0]]
        sampler = ReplicaExchangeSampler(states, simulation, 5)

        # Add reporters to it.

        stateIndices = []
        energies = []
        volumes = []
        conformations = []

        def report(sampler):
            if sampler.currentIteration % 3 == 0:
                stateIndices.append(sampler.replicaStateIndex[:])
                energies.append(sampler.replicaStateEnergy[:])
                volumes.append([sampler.replicaConformation[i].getPeriodicBoxVolume().value_in_unit(nanometer**3) for i in range(len(states))])
                conformations.append(sampler.replicaConformation[:])

        with tempfile.TemporaryDirectory() as directory:
            sampler.reporters.append(ReplicaExchangeReporter(directory, 3, sampler, trajectoryPerState=True, trajectoryPerReplica=True, energy=True, volume=True))
            sampler.reporters.append(report)

            # Generate some output.

            sampler.simulate(15)

            # Delete all objects from the simulation and create a new one, telling it to resume from the files.

            del sampler
            del simulation
            del integrator
            integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.001*picosecond)
            simulation = Simulation(pdb.topology, system, integrator, Platform.getPlatform('Reference'))
            sampler = ReplicaExchangeSampler(states, simulation, 5)
            sampler.reporters.append(ReplicaExchangeReporter(directory, 3, sampler, trajectoryPerState=True, trajectoryPerReplica=True, energy=True, volume=True, resume=True))
            sampler.reporters.append(report)

            # Check that it loaded the checkpoints correctly.

            for i in range(len(states)):
                xml1 = XmlSerializer.serialize(sampler.replicaConformation[i])
                xml2 = open(os.path.join(directory, f'checkpoint_{i}.xml')).read()
                self.assertEqual(xml1, xml2)

            # Generate some more output.

            sampler.simulate(15)

            # Check the log file.

            with open(os.path.join(directory, 'log.csv')) as input:
                lines = input.readlines()[1:]
            for i, line in enumerate(lines):
                fields = [int(x) for x in line.split(',')]
                self.assertEqual(fields[0], 3*(i+1))
                self.assertEqual(fields[1], 15*(i+1))
                for j in range(len(states)):
                    self.assertEqual(stateIndices[i][j], fields[j+2])

            # Check the energy file.

            energy = np.loadtxt(os.path.join(directory, 'energy.csv'), delimiter=',').reshape(-1, len(states), len(states))
            for i in range(10):
                for j in range(len(states)):
                    for k in range(len(states)):
                        self.assertAlmostEqual(energy[i][j][k], energies[i][j][k]/sampler._kT[k])

            # Check the volume file.

            volume = np.loadtxt(os.path.join(directory, 'volume.csv'), delimiter=',').reshape(-1, len(states))
            for i in range(10):
                for j in range(len(states)):
                    self.assertAlmostEqual(volume[i][j], volumes[i][j])

            # Check the trajectory files.

            for i in range(len(states)):
                # Check the per-replica trajectories.

                coords_read, box_read, time, step = read_xtc(os.path.join(directory, f'replica_{i}.xtc').encode('utf-8'))
                self.assertTrue(np.allclose(step, np.linspace(15, 150, 10)))
                self.assertTrue(np.allclose(time, 0.005*np.linspace(3, 30, 10)))
                for j in range(10):
                    conf = conformations[j][i].getPositions().value_in_unit(nanometers)
                    self.assertTrue(np.allclose(conf, coords_read[:,:,j], atol=0.001))

                # Check the per-state trajectories.

                coords_read, box_read, time, step = read_xtc(os.path.join(directory, f'state_{i}.xtc').encode('utf-8'))
                self.assertTrue(np.allclose(step, np.linspace(15, 150, 10)))
                self.assertTrue(np.allclose(time, 0.005*np.linspace(3, 30, 10)))
                for j in range(10):
                    replica = stateIndices[j].index(i)
                    conf = conformations[j][replica].getPositions().value_in_unit(nanometers)
                    self.assertTrue(np.allclose(conf, coords_read[:,:,j], atol=0.001))

            # Creating a new reporter for the same directory should fail.

            with self.assertRaises(ValueError):
                ReplicaExchangeReporter(directory, 3, sampler)
            del sampler
            del simulation
            del integrator


    def testKTFallback(self):
        """Ensure sampler._kT is constructed from the integrator when states omit temperature."""
        system = System()
        system.addParticle(6.0)
        force = CustomExternalForce('0.5*k*x*x')
        force.addGlobalParameter('k', 1.0)
        # apply the force to the single particle
        force.addParticle(0)
        system.addForce(force)

        top = Topology()
        chain = top.addChain()
        res = top.addResidue('RES', chain)
        top.addAtom('C1', Element.getByAtomicNumber(6), res)

        # States differ in a force parameter but do not specify 'temperature'.
        states = [{'k':k*kilojoule_per_mole/(nanometer**2)} for k in [1.0, 2.0, 3.0]]
        integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 1 * femtosecond)
        simulation = Simulation(top, system, integrator, Platform.getPlatform('Reference'))
        simulation.context.setPositions([Vec3(0, 0, 0)])

        repex = ReplicaExchangeSampler(states, simulation, 5)

        # _kT should be a per-state list constructed from integrator.getTemperature()
        expected_kT = MOLAR_GAS_CONSTANT_R * integrator.getTemperature()
        self.assertEqual(len(repex._kT), len(states))
        for kT in repex._kT:
            self.assertEqual(kT, expected_kT)

        # Add reporters to it.
        energies = []

        def report(sampler):
            if sampler.currentIteration % 1 == 0:
                energies.append(sampler.replicaStateEnergy[:])


        with tempfile.TemporaryDirectory() as directory:
            # add reporters that write into the temporary directory and collect runtime data
            repex.reporters.append(ReplicaExchangeReporter(directory, 1, repex, trajectoryPerState=False, energy=True))
            repex.reporters.append(report)

            # Exercise a short simulation to ensure no runtime errors with the fallback pathway.
            repex.simulate(15)

        del simulation
        del integrator

        # Check the energy file make sure we can load it without stripping the units
        energy = np.loadtxt(os.path.join(directory, 'energy.csv'), delimiter=',').reshape(-1, len(states), len(states))
        for i in range(10):
            for j in range(len(states)):
                for k in range(len(states)):
                    self.assertAlmostEqual(energy[i][j][k], energies[i][j][k] / repex._kT[k])
