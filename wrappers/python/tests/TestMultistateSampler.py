from openmm import *
from openmm.unit import *
from openmm.app.internal.multistatesampler import MultistateSampler
import unittest

class TestMultistateSampler(unittest.TestCase):
    def setUp(self):
        # Create a Context for use in tests.  Force groups should have the following energies:
        # 1: 4*k1
        # 2: 2*k2
        # 3: 10
        self.system = System()
        self.system.addParticle(1.0)
        self.system.addParticle(1.0)
        force1 = CustomBondForce('k1*r^2')
        force1.addGlobalParameter('k1', 1.0)
        force1.addBond(0, 1)
        force1.setForceGroup(1)
        self.system.addForce(force1)
        force2 = CustomBondForce('k2*r')
        force2.addGlobalParameter('k2', 1.0)
        force2.addBond(0, 1)
        force2.setForceGroup(2)
        self.system.addForce(force2)
        force3 = CustomExternalForce('10+periodicdistance(x, y, z, 0, 0, 0)')
        force3.addParticle(0)
        force3.setForceGroup(3)
        self.system.addForce(force3)
        self.system.addForce(MonteCarloBarostat(1.0*bar, 300.0*kelvin))
        integrator = LangevinIntegrator(300*kelvin, 1.0/picosecond, 0.004*picoseconds)
        self.context = Context(self.system, integrator, Platform.getPlatform('Reference'))
        self.context.setPositions([Vec3(0, 0, 0), Vec3(0, 2, 0)])

    def validateEnergies(self, sampler, expected):
        """Check that all states have the correct energies when computed in various ways."""
        for i in range(len(sampler.states)):
            energy = sampler.computeEnergy(i)
            self.assertAlmostEqual(expected[i], energy.value_in_unit(kilojoules_per_mole), 5)
        energies = sampler.computeAllEnergies()
        self.assertEqual(len(sampler.states), len(energies))
        for e1, e2 in zip(expected, energies):
            self.assertAlmostEqual(e1, e2.value_in_unit(kilojoules_per_mole), 5)
        relative = sampler.computeRelativeEnergies()
        self.assertEqual(len(sampler.states), len(relative))
        for e1, e2 in zip(expected, relative):
            self.assertAlmostEqual(e1-expected[0], (e2-relative[0]).value_in_unit(kilojoules_per_mole), 6)
        if sampler.groups is None:
            for i in range(len(sampler.states)):
                sampler.applyState(i)
                energy = sampler.context.getState(energy=True).getPotentialEnergy()
                self.assertAlmostEqual(expected[i], energy.value_in_unit(kilojoules_per_mole), 5)

    def testTemperatures(self):
        """Test a set of states that vary only in temperature."""
        states = [{'temperature':t} for t in [300, 350, 400, 450, 500]]
        sampler = MultistateSampler(states, self.context)
        self.assertEqual(1, len(sampler.subsets))
        for i in range(len(states)):
            sampler.applyState(i)
            self.assertAlmostEqual(states[i]['temperature'], self.context.getIntegrator().getTemperature().value_in_unit(kelvin), 5)
            self.assertAlmostEqual(states[i]['temperature'], self.context.getParameter('MonteCarloTemperature'), 5)
        self.validateEnergies(sampler, [16.0]*len(states))

    def testGroups(self):
        """Test a set of states that use different force groups."""
        states = [{'temperature':300, 'groups':{1}},
                  {'temperature':300, 'groups':{1, 2, 3}},
                  {'temperature':500, 'groups':{1:0.5}},
                  {'temperature':500, 'groups':{1:1.0, 2:0.5, 3:0.5}}]
        sampler = MultistateSampler(states, self.context)
        self.assertEqual(1, len(sampler.subsets))
        self.assertEqual(2, len(sampler.groups_of_groups[0]))
        for i in range(len(states)):
            sampler.applyState(i)
            self.assertAlmostEqual(states[i]['temperature'], self.context.getIntegrator().getTemperature().value_in_unit(kelvin), 5)
            self.assertAlmostEqual(states[i]['temperature'], self.context.getParameter('MonteCarloTemperature'), 5)
        self.validateEnergies(sampler, [4.0, 16.0, 2.0, 10.0])

    def testParameters(self):
        """Test a set of states that set parameters to different values."""
        states = [{'k1':1.0, 'k2':1.0},
                  {'k1':1.0, 'k2':2.0},
                  {'k1':2.0, 'k2':1.0}]
        sampler = MultistateSampler(states, self.context)
        self.assertEqual(3, len(sampler.subsets))
        for i in range(len(states)):
            sampler.applyState(i)
            self.assertAlmostEqual(states[i]['k1'], self.context.getParameter('k1'), 5)
            self.assertAlmostEqual(states[i]['k2'], self.context.getParameter('k2'), 5)
        self.validateEnergies(sampler, [16.0, 18.0, 20.0])

    def testParametersAndGroups(self):
        """Test a set of states that differ in both parameters and force groups."""
        states = [{'k1':1.0, 'groups':{1}},
                  {'k1':1.0, 'groups':{1, 2, 3}},
                  {'k1':2.0, 'groups':{1:0.5}},
                  {'k1':2.0, 'groups':{1:1.0, 2:0.5, 3:0.5}}]
        sampler = MultistateSampler(states, self.context)
        self.assertEqual(2, len(sampler.subsets))
        self.assertEqual(2, len(sampler.groups_of_groups[0]))
        self.assertEqual(2, len(sampler.groups_of_groups[1]))
        for i in range(len(states)):
            sampler.applyState(i)
            self.assertAlmostEqual(states[i]['k1'], self.context.getParameter('k1'), 5)
            self.assertAlmostEqual(1.0, self.context.getParameter('k2'), 5)
        self.validateEnergies(sampler, [4.0, 16.0, 4.0, 14.0])

    def testCustomIntegrator(self):
        """Test a set of states that set global variables on a CustomIntegrator."""
        integrator = CustomIntegrator(0.001)
        integrator.addGlobalVariable('temperature', 300.0)
        integrator.addGlobalVariable('friction', 1.0)
        context = Context(self.system, integrator, Platform.getPlatform('Reference'))
        context.setPositions([Vec3(0, 0, 0), Vec3(0, 2, 0)])
        states = [{'temperature':300, 'friction':1.0},
                  {'temperature':300, 'friction':2.0},
                  {'temperature':500, 'friction':1.0},
                  {'temperature':500, 'friction':2.0}]
        sampler = MultistateSampler(states, context)
        self.assertEqual(1, len(sampler.subsets))
        for i in range(len(states)):
            sampler.applyState(i)
            self.assertAlmostEqual(states[i]['temperature'], integrator.getGlobalVariableByName('temperature'), 5)
            self.assertAlmostEqual(states[i]['friction'], integrator.getGlobalVariableByName('friction'), 5)
            self.assertAlmostEqual(states[i]['temperature'], context.getParameter('MonteCarloTemperature'), 5)
        self.validateEnergies(sampler, [16.0]*len(states))
