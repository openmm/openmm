import unittest
import itertools
import simtk.openmm as mm


class TestForceGroups(unittest.TestCase):
    def setUp(self):
        system = mm.System()
        system.addParticle(1.0)

        for i in range(32):
            force = mm.CustomExternalForce(str(i))
            force.addParticle(0, [])
            force.setForceGroup(i)
            system.addForce(force)

        platform = mm.Platform.getPlatformByName('Reference')
        context = mm.Context(system, mm.VerletIntegrator(0), platform)
        context.setPositions([(0,0,0)])
        self.context = context

    def test1(self):
        n = 32
        for (i,j) in itertools.combinations(range(n), 2):
            groups = 1<<i | 1<<j
            e_0 = self.context.getState(getEnergy=True, groups=groups).getPotentialEnergy()._value
            e_1 = self.context.getState(getEnergy=True, groups={i,j}).getPotentialEnergy()._value
            e_ref = i+j
            self.assertEqual(e_0, e_ref)
            self.assertEqual(e_1, e_ref)

    def test2(self):
        with self.assertRaises(TypeError):
            # groups must be an int or set
            self.context.getState(getEnergy=True, groups=(1, 2))

    def test3(self):
        e_0 = self.context.getState(getEnergy=True, groups=-1).getPotentialEnergy()._value
        e_ref = sum(range(32))
        self.assertEqual(e_0, e_ref)


if __name__ == '__main__':
    unittest.main()


