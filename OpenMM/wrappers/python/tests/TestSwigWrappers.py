import unittest
import simtk.openmm as mm


class TestSwigWrappers(unittest.TestCase):
    def test_1(self):
        # This tests for a refcounting bug in the swig wrappers
        # that was previously problematic.
        # See https://github.com/pandegroup/openmm/issues/1214
        for cycle in range(10):
            system = mm.System()
            system.getDefaultPeriodicBoxVectors()

if __name__ == '__main__':
    unittest.main()
