import unittest
import simtk.openmm as mm


class TestBytes(unittest.TestCase):
    def test_createCheckpoint(self):
        # check that the return value of createCheckpoint is of type bytes (non-unicode)
        system = mm.System()
        system.addParticle(1.0)
        mm.Context(system, mm.VerletIntegrator(0))
        chk = mm.Context(system, mm.VerletIntegrator(0)).createCheckpoint()
        assert isinstance(chk, bytes)

if __name__ == '__main__':
    unittest.main()
