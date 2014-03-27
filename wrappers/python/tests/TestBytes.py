import unittest
import simtk.openmm as mm


class TestBytes(unittest.TestCase):
    def test_createCheckpoint(self):
        system = mm.System()
        system.addParticle(1.0)
        refPositions = [(0,0,0)]

        context = mm.Context(system, mm.VerletIntegrator(0))
        context.setPositions(refPositions)
        chk = context.createCheckpoint()
        # check that the return value of createCheckpoint is of type bytes (non-unicode)
        assert isinstance(chk, bytes)

        # set the positions to something random then reload the checkpoint, and
        # make sure that the positions get restored correctly
        context.setPositions([(12345, 12345, 123451)])
        context.loadCheckpoint(chk)
        newPositions = context.getState(getPositions=True).getPositions()._value

        assert newPositions == refPositions

        # try encoding the checkpoint in utf-8. OpenMM should be able to handle this too
        context.setPositions([(12345, 12345, 123451)])
        context.loadCheckpoint(chk.decode('utf-8'))
        newPositions = context.getState(getPositions=True).getPositions()._value

        assert newPositions == refPositions


if __name__ == '__main__':
    unittest.main()

    
