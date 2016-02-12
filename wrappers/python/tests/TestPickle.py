import unittest
from validateConstraints import *
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import simtk.openmm
import simtk.openmm.app.element as elem
import simtk.openmm.app.forcefield as forcefield
import copy
import pickle

class TestPickle(unittest.TestCase):
    """Pickling / deepcopy of OpenMM objects."""

    def setUp(self):
        """Set up the tests by loading the input pdb files and force field
        xml files.

        """
        # alanine dipeptide with explicit water
        self.pdb1 = PDBFile('systems/alanine-dipeptide-explicit.pdb')
        self.forcefield1 = ForceField('amber99sb.xml', 'tip3p.xml')
        self.topology1 = self.pdb1.topology
        self.topology1.setUnitCellDimensions(Vec3(2, 2, 2))

        # alalnine dipeptide with implicit water
        self.pdb2 = PDBFile('systems/alanine-dipeptide-implicit.pdb')
        self.forcefield2 = ForceField('amber99sb.xml', 'amber99_obc.xml')

    def check_copy(self, object, object_copy):
        """Check that an object's copy is an accurate replica."""
        # Check class name is same.
        self.assertEqual(object.__class__.__name__, object_copy.__class__.__name__)
        # Check serialized contents are the same.
        self.assertEqual(XmlSerializer.serialize(object), XmlSerializer.serialize(object_copy))

    def test_deepcopy(self):
        """Test that serialization/deserialization works (via deepcopy)."""

        # Create system, integrator, and state.
        system = self.forcefield1.createSystem(self.pdb1.topology)
        integrator = VerletIntegrator(2*femtosecond)
        context = Context(system, integrator)
        context.setPositions(self.pdb1.positions)
        state = context.getState(getPositions=True, getForces=True, getEnergy=True)

        #
        # Test deepcopy
        #

        self.check_copy(system, copy.deepcopy(system))
        self.check_copy(integrator, copy.deepcopy(integrator))
        self.check_copy(state, copy.deepcopy(state))
        for force_index in range(system.getNumForces()):
            force = system.getForce(force_index)
            force_copy = copy.deepcopy(force)
            self.check_copy(force, force_copy)

        #
        # Test pickle
        #

        self.check_copy(system, pickle.loads(pickle.dumps(system)))
        self.check_copy(integrator, pickle.loads(pickle.dumps(integrator)))
        self.check_copy(state, pickle.loads(pickle.dumps(state)))
        for force_index in range(system.getNumForces()):
            force = system.getForce(force_index)
            force_copy = pickle.loads(pickle.dumps(force))
            self.check_copy(force, force_copy)

if __name__ == '__main__':
    unittest.main()

