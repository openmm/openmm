import pickle
import unittest
from simtk.unit import dalton
from simtk.openmm.app import element

class TestElement(unittest.TestCase):
    def test_immutable(self):
        def modifyElement():
            # this should not be allowed
            element.sulfur.mass = 100*dalton
        self.assertRaises(AttributeError, modifyElement)
    
    def test_pickleable(self):
        newsulfur = pickle.loads(pickle.dumps(element.sulfur))
        # make sure that a new object is not created during the pickle/unpickle
        # cycle
        assert element.sulfur == newsulfur
        assert element.sulfur is newsulfur
        assert id(element.sulfur) == id(newsulfur)
    
    def test_attributes(self):
        assert element.hydrogen.atomic_number == 1
        assert element.hydrogen.symbol == 'H'
        assert element.hydrogen.name == 'hydrogen'
        assert element.hydrogen.mass == 1.007947 * dalton


if __name__ == '__main__':
    unittest.main()

        