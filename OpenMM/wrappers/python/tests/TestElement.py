import pickle
import random
from simtk.unit import dalton, is_quantity
from simtk.openmm.app import element
import unittest

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
        self.assertEqual(element.sulfur, newsulfur)
        self.assertTrue(element.sulfur is newsulfur)
    
    def test_attributes(self):
        self.assertEqual(element.hydrogen.atomic_number, 1)
        self.assertEqual(element.hydrogen.symbol, 'H')
        self.assertEqual(element.hydrogen.name, 'hydrogen')
        self.assertEqual(element.hydrogen.mass, 1.007947 * dalton)

    def test_getByMass(self):
        """ Tests the getByMass method """
        def exhaustive_search(mass):
            """
            Searches through all element symbols and finds the one with the
            smallest mass difference
            """
            min_diff = mass
            closest_element = None
            for elem in sorted(element.Element._elements_by_symbol.values(),
                               key=lambda x:x.mass):
                diff = abs(elem.mass._value - mass)
                if diff < min_diff:
                    min_diff = diff
                    closest_element = elem
            return closest_element

        # Check 500 random numbers between 0 and 200
        for i in range(500):
            mass = random.random() * 200
            elem = element.Element.getByMass(mass)
            self.assertTrue(elem is exhaustive_search(mass))


if __name__ == '__main__':
    unittest.main()
