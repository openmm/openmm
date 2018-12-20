""" Tests the Vec3 object """
from unittest import TestCase
from simtk.openmm import Vec3

class TestVectors(TestCase):
    """ Tests the Vec3 type """

    def testVec3Attributes(self):
        vec1 = Vec3(1, 2, 3)
        self.assertEqual(vec1.x, 1)
        self.assertEqual(vec1.y, 2)
        self.assertEqual(vec1.z, 3)

    def testNegation(self):
        vec1 = Vec3(1, 2, 3)
        vec1_neg = Vec3(-1, -2, -3)
        self.assertEqual(-vec1, vec1_neg)

    def testVec3Equality(self):
        vec1 = Vec3(1, 2, 3)
        vec2 = Vec3(1, 2, 3)
        self.assertEqual(vec1, vec2)

    def testVec3Addition(self):
        vec1 = Vec3(1, 2, 3)
        vec2 = Vec3(4, 5, 6)
        vec2_tup = (4, 5, 6)
        result = Vec3(5, 7, 9)
        self.assertEqual(vec1 + vec2, result)
        self.assertEqual(vec1 + vec2_tup, result)

    def testVec3Subtraction(self):
        vec1 = Vec3(1, 2, 3)
        vec2 = Vec3(3, 2, 1)
        vec2_tup = (3, 2, 1)
        result = Vec3(-2, 0, 2)
        self.assertEqual(vec1 - vec2, result)
        self.assertEqual(vec1 - vec2_tup, result)
        self.assertEqual(vec2_tup - vec1, -result)

    def testVec3Multiplication(self):
        vec1 = Vec3(1, 2, 3)
        factor = 2
        result = Vec3(2, 4, 6)
        self.assertEqual(vec1 * factor, result)
        self.assertEqual(factor * vec1, result)

    def testVec3Division(self):
        vec1 = Vec3(4, 5, 6)
        factor = 2
        result = Vec3(2, 2.5, 3)
        self.assertEqual(vec1 / factor, result)
        with self.assertRaises(TypeError):
            2 / vec1
