import unittest

from openmm.app.internal.amoebaforces import AmoebaOutOfPlaneBendForceBuilder as Builder


class TestAmoebaClassifyAngles(unittest.TestCase):
    def test_partial_opbend_generic_is_three_tuple(self):
        angles = [(0, 1, 2)]
        atom_classes = ["P0", "M", "P2", "P3"]
        bonded = [[1], [0, 2, 3], [1], [1]]
        opbend_types = [("P0", "M", "", ""), ("P3", "M", "", "")]
        in_plane, out_of_plane, generic = Builder.classifyAngles(
            angles, atom_classes, bonded, opbend_types
        )
        self.assertEqual(generic, [(0, 1, 2)])
        self.assertEqual(in_plane, [])
        self.assertEqual(out_of_plane, [])

    def test_full_three_partners_in_plane_four_tuple(self):
        angles = [(0, 1, 2)]
        atom_classes = ["P0", "M", "P2", "P3"]
        bonded = [[1], [0, 2, 3], [1], [1]]
        opbend_types = [("P0", "M", "", ""), ("P2", "M", "", ""), ("P3", "M", "", "")]
        in_plane, out_of_plane, generic = Builder.classifyAngles(
            angles, atom_classes, bonded, opbend_types
        )
        self.assertEqual(generic, [])
        self.assertEqual(len(in_plane[0]), 4)
        self.assertEqual(len(out_of_plane), 3)

    def test_covalency_two_generic_three_tuple(self):
        angles = [(0, 1, 2)]
        atom_classes = ["P0", "M", "P2", "P3"]
        bonded = [[1], [0, 2], [1], []]
        opbend_types = [("P0", "M", "", ""), ("P3", "M", "", "")]
        in_plane, out_of_plane, generic = Builder.classifyAngles(
            angles, atom_classes, bonded, opbend_types
        )
        self.assertEqual(generic, [(0, 1, 2)])
        self.assertEqual(in_plane, [])


if __name__ == "__main__":
    unittest.main()
