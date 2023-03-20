__author__ = "Raul P. Pelaez"
import sys
import unittest
from openmm.app import *
from openmm import *
from openmm.unit import *
import tempfile
import openmm.app.element as elem
import numpy as np

if sys.version_info >= (3, 0):
    from io import StringIO
else:
    from cStringIO import StringIO


class TestXtcFile(unittest.TestCase):
    """Test the XTC file parser
    Uses a test file containing 10 frames of a 4-particle system.
    """

    _file = "systems/test.xtc"
    _n_atoms = 4
    _n_frames = 10

    _time = np.array(
        [
            0.7941250801086426,
            0.15734536945819855,
            0.32645168900489807,
            0.5920419692993164,
            0.9089487791061401,
            0.8441296815872192,
            0.815639317035675,
            0.397926390171051,
            0.5505931973457336,
            0.06827980279922485,
        ]
    )

    _frame = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

    _dimensions = (
        np.array(
            [
                [0.79335505, 0.12560561, 0.91993958],
                [0.72167063, 0.27390003, 0.3431609],
                [0.14339106, 0.03506008, 0.003237],
                [0.55738211, 0.92884773, 0.3487983],
                [0.0369907, 0.61924493, 0.95431805],
                [0.47246876, 0.14322655, 0.29019555],
                [0.84334046, 0.79048061, 0.69617069],
                [0.60610449, 0.45382488, 0.18460917],
                [0.21251759, 0.48684791, 0.14755082],
                [0.95540214, 0.72436297, 0.15664689],
            ]
        ).transpose(1, 0)
        * 0.1
    )

    _coordinates = (
        np.array(
            [
                [
                    [0.22561085, 0.57554418, 0.76663136],
                    [0.95949423, 0.68388247, 0.99240756],
                    [0.09906025, 0.35062689, 0.59255189],
                    [0.4234159, 0.08547893, 0.70661318],
                ],
                [
                    [0.65666175, 0.58463252, 0.34902617],
                    [0.87807858, 0.73946381, 0.17238772],
                    [0.01438851, 0.47097427, 0.30818549],
                    [0.57738721, 0.02329309, 0.55919057],
                ],
                [
                    [0.89687312, 0.0954535, 0.93621975],
                    [0.92152476, 0.36477709, 0.00991035],
                    [0.5440616, 0.00308851, 0.12287426],
                    [0.9326995, 0.0297816, 0.947734],
                ],
                [
                    [0.82742298, 0.74183249, 0.37766394],
                    [0.71317816, 0.48780683, 0.59826714],
                    [0.45111823, 0.74698365, 0.08320894],
                    [0.21948391, 0.40945643, 0.19500691],
                ],
                [
                    [0.37078348, 0.32879841, 0.86849266],
                    [0.17665407, 0.45537287, 0.65787739],
                    [0.50271708, 0.91333556, 0.1114412],
                    [0.18044752, 0.63082135, 0.99525344],
                ],
                [
                    [0.11403982, 0.10707579, 0.20506455],
                    [0.85690928, 0.57547241, 0.74563378],
                    [0.60851145, 0.33486977, 0.14469206],
                    [0.63768005, 0.37128955, 0.76609194],
                ],
                [
                    [0.57743162, 0.92428231, 0.44598326],
                    [0.78896332, 0.31820732, 0.27300602],
                    [0.84454942, 0.98303288, 0.98008746],
                    [0.51634634, 0.4087967, 0.48729214],
                ],
                [
                    [0.00800326, 0.2646499, 0.52253413],
                    [0.73950177, 0.94564372, 0.5809049],
                    [0.3633326, 0.44803208, 0.92811674],
                    [0.08597767, 0.78300983, 0.06905601],
                ],
                [
                    [0.63211572, 0.99625397, 0.87334663],
                    [0.71977252, 0.10461613, 0.03891311],
                    [0.30366307, 0.62271225, 0.72741115],
                    [0.13346681, 0.79658383, 0.27278173],
                ],
                [
                    [0.02055479, 0.53108507, 0.43823507],
                    [0.60833883, 0.20549011, 0.22490236],
                    [0.38676259, 0.23006943, 0.07595193],
                    [0.31911924, 0.40371054, 0.84695876],
                ],
            ]
        ).transpose(1, 2, 0)
        * 0.1
    )

    def test_ReadFileMatchesStoredValues(self):
        """Read a test file and make sure it matches the original."""
        xtc = XTCFile(self._file)
        coords, box, time, frames = xtc.read()
        self.assertEqual(coords.shape, self._coordinates.shape)
        self.assertEqual(box.shape, self._dimensions.shape)
        self.assertEqual(time.shape, (self._n_frames,))
        self.assertEqual(frames.shape, (self._n_frames,))
        self.assertTrue(np.allclose(coords, self._coordinates))
        self.assertTrue(np.allclose(box, self._dimensions))
        self.assertTrue(np.allclose(time, self._time))
        self.assertTrue(np.allclose(frames, self._frame))

    def test_CorrectNumberOfFrames(self):
        xtcf = XTCFile(self._file)
        self.assertEqual(xtcf.getNumberOfFrames(), self._n_frames)

    def test_ReadIndividualFrames(self):
        """Read some frames from a file and make sure it matches the original."""
        xtc_some = XTCFile(self._file)
        coords2, box2, time2, step2 = xtc_some.readFrames([2, 8, 9])
        self.assertTrue(np.allclose(coords2, self._coordinates[:, :, [2, 8, 9]]))
        self.assertTrue(np.allclose(box2, self._dimensions[:, [2, 8, 9]]))
        self.assertTrue(np.allclose(time2, self._time[[2, 8, 9]]))
        self.assertTrue(np.allclose(step2, self._frame[[2, 8, 9]]))

    def test_WriteAndReadFileMatches(self):
        """Write a file, read it back, and make sure it matches the original."""
        xtc1 = XTCFile(self._file)
        coords, box, time, step = xtc1.read()
        #Create a temporary file using python tempfile module
        nframes = time.shape[0]
        with tempfile.NamedTemporaryFile() as temp:
            for t in range(nframes):
                xtc2 = XTCFile(temp.name)
                xtc2.writeFrame(np.array(coords[:,:,t]), np.array(box[:,t]), time[t], step[t])

            xtc3 = XTCFile(temp.name)
            coords2, box2, time2, step2 = xtc3.read()
            self.assertEqual(coords.shape, coords2.shape)
            self.assertEqual(box.shape, box2.shape)
            self.assertEqual(time.shape, time2.shape)
            self.assertEqual(step.shape, step2.shape)
            self.assertTrue((time == time2).all())
            self.assertTrue((step == step2).all())
            self.assertTrue((box == box2).all())
            self.assertTrue((coords == coords2).all())


if __name__ == "__main__":
    unittest.main()
