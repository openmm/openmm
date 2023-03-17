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
    """Test the XTC file parser"""

    def test_WriteFile(self):
        """Write a file, read it back, and make sure it matches the original."""
        xtc1 = XTCFile('systems/dopcdope.xtc')
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
if __name__ == '__main__':
    unittest.main()
