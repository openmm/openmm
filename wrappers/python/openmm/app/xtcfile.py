
from __future__ import print_function, division, absolute_import
__author__ = "Raul P. Pelaez"
__version__ = "1.0"

from openmm.app.internal.xtc_utils import read_xtc, xtc_write_frame
import numpy as np
import os

class XTCFile(object):

    """Parses an XTC file and provides access to the coordinates, box
    vectors, time and step of each  frame.
    This class also provides a method  to write new frames into an XTC
    file from a set of coordinates, box vectors, time and step.
    """

    def __init__(self, filename):
        """Create an XTCFile object.
        Parameters
        ----------
        filename : str
            The name of the XTC file to open.
        """
        self._filename = filename

    def read(self):
        """Reads the XTC file and returns the coordinates, box vectors,
        time and step of each frame.
        Returns
        -------
        coords : np.array
            The coordinates of the atoms in the frame. Shape (n_frames, n_atoms, 3)
        box : np.array
            The box vectors of the frame. Shape (n_frames, 3)
        time : np.array
            The time of the frame. Shape (n_frames,)
        step : np.array
            The step of the frame. Shape (n_frames,)
        """
        if not os.path.isfile(self._filename):
            raise FileNotFoundError(f"The file '{self._filename}' does not exist.")
        _coords, _box, _time, _step = read_xtc(self._filename.encode('utf-8'))
        return _coords, _box, _time, _step

    def writeFrame(self, coords, box, time, step):
        """Writes a new frame into the XTC file.
        Parameters
        ----------
        coords : np.array
            The coordinates of the atoms in the frame. Shape (n_atoms, 3)
        box : np.array
            The box vectors of the frame. Shape (3,)
        time : float
            The time of the frame.
        step : int
            The step of the frame.
        """
        xtc_write_frame(self._filename.encode('utf-8'), np.array(coords), np.array(box), time, step)