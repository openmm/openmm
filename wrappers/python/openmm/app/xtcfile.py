
from __future__ import print_function, division, absolute_import
__author__ = "Raul P. Pelaez"
__version__ = "1.0"

from openmm.app.internal.xtc_utils import read_xtc, xtc_write_frames, xtc_open_for_write, xdrfile_close

import os

class XTCTrajectoryFile(object):
    def __init__(self, filename):
        self._filename = filename

    def read(self):
        if not os.path.isfile(self._filename):
            raise FileNotFoundError(f"The file '{self._filename}' does not exist.")
        self._coords, self._box, self._time, self._step = read_xtc(self._filename.encode('utf-8'))
        return self._coords, self._box, self._time, self._step

    def write_frame(self, coords, box, time, step):
        if not hasattr(self, '_file'):
            self._file = xtc_open_for_write(self._filename.encode('utf-8'))
        xtc_write_frames(self._file, coords, box, time, step)

    def close(self):
        if hasattr(self, '_file'):
            xdrfile_close(self._file)
