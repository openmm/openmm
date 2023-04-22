# MIT License

# Copyright (c) 2023 Accellera

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Contributors: Stefan Doerr, Raul P. Pelaez

import numpy as np
cimport numpy as np
cimport xtclib
from libcpp.string cimport string
ctypedef np.float32_t FLOAT32_t

def get_xtc_nframes(string filename):
    """
    Get the number of frames in a xtc file.
    Parameters
    ----------
    filename: string
        The filename of the xtc file. You need to pass the string with filename.encode("UTF-8") to this function
    Returns
    -------
    nframes: int
        The number of frames in the xtc file
    """
    return xtclib.xtc_nframes(filename)

def get_xtc_natoms(string filename):
    """
    Get the number of atoms in a xtc file.
    Parameters
    ----------
    filename: string
        The filename of the xtc file. You need to pass the string with filename.encode("UTF-8") to this function
    Returns
    -------
    natoms: int
        The number of atoms in the xtc file
    """
    return xtclib.xtc_natoms(filename)

def read_xtc(string filename):
    """
    Reads a xtc file and return its contents.

    Parameters
    ----------
    filename: string
        The filename of the xtc file. You need to pass the string with filename.encode("UTF-8") to this function
    Returns
    -------
    coords: np.ndarray
        The coordinates of the atoms in the xtc file. Shape: (n_atoms, 3, n_frames)
    box: np.ndarray
        The box vectors of the xtc file. Shape: (3, 3, n_frames)
    time: np.ndarray
        The time of each frame. Shape: (n_frames,)
    step: np.ndarray
        The step of each frame. Shape: (n_frames,)
    """
    cdef int natoms = get_xtc_natoms(filename)
    cdef int nframes = get_xtc_nframes(filename)

    cdef FLOAT32_t[:, :, ::1] coords = np.zeros((natoms, 3, nframes), dtype=np.float32)
    cdef FLOAT32_t[:, :, ::1] box = np.zeros((3, 3, nframes), dtype=np.float32)
    cdef FLOAT32_t[::1] time = np.zeros(nframes, dtype=np.float32)
    cdef int[::1] step = np.zeros(nframes, dtype=np.int32)

    xtclib.xtc_read(
        filename,
        &coords[0, 0, 0],
        &box[0, 0, 0],
        &time[0],
        &step[0],
        natoms,
        nframes,
    )
    return np.asarray(coords), np.asarray(box), np.asarray(time), np.asarray(step)

def xtc_write_frame(string filename, float[:, :] coords, float[:, :] box, float time, int step):
    """
    Appends a single frame to a xtc file (if the file does not exist it is created by this function).
    Parameters
    ----------
    filename: string
        The filename of the xtc file. You need to pass the string with filename.encode("UTF-8") to this function
    coords: np.ndarray
        The coordinates of the atoms in the frame. Shape: (n_atoms, 3)
    box: np.ndarray
        The box vectors of the frame. Shape: (3, 3)
    time: float
        The time of the frame
    step: int
        The step of the frame
    """
    cdef int natoms = coords.shape[0]
    cdef int nframes = 1
    xtclib.xtc_write(
        filename,
        natoms,
        nframes,
        &step,
        &time,
        &coords[0, 0],
        &box[0, 0]
    )


def xtc_rewrite_with_new_timestep(string filename_in, string filename_out,
				  int first_step, int interval, float dt):
    """
    Rewrites a trajectory file with a new timestep and starting step number.
    Parameters
    ----------
    filename_in: string
        The filename of the input xtc file. You need to pass the string with filename.encode("UTF-8") to this function
    filename_out: string
        The filename of the output xtc file. You need to pass the string with filename.encode("UTF-8") to this function
    first_step: int
        The first step to be written to the output file
    interval: int
        The interval between steps to be written to the output file
    dt: float
        The timestep of the output file
    """
    xtclib.xtc_rewrite_with_new_timestep(filename_in, filename_out, first_step, interval, dt)
