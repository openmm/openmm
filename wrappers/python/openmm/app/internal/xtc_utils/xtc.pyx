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

import os
import warnings
import numpy as np
cimport numpy as np

cimport xtclib

from libc.stdio cimport SEEK_SET, SEEK_CUR
from libc.math cimport ceil

ctypedef np.npy_int64   int64_t
ctypedef np.npy_float32 float32_t
ctypedef np.int32_t INT32_t
ctypedef np.int64_t INT64_t
ctypedef np.uint32_t UINT32_t
ctypedef np.float32_t FLOAT32_t
ctypedef np.float64_t FLOAT64_t

import cython

def get_xtc_nframes(const char* filename):
    """ You need to pass the string with filename.encode("UTF-8") to this function """
    return xtclib.xtc_nframes(filename)

def get_xtc_natoms(const char* filename):
    """ You need to pass the string with filename.encode("UTF-8") to this function """
    return xtclib.xtc_natoms(filename)

def read_xtc(const char* filename):
    """ You need to pass the string with filename.encode("UTF-8") to this function """
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

def xtc_write_frame(const char * filename, float[:, :] coords, float[:, :] box, float time, int step):
    """ You need to pass the string with filename.encode("UTF-8") to this function """
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


def xtc_rewrite_with_new_timestep(const char* filename_in, const char* filename_out,
				  int first_step, int interval, float dt):
    xtclib.xtc_rewrite_with_new_timestep(filename_in, filename_out, first_step, interval, dt)
