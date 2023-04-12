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

cimport numpy as np
ctypedef np.npy_int64 int64_t
ctypedef np.npy_float32 float32_t

cdef extern from "include/xtc.h":
    int xtc_nframes(const char *filename)
    int xtc_natoms(const char *filename)
    void xtc_read(const char *filename, float *coords_arr, float *box_arr, float *time_arr, int *step_arr, int natoms, int nframes)
    void xtc_write(const char *filename, int natoms, int nframes, int *step, float *timex, float *pos, float *box)
