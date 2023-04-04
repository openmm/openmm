/*
MIT License

Copyright (c) 2023 Accellera

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Contributors: Stefan Doerr, Raul P. Pelaez
*/
#ifdef PLATFORM_Linux
    #if defined(__i386__) || defined(__x86_64__)
__asm__(".symver memcpy,memcpy@GLIBC_2.2.5");
    #endif
#endif
#ifndef XTC
    #define XTC 1
    #include "xdrfile.h"
    #ifdef __cplusplus
extern "C" {
    #endif

// Reads a nframes from a trajectory file starting at "frame". Opens the
// file, reads the frame and closes the file. Puts the read frame into
// the arrays, at position fidx.
void xtc_read_frame(char* filename, float* coords_arr, float* box_arr, float* time_arr, int* step_arr, int natoms, int frame, int nframes, int fidx);

// Writes a full trajectory to a file. Opens the file with append, write nframes and closes
// the file.
int xtc_write(char* filename, int natoms, int nframes, int* step, float* timex, float* pos, float* box);

struct XTC_frame {
    float box[9];
    int natoms;
    unsigned long step;
    double time;
    float* pos;
}; // XTC_frame;

// Read all frames from an XTC file,  returns a pointer to an array of
// XTC_frame structs for each frame in the file. Argument pointers are
// filled with  the number  of atoms, number  of frames,  timestep and
// stepsize found in the file.
struct XTC_frame* xtc_read(char* filename, int* natoms, int* Nframes, double* dt, int* dstep);

// Reads nframes from  an xtc file, fills the arrays  in the input with
// the data. Each array must be allocated with the correct size.
void xtc_read_new(char* filename, float* coords_arr, float* box_arr, float* time_arr, int* step_arr, int natoms, int nframes);

// Get the number of frames in a trajectory file
int xtc_nframes(char* filename);

// Get the number of atoms in a trajectory file
int xtc_natoms(char* filename);

// Open an XTC file for writing as append, returning an XDRFILE pointer
XDRFILE* xtc_open_for_write(char* filename);

    #ifdef __cplusplus
}
    #endif

#endif
