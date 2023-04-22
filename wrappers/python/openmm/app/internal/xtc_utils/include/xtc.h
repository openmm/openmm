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
#ifndef XTC
#define XTC
#include "xdrfile.h"
#include<string>
// Get the number of frames in a trajectory file
int xtc_nframes(std::string filename);

// Get the number of atoms in a trajectory file
int xtc_natoms(std::string filename);

// Reads nframes from  an xtc file, fills the arrays  in the input with
// the data. Each array must be allocated with the correct size.
void xtc_read(std::string filename, float* coords_arr, float* box_arr, float* time_arr, int* step_arr, int natoms, int nframes);

// Appends a trajectory to a file.
void xtc_write(std::string filename, int natoms, int nframes, int* step, float* timex, float* pos, float* box);

// Rewrites a trajectory file with a new timestep and starting step number.
// Useful when the step number is larger than 2^32.
void xtc_rewrite_with_new_timestep(std::string filename_in, std::string filename_out, int first_step, int interval, float dt);
#endif
