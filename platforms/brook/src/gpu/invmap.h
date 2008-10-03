#ifndef __INVMAP_H__
#define __INVMAP_H__

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Mark Friedrichs, Chris Bruns                       *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/*
 * For each atom, calculates the positions at which it's
 * forces are to be picked up from and stores the position
 * in the appropriate index.
 *
 * Input: number of dihedrals, the atom indices, and a flag indicating
 *        whether we're doing i(0), j(1), k(2) or l(3)
 * Output: an array of counts per atom
 *         arrays of inversemaps
 *         nimaps - the number of invmaps actually used.
 * 
 * */
int
gpuCalcInvMap( 
		int posflag,  //0-niatoms-1
		int niatoms,  //3 for angles, 4 for torsions, impropers
		int nints,    //number of interactions
		int natoms,   //number of atoms
		int *atoms, //gromacs interaction list
		int nmaps,      //maximum number of inverse maps
      int counts[],   //output counts of how many places each atom occurs
		float4 *invmaps[], //output array of nmaps inverse maps
		int *nimaps        //output max number of inverse maps actually used
		);

void
gpuPrintInvMaps( int nmaps, int natoms, int counts[], float4 *invmap[],	FILE* logFile );

/* We are still plagued by kernel call overheads. This is for a big fat
 * merged inverse gather kernel:
 * Since we have 32 bit floats, we have 23 bits of mantissa or the largest
 * integer we can represent is 2^23. So it should be quite safe to add 
 * 100000 * n to the index where n is the stream in which we should do the
 * lookup. This assumes that nints < 100000, preferably nints << 100000
 * which should always be true
 * */
int
gpuCalcInvMap_merged( 
		int nints,    //number of interactions
		int natoms,   //number of atoms
		int *atoms,   //ijkl,ijkl,ijkl...
		int nmaps,      //maximum number of inverse maps
        int counts[],   //output counts of how many places each atom occurs
		float4 *invmaps[], //output array of nmaps inverse maps
		int *nimaps        //output max number of inverse maps actually used
		);

/* Repacks the invmap streams for more efficient access in the
 * merged inverse gather kernel
 *
 * buf should be nimaps * natoms large.
 * */
int
gpuRepackInvMap_merged( int natoms, int nmaps, int *counts, 
		float4 *invmaps[], float4 *buf );

#endif //__INVMAP_H__
