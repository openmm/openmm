#ifndef __KFORCE_H__
#define __KFORCE_H__

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

//Define a generic force kernel prototype
//To make switching easier to look at
//This should match the kernels from kforce.br
typedef void (*gpuNBForceFunction)(
		const float natoms,
		const float  dupfac,
		const float  strheight,
		const float  strwidth,
		const float  jstrwidth,
		const float  fstrwidth,
		const float  epsfac,
		const float4 params,
		::brook::stream posq,
		::brook::stream isigeps,
		::brook::stream sigma,
		::brook::stream epsilon,
		::brook::stream excl,
		::brook::stream outforce1,
		::brook::stream outforce2 
		);

void  knbforce_CDLJ (const float  natoms,
		const float  dupfac,
		const float  strheight,
		const float  strwidth,
		const float  jstrwidth,
		const float  fstrwidth,
		const float  epsfac,
		::brook::stream posq,
		::brook::stream isigeps,
		::brook::stream sigma,
		::brook::stream epsilon,
		::brook::stream excl,
		::brook::stream outforce1,
		::brook::stream outforce2); 

void  knbforce_CDLJ4(
      const float  natoms,
      const float  nAtomsCeiling,
		const float  dupfac,
		const float  strheight,
		const float  strwidth,
		const float  jstrwidth,
		const float  fstrwidth,
		const float  epsfac,
		::brook::stream posq,
		::brook::stream charge,
		::brook::stream isigeps,
		::brook::stream sigma,
		::brook::stream epsilon,
		::brook::stream excl,
		::brook::stream outforce1,
		::brook::stream outforce2, 
		::brook::stream outforce3, 
		::brook::stream outforce4); 

void  kmerge_partial_forces (const float  repfac,
		const float  atomStrWidth,
		const float  pforceStrWidth,
		const float  natoms,
		::brook::stream pforce1,
		::brook::stream pforce2,
		::brook::stream force); 

void kMergeFloat3_4( 
      float repfac, 
      float atomStreamWidth, 
      float pStreamWidth,
      float natoms,
      float roundNatoms,
      float iUnroll,
		::brook::stream pforce1,
		::brook::stream pforce2,
		::brook::stream pforce3,
		::brook::stream pforce4,
		::brook::stream force );


void kMergeFloat3_4_nobranch( 
      float repfac, 
      float atomStreamWidth, 
      float pStreamWidth,
      float natoms,
      float roundNatoms,
      float iUnroll,
		::brook::stream pforce1,
		::brook::stream pforce2,
		::brook::stream pforce3,
		::brook::stream pforce4,
		::brook::stream force );

typedef void ( *gpuNBForceFunction14 )
	(const float  natoms,
		const float  strwidth,
		const float  epsfac,
		const float4 params,
		::brook::stream atoms,
		::brook::stream posq,
		::brook::stream sigeps,
		::brook::stream fi,
		::brook::stream fj); 


void  kforce14_CDLJ
	(const float  natoms,
		const float  strwidth,
		const float  epsfac,
		const float4 params,
		::brook::stream atoms,
		::brook::stream posq,
		::brook::stream sigeps,
		::brook::stream fi,
		::brook::stream fj); 

void  kforce14_LDLJ 
	(const float  natoms,
		const float  strwidth,
		const float  epsfac,
		const float4 params,
		::brook::stream atoms,
		::brook::stream posq,
		::brook::stream sigeps,
		::brook::stream fi,
		::brook::stream fj); 

void  kforce14_SFDLJ 
	(const float  natoms,
		const float  strwidth,
		const float  epsfac,
		const float4 params,
		::brook::stream atoms,
		::brook::stream posq,
		::brook::stream sigeps,
		::brook::stream fi,
		::brook::stream fj); 

void  kforce14_SHEFFIELD
	(const float  natoms,
		const float  strwidth,
		const float  epsfac,
		const float4 params,
		::brook::stream atoms,
		::brook::stream posq,
		::brook::stream sigeps,
		::brook::stream fi,
		::brook::stream fj); 

typedef void (*gpuBondedFunction)(
		const float  epsfac,
		const float  xstrwidth,
		const float4 params,
		::brook::stream posq,
		::brook::stream atoms,
		::brook::stream parm0,
		::brook::stream parm1,
		::brook::stream parm2,
		::brook::stream parm3,
		::brook::stream parm4,
		::brook::stream fi,
		::brook::stream fj,
		::brook::stream fk,
		::brook::stream fl);

void  kbonded_CDLJ (const float  epsfac,
		const float  xstrwidth,
		::brook::stream posq,
		::brook::stream charge,
		::brook::stream atoms,
		::brook::stream parm0,
		::brook::stream parm1,
		::brook::stream parm2,
		::brook::stream parm3,
		::brook::stream parm4,
		::brook::stream fi,
		::brook::stream fj,
		::brook::stream fk,
		::brook::stream fl); 

void  kbonded_CDLJDebug (const float  epsfac,
		const float  xstrwidth,
		const float4 params,
		::brook::stream posq,
		::brook::stream atoms,
		::brook::stream parm0,
		::brook::stream parm1,
		::brook::stream parm2,
		::brook::stream parm3,
		::brook::stream parm4,
		::brook::stream fi,
		::brook::stream fj,
		::brook::stream fk,
		::brook::stream fl); 

void  kbonded_LDLJ (const float  epsfac,
		const float  xstrwidth,
		const float4 params,
		::brook::stream posq,
		::brook::stream atoms,
		::brook::stream parm0,
		::brook::stream parm1,
		::brook::stream parm2,
		::brook::stream parm3,
		::brook::stream parm4,
		::brook::stream fi,
		::brook::stream fj,
		::brook::stream fk,
		::brook::stream fl); 

void  kbonded_SFDLJ (const float  epsfac,
		const float  xstrwidth,
		const float4 params,
		::brook::stream posq,
		::brook::stream atoms,
		::brook::stream parm0,
		::brook::stream parm1,
		::brook::stream parm2,
		::brook::stream parm3,
		::brook::stream parm4,
		::brook::stream fi,
		::brook::stream fj,
		::brook::stream fk,
		::brook::stream fl); 

void  kbonded_SHEFFIELD (const float  epsfac,
		const float  xstrwidth,
		const float4 params,
		::brook::stream posq,
		::brook::stream atoms,
		::brook::stream parm0,
		::brook::stream parm1,
		::brook::stream parm2,
		::brook::stream parm3,
		::brook::stream parm4,
		::brook::stream fi,
		::brook::stream fj,
		::brook::stream fk,
		::brook::stream fl); 

void kAddForces3_4( 
		const float  conversion,
      ::brook::stream force1, 
      ::brook::stream force2, 
      ::brook::stream outForce );
 
#endif //__KFORCE_H__
