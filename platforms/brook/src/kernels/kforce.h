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
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Mark Friedrichs, Mike Houston                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
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
