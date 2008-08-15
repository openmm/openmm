#ifndef __KFORCE_H__
#define __KFORCE_H__

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
		const float4 params,
		::brook::stream posq,
		::brook::stream isigeps,
		::brook::stream sigma,
		::brook::stream epsilon,
		::brook::stream excl,
		::brook::stream outforce1,
		::brook::stream outforce2); 

void  knbforce_CDLJ2(const float  natoms,
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
		const float4 params,
		::brook::stream posq,
		::brook::stream isigeps,
		::brook::stream sigma,
		::brook::stream epsilon,
		::brook::stream excl,
		::brook::stream outforce1,
		::brook::stream outforce2, 
		::brook::stream outforce3, 
		::brook::stream outforce4); 

void  knbforce_CDLJ4Debug(
      const float  natoms,
      const float  nAtomsCeiling,
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
		::brook::stream outforce2, 
		::brook::stream outforce3, 
		::brook::stream outforce4); 

void  knbforce_CDLJ_1(
      const float  natoms,
      const float  nAtomsCeiling,
		const float  strwidth,
		const float  epsfac,
		::brook::stream posq,
		::brook::stream isigeps,
		::brook::stream excl,
		::brook::stream outforce ); 

void  knbforce_CDLJ4NoEx(
      const float  natoms,
      const float  nAtomsCeiling,
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
		::brook::stream outforce2, 
		::brook::stream outforce3, 
		::brook::stream outforce4); 

void  knbforce_LDLJ (const float  natoms,
		const float  dupfac,
		const float  strheight,
		const float  strwidth,
		const float  jstrwidth,
		const float  fstrwidth,
		const float  epsfac,
		const float4 params,
		const __BRTIter& a_iatom2,
		::brook::stream posq,
		::brook::stream isigeps,
		::brook::stream sigma,
		::brook::stream epsilon,
		::brook::stream excl,
		::brook::stream outforce1,
		::brook::stream outforce2); 

void  knbforce_SFDLJ (const float  natoms,
		const float  dupfac,
		const float  strheight,
		const float  strwidth,
		const float  jstrwidth,
		const float  fstrwidth,
		const float  epsfac,
		const float4 params,
		const __BRTIter& a_iatom2,
		::brook::stream posq,
		::brook::stream isigeps,
		::brook::stream sigma,
		::brook::stream epsilon,
		::brook::stream excl,
		::brook::stream outforce1,
		::brook::stream outforce2); 

void  knbforce_SHEFFIELD (const float  natoms,
		const float  dupfac,
		const float  strheight,
		const float  strwidth,
		const float  jstrwidth,
		const float  fstrwidth,
		const float  epsfac,
		const float4 params,
		const __BRTIter& a_iatom2,
		::brook::stream posq,
		::brook::stream isigeps,
		::brook::stream sigma,
		::brook::stream epsilon,
		::brook::stream excl,
		::brook::stream outforce1,
		::brook::stream outforce2); 

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
