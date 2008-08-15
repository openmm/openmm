#ifndef __INVMAP_H__
#define __INVMAP_H__

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
