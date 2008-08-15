
/****************************************************************
* This file is part of the gpu acceleration library for gromacs.
* Author: V. Vishal
* Copyright (C) Pande Group, Stanford, 2006
*****************************************************************/
#include <stdio.h>
#include <brook/brook.hpp>
// #include "typedefs.h"
#include "invmap.h"


/*
 * Helper functions for building inverse maps for 
 * torsions, impropers and angles.
 * 
 * */

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
		)
{
	int i, j;
	int atom;
	int mapnum, mapcomp;
	
	for ( i = 0; i < natoms; i++ )
		counts[i] = 0;
	
	for ( i = 0; i < nmaps; i++ ) {
		for ( j = 0; j < natoms; j++ ) {
			invmaps[i][j] = float4( -1.0, -1.0, -1.0, -1.0 );
		}
	}
	
	//This will hold the number of imaps actually used
	*nimaps = -1;

	//Now note down the positions where each atom occurs
	for ( i = 0; i < nints; i++ ) {
		//This is our atom
		atom = atoms[ (niatoms + 1) * i + posflag + 1 ];

		//Special for merged bondeds
		if ( atom == -1 ) {
			continue;
		}

		//Check to make sure we're inside the limits
		if ( counts[atom] > nmaps * 4 ) {
			printf( "Atom %d has too many proper dihedrals(%d, max %d)\n",
			         atom, counts[atom], nmaps * 4 );
			return 0;
		}
		
		//Which invmap will this go into
		mapnum = counts[atom] / 4;

		if ( mapnum > *nimaps )
			*nimaps = mapnum;

		//Which component will it be
		mapcomp = counts[atom] % 4;

		//Set it
		//This is silly, but otherwise I have to declare it as float*
		//and things get even more confusing. :)
		switch (mapcomp) {
			case 0: invmaps[mapnum][atom].x = (float) i; break;
			case 1: invmaps[mapnum][atom].y = (float) i; break;
			case 2: invmaps[mapnum][atom].z = (float) i; break;
			case 3: invmaps[mapnum][atom].w = (float) i; break;
		}
		
		counts[atom]++;
	}

	(*nimaps)++;
	return 1;	
}


void
gpuPrintInvMaps( int nmaps, int natoms, int counts[], float4 *invmap[], FILE* logFile )
{
	int i;
	int j;
	for ( i = 0; i < natoms; i++ ) {
		fprintf( logFile, "%d %d ", i, counts[i] );
		for ( j = 0; j < nmaps; j++ ) {
			fprintf( logFile, "%6.0f %6.0f %6.0f %6.0f", invmap[j][i].x, invmap[j][i].y, 
			         invmap[j][i].z, invmap[j][i].w );
		}
		fprintf( logFile, "\n");
	}
}

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
		)
{
	int i, j;
	int atom;
	int mapnum, mapcomp;
	int pos;
	
	for ( i = 0; i < natoms; i++ )
		counts[i] = 0;

	for ( i = 0; i < nmaps; i++ ) {
		for ( j = 0; j < natoms; j++ ) {
			invmaps[i][j] = float4( -1.0, -1.0, -1.0, -1.0 );
		}
	}

	//This will hold the number of imaps actually used
	*nimaps = -1;

	//For each atom
	for ( i = 0; i < nints; i++ ) {
		for ( j = 0; j < 4; j++ ) {
			
			atom = atoms[ i * 4 + j ];
			
			if ( atom == -1 ) {
				//Nothing to be done for this atom, go to next
				continue;
			}
			
			//Which map
			mapnum = counts[ atom ] / 4;
			
			//Make sure we have space
			if ( mapnum >= nmaps ) {
				printf( "Atom %d has too many bondeds(%d, max %d)\n",
						 atom, counts[atom], nmaps * 4 );
				return 0;
			}
				
			if ( mapnum > *nimaps ) {
				*nimaps = mapnum;
			}

			//Which component
			mapcomp = counts[ atom ] % 4;
			
			//Encode target stream and position
			pos = 100000 * j + i;

			switch ( mapcomp ) {
				case 0: invmaps[mapnum][atom].x = (float) pos; break;
				case 1: invmaps[mapnum][atom].y = (float) pos; break;
				case 2: invmaps[mapnum][atom].z = (float) pos; break;
				case 3: invmaps[mapnum][atom].w = (float) pos; break;
			}

			counts[ atom ]++;

		}
	}
	
	(*nimaps)++;
	return 1;
}

/* Repacks the invmap streams for more efficient access in the
 * merged inverse gather kernel
 *
 * buf should be nimaps * natoms large.
 * */
int
gpuRepackInvMap_merged( int natoms, int nmaps, int *counts, 
		float4 *invmaps[], float4 *buf )
{
	int i, j;
	int nmaps_i;

	for ( i = 0; i < natoms; i++ ) {
		for ( j = 0; j < nmaps; j++ ) {
			buf[ i + j*natoms ] = float4( -1.0f, -1.0f, -1.0f, -1.0f );
		}
	}
	
	for ( i = 0; i < natoms; i++ ) {
		
		nmaps_i = counts[i] / 4;
		if ( counts[i] % 4 ) 
			nmaps_i++;
		
		for ( j = 0; j < nmaps_i; j++ ) {
			buf[ i + j * natoms ] = invmaps[j][i];
		}
	}
	return 1;
}
