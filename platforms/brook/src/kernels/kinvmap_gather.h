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

void  kinvmap_gather (const float  strwidth,
		::brook::stream invmap,
		::brook::stream forces,
		::brook::stream inforce,
		::brook::stream outforce); 

void  kinvmap_gather2 (const float  strwidth,
		::brook::stream invmap1,
		::brook::stream invmap2,
		::brook::stream forces,
		::brook::stream inforce,
		::brook::stream outforce); 

void  kinvmap_gather3 (const float  strwidth,
		::brook::stream invmap1,
		::brook::stream invmap2,
		::brook::stream invmap3,
		::brook::stream forces,
		::brook::stream inforce,
		::brook::stream outforce);

void  kinvmap_gather4 (const float  strwidth,
		::brook::stream invmap1,
		::brook::stream invmap2,
		::brook::stream invmap3,
		::brook::stream invmap4,
		::brook::stream forces,
		::brook::stream inforce,
		::brook::stream outforce);

void  kinvmap_gather5 (const float  strwidth,
		::brook::stream invmap1,
		::brook::stream invmap2,
		::brook::stream invmap3,
		::brook::stream invmap4,
		::brook::stream invmap5,
		::brook::stream forces,
		::brook::stream inforce,
		::brook::stream outforce);

void  kinvmap_gather6 (const float  strwidth,
		::brook::stream invmap1,
		::brook::stream invmap2,
		::brook::stream invmap3,
		::brook::stream invmap4,
		::brook::stream invmap5,
		::brook::stream invmap6,
		::brook::stream forces,
		::brook::stream inforce,
		::brook::stream outforce); 


void  kinvmap_gather2_2 (const float  strwidth,
		::brook::stream invmap3_1,
		::brook::stream invmap3_2,
		::brook::stream forces3,
		::brook::stream invmap4_1,
		::brook::stream invmap4_2,
		::brook::stream forces4,
		::brook::stream inforce,
		::brook::stream outforce);

void  kinvmap_gather2_1 (const float  strwidth,
		::brook::stream invmap3_1,
		::brook::stream invmap3_2,
		::brook::stream forces3,
		::brook::stream invmap4_1,
		::brook::stream forces4,
		::brook::stream inforce,
		::brook::stream outforce);

void  kinvmap_gather2_3 (const float  strwidth,
		::brook::stream invmap3_1,
		::brook::stream invmap3_2,
		::brook::stream forces3,
		::brook::stream invmap4_1,
		::brook::stream invmap4_2,
		::brook::stream invmap4_3,
		::brook::stream forces4,
		::brook::stream inforce,
		::brook::stream outforce);

void  kinvmap_gather1_2 (const float  strwidth,
		::brook::stream invmap3_1,
		::brook::stream forces3,
		::brook::stream invmap4_1,
		::brook::stream invmap4_2,
		::brook::stream forces4,
		::brook::stream inforce,
		::brook::stream outforce);

void  kinvmap_gather2_4 (const float  strwidth,
		::brook::stream invmap3_1,
		::brook::stream invmap3_2,
		::brook::stream forces3,
		::brook::stream invmap4_1,
		::brook::stream invmap4_2,
		::brook::stream invmap4_3,
		::brook::stream invmap4_4,
		::brook::stream forces4,
		::brook::stream inforce,
		::brook::stream outforce);

void  kinvmap_gather2_5 (const float  strwidth,
		::brook::stream invmap3_1,
		::brook::stream invmap3_2,
		::brook::stream forces3,
		::brook::stream invmap4_1,
		::brook::stream invmap4_2,
		::brook::stream invmap4_3,
		::brook::stream invmap4_4,
		::brook::stream invmap4_5,
		::brook::stream forces4,
		::brook::stream inforce,
		::brook::stream outforce);

void  kinvmap_gather3_2 (const float  strwidth,
		::brook::stream invmap3_1,
		::brook::stream invmap3_2,
		::brook::stream invmap3_3,
		::brook::stream forces3,
		::brook::stream invmap4_1,
		::brook::stream invmap4_2,
		::brook::stream forces4,
		::brook::stream inforce,
		::brook::stream outforce);

void  kinvmap_gather3_1 (const float  strwidth,
		::brook::stream invmap3_1,
		::brook::stream invmap3_2,
		::brook::stream invmap3_3,
		::brook::stream forces3,
		::brook::stream invmap4_1,
		::brook::stream forces4,
		::brook::stream inforce,
		::brook::stream outforce);

void  kinvmap_gather3_3 (const float  strwidth,
		::brook::stream invmap3_1,
		::brook::stream invmap3_2,
		::brook::stream invmap3_3,
		::brook::stream forces3,
		::brook::stream invmap4_1,
		::brook::stream invmap4_2,
		::brook::stream invmap4_3,
		::brook::stream forces4,
		::brook::stream inforce,
		::brook::stream outforce);

void  kinvmap_gather3_4 (const float  strwidth,
		::brook::stream invmap3_1,
		::brook::stream invmap3_2,
		::brook::stream invmap3_3,
		::brook::stream forces3,
		::brook::stream invmap4_1,
		::brook::stream invmap4_2,
		::brook::stream invmap4_3,
		::brook::stream invmap4_4,
		::brook::stream forces4,
		::brook::stream inforce,
		::brook::stream outforce);

void  kinvmap_gather3_5 (const float  strwidth,
		::brook::stream invmap3_1,
		::brook::stream invmap3_2,
		::brook::stream invmap3_3,
		::brook::stream forces3,
		::brook::stream invmap5_1,
		::brook::stream invmap5_2,
		::brook::stream invmap5_3,
		::brook::stream invmap5_4,
		::brook::stream invmap5_5,
		::brook::stream forces5,
		::brook::stream inforce,
		::brook::stream outforce);


void  kinvmap_gather1_1 (const float  strwidth,
		::brook::stream invmap3_1,
		::brook::stream forces3,
		::brook::stream invmap5_1,
		::brook::stream forces5,
		::brook::stream inforce,
		::brook::stream outforce);


void  kinvmap_gather5_2 (const float  strwidth,
		::brook::stream invmap5_1,
		::brook::stream invmap5_2,
		::brook::stream invmap5_3,
		::brook::stream invmap5_4,
		::brook::stream invmap5_5,
		::brook::stream forces5,
		::brook::stream invmap2_1,
		::brook::stream invmap2_2,
		::brook::stream forces2,
		::brook::stream inforce,
		::brook::stream outforce);


void  kinvmap_gather5_3 (const float  strwidth,
		::brook::stream invmap5_1,
		::brook::stream invmap5_2,
		::brook::stream invmap5_3,
		::brook::stream invmap5_4,
		::brook::stream invmap5_5,
		::brook::stream forces5,
		::brook::stream invmap2_1,
		::brook::stream invmap2_2,
		::brook::stream invmap2_3,
		::brook::stream forces2,
		::brook::stream inforce,
		::brook::stream outforce);


void  kinvmap_gather4_3 (const float  strwidth,
		::brook::stream invmap5_1,
		::brook::stream invmap5_2,
		::brook::stream invmap5_3,
		::brook::stream invmap5_4,
		::brook::stream forces5,
		::brook::stream invmap2_1,
		::brook::stream invmap2_2,
		::brook::stream invmap2_3,
		::brook::stream forces2,
		::brook::stream inforce,
		::brook::stream outforce);


void  kinvmap_gather4_2 (const float  strwidth,
		::brook::stream invmap4_1,
		::brook::stream invmap4_2,
		::brook::stream invmap4_3,
		::brook::stream invmap4_4,
		::brook::stream forces4,
		::brook::stream invmap2_1,
		::brook::stream invmap2_2,
		::brook::stream forces2,
		::brook::stream inforce,
		::brook::stream outforce);


void  kinvmap_gather_merged (const float  natoms,
		const float  strwidth,
		const float  istrwidth,
		const float  fstrwidth,
		::brook::stream nimap,
		::brook::stream invmap,
		::brook::stream fi,
		::brook::stream fj,
		::brook::stream fk,
		::brook::stream fl,
		::brook::stream inforce,
		::brook::stream outforce); 

void  kinvmap_gather_merged9 (const float  natoms,
		const float  strwidth,
		::brook::stream invmap0,
		::brook::stream invmap1,
		::brook::stream invmap2,
		::brook::stream invmap3,
		::brook::stream invmap4,
		::brook::stream invmap5,
		::brook::stream invmap6,
		::brook::stream invmap7,
		::brook::stream invmap8,
		::brook::stream fi,
		::brook::stream fj,
		::brook::stream fk,
		::brook::stream fl,
		::brook::stream inforce,
		::brook::stream outforce); 

void  kinvmap_gather_merged5 (const float  natoms,
		const float  strwidth,
		::brook::stream invmap0,
		::brook::stream invmap1,
		::brook::stream invmap2,
		::brook::stream invmap3,
		::brook::stream invmap4,
		::brook::stream fi,
		::brook::stream fj,
		::brook::stream fk,
		::brook::stream fl,
		::brook::stream inforce,
		::brook::stream outforce); 
