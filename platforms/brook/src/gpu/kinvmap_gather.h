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
