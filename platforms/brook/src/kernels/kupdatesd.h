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

void  kupdate_sd1 (
		const float  xstrwidth,
		const float  gstrwidth,
		const float  goffset,
		const float  cem,
		const float  pc1,
		const float  pc2,
		const float  pc3,
		::brook::stream sdpc,
		::brook::stream fgauss,
		::brook::stream sd2X,
		::brook::stream posq,
		::brook::stream f,
		::brook::stream v,
		::brook::stream invmass,
		::brook::stream sd1V,
		::brook::stream vnew,
		::brook::stream posqp); 

void  kupdate_sd2 (
		const float  xstrwidth,
		const float  gstrwidth,
		const float  goffset,
		const float  pc1,
		const float  pc2,
		::brook::stream sdpc,
		::brook::stream fgauss,
		::brook::stream sd1V,
		::brook::stream posq,
		::brook::stream posqp,
		::brook::stream vnew,
		::brook::stream sd2X,
		::brook::stream v,
		::brook::stream posqp2); 

void  kpermute_vectors (const float  gstrwidth,
		::brook::stream perm,
		::brook::stream gvin,
		::brook::stream gvout); 

void  kupdate_sd2_fix1 (const float  xstrwidth,
		const float  gstrwidth,
		const float  goffset,
		const float  pc1,
		const float  pc2,
		::brook::stream sdpc,
		::brook::stream fgauss,
		::brook::stream sd1V,
		::brook::stream posq,
		::brook::stream posqp,
		::brook::stream vnew,
		::brook::stream sd2X,
		::brook::stream v,
		::brook::stream posqp2); 

void  kupdate_sd1_fix1 (const float  xstrwidth,
		const float  gstrwidth,
		const float  goffset,
		const float  cem,
		const float  pc1,
		const float  pc2,
		const float  pc3,
		::brook::stream sdpc,
		::brook::stream fgauss,
		::brook::stream sd2X,
		::brook::stream posq,
		::brook::stream f,
		::brook::stream v,
		::brook::stream invmass,
		::brook::stream sd1V,
		::brook::stream vnew,
		::brook::stream posqp);

void  kupdate_sd2_fix1_FixedRV(const float  xstrwidth,
		const float  gstrwidth,
		const float  goffset,
		const float  pc1,
		const float  pc2,
		::brook::stream sdpc,
		::brook::stream fgauss,
		::brook::stream sd1V,
		::brook::stream posq,
		::brook::stream posqp,
		::brook::stream vnew,
		::brook::stream sd2X,
		::brook::stream v,
		::brook::stream posqp2); 

void  kupdate_sd1_fix1_FixedRV(const float  xstrwidth,
		const float  gstrwidth,
		const float  goffset,
		const float  cem,
		const float  pc1,
		const float  pc2,
		const float  pc3,
		::brook::stream sdpc,
		::brook::stream fgauss,
		::brook::stream sd2X,
		::brook::stream posq,
		::brook::stream f,
		::brook::stream v,
		::brook::stream invmass,
		::brook::stream sd1V,
		::brook::stream vnew,
		::brook::stream posqp);
