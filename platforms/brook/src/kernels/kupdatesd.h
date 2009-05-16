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
