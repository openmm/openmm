
/* Portions copyright (c) 2006 Stanford University and Simbios.
 * Contributors: Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef __CpuGrycukInterface_H__
#define __CpuGrycukInterface_H__

#ifdef __cplusplus
#define externC extern "C"
#else
#define externC extern
#endif

#include "SimTKOpenMMRealType.h"
#include <stdio.h>

/**---------------------------------------------------------------------------------------

   Delete the Grycuk associated object(s)

   @return 0 if static CpuGrycuk object was set; else return -1

   --------------------------------------------------------------------------------------- */

externC int cpuDeleteGrycukParameters( void );

/**---------------------------------------------------------------------------------------

	Setup for Grycuk calculations from Gromacs

   @param numberOfAtoms            number of atoms
   @param grycukScaleFactors          array of Grycuk scale factors (one entry each atom)
   @param atomicRadii              atomic radii in Angstrom (one entry each atom)
   @param includeAceApproximation  if true, then include nonpolar 
                                   ACE term in calculations
   @param soluteDielectric         solute dielectric
   @param solventDielectric        solvent dielectric
   @param log                      log reference -- if NULL, then errors/warnings
                                   output to stderr

   The method creates a CpuGrycuk instance -- currently the Grycuk type II model is the
   default (see paper). If the Grycuk type I model is desired change

      GrycukParameters* grycukParameters  = new GrycukParameters( numberOfAtoms, GrycukParameters::GrycukTypeII );
   to
      GrycukParameters* grycukParameters  = new GrycukParameters( numberOfAtoms, GrycukParameters::GrycukTypeI  );

   The created object is a static member of the class CpuGrycuk; 
   when the force routine, cpuCalculateGrycukForces(), is called, 
   the static object is used to compute the forces and energy 

   @return 0

   --------------------------------------------------------------------------------------- */

externC 
int cpuSetGrycukParameters( int numberOfAtoms, RealOpenMM* atomicRadii, RealOpenMM* grycukScaleFactors,
                            int includeAceApproximation, RealOpenMM soluteDielectric,
                             RealOpenMM solventDielectric, FILE* log );

/**---------------------------------------------------------------------------------------

   Get Grycuk scale factors given masses

   @param numberOfAtoms number of atoms
   @param masses        input masses 
   @param scaleFactors  output atomic numbers

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

externC int getGrycukScaleFactorsGivenAtomMasses( int numberOfAtoms, const RealOpenMM* masses,
                                               RealOpenMM* scaleFactors );

/**---------------------------------------------------------------------------------------

   Get Grycuk scale factors given atomic numbers

   @param numberOfAtoms number of atoms
   @param atomicNumber  input atomic number for each atom
   @param scaleFactors  output atomic numbers

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

externC int getGrycukScaleFactors( int numberOfAtoms, const int* atomicNumber,
                                    RealOpenMM* scaleFactors );

#undef externC

#endif
