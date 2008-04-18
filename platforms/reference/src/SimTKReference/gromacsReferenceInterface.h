
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

#ifndef __GromacsReferenceInterface_H__
#define __GromacsReferenceInterface_H__

#ifdef __cplusplus
#define externC extern "C"
#else
#define externC extern
#endif

#include "../SimTKUtilities/SimTKOpenMMCommon.h"

/**---------------------------------------------------------------------------------------

   Calculate forces for given configuration and topology

   @param top                      Gromacs t_topology data struct
   @param gromacAtomCoordinates    atom configuration
   @param partialChargesIn         array of partial charges
   @param fr                       Gromac's force_record?????
   @param log                      log reference (stdlog in md.c)

   @return 0

   --------------------------------------------------------------------------------------- */

externC
int gromacsReferenceForce( const t_topology* top, const rvec* gromacAtomCoordinates,
                           const t_mdatoms *md, t_forcerec* fr,
                           char* baseFileName, FILE* log );

externC
int gromacsReferenceSdUpdate( int numberOfAtoms, const t_topology* top, const rvec* gromacAtomCoordinates,
                              const t_mdatoms *md, t_forcerec* fr,
                              const rvec* gromacsVelocities, const rvec* gromacsForces,
                              real* masses, real deltaT, real tau, real temperature,
                              char* baseFileName, FILE* log );

/**---------------------------------------------------------------------------------------

   Calculate forces for given configuration and topology
   Used in self-test; parameter arrays, ... only allocated once (static); freed
   if debug is -1

   @param top                      Gromacs t_topology data struct
   @param gromacAtomCoordinates    atom configuration
   @param md                       Gromacs t_mdatoms data struct
   @param fr                       Gromacs t_forcerec data struct
   @param includeNonBonded         include nonbonded ixn
   @param baseFileName             Base file name
   @param debug                    debug flag (== -1, free arrays and return)
   @param log                      log reference (stdlog in md.c)

   @return 0

   --------------------------------------------------------------------------------------- */

externC
int gromacsReferenceForceSelfTest( 
                                   const t_topology* top, const rvec* gromacAtomCoordinates,
                                   const t_mdatoms *md, t_forcerec* fr,
                                   rvec* forces, int includeNonBonded, char* baseFileName,
                                   int debug, FILE* log );

/**---------------------------------------------------------------------------------------

   Remove linear momentum from velocities

   @param numberOfAtoms            number of atoms
   @param gromacsMasses            masses
   @param gromacsVelocities        velocities
   @param baseFileName             base file name
   @param log                      log reference (stdlog in md.c)

   @return 0

   --------------------------------------------------------------------------------------- */

externC
int gromacsReferenceRemoveLinearMomentum( int numberOfAtoms, const real* gromacsMasses,
                                          rvec* gromacsVelocities, char* baseFileName, FILE* log );

#endif
