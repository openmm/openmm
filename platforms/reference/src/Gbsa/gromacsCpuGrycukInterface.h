
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

#ifndef __GromacsCpuGrycukInterface_H__
#define __GromacsCpuGrycukInterface_H__

#ifdef __cplusplus
#define externC extern "C"
#else
#define externC extern
#endif

/**---------------------------------------------------------------------------------------

   Setup for Grycuk calculations

   @param top                      Gromacs t_topology (as in md.c)
   @param log                      log reference (stdlog in md.c)
   @param includeAceApproximation  if true, then include nonpolar 
                                   ACE term in calculataions
   @param soluteDielectric         solute dielectric
   @param solventDielectric        solvent dielectric

   Method creates a CpuGrycuk instance w/ Grycuk parameters set

   The created object is a static member of the
   class CpuGrycuk that is referenced to calculate the  implicit solvation
   forces and energy 

  @return 0

   --------------------------------------------------------------------------------------- */

externC 
int gromacsCpuInitialSetupGrycuk( const t_topology* top, FILE* log, int includeAceApproximation,
                                  float soluteDielectric, float solventDielectric );

/**---------------------------------------------------------------------------------------

   Calculate Grycuk forces and energy

   @param atomCoordinates   Gromacs atom coordinates ('x' in md.c)
   @param partialCharges    Gromacs charges ('mdatoms->chargeA' in md.c)
   @param forces            output forces in kJ/mol.A; the computed forces
                            are added to the entries in the array

   Function calls a static method in CpuGrycuk class to calculate forces/energy

   @return 0

   --------------------------------------------------------------------------------------- */

externC
int gromacsCpuCalculateGrycukForces( const rvec* atomCoordinates, const float* partialChargesIn, 
                                     rvec* forces );

#endif
