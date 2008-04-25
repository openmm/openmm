
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

#include "CpuGbsa.h"

// 'gmxGbsa.h' has declarations for cpuSetGbsaParameters() and
// and cpuCalculateGbsaForces()

#include "gmxGbsa.h"

/**---------------------------------------------------------------------------------------

	Example calling sequence:

		float** forces = NULL;
      float energy;
      cpuSetGbsaParameters( mdatoms->nr, top, stdlog );

      // if 'forces' array is initially null, then memory will
      // be allocated by cpuCalculateGbsaForces()

      // format for 'forces' array:
		//		forces[3][numberAtoms]

      forces = cpuCalculateGbsaForces( x, mdatoms->chargeA, &energy,
                                      forces, stdlog );

   --------------------------------------------------------------------------------------- */

/**---------------------------------------------------------------------------------------

	Setup for Gbsa calculations

   @param numberOfAtoms            number of atoms
   @param top                      Gromacs t_topology (as in md.c)
   @param log                      log reference (stdlog in md.c)
   @param includeAceApproximation  if true, then include nonpolar 
                                   ACE term in calculataions
   @param innerDielectric          nonsolvent dielectric
   @param solventDielectric        solvent dielectric

   function creates a CpuGbsa instance; the created object is a static member of the
   class CpuGbsa; when the force routine (a static method) is called, the object is
   used to compute the forces and energy 

  @return 0

   --------------------------------------------------------------------------------------- */

extern "C" int 
cpuSetGbsaParameters( int numberOfAtoms, const t_topology* top, FILE* log,
                      bool includeAceApproximation,
                      float innerDielectric, float solventDielectric ){

   // ---------------------------------------------------------------------------------------

   static bool printGbsaParameters           = false;
   static const char* methodName             = "\ncpuSetGbsaParameters: ";

   // static bool debug                         = true;
   static bool debug                         = false;

   // ---------------------------------------------------------------------------------------
   
   // toggle whether logging enabled

   FILE* cpuGbsaLog               = debug ? log : NULL;
   CpuGbsa* cpuGbsa               = new CpuGbsa( numberOfAtoms, top, cpuGbsaLog );

   // include/not include ACE approximation (nonpolar solvation)

   cpuGbsa->setIncludeAceApproximation( includeAceApproximation );
   // cpuGbsa->setIncludeAceApproximation( false );
   // cpuGbsa->setIncludeAceApproximation( true );

   // default values are      1.0/1.0      -> gives values in 10.0*kcal/A (force) and 100*kcal/A (energy)
   //      if values are 0.41868/0.41868   -> gives values in Gromacs units kJ/nm

   if( !debug ){
      cpuGbsa->setForceConversionFactor(  0.41868f );
      cpuGbsa->setEnergyConversionFactor( 0.41868f );
   }

   // dielectrics

   cpuGbsa->getGbsaParameters()->setSolventDielectric( solventDielectric );
   cpuGbsa->getGbsaParameters()->setInnerDielectric( innerDielectric );

   // set static member for subsequent calls to calculate forces/energy 

   CpuGbsa::setCpuGbsa( cpuGbsa );

   // ---------------------------------------------------------------------------------------

   if( log != NULL ){
      cpuGbsa->logState( methodName, log );
      (void) fprintf( log, "%s done", methodName );
      (void) fflush( log );
   }

   // ---------------------------------------------------------------------------------------

   return 0;
}

/**---------------------------------------------------------------------------------------

   Calculate Gbsa forces and energy \n

   @param atomCoordinates   Gromacs atom coordinates ('x' in md.c)
   @param partialCharges    Gromacs charges ('mdatoms->chargeA' in md.c)
   @param energy            output energy
   @param forces            output forces (format currently forces[3][numberOfAtoms])
                            if 'forces' reference is NULL upon input, then \n
                            memory is allocated for array; responsibilty of callee  \n
                            to free \n
                              delete forces[0]; \n
                              delete forces; \n
                            Note memory is allocated as single block of size \n
                            3*numberOfAtoms*sizeof( float ), so only free \n
                            'forces[0]' \n
   @param log               log reference (stdlog in md.c)

   Function calls a static method in CpuGbsa class to calculate forces/energy

   Logging can be toggled by setting cpuGbsaLog to input 'log' or NULL

   @return force array; energy is returned in *energy

   --------------------------------------------------------------------------------------- */

extern "C" float** 
cpuCalculateGbsaForces( const rvec *atomCoordinates, const float* partialCharges, 
                        float** forces, FILE* log ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName             = "\ncpuCalculateGbsaForces: ";

   // ---------------------------------------------------------------------------------------

   // FILE* cpuGbsaLog               = NULL;
   FILE* cpuGbsaLog               = log;

   // (void) fprintf( log, "%s start", methodName );
   // (void) fflush( log );

   forces  = CpuGbsa::computeGbsaForces( atomCoordinates, partialCharges, forces, cpuGbsaLog );
//   *energy = CpuGbsa::getCpuGbsa()->getEnergy(); 

   // (void) fprintf( log, "%s done E=%.4f", methodName, *energy );
   // (void) fflush( log );

   return forces;

}

/**---------------------------------------------------------------------------------------

   Retrieve the calculated Gbsa energy from the static class member

   @return the calculated Gbsa energy from the static class member

   --------------------------------------------------------------------------------------- */

extern "C" float cpuGetEnergy( void ){
   return CpuGbsa::getCpuGbsa()->getEnergy();
}

/**---------------------------------------------------------------------------------------

   Delete the Gbsa associated object(s)

   @return 0 if static CpuGbsa object was set; else return -1

   --------------------------------------------------------------------------------------- */

extern "C" int cpuDeleteGbsaParameters( void ){
   return CpuGbsa::deleteCpuGbsa();
}

