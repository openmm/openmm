
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

#include <string.h>
#include <sstream>

#include "../SimTKUtilities/SimTKOpenMMCommon.h"
#include "../SimTKUtilities/SimTKOpenMMLog.h"
#include "../SimTKUtilities/SimTKOpenMMUtilities.h"
#include "ReferenceForce.h"

/**---------------------------------------------------------------------------------------

   ReferenceForce constructor

   --------------------------------------------------------------------------------------- */

ReferenceForce::ReferenceForce( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceForce::ReferenceForce";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   ReferenceForce destructor

   --------------------------------------------------------------------------------------- */

ReferenceForce::~ReferenceForce( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceForce::~ReferenceForce";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   Given two coordinates on a periodic lattice, return the difference between them.

   --------------------------------------------------------------------------------------- */

RealOpenMM ReferenceForce::periodicDifference(RealOpenMM val1, RealOpenMM val2, RealOpenMM period) {
    RealOpenMM diff = val1-val2;
    RealOpenMM base = floor(diff/period+0.5)*period;
    return diff-base;

}


/**---------------------------------------------------------------------------------------

   Get deltaR and distance and distance**2 between atomI and atomJ (static method)
   deltaR: j - i

   @param atomCoordinatesI    atom i coordinates
   @param atomCoordinatesI    atom j coordinates
   @param deltaR              deltaX, deltaY, deltaZ, R2, R upon return

   @return ReferenceForce::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceForce::getDeltaR( const RealOpenMM* atomCoordinatesI, const RealOpenMM* atomCoordinatesJ,
                               RealOpenMM* deltaR ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nReferenceForce::getDeltaR";

   // ---------------------------------------------------------------------------------------

   deltaR[XIndex]    = atomCoordinatesJ[0] - atomCoordinatesI[0];
   deltaR[YIndex]    = atomCoordinatesJ[1] - atomCoordinatesI[1];
   deltaR[ZIndex]    = atomCoordinatesJ[2] - atomCoordinatesI[2];

   deltaR[R2Index]   = DOT3( deltaR, deltaR );
   deltaR[RIndex]    = (RealOpenMM) SQRT( deltaR[R2Index] );

   return ReferenceForce::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Get deltaR and distance and distance**2 between atomI and atomJ, assuming periodic
   boundary conditions (static method); deltaR: j - i

   @param atomCoordinatesI    atom i coordinates
   @param atomCoordinatesI    atom j coordinates
   @param boxSize             X, Y, and Z sizes of the periodic box
   @param deltaR              deltaX, deltaY, deltaZ, R2, R upon return

   @return ReferenceForce::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceForce::getDeltaRPeriodic( const RealOpenMM* atomCoordinatesI, const RealOpenMM* atomCoordinatesJ,
                               const RealOpenMM* boxSize, RealOpenMM* deltaR ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nReferenceForce::getDeltaR";

   // ---------------------------------------------------------------------------------------

   deltaR[XIndex]    = periodicDifference(atomCoordinatesJ[0], atomCoordinatesI[0], boxSize[0]);
   deltaR[YIndex]    = periodicDifference(atomCoordinatesJ[1], atomCoordinatesI[1], boxSize[1]);
   deltaR[ZIndex]    = periodicDifference(atomCoordinatesJ[2], atomCoordinatesI[2], boxSize[2]);

   deltaR[R2Index]   = DOT3( deltaR, deltaR );
   deltaR[RIndex]    = (RealOpenMM) SQRT( deltaR[R2Index] );

   return ReferenceForce::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Get deltaR between atomI and atomJ (static method); deltaR: j - i

   @param atomCoordinatesI    atom i coordinates
   @param atomCoordinatesI    atom j coordinates
   @param deltaR              deltaX, deltaY, deltaZ upon return

   @return ReferenceForce::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceForce::getDeltaROnly( const RealOpenMM* atomCoordinatesI,
                                   const RealOpenMM* atomCoordinatesJ,
                                   RealOpenMM* deltaR ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nReferenceForce::getDeltaR";

   // ---------------------------------------------------------------------------------------

   deltaR[XIndex]    = atomCoordinatesJ[0] - atomCoordinatesI[0];
   deltaR[YIndex]    = atomCoordinatesJ[1] - atomCoordinatesI[1];
   deltaR[ZIndex]    = atomCoordinatesJ[2] - atomCoordinatesI[2];

   return ReferenceForce::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Write coordinates, energies and forces (Simbios)

   @param numberOfAtoms       number of atoms
   @param atomsPerBond        atoms/bond
   @param atomCoordinates     atomic coordinates
   @param forces              forces
   @param energies            energies (optional)
   @param resultsFileName     output file name
   @param lengthConversion    length conversion (optional defaults to 1)
   @param energyConversion    energy conversion (optional defaults to 1)
   @param forceConversion     force  conversion (optional defaults to 1)


   @return ReferenceForce::DefaultReturn unless
           file cannot be opened
           in which case return ReferenceForce::ErrorReturn

   --------------------------------------------------------------------------------------- */

int ReferenceForce::writeForces( int numberOfAtoms, int atomsPerBond,
                                 RealOpenMM** atomCoordinates, RealOpenMM** forces,
                                 RealOpenMM* energies, const std::string& resultsFileName,
                                 RealOpenMM lengthConversion  = 1.0,
                                 RealOpenMM forceConversion   = 1.0,
                                 RealOpenMM energyConversion  = 1.0 ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName  = "\nReferenceForce::writeBornEnergyForces";
   static const RealOpenMM zero   = 0.0;

   // ---------------------------------------------------------------------------------------

   // open file -- return if unsuccessful

   FILE* resultsFile = NULL;
#ifdef WIN32
   fopen_s( &resultsFile, resultsFileName.c_str(), "w" );
#else
   resultsFile = fopen( resultsFileName.c_str(), "w" );
#endif

   // ---------------------------------------------------------------------------------------

   // diagnostics

   std::stringstream message;
   message << methodName;
   if( resultsFile != NULL ){
      std::stringstream message;
      message << methodName;
      message << " Opened file=<" << resultsFileName << ">.";
      SimTKOpenMMLog::printMessage( message );
   } else {
      std::stringstream message;
      message << methodName;
      message << "  could not open file=<" << resultsFileName << "> -- abort output.";
      SimTKOpenMMLog::printMessage( message );
      return ReferenceForce::ErrorReturn;
   }

   // total energy and normalize by number of atoms/bond

   RealOpenMM totalEnergy = zero;
   if( energies ){

      for( int ii = 1; ii < numberOfAtoms; ii++ ){
         totalEnergy += energies[ii];
      }
      if( atomsPerBond > 0 ){
         totalEnergy /= (RealOpenMM) atomsPerBond;
      }
   }
 
   // header

   (void) fprintf( resultsFile, "# %d atoms E=%.7e  format: coords(3) forces(3) energies\n",
                   numberOfAtoms, totalEnergy );

   // output

   if( forces != NULL && atomCoordinates != NULL ){

      for( int ii = 0; ii < numberOfAtoms; ii++ ){
            (void) fprintf( resultsFile, "%5d %15.7e %15.7e %15.7e  %15.7e %15.7e %15.7e", ii,
                            lengthConversion*atomCoordinates[ii][0],
                            lengthConversion*atomCoordinates[ii][1], 
                            lengthConversion*atomCoordinates[ii][2],
                            forceConversion*forces[ii][0],
                            forceConversion*forces[ii][1],
                            forceConversion*forces[ii][2] 
                          );
         if( energies != NULL ){
            (void) fprintf( resultsFile, "  %15.7e", energyConversion*energies[ii] );
         }
         (void) fprintf( resultsFile, "\n" );
      }
   }
   (void) fclose( resultsFile );

   return ReferenceForce::DefaultReturn;

}
