
/* Portions copyright (c) 2006-2009 Stanford University and Simbios.
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

#include "SimTKOpenMMCommon.h"
#include "SimTKOpenMMUtilities.h"
#include "SimTKOpenMMLog.h"
#include "ReferenceShakeAlgorithm.h"
#include "ReferenceDynamics.h"
#include "openmm/OpenMMException.h"

using std::vector;
using OpenMM::RealVec;

/**---------------------------------------------------------------------------------------

   ReferenceShakeAlgorithm constructor

   @param numberOfConstraints      number of constraints
   @param atomIndices              block of atom indices
   @param shakeParameters          Shake parameters
   @param tolerance                constraint tolerance

   --------------------------------------------------------------------------------------- */

ReferenceShakeAlgorithm::ReferenceShakeAlgorithm( int numberOfConstraints,
                                                  int** atomIndices,
                                                  RealOpenMM* distance,
                                                  RealOpenMM tolerance){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceShakeAlgorithm::ReferenceShakeAlgorithm";

   static const RealOpenMM zero        =  0.0;
   static const RealOpenMM one         =  1.0;
   static const RealOpenMM two         =  2.0;
   static const RealOpenMM three       =  3.0;
   static const RealOpenMM oneM        = -1.0;

   static const int threeI             =  3;

   // ---------------------------------------------------------------------------------------

   _numberOfConstraints        = numberOfConstraints;
   _atomIndices                = atomIndices;
   _distance                   = distance;

   _maximumNumberOfIterations  = 150;
   _tolerance                  = tolerance;
   _hasInitializedMasses       = false;

   // work arrays

   if (_numberOfConstraints > 0) {
       _r_ij                       = SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray( numberOfConstraints, threeI, NULL, 
                                                                                        1, zero, "r_ij" );
    
       _d_ij2                      = SimTKOpenMMUtilities::allocateOneDRealOpenMMArray( numberOfConstraints, NULL, 1, zero, "dij_2" );
       _distanceTolerance          = SimTKOpenMMUtilities::allocateOneDRealOpenMMArray( numberOfConstraints, NULL, 1, zero, "distanceTolerance" );
       _reducedMasses              = SimTKOpenMMUtilities::allocateOneDRealOpenMMArray( numberOfConstraints, NULL, 1, zero, "reducedMasses" );
   }
}

/**---------------------------------------------------------------------------------------

   ReferenceShakeAlgorithm destructor

   --------------------------------------------------------------------------------------- */

ReferenceShakeAlgorithm::~ReferenceShakeAlgorithm( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceShakeAlgorithm::~ReferenceShakeAlgorithm";

   // ---------------------------------------------------------------------------------------

    if (_numberOfConstraints > 0) {
       SimTKOpenMMUtilities::freeTwoDRealOpenMMArray( _r_ij,  "r_ij" );
    
       SimTKOpenMMUtilities::freeOneDRealOpenMMArray( _d_ij2, "d_ij2" );
       SimTKOpenMMUtilities::freeOneDRealOpenMMArray( _distanceTolerance, "distanceTolerance" );
       SimTKOpenMMUtilities::freeOneDRealOpenMMArray( _reducedMasses, "reducedMasses" );
    }
}

/**---------------------------------------------------------------------------------------

   Get number of constraints

   @return number of constraints

   --------------------------------------------------------------------------------------- */
   
int ReferenceShakeAlgorithm::getNumberOfConstraints( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceShakeAlgorithm::getNumberOfConstraints";

   // ---------------------------------------------------------------------------------------

   return _numberOfConstraints;
}
   
/**---------------------------------------------------------------------------------------

   Get maximum number of iterations

   @return maximum number of iterations

   --------------------------------------------------------------------------------------- */
   
int ReferenceShakeAlgorithm::getMaximumNumberOfIterations( void ) const {
   
   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceShakeAlgorithm::getMaximumNumberOfIterations";

   // ---------------------------------------------------------------------------------------

   return _maximumNumberOfIterations;
}

/**---------------------------------------------------------------------------------------

   Set maximum number of iterations

   @param maximumNumberOfIterations   new maximum number of iterations

   --------------------------------------------------------------------------------------- */
   
void ReferenceShakeAlgorithm::setMaximumNumberOfIterations( int maximumNumberOfIterations ){
   
   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceShakeAlgorithm::setMaximumNumberOfIterations";

   // ---------------------------------------------------------------------------------------

   _maximumNumberOfIterations = maximumNumberOfIterations;
}

/**---------------------------------------------------------------------------------------

   Get tolerance

   @return tolerance

   --------------------------------------------------------------------------------------- */

RealOpenMM ReferenceShakeAlgorithm::getTolerance( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceShakeAlgorithm::getTolerance";

   // ---------------------------------------------------------------------------------------

   return _tolerance;
}

/**---------------------------------------------------------------------------------------

   Set tolerance

   @param tolerance new tolerance

   --------------------------------------------------------------------------------------- */

void ReferenceShakeAlgorithm::setTolerance( RealOpenMM tolerance ){


   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceShakeAlgorithm::setTolerance";

   // ---------------------------------------------------------------------------------------

   _tolerance = tolerance;;
}

/**---------------------------------------------------------------------------------------

   Apply Shake algorithm

   @param numberOfAtoms    number of atoms
   @param atomCoordinates  atom coordinates
   @param atomCoordinatesP atom coordinates prime
   @param inverseMasses    1/mass

   @return SimTKOpenMMCommon::DefaultReturn if converge; else
    return SimTKOpenMMCommon::ErrorReturn

   --------------------------------------------------------------------------------------- */

int ReferenceShakeAlgorithm::apply( int numberOfAtoms, vector<RealVec>& atomCoordinates,
                                         vector<RealVec>& atomCoordinatesP,
                                         vector<RealOpenMM>& inverseMasses ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nReferenceShakeAlgorithm::apply";

   static const RealOpenMM zero        =  0.0;
   static const RealOpenMM one         =  1.0;
   static const RealOpenMM two         =  2.0;
   static const RealOpenMM three       =  3.0;
   static const RealOpenMM oneM        = -1.0;
   static const RealOpenMM half        =  0.5;

   static const RealOpenMM epsilon6    = (RealOpenMM) 1.0e-06;

   // ---------------------------------------------------------------------------------------

   int numberOfConstraints = getNumberOfConstraints();

   // temp arrays

   RealOpenMM** r_ij                   = _r_ij;
   RealOpenMM* d_ij2                   = _d_ij2;
   RealOpenMM* distanceTolerance       = _distanceTolerance;
   RealOpenMM* reducedMasses           = _reducedMasses;

   // calculate reduced masses on 1st pass

   if( !_hasInitializedMasses ){
      _hasInitializedMasses = true;
      for( int ii = 0; ii < _numberOfConstraints; ii++ ){
         int atomI          = _atomIndices[ii][0];
         int atomJ          = _atomIndices[ii][1];
         reducedMasses[ii]  = half/( inverseMasses[atomI] + inverseMasses[atomJ] );
      }
   }

   // setup: r_ij for each (i,j) constraint

   RealOpenMM tolerance     = getTolerance();
              tolerance    *= two;
   for( int ii = 0; ii < _numberOfConstraints; ii++ ){

      int atomI   = _atomIndices[ii][0];
      int atomJ   = _atomIndices[ii][1];
      for( int jj = 0; jj < 3; jj++ ){
         r_ij[ii][jj]   = atomCoordinates[atomI][jj]  - atomCoordinates[atomJ][jj];
      }
      d_ij2[ii]              = DOT3( r_ij[ii], r_ij[ii] );
      distanceTolerance[ii]  = d_ij2[ii]*tolerance;
      if( distanceTolerance[ii] > zero ){
         distanceTolerance[ii] = one/distanceTolerance[ii];
      }
   }

   // main loop

   int done                 = 0;
   int iterations           = 0;
   int numberConverged      = 0;
   while( !done && iterations++ < getMaximumNumberOfIterations() ){
      numberConverged  = 0; 
      for( int ii = 0; ii < _numberOfConstraints; ii++ ){

         int atomI   = _atomIndices[ii][0];
         int atomJ   = _atomIndices[ii][1];

         RealOpenMM rp_ij[3];
         for( int jj = 0; jj < 3; jj++ ){
            rp_ij[jj] = atomCoordinatesP[atomI][jj] - atomCoordinatesP[atomJ][jj]; 
         }

         RealOpenMM rp2  = DOT3( rp_ij, rp_ij );
         RealOpenMM dist2= _distance[ii]*_distance[ii];
         RealOpenMM diff = dist2 - rp2;
         int iconv       = (int) ( FABS( diff )*distanceTolerance[ii] );
         if( iconv ){
            RealOpenMM rrpr  = DOT3(  rp_ij, r_ij[ii] );
            RealOpenMM acor;
            if( rrpr <  d_ij2[ii]*epsilon6 ){
               std::stringstream message;
               message << iterations << " Error: sign of rrpr < 0?\n";
               SimTKOpenMMLog::printMessage( message );
            } else {
               acor  = reducedMasses[ii]*diff/rrpr;
               for( int jj = 0; jj < 3; jj++ ){
                  RealOpenMM dr                = acor*r_ij[ii][jj];
                  atomCoordinatesP[atomI][jj] += inverseMasses[atomI]*dr;
                  atomCoordinatesP[atomJ][jj] -= inverseMasses[atomJ]*dr;
               }
            }
         } else {
            numberConverged++;
         }
      }
      if( numberConverged == _numberOfConstraints ){
         done = true;
      }
   }

   return (done ? SimTKOpenMMCommon::DefaultReturn : SimTKOpenMMCommon::ErrorReturn);

}

/**---------------------------------------------------------------------------------------

   Apply constraint algorithm to velocities.

   @param numberOfAtoms    number of atoms
   @param atomCoordinates  atom coordinates
   @param velocities       atom velocities
   @param inverseMasses    1/mass

   @return SimTKOpenMMCommon::DefaultReturn if converge; else
    return SimTKOpenMMCommon::ErrorReturn

   --------------------------------------------------------------------------------------- */

int ReferenceShakeAlgorithm::applyToVelocities(int numberOfAtoms, std::vector<OpenMM::RealVec>& atomCoordinates,
               std::vector<OpenMM::RealVec>& velocities, std::vector<RealOpenMM>& inverseMasses) {
    throw OpenMM::OpenMMException("applyToVelocities is not implemented");
}

/**---------------------------------------------------------------------------------------

   Report any violated constriants

   @param numberOfAtoms    number of atoms
   @param atomCoordinates  atom coordinates
   @param message          report

   @return number of violated constraints

   --------------------------------------------------------------------------------------- */

int ReferenceShakeAlgorithm::reportShake( int numberOfAtoms, vector<RealVec>& atomCoordinates,
                                          std::stringstream& message ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nReferenceShakeAlgorithm::reportShake";

   static const RealOpenMM zero        =  0.0;
   static const RealOpenMM one         =  1.0;
   static const RealOpenMM two         =  2.0;
   static const RealOpenMM three       =  3.0;
   static const RealOpenMM oneM        = -1.0;
   static const RealOpenMM half        =  0.5;

   // ---------------------------------------------------------------------------------------

   int numberOfConstraints = getNumberOfConstraints();

   // loop over constraints calculating distance and comparing to
   // expected distance -- report any contraints that are violated

   int numberConverged  = 0; 
   RealOpenMM tolerance = getTolerance();
   for( int ii = 0; ii < _numberOfConstraints; ii++ ){

      int atomI   = _atomIndices[ii][0];
      int atomJ   = _atomIndices[ii][1];

      RealOpenMM rp2  = zero;
      for( int jj = 0; jj < 3; jj++ ){
         rp2 += (atomCoordinates[atomI][jj] - atomCoordinates[atomJ][jj])*(atomCoordinates[atomI][jj] - atomCoordinates[atomJ][jj]); 
      }
      RealOpenMM diff = FABS( rp2 - (_distance[ii]*_distance[ii]) );
      if( diff > tolerance ){
         message << ii << " constraint violated: " << atomI << " " << atomJ << "] d=" << SQRT( rp2 ) << " " << rp2 << " d0=" << _distance[ii];
         message << " diff=" << diff;
         message << " [" << atomCoordinates[atomI][0] << " " << atomCoordinates[atomI][1] << " " << atomCoordinates[atomI][2] << "] ";
         message << " [" << atomCoordinates[atomJ][0] << " " << atomCoordinates[atomJ][1] << " " << atomCoordinates[atomJ][2] << "] ";
         message << "\n";
      } else {
         numberConverged++;
      }
   }
          
   return (numberOfConstraints-numberConverged);

}
