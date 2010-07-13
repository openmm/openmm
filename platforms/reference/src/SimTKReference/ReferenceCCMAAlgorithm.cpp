
/* Portions copyright (c) 2006-2009 Stanford University and Simbios.
 * Contributors: Peter Eastman, Pande Group
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
#include "../SimTKUtilities/SimTKOpenMMUtilities.h"
#include "../SimTKUtilities/SimTKOpenMMLog.h"
#include "ReferenceCCMAAlgorithm.h"
#include "ReferenceDynamics.h"
#include "quern.h"
#include "openmm/Vec3.h"
#include <map>

using std::map;
using std::pair;
using std::vector;
using std::set;
using OpenMM::Vec3;

/**---------------------------------------------------------------------------------------

   ReferenceCCMAAlgorithm constructor

         @param numberOfAtoms    number of atoms
         @param numberOfConstraints      number of constraints
         @param atomIndices              atom indices for contraints
         @param distance                 distances for constraints
         @param masses                   atom masses
         @param angles                   angle force field terms
         @param tolerance                constraint tolerance

   --------------------------------------------------------------------------------------- */

ReferenceCCMAAlgorithm::ReferenceCCMAAlgorithm( int numberOfAtoms,
                                                  int numberOfConstraints,
                                                  int** atomIndices,
                                                  RealOpenMM* distance,
                                                  RealOpenMM* masses,
                                                  vector<AngleInfo>& angles,
                                                  RealOpenMM tolerance){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceCCMAAlgorithm::ReferenceCCMAAlgorithm";

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
   if (numberOfConstraints > 0)
   {
       // Compute the constraint coupling matrix

       vector<vector<int> > atomAngles(numberOfAtoms);
       for (int i = 0; i < (int) angles.size(); i++)
           atomAngles[angles[i].atom2].push_back(i);
       vector<vector<pair<int, double> > > matrix(numberOfConstraints);
       for (int j = 0; j < numberOfConstraints; j++) {
           for (int k = 0; k < numberOfConstraints; k++) {
               if (j == k) {
                   matrix[j].push_back(pair<int, double>(j, 1.0));
                   continue;
               }
               double scale;
               int atomj0 = _atomIndices[j][0];
               int atomj1 = _atomIndices[j][1];
               int atomk0 = _atomIndices[k][0];
               int atomk1 = _atomIndices[k][1];
               RealOpenMM invMass0 = one/masses[atomj0];
               RealOpenMM invMass1 = one/masses[atomj1];
               int atoma, atomb, atomc;
               if (atomj0 == atomk0) {
                   atoma = atomj1;
                   atomb = atomj0;
                   atomc = atomk1;
                   scale = invMass0/(invMass0+invMass1);
               }
               else if (atomj1 == atomk1) {
                   atoma = atomj0;
                   atomb = atomj1;
                   atomc = atomk0;
                   scale = invMass1/(invMass0+invMass1);
               }
               else if (atomj0 == atomk1) {
                   atoma = atomj1;
                   atomb = atomj0;
                   atomc = atomk0;
                   scale = invMass0/(invMass0+invMass1);
               }
               else if (atomj1 == atomk0) {
                   atoma = atomj0;
                   atomb = atomj1;
                   atomc = atomk1;
                   scale = invMass1/(invMass0+invMass1);
               }
               else
                   continue; // These constraints are not connected.

               // Look for a third constraint forming a triangle with these two.

               bool foundConstraint = false;
               for (int other = 0; other < numberOfConstraints; other++) {
                   if ((_atomIndices[other][0] == atoma && _atomIndices[other][1] == atomc) || (_atomIndices[other][0] == atomc && _atomIndices[other][1] == atoma)) {
                       double d1 = _distance[j];
                       double d2 = _distance[k];
                       double d3 = _distance[other];
                       matrix[j].push_back(pair<int, double>(k, scale*(d1*d1+d2*d2-d3*d3)/(2.0*d1*d2)));
                       foundConstraint = true;
                       break;
                   }
               }
               if (!foundConstraint) {
                   // We didn't find one, so look for an angle force field term.

                   const vector<int>& angleCandidates = atomAngles[atomb];
                   for (vector<int>::const_iterator iter = angleCandidates.begin(); iter != angleCandidates.end(); iter++) {
                       const AngleInfo& angle = angles[*iter];
                       if ((angle.atom1 == atoma && angle.atom3 == atomc) || (angle.atom3 == atoma && angle.atom1 == atomc)) {
                           matrix[j].push_back(pair<int, double>(k, scale*cos(angle.angle)));
                           break;
                       }
                   }
               }
           }
       }

       // Invert it using QR.

       vector<int> matrixRowStart;
       vector<int> matrixColIndex;
       vector<double> matrixValue;
       for (int i = 0; i < numberOfConstraints; i++) {
           matrixRowStart.push_back(matrixValue.size());
           for (int j = 0; j < (int) matrix[i].size(); j++) {
               pair<int, double> element = matrix[i][j];
               matrixColIndex.push_back(element.first);
               matrixValue.push_back(element.second);
           }
       }
       matrixRowStart.push_back(matrixValue.size());
       int *qRowStart, *qColIndex, *rRowStart, *rColIndex;
       double *qValue, *rValue;
       QUERN_compute_qr(numberOfConstraints, numberOfConstraints, &matrixRowStart[0], &matrixColIndex[0], &matrixValue[0], NULL,
               &qRowStart, &qColIndex, &qValue, &rRowStart, &rColIndex, &rValue);
       vector<double> rhs(numberOfConstraints);
       _matrix.resize(numberOfConstraints);
       for (int i = 0; i < numberOfConstraints; i++) {
           // Extract column i of the inverse matrix.

           for (int j = 0; j < numberOfConstraints; j++)
               rhs[j] = (i == j ? 1.0 : 0.0);
           QUERN_multiply_with_q_transpose(numberOfConstraints, qRowStart, qColIndex, qValue, &rhs[0]);
           QUERN_solve_with_r(numberOfConstraints, rRowStart, rColIndex, rValue, &rhs[0], &rhs[0]);
           for (int j = 0; j < numberOfConstraints; j++) {
               double value = rhs[j]*_distance[i]/_distance[j];
               if (FABS((RealOpenMM)value) > 0.02)
                   _matrix[j].push_back(pair<int, RealOpenMM>(i, (RealOpenMM) value));
           }
       }
       QUERN_free_result(qRowStart, qColIndex, qValue);
       QUERN_free_result(rRowStart, rColIndex, rValue);
   }
}

/**---------------------------------------------------------------------------------------

   ReferenceCCMAAlgorithm destructor

   --------------------------------------------------------------------------------------- */

ReferenceCCMAAlgorithm::~ReferenceCCMAAlgorithm( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceCCMAAlgorithm::~ReferenceCCMAAlgorithm";

   // ---------------------------------------------------------------------------------------

    if (_numberOfConstraints > 0) {
       SimTKOpenMMUtilities::freeTwoDRealOpenMMArray( _r_ij,  "r_ij" );

       SimTKOpenMMUtilities::freeOneDRealOpenMMArray( _d_ij2, "d_ij2" );
       SimTKOpenMMUtilities::freeOneDRealOpenMMArray( _distanceTolerance, "distanceTolerance" );
       SimTKOpenMMUtilities::freeOneDRealOpenMMArray( _reducedMasses, "reducedMasses" );
    }
//    for (unsigned int i = 0; i < _matrices.size(); i++)
//        SimTKOpenMMUtilities::freeTwoDRealOpenMMArray(_matrices[i], "");
}

/**---------------------------------------------------------------------------------------

   Get number of constraints

   @return number of constraints

   --------------------------------------------------------------------------------------- */

int ReferenceCCMAAlgorithm::getNumberOfConstraints( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceCCMAAlgorithm::getNumberOfConstraints";

   // ---------------------------------------------------------------------------------------

   return _numberOfConstraints;
}

/**---------------------------------------------------------------------------------------

   Get maximum number of iterations

   @return maximum number of iterations

   --------------------------------------------------------------------------------------- */

int ReferenceCCMAAlgorithm::getMaximumNumberOfIterations( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceCCMAAlgorithm::getMaximumNumberOfIterations";

   // ---------------------------------------------------------------------------------------

   return _maximumNumberOfIterations;
}

/**---------------------------------------------------------------------------------------

   Set maximum number of iterations

   @param maximumNumberOfIterations   new maximum number of iterations

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceCCMAAlgorithm::setMaximumNumberOfIterations( int maximumNumberOfIterations ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceCCMAAlgorithm::setMaximumNumberOfIterations";

   // ---------------------------------------------------------------------------------------

   _maximumNumberOfIterations = maximumNumberOfIterations;

   return ReferenceDynamics::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Get tolerance

   @return tolerance

   --------------------------------------------------------------------------------------- */

RealOpenMM ReferenceCCMAAlgorithm::getTolerance( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceCCMAAlgorithm::getTolerance";

   // ---------------------------------------------------------------------------------------

   return _tolerance;
}

/**---------------------------------------------------------------------------------------

   Set tolerance

   @param tolerance new tolerance

   --------------------------------------------------------------------------------------- */

void ReferenceCCMAAlgorithm::setTolerance( RealOpenMM tolerance ){


   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceCCMAAlgorithm::setTolerance";

   // ---------------------------------------------------------------------------------------

   _tolerance = tolerance;;
}

/**---------------------------------------------------------------------------------------

   Apply CCMA algorithm

   @param numberOfAtoms    number of atoms
   @param atomCoordinates  atom coordinates
   @param atomCoordinatesP atom coordinates prime
   @param inverseMasses    1/mass

   @return ReferenceDynamics::DefaultReturn if converge; else
    return ReferenceDynamics::ErrorReturn

   --------------------------------------------------------------------------------------- */

int ReferenceCCMAAlgorithm::apply( int numberOfAtoms, RealOpenMM** atomCoordinates,
                                         RealOpenMM** atomCoordinatesP,
                                         RealOpenMM* inverseMasses ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nReferenceCCMAAlgorithm::apply";

   static const RealOpenMM zero        =  0.0;
   static const RealOpenMM one         =  1.0;
   static const RealOpenMM two         =  2.0;
   static const RealOpenMM half        =  0.5;

   static const RealOpenMM epsilon6    = (RealOpenMM) 1.0e-06;

   static int debug                    = 0;

   // ---------------------------------------------------------------------------------------

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
         r_ij[ii][jj] = atomCoordinates[atomI][jj] - atomCoordinates[atomJ][jj];
      }
      d_ij2[ii]              = DOT3( r_ij[ii], r_ij[ii] );
      distanceTolerance[ii]  = d_ij2[ii]*tolerance;
      if( distanceTolerance[ii] > zero ){
         distanceTolerance[ii] = one/distanceTolerance[ii];
      }
   }
   RealOpenMM lowerTol = one-two*getTolerance()+getTolerance()*getTolerance();
   RealOpenMM upperTol = one+two*getTolerance()+getTolerance()*getTolerance();

   // main loop

   int iterations           = 0;
   int numberConverged      = 0;
   vector<RealOpenMM> constraintDelta(_numberOfConstraints);
   vector<RealOpenMM> tempDelta(_numberOfConstraints);
   while (iterations < getMaximumNumberOfIterations()) {
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
         constraintDelta[ii] = zero;
         RealOpenMM rrpr  = DOT3(  rp_ij, r_ij[ii] );
         if( rrpr <  d_ij2[ii]*epsilon6 ){
             std::stringstream message;
             message << iterations <<" "<<atomI<<" "<<atomJ<< " Error: sign of rrpr < 0?\n";
             SimTKOpenMMLog::printMessage( message );
         } else {
             constraintDelta[ii] = reducedMasses[ii]*diff/rrpr;
         }
         if (rp2 >= lowerTol*dist2 && rp2 <= upperTol*dist2) {
            numberConverged++;
         }
      }
      if( numberConverged == _numberOfConstraints )
         break;
      iterations++;

      if (_matrix.size() > 0) {
          for (int i = 0; i < _numberOfConstraints; i++) {
              RealOpenMM sum = 0.0;
              for (int j = 0; j < (int) _matrix[i].size(); j++) {
                  pair<int, RealOpenMM> element = _matrix[i][j];
                  sum += element.second*constraintDelta[element.first];
              }
              tempDelta[i] = sum;
          }
          constraintDelta = tempDelta;
      }
      for( int ii = 0; ii < _numberOfConstraints; ii++ ){

         int atomI   = _atomIndices[ii][0];
         int atomJ   = _atomIndices[ii][1];
         for( int jj = 0; jj < 3; jj++ ){
            RealOpenMM dr                = constraintDelta[ii]*r_ij[ii][jj];
            atomCoordinatesP[atomI][jj] += inverseMasses[atomI]*dr;
            atomCoordinatesP[atomJ][jj] -= inverseMasses[atomJ]*dr;
         }
      }
   }

   // diagnostics

   if( debug || numberConverged < _numberOfConstraints ){
      std::stringstream message;
      message << methodName;
      message << " iterations=" << iterations << " no. converged=" << numberConverged << " out of " << _numberOfConstraints;
      if( numberConverged == _numberOfConstraints ){
         message << " SUCCESS";
      } else {
         message << " FAILED";
      }
      message << "\n";
      int errors = reportCCMA( numberOfAtoms, atomCoordinatesP, message );
      if( !errors ){
         message << "*** no errors recorded in explicit check ***";
      }
      message << "\n";
      SimTKOpenMMLog::printMessage( message );
   }

   return (numberConverged == _numberOfConstraints ? ReferenceDynamics::DefaultReturn : ReferenceDynamics::ErrorReturn);

}

/**---------------------------------------------------------------------------------------

   Report any violated constriants

   @param numberOfAtoms    number of atoms
   @param atomCoordinates  atom coordinates
   @param message          report

   @return number of violated constraints

   --------------------------------------------------------------------------------------- */

int ReferenceCCMAAlgorithm::reportCCMA( int numberOfAtoms, RealOpenMM** atomCoordinates,
                                          std::stringstream& message ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nReferenceCCMAAlgorithm::reportCCMA";

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

