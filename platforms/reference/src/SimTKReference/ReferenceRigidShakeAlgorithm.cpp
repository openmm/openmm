
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
#include "ReferenceRigidShakeAlgorithm.h"
#include "ReferenceDynamics.h"
#include "jama_svd.h"
#include "openmm/Vec3.h"
#include <map>

using std::map;
using std::vector;
using std::set;
using OpenMM::Vec3;
using TNT::Array2D;
using JAMA::SVD;

/**---------------------------------------------------------------------------------------

   ReferenceQShakeAlgorithm constructor

   @param numberOfAtoms            number of atoms
   @param numberOfConstraints      number of constraints
   @param atomIndices              block of atom indices
   @param shakeParameters          Shake parameters
   @param tolerance                constraint tolerance

   --------------------------------------------------------------------------------------- */

ReferenceRigidShakeAlgorithm::ReferenceRigidShakeAlgorithm( int numberOfAtoms,
                                                  int numberOfConstraints,
                                                  int** atomIndices,
                                                  RealOpenMM* distance,
                                                  RealOpenMM tolerance){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceQShakeAlgorithm::ReferenceQShakeAlgorithm";

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

   // Identify rigid clusters.

   vector<map<int, int> > atomConstraints(numberOfAtoms);
   for (int i = 0; i < numberOfConstraints; i++) {
       atomConstraints[atomIndices[i][0]][atomIndices[i][1]] = i;
       atomConstraints[atomIndices[i][1]][atomIndices[i][0]] = i;
   }
   for (int i = 0; i < numberOfAtoms; i++) {
       if (atomConstraints[i].size() < 2)
           continue;

       // Begin by looking for a triangle this atom is part of.

       set<int> atoms;
       atoms.insert(i);
       for (map<int, int>::const_iterator atom1 = atomConstraints[i].begin(); atom1 != atomConstraints[i].end() && atoms.size() == 1; ++atom1) {
           for (map<int, int>::const_iterator atom2 = atomConstraints[atom1->first].begin(); atom2 != atomConstraints[atom1->first].end(); ++atom2) {
               if (atomConstraints[i].count(atom2->first) != 0) {
                   atoms.insert(atom1->first);
                   atoms.insert(atom2->first);
                   break;
               }
           }
       }
       if (atoms.size() == 1)
           continue;

       // We have three atoms that are part of a cluster, so look for other atoms we can add.

       bool done = false;
       while (!done) {
           done = true;
           for (set<int>::const_iterator atom1 = atoms.begin(); atom1 != atoms.end(); ++atom1) {
               for (map<int, int>::const_iterator atom2 = atomConstraints[*atom1].begin(); atom2 != atomConstraints[*atom1].end(); ++atom2) {
                   if (atoms.find(atom2->first) != atoms.end())
                       continue; // This atom is already in the cluster.

                   // See if this atom is linked to three other atoms in the cluster.

                   int linkCount = 0;
                   for (map<int, int>::const_iterator atom3 = atomConstraints[atom2->first].begin(); atom3 != atomConstraints[atom2->first].end(); ++atom3)
                       if (atoms.find(atom3->first) != atoms.end())
                           linkCount++;
                   if (linkCount > 2) {
                       atoms.insert(atom2->first);
                       done = false;
                   }
               }
           }
       }

       // Record the cluster.

       vector<int> constraints;
       for (set<int>::const_iterator atom1 = atoms.begin(); atom1 != atoms.end(); ++atom1) {
           for (map<int, int>::const_iterator atom2 = atomConstraints[*atom1].begin(); atom2 != atomConstraints[*atom1].end(); ++atom2) {
               if (*atom1 < atom2->first && atoms.find(atom2->first) != atoms.end())
                   constraints.push_back(atom2->second);
           }
       }
       _rigidClusters.push_back(constraints);
       for (set<int>::const_iterator atom1 = atoms.begin(); atom1 != atoms.end(); ++atom1) {
           for (map<int, int>::const_iterator atom2 = atomConstraints[*atom1].begin(); atom2 != atomConstraints[*atom1].end(); ++atom2)
               atomConstraints[atom2->first].erase(*atom1);
           atomConstraints[*atom1].clear();
       }
       _matrices.push_back(SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray(constraints.size(), constraints.size(), NULL, false, 0.0, ""));
   }
}

/**---------------------------------------------------------------------------------------

   ReferenceQShakeAlgorithm destructor

   --------------------------------------------------------------------------------------- */

ReferenceRigidShakeAlgorithm::~ReferenceRigidShakeAlgorithm( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceQShakeAlgorithm::~ReferenceQShakeAlgorithm";

   // ---------------------------------------------------------------------------------------

    if (_numberOfConstraints > 0) {
       SimTKOpenMMUtilities::freeTwoDRealOpenMMArray( _r_ij,  "r_ij" );

       SimTKOpenMMUtilities::freeOneDRealOpenMMArray( _d_ij2, "d_ij2" );
       SimTKOpenMMUtilities::freeOneDRealOpenMMArray( _distanceTolerance, "distanceTolerance" );
       SimTKOpenMMUtilities::freeOneDRealOpenMMArray( _reducedMasses, "reducedMasses" );
    }
    for (int i = 0; i < _matrices.size(); i++)
        SimTKOpenMMUtilities::freeTwoDRealOpenMMArray(_matrices[i], "");
}

/**---------------------------------------------------------------------------------------

   Get number of constraints

   @return number of constraints

   --------------------------------------------------------------------------------------- */

int ReferenceRigidShakeAlgorithm::getNumberOfConstraints( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceQShakeAlgorithm::getNumberOfConstraints";

   // ---------------------------------------------------------------------------------------

   return _numberOfConstraints;
}

/**---------------------------------------------------------------------------------------

   Get maximum number of iterations

   @return maximum number of iterations

   --------------------------------------------------------------------------------------- */

int ReferenceRigidShakeAlgorithm::getMaximumNumberOfIterations( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceQShakeAlgorithm::getMaximumNumberOfIterations";

   // ---------------------------------------------------------------------------------------

   return _maximumNumberOfIterations;
}

/**---------------------------------------------------------------------------------------

   Set maximum number of iterations

   @param maximumNumberOfIterations   new maximum number of iterations

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceRigidShakeAlgorithm::setMaximumNumberOfIterations( int maximumNumberOfIterations ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceQShakeAlgorithm::setMaximumNumberOfIterations";

   // ---------------------------------------------------------------------------------------

   _maximumNumberOfIterations = maximumNumberOfIterations;

   return ReferenceDynamics::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Get tolerance

   @return tolerance

   --------------------------------------------------------------------------------------- */

RealOpenMM ReferenceRigidShakeAlgorithm::getTolerance( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceQShakeAlgorithm::getTolerance";

   // ---------------------------------------------------------------------------------------

   return _tolerance;
}

/**---------------------------------------------------------------------------------------

   Set tolerance

   @param tolerance new tolerance

   @return tolerance

   --------------------------------------------------------------------------------------- */

int ReferenceRigidShakeAlgorithm::setTolerance( RealOpenMM tolerance ){


   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceQShakeAlgorithm::setTolerance";

   // ---------------------------------------------------------------------------------------

   _tolerance = tolerance;;

   return ReferenceDynamics::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Apply Shake algorithm

   @param numberOfAtoms    number of atoms
   @param atomCoordinates  atom coordinates
   @param atomCoordinatesP atom coordinates prime
   @param inverseMasses    1/mass

   @return ReferenceDynamics::DefaultReturn if converge; else
    return ReferenceDynamics::ErrorReturn

   --------------------------------------------------------------------------------------- */

int ReferenceRigidShakeAlgorithm::apply( int numberOfAtoms, RealOpenMM** atomCoordinates,
                                         RealOpenMM** atomCoordinatesP,
                                         RealOpenMM* inverseMasses ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nReferenceQShakeAlgorithm::apply";

   static const RealOpenMM zero        =  0.0;
   static const RealOpenMM one         =  1.0;
   static const RealOpenMM two         =  2.0;
   static const RealOpenMM three       =  3.0;
   static const RealOpenMM oneM        = -1.0;
   static const RealOpenMM half        =  0.5;

   static const RealOpenMM epsilon6    = (RealOpenMM) 1.0e-06;

   static int debug                    = 0;

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
      vector<double> temp;
      for (int i = 0; i < _rigidClusters.size(); i++) {
          // Compute the constraint coupling matrix for this cluster.

          const vector<int>& cluster = _rigidClusters[i];
          int size = cluster.size();
          vector<Vec3> r(size);
          for (int j = 0; j < cluster.size(); j++) {
              int atom1 = _atomIndices[cluster[j]][0];
              int atom2 = _atomIndices[cluster[j]][1];
              r[j] = Vec3(atomCoordinates[atom1][0]-atomCoordinates[atom2][0],
                          atomCoordinates[atom1][1]-atomCoordinates[atom2][1],
                          atomCoordinates[atom1][2]-atomCoordinates[atom2][2]);
              double invLength = 1.0/sqrt(r[j].dot(r[j]));
              r[j][0] *= invLength;
              r[j][1] *= invLength;
              r[j][2] *= invLength;
          }
          Array2D<double> matrix(size, size);
          for (int j = 0; j < size; j++) {
              for (int k = 0; k < size; k++) {
                  double dot;
                  int atomj0 = _atomIndices[cluster[j]][0];
                  int atomj1 = _atomIndices[cluster[j]][1];
                  int atomk0 = _atomIndices[cluster[k]][0];
                  int atomk1 = _atomIndices[cluster[k]][1];
                  if (atomj0 == atomk0)
                      dot = r[j].dot(r[k])*inverseMasses[atomj0]/(inverseMasses[atomj0]+inverseMasses[atomj1]);
                  else if (atomj1 == atomk1)
                      dot = r[j].dot(r[k])*inverseMasses[atomj1]/(inverseMasses[atomj0]+inverseMasses[atomj1]);
                  else if (atomj0 == atomk1)
                      dot = -r[j].dot(r[k])*inverseMasses[atomj0]/(inverseMasses[atomj0]+inverseMasses[atomj1]);
                  else if (atomj1 == atomk0)
                      dot = -r[j].dot(r[k])*inverseMasses[atomj1]/(inverseMasses[atomj0]+inverseMasses[atomj1]);
                  else
                      dot = 0.0;
                  matrix[j][k] = dot;
              }
              matrix[j][j] = 1.0;
          }

          // Invert it using SVD.

          Array2D<double> u, v;
          Array1D<double> w;
          SVD<double> svd(matrix);
          svd.getU(u);
          svd.getV(v);
          svd.getSingularValues(w);
          double singularValueCutoff = 0.01*w[0];
          for (int j = 0; j < size; j++)
              w[j] = (w[j] < singularValueCutoff ? 0.0 : 1.0/w[j]);
          if (temp.size() < size)
              temp.resize(size);
          for (int j = 0; j < size; j++) {
              for (int k = 0; k < size; k++) {
                  matrix[j][k] = 0.0;
                  for (int m = 0; m < size; m++)
                      matrix[j][k] += v[j][m]*w[m]*u[k][m];
              }
          }

          // Record the inverted matrix.

          for (int j = 0; j < size; j++)
              for (int k = 0; k < size; k++)
                  _matrices[i][j][k] = matrix[j][k];
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

   int done                 = 0;
   int iterations           = 0;
   int numberConverged      = 0;
   vector<RealOpenMM> constraintForce(_numberOfConstraints);
   vector<RealOpenMM> tempForce(10);
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
         constraintForce[ii] = zero;
         RealOpenMM rrpr  = DOT3(  rp_ij, r_ij[ii] );
         RealOpenMM acor;
         if( rrpr <  d_ij2[ii]*epsilon6 ){
             std::stringstream message;
             message << iterations <<" "<<atomI<<" "<<atomJ<< " Error: sign of rrpr < 0?\n";
             SimTKOpenMMLog::printMessage( message );
         } else {
             constraintForce[ii] = reducedMasses[ii]*diff/rrpr;
         }
         if (rp2 >= lowerTol*dist2 && rp2 <= upperTol*dist2) {
            numberConverged++;
         }
      }
      if( numberConverged == _numberOfConstraints ){
         done = true;
      }
      for (int i = 0; i < _rigidClusters.size(); i++) {
          const vector<int>& cluster = _rigidClusters[i];
          RealOpenMM** matrix = _matrices[i];
          int size = cluster.size();
          if (size > tempForce.size())
              tempForce.resize(size);
          for (int j = 0; j < size; j++) {
              tempForce[j] = zero;
              for (int k = 0; k < size; k++)
                  tempForce[j] += matrix[j][k]*constraintForce[cluster[k]]*_distance[cluster[k]]/_distance[cluster[j]];
          }
          for (int j = 0; j < size; j++)
              constraintForce[cluster[j]] = tempForce[j];

      }
RealOpenMM damping = one;//(RealOpenMM) (iterations%2 == 0 ? 0.5 : 1.0);
      for( int ii = 0; ii < _numberOfConstraints; ii++ ){

         int atomI   = _atomIndices[ii][0];
         int atomJ   = _atomIndices[ii][1];
         for( int jj = 0; jj < 3; jj++ ){
            RealOpenMM dr                = constraintForce[ii]*r_ij[ii][jj]*damping;
            atomCoordinatesP[atomI][jj] += inverseMasses[atomI]*dr;
            atomCoordinatesP[atomJ][jj] -= inverseMasses[atomJ]*dr;
         }
      }
   }
   static int sum = 0;
   static int count = 0;
   sum += iterations;
   count++;
   if (count == 100) {
       printf("%d iterations\n", sum);
       sum = 0;
       count = 0;
   }

   // diagnostics

   if( debug || !done ){
      std::stringstream message;
      message << methodName;
      message << " iterations=" << iterations << " no. converged=" << numberConverged << " out of " << _numberOfConstraints;
      if( done ){
         message << " SUCCESS";
      } else {
         message << " FAILED";
      }
      message << "\n";
      int errors = reportShake( numberOfAtoms, atomCoordinatesP, message );
      if( !errors ){
         message << "*** no errors recorded in explicit check ***";
      }
      message << "\n";
      SimTKOpenMMLog::printMessage( message );
   }

   return (done ? ReferenceDynamics::DefaultReturn : ReferenceDynamics::ErrorReturn);

}

/**---------------------------------------------------------------------------------------

   Report any violated constriants

   @param numberOfAtoms    number of atoms
   @param atomCoordinates  atom coordinates
   @param message          report

   @return number of violated constraints

   --------------------------------------------------------------------------------------- */

int ReferenceRigidShakeAlgorithm::reportShake( int numberOfAtoms, RealOpenMM** atomCoordinates,
                                          std::stringstream& message ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nReferenceQShakeAlgorithm::reportShake";

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

