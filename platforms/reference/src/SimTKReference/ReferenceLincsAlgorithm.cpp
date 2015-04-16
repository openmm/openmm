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

#include "SimTKOpenMMUtilities.h"
#include "ReferenceLincsAlgorithm.h"
#include "ReferenceDynamics.h"
#include "openmm/OpenMMException.h"

using std::vector;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceLincsAlgorithm constructor

   @param numberOfConstraints      number of constraints
   @param atomIndices              block of atom indices
   @param distance                 distances for constraints

   --------------------------------------------------------------------------------------- */

ReferenceLincsAlgorithm::ReferenceLincsAlgorithm(int numberOfConstraints,
                                                 int** atomIndices,
                                                 RealOpenMM* distance) {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceLincsAlgorithm::ReferenceLincsAlgorithm";

   // ---------------------------------------------------------------------------------------

   _numberOfConstraints        = numberOfConstraints;
   _atomIndices                = atomIndices;
   _distance                   = distance;

   _numTerms                   = 4;
   _hasInitialized             = false;
}

/**---------------------------------------------------------------------------------------

   Get number of constraints

   @return number of constraints

   --------------------------------------------------------------------------------------- */

int ReferenceLincsAlgorithm::getNumberOfConstraints() const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceLincsAlgorithm::getNumberOfConstraints";

   // ---------------------------------------------------------------------------------------

   return _numberOfConstraints;
}

/**---------------------------------------------------------------------------------------

   Get the number of terms to use in the series expansion

   @return terms

   --------------------------------------------------------------------------------------- */

int ReferenceLincsAlgorithm::getNumTerms() const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceLincsAlgorithm::getNumTerms";

   // ---------------------------------------------------------------------------------------

   return _numTerms;
}

/**---------------------------------------------------------------------------------------

   Set the number of terms to use in the series expansion

   --------------------------------------------------------------------------------------- */

void ReferenceLincsAlgorithm::setNumTerms(int terms) {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceLincsAlgorithm::setNumTerms";

   // ---------------------------------------------------------------------------------------

   _numTerms = terms;
}


/**---------------------------------------------------------------------------------------

   Initialize internal data structures.

   @param numberOfAtoms    number of atoms
   @param inverseMasses    1/mass

   --------------------------------------------------------------------------------------- */

void ReferenceLincsAlgorithm::initialize(int numberOfAtoms, vector<RealOpenMM>& inverseMasses) {
    static const RealOpenMM one = 1.0;
    _hasInitialized = true;
    vector<vector<int> > atomConstraints(numberOfAtoms);
    for (int constraint = 0; constraint < _numberOfConstraints; constraint++) {
        atomConstraints[_atomIndices[constraint][0]].push_back(constraint);
        atomConstraints[_atomIndices[constraint][1]].push_back(constraint);
    }
    _linkedConstraints.resize(_numberOfConstraints);
    for (int atom = 0; atom < numberOfAtoms; atom++) {
        for (int i = 0; i < (int)atomConstraints[atom].size(); i++)
            for (int j = 0; j < i; j++) {
                int c1 = atomConstraints[atom][i];
                int c2 = atomConstraints[atom][j];
                _linkedConstraints[c1].push_back(c2);
                _linkedConstraints[c2].push_back(c1);
            }
    }
    _sMatrix.resize(_numberOfConstraints);
    for (int constraint = 0; constraint < _numberOfConstraints; constraint++)
        _sMatrix[constraint] = one/SQRT(inverseMasses[_atomIndices[constraint][0]]+inverseMasses[_atomIndices[constraint][1]]);
    _couplingMatrix.resize(_numberOfConstraints);
    for (int constraint = 0; constraint < _numberOfConstraints; constraint++)
        _couplingMatrix[constraint].resize(_linkedConstraints[constraint].size());
    _constraintDir.resize(_numberOfConstraints);
    _rhs1.resize(_numberOfConstraints);
    _rhs2.resize(_numberOfConstraints);
    _solution.resize(_numberOfConstraints);
}

/**---------------------------------------------------------------------------------------

 Solve the matrix equation

 --------------------------------------------------------------------------------------- */

void ReferenceLincsAlgorithm::solveMatrix() {
    static const RealOpenMM zero = 0.0;
    for (int iteration = 0; iteration < _numTerms; iteration++) {
        vector<RealOpenMM>& rhs1 = (iteration%2 == 0 ? _rhs1 : _rhs2);
        vector<RealOpenMM>& rhs2 = (iteration%2 == 0 ? _rhs2 : _rhs1);
        for (int c1 = 0; c1 < _numberOfConstraints; c1++) {
            rhs2[c1] = zero;
            for (int j = 0; j < (int)_linkedConstraints[c1].size(); j++) {
                int c2 = _linkedConstraints[c1][j];
                rhs2[c1] += _couplingMatrix[c1][j]*rhs1[c2];
            }
            _solution[c1] += rhs2[c1];
        }
    }
}

/**---------------------------------------------------------------------------------------

 Update the atom position based on the solution to the matrix equation.

 @param numberOfAtoms    number of atoms
 @param atomCoordinates  atom coordinates
 @param inverseMasses    1/mass

 --------------------------------------------------------------------------------------- */

void ReferenceLincsAlgorithm::updateAtomPositions(int numberOfAtoms, vector<RealVec>& atomCoordinates, vector<RealOpenMM>& inverseMasses) {
    for (int i = 0; i < _numberOfConstraints; i++) {
        RealVec delta(_sMatrix[i]*_solution[i]*_constraintDir[i][0],
                   _sMatrix[i]*_solution[i]*_constraintDir[i][1],
                   _sMatrix[i]*_solution[i]*_constraintDir[i][2]);
        int atom1 = _atomIndices[i][0];
        int atom2 = _atomIndices[i][1];
        atomCoordinates[atom1][0] -= (RealOpenMM)(inverseMasses[atom1]*delta[0]);
        atomCoordinates[atom1][1] -= (RealOpenMM)(inverseMasses[atom1]*delta[1]);
        atomCoordinates[atom1][2] -= (RealOpenMM)(inverseMasses[atom1]*delta[2]);
        atomCoordinates[atom2][0] += (RealOpenMM)(inverseMasses[atom2]*delta[0]);
        atomCoordinates[atom2][1] += (RealOpenMM)(inverseMasses[atom2]*delta[1]);
        atomCoordinates[atom2][2] += (RealOpenMM)(inverseMasses[atom2]*delta[2]);
    }
}

/**---------------------------------------------------------------------------------------

   Apply Lincs algorithm

   @param numberOfAtoms    number of atoms
   @param atomCoordinates  atom coordinates
   @param atomCoordinatesP atom coordinates prime
   @param inverseMasses    1/mass

   --------------------------------------------------------------------------------------- */

void ReferenceLincsAlgorithm::apply(int numberOfAtoms, vector<RealVec>& atomCoordinates,
                                         vector<RealVec>& atomCoordinatesP,
                                         vector<RealOpenMM>& inverseMasses) {

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nReferenceLincsAlgorithm::apply";

   static const RealOpenMM zero        =  0.0;
   static const RealOpenMM one         =  1.0;
   static const RealOpenMM two         =  2.0;

   // ---------------------------------------------------------------------------------------

   if (_numberOfConstraints == 0)
       return;

   if (!_hasInitialized)
       initialize(numberOfAtoms, inverseMasses);

   // Calculate the direction of each constraint, along with the initial RHS and solution vectors.

   for (int i = 0; i < _numberOfConstraints; i++) {
       int atom1 = _atomIndices[i][0];
       int atom2 = _atomIndices[i][1];
       _constraintDir[i] = RealVec(atomCoordinatesP[atom1][0]-atomCoordinatesP[atom2][0],
                                atomCoordinatesP[atom1][1]-atomCoordinatesP[atom2][1],
                                atomCoordinatesP[atom1][2]-atomCoordinatesP[atom2][2]);
       RealOpenMM invLength = (RealOpenMM)(1/SQRT((RealOpenMM)_constraintDir[i].dot(_constraintDir[i])));
       _constraintDir[i][0] *= invLength;
       _constraintDir[i][1] *= invLength;
       _constraintDir[i][2] *= invLength;
       _rhs1[i] = _solution[i] = _sMatrix[i]*(one/invLength-_distance[i]);
   }

   // Build the coupling matrix.

   for (int c1 = 0; c1 < (int)_couplingMatrix.size(); c1++) {
       RealVec& dir1 = _constraintDir[c1];
       for (int j = 0; j < (int)_couplingMatrix[c1].size(); j++) {
           int c2 = _linkedConstraints[c1][j];
           RealVec& dir2 = _constraintDir[c2];
           if (_atomIndices[c1][0] == _atomIndices[c2][0] || _atomIndices[c1][1] == _atomIndices[c2][1])
               _couplingMatrix[c1][j] = (RealOpenMM)(-inverseMasses[_atomIndices[c1][0]]*_sMatrix[c1]*dir1.dot(dir2)*_sMatrix[c2]);
           else
               _couplingMatrix[c1][j] = (RealOpenMM)(inverseMasses[_atomIndices[c1][1]]*_sMatrix[c1]*dir1.dot(dir2)*_sMatrix[c2]);
       }
   }

   // Solve the matrix equation and update the positions.

   solveMatrix();
   updateAtomPositions(numberOfAtoms, atomCoordinatesP, inverseMasses);

   // Correct for rotational lengthening.

   for (int i = 0; i < _numberOfConstraints; i++) {
       int atom1 = _atomIndices[i][0];
       int atom2 = _atomIndices[i][1];
       RealVec delta(atomCoordinatesP[atom1][0]-atomCoordinatesP[atom2][0],
                  atomCoordinatesP[atom1][1]-atomCoordinatesP[atom2][1],
                  atomCoordinatesP[atom1][2]-atomCoordinatesP[atom2][2]);
       RealOpenMM p2 = (RealOpenMM)(two*_distance[i]*_distance[i]-delta.dot(delta));
       if (p2 < zero)
           p2 = zero;
       _rhs1[i] = _solution[i] = _sMatrix[i]*(_distance[i]-SQRT(p2));
   }
   solveMatrix();
   updateAtomPositions(numberOfAtoms, atomCoordinatesP, inverseMasses);
}

/**---------------------------------------------------------------------------------------

   Apply constraint algorithm to velocities.

   @param numberOfAtoms    number of atoms
   @param atomCoordinates  atom coordinates
   @param velocities       atom velocities
   @param inverseMasses    1/mass

   --------------------------------------------------------------------------------------- */

void ReferenceLincsAlgorithm::applyToVelocities(int numberOfAtoms, std::vector<OpenMM::RealVec>& atomCoordinates,
               std::vector<OpenMM::RealVec>& velocities, std::vector<RealOpenMM>& inverseMasses) {
    throw OpenMM::OpenMMException("applyToVelocities is not implemented");
}
