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

#ifndef __ReferenceLincsAlgorithm_H__
#define __ReferenceLincsAlgorithm_H__

#include "ReferenceConstraintAlgorithm.h"
#include "SimTKOpenMMRealType.h"
#include <vector>

namespace OpenMM {

class ReferenceLincsAlgorithm : public ReferenceConstraintAlgorithm {

   protected:

      int _numTerms;
      int _numberOfConstraints;
      int** _atomIndices;
      double* _distance;

      bool _hasInitialized;
      std::vector<std::vector<int> > _linkedConstraints;
      std::vector<double> _sMatrix;
      std::vector<double> _rhs1;
      std::vector<double> _rhs2;
      std::vector<double> _solution;
      std::vector<std::vector<double> > _couplingMatrix;
      std::vector<OpenMM::Vec3> _constraintDir;


      /**---------------------------------------------------------------------------------------

         Set the number of terms to use in the series expansion

         @param numberOfAtoms    number of atoms
         @param inverseMasses    1/mass

         --------------------------------------------------------------------------------------- */

      void initialize(int numberOfAtoms, std::vector<double>& inverseMasses);

      /**---------------------------------------------------------------------------------------

         Solve the matrix equation

         --------------------------------------------------------------------------------------- */

      void solveMatrix();

      /**---------------------------------------------------------------------------------------

         Update the atom position based on the solution to the matrix equation.

         @param numberOfAtoms    number of atoms
         @param atomCoordinates  atom coordinates
         @param inverseMasses    1/mass

         --------------------------------------------------------------------------------------- */

      void updateAtomPositions(int numberOfAtoms, std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<double>& inverseMasses);

   public:

      /**---------------------------------------------------------------------------------------

         ReferenceLincsAlgorithm constructor

         @param numberOfConstraints      number of constraints
         @param atomIndices              atom indices for contraints
         @param distance                 distances for constraints

         --------------------------------------------------------------------------------------- */

      ReferenceLincsAlgorithm(int numberOfConstraints, int** atomIndices, double* distance);

      /**---------------------------------------------------------------------------------------

         Get number of constraints

         @return number of constraints

         --------------------------------------------------------------------------------------- */

      int getNumberOfConstraints() const;

      /**---------------------------------------------------------------------------------------

         Get the number of terms to use in the series expansion

         @return terms

         --------------------------------------------------------------------------------------- */

      int getNumTerms() const;

      /**---------------------------------------------------------------------------------------

         Set the number of terms to use in the series expansion

         --------------------------------------------------------------------------------------- */

      void setNumTerms(int terms);

      /**---------------------------------------------------------------------------------------

         Apply Lincs algorithm

         @param numberOfAtoms    number of atoms
         @param atomCoordinates  atom coordinates
         @param atomCoordinatesP atom coordinates prime
         @param inverseMasses    1/mass

         --------------------------------------------------------------------------------------- */

      void apply(int numberOfAtoms, std::vector<OpenMM::Vec3>& atomCoordinates,
                       std::vector<OpenMM::Vec3>& atomCoordinatesP, std::vector<double>& inverseMasses);

      /**---------------------------------------------------------------------------------------

         Apply constraint algorithm to velocities.

         @param numberOfAtoms    number of atoms
         @param atomCoordinates  atom coordinates
         @param velocities       atom velocities
         @param inverseMasses    1/mass

         --------------------------------------------------------------------------------------- */

      void applyToVelocities(int numberOfAtoms, std::vector<OpenMM::Vec3>& atomCoordinates,
                     std::vector<OpenMM::Vec3>& velocities, std::vector<double>& inverseMasses);
};

} // namespace OpenMM

#endif // __ReferenceLincsAlgorithm_H__
