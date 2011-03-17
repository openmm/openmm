
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

#ifndef __ReferenceCCMAAlgorithm_H__
#define __ReferenceCCMAAlgorithm_H__

#include "ReferenceConstraintAlgorithm.h"
#include <utility>
#include <vector>
#include <set>

// ---------------------------------------------------------------------------------------

class OPENMM_EXPORT ReferenceCCMAAlgorithm : public ReferenceConstraintAlgorithm {

   protected:

      int _maximumNumberOfIterations;
      RealOpenMM _tolerance;

      int _numberOfConstraints;
      int** _atomIndices;
      RealOpenMM* _distance;

      RealOpenMM** _r_ij;
      RealOpenMM* _d_ij2;
      RealOpenMM* _distanceTolerance;
      RealOpenMM* _reducedMasses;
      bool _hasInitializedMasses;
      std::vector<std::vector<std::pair<int, RealOpenMM> > > _matrix;

   public:
      class AngleInfo;

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

      ReferenceCCMAAlgorithm( int numberOfAtoms, int numberOfConstraints, int** atomIndices, RealOpenMM* distance, std::vector<RealOpenMM>& masses, std::vector<AngleInfo>& angles, RealOpenMM tolerance );

      /**---------------------------------------------------------------------------------------

         Destructor

         --------------------------------------------------------------------------------------- */

       ~ReferenceCCMAAlgorithm( );

      /**---------------------------------------------------------------------------------------

         Get number of constraints

         @return number of constraints

         --------------------------------------------------------------------------------------- */

      int getNumberOfConstraints( void ) const;

      /**---------------------------------------------------------------------------------------

         Get maximum number of iterations

         @return maximum number of iterations

         --------------------------------------------------------------------------------------- */

      int getMaximumNumberOfIterations( void ) const;

      /**---------------------------------------------------------------------------------------

         Set maximum number of iterations

         @param maximumNumberOfIterations   new maximum number of iterations

         --------------------------------------------------------------------------------------- */

      void setMaximumNumberOfIterations( int maximumNumberOfIterations );

      /**---------------------------------------------------------------------------------------

         Get tolerance

         @return tolerance

         --------------------------------------------------------------------------------------- */

      RealOpenMM getTolerance( void ) const;

      /**---------------------------------------------------------------------------------------

         Set tolerance

         @param tolerance new tolerance

         --------------------------------------------------------------------------------------- */

      void setTolerance( RealOpenMM tolerance );

      /**---------------------------------------------------------------------------------------

         Apply CCMA algorithm

         @param numberOfAtoms    number of atoms
         @param atomCoordinates  atom coordinates
         @param atomCoordinatesP atom coordinates prime
         @param inverseMasses    1/mass

         @return SimTKOpenMMCommon::DefaultReturn if converge; else
          return SimTKOpenMMCommon::ErrorReturn

         --------------------------------------------------------------------------------------- */

      int apply( int numberOfAtoms, std::vector<OpenMM::RealVec>& atomCoordinates,
                       std::vector<OpenMM::RealVec>& atomCoordinatesP, std::vector<RealOpenMM>& inverseMasses );

      /**---------------------------------------------------------------------------------------

         Report any violated constraints

         @param numberOfAtoms    number of atoms
         @param atomCoordinates  atom coordinates
         @param message          report

         @return number of violated constraints

         --------------------------------------------------------------------------------------- */

      int reportCCMA( int numberOfAtoms, std::vector<OpenMM::RealVec>& atomCoordinates, std::stringstream& message );
};

class ReferenceCCMAAlgorithm::AngleInfo
{
public:
    int atom1, atom2, atom3;
    RealOpenMM angle;
    AngleInfo(int atom1, int atom2, int atom3, RealOpenMM angle) :
        atom1(atom1), atom2(atom2), atom3(atom3), angle(angle)
    {
    }
};

// ---------------------------------------------------------------------------------------

#endif // __ReferenceCCMAAlgorithm_H__
