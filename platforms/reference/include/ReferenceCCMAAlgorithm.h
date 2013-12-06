
/* Portions copyright (c) 2006-2013 Stanford University and Simbios.
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

    int _numberOfConstraints;
    std::vector<std::pair<int, int> > _atomIndices;
    std::vector<RealOpenMM> _distance;

    std::vector<OpenMM::RealVec> _r_ij;
    RealOpenMM* _d_ij2;
    RealOpenMM* _distanceTolerance;
    RealOpenMM* _reducedMasses;
    bool _hasInitializedMasses;
    std::vector<std::vector<std::pair<int, RealOpenMM> > > _matrix;

private:

    void applyConstraints(std::vector<OpenMM::RealVec>& atomCoordinates,
                       std::vector<OpenMM::RealVec>& atomCoordinatesP, std::vector<RealOpenMM>& inverseMasses, bool constrainingVelocities, RealOpenMM tolerance);
          
public:
    class AngleInfo;

    /**
     * Create a ReferenceCCMAAlgorithm object.
     * 
     * @param numberOfAtoms            the number of atoms in the system
     * @param numberOfConstraints      the number of constraints
     * @param atomIndices              atom indices for contraints
     * @param distance                 distances for constraints
     * @param masses                   atom masses
     * @param angles                   angle force field terms
     */
    ReferenceCCMAAlgorithm(int numberOfAtoms, int numberOfConstraints, const std::vector<std::pair<int, int> >& atomIndices, const std::vector<RealOpenMM>& distance, std::vector<RealOpenMM>& masses, std::vector<AngleInfo>& angles);

    ~ReferenceCCMAAlgorithm( );

    /**
     * Get the number of constraints.
     */
    int getNumberOfConstraints( void ) const;

    /**
     * Get the maximum number of iterations to perform.
     */
    int getMaximumNumberOfIterations( void ) const;

    /**
     * Set the maximum number of iterations to perform.
     */
    void setMaximumNumberOfIterations( int maximumNumberOfIterations );

    /**
     * Apply the constraint algorithm.
     * 
     * @param atomCoordinates  the original atom coordinates
     * @param atomCoordinatesP the new atom coordinates
     * @param inverseMasses    1/mass
     * @param tolerance        the constraint tolerance
     */
    void apply(std::vector<OpenMM::RealVec>& atomCoordinates,
                       std::vector<OpenMM::RealVec>& atomCoordinatesP, std::vector<RealOpenMM>& inverseMasses, RealOpenMM tolerance);

    /**
     * Apply the constraint algorithm to velocities.
     * 
     * @param atomCoordinates  the atom coordinates
     * @param atomCoordinatesP the velocities to modify
     * @param inverseMasses    1/mass
     * @param tolerance        the constraint tolerance
     */
    void applyToVelocities(std::vector<OpenMM::RealVec>& atomCoordinates,
                     std::vector<OpenMM::RealVec>& velocities, std::vector<RealOpenMM>& inverseMasses, RealOpenMM tolerance);
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
