#ifndef OPENMM_REFERENCESETTLEALGORITHM_H_
#define OPENMM_REFERENCESETTLEALGORITHM_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "ReferenceConstraintAlgorithm.h"
#include <vector>

namespace OpenMM {

/**
 * This implements the SETTLE algorithm for constraining water molecules.
 */

class OPENMM_EXPORT ReferenceSETTLEAlgorithm : public ReferenceConstraintAlgorithm {
public:
    ReferenceSETTLEAlgorithm(const std::vector<int>& atom1, const std::vector<int>& atom2, const std::vector<int>& atom3,
            const std::vector<RealOpenMM>& distance1, const std::vector<RealOpenMM>& distance2, std::vector<RealOpenMM>& masses, RealOpenMM tolerance);
    /**
     * Get the constraint tolerance.
     */
    RealOpenMM getTolerance() const;

    /**
     * Set the constraint tolerance.
     */
    void setTolerance(RealOpenMM tolerance);

    /**
     * Apply the constraint algorithm.
     * 
     * @param numberOfAtoms    number of atoms
     * @param atomCoordinates  the original atom coordinates
     * @param atomCoordinatesP the new atom coordinates
     * @param inverseMasses    1/mass
     */
    int apply(int numberOfAtoms, std::vector<OpenMM::RealVec>& atomCoordinates, std::vector<OpenMM::RealVec>& atomCoordinatesP, std::vector<RealOpenMM>& inverseMasses);

    /**
     * Apply the constraint algorithm to velocities.
     * 
     * @param numberOfAtoms    number of atoms
     * @param atomCoordinates  the atom coordinates
     * @param atomCoordinatesP the velocities to modify
     * @param inverseMasses    1/mass
     */
    int applyToVelocities(int numberOfAtoms, std::vector<OpenMM::RealVec>& atomCoordinates, std::vector<OpenMM::RealVec>& velocities, std::vector<RealOpenMM>& inverseMasses);
private:
    RealOpenMM tolerance;
    std::vector<int> atom1;
    std::vector<int> atom2;
    std::vector<int> atom3;
    std::vector<RealOpenMM> distance1;
    std::vector<RealOpenMM> distance2;
    std::vector<RealOpenMM> masses;
};

} // namespace OpenMM

#endif /*OPENMM_REFERENCESETTLEALGORITHM_H_*/
