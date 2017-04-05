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
            const std::vector<double>& distance1, const std::vector<double>& distance2, std::vector<double>& masses);

    /**
     * Get the number of clusters.
     */
    int getNumClusters() const;
    /**
     * Get the parameters describing one cluster.
     * 
     * @param index       the index of the cluster to get
     * @param atom1       the index of the first atom in the cluster
     * @param atom2       the index of the second atom in the cluster
     * @param atom3       the index of the third atom in the cluster
     * @param distance1   the distance between atoms 1 and 2
     * @param distance2   the distance between atoms 2 and 3
     */
    void getClusterParameters(int index, int& atom1, int& atom2, int& atom3, double& distance1, double& distance2) const;
    /**
     * Apply the constraint algorithm.
     * 
     * @param atomCoordinates  the original atom coordinates
     * @param atomCoordinatesP the new atom coordinates
     * @param inverseMasses    1/mass
     * @param tolerance        the constraint tolerance
     */
    void apply(std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<OpenMM::Vec3>& atomCoordinatesP, std::vector<double>& inverseMasses, double tolerance);

    /**
     * Apply the constraint algorithm to velocities.
     * 
     * @param atomCoordinates  the atom coordinates
     * @param atomCoordinatesP the velocities to modify
     * @param inverseMasses    1/mass
     * @param tolerance        the constraint tolerance
     */
    void applyToVelocities(std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<OpenMM::Vec3>& velocities, std::vector<double>& inverseMasses, double tolerance);
private:
    std::vector<int> atom1;
    std::vector<int> atom2;
    std::vector<int> atom3;
    std::vector<double> distance1;
    std::vector<double> distance2;
    std::vector<double> masses;
};

} // namespace OpenMM

#endif /*OPENMM_REFERENCESETTLEALGORITHM_H_*/
