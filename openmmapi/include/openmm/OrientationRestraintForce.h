#ifndef OPENMM_ORIENTATIONRESTRAINTFORCE_H_
#define OPENMM_ORIENTATIONRESTRAINTFORCE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
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

#include "Force.h"
#include "Vec3.h"
#include <vector>
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This is an approximately harmonic restraint force for keeping a group of particles
 * in a particular orientation.  You specify a set of particles and reference positions
 * for them.  It finds the rigid transformation that optimally aligns the particles
 * to the reference positions, consisting of a translation followed by a rotation.
 * The interaction energy is 2*k*sin^2(theta/2), where theta is the rotation angle
 * and k is a user defined force constant.  When the particles are close to the target
 * orientation, this approximately equals the harmonic form (k/2)*theta^2.  For larger
 * angles, the force is reduced from the harmonic form, and it goes to zero at theta=pi/2.
 * This avoids the discontinuity in the force that would be present if the force were
 * nonzero at pi/2.
 */

class OPENMM_EXPORT OrientationRestraintForce : public Force {
public:
    /**
     * Create an OrientationRestraintForce.
     *
     * @param k                   the force constant, measured in kJ/mol
     * @param referencePositions  the reference positions to compute the rotation
     *                            from.  The length of this vector must equal the
     *                            number of particles in the system, even if not
     *                            all particles are used in computing the rotation.
     * @param particles           the indices of the particles to use when computing
     *                            the rotation.  If this is empty (the default), all
     *                            particles in the system will be used.
     */
    explicit OrientationRestraintForce(double k, const std::vector<Vec3>& referencePositions,
                       const std::vector<int>& particles=std::vector<int>());
    /**
     * Get the force constant.
     * 
     * @returns the force constant in kJ/mol
     */
    double getK() const {
        return k;
    }
    /**
     * Set the force constant.
     * 
     * @param k    the force constant, measured in kJ/mol
     */
    void setK(double k);
    /**
     * Get the reference positions to compute the rotation from.
     */
    const std::vector<Vec3>& getReferencePositions() const {
        return referencePositions;
    }
    /**
     * Set the reference positions to compute the rotation from.
     */
    void setReferencePositions(const std::vector<Vec3>& positions);
    /**
     * Get the indices of the particles to use when computing the rotation.  If this
     * is empty, all particles in the system will be used.
     */
    const std::vector<int>& getParticles() const {
        return particles;
    }
    /**
     * Set the indices of the particles to use when computing the rotation.  If this
     * is empty, all particles in the system will be used.
     */
    void setParticles(const std::vector<int>& particles);
    /**
     * Update the reference positions and particle indices in a Context to match those stored
     * in this Force object.  This method provides an efficient method to update certain parameters
     * in an existing Context without needing to reinitialize it.  Simply call setReferencePositions()
     * setParticles(), and setK() to modify this object's parameters, then call updateParametersInContext()
     * to copy them over to the Context.
     */
    void updateParametersInContext(Context& context);
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * @returns true if force uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const {
        return false;
    }
protected:
    ForceImpl* createImpl() const;
private:
    double k;
    std::vector<Vec3> referencePositions;
    std::vector<int> particles;
};

} // namespace OpenMM

#endif /*OPENMM_ORIENTATIONRESTRAINTFORCE_H_*/
