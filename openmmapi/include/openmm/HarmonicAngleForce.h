#ifndef OPENMM_HARMONICANGLEFORCE_H_
#define OPENMM_HARMONICANGLEFORCE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
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
#include <map>
#include <vector>
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This class implements an interaction between groups of three particles that varies harmonically with the angle
 * between them.  To use it, create a HarmonicAngleForce object then call addAngle() once for each angle.  After
 * an angle has been added, you can modify its force field parameters by calling setAngleParameters().  This will
 * have no effect on Contexts that already exist unless you call updateParametersInContext().
 */

class OPENMM_EXPORT HarmonicAngleForce : public Force {
public:
    /**
     * Create a HarmonicAngleForce.
     */
    HarmonicAngleForce();
    /**
     * Get the number of harmonic bond angle terms in the potential function
     */
    int getNumAngles() const {
        return angles.size();
    }
    /**
     * Add an angle term to the force field.
     *
     * @param particle1 the index of the first particle forming the angle
     * @param particle2 the index of the second particle forming the angle
     * @param particle3 the index of the third particle forming the angle
     * @param angle     the equilibrium angle, measured in radians
     * @param k         the harmonic force constant for the angle, measured in kJ/mol/radian^2
     * @return the index of the angle that was added
     */
    int addAngle(int particle1, int particle2, int particle3, double angle, double k);
    /**
     * Get the force field parameters for an angle term.
     *
     * @param      index     the index of the angle for which to get parameters
     * @param[out] particle1 the index of the first particle forming the angle
     * @param[out] particle2 the index of the second particle forming the angle
     * @param[out] particle3 the index of the third particle forming the angle
     * @param[out] angle     the equilibrium angle, measured in radians
     * @param[out] k         the harmonic force constant for the angle, measured in kJ/mol/radian^2
     */
    void getAngleParameters(int index, int& particle1, int& particle2, int& particle3, double& angle, double& k) const;
    /**
     * Set the force field parameters for an angle term.
     *
     * @param index     the index of the angle for which to set parameters
     * @param particle1 the index of the first particle forming the angle
     * @param particle2 the index of the second particle forming the angle
     * @param particle3 the index of the third particle forming the angle
     * @param angle     the equilibrium angle, measured in radians
     * @param k         the harmonic force constant for the angle, measured in kJ/mol/radian^2
     */
    void setAngleParameters(int index, int particle1, int particle2, int particle3, double angle, double k);
    /**
     * Update the per-angle parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setAngleParameters() to modify this object's parameters, then call updateParametersInContext()
     * to copy them over to the Context.
     *
     * The only information this method updates is the values of per-angle parameters.  The set of particles involved
     * in a angle cannot be changed, nor can new angles be added.
     */
    void updateParametersInContext(Context& context);
    /**
     * Set whether this force should apply periodic boundary conditions when calculating displacements.
     * Usually this is not appropriate for bonded forces, but there are situations when it can be useful.
     */
    void setUsesPeriodicBoundaryConditions(bool periodic);
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * @returns true if force uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const;protected:
    ForceImpl* createImpl() const;
private:
    class AngleInfo;
    std::vector<AngleInfo> angles;
    bool usePeriodic;
};

/**
 * This is an internal class used to record information about an angle.
 * @private
 */
class HarmonicAngleForce::AngleInfo {
public:
    int particle1, particle2, particle3;
    double angle, k;
    AngleInfo() {
        particle1 = particle2 = particle3 = -1;
        angle = k = 0.0;
    }
    AngleInfo(int particle1, int particle2, int particle3, double angle, double k) :
        particle1(particle1), particle2(particle2), particle3(particle3), angle(angle), k(k) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_HARMONICANGLEFORCE_H_*/
