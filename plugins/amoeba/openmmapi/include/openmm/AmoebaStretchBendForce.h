#ifndef OPENMM_AMOEBA_STRETCH_BEND_FORCE_H_
#define OPENMM_AMOEBA_STRETCH_BEND_FORCE_H_

/* -------------------------------------------------------------------------- *
 *                              OpenMMAmoeba                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
 * Authors: Mark Friedrichs, Peter Eastman                                    *
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

#include "openmm/Force.h"
#include "internal/windowsExportAmoeba.h"
#include <vector>

namespace OpenMM {

/**
 * This class implements the Amoeba stretch-bend interaction.
 * 
 * To use it, create a StretchBendForce object then call addStretchBend() once for each stretch-bend.  After
 * a stretch-bend has been added, you can modify its force field parameters by calling setStretchBendParameters().
 * This will have no effect on Contexts that already exist unless you call updateParametersInContext().
 */

class OPENMM_EXPORT_AMOEBA AmoebaStretchBendForce : public Force {

public:

    /**
     * Create an AmoebaStretchBendForce.
     */
    AmoebaStretchBendForce();

    /**
     * Get the number of stretch-bend terms in the potential function
     */
    int getNumStretchBends() const {
        return stretchBends.size();
    }

    /**
     * Add a stretch-bend term to the force field.
     *
     * @param particle1     the index of the first particle connected by the stretch-bend
     * @param particle2     the index of the second particle connected by the stretch-bend
     * @param particle3     the index of the third particle connected by the stretch-bend
     * @param lengthAB      the equilibrium length of the stretch-bend in bond ab [particle1, particle2], measured in nm
     * @param lengthCB      the equilibrium length of the stretch-bend in bond cb [particle3, particle2], measured in nm
     * @param angle         the equilibrium angle in radians
     * @param k1            the force constant of the product of bond ab and angle a-b-c
     * @param k2            the force constant of the product of bond bc and angle a-b-c (optional, default is the same as k1)
     * @return the index of the stretch-bend that was added
     */
    int addStretchBend(int particle1, int particle2, int particle3, double lengthAB,  double lengthCB, double angle,
                       double k1, double k2);

    /**
     * Get the force field parameters for a stretch-bend term.
     * 
     * @param index              the index of the stretch-bend for which to get parameters
     * @param[out] particle1     the index of the first particle connected by the stretch-bend
     * @param[out] particle2     the index of the second particle connected by the stretch-bend
     * @param[out] particle3     the index of the third particle connected by the stretch-bend
     * @param[out] lengthAB      the equilibrium length of the stretch-bend in bond ab [particle1, particle2], measured in nm
     * @param[out] lengthCB      the equilibrium length of the stretch-bend in bond cb [particle3, particle2], measured in nm
     * @param[out] angle         the equilibrium angle in radians
     * @param[out] k1            the force constant of the product of bond ab and angle a-b-c
     * @param[out] k2            the force constant of the product of bond bc and angle a-b-c
     */
    void getStretchBendParameters(int index, int& particle1, int& particle2, int& particle3, double& lengthAB,
                                  double& lengthCB, double& angle, double& k1, double& k2) const;

    /**
     * Set the force field parameters for a stretch-bend term.
     * 
     * @param index         the index of the stretch-bend for which to set parameters
     * @param particle1     the index of the first particle connected by the stretch-bend
     * @param particle2     the index of the second particle connected by the stretch-bend
     * @param particle3     the index of the third particle connected by the stretch-bend
     * @param lengthAB      the equilibrium length of the stretch-bend in bond ab [particle1, particle2], measured in nm
     * @param lengthCB      the equilibrium length of the stretch-bend in bond cb [particle3, particle2], measured in nm
     * @param angle         the equilibrium angle in radians
     * @param k1            the force constant of the product of bond ab and angle a-b-c
     * @param k2            the force constant of the product of bond bc and angle a-b-c (optional, default is the same as k1)
     */
    void setStretchBendParameters(int index, int particle1, int particle2, int particle3, 
                                  double lengthAB,  double lengthCB, double angle, double k1, double k2);
    /**
     * Update the per-stretch-bend term parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setStretchBendParameters() to modify this object's parameters, then call updateParametersInContext()
     * to copy them over to the Context.
     * 
     * The only information this method updates is the values of per-stretch-bend term parameters.  The set of particles involved
     * in a term cannot be changed, nor can new terms be added.
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
    bool usesPeriodicBoundaryConditions() const;
protected:
    ForceImpl* createImpl() const;
private:
    class StretchBendInfo;
    std::vector<StretchBendInfo> stretchBends;
    bool usePeriodic;
};

/**
 * This is an internal class used to record information about a stretch-bend.
 * @private
 */
class AmoebaStretchBendForce::StretchBendInfo {
public:
    int particle1, particle2, particle3;
    double lengthAB, lengthCB, angle, k1, k2;
    StretchBendInfo() {
        particle1 = particle2  = particle3 = -1;
        lengthAB  = lengthCB   = angle     = k1 = k2 = 0.0;
    }
    StretchBendInfo(int particle1, int particle2, int particle3, 
                    double lengthAB,  double lengthCB, double angle, double k1, double k2) :
                    particle1(particle1), particle2(particle2), particle3(particle3), 
                    lengthAB(lengthAB), lengthCB(lengthCB), angle(angle), k1(k1), k2(k2) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_AMOEBA_STRETCH_BEND_FORCE_H_*/
