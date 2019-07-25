#ifndef OPENMM_AMOEBA_OUT_OF_PLANE_BEND_FORCE_H_
#define OPENMM_AMOEBA_OUT_OF_PLANE_BEND_FORCE_H_

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
 * This class implements the Amoeba out-of-plane bend interaction.
 * 
 * To use it, create an OutOfPlaneBendForce object then call addOutOfPlaneBend() once for each outOfPlaneBend.  After
 * an out-of-plane bend has been added, you can modify its force field parameters by calling setOutOfPlaneBendParameters().
 * This will have no effect on Contexts that already exist unless you call updateParametersInContext().
 */

class OPENMM_EXPORT_AMOEBA AmoebaOutOfPlaneBendForce : public Force {

public:

    /**
     * Create an AmoebaOutOfPlaneBendForce.
     */
    AmoebaOutOfPlaneBendForce();

    /**
     * Get the number of out-of-plane bend terms in the potential function
     */
    int getNumOutOfPlaneBends() const {
        return outOfPlaneBends.size();
    }

    /** 
     * Set the global cubic term
     * 
     * @param cubicK        the cubic force constant for the angle
     */
    void setAmoebaGlobalOutOfPlaneBendCubic(double cubicK);

    /** 
     * Get the global cubic term
     * 
     * @return global cubicK term
     */
    double getAmoebaGlobalOutOfPlaneBendCubic() const;

    /** 
     * Set the global cubic term
     * 
     * @param quarticK       the quartic force constant for the angle
     */
    void setAmoebaGlobalOutOfPlaneBendQuartic(double quarticK);

    /** 
     * Get the global quartic term
     * 
     * @return global  quartic term
     */
    double getAmoebaGlobalOutOfPlaneBendQuartic() const;

    /** 
     * Set the global pentic term
     * 
     * @param penticK the pentic force constant for the angle
     */
    void setAmoebaGlobalOutOfPlaneBendPentic(double penticK);

    /** 
     * Get the global pentic term
     * 
     * @return global penticK term
     */
    double getAmoebaGlobalOutOfPlaneBendPentic() const;

    /** 
     * Set the global sextic term
     * 
     * @param sexticK       the sextic force constant for the angle
     */
    void setAmoebaGlobalOutOfPlaneBendSextic(double sexticK);

    /** 
     * Get the global sextic term
     * 
     * @return global sexticK term
     */
    double getAmoebaGlobalOutOfPlaneBendSextic() const;

    /**
     * Add an out-of-plane bend term to the force field.
     *
     * @param particle1     the index of the first particle connected by the outOfPlaneBend
     * @param particle2     the index of the second particle connected by the outOfPlaneBend
     * @param particle3     the index of the third particle connected by the outOfPlaneBend
     * @param particle4     the index of the fourth particle connected by the outOfPlaneBend
     * @param k             the force constant for the out-of-plane bend
     * @return the index of the out-of-plane bend that was added
     */
    int addOutOfPlaneBend(int particle1, int particle2, int particle3, int particle4, double k);

    /**
     * Get the force field parameters for an out-of-plane bend term.
     * 
     * @param index              the index of the outOfPlaneBend for which to get parameters
     * @param[out] particle1     the index of the first particle connected by the outOfPlaneBend
     * @param[out] particle2     the index of the second particle connected by the outOfPlaneBend
     * @param[out] particle3     the index of the third particle connected by the outOfPlaneBend
     * @param[out] particle4     the index of the fourth particle connected by the outOfPlaneBend
     * @param[out] k             the force constant for the out-of-plane bend
     */
    void getOutOfPlaneBendParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4, double& k) const;

    /**
     * Set the force field parameters for an out-of-plane bend term.
     * 
     * @param index         the index of the outOfPlaneBend for which to set parameters
     * @param particle1     the index of the first particle connected by the outOfPlaneBend
     * @param particle2     the index of the second particle connected by the outOfPlaneBend
     * @param particle3     the index of the third particle connected by the outOfPlaneBend
     * @param particle4     the index of the fourth particle connected by the outOfPlaneBend
     * @param k             the force constant for the out-of-plane bend
     */
    void setOutOfPlaneBendParameters(int index, int particle1, int particle2, int particle3, int particle4, double k);
    /**
     * Update the per-bend term parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setOutOfPlaneBendParameters() to modify this object's parameters, then call updateParametersInContext()
     * to copy them over to the Context.
     * 
     * The only information this method updates is the values of per-bend term parameters.  The set of particles involved
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
    double _globalCubicK, _globalQuarticK, _globalPenticK, _globalSexticK;
private:
    class OutOfPlaneBendInfo;
    std::vector<OutOfPlaneBendInfo> outOfPlaneBends;
    bool usePeriodic;
};

/**
 * This is an internal class used to record information about a bend.
 * @private
 */
class AmoebaOutOfPlaneBendForce::OutOfPlaneBendInfo {
public:
    int particle1, particle2, particle3, particle4;
    double k;
    OutOfPlaneBendInfo() {
        particle1 = particle2  = particle3 = particle4 = -1;
        k   = 0.0;
    }
    OutOfPlaneBendInfo(int particle1, int particle2, int particle3, int particle4,
                           double k) :
                    particle1(particle1), particle2(particle2), particle3(particle3), particle4(particle4), k(k) {
     
    }
};

} // namespace OpenMM

#endif /*OPENMM_AMOEBA_OUT_OF_PLANE_BEND_FORCE_H_*/
