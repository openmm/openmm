#ifndef OPENMM_AMOEBA_GK_FORCE_FIELD_IMPL_H_
#define OPENMM_AMOEBA_GK_FORCE_FIELD_IMPL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "openmm/internal/ForceImpl.h"
#include "openmm/AmoebaGeneralizedKirkwoodForce.h"
#include "openmm/Kernel.h"
#include <string>

namespace OpenMM {

/**
 * This is the internal implementation of AmoebaGeneralizedKirkwoodForce.
 */

class OPENMM_EXPORT_AMOEBA AmoebaGeneralizedKirkwoodForceImpl : public ForceImpl {
public:
    AmoebaGeneralizedKirkwoodForceImpl(const AmoebaGeneralizedKirkwoodForce& owner);
    void initialize(ContextImpl& context);
    const AmoebaGeneralizedKirkwoodForce& getOwner() const {
        return owner;
    }
    void updateContextState(ContextImpl& context, bool& forcesInvalid) {
        // This force field doesn't update the state directly.
    }
    double calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups);
    std::map<std::string, double> getDefaultParameters() {
        return std::map<std::string, double>(); // This force field doesn't define any parameters.
    }
    std::vector<std::string> getKernelNames();
    Kernel& getKernel() {
        return kernel;
    }
    void updateParametersInContext(ContextImpl& context);

    /**
     * Compute a "neck" descreening contribution.
     * <p>
     * Using the tabulated Aij and Bij values, compute the neck integral
     * using Equations 13 and 14 from Aguilar, Shadrach and Onufriev (10.1021/ct100392h).
     *
     * @param r separation distance.
     * @param radius Radius of the current atom.
     * @param radiusK Radius of the descreening atom.
     * @param sneck The neck scale factor.
     * @return The integral value.
     */
    static double neckDescreen(double r, double radius, double radiusK, double sneck);

    /**
     * Get Neck Aij and Bij constants.
     *
     * Based on the radii of two input atoms ("radius" for the atom being descreened
     * and "radiusK" for the descreening atom), the neck parameters Aij and Bij are
     * interpolated from tabulated values. The definitions of Aij and Bij are defined by
     * Eq. 11 from Corrigan et al. (10.1063/5.0158914) while Figure 2 describes
     * a neck descreening region.
     *
     * @param radius Radius of the current atom.
     * @param radiusK Radius of the descreening atom.
     * @param aij The Aij neck constant.
     * @param bij the Bij neck constant.
     */
    static void getNeckConstants(double radius, double radiusK, double &aij, double &bij);

     /**
      * The array of neck tabulated radii values.
      */
    static const std::vector<float>& getNeckRadii();

    const static int NUM_NECK_RADII = 45;

    /**
     * The tabulated Aij parameters.
     */
    const static float (&getAij())[NUM_NECK_RADII][NUM_NECK_RADII];

     /**
      * The tabulated Bij parameters.
      */
     const static float (&getBij())[NUM_NECK_RADII][NUM_NECK_RADII];

private:
    const AmoebaGeneralizedKirkwoodForce& owner;
    Kernel kernel;
};

} // namespace OpenMM

#endif /*OPENMM_AMOEBA_GBSA_OBC_FORCE_FIELD_IMPL_H_*/
