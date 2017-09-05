#ifndef OPENMM_AMOEBA_WCA_DISPERSION_FORCE_IMPL_H_
#define OPENMM_AMOEBA_WCA_DISPERSION_FORCE_IMPL_H_

/* -------------------------------------------------------------------------- *
 *                                OpenMMAmoeba                                *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors:                                                                   *
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
#include "openmm/AmoebaWcaDispersionForce.h"
#include "openmm/Kernel.h"
#include <utility>
#include <set>
#include <string>

namespace OpenMM {

/**
 * This is the internal implementation of AmoebaWcaDispersionForce.
 */

class OPENMM_EXPORT_AMOEBA AmoebaWcaDispersionForceImpl : public ForceImpl {
public:
    AmoebaWcaDispersionForceImpl(const AmoebaWcaDispersionForce& owner);
    ~AmoebaWcaDispersionForceImpl();
    void initialize(ContextImpl& context);
    const AmoebaWcaDispersionForce& getOwner() const {
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

    /** 
     * Get the maximum dispersion energy for a particle
     * 
     * @param force               AmoebaWcaDispersionForce reference
     * @param particleIndex       the particle index
     * @param maxDispersionEnergy maximum dispersion energy
     */
    static void getMaximumDispersionEnergy(const AmoebaWcaDispersionForce& force, int particleIndex, double& maxDispersionEnergy);

    /** 
     * Get the total maximum dispersion energy
     * 
     * @param force               AmoebaWcaDispersionForce reference
     *
     * @return total maximum dispersion energy for the system
     */
    static double getTotalMaximumDispersionEnergy(const AmoebaWcaDispersionForce& force);
    void updateParametersInContext(ContextImpl& context);

private:
    const AmoebaWcaDispersionForce& owner;
    Kernel kernel;
};

} // namespace OpenMM

#endif /*OPENMM_AMOEBA_WCA_DISPERSION_FORCE_IMPL_H_*/
