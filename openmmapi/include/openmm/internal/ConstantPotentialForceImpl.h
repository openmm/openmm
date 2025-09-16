#ifndef OPENMM_CONSTANTPOTENTIALFORCEIMPL_H_
#define OPENMM_CONSTANTPOTENTIALFORCEIMPL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2008-2025 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Evan Pretti                                        *
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

#include "ForceImpl.h"
#include "openmm/ConstantPotentialForce.h"
#include "openmm/Kernel.h"
#include <utility>
#include <set>
#include <string>

namespace OpenMM {

class System;

/**
 * This is the internal implementation of ConstantPotentialForce.
 */

class OPENMM_EXPORT ConstantPotentialForceImpl : public ForceImpl {
public:
    ConstantPotentialForceImpl(const ConstantPotentialForce& owner);
    ~ConstantPotentialForceImpl();
    void initialize(ContextImpl& context);
    const ConstantPotentialForce& getOwner() const {
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
    void updateParametersInContext(ContextImpl& context, int firstParticle, int lastParticle, int firstException, int lastException, int firstElectrode, int lastElectrode);
    void getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const;
    void getCharges(ContextImpl& context, std::vector<double>& charges);
    /**
     * This is a utility routine that calculates the values to use for alpha and
     * grid size when using Particle Mesh Ewald.
     */
    static void calcPMEParameters(const System& system, const ConstantPotentialForce& force, double& alpha, int& xsize, int& ysize, int& zsize);
private:
    const ConstantPotentialForce& owner;
    Kernel kernel;
};

} // namespace OpenMM

#endif /*OPENMM_CONSTANTPOTENTIALFORCEIMPL_H_*/
