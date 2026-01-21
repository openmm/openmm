#ifndef OPENMM_BUSSITHERMOSTATIMPL_H_
#define OPENMM_BUSSITHERMOSTATIMPL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
 * Authors: Muhammad Hasyim                                                   *
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
#include "openmm/BussiThermostat.h"
#include "openmm/Kernel.h"
#include <string>
#include <vector>

namespace OpenMM {

class System;

/**
 * This is the internal implementation of BussiThermostat.
 */
class OPENMM_EXPORT BussiThermostatImpl : public ForceImpl {
public:
    BussiThermostatImpl(const BussiThermostat& owner);
    void initialize(ContextImpl& context);
    const BussiThermostat& getOwner() const {
        return owner;
    }
    void updateContextState(ContextImpl& context, bool& forcesInvalid);
    double calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
        // This force doesn't apply forces to particles (only rescales velocities).
        return 0.0;
    }
    std::map<std::string, double> getDefaultParameters();
    std::vector<std::string> getKernelNames();
    /**
     * Get the list of particle indices that this thermostat applies to.
     * If the thermostat is set to apply to all particles, this returns all 
     * non-zero-mass particles.
     */
    const std::vector<int>& getParticleIndices() const {
        return particleIndices;
    }
private:
    const BussiThermostat& owner;
    Kernel kernel;
    std::vector<int> particleIndices;
    int stepCount;
};

} // namespace OpenMM

#endif /*OPENMM_BUSSITHERMOSTATIMPL_H_*/
