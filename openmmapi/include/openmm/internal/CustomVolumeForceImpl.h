#ifndef OPENMM_CUSTOMVOLUMEFORCEIMPL_H_
#define OPENMM_CUSTOMVOLUMEFORCEIMPL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
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

#include "CustomCPPForceImpl.h"
#include "openmm/CustomVolumeForce.h"
#include "lepton/CompiledExpression.h"
#include <string>

namespace OpenMM {

/**
 * This is the internal implementation of CustomVolumeForce.
 */

class CustomVolumeForceImpl : public CustomCPPForceImpl {
public:
    CustomVolumeForceImpl(const CustomVolumeForce& owner);
    void initialize(ContextImpl& context);
    const CustomVolumeForce& getOwner() const {
        return owner;
    }
    double computeForce(ContextImpl& context, const std::vector<Vec3>& positions, std::vector<Vec3>& forces);
    std::map<std::string, double> getDefaultParameters();
private:
    const CustomVolumeForce& owner;
    std::map<std::string, double> defaultParameters;
    Lepton::CompiledExpression energyExpression;
    std::vector<std::string> globalParameterNames;
    std::vector<double> globalValues;
    Vec3 a, b, c;
    double volume;
};

} // namespace OpenMM

#endif /*OPENMM_CUSTOMVOLUMEFORCEIMPL_H_*/
