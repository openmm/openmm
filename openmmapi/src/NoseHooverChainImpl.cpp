/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2010 Stanford University and the Authors.      *
 * Authors: Andreas Krämer and Andrew C. Simmonett                            *
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

#include "openmm/internal/NoseHooverChainImpl.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/Integrator.h"
#include "openmm/System.h"
#include "openmm/kernels.h"
#include "SimTKOpenMMRealType.h"
#include "sfmt/SFMT.h"
#include <string>
#include <vector>
#include <iostream>

using namespace OpenMM;
using std::vector;

NoseHooverChainImpl::NoseHooverChainImpl(const NoseHooverChain& owner) : owner(owner) {
}

void NoseHooverChainImpl::initialize(ContextImpl& context) {
}

void NoseHooverChainImpl::updateContextState(ContextImpl& context, bool& forcesInvalid) {
}

double NoseHooverChainImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    return 0;
}

std::map<std::string, double> NoseHooverChainImpl::getDefaultParameters() {
    std::map<std::string, double> parameters;
    const auto &owner = getOwner();
    int chainLength = owner.getDefaultChainLength();
    for(int i = 0; i < chainLength; ++i) {
        parameters[owner.Position(i)] = 0;
        parameters[owner.Velocity(i)] = 0;
    }
    return parameters;
}

std::vector<std::string> NoseHooverChainImpl::getKernelNames() {
    return std::vector<std::string>();
}
