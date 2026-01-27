/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                             *
 * Portions copyright (c) 2025 Stanford University and the Authors.            *
 * Authors: Muhammad Hasyim                                                    *
 * Contributors:                                                               *
 *                                                                             *
 * Permission is hereby granted, free of charge, to any person obtaining a     *
 * copy of this software and associated documentation files (the "Software"),  *
 * to deal in the Software without restriction, including without limitation   *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,    *
 * and/or sell copies of the Software, and to permit persons to whom the       *
 * Software is furnished to do so, subject to the following conditions:        *
 *                                                                             *
 * The above copyright notice and this permission notice shall be included in  *
 * all copies or substantial portions of the Software.                         *
 *                                                                             *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL     *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,     *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR       *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE   *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                      *
 * -------------------------------------------------------------------------- */

#include "openmm/internal/UMAForceImpl.h"
#include "openmm/UMAKernels.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"
#include <set>

using namespace OpenMM;
using std::set;
using std::string;
using std::vector;

UMAForceImpl::UMAForceImpl(const UMAForce& owner) : owner(owner) {
}

UMAForceImpl::~UMAForceImpl() {
}

void UMAForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcUMAForceKernel::Name(), context);
    kernel.getAs<CalcUMAForceKernel>().initialize(context.getSystem(), owner);
}

double UMAForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0)
        return kernel.getAs<CalcUMAForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

vector<string> UMAForceImpl::getKernelNames() {
    vector<string> names;
    names.push_back(CalcUMAForceKernel::Name());
    return names;
}

void UMAForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcUMAForceKernel>().copyParametersToContext(context, owner);
    context.systemChanged();
}

set<int> UMAForceImpl::getAtomIndices() const {
    const vector<int>& subset = owner.getAtomSubset();
    if (subset.empty())
        return set<int>(); // Empty set means all atoms
    return set<int>(subset.begin(), subset.end());
}
