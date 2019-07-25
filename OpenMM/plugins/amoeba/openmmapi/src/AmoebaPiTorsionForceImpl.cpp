/* -------------------------------------------------------------------------- *
 *                               OpenMMAmoeba                                 *
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

#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/AmoebaPiTorsionForceImpl.h"
#include "openmm/amoebaKernels.h"

using namespace OpenMM;

using std::pair;
using std::vector;
using std::set;

AmoebaPiTorsionForceImpl::AmoebaPiTorsionForceImpl(const AmoebaPiTorsionForce& owner) : owner(owner) {
}

AmoebaPiTorsionForceImpl::~AmoebaPiTorsionForceImpl() {
}

void AmoebaPiTorsionForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcAmoebaPiTorsionForceKernel::Name(), context);
    kernel.getAs<CalcAmoebaPiTorsionForceKernel>().initialize(context.getSystem(), owner);
}

double AmoebaPiTorsionForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0)
        return kernel.getAs<CalcAmoebaPiTorsionForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

std::vector<std::string> AmoebaPiTorsionForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcAmoebaPiTorsionForceKernel::Name());
    return names;
}

void AmoebaPiTorsionForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcAmoebaPiTorsionForceKernel>().copyParametersToContext(context, owner);
    context.systemChanged();
}
