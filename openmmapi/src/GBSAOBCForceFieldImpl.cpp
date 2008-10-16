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

#include "internal/GBSAOBCForceFieldImpl.h"
#include "internal/OpenMMContextImpl.h"
#include "kernels.h"
#include <vector>

using namespace OpenMM;
using std::vector;

GBSAOBCForceFieldImpl::GBSAOBCForceFieldImpl(GBSAOBCForceField& owner) : owner(owner) {
}

void GBSAOBCForceFieldImpl::initialize(OpenMMContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcGBSAOBCForceFieldKernel::Name(), context);
    dynamic_cast<CalcGBSAOBCForceFieldKernel&>(kernel.getImpl()).initialize(context.getSystem(), owner);
}

void GBSAOBCForceFieldImpl::calcForces(OpenMMContextImpl& context, Stream& forces) {
    dynamic_cast<CalcGBSAOBCForceFieldKernel&>(kernel.getImpl()).executeForces(context);
}

double GBSAOBCForceFieldImpl::calcEnergy(OpenMMContextImpl& context) {
    return dynamic_cast<CalcGBSAOBCForceFieldKernel&>(kernel.getImpl()).executeEnergy(context);
}

std::vector<std::string> GBSAOBCForceFieldImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcGBSAOBCForceFieldKernel::Name());
    return names;
}
