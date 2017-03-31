/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
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

#include "openmm/Context.h"
#include "openmm/Force.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/ForceImpl.h"
#include <vector>

using namespace OpenMM;
using namespace std;

int Force::getForceGroup() const {
    return forceGroup;
}

void Force::setForceGroup(int group) {
    if (group < 0 || group > 31)
        throw OpenMMException("Force group must be between 0 and 31");
    forceGroup = group;
}

bool Force::usesPeriodicBoundaryConditions() const {
    throw OpenMMException("usesPeriodicBoundaryConditions is not implemented");
}

ForceImpl& Force::getImplInContext(Context& context) {
    for (auto impl : context.getImpl().getForceImpls())
        if (&impl->getOwner() == this)
            return *impl;
    throw OpenMMException("getImplInContext: This Force is not present in the Context");
}

const ForceImpl& Force::getImplInContext(const Context& context) const {
    for (auto impl : context.getImpl().getForceImpls())
        if (&impl->getOwner() == this)
            return *impl;
    throw OpenMMException("getImplInContext: This Force is not present in the Context");
}

ContextImpl& Force::getContextImpl(Context& context) {
    return context.getImpl();
}
