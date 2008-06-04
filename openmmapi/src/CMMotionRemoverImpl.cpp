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

#include "internal/CMMotionRemoverImpl.h"
#include "internal/OpenMMContextImpl.h"
#include "Integrator.h"
#include "System.h"
#include "kernels.h"
#include <vector>

using namespace OpenMM;
using std::vector;

CMMotionRemoverImpl::CMMotionRemoverImpl(CMMotionRemover& owner) : owner(owner) {
}

void CMMotionRemoverImpl::initialize(OpenMMContextImpl& context) {
    kernel = context.getPlatform().createKernel(RemoveCMMotionKernel::Name(), context);
    const System& system = context.getSystem();
    vector<double> masses(system.getNumAtoms());
    for (int i = 0; i < system.getNumAtoms(); ++i)
        masses[i] = system.getAtomMass(i);
    dynamic_cast<RemoveCMMotionKernel&>(kernel.getImpl()).initialize(masses);
}

void CMMotionRemoverImpl::updateContextState(OpenMMContextImpl& context) {
    dynamic_cast<RemoveCMMotionKernel&>(kernel.getImpl()).execute(context.getVelocities());
}

std::vector<std::string> CMMotionRemoverImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(RemoveCMMotionKernel::Name());
    return names;
}
