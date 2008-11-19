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

#include "ReferenceKernelFactory.h"
#include "ReferenceKernels.h"
#include "OpenMMException.h"

using namespace OpenMM;

KernelImpl* ReferenceKernelFactory::createKernelImpl(std::string name, const Platform& platform, OpenMMContextImpl& context) const {
    if (name == InitializeForcesKernel::Name())
        return new ReferenceInitializeForcesKernel(name, platform);
    else if (name == CalcNonbondedForceKernel::Name())
        return new ReferenceCalcNonbondedForceKernel(name, platform);
    else if (name == CalcHarmonicBondForceKernel::Name())
        return new ReferenceCalcHarmonicBondForceKernel(name, platform);
    else if (name == CalcHarmonicAngleForceKernel::Name())
        return new ReferenceCalcHarmonicAngleForceKernel(name, platform);
    else if (name == CalcHarmonicAngleForceKernel::Name())
        return new ReferenceCalcHarmonicAngleForceKernel(name, platform);
    else if (name == CalcPeriodicTorsionForceKernel::Name())
        return new ReferenceCalcPeriodicTorsionForceKernel(name, platform);
    else if (name == CalcRBTorsionForceKernel::Name())
        return new ReferenceCalcRBTorsionForceKernel(name, platform);
    else if (name == CalcGBSAOBCForceKernel::Name())
        return new ReferenceCalcGBSAOBCForceKernel(name, platform);
    else if (name == IntegrateVerletStepKernel::Name())
        return new ReferenceIntegrateVerletStepKernel(name, platform);
    else if (name == IntegrateLangevinStepKernel::Name())
        return new ReferenceIntegrateLangevinStepKernel(name, platform);
    else if (name == IntegrateBrownianStepKernel::Name())
        return new ReferenceIntegrateBrownianStepKernel(name, platform);
    else if (name == ApplyAndersenThermostatKernel::Name())
        return new ReferenceApplyAndersenThermostatKernel(name, platform);
    else if (name == CalcKineticEnergyKernel::Name())
        return new ReferenceCalcKineticEnergyKernel(name, platform);
    else if (name == RemoveCMMotionKernel::Name())
        return new ReferenceRemoveCMMotionKernel(name, platform);
    else {
        throw OpenMMException( (std::string("Tried to create kernel with illegal kernel name '") + name + "'").c_str() );
    }
}
