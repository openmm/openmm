/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2016 Stanford University and the Authors.      *
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

#include "CpuKernelFactory.h"
#include "CpuKernels.h"
#include "CpuPlatform.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

KernelImpl* CpuKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    CpuPlatform::PlatformData& data = CpuPlatform::getPlatformData(context);
    if (name == CalcForcesAndEnergyKernel::Name())
        return new CpuCalcForcesAndEnergyKernel(name, platform, data, context);
    if (name == CalcHarmonicAngleForceKernel::Name())
        return new CpuCalcHarmonicAngleForceKernel(name, platform, data);
    if (name == CalcPeriodicTorsionForceKernel::Name())
        return new CpuCalcPeriodicTorsionForceKernel(name, platform, data);
    if (name == CalcRBTorsionForceKernel::Name())
        return new CpuCalcRBTorsionForceKernel(name, platform, data);
    if (name == CalcNonbondedForceKernel::Name())
        return new CpuCalcNonbondedForceKernel(name, platform, data);
    if (name == CalcCustomNonbondedForceKernel::Name())
        return new CpuCalcCustomNonbondedForceKernel(name, platform, data);
    if (name == CalcCustomManyParticleForceKernel::Name())
        return new CpuCalcCustomManyParticleForceKernel(name, platform, data);
    if (name == CalcGBSAOBCForceKernel::Name())
        return new CpuCalcGBSAOBCForceKernel(name, platform, data);
    if (name == CalcCustomGBForceKernel::Name())
        return new CpuCalcCustomGBForceKernel(name, platform, data);
    if (name == CalcGayBerneForceKernel::Name())
        return new CpuCalcGayBerneForceKernel(name, platform, data);
    if (name == IntegrateLangevinStepKernel::Name())
        return new CpuIntegrateLangevinStepKernel(name, platform, data);
    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '") + name + "'").c_str());
}
