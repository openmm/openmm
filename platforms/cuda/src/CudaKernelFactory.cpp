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

#include "CudaKernelFactory.h"
#include "CudaKernels.h"
#include "internal/OpenMMContextImpl.h"

using namespace OpenMM;

KernelImpl* CudaKernelFactory::createKernelImpl(std::string name, const Platform& platform, OpenMMContextImpl& context) const {
    _gpuContext* gpu = static_cast<_gpuContext*>(context.getPlatformData());
    if (name == CalcStandardMMForceFieldKernel::Name())
        return new CudaCalcStandardMMForceFieldKernel(name, platform, gpu);
//    if (name == CalcGBSAOBCForceFieldKernel::Name())
//        return new CudaCalcGBSAOBCForceFieldKernel(name, platform);
//    if (name == IntegrateVerletStepKernel::Name())
//        return new CudaIntegrateVerletStepKernel(name, platform);
//    if (name == IntegrateLangevinStepKernel::Name())
//        return new CudaIntegrateLangevinStepKernel(name, platform);
//    if (name == IntegrateBrownianStepKernel::Name())
//        return new CudaIntegrateBrownianStepKernel(name, platform);
//    if (name == ApplyAndersenThermostatKernel::Name())
//        return new CudaApplyAndersenThermostatKernel(name, platform);
//    if (name == CalcKineticEnergyKernel::Name())
//        return new CudaCalcKineticEnergyKernel(name, platform);
//    if (name == RemoveCMMotionKernel::Name())
//        return new CudaRemoveCMMotionKernel(name, platform);
}
