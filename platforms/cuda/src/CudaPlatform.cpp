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
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "CudaPlatform.h"
#include "CudaKernelFactory.h"
#include "CudaKernels.h"
#include "openmm/PluginInitializer.h"
#include "openmm/internal/OpenMMContextImpl.h"
#include "kernels/gputypes.h"
#include "openmm/System.h"

using namespace OpenMM;

extern "C" void initOpenMMPlugin() {
    if (gpuIsAvailable())
        Platform::registerPlatform(new CudaPlatform());
}

CudaPlatform::CudaPlatform() {
    CudaKernelFactory* factory = new CudaKernelFactory();
    registerKernelFactory(InitializeForcesKernel::Name(), factory);
    registerKernelFactory(UpdateTimeKernel::Name(), factory);
    registerKernelFactory(CalcHarmonicBondForceKernel::Name(), factory);
    registerKernelFactory(CalcHarmonicAngleForceKernel::Name(), factory);
    registerKernelFactory(CalcPeriodicTorsionForceKernel::Name(), factory);
    registerKernelFactory(CalcRBTorsionForceKernel::Name(), factory);
    registerKernelFactory(CalcNonbondedForceKernel::Name(), factory);
    registerKernelFactory(CalcGBSAOBCForceKernel::Name(), factory);
    registerKernelFactory(IntegrateVerletStepKernel::Name(), factory);
    registerKernelFactory(IntegrateLangevinStepKernel::Name(), factory);
    registerKernelFactory(IntegrateBrownianStepKernel::Name(), factory);
    registerKernelFactory(IntegrateVariableVerletStepKernel::Name(), factory);
    registerKernelFactory(ApplyAndersenThermostatKernel::Name(), factory);
    registerKernelFactory(CalcKineticEnergyKernel::Name(), factory);
    registerKernelFactory(RemoveCMMotionKernel::Name(), factory);
}

bool CudaPlatform::supportsDoublePrecision() const {
    return false;
}

const StreamFactory& CudaPlatform::getDefaultStreamFactory() const {
    return defaultStreamFactory;
}

void CudaPlatform::contextCreated(OpenMMContextImpl& context) const {
    int numParticles = context.getSystem().getNumParticles();
    _gpuContext* gpu = (_gpuContext*) gpuInit(numParticles);
    context.setPlatformData(new PlatformData(gpu));
}

void CudaPlatform::contextDestroyed(OpenMMContextImpl& context) const {
    PlatformData* data = reinterpret_cast<PlatformData*>(context.getPlatformData());
    gpuShutDown(data->gpu);
    delete data;
}
