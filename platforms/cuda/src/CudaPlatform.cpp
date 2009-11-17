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
#include "openmm/internal/ContextImpl.h"
#include "kernels/gputypes.h"
#include "openmm/Context.h"
#include "openmm/System.h"
#include <sstream>

using namespace OpenMM;
using std::map;
using std::string;
using std::stringstream;

extern "C" void initOpenMMPlugin() {
    if (gpuIsAvailable())
        Platform::registerPlatform(new CudaPlatform());
}

CudaPlatform::CudaPlatform() {
    CudaKernelFactory* factory = new CudaKernelFactory();
    registerKernelFactory(CalcForcesAndEnergyKernel::Name(), factory);
    registerKernelFactory(UpdateStateDataKernel::Name(), factory);
    registerKernelFactory(CalcHarmonicBondForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomBondForceKernel::Name(), factory);
    registerKernelFactory(CalcHarmonicAngleForceKernel::Name(), factory);
    registerKernelFactory(CalcPeriodicTorsionForceKernel::Name(), factory);
    registerKernelFactory(CalcRBTorsionForceKernel::Name(), factory);
    registerKernelFactory(CalcNonbondedForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomNonbondedForceKernel::Name(), factory);
    registerKernelFactory(CalcGBSAOBCForceKernel::Name(), factory);
    registerKernelFactory(CalcGBVIForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomExternalForceKernel::Name(), factory);
    registerKernelFactory(IntegrateVerletStepKernel::Name(), factory);
    registerKernelFactory(IntegrateLangevinStepKernel::Name(), factory);
    registerKernelFactory(IntegrateBrownianStepKernel::Name(), factory);
    registerKernelFactory(IntegrateVariableVerletStepKernel::Name(), factory);
    registerKernelFactory(IntegrateVariableLangevinStepKernel::Name(), factory);
    registerKernelFactory(ApplyAndersenThermostatKernel::Name(), factory);
    registerKernelFactory(CalcKineticEnergyKernel::Name(), factory);
    registerKernelFactory(RemoveCMMotionKernel::Name(), factory);
    platformProperties.push_back(CudaDevice());
    platformProperties.push_back(CudaUseBlockingSync());
    setPropertyDefaultValue(CudaDevice(), "0");
    setPropertyDefaultValue(CudaUseBlockingSync(), "true");
}

bool CudaPlatform::supportsDoublePrecision() const {
    return false;
}

const string& CudaPlatform::getPropertyValue(const Context& context, const string& property) const {
    const ContextImpl& impl = getContextImpl(context);
    const PlatformData* data = reinterpret_cast<const PlatformData*>(impl.getPlatformData());
    map<string, string>::const_iterator value = data->propertyValues.find(property);
    if (value != data->propertyValues.end())
        return value->second;
    return Platform::getPropertyValue(context, property);
}

void CudaPlatform::setPropertyValue(Context& context, const string& property, const string& value) const {
}

void CudaPlatform::contextCreated(ContextImpl& context) const {
    unsigned int device = 0;
    const string& devicePropValue = getPropertyDefaultValue(CudaDevice());
    if (devicePropValue.length() > 0)
        stringstream(devicePropValue) >> device;
    int numParticles = context.getSystem().getNumParticles();
    _gpuContext* gpu = (_gpuContext*) gpuInit(numParticles, device, getPropertyDefaultValue(CudaUseBlockingSync()) == "true");
    context.setPlatformData(new PlatformData(gpu));
}

void CudaPlatform::contextDestroyed(ContextImpl& context) const {
    PlatformData* data = reinterpret_cast<PlatformData*>(context.getPlatformData());
    gpuShutDown(data->gpu);
    delete data;
}

CudaPlatform::PlatformData::PlatformData(_gpuContext* gpu) : gpu(gpu), removeCM(false), nonbondedMethod(0), customNonbondedMethod(0), hasBonds(false), hasAngles(false),
        hasPeriodicTorsions(false), hasRB(false), hasNonbonded(false), hasCustomNonbonded(false), stepCount(0), computeForceCount(0), time(0.0),
        ewaldSelfEnergy(0.0) {
    stringstream device;
    device << gpu->device;
    propertyValues[CudaPlatform::CudaDevice()] = device.str();
    propertyValues[CudaPlatform::CudaUseBlockingSync()] = (gpu->useBlockingSync ? "true" : "false");
}
