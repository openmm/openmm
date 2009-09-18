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

#include "OpenCLContext.h"
#include "OpenCLPlatform.h"
#include "OpenCLKernelFactory.h"
#include "OpenCLKernels.h"
#include "openmm/PluginInitializer.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/Context.h"
#include "openmm/System.h"
#include <sstream>

using namespace OpenMM;
using std::map;
using std::string;
using std::stringstream;

extern "C" void initOpenMMPlugin() {
    Platform::registerPlatform(new OpenCLPlatform());
}

OpenCLPlatform::OpenCLPlatform() {
    OpenCLKernelFactory* factory = new OpenCLKernelFactory();
    registerKernelFactory(InitializeForcesKernel::Name(), factory);
    registerKernelFactory(UpdateTimeKernel::Name(), factory);
//    registerKernelFactory(CalcHarmonicBondForceKernel::Name(), factory);
//    registerKernelFactory(CalcHarmonicAngleForceKernel::Name(), factory);
//    registerKernelFactory(CalcPeriodicTorsionForceKernel::Name(), factory);
//    registerKernelFactory(CalcRBTorsionForceKernel::Name(), factory);
//    registerKernelFactory(CalcNonbondedForceKernel::Name(), factory);
//    registerKernelFactory(CalcCustomNonbondedForceKernel::Name(), factory);
//    registerKernelFactory(CalcGBSAOBCForceKernel::Name(), factory);
    registerKernelFactory(IntegrateVerletStepKernel::Name(), factory);
//    registerKernelFactory(IntegrateLangevinStepKernel::Name(), factory);
//    registerKernelFactory(IntegrateBrownianStepKernel::Name(), factory);
//    registerKernelFactory(IntegrateVariableVerletStepKernel::Name(), factory);
//    registerKernelFactory(IntegrateVariableLangevinStepKernel::Name(), factory);
//    registerKernelFactory(ApplyAndersenThermostatKernel::Name(), factory);
    registerKernelFactory(CalcKineticEnergyKernel::Name(), factory);
//    registerKernelFactory(RemoveCMMotionKernel::Name(), factory);
    platformProperties.push_back(OpenCLPlatformIndex());
    platformProperties.push_back(OpenCLDeviceIndex());
    setPropertyDefaultValue(OpenCLPlatformIndex(), "0");
    setPropertyDefaultValue(OpenCLDeviceIndex(), "0");
}

bool OpenCLPlatform::supportsDoublePrecision() const {
    return false;
}

const string& OpenCLPlatform::getPropertyValue(const Context& context, const string& property) const {
    const ContextImpl& impl = getContextImpl(context);
    const PlatformData* data = reinterpret_cast<const PlatformData*>(impl.getPlatformData());
    map<string, string>::const_iterator value = data->propertyValues.find(property);
    if (value != data->propertyValues.end())
        return value->second;
    return Platform::getPropertyValue(context, property);
}

void OpenCLPlatform::setPropertyValue(Context& context, const string& property, const string& value) const {
}

const StreamFactory& OpenCLPlatform::getDefaultStreamFactory() const {
    return defaultStreamFactory;
}

void OpenCLPlatform::contextCreated(ContextImpl& context) const {
    unsigned int platformIndex = 0;
    const string& platformPropValue = getPropertyDefaultValue(OpenCLPlatformIndex());
    if (platformPropValue.length() > 0)
        stringstream(platformPropValue) >> platformIndex;
    unsigned int deviceIndex = 0;
    const string& devicePropValue = getPropertyDefaultValue(OpenCLDeviceIndex());
    if (devicePropValue.length() > 0)
        stringstream(devicePropValue) >> deviceIndex;
    int numParticles = context.getSystem().getNumParticles();
    context.setPlatformData(new PlatformData(numParticles, platformIndex, deviceIndex));
}

void OpenCLPlatform::contextDestroyed(ContextImpl& context) const {
    PlatformData* data = reinterpret_cast<PlatformData*>(context.getPlatformData());
//    gpuShutDown(data->gpu);
    delete data;
}

OpenCLPlatform::PlatformData::PlatformData(int numParticles, int platformIndex, int deviceIndex) : removeCM(false), stepCount(0), computeForceCount(0), time(0.0)  {
    context = new OpenCLContext(numParticles, platformIndex, deviceIndex);
    stringstream platform;
//    device << gpu->platform;
    stringstream device;
//    device << gpu->device;
    propertyValues[OpenCLPlatform::OpenCLPlatformIndex()] = platform.str();
    propertyValues[OpenCLPlatform::OpenCLDeviceIndex()] = device.str();
}
