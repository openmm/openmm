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
#include "openmm/internal/ContextImpl.h"
#include "openmm/Context.h"
#include "openmm/System.h"
#include <algorithm>
#include <cctype>
#include <sstream>
#ifdef __APPLE__
#include "sys/sysctl.h"
#endif


using namespace OpenMM;
using namespace std;

#ifdef OPENMM_OPENCL_BUILDING_STATIC_LIBRARY
extern "C" void registerOpenCLPlatform() {
    if (OpenCLPlatform::isPlatformSupported())
        Platform::registerPlatform(new OpenCLPlatform());
}
#else
extern "C" OPENMM_EXPORT_OPENCL void registerPlatforms() {
    if (OpenCLPlatform::isPlatformSupported())
        Platform::registerPlatform(new OpenCLPlatform());
}
#endif

OpenCLPlatform::OpenCLPlatform() {
    OpenCLKernelFactory* factory = new OpenCLKernelFactory();
    registerKernelFactory(CalcForcesAndEnergyKernel::Name(), factory);
    registerKernelFactory(UpdateStateDataKernel::Name(), factory);
    registerKernelFactory(ApplyConstraintsKernel::Name(), factory);
    registerKernelFactory(VirtualSitesKernel::Name(), factory);
    registerKernelFactory(CalcHarmonicBondForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomBondForceKernel::Name(), factory);
    registerKernelFactory(CalcHarmonicAngleForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomAngleForceKernel::Name(), factory);
    registerKernelFactory(CalcPeriodicTorsionForceKernel::Name(), factory);
    registerKernelFactory(CalcRBTorsionForceKernel::Name(), factory);
    registerKernelFactory(CalcCMAPTorsionForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomTorsionForceKernel::Name(), factory);
    registerKernelFactory(CalcNonbondedForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomNonbondedForceKernel::Name(), factory);
    registerKernelFactory(CalcGBSAOBCForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomGBForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomExternalForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomHbondForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomCompoundBondForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomManyParticleForceKernel::Name(), factory);
    registerKernelFactory(IntegrateVerletStepKernel::Name(), factory);
    registerKernelFactory(IntegrateLangevinStepKernel::Name(), factory);
    registerKernelFactory(IntegrateBrownianStepKernel::Name(), factory);
    registerKernelFactory(IntegrateVariableVerletStepKernel::Name(), factory);
    registerKernelFactory(IntegrateVariableLangevinStepKernel::Name(), factory);
    registerKernelFactory(IntegrateCustomStepKernel::Name(), factory);
    registerKernelFactory(ApplyAndersenThermostatKernel::Name(), factory);
    registerKernelFactory(ApplyMonteCarloBarostatKernel::Name(), factory);
    registerKernelFactory(RemoveCMMotionKernel::Name(), factory);
    platformProperties.push_back(OpenCLDeviceIndex());
    platformProperties.push_back(OpenCLDeviceName());
    platformProperties.push_back(OpenCLPlatformIndex());
    platformProperties.push_back(OpenCLPlatformName());
    platformProperties.push_back(OpenCLPrecision());
    platformProperties.push_back(OpenCLUseCpuPme());
    setPropertyDefaultValue(OpenCLDeviceIndex(), "");
    setPropertyDefaultValue(OpenCLDeviceName(), "");
    setPropertyDefaultValue(OpenCLPlatformIndex(), "");
    setPropertyDefaultValue(OpenCLPlatformName(), "");
    setPropertyDefaultValue(OpenCLPrecision(), "single");
    setPropertyDefaultValue(OpenCLUseCpuPme(), "false");
}

double OpenCLPlatform::getSpeed() const {
    return 50;
}

bool OpenCLPlatform::supportsDoublePrecision() const {
    return true;
}

bool OpenCLPlatform::isPlatformSupported() {
    // Return false for OpenCL implementations that are known
    // to be buggy (Apple OSX since 10.7.5)

#ifdef __APPLE__
    char str[256];
    size_t size = sizeof(str);
    int ret = sysctlbyname("kern.osrelease", str, &size, NULL, 0);
    if (ret != 0)
        return false;

    int major, minor, micro;
    if (sscanf(str, "%d.%d.%d", &major, &minor, &micro) != 3)
        return false;

    if ((major > 11) || (major == 11 && minor > 4) || (major == 11 && minor == 4 && micro >= 2))
        // 11.4.2 is the darwin release corresponding to OSX 10.7.5, which is the
        // point at which a number of serious bugs were introduced into the
        // Apple OpenCL libraries, resulting in catistrophically incorrect MD simulations
        // (see https://github.com/SimTk/openmm/issues/395 for example). Once a fix is released,
        // this version check should be updated.
        return false;
#endif

    return true;
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

void OpenCLPlatform::contextCreated(ContextImpl& context, const map<string, string>& properties) const {
    const string& platformPropValue = (properties.find(OpenCLPlatformIndex()) == properties.end() ?
            getPropertyDefaultValue(OpenCLPlatformIndex()) : properties.find(OpenCLPlatformIndex())->second);
    const string& devicePropValue = (properties.find(OpenCLDeviceIndex()) == properties.end() ?
            getPropertyDefaultValue(OpenCLDeviceIndex()) : properties.find(OpenCLDeviceIndex())->second);
    string precisionPropValue = (properties.find(OpenCLPrecision()) == properties.end() ?
            getPropertyDefaultValue(OpenCLPrecision()) : properties.find(OpenCLPrecision())->second);
    string cpuPmePropValue = (properties.find(OpenCLUseCpuPme()) == properties.end() ?
            getPropertyDefaultValue(OpenCLUseCpuPme()) : properties.find(OpenCLUseCpuPme())->second);
    transform(precisionPropValue.begin(), precisionPropValue.end(), precisionPropValue.begin(), ::tolower);
    transform(cpuPmePropValue.begin(), cpuPmePropValue.end(), cpuPmePropValue.begin(), ::tolower);
    vector<string> pmeKernelName;
    pmeKernelName.push_back(CalcPmeReciprocalForceKernel::Name());
    if (!supportsKernels(pmeKernelName))
        cpuPmePropValue = "false";
    context.setPlatformData(new PlatformData(context.getSystem(), platformPropValue, devicePropValue, precisionPropValue, cpuPmePropValue));
}

void OpenCLPlatform::contextDestroyed(ContextImpl& context) const {
    PlatformData* data = reinterpret_cast<PlatformData*>(context.getPlatformData());
    delete data;
}

OpenCLPlatform::PlatformData::PlatformData(const System& system, const string& platformPropValue, const string& deviceIndexProperty,
        const string& precisionProperty, const string& cpuPmeProperty) : removeCM(false), stepCount(0), computeForceCount(0), time(0.0)  {
    int platformIndex = -1;
    if (platformPropValue.length() > 0)
        stringstream(platformPropValue) >> platformIndex;
    vector<string> devices;
    size_t searchPos = 0, nextPos;
    while ((nextPos = deviceIndexProperty.find_first_of(", ", searchPos)) != string::npos) {
        devices.push_back(deviceIndexProperty.substr(searchPos, nextPos-searchPos));
        searchPos = nextPos+1;
    }
    devices.push_back(deviceIndexProperty.substr(searchPos));
    try {
        for (int i = 0; i < (int) devices.size(); i++) {
            if (devices[i].length() > 0) {
                unsigned int deviceIndex;
                stringstream(devices[i]) >> deviceIndex;
                contexts.push_back(new OpenCLContext(system, platformIndex, deviceIndex, precisionProperty, *this));
            }
        }
        if (contexts.size() == 0)
            contexts.push_back(new OpenCLContext(system, platformIndex, -1, precisionProperty, *this));
    }
    catch (...) {
        // If an exception was thrown, do our best to clean up memory.
        
        for (int i = 0; i < (int) contexts.size(); i++)
            delete contexts[i];
        throw;
    }
    stringstream deviceIndex, deviceName;
    for (int i = 0; i < (int) contexts.size(); i++) {
        if (i > 0) {
            deviceIndex << ',';
            deviceName << ',';
        }
        deviceIndex << contexts[i]->getDeviceIndex();
        deviceName << contexts[i]->getDevice().getInfo<CL_DEVICE_NAME>();
    }
    platformIndex = contexts[0]->getPlatformIndex();

    useCpuPme = (cpuPmeProperty == "true" && !contexts[0]->getUseDoublePrecision());
    propertyValues[OpenCLPlatform::OpenCLDeviceIndex()] = deviceIndex.str();
    propertyValues[OpenCLPlatform::OpenCLDeviceName()] = deviceName.str();
    propertyValues[OpenCLPlatform::OpenCLPlatformIndex()] = contexts[0]->intToString(platformIndex);
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);
    propertyValues[OpenCLPlatform::OpenCLPlatformName()] = platforms[platformIndex].getInfo<CL_PLATFORM_NAME>();
    propertyValues[OpenCLPlatform::OpenCLPrecision()] = precisionProperty;
    propertyValues[OpenCLPlatform::OpenCLUseCpuPme()] = useCpuPme ? "true" : "false";
    contextEnergy.resize(contexts.size());
}

OpenCLPlatform::PlatformData::~PlatformData() {
    for (int i = 0; i < (int) contexts.size(); i++)
        delete contexts[i];
}

void OpenCLPlatform::PlatformData::initializeContexts(const System& system) {
    for (int i = 0; i < (int) contexts.size(); i++)
        contexts[i]->initialize();
}

void OpenCLPlatform::PlatformData::syncContexts() {
    for (int i = 0; i < (int) contexts.size(); i++)
        contexts[i]->getWorkThread().flush();
}
