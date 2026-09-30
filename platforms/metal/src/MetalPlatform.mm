/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM OpenCL Platform.                                    *
 * Source: platforms/opencl/src/OpenCLPlatform.cpp                         *
 *                                                                            *
 * Original OpenCL Platform code:                                             *
 * Portions copyright (c) 2008-2026 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 *                                                                            *
 * Metal Platform code:                                                       *
 * Portions copyright (c) 2026 Chun-Chi Hung.                                 *
 * Authors: Chun-Chi Hung                                                     *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the               *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.       *
 * -------------------------------------------------------------------------- */

#include "MetalPlatform.h"
#include "MetalContext.h"
#include "MetalKernelFactory.h"
#include "MetalKernels.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#import <Metal/Metal.h>
#include <algorithm>
#include <cctype>

using namespace OpenMM;
using namespace std;

#ifdef OPENMM_COMMON_BUILDING_STATIC_LIBRARY
extern "C" void registerMetalPlatform() {
#else
extern "C" OPENMM_EXPORT_COMMON void registerPlatforms() {
#endif
    if (MetalPlatform::isPlatformSupported())
        Platform::registerPlatform(new MetalPlatform());
}

namespace {

/** @brief Normalize textual enum properties without signed-char ctype ambiguity. */
string lowercase(string value) {
    transform(value.begin(), value.end(), value.begin(), [](unsigned char c) { return static_cast<char>(tolower(c)); });
    return value;
}

/** @brief Return the concrete properties of a supported native device. */
map<string, string> deviceProperties(id<MTLDevice> device) {
    return {{MetalPlatform::MetalDeviceIndex(), "0"},
            {MetalPlatform::MetalDeviceName(), string(device.name.UTF8String)},
            {MetalPlatform::MetalPrecision(), "single"},
            {MetalPlatform::MetalUseCpuPme(), "false"}};
}

} // namespace

MetalPlatform::MetalPlatform(bool floatingAccumulators) : floatingAccumulators(floatingAccumulators) {
    MetalKernelFactory* factory = new MetalKernelFactory();
    registerKernelFactory(CalcForcesAndEnergyKernel::Name(), factory);
    registerKernelFactory(UpdateStateDataKernel::Name(), factory);
    registerKernelFactory(ApplyConstraintsKernel::Name(), factory);
    registerKernelFactory(VirtualSitesKernel::Name(), factory);
    registerKernelFactory(MinimizeKernel::Name(), factory);
    registerKernelFactory(CalcHarmonicBondForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomBondForceKernel::Name(), factory);
    registerKernelFactory(CalcHarmonicAngleForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomAngleForceKernel::Name(), factory);
    registerKernelFactory(CalcPeriodicTorsionForceKernel::Name(), factory);
    registerKernelFactory(CalcRBTorsionForceKernel::Name(), factory);
    registerKernelFactory(CalcCMAPTorsionForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomTorsionForceKernel::Name(), factory);
    registerKernelFactory(CalcNonbondedForceKernel::Name(), factory);
    registerKernelFactory(CalcConstantPotentialForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomNonbondedForceKernel::Name(), factory);
    registerKernelFactory(CalcGBSAOBCForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomGBForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomExternalForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomHbondForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomCentroidBondForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomCompoundBondForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomCVForceKernel::Name(), factory);
    registerKernelFactory(CalcATMForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomCPPForceKernel::Name(), factory);
    registerKernelFactory(CalcOrientationRestraintForceKernel::Name(), factory);
    registerKernelFactory(CalcPythonForceKernel::Name(), factory);
    registerKernelFactory(CalcRGForceKernel::Name(), factory);
    registerKernelFactory(CalcRMSDForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomManyParticleForceKernel::Name(), factory);
    registerKernelFactory(CalcGayBerneForceKernel::Name(), factory);
    registerKernelFactory(CalcLCPOForceKernel::Name(), factory);
    registerKernelFactory(IntegrateVerletStepKernel::Name(), factory);
    registerKernelFactory(IntegrateLangevinMiddleStepKernel::Name(), factory);
    registerKernelFactory(IntegrateBrownianStepKernel::Name(), factory);
    registerKernelFactory(IntegrateVariableVerletStepKernel::Name(), factory);
    registerKernelFactory(IntegrateVariableLangevinStepKernel::Name(), factory);
    registerKernelFactory(IntegrateCustomStepKernel::Name(), factory);
    registerKernelFactory(IntegrateDPDStepKernel::Name(), factory);
    registerKernelFactory(IntegrateQTBStepKernel::Name(), factory);
    registerKernelFactory(ApplyAndersenThermostatKernel::Name(), factory);
    registerKernelFactory(IntegrateNoseHooverStepKernel::Name(), factory);
    registerKernelFactory(ApplyMonteCarloBarostatKernel::Name(), factory);
    registerKernelFactory(RemoveCMMotionKernel::Name(), factory);
    platformProperties.push_back(MetalDeviceIndex());
    platformProperties.push_back(MetalDeviceName());
    platformProperties.push_back(MetalPrecision());
    platformProperties.push_back(MetalUseCpuPme());
    setPropertyDefaultValue(MetalDeviceIndex(), "");
    setPropertyDefaultValue(MetalDeviceName(), "");
    setPropertyDefaultValue(MetalPrecision(), "single");
    setPropertyDefaultValue(MetalUseCpuPme(), "false");
}

const string& MetalPlatform::getName() const {
    static const string name = "Metal";
    return name;
}

double MetalPlatform::getSpeed() const {
    return 50;
}

bool MetalPlatform::supportsDoublePrecision() const {
    return false;
}

bool MetalPlatform::isPlatformSupported() {
    @autoreleasepool {
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        return device != nil && [device supportsFamily:MTLGPUFamilyApple7];
    }
}

const string& MetalPlatform::MetalDeviceIndex() {
    static const string property = "DeviceIndex";
    return property;
}

const string& MetalPlatform::MetalDeviceName() {
    static const string property = "DeviceName";
    return property;
}

const string& MetalPlatform::MetalPrecision() {
    static const string property = "Precision";
    return property;
}

const string& MetalPlatform::MetalUseCpuPme() {
    static const string property = "UseCpuPme";
    return property;
}

const string& MetalPlatform::getPropertyValue(const Context& context, const string& property) const {
    const PlatformData& data = *static_cast<const PlatformData*>(getContextImpl(context).getPlatformData());
    auto value = data.propertyValues.find(property);
    if (value != data.propertyValues.end())
        return value->second;
    return Platform::getPropertyValue(context, property);
}

void MetalPlatform::setPropertyValue(Context& context, const string& property, const string& value) const {
    if (getPropertyValue(context, property) != value)
        throw OpenMMException("The Metal Platform property '"+property+"' cannot be changed after Context creation");
}

vector<map<string, string> > MetalPlatform::getDevices(const map<string, string>& filters) const {
    @autoreleasepool {
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        if (device == nil || ![device supportsFamily:MTLGPUFamilyApple7])
            return {};
        map<string, string> properties = deviceProperties(device);
        for (const auto& filter : filters) {
            auto value = properties.find(filter.first);
            if (value == properties.end())
                return {};
            string requested = filter.second;
            if (filter.first == MetalPrecision() || filter.first == MetalUseCpuPme())
                requested = lowercase(requested);
            // Empty device selectors choose the default; other values are exact filters.
            if ((filter.first == MetalDeviceIndex() || filter.first == MetalDeviceName()) && requested.empty())
                continue;
            if (requested != value->second)
                return {};
        }
        return {properties};
    }
}

void MetalPlatform::contextCreated(ContextImpl& context, const map<string, string>& properties) const {
    map<string, string> selected;
    for (const string& property : platformProperties)
        selected[property] = getPropertyDefaultValue(property);
    for (const auto& property : properties) {
        if (selected.find(property.first) == selected.end())
            throw OpenMMException("Unknown Metal Platform property '"+property.first+"'");
        selected[property.first] = property.second;
    }
    selected[MetalPrecision()] = lowercase(selected[MetalPrecision()]);
    selected[MetalUseCpuPme()] = lowercase(selected[MetalUseCpuPme()]);
    if (selected[MetalPrecision()] != "single")
        throw OpenMMException("The Metal Platform supports only single precision");
    if (selected[MetalUseCpuPme()] != "false")
        throw OpenMMException("The Metal Platform does not support CPU PME or CPU fallback");
    if (!selected[MetalDeviceIndex()].empty() && selected[MetalDeviceIndex()] != "0")
        throw OpenMMException("The Metal Platform supports only the default GPU (DeviceIndex 0)");
    if (getDevices(selected).empty())
        throw OpenMMException("No supported Apple silicon GPU matches the requested Metal properties");
    context.setPlatformData(new PlatformData(context, nullptr, floatingAccumulators));
}

void MetalPlatform::linkedContextCreated(ContextImpl& context, ContextImpl& originalContext) const {
    const PlatformData& original = *static_cast<const PlatformData*>(originalContext.getPlatformData());
    context.setPlatformData(new PlatformData(context, original.computeContext.get(),
            original.computeContext->getUseFloatingPointAccumulators()));
}

void MetalPlatform::contextDestroyed(ContextImpl& context) const {
    delete static_cast<PlatformData*>(context.getPlatformData());
}

MetalPlatform::PlatformData::PlatformData(ContextImpl& context, MetalContext* linked, bool floatingAccumulators) : context(&context),
        computeContext(new MetalContext(context.getSystem(), &context, linked, floatingAccumulators)) {
    @autoreleasepool {
        id<MTLDevice> device = (__bridge id<MTLDevice>) computeContext->getDevice();
        propertyValues = deviceProperties(device);
    }
}

MetalPlatform::PlatformData::~PlatformData() {
}
