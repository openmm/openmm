/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
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

#include "CudaContext.h"
#include "CudaExpressionUtilities.h"
#include "CudaPlatform.h"
#include "CudaKernelFactory.h"
#include "CudaKernels.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/Context.h"
#include "openmm/System.h"
#include <algorithm>
#include <cctype>
#include <sstream>
#include <cstdio>
#ifdef _MSC_VER
    #include <Windows.h>
#endif
using namespace OpenMM;
using namespace std;

extern "C" OPENMM_EXPORT void registerPlatforms() {
    Platform::registerPlatform(new CudaPlatform());
}

CudaPlatform::CudaPlatform() {
    CudaKernelFactory* factory = new CudaKernelFactory();
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
    registerKernelFactory(IntegrateVerletStepKernel::Name(), factory);
    registerKernelFactory(IntegrateLangevinStepKernel::Name(), factory);
    registerKernelFactory(IntegrateBrownianStepKernel::Name(), factory);
    registerKernelFactory(IntegrateVariableVerletStepKernel::Name(), factory);
    registerKernelFactory(IntegrateVariableLangevinStepKernel::Name(), factory);
    registerKernelFactory(IntegrateCustomStepKernel::Name(), factory);
    registerKernelFactory(ApplyAndersenThermostatKernel::Name(), factory);
    registerKernelFactory(ApplyMonteCarloBarostatKernel::Name(), factory);
    registerKernelFactory(CalcKineticEnergyKernel::Name(), factory);
    registerKernelFactory(RemoveCMMotionKernel::Name(), factory);
    platformProperties.push_back(CudaDeviceIndex());
    platformProperties.push_back(CudaUseBlockingSync());
    platformProperties.push_back(CudaPrecision());
    platformProperties.push_back(CudaCompiler());
    platformProperties.push_back(CudaTempDirectory());
    setPropertyDefaultValue(CudaDeviceIndex(), "");
    setPropertyDefaultValue(CudaUseBlockingSync(), "true");
    setPropertyDefaultValue(CudaPrecision(), "single");
#ifdef _MSC_VER
    char* bindir = getenv("CUDA_BIN_PATH");
	string nvcc = (bindir == NULL ? "nvcc.exe" : string(bindir)+"\\nvcc.exe");
    int length = GetShortPathName(nvcc.c_str(), NULL, 0);
	if (length > 0) {
		vector<char> shortName(length);
	    GetShortPathName(nvcc.c_str(), &shortName[0], length);
		nvcc = string(&shortName[0]);
	}
    setPropertyDefaultValue(CudaCompiler(), nvcc);
    setPropertyDefaultValue(CudaTempDirectory(), string(getenv("TEMP")));
#else
    setPropertyDefaultValue(CudaCompiler(), "/usr/local/cuda/bin/nvcc");
    char* tmpdir = getenv("TMPDIR");
    if (tmpdir == NULL)
        tmpdir = P_tmpdir;
    setPropertyDefaultValue(CudaTempDirectory(), string(tmpdir));
#endif
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

void CudaPlatform::contextCreated(ContextImpl& context, const map<string, string>& properties) const {
    const string& devicePropValue = (properties.find(CudaDeviceIndex()) == properties.end() ?
            getPropertyDefaultValue(CudaDeviceIndex()) : properties.find(CudaDeviceIndex())->second);
    string blockingPropValue = (properties.find(CudaUseBlockingSync()) == properties.end() ?
            getPropertyDefaultValue(CudaUseBlockingSync()) : properties.find(CudaUseBlockingSync())->second);
    string precisionPropValue = (properties.find(CudaPrecision()) == properties.end() ?
            getPropertyDefaultValue(CudaPrecision()) : properties.find(CudaPrecision())->second);
    const string& compilerPropValue = (properties.find(CudaCompiler()) == properties.end() ?
            getPropertyDefaultValue(CudaCompiler()) : properties.find(CudaCompiler())->second);
    const string& tempPropValue = (properties.find(CudaTempDirectory()) == properties.end() ?
            getPropertyDefaultValue(CudaTempDirectory()) : properties.find(CudaTempDirectory())->second);
    transform(blockingPropValue.begin(), blockingPropValue.end(), blockingPropValue.begin(), ::tolower);
    transform(precisionPropValue.begin(), precisionPropValue.end(), precisionPropValue.begin(), ::tolower);
    context.setPlatformData(new PlatformData(context.getSystem(), devicePropValue, blockingPropValue, precisionPropValue, compilerPropValue, tempPropValue));
}

void CudaPlatform::contextDestroyed(ContextImpl& context) const {
    PlatformData* data = reinterpret_cast<PlatformData*>(context.getPlatformData());
    delete data;
}

CudaPlatform::PlatformData::PlatformData(const System& system, const string& deviceIndexProperty, const string& blockingProperty, const string& precisionProperty,
            const string& compilerProperty, const string& tempProperty) : removeCM(false), stepCount(0), computeForceCount(0), time(0.0)  {
    bool blocking = (blockingProperty == "true");
    vector<string> devices;
    size_t searchPos = 0, nextPos;
    while ((nextPos = deviceIndexProperty.find_first_of(", ", searchPos)) != string::npos) {
        devices.push_back(deviceIndexProperty.substr(searchPos, nextPos-searchPos));
        searchPos = nextPos+1;
    }
    devices.push_back(deviceIndexProperty.substr(searchPos));
    for (int i = 0; i < (int) devices.size(); i++) {
        if (devices[i].length() > 0) {
            unsigned int deviceIndex;
            stringstream(devices[i]) >> deviceIndex;
            contexts.push_back(new CudaContext(system, deviceIndex, blocking, precisionProperty, compilerProperty, tempProperty, *this));
        }
    }
    if (contexts.size() == 0)
        contexts.push_back(new CudaContext(system, -1, blocking, precisionProperty, compilerProperty, tempProperty, *this));
    stringstream device;
    for (int i = 0; i < (int) contexts.size(); i++) {
        if (i > 0)
            device << ',';
        device << contexts[i]->getDeviceIndex();
    }
    propertyValues[CudaPlatform::CudaDeviceIndex()] = device.str();
    propertyValues[CudaPlatform::CudaUseBlockingSync()] = blocking ? "true" : "false";
    propertyValues[CudaPlatform::CudaPrecision()] = precisionProperty;
    propertyValues[CudaPlatform::CudaCompiler()] = compilerProperty;
    propertyValues[CudaPlatform::CudaTempDirectory()] = tempProperty;
    contextEnergy.resize(contexts.size());
}

CudaPlatform::PlatformData::~PlatformData() {
    for (int i = 0; i < (int) contexts.size(); i++)
        delete contexts[i];
}

void CudaPlatform::PlatformData::initializeContexts(const System& system) {
    for (int i = 0; i < (int) contexts.size(); i++)
        contexts[i]->initialize();
}

void CudaPlatform::PlatformData::syncContexts() {
    for (int i = 0; i < (int) contexts.size(); i++)
        contexts[i]->getWorkThread().flush();
}
