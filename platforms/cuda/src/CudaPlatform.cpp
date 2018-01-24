/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2017 Stanford University and the Authors.      *
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
#include "openmm/Context.h"
#include "openmm/System.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/hardware.h"
#include <algorithm>
#include <cctype>
#include <sstream>
#include <cstdio>
#ifdef _MSC_VER
    #include <Windows.h>
#endif
using namespace OpenMM;
using namespace std;

#define CHECK_RESULT(result, prefix) \
    if (result != CUDA_SUCCESS) { \
        std::stringstream m; \
        m<<prefix<<": "<<CudaContext::getErrorString(result)<<" ("<<result<<")"<<" at "<<__FILE__<<":"<<__LINE__; \
        throw OpenMMException(m.str());\
    }


#ifdef OPENMM_CUDA_BUILDING_STATIC_LIBRARY
extern "C" void registerCudaPlatform() {
    Platform::registerPlatform(new CudaPlatform());
}
#else
extern "C" OPENMM_EXPORT_CUDA void registerPlatforms() {
    Platform::registerPlatform(new CudaPlatform());
}
#endif

CudaPlatform::CudaPlatform() {
    deprecatedPropertyReplacements["CudaDeviceIndex"] = CudaDeviceIndex();
    deprecatedPropertyReplacements["CudaDeviceName"] = CudaDeviceName();
    deprecatedPropertyReplacements["CudaUseBlockingSync"] = CudaUseBlockingSync();
    deprecatedPropertyReplacements["CudaPrecision"] = CudaPrecision();
    deprecatedPropertyReplacements["CudaUseCpuPme"] = CudaUseCpuPme();
    deprecatedPropertyReplacements["CudaTempDirectory"] = CudaTempDirectory();
    deprecatedPropertyReplacements["CudaDisablePmeStream"] = CudaDisablePmeStream();
    deprecatedPropertyReplacements["CudaDeterministicForces"] = CudaDeterministicForces();
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
    registerKernelFactory(CalcCustomCentroidBondForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomCompoundBondForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomCVForceKernel::Name(), factory);
    registerKernelFactory(CalcRMSDForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomManyParticleForceKernel::Name(), factory);
    registerKernelFactory(CalcGayBerneForceKernel::Name(), factory);
    registerKernelFactory(IntegrateVerletStepKernel::Name(), factory);
    registerKernelFactory(IntegrateLangevinStepKernel::Name(), factory);
    registerKernelFactory(IntegrateBrownianStepKernel::Name(), factory);
    registerKernelFactory(IntegrateVariableVerletStepKernel::Name(), factory);
    registerKernelFactory(IntegrateVariableLangevinStepKernel::Name(), factory);
    registerKernelFactory(IntegrateCustomStepKernel::Name(), factory);
    registerKernelFactory(ApplyAndersenThermostatKernel::Name(), factory);
    registerKernelFactory(ApplyMonteCarloBarostatKernel::Name(), factory);
    registerKernelFactory(RemoveCMMotionKernel::Name(), factory);
    platformProperties.push_back(CudaDeviceIndex());
    platformProperties.push_back(CudaDeviceName());
    platformProperties.push_back(CudaUseBlockingSync());
    platformProperties.push_back(CudaPrecision());
    platformProperties.push_back(CudaUseCpuPme());
    platformProperties.push_back(CudaCompiler());
    platformProperties.push_back(CudaTempDirectory());
    platformProperties.push_back(CudaHostCompiler());
    platformProperties.push_back(CudaDisablePmeStream());
    platformProperties.push_back(CudaDeterministicForces());
    setPropertyDefaultValue(CudaDeviceIndex(), "");
    setPropertyDefaultValue(CudaDeviceName(), "");
    setPropertyDefaultValue(CudaUseBlockingSync(), "true");
    setPropertyDefaultValue(CudaPrecision(), "single");
    setPropertyDefaultValue(CudaUseCpuPme(), "false");
    setPropertyDefaultValue(CudaDisablePmeStream(), "false");
    setPropertyDefaultValue(CudaDeterministicForces(), "false");
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
    char* compiler = getenv("OPENMM_CUDA_COMPILER");
    string nvcc = (compiler == NULL ? "/usr/local/cuda/bin/nvcc" : string(compiler));
    setPropertyDefaultValue(CudaCompiler(), nvcc);
    char* tmpdir = getenv("TMPDIR");
    string tmp = (tmpdir == NULL ? string(P_tmpdir) : string(tmpdir));
    setPropertyDefaultValue(CudaTempDirectory(), tmp);
#endif
    char* hostCompiler = getenv("CUDA_HOST_COMPILER");
    setPropertyDefaultValue(CudaHostCompiler(), (hostCompiler == NULL ? "" : string(hostCompiler)));
}

double CudaPlatform::getSpeed() const {
    return 100;
}

bool CudaPlatform::supportsDoublePrecision() const {
    return true;
}

const string& CudaPlatform::getPropertyValue(const Context& context, const string& property) const {
    const ContextImpl& impl = getContextImpl(context);
    const PlatformData* data = reinterpret_cast<const PlatformData*>(impl.getPlatformData());
    string propertyName = property;
    if (deprecatedPropertyReplacements.find(property) != deprecatedPropertyReplacements.end())
        propertyName = deprecatedPropertyReplacements.find(property)->second;
    map<string, string>::const_iterator value = data->propertyValues.find(propertyName);
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
    string cpuPmePropValue = (properties.find(CudaUseCpuPme()) == properties.end() ?
            getPropertyDefaultValue(CudaUseCpuPme()) : properties.find(CudaUseCpuPme())->second);
    const string& compilerPropValue = (properties.find(CudaCompiler()) == properties.end() ?
            getPropertyDefaultValue(CudaCompiler()) : properties.find(CudaCompiler())->second);
    const string& tempPropValue = (properties.find(CudaTempDirectory()) == properties.end() ?
            getPropertyDefaultValue(CudaTempDirectory()) : properties.find(CudaTempDirectory())->second);
    const string& hostCompilerPropValue = (properties.find(CudaHostCompiler()) == properties.end() ?
            getPropertyDefaultValue(CudaHostCompiler()) : properties.find(CudaHostCompiler())->second);
    string pmeStreamPropValue = (properties.find(CudaDisablePmeStream()) == properties.end() ?
            getPropertyDefaultValue(CudaDisablePmeStream()) : properties.find(CudaDisablePmeStream())->second);
    string deterministicForcesValue = (properties.find(CudaDeterministicForces()) == properties.end() ?
            getPropertyDefaultValue(CudaDeterministicForces()) : properties.find(CudaDeterministicForces())->second);
    transform(blockingPropValue.begin(), blockingPropValue.end(), blockingPropValue.begin(), ::tolower);
    transform(precisionPropValue.begin(), precisionPropValue.end(), precisionPropValue.begin(), ::tolower);
    transform(cpuPmePropValue.begin(), cpuPmePropValue.end(), cpuPmePropValue.begin(), ::tolower);
    transform(pmeStreamPropValue.begin(), pmeStreamPropValue.end(), pmeStreamPropValue.begin(), ::tolower);
    transform(deterministicForcesValue.begin(), deterministicForcesValue.end(), deterministicForcesValue.begin(), ::tolower);
    vector<string> pmeKernelName;
    pmeKernelName.push_back(CalcPmeReciprocalForceKernel::Name());
    if (!supportsKernels(pmeKernelName))
        cpuPmePropValue = "false";
    int threads = getNumProcessors();
    char* threadsEnv = getenv("OPENMM_CPU_THREADS");
    if (threadsEnv != NULL)
        stringstream(threadsEnv) >> threads;
    context.setPlatformData(new PlatformData(&context, context.getSystem(), devicePropValue, blockingPropValue, precisionPropValue, cpuPmePropValue, compilerPropValue, tempPropValue,
            hostCompilerPropValue, pmeStreamPropValue, deterministicForcesValue, threads, NULL));
}

void CudaPlatform::linkedContextCreated(ContextImpl& context, ContextImpl& originalContext) const {
    Platform& platform = originalContext.getPlatform();
    string devicePropValue = platform.getPropertyValue(originalContext.getOwner(), CudaDeviceIndex());
    string blockingPropValue = platform.getPropertyValue(originalContext.getOwner(), CudaUseBlockingSync());
    string precisionPropValue = platform.getPropertyValue(originalContext.getOwner(), CudaPrecision());
    string cpuPmePropValue = platform.getPropertyValue(originalContext.getOwner(), CudaUseCpuPme());
    string compilerPropValue = platform.getPropertyValue(originalContext.getOwner(), CudaCompiler());
    string tempPropValue = platform.getPropertyValue(originalContext.getOwner(), CudaTempDirectory());
    string hostCompilerPropValue = platform.getPropertyValue(originalContext.getOwner(), CudaHostCompiler());
    string pmeStreamPropValue = platform.getPropertyValue(originalContext.getOwner(), CudaDisablePmeStream());
    string deterministicForcesValue = platform.getPropertyValue(originalContext.getOwner(), CudaDeterministicForces());
    int threads = reinterpret_cast<PlatformData*>(originalContext.getPlatformData())->threads.getNumThreads();
    context.setPlatformData(new PlatformData(&context, context.getSystem(), devicePropValue, blockingPropValue, precisionPropValue, cpuPmePropValue, compilerPropValue, tempPropValue,
            hostCompilerPropValue, pmeStreamPropValue, deterministicForcesValue, threads, &originalContext));
}

void CudaPlatform::contextDestroyed(ContextImpl& context) const {
    PlatformData* data = reinterpret_cast<PlatformData*>(context.getPlatformData());
    delete data;
}

CudaPlatform::PlatformData::PlatformData(ContextImpl* context, const System& system, const string& deviceIndexProperty, const string& blockingProperty, const string& precisionProperty,
            const string& cpuPmeProperty, const string& compilerProperty, const string& tempProperty, const string& hostCompilerProperty, const string& pmeStreamProperty,
            const string& deterministicForcesProperty, int numThreads, ContextImpl* originalContext) :
                context(context), removeCM(false), stepCount(0), computeForceCount(0), time(0.0), hasInitializedContexts(false), threads(numThreads) {
    bool blocking = (blockingProperty == "true");
    vector<string> devices;
    size_t searchPos = 0, nextPos;
    while ((nextPos = deviceIndexProperty.find_first_of(", ", searchPos)) != string::npos) {
        devices.push_back(deviceIndexProperty.substr(searchPos, nextPos-searchPos));
        searchPos = nextPos+1;
    }
    devices.push_back(deviceIndexProperty.substr(searchPos));
    PlatformData* originalData = NULL;
    if (originalContext != NULL)
        originalData = reinterpret_cast<PlatformData*>(originalContext->getPlatformData());
    try {
        for (int i = 0; i < (int) devices.size(); i++) {
            if (devices[i].length() > 0) {
                int deviceIndex;
                stringstream(devices[i]) >> deviceIndex;
                contexts.push_back(new CudaContext(system, deviceIndex, blocking, precisionProperty, compilerProperty, tempProperty, hostCompilerProperty, *this, (originalData == NULL ? NULL : originalData->contexts[i])));
            }
        }
        if (contexts.size() == 0)
            contexts.push_back(new CudaContext(system, -1, blocking, precisionProperty, compilerProperty, tempProperty, hostCompilerProperty, *this, (originalData == NULL ? NULL : originalData->contexts[0])));
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
        char name[1000];
        CHECK_RESULT(cuDeviceGetName(name, 1000, contexts[i]->getDevice()), "Error querying device name");
        deviceName << name;
    }
    size_t printfsize;
    cuCtxGetLimit(&printfsize, CU_LIMIT_PRINTF_FIFO_SIZE);
    cuCtxSetLimit(CU_LIMIT_PRINTF_FIFO_SIZE, 10*printfsize);

    useCpuPme = (cpuPmeProperty == "true" && !contexts[0]->getUseDoublePrecision());
    disablePmeStream = (pmeStreamProperty == "true");
    deterministicForces = (deterministicForcesProperty == "true");
    propertyValues[CudaPlatform::CudaDeviceIndex()] = deviceIndex.str();
    propertyValues[CudaPlatform::CudaDeviceName()] = deviceName.str();
    propertyValues[CudaPlatform::CudaUseBlockingSync()] = blocking ? "true" : "false";
    propertyValues[CudaPlatform::CudaPrecision()] = precisionProperty;
    propertyValues[CudaPlatform::CudaUseCpuPme()] = useCpuPme ? "true" : "false";
    propertyValues[CudaPlatform::CudaCompiler()] = compilerProperty;
    propertyValues[CudaPlatform::CudaTempDirectory()] = tempProperty;
    propertyValues[CudaPlatform::CudaHostCompiler()] = hostCompilerProperty;
    propertyValues[CudaPlatform::CudaDisablePmeStream()] = disablePmeStream ? "true" : "false";
    propertyValues[CudaPlatform::CudaDeterministicForces()] = deterministicForces ? "true" : "false";
    contextEnergy.resize(contexts.size());
    
    // Determine whether peer-to-peer copying is supported, and enable it if so.
    
    peerAccessSupported = true;
    for (int i = 1; i < contexts.size(); i++) {
        int canAccess;
        cuDeviceCanAccessPeer(&canAccess, contexts[i]->getDevice(), contexts[0]->getDevice());
        if (!canAccess) {
            peerAccessSupported = false;
            break;
        }
    }
}

CudaPlatform::PlatformData::~PlatformData() {
    for (int i = 0; i < (int) contexts.size(); i++)
        delete contexts[i];
}

void CudaPlatform::PlatformData::initializeContexts(const System& system) {
    if (hasInitializedContexts)
        return;
    for (int i = 0; i < (int) contexts.size(); i++)
        contexts[i]->initialize();
    hasInitializedContexts = true;
}

void CudaPlatform::PlatformData::syncContexts() {
    for (int i = 0; i < (int) contexts.size(); i++)
        contexts[i]->getWorkThread().flush();
}
