/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2016 Stanford University and the Authors.      *
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

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include <cmath>
#include "CudaContext.h"
#include "CudaArray.h"
#include "CudaBondedUtilities.h"
#include "CudaForceInfo.h"
#include "CudaIntegrationUtilities.h"
#include "CudaKernels.h"
#include "CudaKernelSources.h"
#include "CudaNonbondedUtilities.h"
#include "SHA1.h"
#include "hilbert.h"
#include "openmm/OpenMMException.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/VirtualSite.h"
#include "CudaExpressionUtilities.h"
#include "openmm/internal/ContextImpl.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <typeinfo>
#include <sys/stat.h>
#include <cudaProfiler.h>
#ifndef WIN32
  #include <unistd.h>
#endif


#define CHECK_RESULT(result) CHECK_RESULT2(result, errorMessage);
#define CHECK_RESULT2(result, prefix) \
    if (result != CUDA_SUCCESS) { \
        std::stringstream m; \
        m<<prefix<<": "<<getErrorString(result)<<" ("<<result<<")"<<" at "<<__FILE__<<":"<<__LINE__; \
        throw OpenMMException(m.str());\
    }

using namespace OpenMM;
using namespace std;

const int CudaContext::ThreadBlockSize = 64;
const int CudaContext::TileSize = sizeof(tileflags)*8;
bool CudaContext::hasInitializedCuda = false;

#ifdef WIN32
#include <Windows.h>
static int executeInWindows(const string &command) {
    // COMSPEC is an env variable pointing to full dir of cmd.exe
    // it always defined on pretty much all Windows OSes
    string fullcommand = getenv("COMSPEC") + string(" /C ") + command;
    STARTUPINFO si;
    PROCESS_INFORMATION pi;
    ZeroMemory( &si, sizeof(si) );
    si.cb = sizeof(si);
    ZeroMemory( &pi, sizeof(pi) );
    vector<char> args(std::max(1000, (int) fullcommand.size()+1));
    strcpy(&args[0], fullcommand.c_str());
    si.dwFlags = STARTF_USESHOWWINDOW;
    si.wShowWindow = SW_HIDE;
    if (!CreateProcess(NULL, &args[0], NULL, NULL, FALSE, 0, NULL, NULL, &si, &pi)) {
        return -1;
    }
    WaitForSingleObject(pi.hProcess, INFINITE);
    DWORD exitCode = -1;
    if(!GetExitCodeProcess(pi.hProcess, &exitCode)) {
        throw(OpenMMException("Could not get nvcc.exe's exit code\n"));
    } else {
        if(exitCode == 0)
            return 0;
        else
            return -1;
    }
}
#endif

CudaContext::CudaContext(const System& system, int deviceIndex, bool useBlockingSync, const string& precision, const string& compiler,
        const string& tempDir, const std::string& hostCompiler, CudaPlatform::PlatformData& platformData) : system(system), currentStream(0),
        time(0.0), platformData(platformData), stepCount(0), computeForceCount(0), stepsSinceReorder(99999), contextIsValid(false), atomsWereReordered(false), hasCompilerKernel(false), isNvccAvailable(false),
        pinnedBuffer(NULL), posq(NULL), posqCorrection(NULL), velm(NULL), force(NULL), energyBuffer(NULL), energyParamDerivBuffer(NULL), atomIndexDevice(NULL), integration(NULL), expression(NULL), bonded(NULL), nonbonded(NULL), thread(NULL) {
    // Determine what compiler to use.
    
    this->compiler = "\""+compiler+"\"";
    if (platformData.context != NULL) {
        try {
            compilerKernel = platformData.context->getPlatform().createKernel(CudaCompilerKernel::Name(), *platformData.context);
            hasCompilerKernel = true;
        }
        catch (...) {
            // The runtime compiler plugin isn't available.
        }
    }
#ifdef WIN32
    string testCompilerCommand = this->compiler+" --version > nul 2> nul";
    int res = executeInWindows(testCompilerCommand.c_str());
#else
    string testCompilerCommand = this->compiler+" --version > /dev/null 2> /dev/null";
    int res = std::system(testCompilerCommand.c_str());
#endif
    struct stat info;
    isNvccAvailable = (res == 0 && stat(tempDir.c_str(), &info) == 0);
    int cudaDriverVersion;
    cuDriverGetVersion(&cudaDriverVersion);
    static bool hasShownNvccWarning = false;
    if (hasCompilerKernel && !isNvccAvailable && !hasShownNvccWarning && cudaDriverVersion < 8000) {
        hasShownNvccWarning = true;
        printf("Could not find nvcc.  Using runtime compiler, which may produce slower performance.  ");
#ifdef WIN32
        printf("Set CUDA_BIN_PATH to specify where nvcc is located.\n");
#else
        printf("Set OPENMM_CUDA_COMPILER to specify where nvcc is located.\n");
#endif
    }
    if (hostCompiler.size() > 0)
        this->compiler = compiler+" --compiler-bindir "+hostCompiler;
    if (!hasInitializedCuda) {
        CHECK_RESULT2(cuInit(0), "Error initializing CUDA");
        hasInitializedCuda = true;
    }
    if (precision == "single") {
        useDoublePrecision = false;
        useMixedPrecision = false;
    }
    else if (precision == "mixed") {
        useDoublePrecision = false;
        useMixedPrecision = true;
    }
    else if (precision == "double") {
        useDoublePrecision = true;
        useMixedPrecision = false;
    }
    else
        throw OpenMMException("Illegal value for Precision: "+precision);
    char* cacheVariable = getenv("OPENMM_CACHE_DIR");
    cacheDir = (cacheVariable == NULL ? tempDir : string(cacheVariable));
#ifdef WIN32
    this->tempDir = tempDir+"\\";
    cacheDir = cacheDir+"\\";
#else
    this->tempDir = tempDir+"/";
    cacheDir = cacheDir+"/";
#endif
    contextIndex = platformData.contexts.size();
    int numDevices;
    string errorMessage = "Error initializing Context";
    CHECK_RESULT(cuDeviceGetCount(&numDevices));
    if (deviceIndex < -1 || deviceIndex >= numDevices)
        throw OpenMMException("Illegal value for DeviceIndex: "+intToString(deviceIndex));

    vector<int> devicePrecedence;
    if (deviceIndex == -1) {
        devicePrecedence = getDevicePrecedence();
    } else {
        devicePrecedence.push_back(deviceIndex);
    }

    this->deviceIndex = -1;
    for (int i = 0; i < static_cast<int>(devicePrecedence.size()); i++) {
        int trialDeviceIndex = devicePrecedence[i];
        CHECK_RESULT(cuDeviceGet(&device, trialDeviceIndex));
        defaultOptimizationOptions = "--use_fast_math";
        unsigned int flags = CU_CTX_MAP_HOST;
        if (useBlockingSync)
            flags += CU_CTX_SCHED_BLOCKING_SYNC;
        else
            flags += CU_CTX_SCHED_SPIN;

        if (cuCtxCreate(&context, flags, device) == CUDA_SUCCESS) {
            this->deviceIndex = trialDeviceIndex;
            break;
        }
    }
    if (this->deviceIndex == -1)
        if (deviceIndex != -1)
            throw OpenMMException("The requested CUDA device could not be loaded");
        else
            throw OpenMMException("No compatible CUDA device is available");

    int major, minor;
    CHECK_RESULT(cuDeviceComputeCapability(&major, &minor, device));
    int numThreadBlocksPerComputeUnit = (major >= 6 ? 4 : 6);
    if (cudaDriverVersion < 7000) {
        // This is a workaround to support GTX 980 with CUDA 6.5.  It reports
        // its compute capability as 5.2, but the compiler doesn't support
        // anything beyond 5.0.
        if (major == 5)
            minor = 0;
    }
    if (cudaDriverVersion < 8000) {
        // This is a workaround to support Pascal with CUDA 7.5.  It reports
        // its compute capability as 6.x, but the compiler doesn't support
        // anything beyond 5.3.
        if (major == 6) {
            major = 5;
            minor = 3;
        }
    }
    gpuArchitecture = intToString(major)+intToString(minor);
    computeCapability = major+0.1*minor;

    contextIsValid = true;
    CHECK_RESULT(cuCtxSetCacheConfig(CU_FUNC_CACHE_PREFER_SHARED));
    if (contextIndex > 0) {
        int canAccess;
        cuDeviceCanAccessPeer(&canAccess, getDevice(), platformData.contexts[0]->getDevice());
        if (canAccess) {
            platformData.contexts[0]->setAsCurrent();
            CHECK_RESULT(cuCtxEnablePeerAccess(getContext(), 0));
            setAsCurrent();
            CHECK_RESULT(cuCtxEnablePeerAccess(platformData.contexts[0]->getContext(), 0));
        }
    }
    numAtoms = system.getNumParticles();
    paddedNumAtoms = TileSize*((numAtoms+TileSize-1)/TileSize);
    numAtomBlocks = (paddedNumAtoms+(TileSize-1))/TileSize;
    int multiprocessors;
    CHECK_RESULT(cuDeviceGetAttribute(&multiprocessors, CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT, device));
    numThreadBlocks = numThreadBlocksPerComputeUnit*multiprocessors;
    if (useDoublePrecision) {
        posq = CudaArray::create<double4>(*this, paddedNumAtoms, "posq");
        velm = CudaArray::create<double4>(*this, paddedNumAtoms, "velm");
        compilationDefines["USE_DOUBLE_PRECISION"] = "1";
        compilationDefines["make_real2"] = "make_double2";
        compilationDefines["make_real3"] = "make_double3";
        compilationDefines["make_real4"] = "make_double4";
        compilationDefines["make_mixed2"] = "make_double2";
        compilationDefines["make_mixed3"] = "make_double3";
        compilationDefines["make_mixed4"] = "make_double4";
    }
    else if (useMixedPrecision) {
        posq = CudaArray::create<float4>(*this, paddedNumAtoms, "posq");
        posqCorrection = CudaArray::create<float4>(*this, paddedNumAtoms, "posqCorrection");
        velm = CudaArray::create<double4>(*this, paddedNumAtoms, "velm");
        compilationDefines["USE_MIXED_PRECISION"] = "1";
        compilationDefines["make_real2"] = "make_float2";
        compilationDefines["make_real3"] = "make_float3";
        compilationDefines["make_real4"] = "make_float4";
        compilationDefines["make_mixed2"] = "make_double2";
        compilationDefines["make_mixed3"] = "make_double3";
        compilationDefines["make_mixed4"] = "make_double4";
    }
    else {
        posq = CudaArray::create<float4>(*this, paddedNumAtoms, "posq");
        velm = CudaArray::create<float4>(*this, paddedNumAtoms, "velm");
        compilationDefines["make_real2"] = "make_float2";
        compilationDefines["make_real3"] = "make_float3";
        compilationDefines["make_real4"] = "make_float4";
        compilationDefines["make_mixed2"] = "make_float2";
        compilationDefines["make_mixed3"] = "make_float3";
        compilationDefines["make_mixed4"] = "make_float4";
    }
    posCellOffsets.resize(paddedNumAtoms, make_int4(0, 0, 0, 0));

    // Create utility kernels that are used in multiple places.

    CUmodule utilities = createModule(CudaKernelSources::vectorOps+CudaKernelSources::utilities);
    clearBufferKernel = getKernel(utilities, "clearBuffer");
    clearTwoBuffersKernel = getKernel(utilities, "clearTwoBuffers");
    clearThreeBuffersKernel = getKernel(utilities, "clearThreeBuffers");
    clearFourBuffersKernel = getKernel(utilities, "clearFourBuffers");
    clearFiveBuffersKernel = getKernel(utilities, "clearFiveBuffers");
    clearSixBuffersKernel = getKernel(utilities, "clearSixBuffers");

    // Set defines based on the requested precision.

    compilationDefines["SQRT"] = useDoublePrecision ? "sqrt" : "sqrtf";
    compilationDefines["RSQRT"] = useDoublePrecision ? "rsqrt" : "rsqrtf";
    compilationDefines["RECIP"] = useDoublePrecision ? "1.0/" : "1.0f/";
    compilationDefines["EXP"] = useDoublePrecision ? "exp" : "expf";
    compilationDefines["LOG"] = useDoublePrecision ? "log" : "logf";
    compilationDefines["POW"] = useDoublePrecision ? "pow" : "powf";
    compilationDefines["COS"] = useDoublePrecision ? "cos" : "cosf";
    compilationDefines["SIN"] = useDoublePrecision ? "sin" : "sinf";
    compilationDefines["TAN"] = useDoublePrecision ? "tan" : "tanf";
    compilationDefines["ACOS"] = useDoublePrecision ? "acos" : "acosf";
    compilationDefines["ASIN"] = useDoublePrecision ? "asin" : "asinf";
    compilationDefines["ATAN"] = useDoublePrecision ? "atan" : "atanf";
    compilationDefines["ERF"] = useDoublePrecision ? "erf" : "erff";
    compilationDefines["ERFC"] = useDoublePrecision ? "erfc" : "erfcf";

    // Set defines for applying periodic boundary conditions.

    Vec3 boxVectors[3];
    system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    boxIsTriclinic = (boxVectors[0][1] != 0.0 || boxVectors[0][2] != 0.0 ||
                      boxVectors[1][0] != 0.0 || boxVectors[1][2] != 0.0 ||
                      boxVectors[2][0] != 0.0 || boxVectors[2][1] != 0.0);
    if (boxIsTriclinic) {
        compilationDefines["APPLY_PERIODIC_TO_DELTA(delta)"] =
            "{"
            "real scale3 = floor(delta.z*invPeriodicBoxSize.z+0.5f); \\\n"
            "delta.x -= scale3*periodicBoxVecZ.x; \\\n"
            "delta.y -= scale3*periodicBoxVecZ.y; \\\n"
            "delta.z -= scale3*periodicBoxVecZ.z; \\\n"
            "real scale2 = floor(delta.y*invPeriodicBoxSize.y+0.5f); \\\n"
            "delta.x -= scale2*periodicBoxVecY.x; \\\n"
            "delta.y -= scale2*periodicBoxVecY.y; \\\n"
            "real scale1 = floor(delta.x*invPeriodicBoxSize.x+0.5f); \\\n"
            "delta.x -= scale1*periodicBoxVecX.x;}";
        compilationDefines["APPLY_PERIODIC_TO_POS(pos)"] =
            "{"
            "real scale3 = floor(pos.z*invPeriodicBoxSize.z); \\\n"
            "pos.x -= scale3*periodicBoxVecZ.x; \\\n"
            "pos.y -= scale3*periodicBoxVecZ.y; \\\n"
            "pos.z -= scale3*periodicBoxVecZ.z; \\\n"
            "real scale2 = floor(pos.y*invPeriodicBoxSize.y); \\\n"
            "pos.x -= scale2*periodicBoxVecY.x; \\\n"
            "pos.y -= scale2*periodicBoxVecY.y; \\\n"
            "real scale1 = floor(pos.x*invPeriodicBoxSize.x); \\\n"
            "pos.x -= scale1*periodicBoxVecX.x;}";
        compilationDefines["APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, center)"] =
            "{"
            "real scale3 = floor((pos.z-center.z)*invPeriodicBoxSize.z+0.5f); \\\n"
            "pos.x -= scale3*periodicBoxVecZ.x; \\\n"
            "pos.y -= scale3*periodicBoxVecZ.y; \\\n"
            "pos.z -= scale3*periodicBoxVecZ.z; \\\n"
            "real scale2 = floor((pos.y-center.y)*invPeriodicBoxSize.y+0.5f); \\\n"
            "pos.x -= scale2*periodicBoxVecY.x; \\\n"
            "pos.y -= scale2*periodicBoxVecY.y; \\\n"
            "real scale1 = floor((pos.x-center.x)*invPeriodicBoxSize.x+0.5f); \\\n"
            "pos.x -= scale1*periodicBoxVecX.x;}";
    }
    else {
        compilationDefines["APPLY_PERIODIC_TO_DELTA(delta)"] =
            "{"
            "delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x; \\\n"
            "delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y; \\\n"
            "delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;}";
        compilationDefines["APPLY_PERIODIC_TO_POS(pos)"] =
            "{"
            "pos.x -= floor(pos.x*invPeriodicBoxSize.x)*periodicBoxSize.x; \\\n"
            "pos.y -= floor(pos.y*invPeriodicBoxSize.y)*periodicBoxSize.y; \\\n"
            "pos.z -= floor(pos.z*invPeriodicBoxSize.z)*periodicBoxSize.z;}";
        compilationDefines["APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, center)"] =
            "{"
            "pos.x -= floor((pos.x-center.x)*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x; \\\n"
            "pos.y -= floor((pos.y-center.y)*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y; \\\n"
            "pos.z -= floor((pos.z-center.z)*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;}";
    }

    // Create the work thread used for parallelization when running on multiple devices.

    thread = new WorkThread();

    // Create utilities objects.

    bonded = new CudaBondedUtilities(*this);
    nonbonded = new CudaNonbondedUtilities(*this);
    integration = new CudaIntegrationUtilities(*this, system);
    expression = new CudaExpressionUtilities(*this);
}

CudaContext::~CudaContext() {
    setAsCurrent();
    for (int i = 0; i < (int) forces.size(); i++)
        delete forces[i];
    for (int i = 0; i < (int) reorderListeners.size(); i++)
        delete reorderListeners[i];
    for (int i = 0; i < (int) preComputations.size(); i++)
        delete preComputations[i];
    for (int i = 0; i < (int) postComputations.size(); i++)
        delete postComputations[i];
    if (pinnedBuffer != NULL)
        cuMemFreeHost(pinnedBuffer);
    if (posq != NULL)
        delete posq;
    if (posqCorrection != NULL)
        delete posqCorrection;
    if (velm != NULL)
        delete velm;
    if (force != NULL)
        delete force;
    if (energyBuffer != NULL)
        delete energyBuffer;
    if (energyParamDerivBuffer != NULL)
        delete energyParamDerivBuffer;
    if (atomIndexDevice != NULL)
        delete atomIndexDevice;
    if (integration != NULL)
        delete integration;
    if (expression != NULL)
        delete expression;
    if (bonded != NULL)
        delete bonded;
    if (nonbonded != NULL)
        delete nonbonded;
    if (thread != NULL)
        delete thread;
    string errorMessage = "Error deleting Context";
    if (contextIsValid) {
        cuProfilerStop();
        CHECK_RESULT(cuCtxDestroy(context));
    }
    contextIsValid = false;
}

void CudaContext::initialize() {
    cuCtxSetCurrent(context);
    string errorMessage = "Error initializing Context";
    int numEnergyBuffers = max(numThreadBlocks*ThreadBlockSize, nonbonded->getNumEnergyBuffers());
    if (useDoublePrecision) {
        energyBuffer = CudaArray::create<double>(*this, numEnergyBuffers, "energyBuffer");
        int pinnedBufferSize = max(paddedNumAtoms*4, numEnergyBuffers);
        CHECK_RESULT(cuMemHostAlloc(&pinnedBuffer, pinnedBufferSize*sizeof(double), 0));
    }
    else if (useMixedPrecision) {
        energyBuffer = CudaArray::create<double>(*this, numEnergyBuffers, "energyBuffer");
        int pinnedBufferSize = max(paddedNumAtoms*4, numEnergyBuffers);
        CHECK_RESULT(cuMemHostAlloc(&pinnedBuffer, pinnedBufferSize*sizeof(double), 0));
    }
    else {
        energyBuffer = CudaArray::create<float>(*this, numEnergyBuffers, "energyBuffer");
        int pinnedBufferSize = max(paddedNumAtoms*6, numEnergyBuffers);
        CHECK_RESULT(cuMemHostAlloc(&pinnedBuffer, pinnedBufferSize*sizeof(float), 0));
    }
    for (int i = 0; i < numAtoms; i++) {
        double mass = system.getParticleMass(i);
        if (useDoublePrecision || useMixedPrecision)
            ((double4*) pinnedBuffer)[i] = make_double4(0.0, 0.0, 0.0, mass == 0.0 ? 0.0 : 1.0/mass);
        else
            ((float4*) pinnedBuffer)[i] = make_float4(0.0f, 0.0f, 0.0f, mass == 0.0 ? 0.0f : (float) (1.0/mass));
    }
    velm->upload(pinnedBuffer);
    bonded->initialize(system);
    force = CudaArray::create<long long>(*this, paddedNumAtoms*3, "force");
    addAutoclearBuffer(force->getDevicePointer(), force->getSize()*force->getElementSize());
    addAutoclearBuffer(energyBuffer->getDevicePointer(), energyBuffer->getSize()*energyBuffer->getElementSize());
    int numEnergyParamDerivs = energyParamDerivNames.size();
    if (numEnergyParamDerivs > 0) {
        if (useDoublePrecision || useMixedPrecision)
            energyParamDerivBuffer = CudaArray::create<double>(*this, numEnergyParamDerivs*numEnergyBuffers, "energyParamDerivBuffer");
        else
            energyParamDerivBuffer = CudaArray::create<float>(*this, numEnergyParamDerivs*numEnergyBuffers, "energyParamDerivBuffer");
        addAutoclearBuffer(*energyParamDerivBuffer);
    }
    atomIndexDevice = CudaArray::create<int>(*this, paddedNumAtoms, "atomIndex");
    atomIndex.resize(paddedNumAtoms);
    for (int i = 0; i < paddedNumAtoms; ++i)
        atomIndex[i] = i;
    atomIndexDevice->upload(atomIndex);
    findMoleculeGroups();
    nonbonded->initialize(system);
}

void CudaContext::addForce(CudaForceInfo* force) {
    forces.push_back(force);
}

void CudaContext::setAsCurrent() {
    if (contextIsValid)
        cuCtxSetCurrent(context);
}

string CudaContext::replaceStrings(const string& input, const std::map<std::string, std::string>& replacements) const {
    static set<char> symbolChars;
    if (symbolChars.size() == 0) {
        symbolChars.insert('_');
        for (char c = 'a'; c <= 'z'; c++)
            symbolChars.insert(c);
        for (char c = 'A'; c <= 'Z'; c++)
            symbolChars.insert(c);
        for (char c = '0'; c <= '9'; c++)
            symbolChars.insert(c);
    }
    string result = input;
    for (map<string, string>::const_iterator iter = replacements.begin(); iter != replacements.end(); iter++) {
        int index = 0;
        int size = iter->first.size();
        do {
            index = result.find(iter->first, index);
            if (index != result.npos) {
                if ((index == 0 || symbolChars.find(result[index-1]) == symbolChars.end()) && (index == result.size()-size || symbolChars.find(result[index+size]) == symbolChars.end())) {
                    // We have found a complete symbol, not part of a longer symbol.

                    result.replace(index, size, iter->second);
                    index += iter->second.size();
                }
                else
                    index++;
            }
        } while (index != result.npos);
    }
    return result;
}

CUmodule CudaContext::createModule(const string source, const char* optimizationFlags) {
    return createModule(source, map<string, string>(), optimizationFlags);
}

CUmodule CudaContext::createModule(const string source, const map<string, string>& defines, const char* optimizationFlags) {
    string bits = intToString(8*sizeof(void*));
    string options = (optimizationFlags == NULL ? defaultOptimizationOptions : string(optimizationFlags));
    stringstream src;
    if (!options.empty())
        src << "// Compilation Options: " << options << endl << endl;
    for (map<string, string>::const_iterator iter = compilationDefines.begin(); iter != compilationDefines.end(); ++iter) {
        src << "#define " << iter->first;
        if (!iter->second.empty())
            src << " " << iter->second;
        src << endl;
    }
    if (!compilationDefines.empty())
        src << endl;
    if (useDoublePrecision) {
        src << "typedef double real;\n";
        src << "typedef double2 real2;\n";
        src << "typedef double3 real3;\n";
        src << "typedef double4 real4;\n";
    }
    else {
        src << "typedef float real;\n";
        src << "typedef float2 real2;\n";
        src << "typedef float3 real3;\n";
        src << "typedef float4 real4;\n";
    }
    if (useDoublePrecision || useMixedPrecision) {
        src << "typedef double mixed;\n";
        src << "typedef double2 mixed2;\n";
        src << "typedef double3 mixed3;\n";
        src << "typedef double4 mixed4;\n";
    }
    else {
        src << "typedef float mixed;\n";
        src << "typedef float2 mixed2;\n";
        src << "typedef float3 mixed3;\n";
        src << "typedef float4 mixed4;\n";
    }
    src << "typedef unsigned int tileflags;\n";
    for (map<string, string>::const_iterator iter = defines.begin(); iter != defines.end(); ++iter) {
        src << "#define " << iter->first;
        if (!iter->second.empty())
            src << " " << iter->second;
        src << endl;
    }
    if (!defines.empty())
        src << endl;
    src << source << endl;

    // See whether we already have PTX for this kernel cached.

    CSHA1 sha1;
    sha1.Update((const UINT_8*) src.str().c_str(), src.str().size());
    sha1.Final();
    UINT_8 hash[20];
    sha1.GetHash(hash);
    stringstream cacheFile;
    cacheFile << cacheDir;
    cacheFile.flags(ios::hex);
    for (int i = 0; i < 20; i++)
        cacheFile << setw(2) << setfill('0') << (int) hash[i];
    cacheFile << '_' << gpuArchitecture << '_' << bits;
    CUmodule module;
    if (cuModuleLoad(&module, cacheFile.str().c_str()) == CUDA_SUCCESS)
        return module;

    // Select names for the various temporary files.

    stringstream tempFileName;
    tempFileName << "openmmTempKernel" << this; // Include a pointer to this context as part of the filename to avoid collisions.
#ifdef WIN32
    tempFileName << "_" << GetCurrentProcessId();
#else
    tempFileName << "_" << getpid();
#endif
    string inputFile = (tempDir+tempFileName.str()+".cu");
    string outputFile = (tempDir+tempFileName.str()+".ptx");
    string logFile = (tempDir+tempFileName.str()+".log");
    int res = 0;

    // If the runtime compiler plugin is available, use it.

    if (hasCompilerKernel && !isNvccAvailable) {
        string ptx = compilerKernel.getAs<CudaCompilerKernel>().createModule(src.str(), "-arch=compute_"+gpuArchitecture+" "+options, *this);

        // If possible, write the PTX out to a temporary file so we can cache it for later use.

        bool wroteCache = false;
        try {
            ofstream out(outputFile.c_str());
            out << ptx;
            out.close();
            if (!out.fail())
                wroteCache = true;
        }
        catch (...) {
            // Ignore.
        }
        if (!wroteCache) {
            // An error occurred.  Possibly we don't have permission to write to the temp directory.  Just try to load the module directly.

            CHECK_RESULT2(cuModuleLoadDataEx(&module, &ptx[0], 0, NULL, NULL), "Error loading CUDA module");
            return module;
        }
    }
    else {
        // Write out the source to a temporary file.

        ofstream out(inputFile.c_str());
        out << src.str();
        out.close();
#ifdef WIN32
#ifdef _DEBUG
        string command = compiler+" --ptx -G -g --machine "+bits+" -arch=sm_"+gpuArchitecture+" -o "+outputFile+" "+options+" "+inputFile+" 2> "+logFile;
#else
        string command = compiler+" --ptx -lineinfo --machine "+bits+" -arch=sm_"+gpuArchitecture+" -o "+outputFile+" "+options+" "+inputFile+" 2> "+logFile;
#endif
        int res = executeInWindows(command);
#else
        string command = compiler+" --ptx --machine "+bits+" -arch=sm_"+gpuArchitecture+" -o \""+outputFile+"\" "+options+" \""+inputFile+"\" 2> \""+logFile+"\"";
        res = std::system(command.c_str());
#endif
    }
    try {
        if (res != 0) {
            // Load the error log.

            stringstream error;
            error << "Error launching CUDA compiler: " << res;
            ifstream log(logFile.c_str());
            if (log.is_open()) {
                string line;
                while (!log.eof()) {
                    getline(log, line);
                    error << '\n' << line;
                }
                log.close();
            }
            throw OpenMMException(error.str());
        }
        CUresult result = cuModuleLoad(&module, outputFile.c_str());
        if (result != CUDA_SUCCESS) {
            std::stringstream m;
            m<<"Error loading CUDA module: "<<getErrorString(result)<<" ("<<result<<")";
            throw OpenMMException(m.str());
        }
        remove(inputFile.c_str());
        if (rename(outputFile.c_str(), cacheFile.str().c_str()) != 0)
            remove(outputFile.c_str());
        remove(logFile.c_str());
        return module;
    }
    catch (...) {
        remove(inputFile.c_str());
        remove(outputFile.c_str());
        remove(logFile.c_str());
        throw;
    }
}

CUfunction CudaContext::getKernel(CUmodule& module, const string& name) {
    CUfunction function;
    CUresult result = cuModuleGetFunction(&function, module, name.c_str());
    if (result != CUDA_SUCCESS) {
        std::stringstream m;
        m<<"Error creating kernel "<<name<<": "<<getErrorString(result)<<" ("<<result<<")";
        throw OpenMMException(m.str());
    }
    return function;
}

CUstream CudaContext::getCurrentStream() {
    return currentStream;
}

void CudaContext::setCurrentStream(CUstream stream) {
    currentStream = stream;
}

void CudaContext::restoreDefaultStream() {
    setCurrentStream(0);
}

string CudaContext::doubleToString(double value) const {
    stringstream s;
    s.precision(useDoublePrecision ? 16 : 8);
    s << scientific << value;
    if (!useDoublePrecision)
        s << "f";
    return s.str();
}

string CudaContext::intToString(int value) const {
    stringstream s;
    s << value;
    return s.str();
}

std::string CudaContext::getErrorString(CUresult result) {
    switch (result) {
        case CUDA_SUCCESS: return "CUDA_SUCCESS";
        case CUDA_ERROR_INVALID_VALUE: return "CUDA_ERROR_INVALID_VALUE";
        case CUDA_ERROR_OUT_OF_MEMORY: return "CUDA_ERROR_OUT_OF_MEMORY";
        case CUDA_ERROR_NOT_INITIALIZED: return "CUDA_ERROR_NOT_INITIALIZED";
        case CUDA_ERROR_DEINITIALIZED: return "CUDA_ERROR_DEINITIALIZED";
        case CUDA_ERROR_PROFILER_DISABLED: return "CUDA_ERROR_PROFILER_DISABLED";
        case CUDA_ERROR_PROFILER_NOT_INITIALIZED: return "CUDA_ERROR_PROFILER_NOT_INITIALIZED";
        case CUDA_ERROR_PROFILER_ALREADY_STARTED: return "CUDA_ERROR_PROFILER_ALREADY_STARTED";
        case CUDA_ERROR_PROFILER_ALREADY_STOPPED: return "CUDA_ERROR_PROFILER_ALREADY_STOPPED";
        case CUDA_ERROR_NO_DEVICE: return "CUDA_ERROR_NO_DEVICE";
        case CUDA_ERROR_INVALID_DEVICE: return "CUDA_ERROR_INVALID_DEVICE";
        case CUDA_ERROR_INVALID_IMAGE: return "CUDA_ERROR_INVALID_IMAGE";
        case CUDA_ERROR_INVALID_CONTEXT: return "CUDA_ERROR_INVALID_CONTEXT";
        case CUDA_ERROR_CONTEXT_ALREADY_CURRENT: return "CUDA_ERROR_CONTEXT_ALREADY_CURRENT";
        case CUDA_ERROR_MAP_FAILED: return "CUDA_ERROR_MAP_FAILED";
        case CUDA_ERROR_UNMAP_FAILED: return "CUDA_ERROR_UNMAP_FAILED";
        case CUDA_ERROR_ARRAY_IS_MAPPED: return "CUDA_ERROR_ARRAY_IS_MAPPED";
        case CUDA_ERROR_ALREADY_MAPPED: return "CUDA_ERROR_ALREADY_MAPPED";
        case CUDA_ERROR_NO_BINARY_FOR_GPU: return "CUDA_ERROR_NO_BINARY_FOR_GPU";
        case CUDA_ERROR_ALREADY_ACQUIRED: return "CUDA_ERROR_ALREADY_ACQUIRED";
        case CUDA_ERROR_NOT_MAPPED: return "CUDA_ERROR_NOT_MAPPED";
        case CUDA_ERROR_NOT_MAPPED_AS_ARRAY: return "CUDA_ERROR_NOT_MAPPED_AS_ARRAY";
        case CUDA_ERROR_NOT_MAPPED_AS_POINTER: return "CUDA_ERROR_NOT_MAPPED_AS_POINTER";
        case CUDA_ERROR_ECC_UNCORRECTABLE: return "CUDA_ERROR_ECC_UNCORRECTABLE";
        case CUDA_ERROR_UNSUPPORTED_LIMIT: return "CUDA_ERROR_UNSUPPORTED_LIMIT";
        case CUDA_ERROR_CONTEXT_ALREADY_IN_USE: return "CUDA_ERROR_CONTEXT_ALREADY_IN_USE";
        case CUDA_ERROR_PEER_ACCESS_UNSUPPORTED: return "CUDA_ERROR_PEER_ACCESS_UNSUPPORTED";
        case CUDA_ERROR_INVALID_SOURCE: return "CUDA_ERROR_INVALID_SOURCE";
        case CUDA_ERROR_FILE_NOT_FOUND: return "CUDA_ERROR_FILE_NOT_FOUND";
        case CUDA_ERROR_SHARED_OBJECT_SYMBOL_NOT_FOUND: return "CUDA_ERROR_SHARED_OBJECT_SYMBOL_NOT_FOUND";
        case CUDA_ERROR_SHARED_OBJECT_INIT_FAILED: return "CUDA_ERROR_SHARED_OBJECT_INIT_FAILED";
        case CUDA_ERROR_OPERATING_SYSTEM: return "CUDA_ERROR_OPERATING_SYSTEM";
        case CUDA_ERROR_INVALID_HANDLE: return "CUDA_ERROR_INVALID_HANDLE";
        case CUDA_ERROR_NOT_FOUND: return "CUDA_ERROR_NOT_FOUND";
        case CUDA_ERROR_NOT_READY: return "CUDA_ERROR_NOT_READY";
        case CUDA_ERROR_LAUNCH_FAILED: return "CUDA_ERROR_LAUNCH_FAILED";
        case CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES: return "CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES";
        case CUDA_ERROR_LAUNCH_TIMEOUT: return "CUDA_ERROR_LAUNCH_TIMEOUT";
        case CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING: return "CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING";
        case CUDA_ERROR_PEER_ACCESS_ALREADY_ENABLED: return "CUDA_ERROR_PEER_ACCESS_ALREADY_ENABLED";
        case CUDA_ERROR_PEER_ACCESS_NOT_ENABLED: return "CUDA_ERROR_PEER_ACCESS_NOT_ENABLED";
        case CUDA_ERROR_PRIMARY_CONTEXT_ACTIVE: return "CUDA_ERROR_PRIMARY_CONTEXT_ACTIVE";
        case CUDA_ERROR_CONTEXT_IS_DESTROYED: return "CUDA_ERROR_CONTEXT_IS_DESTROYED";
        case CUDA_ERROR_ASSERT: return "CUDA_ERROR_ASSERT";
        case CUDA_ERROR_TOO_MANY_PEERS: return "CUDA_ERROR_TOO_MANY_PEERS";
        case CUDA_ERROR_HOST_MEMORY_ALREADY_REGISTERED: return "CUDA_ERROR_HOST_MEMORY_ALREADY_REGISTERED";
        case CUDA_ERROR_HOST_MEMORY_NOT_REGISTERED: return "CUDA_ERROR_HOST_MEMORY_NOT_REGISTERED";
        case CUDA_ERROR_NOT_PERMITTED: return "CUDA_ERROR_NOT_PERMITTED";
        case CUDA_ERROR_NOT_SUPPORTED: return "CUDA_ERROR_NOT_SUPPORTED";
        case CUDA_ERROR_UNKNOWN: return "CUDA_ERROR_UNKNOWN";
    }
    return "CUDA error";
}

void CudaContext::executeKernel(CUfunction kernel, void** arguments, int threads, int blockSize, unsigned int sharedSize) {
    if (blockSize == -1)
        blockSize = ThreadBlockSize;
    int gridSize = std::min((threads+blockSize-1)/blockSize, numThreadBlocks);
    CUresult result = cuLaunchKernel(kernel, gridSize, 1, 1, blockSize, 1, 1, sharedSize, currentStream, arguments, NULL);
    if (result != CUDA_SUCCESS) {
        stringstream str;
        str<<"Error invoking kernel: "<<getErrorString(result)<<" ("<<result<<")";
        throw OpenMMException(str.str());
    }
}

int CudaContext::computeThreadBlockSize(double memory, bool preferShared) const {
    int maxShared = 16*1024;
    if (computeCapability >= 2.0 && preferShared)
        maxShared = 48*1024;
    int max = (int) (maxShared/memory);
    if (max < 64)
        return 32;
    int threads = 64;
    while (threads+64 < max)
        threads += 64;
    return threads;
}

void CudaContext::clearBuffer(CudaArray& array) {
    clearBuffer(array.getDevicePointer(), array.getSize()*array.getElementSize());
}

void CudaContext::clearBuffer(CUdeviceptr memory, int size) {
    int words = size/4;
    void* args[] = {&memory, &words};
    executeKernel(clearBufferKernel, args, words, 128);
}

void CudaContext::addAutoclearBuffer(CudaArray& array) {
    addAutoclearBuffer(array.getDevicePointer(), array.getSize()*array.getElementSize());
}

void CudaContext::addAutoclearBuffer(CUdeviceptr memory, int size) {
    autoclearBuffers.push_back(memory);
    autoclearBufferSizes.push_back(size/4);
}

void CudaContext::clearAutoclearBuffers() {
    int base = 0;
    int total = autoclearBufferSizes.size();
    while (total-base >= 6) {
        void* args[] = {&autoclearBuffers[base], &autoclearBufferSizes[base],
                        &autoclearBuffers[base+1], &autoclearBufferSizes[base+1],
                        &autoclearBuffers[base+2], &autoclearBufferSizes[base+2],
                        &autoclearBuffers[base+3], &autoclearBufferSizes[base+3],
                        &autoclearBuffers[base+4], &autoclearBufferSizes[base+4],
                        &autoclearBuffers[base+5], &autoclearBufferSizes[base+5]};
        executeKernel(clearSixBuffersKernel, args, max(max(max(max(max(autoclearBufferSizes[base], autoclearBufferSizes[base+1]), autoclearBufferSizes[base+2]), autoclearBufferSizes[base+3]), autoclearBufferSizes[base+4]), autoclearBufferSizes[base+5]), 128);
        base += 6;
    }
    if (total-base == 5) {
        void* args[] = {&autoclearBuffers[base], &autoclearBufferSizes[base],
                        &autoclearBuffers[base+1], &autoclearBufferSizes[base+1],
                        &autoclearBuffers[base+2], &autoclearBufferSizes[base+2],
                        &autoclearBuffers[base+3], &autoclearBufferSizes[base+3],
                        &autoclearBuffers[base+4], &autoclearBufferSizes[base+4]};
        executeKernel(clearFiveBuffersKernel, args, max(max(max(max(autoclearBufferSizes[base], autoclearBufferSizes[base+1]), autoclearBufferSizes[base+2]), autoclearBufferSizes[base+3]), autoclearBufferSizes[base+4]), 128);
    }
    else if (total-base == 4) {
        void* args[] = {&autoclearBuffers[base], &autoclearBufferSizes[base],
                        &autoclearBuffers[base+1], &autoclearBufferSizes[base+1],
                        &autoclearBuffers[base+2], &autoclearBufferSizes[base+2],
                        &autoclearBuffers[base+3], &autoclearBufferSizes[base+3]};
        executeKernel(clearFourBuffersKernel, args, max(max(max(autoclearBufferSizes[base], autoclearBufferSizes[base+1]), autoclearBufferSizes[base+2]), autoclearBufferSizes[base+3]), 128);
    }
    else if (total-base == 3) {
        void* args[] = {&autoclearBuffers[base], &autoclearBufferSizes[base],
                        &autoclearBuffers[base+1], &autoclearBufferSizes[base+1],
                        &autoclearBuffers[base+2], &autoclearBufferSizes[base+2]};
        executeKernel(clearThreeBuffersKernel, args, max(max(autoclearBufferSizes[base], autoclearBufferSizes[base+1]), autoclearBufferSizes[base+2]), 128);
    }
    else if (total-base == 2) {
        void* args[] = {&autoclearBuffers[base], &autoclearBufferSizes[base],
                        &autoclearBuffers[base+1], &autoclearBufferSizes[base+1]};
        executeKernel(clearTwoBuffersKernel, args, max(autoclearBufferSizes[base], autoclearBufferSizes[base+1]), 128);
    }
    else if (total-base == 1) {
        clearBuffer(autoclearBuffers[base], autoclearBufferSizes[base]*4);
    }
}

/**
 * This class ensures that atom reordering doesn't break virtual sites.
 */
class CudaContext::VirtualSiteInfo : public CudaForceInfo {
public:
    VirtualSiteInfo(const System& system) {
        for (int i = 0; i < system.getNumParticles(); i++) {
            if (system.isVirtualSite(i)) {
                siteTypes.push_back(&typeid(system.getVirtualSite(i)));
                vector<int> particles;
                particles.push_back(i);
                for (int j = 0; j < system.getVirtualSite(i).getNumParticles(); j++)
                    particles.push_back(system.getVirtualSite(i).getParticle(j));
                siteParticles.push_back(particles);
                vector<double> weights;
                if (dynamic_cast<const TwoParticleAverageSite*>(&system.getVirtualSite(i)) != NULL) {
                    // A two particle average.

                    const TwoParticleAverageSite& site = dynamic_cast<const TwoParticleAverageSite&>(system.getVirtualSite(i));
                    weights.push_back(site.getWeight(0));
                    weights.push_back(site.getWeight(1));
                }
                else if (dynamic_cast<const ThreeParticleAverageSite*>(&system.getVirtualSite(i)) != NULL) {
                    // A three particle average.

                    const ThreeParticleAverageSite& site = dynamic_cast<const ThreeParticleAverageSite&>(system.getVirtualSite(i));
                    weights.push_back(site.getWeight(0));
                    weights.push_back(site.getWeight(1));
                    weights.push_back(site.getWeight(2));
                }
                else if (dynamic_cast<const OutOfPlaneSite*>(&system.getVirtualSite(i)) != NULL) {
                    // An out of plane site.

                    const OutOfPlaneSite& site = dynamic_cast<const OutOfPlaneSite&>(system.getVirtualSite(i));
                    weights.push_back(site.getWeight12());
                    weights.push_back(site.getWeight13());
                    weights.push_back(site.getWeightCross());
                }
                siteWeights.push_back(weights);
            }
        }
    }
    int getNumParticleGroups() {
        return siteTypes.size();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
        particles = siteParticles[index];
    }
    bool areGroupsIdentical(int group1, int group2) {
        if (siteTypes[group1] != siteTypes[group2])
            return false;
        int numParticles = siteWeights[group1].size();
        if (siteWeights[group2].size() != numParticles)
            return false;
        for (int i = 0; i < numParticles; i++)
            if (siteWeights[group1][i] != siteWeights[group2][i])
                return false;
        return true;
    }
private:
    vector<const type_info*> siteTypes;
    vector<vector<int> > siteParticles;
    vector<vector<double> > siteWeights;
};

void CudaContext::findMoleculeGroups() {
    // The first time this is called, we need to identify all the molecules in the system.

    if (moleculeGroups.size() == 0) {
        // Add a ForceInfo that makes sure reordering doesn't break virtual sites.

        addForce(new VirtualSiteInfo(system));

        // First make a list of every other atom to which each atom is connect by a constraint or force group.

        vector<vector<int> > atomBonds(system.getNumParticles());
        for (int i = 0; i < system.getNumConstraints(); i++) {
            int particle1, particle2;
            double distance;
            system.getConstraintParameters(i, particle1, particle2, distance);
            atomBonds[particle1].push_back(particle2);
            atomBonds[particle2].push_back(particle1);
        }
        for (int i = 0; i < (int) forces.size(); i++) {
            for (int j = 0; j < forces[i]->getNumParticleGroups(); j++) {
                vector<int> particles;
                forces[i]->getParticlesInGroup(j, particles);
                for (int k = 0; k < (int) particles.size(); k++)
                    for (int m = 0; m < (int) particles.size(); m++)
                        if (k != m)
                            atomBonds[particles[k]].push_back(particles[m]);
            }
        }

        // Now identify atoms by which molecule they belong to.

        vector<vector<int> > atomIndices = ContextImpl::findMolecules(numAtoms, atomBonds);
        int numMolecules = atomIndices.size();
        vector<int> atomMolecule(numAtoms);
        for (int i = 0; i < (int) atomIndices.size(); i++)
            for (int j = 0; j < (int) atomIndices[i].size(); j++)
                atomMolecule[atomIndices[i][j]] = i;

        // Construct a description of each molecule.

        molecules.resize(numMolecules);
        for (int i = 0; i < numMolecules; i++) {
            molecules[i].atoms = atomIndices[i];
            molecules[i].groups.resize(forces.size());
        }
        for (int i = 0; i < system.getNumConstraints(); i++) {
            int particle1, particle2;
            double distance;
            system.getConstraintParameters(i, particle1, particle2, distance);
            molecules[atomMolecule[particle1]].constraints.push_back(i);
        }
        for (int i = 0; i < (int) forces.size(); i++)
            for (int j = 0; j < forces[i]->getNumParticleGroups(); j++) {
                vector<int> particles;
                forces[i]->getParticlesInGroup(j, particles);
                if (particles.size() > 0)
                    molecules[atomMolecule[particles[0]]].groups[i].push_back(j);
            }
    }

    // Sort them into groups of identical molecules.

    vector<Molecule> uniqueMolecules;
    vector<vector<int> > moleculeInstances;
    vector<vector<int> > moleculeOffsets;
    for (int molIndex = 0; molIndex < (int) molecules.size(); molIndex++) {
        Molecule& mol = molecules[molIndex];

        // See if it is identical to another molecule.

        bool isNew = true;
        for (int j = 0; j < (int) uniqueMolecules.size() && isNew; j++) {
            Molecule& mol2 = uniqueMolecules[j];
            bool identical = (mol.atoms.size() == mol2.atoms.size() && mol.constraints.size() == mol2.constraints.size());

            // See if the atoms are identical.

            int atomOffset = mol2.atoms[0]-mol.atoms[0];
            for (int i = 0; i < (int) mol.atoms.size() && identical; i++) {
                if (mol.atoms[i] != mol2.atoms[i]-atomOffset || system.getParticleMass(mol.atoms[i]) != system.getParticleMass(mol2.atoms[i]))
                    identical = false;
                for (int k = 0; k < (int) forces.size(); k++)
                    if (!forces[k]->areParticlesIdentical(mol.atoms[i], mol2.atoms[i]))
                        identical = false;
            }

            // See if the constraints are identical.

            for (int i = 0; i < (int) mol.constraints.size() && identical; i++) {
                int c1particle1, c1particle2, c2particle1, c2particle2;
                double distance1, distance2;
                system.getConstraintParameters(mol.constraints[i], c1particle1, c1particle2, distance1);
                system.getConstraintParameters(mol2.constraints[i], c2particle1, c2particle2, distance2);
                if (c1particle1 != c2particle1-atomOffset || c1particle2 != c2particle2-atomOffset || distance1 != distance2)
                    identical = false;
            }

            // See if the force groups are identical.

            for (int i = 0; i < (int) forces.size() && identical; i++) {
                if (mol.groups[i].size() != mol2.groups[i].size())
                    identical = false;
                for (int k = 0; k < (int) mol.groups[i].size() && identical; k++)
                    if (!forces[i]->areGroupsIdentical(mol.groups[i][k], mol2.groups[i][k]))
                        identical = false;
            }
            if (identical) {
                moleculeInstances[j].push_back(molIndex);
                moleculeOffsets[j].push_back(mol.atoms[0]);
                isNew = false;
            }
        }
        if (isNew) {
            uniqueMolecules.push_back(mol);
            moleculeInstances.push_back(vector<int>());
            moleculeInstances[moleculeInstances.size()-1].push_back(molIndex);
            moleculeOffsets.push_back(vector<int>());
            moleculeOffsets[moleculeOffsets.size()-1].push_back(mol.atoms[0]);
        }
    }
    moleculeGroups.resize(moleculeInstances.size());
    for (int i = 0; i < (int) moleculeInstances.size(); i++)
    {
        moleculeGroups[i].instances = moleculeInstances[i];
        moleculeGroups[i].offsets = moleculeOffsets[i];
        vector<int>& atoms = uniqueMolecules[i].atoms;
        moleculeGroups[i].atoms.resize(atoms.size());
        for (int j = 0; j < (int) atoms.size(); j++)
            moleculeGroups[i].atoms[j] = atoms[j]-atoms[0];
    }
}

void CudaContext::invalidateMolecules() {
    if (numAtoms == 0 || nonbonded == NULL || !nonbonded->getUseCutoff())
        return;
    bool valid = true;
    for (int group = 0; valid && group < (int) moleculeGroups.size(); group++) {
        MoleculeGroup& mol = moleculeGroups[group];
        vector<int>& instances = mol.instances;
        vector<int>& offsets = mol.offsets;
        vector<int>& atoms = mol.atoms;
        int numMolecules = instances.size();
        Molecule& m1 = molecules[instances[0]];
        int offset1 = offsets[0];
        for (int j = 1; valid && j < numMolecules; j++) {
            // See if the atoms are identical.

            Molecule& m2 = molecules[instances[j]];
            int offset2 = offsets[j];
            for (int i = 0; i < (int) atoms.size() && valid; i++) {
                for (int k = 0; k < (int) forces.size(); k++)
                    if (!forces[k]->areParticlesIdentical(atoms[i]+offset1, atoms[i]+offset2))
                        valid = false;
            }

            // See if the force groups are identical.

            for (int i = 0; i < (int) forces.size() && valid; i++) {
                for (int k = 0; k < (int) m1.groups[i].size() && valid; k++)
                    if (!forces[i]->areGroupsIdentical(m1.groups[i][k], m2.groups[i][k]))
                        valid = false;
            }
        }
    }
    if (valid)
        return;

    // The list of which molecules are identical is no longer valid.  We need to restore the
    // atoms to their original order, rebuild the list of identical molecules, and sort them
    // again.

    vector<int4> newCellOffsets(numAtoms);
    if (useDoublePrecision) {
        vector<double4> oldPosq(paddedNumAtoms);
        vector<double4> newPosq(paddedNumAtoms, make_double4(0, 0, 0, 0));
        vector<double4> oldVelm(paddedNumAtoms);
        vector<double4> newVelm(paddedNumAtoms, make_double4(0, 0, 0, 0));
        posq->download(oldPosq);
        velm->download(oldVelm);
        for (int i = 0; i < numAtoms; i++) {
            int index = atomIndex[i];
            newPosq[index] = oldPosq[i];
            newVelm[index] = oldVelm[i];
            newCellOffsets[index] = posCellOffsets[i];
        }
        posq->upload(newPosq);
        velm->upload(newVelm);
    }
    else if (useMixedPrecision) {
        vector<float4> oldPosq(paddedNumAtoms);
        vector<float4> newPosq(paddedNumAtoms, make_float4(0, 0, 0, 0));
        vector<float4> oldPosqCorrection(paddedNumAtoms);
        vector<float4> newPosqCorrection(paddedNumAtoms, make_float4(0, 0, 0, 0));
        vector<double4> oldVelm(paddedNumAtoms);
        vector<double4> newVelm(paddedNumAtoms, make_double4(0, 0, 0, 0));
        posq->download(oldPosq);
        velm->download(oldVelm);
        for (int i = 0; i < numAtoms; i++) {
            int index = atomIndex[i];
            newPosq[index] = oldPosq[i];
            newPosqCorrection[index] = oldPosqCorrection[i];
            newVelm[index] = oldVelm[i];
            newCellOffsets[index] = posCellOffsets[i];
        }
        posq->upload(newPosq);
        posqCorrection->upload(newPosqCorrection);
        velm->upload(newVelm);
    }
    else {
        vector<float4> oldPosq(paddedNumAtoms);
        vector<float4> newPosq(paddedNumAtoms, make_float4(0, 0, 0, 0));
        vector<float4> oldVelm(paddedNumAtoms);
        vector<float4> newVelm(paddedNumAtoms, make_float4(0, 0, 0, 0));
        posq->download(oldPosq);
        velm->download(oldVelm);
        for (int i = 0; i < numAtoms; i++) {
            int index = atomIndex[i];
            newPosq[index] = oldPosq[i];
            newVelm[index] = oldVelm[i];
            newCellOffsets[index] = posCellOffsets[i];
        }
        posq->upload(newPosq);
        velm->upload(newVelm);
    }
    for (int i = 0; i < numAtoms; i++) {
        atomIndex[i] = i;
        posCellOffsets[i] = newCellOffsets[i];
    }
    atomIndexDevice->upload(atomIndex);
    findMoleculeGroups();
    for (int i = 0; i < (int) reorderListeners.size(); i++)
        reorderListeners[i]->execute();
    reorderAtoms();
}

void CudaContext::reorderAtoms() {
    atomsWereReordered = false;
    if (numAtoms == 0 || nonbonded == NULL || !nonbonded->getUseCutoff() || stepsSinceReorder < 250) {
        stepsSinceReorder++;
        return;
    }
    atomsWereReordered = true;
    stepsSinceReorder = 0;
    if (useDoublePrecision)
        reorderAtomsImpl<double, double4, double, double4>();
    else if (useMixedPrecision)
        reorderAtomsImpl<float, float4, double, double4>();
    else
        reorderAtomsImpl<float, float4, float, float4>();
}

template <class Real, class Real4, class Mixed, class Mixed4>
void CudaContext::reorderAtomsImpl() {
    // Find the range of positions and the number of bins along each axis.

    Real4 padding = {0, 0, 0, 0};
    vector<Real4> oldPosq(paddedNumAtoms, padding);
    vector<Real4> oldPosqCorrection(paddedNumAtoms, padding);
    Mixed4 paddingMixed = {0, 0, 0, 0};
    vector<Mixed4> oldVelm(paddedNumAtoms, paddingMixed);
    posq->download(oldPosq);
    velm->download(oldVelm);
    if (useMixedPrecision)
        posqCorrection->download(oldPosqCorrection);
    Real minx = oldPosq[0].x, maxx = oldPosq[0].x;
    Real miny = oldPosq[0].y, maxy = oldPosq[0].y;
    Real minz = oldPosq[0].z, maxz = oldPosq[0].z;
    if (nonbonded->getUsePeriodic()) {
        minx = miny = minz = 0.0;
        maxx = periodicBoxSize.x;
        maxy = periodicBoxSize.y;
        maxz = periodicBoxSize.z;
    }
    else {
        for (int i = 1; i < numAtoms; i++) {
            const Real4& pos = oldPosq[i];
            minx = min(minx, pos.x);
            maxx = max(maxx, pos.x);
            miny = min(miny, pos.y);
            maxy = max(maxy, pos.y);
            minz = min(minz, pos.z);
            maxz = max(maxz, pos.z);
        }
    }

    // Loop over each group of identical molecules and reorder them.

    vector<int> originalIndex(numAtoms);
    vector<Real4> newPosq(paddedNumAtoms);
    vector<Real4> newPosqCorrection(paddedNumAtoms);
    vector<Mixed4> newVelm(paddedNumAtoms);
    vector<int4> newCellOffsets(numAtoms);
    for (int group = 0; group < (int) moleculeGroups.size(); group++) {
        // Find the center of each molecule.

        MoleculeGroup& mol = moleculeGroups[group];
        int numMolecules = mol.offsets.size();
        vector<int>& atoms = mol.atoms;
        vector<Real4> molPos(numMolecules);
        Real invNumAtoms = (Real) (1.0/atoms.size());
        for (int i = 0; i < numMolecules; i++) {
            molPos[i].x = 0.0f;
            molPos[i].y = 0.0f;
            molPos[i].z = 0.0f;
            for (int j = 0; j < (int)atoms.size(); j++) {
                int atom = atoms[j]+mol.offsets[i];
                const Real4& pos = oldPosq[atom];
                molPos[i].x += pos.x;
                molPos[i].y += pos.y;
                molPos[i].z += pos.z;
            }
            molPos[i].x *= invNumAtoms;
            molPos[i].y *= invNumAtoms;
            molPos[i].z *= invNumAtoms;
            if (molPos[i].x != molPos[i].x)
                throw OpenMMException("Particle coordinate is nan");
        }
        if (nonbonded->getUsePeriodic()) {
            // Move each molecule position into the same box.

            for (int i = 0; i < numMolecules; i++) {
                Real4 center = molPos[i];
                int zcell = (int) floor(center.z*invPeriodicBoxSize.z);
                center.x -= zcell*periodicBoxVecZ.x;
                center.y -= zcell*periodicBoxVecZ.y;
                center.z -= zcell*periodicBoxVecZ.z;
                int ycell = (int) floor(center.y*invPeriodicBoxSize.y);
                center.x -= ycell*periodicBoxVecY.x;
                center.y -= ycell*periodicBoxVecY.y;
                int xcell = (int) floor(center.x*invPeriodicBoxSize.x);
                center.x -= xcell*periodicBoxVecX.x;
                if (xcell != 0 || ycell != 0 || zcell != 0) {
                    Real dx = molPos[i].x-center.x;
                    Real dy = molPos[i].y-center.y;
                    Real dz = molPos[i].z-center.z;
                    molPos[i] = center;
                    for (int j = 0; j < (int) atoms.size(); j++) {
                        int atom = atoms[j]+mol.offsets[i];
                        Real4 p = oldPosq[atom];
                        p.x -= dx;
                        p.y -= dy;
                        p.z -= dz;
                        oldPosq[atom] = p;
                        posCellOffsets[atom].x -= xcell;
                        posCellOffsets[atom].y -= ycell;
                        posCellOffsets[atom].z -= zcell;
                    }
                }
            }
        }

        // Select a bin for each molecule, then sort them by bin.

        bool useHilbert = (numMolecules > 5000 || atoms.size() > 8); // For small systems, a simple zigzag curve works better than a Hilbert curve.
        Real binWidth;
        if (useHilbert)
            binWidth = (Real) (max(max(maxx-minx, maxy-miny), maxz-minz)/255.0);
        else
            binWidth = (Real) (0.2*nonbonded->getMaxCutoffDistance());
        Real invBinWidth = (Real) (1.0/binWidth);
        int xbins = 1 + (int) ((maxx-minx)*invBinWidth);
        int ybins = 1 + (int) ((maxy-miny)*invBinWidth);
        vector<pair<int, int> > molBins(numMolecules);
        bitmask_t coords[3];
        for (int i = 0; i < numMolecules; i++) {
            int x = (int) ((molPos[i].x-minx)*invBinWidth);
            int y = (int) ((molPos[i].y-miny)*invBinWidth);
            int z = (int) ((molPos[i].z-minz)*invBinWidth);
            int bin;
            if (useHilbert) {
                coords[0] = x;
                coords[1] = y;
                coords[2] = z;
                bin = (int) hilbert_c2i(3, 8, coords);
            }
            else {
                int yodd = y&1;
                int zodd = z&1;
                bin = z*xbins*ybins;
                bin += (zodd ? ybins-y : y)*xbins;
                bin += (yodd ? xbins-x : x);
            }
            molBins[i] = pair<int, int>(bin, i);
        }
        sort(molBins.begin(), molBins.end());

        // Reorder the atoms.

        for (int i = 0; i < numMolecules; i++) {
            for (int j = 0; j < (int)atoms.size(); j++) {
                int oldIndex = mol.offsets[molBins[i].second]+atoms[j];
                int newIndex = mol.offsets[i]+atoms[j];
                originalIndex[newIndex] = atomIndex[oldIndex];
                newPosq[newIndex] = oldPosq[oldIndex];
                if (useMixedPrecision)
                    newPosqCorrection[newIndex] = oldPosqCorrection[oldIndex];
                newVelm[newIndex] = oldVelm[oldIndex];
                newCellOffsets[newIndex] = posCellOffsets[oldIndex];
            }
        }
    }

    // Update the streams.

    for (int i = 0; i < numAtoms; i++) {
        atomIndex[i] = originalIndex[i];
        posCellOffsets[i] = newCellOffsets[i];
    }
    posq->upload(newPosq);
    if (useMixedPrecision)
        posqCorrection->upload(newPosqCorrection);
    velm->upload(newVelm);
    atomIndexDevice->upload(atomIndex);
    for (int i = 0; i < (int) reorderListeners.size(); i++)
        reorderListeners[i]->execute();
}

void CudaContext::addReorderListener(ReorderListener* listener) {
    reorderListeners.push_back(listener);
}

void CudaContext::addPreComputation(ForcePreComputation* computation) {
    preComputations.push_back(computation);
}

void CudaContext::addPostComputation(ForcePostComputation* computation) {
    postComputations.push_back(computation);
}

void CudaContext::addEnergyParameterDerivative(const string& param) {
    // See if this parameter has already been registered.
    
    for (int i = 0; i < energyParamDerivNames.size(); i++)
        if (param == energyParamDerivNames[i])
            return;
    energyParamDerivNames.push_back(param);
}

struct CudaContext::WorkThread::ThreadData {
    ThreadData(std::queue<CudaContext::WorkTask*>& tasks, bool& waiting,  bool& finished,
            pthread_mutex_t& queueLock, pthread_cond_t& waitForTaskCondition, pthread_cond_t& queueEmptyCondition) :
        tasks(tasks), waiting(waiting), finished(finished), queueLock(queueLock),
        waitForTaskCondition(waitForTaskCondition), queueEmptyCondition(queueEmptyCondition) {
    }
    std::queue<CudaContext::WorkTask*>& tasks;
    bool& waiting;
    bool& finished;
    pthread_mutex_t& queueLock;
    pthread_cond_t& waitForTaskCondition;
    pthread_cond_t& queueEmptyCondition;
};

static void* threadBody(void* args) {
    CudaContext::WorkThread::ThreadData& data = *reinterpret_cast<CudaContext::WorkThread::ThreadData*>(args);
    while (!data.finished || data.tasks.size() > 0) {
        pthread_mutex_lock(&data.queueLock);
        while (data.tasks.empty() && !data.finished) {
            data.waiting = true;
            pthread_cond_signal(&data.queueEmptyCondition);
            pthread_cond_wait(&data.waitForTaskCondition, &data.queueLock);
        }
        CudaContext::WorkTask* task = NULL;
        if (!data.tasks.empty()) {
            data.waiting = false;
            task = data.tasks.front();
            data.tasks.pop();
        }
        pthread_mutex_unlock(&data.queueLock);
        if (task != NULL) {
            task->execute();
            delete task;
        }
    }
    data.waiting = true;
    pthread_cond_signal(&data.queueEmptyCondition);
    delete &data;
    return 0;
}

CudaContext::WorkThread::WorkThread() : waiting(true), finished(false) {
    pthread_mutex_init(&queueLock, NULL);
    pthread_cond_init(&waitForTaskCondition, NULL);
    pthread_cond_init(&queueEmptyCondition, NULL);
    ThreadData* data = new ThreadData(tasks, waiting, finished, queueLock, waitForTaskCondition, queueEmptyCondition);
    pthread_create(&thread, NULL, threadBody, data);
}

CudaContext::WorkThread::~WorkThread() {
    pthread_mutex_lock(&queueLock);
    finished = true;
    pthread_cond_broadcast(&waitForTaskCondition);
    pthread_mutex_unlock(&queueLock);
    pthread_join(thread, NULL);
    pthread_mutex_destroy(&queueLock);
    pthread_cond_destroy(&waitForTaskCondition);
    pthread_cond_destroy(&queueEmptyCondition);
}

void CudaContext::WorkThread::addTask(CudaContext::WorkTask* task) {
    pthread_mutex_lock(&queueLock);
    tasks.push(task);
    waiting = false;
    pthread_cond_signal(&waitForTaskCondition);
    pthread_mutex_unlock(&queueLock);
}

bool CudaContext::WorkThread::isWaiting() {
    return waiting;
}

bool CudaContext::WorkThread::isFinished() {
    return finished;
}

void CudaContext::WorkThread::flush() {
    pthread_mutex_lock(&queueLock);
    while (!waiting)
       pthread_cond_wait(&queueEmptyCondition, &queueLock);
    pthread_mutex_unlock(&queueLock);
}


vector<int> CudaContext::getDevicePrecedence() {
    int numDevices;
    CUdevice thisDevice;
    string errorMessage = "Error initializing Context";
    vector<pair<pair<int, int>, int> > devices;

    CHECK_RESULT(cuDeviceGetCount(&numDevices));
    for (int i = 0; i < numDevices; i++) {
        CHECK_RESULT(cuDeviceGet(&thisDevice, i));
        int major, minor, clock, multiprocessors, speed;
        CHECK_RESULT(cuDeviceComputeCapability(&major, &minor, thisDevice));
        if (major == 1 && minor < 2)
            continue;

        if ((useDoublePrecision || useMixedPrecision) && (major+0.1*minor < 1.3))
            continue;

        CHECK_RESULT(cuDeviceGetAttribute(&clock, CU_DEVICE_ATTRIBUTE_CLOCK_RATE, thisDevice));
        CHECK_RESULT(cuDeviceGetAttribute(&multiprocessors, CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT, thisDevice));
        speed = clock*multiprocessors;
        pair<int, int> deviceProperties = std::make_pair(major, speed);
        devices.push_back(std::make_pair(deviceProperties, -i));
    }

    // sort first by compute capability (higher is better), then speed
    // (higher is better), and finally device index (lower is better)
    std::sort(devices.begin(), devices.end());
    std::reverse(devices.begin(), devices.end());

    vector<int> precedence;
    for (int i = 0; i < static_cast<int>(devices.size()); i++) {
        precedence.push_back(-devices[i].second);
    }

    return precedence;
}
