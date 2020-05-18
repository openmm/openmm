/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2019 Stanford University and the Authors.      *
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
#include "CudaEvent.h"
#include "CudaIntegrationUtilities.h"
#include "CudaKernels.h"
#include "CudaKernelSources.h"
#include "CudaNonbondedUtilities.h"
#include "CudaProgram.h"
#include "openmm/common/ComputeArray.h"
#include "SHA1.h"
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
        const string& tempDir, const std::string& hostCompiler, CudaPlatform::PlatformData& platformData, CudaContext* originalContext) : ComputeContext(system), currentStream(0),
        platformData(platformData), contextIsValid(false), hasAssignedPosqCharges(false),
        hasCompilerKernel(false), isNvccAvailable(false), pinnedBuffer(NULL), integration(NULL), expression(NULL), bonded(NULL), nonbonded(NULL) {
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
    string errorMessage = "Error initializing Context";
    if (originalContext == NULL) {
        isLinkedContext = false;
        int numDevices;
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
        if (this->deviceIndex == -1) {
            if (deviceIndex != -1)
                throw OpenMMException("The requested CUDA device could not be loaded");
            else
                throw OpenMMException("No compatible CUDA device is available");
        }
    }
    else {
        isLinkedContext = true;
        context = originalContext->context;
        this->deviceIndex = originalContext->deviceIndex;
        this->device = originalContext->device;
    }

    int major, minor;
    CHECK_RESULT(cuDeviceGetAttribute(&major, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, device));
    CHECK_RESULT(cuDeviceGetAttribute(&minor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, device));
    int numThreadBlocksPerComputeUnit = (major == 6 ? 4 : 6);
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
    if (cudaDriverVersion >= 9000) {
        compilationDefines["SYNC_WARPS"] = "__syncwarp();";
        compilationDefines["SHFL(var, srcLane)"] = "__shfl_sync(0xffffffff, var, srcLane);";
        compilationDefines["BALLOT(var)"] = "__ballot_sync(0xffffffff, var);";
    }
    else {
        compilationDefines["SYNC_WARPS"] = "";
        compilationDefines["SHFL(var, srcLane)"] = "__shfl(var, srcLane);";
        compilationDefines["BALLOT(var)"] = "__ballot(var);";
    }
    if (useDoublePrecision) {
        posq.initialize<double4>(*this, paddedNumAtoms, "posq");
        velm.initialize<double4>(*this, paddedNumAtoms, "velm");
        compilationDefines["USE_DOUBLE_PRECISION"] = "1";
        compilationDefines["make_real2"] = "make_double2";
        compilationDefines["make_real3"] = "make_double3";
        compilationDefines["make_real4"] = "make_double4";
        compilationDefines["make_mixed2"] = "make_double2";
        compilationDefines["make_mixed3"] = "make_double3";
        compilationDefines["make_mixed4"] = "make_double4";
    }
    else if (useMixedPrecision) {
        posq.initialize<float4>(*this, paddedNumAtoms, "posq");
        posqCorrection.initialize<float4>(*this, paddedNumAtoms, "posqCorrection");
        velm.initialize<double4>(*this, paddedNumAtoms, "velm");
        compilationDefines["USE_MIXED_PRECISION"] = "1";
        compilationDefines["make_real2"] = "make_float2";
        compilationDefines["make_real3"] = "make_float3";
        compilationDefines["make_real4"] = "make_float4";
        compilationDefines["make_mixed2"] = "make_double2";
        compilationDefines["make_mixed3"] = "make_double3";
        compilationDefines["make_mixed4"] = "make_double4";
    }
    else {
        posq.initialize<float4>(*this, paddedNumAtoms, "posq");
        velm.initialize<float4>(*this, paddedNumAtoms, "velm");
        compilationDefines["make_real2"] = "make_float2";
        compilationDefines["make_real3"] = "make_float3";
        compilationDefines["make_real4"] = "make_float4";
        compilationDefines["make_mixed2"] = "make_float2";
        compilationDefines["make_mixed3"] = "make_float3";
        compilationDefines["make_mixed4"] = "make_float4";
    }
    force.initialize<long long>(*this, paddedNumAtoms*3, "force");
    posCellOffsets.resize(paddedNumAtoms, mm_int4(0, 0, 0, 0));
    atomIndexDevice.initialize<int>(*this, paddedNumAtoms, "atomIndex");
    atomIndex.resize(paddedNumAtoms);
    for (int i = 0; i < paddedNumAtoms; ++i)
        atomIndex[i] = i;
    atomIndexDevice.upload(atomIndex);

    // Create utility kernels that are used in multiple places.

    CUmodule utilities = createModule(CudaKernelSources::vectorOps+CudaKernelSources::utilities);
    clearBufferKernel = getKernel(utilities, "clearBuffer");
    clearTwoBuffersKernel = getKernel(utilities, "clearTwoBuffers");
    clearThreeBuffersKernel = getKernel(utilities, "clearThreeBuffers");
    clearFourBuffersKernel = getKernel(utilities, "clearFourBuffers");
    clearFiveBuffersKernel = getKernel(utilities, "clearFiveBuffers");
    clearSixBuffersKernel = getKernel(utilities, "clearSixBuffers");
    reduceEnergyKernel = getKernel(utilities, "reduceEnergy");
    setChargesKernel = getKernel(utilities, "setCharges");

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

    // Create utilities objects.

    bonded = new CudaBondedUtilities(*this);
    nonbonded = new CudaNonbondedUtilities(*this);
    integration = new CudaIntegrationUtilities(*this, system);
    expression = new CudaExpressionUtilities(*this);
}

CudaContext::~CudaContext() {
    setAsCurrent();
    for (auto force : forces)
        delete force;
    for (auto listener : reorderListeners)
        delete listener;
    for (auto computation : preComputations)
        delete computation;
    for (auto computation : postComputations)
        delete computation;
    if (pinnedBuffer != NULL)
        cuMemFreeHost(pinnedBuffer);
    if (integration != NULL)
        delete integration;
    if (expression != NULL)
        delete expression;
    if (bonded != NULL)
        delete bonded;
    if (nonbonded != NULL)
        delete nonbonded;
    string errorMessage = "Error deleting Context";
    if (contextIsValid && !isLinkedContext) {
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
        energyBuffer.initialize<double>(*this, numEnergyBuffers, "energyBuffer");
        energySum.initialize<double>(*this, 1, "energySum");
        int pinnedBufferSize = max(paddedNumAtoms*4, numEnergyBuffers);
        CHECK_RESULT(cuMemHostAlloc(&pinnedBuffer, pinnedBufferSize*sizeof(double), 0));
    }
    else if (useMixedPrecision) {
        energyBuffer.initialize<double>(*this, numEnergyBuffers, "energyBuffer");
        energySum.initialize<double>(*this, 1, "energySum");
        int pinnedBufferSize = max(paddedNumAtoms*4, numEnergyBuffers);
        CHECK_RESULT(cuMemHostAlloc(&pinnedBuffer, pinnedBufferSize*sizeof(double), 0));
    }
    else {
        energyBuffer.initialize<float>(*this, numEnergyBuffers, "energyBuffer");
        energySum.initialize<float>(*this, 1, "energySum");
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
    velm.upload(pinnedBuffer);
    bonded->initialize(system);
    addAutoclearBuffer(force.getDevicePointer(), force.getSize()*force.getElementSize());
    addAutoclearBuffer(energyBuffer.getDevicePointer(), energyBuffer.getSize()*energyBuffer.getElementSize());
    int numEnergyParamDerivs = energyParamDerivNames.size();
    if (numEnergyParamDerivs > 0) {
        if (useDoublePrecision || useMixedPrecision)
            energyParamDerivBuffer.initialize<double>(*this, numEnergyParamDerivs*numEnergyBuffers, "energyParamDerivBuffer");
        else
            energyParamDerivBuffer.initialize<float>(*this, numEnergyParamDerivs*numEnergyBuffers, "energyParamDerivBuffer");
        addAutoclearBuffer(energyParamDerivBuffer);
    }
    findMoleculeGroups();
    nonbonded->initialize(system);
}

void CudaContext::initializeContexts() {
    getPlatformData().initializeContexts(system);
}

void CudaContext::setAsCurrent() {
    if (contextIsValid)
        cuCtxSetCurrent(context);
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
    for (auto& pair : compilationDefines) {
        src << "#define " << pair.first;
        if (!pair.second.empty())
            src << " " << pair.second;
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
    src << CudaKernelSources::common << endl;
    for (auto& pair : defines) {
        src << "#define " << pair.first;
        if (!pair.second.empty())
            src << " " << pair.second;
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

    if (hasCompilerKernel) {
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

CudaArray* CudaContext::createArray() {
    return new CudaArray();
}

ComputeEvent CudaContext::createEvent() {
    return shared_ptr<ComputeEventImpl>(new CudaEvent(*this));
}

ComputeProgram CudaContext::compileProgram(const std::string source, const std::map<std::string, std::string>& defines) {
    CUmodule module = createModule(CudaKernelSources::vectorOps+source, defines);
    return shared_ptr<ComputeProgramImpl>(new CudaProgram(*this, module));
}

CudaArray& CudaContext::unwrap(ArrayInterface& array) const {
    CudaArray* cuarray;
    ComputeArray* wrapper = dynamic_cast<ComputeArray*>(&array);
    if (wrapper != NULL)
        cuarray = dynamic_cast<CudaArray*>(&wrapper->getArray());
    else
        cuarray = dynamic_cast<CudaArray*>(&array);
    if (cuarray == NULL)
        throw OpenMMException("Array argument is not an CudaArray");
    return *cuarray;
}

std::string CudaContext::getErrorString(CUresult result) {
    const char* message;
    if (cuGetErrorName(result, &message) == CUDA_SUCCESS)
        return string(message);
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

void CudaContext::clearBuffer(ArrayInterface& array) {
    clearBuffer(unwrap(array).getDevicePointer(), array.getSize()*array.getElementSize());
}

void CudaContext::clearBuffer(CUdeviceptr memory, int size) {
    int words = size/4;
    void* args[] = {&memory, &words};
    executeKernel(clearBufferKernel, args, words, 128);
}

void CudaContext::addAutoclearBuffer(ArrayInterface& array) {
    addAutoclearBuffer(unwrap(array).getDevicePointer(), array.getSize()*array.getElementSize());
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

double CudaContext::reduceEnergy() {
    int bufferSize = energyBuffer.getSize();
    int workGroupSize  = 512;
    void* args[] = {&energyBuffer.getDevicePointer(), &energySum.getDevicePointer(), &bufferSize, &workGroupSize};
    executeKernel(reduceEnergyKernel, args, workGroupSize, workGroupSize, workGroupSize*energyBuffer.getElementSize());
    energySum.download(pinnedBuffer);
    if (getUseDoublePrecision() || getUseMixedPrecision())
        return *((double*) pinnedBuffer);
    else
        return *((float*) pinnedBuffer);
}

void CudaContext::setCharges(const vector<double>& charges) {
    if (!chargeBuffer.isInitialized())
        chargeBuffer.initialize(*this, numAtoms, useDoublePrecision ? sizeof(double) : sizeof(float), "chargeBuffer");
    vector<double> c(numAtoms);
    for (int i = 0; i < numAtoms; i++)
        c[i] = charges[i];
    chargeBuffer.upload(c, true);
    void* args[] = {&chargeBuffer.getDevicePointer(), &posq.getDevicePointer(), &atomIndexDevice.getDevicePointer(), &numAtoms};
    executeKernel(setChargesKernel, args, numAtoms);
}

bool CudaContext::requestPosqCharges() {
    bool allow = !hasAssignedPosqCharges;
    hasAssignedPosqCharges = true;
    return allow;
}

void CudaContext::addEnergyParameterDerivative(const string& param) {
    // See if this parameter has already been registered.
    
    for (int i = 0; i < energyParamDerivNames.size(); i++)
        if (param == energyParamDerivNames[i])
            return;
    energyParamDerivNames.push_back(param);
}

void CudaContext::flushQueue() {
    cuStreamSynchronize(getCurrentStream());
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
        CHECK_RESULT(cuDeviceGetAttribute(&major, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, thisDevice));
        CHECK_RESULT(cuDeviceGetAttribute(&minor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, thisDevice));
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
