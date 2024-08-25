/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2019 Stanford University and the Authors.      *
 * Portions copyright (c) 2020 Advanced Micro Devices, Inc.                   *
 * Authors: Peter Eastman, Nicholas Curtis                                    *
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
  #error "Windows unsupported for HIP platform"
#endif
#include <cmath>
#include "HipContext.h"
#include "HipArray.h"
#include "HipBondedUtilities.h"
#include "HipEvent.h"
#include "HipIntegrationUtilities.h"
#include "HipKernels.h"
#include "HipKernelSources.h"
#include "HipNonbondedUtilities.h"
#include "HipProgram.h"
#include "openmm/common/ComputeArray.h"
#include "SHA1.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/VirtualSite.h"
#include "HipExpressionUtilities.h"
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
#include <unistd.h>


#define CHECK_RESULT(result) CHECK_RESULT2(result, errorMessage);
#define CHECK_RESULT2(result, prefix) \
    if (result != hipSuccess) { \
        std::stringstream m; \
        m<<prefix<<": "<<getErrorString(result)<<" ("<<result<<")"<<" at "<<__FILE__<<":"<<__LINE__; \
        throw OpenMMException(m.str());\
    }

using namespace OpenMM;
using namespace std;

const int HipContext::ThreadBlockSize = 64;
const int HipContext::TileSize = 32;
bool HipContext::hasInitializedHip = false;


HipContext::HipContext(const System& system, int deviceIndex, bool useBlockingSync, const string& precision, const string& compiler,
        const string& tempDir, const std::string& hostCompiler, HipPlatform::PlatformData& platformData, HipContext* originalContext) : ComputeContext(system), currentStream(0),
        platformData(platformData), contextIsValid(false), hasAssignedPosqCharges(false),
        hasCompilerKernel(false), isHipccAvailable(false), pinnedBuffer(NULL), integration(NULL), expression(NULL), bonded(NULL), nonbonded(NULL),
        supportsHardwareFloatGlobalAtomicAdd(false) {
    // Determine what compiler to use.

    this->compiler = "\""+compiler+"\"";
    if (platformData.context != NULL) {
        try {
            compilerKernel = platformData.context->getPlatform().createKernel(HipCompilerKernel::Name(), *platformData.context);
            hasCompilerKernel = true;
        }
        catch (...) {
            // The runtime compiler plugin isn't available.
        }
    }
    string testCompilerCommand = this->compiler+" --version > /dev/null 2> /dev/null";
    int res = std::system(testCompilerCommand.c_str());
    struct stat info;
    isHipccAvailable = (res == 0 && stat(tempDir.c_str(), &info) == 0);
    if (!hasInitializedHip) {
        CHECK_RESULT2(hipInit(0), "Error initializing HIP");
        hasInitializedHip = true;
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
    this->tempDir = tempDir+"/";
    cacheDir = cacheDir+"/";
    contextIndex = platformData.contexts.size();
    string errorMessage = "Error initializing Context";
    if (originalContext == NULL) {
        isLinkedContext = false;
        int numDevices;
        CHECK_RESULT(hipGetDeviceCount(&numDevices));
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
            CHECK_RESULT(hipDeviceGet(&device, trialDeviceIndex));
            // try setting device
            if (hipSetDevice(device) == hipSuccess) {
                // and set flags
                unsigned int flags = hipDeviceMapHost;
                if (useBlockingSync)
                    flags += hipDeviceScheduleBlockingSync;
                else
                    flags += hipDeviceScheduleSpin;

                if (hipSetDeviceFlags(flags) == hipSuccess) {
                    this->deviceIndex = trialDeviceIndex;
                    break;
                }
            }

        }
        if (this->deviceIndex == -1) {
            if (deviceIndex != -1)
                throw OpenMMException("The requested HIP device could not be loaded");
            else
                throw OpenMMException("No compatible HIP device is available");
        }
    }
    else {
        isLinkedContext = true;
        this->deviceIndex = originalContext->deviceIndex;
        this->device = originalContext->device;
    }

    hipDeviceProp_t props;
    CHECK_RESULT(hipGetDeviceProperties(&props, device));

    // set device properties
    this->simdWidth = props.warpSize;
    this->sharedMemPerBlock = props.sharedMemPerBlock;

    gpuArchitecture = props.gcnArchName;
    // HIP-TODO: find a good value here
    int numThreadBlocksPerComputeUnit = 6;

    // GPUs starting from CDNA1 and RDNA3 support atomic add for floats (global_atomic_add_f32),
    // which can be used in PME. Older GPUs use fixed point charge spreading instead.
    this->supportsHardwareFloatGlobalAtomicAdd = true;
    if (gpuArchitecture.find("gfx900") == 0 ||
        gpuArchitecture.find("gfx906") == 0 ||
        gpuArchitecture.find("gfx10") == 0) {
        this->supportsHardwareFloatGlobalAtomicAdd = false;
    }

    contextIsValid = true;
    if (contextIndex > 0) {
        int canAccess;
        CHECK_RESULT(hipDeviceCanAccessPeer(&canAccess, getDevice(), platformData.contexts[0]->getDevice()));
        if (canAccess) {
            platformData.contexts[0]->setAsCurrent();
            CHECK_RESULT(hipDeviceEnablePeerAccess(getDevice(), 0));
            setAsCurrent();
            CHECK_RESULT(hipDeviceEnablePeerAccess(platformData.contexts[0]->getDevice(), 0));
        }
    }
    numAtoms = system.getNumParticles();
    paddedNumAtoms = TileSize*((numAtoms+TileSize-1)/TileSize);
    numAtomBlocks = (paddedNumAtoms+(TileSize-1))/TileSize;
    CHECK_RESULT(hipDeviceGetAttribute(&multiprocessors, hipDeviceAttributeMultiprocessorCount, device));
    // For RDNA GPUs hipDeviceAttributeMultiprocessorCount means WGP (work-group processors, two compute units), not CUs.
    if (simdWidth == 32)
        multiprocessors *= 2;
    numThreadBlocks = numThreadBlocksPerComputeUnit*multiprocessors;

    compilationDefines["USE_HIP"] = "1";
    if (simdWidth == 32)
        compilationDefines["AMD_RDNA"] = "1";
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

    hipModule_t utilities = createModule(HipKernelSources::vectorOps+HipKernelSources::utilities);
    clearBufferKernel = getKernel(utilities, "clearBuffer");
    clearTwoBuffersKernel = getKernel(utilities, "clearTwoBuffers");
    clearThreeBuffersKernel = getKernel(utilities, "clearThreeBuffers");
    clearFourBuffersKernel = getKernel(utilities, "clearFourBuffers");
    clearFiveBuffersKernel = getKernel(utilities, "clearFiveBuffers");
    clearSixBuffersKernel = getKernel(utilities, "clearSixBuffers");
    reduceEnergyKernel = getKernel(utilities, "reduceEnergy");
    setChargesKernel = getKernel(utilities, "setCharges");

    // Set defines based on the requested precision.

    compilationDefines["SQRT"] = useDoublePrecision ? "sqrt" : "__fsqrt_rn";
    compilationDefines["RSQRT"] = useDoublePrecision ? "rsqrt" : "__frsqrt_rn";
    compilationDefines["RECIP(x)"] = useDoublePrecision ? "(1.0/(x))" : "(1.0f/(x))";
    compilationDefines["EXP"] = useDoublePrecision ? "exp" : "__expf";
    compilationDefines["LOG"] = useDoublePrecision ? "log" : "__logf";
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

    bonded = new HipBondedUtilities(*this);
    nonbonded = new HipNonbondedUtilities(*this);
    integration = new HipIntegrationUtilities(*this, system);
    expression = new HipExpressionUtilities(*this);
}

HipContext::~HipContext() {
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
        hipHostFree(pinnedBuffer);
    if (integration != NULL)
        delete integration;
    if (expression != NULL)
        delete expression;
    if (bonded != NULL)
        delete bonded;
    if (nonbonded != NULL)
        delete nonbonded;
    contextIsValid = false;
}

void HipContext::initialize() {
    hipSetDevice(device);
    string errorMessage = "Error initializing Context";
    int numEnergyBuffers = max(numThreadBlocks*ThreadBlockSize, nonbonded->getNumEnergyBuffers());
    if (useDoublePrecision) {
        energyBuffer.initialize<double>(*this, numEnergyBuffers, "energyBuffer");
        energySum.initialize<double>(*this, 1, "energySum");
        int pinnedBufferSize = max(paddedNumAtoms*4, numEnergyBuffers);
        CHECK_RESULT(hipHostMalloc(&pinnedBuffer, pinnedBufferSize*sizeof(double), 0));
    }
    else if (useMixedPrecision) {
        energyBuffer.initialize<double>(*this, numEnergyBuffers, "energyBuffer");
        energySum.initialize<double>(*this, 1, "energySum");
        int pinnedBufferSize = max(paddedNumAtoms*4, numEnergyBuffers);
        CHECK_RESULT(hipHostMalloc(&pinnedBuffer, pinnedBufferSize*sizeof(double), 0));
    }
    else {
        energyBuffer.initialize<float>(*this, numEnergyBuffers, "energyBuffer");
        energySum.initialize<float>(*this, 1, "energySum");
        int pinnedBufferSize = max(paddedNumAtoms*6, numEnergyBuffers);
        CHECK_RESULT(hipHostMalloc(&pinnedBuffer, pinnedBufferSize*sizeof(float), 0));
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

void HipContext::initializeContexts() {
    getPlatformData().initializeContexts(system);
}

void HipContext::setAsCurrent() {
    if (contextIsValid)
        hipSetDevice(device);
}

hipModule_t HipContext::createModule(const string source) {
    return createModule(source, map<string, string>());
}

hipModule_t HipContext::createModule(const string source, const map<string, string>& defines) {
    const char* saveTempsEnv = getenv("OPENMM_SAVE_TEMPS");
    bool saveTemps = saveTempsEnv != nullptr;
    string bits = intToString(8*sizeof(void*));
    string options = "-ffast-math -munsafe-fp-atomics -Wall";
    if (gpuArchitecture.find("gfx90a") == 0 ||
        gpuArchitecture.find("gfx94") == 0) {
        // HIP-TODO: Remove it when the compiler does a better job
        // Disable SLP vectorization as it may generate unoptimal packed math instructions on
        // >=MI200 (gfx90a, gfx942): more v_mov, higher register usage etc.
        options += " -fno-slp-vectorize";
    }
    if (getMaxThreadBlockSize() < 1024) {
        options += " --gpu-max-threads-per-block=" + std::to_string(getMaxThreadBlockSize());
    }
    stringstream src;
    if (!options.empty())
        src << "// Compilation Options: " << options << endl << endl;
    for (auto& pair : compilationDefines) {
        // Query defines to avoid duplicate variables
        if (defines.find(pair.first) == defines.end()) {
            src << "#define " << pair.first;
            if (!pair.second.empty())
                src << " " << pair.second;
            src << endl;
        }
    }
    if (!compilationDefines.empty())
        src << endl;

    // include the main header for built-in variables (threadIdx etc.) and functions
    src << "#include \"hip/hip_runtime.h\"\n";

    // include the vector types
    src << "#include \"hip/hip_vector_types.h\"\n";
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
    src << HipKernelSources::common << endl;
    for (auto& pair : defines) {
        src << "#define " << pair.first;
        if (!pair.second.empty())
            src << " " << pair.second;
        src << endl;
    }
    if (!defines.empty())
        src << endl;
    src << HipKernelSources::intrinsics << endl;
    src << source << endl;

    // See whether we already have PTX for this kernel cached.

    CSHA1 sha1;
    sha1.Update((const UINT_8*) src.str().c_str(), src.str().size());
    sha1.Final();
    UINT_8 hash[20];
    sha1.GetHash(hash);
    stringstream cacheHash;
    cacheHash.flags(ios::hex);
    for (int i = 0; i < 20; i++)
        cacheHash << setw(2) << setfill('0') << (int) hash[i];
    stringstream cacheFile;
    cacheFile << cacheDir << cacheHash.str() << '_' << gpuArchitecture << '_' << bits;
    hipModule_t module;
    if (hipModuleLoad(&module, cacheFile.str().c_str()) == hipSuccess)
        return module;

    // Select names for the various temporary files.

    stringstream tempFileName;
    if (saveTemps) {
        tempFileName << saveTempsEnv;
        const char* saveTempsPrefixEnv = getenv("OPENMM_SAVE_TEMPS_PREFIX");
        if (saveTempsPrefixEnv) {
            tempFileName << saveTempsPrefixEnv;
        }
        tempFileName << cacheHash.str();
    }
    else {
        tempFileName << tempDir;
        tempFileName << "openmmTempKernel" << this; // Include a pointer to this context as part of the filename to avoid collisions.
        tempFileName << "_" << getpid();
    }
    string inputFile = (tempFileName.str()+".hip.cpp");
    string outputFile = (tempFileName.str()+".hsaco");
    string logFile = (tempFileName.str()+".log");
    int res = 0;

    // If the runtime compiler plugin is available, use it.

    if (hasCompilerKernel) {
        string ptx = compilerKernel.getAs<HipCompilerKernel>().createModule(src.str(), options, *this);

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

            CHECK_RESULT2(hipModuleLoadDataEx(&module, &ptx[0], 0, NULL, NULL), "Error loading HIP module");
            return module;
        }
    }
    else {
        // Write out the source to a temporary file.

        ofstream out(inputFile.c_str());
        out << src.str();
        out.close();
        string command = compiler + " --genco --amdgpu-target=" + gpuArchitecture + " " + options + (saveTemps ? " -save-temps=obj" : "") +" -o \""+outputFile+"\" " + " \""+inputFile+"\" 2> \""+logFile+"\"";
        res = std::system(command.c_str());
    }
    try {
        if (res != 0) {
            // Load the error log.

            stringstream error;
            error << "Error launching HIP compiler: " << res;
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
        hipError_t result = hipModuleLoad(&module, outputFile.c_str());
        if (result != hipSuccess) {
            std::stringstream m;
            m<<"Error loading HIP module: "<<getErrorString(result)<<" ("<<result<<")";
            throw OpenMMException(m.str());
        }
        if (!saveTemps) {
            remove(inputFile.c_str());
            remove(logFile.c_str());
        }
        if (rename(outputFile.c_str(), cacheFile.str().c_str()) != 0 && !saveTemps)
            remove(outputFile.c_str());
        return module;
    }
    catch (...) {
        if (!saveTemps) {
            remove(inputFile.c_str());
            remove(outputFile.c_str());
            remove(logFile.c_str());
        }
        throw;
    }
}

hipFunction_t HipContext::getKernel(hipModule_t& module, const string& name) {
    hipFunction_t function;
    hipError_t result = hipModuleGetFunction(&function, module, name.c_str());
    if (result != hipSuccess) {
        std::stringstream m;
        m<<"Error creating kernel "<<name<<": "<<getErrorString(result)<<" ("<<result<<")";
        throw OpenMMException(m.str());
    }
    return function;
}

hipStream_t HipContext::getCurrentStream() {
    return currentStream;
}

void HipContext::setCurrentStream(hipStream_t stream) {
    currentStream = stream;
}

void HipContext::restoreDefaultStream() {
    setCurrentStream(0);
}

HipArray* HipContext::createArray() {
    return new HipArray();
}

ComputeEvent HipContext::createEvent() {
    return shared_ptr<ComputeEventImpl>(new HipEvent(*this));
}

ComputeProgram HipContext::compileProgram(const std::string source, const std::map<std::string, std::string>& defines) {
    hipModule_t module = createModule(HipKernelSources::vectorOps+source, defines);
    return shared_ptr<ComputeProgramImpl>(new HipProgram(*this, module));
}

HipArray& HipContext::unwrap(ArrayInterface& array) const {
    HipArray* cuarray;
    ComputeArray* wrapper = dynamic_cast<ComputeArray*>(&array);
    if (wrapper != NULL)
        cuarray = dynamic_cast<HipArray*>(&wrapper->getArray());
    else
        cuarray = dynamic_cast<HipArray*>(&array);
    if (cuarray == NULL)
        throw OpenMMException("Array argument is not an HipArray");
    return *cuarray;
}

std::string HipContext::getErrorString(hipError_t result) {
    return string(hipGetErrorName(result));
}

void HipContext::executeKernel(hipFunction_t kernel, void** arguments, int threads, int blockSize, unsigned int sharedSize) {
    if (blockSize == -1)
        blockSize = ThreadBlockSize;
    int gridSize = std::min((threads+blockSize-1)/blockSize, numThreadBlocks);
    hipError_t result = hipModuleLaunchKernel(kernel, gridSize, 1, 1, blockSize, 1, 1, sharedSize, currentStream, arguments, NULL);
    if (result != hipSuccess) {
        stringstream str;
        str<<"Error invoking kernel: "<<getErrorString(result)<<" ("<<result<<")";
        throw OpenMMException(str.str());
    }
}

void HipContext::executeKernelFlat(hipFunction_t kernel, void** arguments, int threads, int blockSize, unsigned int sharedSize) {
    if (blockSize == -1)
        blockSize = ThreadBlockSize;
    int gridSize = (threads+blockSize-1)/blockSize;
    hipError_t result = hipModuleLaunchKernel(kernel, gridSize, 1, 1, blockSize, 1, 1, sharedSize, currentStream, arguments, NULL);
    if (result != hipSuccess) {
        stringstream str;
        str<<"Error invoking kernel: "<<getErrorString(result)<<" ("<<result<<")";
        throw OpenMMException(str.str());
    }
}

int HipContext::computeThreadBlockSize(double memory) const {
    int maxShared = this->sharedMemPerBlock;
    int max = (int) (maxShared/memory);
    if (max < HipContext::ThreadBlockSize) {
        throw OpenMMException("Too much shared memory requested!");
    }
    int threads = this->simdWidth;
    while (threads+this->simdWidth < max)
        threads += this->simdWidth;
    return threads;
}

void HipContext::clearBuffer(ArrayInterface& array) {
    clearBuffer(unwrap(array).getDevicePointer(), array.getSize()*array.getElementSize());
}

void HipContext::clearBuffer(hipDeviceptr_t memory, int size) {
    int words = size/4;
    void* args[] = {&memory, &words};
    executeKernel(clearBufferKernel, args, words, 4 * this->simdWidth);
}

void HipContext::addAutoclearBuffer(ArrayInterface& array) {
    addAutoclearBuffer(unwrap(array).getDevicePointer(), array.getSize()*array.getElementSize());
}

void HipContext::addAutoclearBuffer(hipDeviceptr_t memory, int size) {
    autoclearBuffers.push_back(memory);
    autoclearBufferSizes.push_back(size/4);
}

void HipContext::clearAutoclearBuffers() {

    int preferredTBSize = this->simdWidth * 4;
    int base = 0;
    int total = autoclearBufferSizes.size();
    while (total-base >= 6) {
        void* args[] = {&autoclearBuffers[base], &autoclearBufferSizes[base],
                        &autoclearBuffers[base+1], &autoclearBufferSizes[base+1],
                        &autoclearBuffers[base+2], &autoclearBufferSizes[base+2],
                        &autoclearBuffers[base+3], &autoclearBufferSizes[base+3],
                        &autoclearBuffers[base+4], &autoclearBufferSizes[base+4],
                        &autoclearBuffers[base+5], &autoclearBufferSizes[base+5]};
        executeKernel(clearSixBuffersKernel, args, max(max(max(max(max(autoclearBufferSizes[base], autoclearBufferSizes[base+1]), autoclearBufferSizes[base+2]), autoclearBufferSizes[base+3]), autoclearBufferSizes[base+4]), autoclearBufferSizes[base+5]), preferredTBSize);
        base += 6;
    }
    if (total-base == 5) {
        void* args[] = {&autoclearBuffers[base], &autoclearBufferSizes[base],
                        &autoclearBuffers[base+1], &autoclearBufferSizes[base+1],
                        &autoclearBuffers[base+2], &autoclearBufferSizes[base+2],
                        &autoclearBuffers[base+3], &autoclearBufferSizes[base+3],
                        &autoclearBuffers[base+4], &autoclearBufferSizes[base+4]};
        executeKernel(clearFiveBuffersKernel, args, max(max(max(max(autoclearBufferSizes[base], autoclearBufferSizes[base+1]), autoclearBufferSizes[base+2]), autoclearBufferSizes[base+3]), autoclearBufferSizes[base+4]), preferredTBSize);
    }
    else if (total-base == 4) {
        void* args[] = {&autoclearBuffers[base], &autoclearBufferSizes[base],
                        &autoclearBuffers[base+1], &autoclearBufferSizes[base+1],
                        &autoclearBuffers[base+2], &autoclearBufferSizes[base+2],
                        &autoclearBuffers[base+3], &autoclearBufferSizes[base+3]};
        executeKernel(clearFourBuffersKernel, args, max(max(max(autoclearBufferSizes[base], autoclearBufferSizes[base+1]), autoclearBufferSizes[base+2]), autoclearBufferSizes[base+3]), preferredTBSize);
    }
    else if (total-base == 3) {
        void* args[] = {&autoclearBuffers[base], &autoclearBufferSizes[base],
                        &autoclearBuffers[base+1], &autoclearBufferSizes[base+1],
                        &autoclearBuffers[base+2], &autoclearBufferSizes[base+2]};
        executeKernel(clearThreeBuffersKernel, args, max(max(autoclearBufferSizes[base], autoclearBufferSizes[base+1]), autoclearBufferSizes[base+2]), preferredTBSize);
    }
    else if (total-base == 2) {
        void* args[] = {&autoclearBuffers[base], &autoclearBufferSizes[base],
                        &autoclearBuffers[base+1], &autoclearBufferSizes[base+1]};
        executeKernel(clearTwoBuffersKernel, args, max(autoclearBufferSizes[base], autoclearBufferSizes[base+1]), preferredTBSize);
    }
    else if (total-base == 1) {
        clearBuffer(autoclearBuffers[base], autoclearBufferSizes[base]*4);
    }
}

double HipContext::reduceEnergy() {
    int bufferSize = energyBuffer.getSize();
    int workGroupSize = getMaxThreadBlockSize();
    void* args[] = {&energyBuffer.getDevicePointer(), &energySum.getDevicePointer(), &bufferSize, &workGroupSize};
    executeKernel(reduceEnergyKernel, args, workGroupSize, workGroupSize, workGroupSize*energyBuffer.getElementSize());
    energySum.download(pinnedBuffer);
    if (getUseDoublePrecision() || getUseMixedPrecision())
        return *((double*) pinnedBuffer);
    else
        return *((float*) pinnedBuffer);
}

void HipContext::setCharges(const vector<double>& charges) {
    if (!chargeBuffer.isInitialized())
        chargeBuffer.initialize(*this, numAtoms, useDoublePrecision ? sizeof(double) : sizeof(float), "chargeBuffer");
    vector<double> c(numAtoms);
    for (int i = 0; i < numAtoms; i++)
        c[i] = charges[i];
    chargeBuffer.upload(c, true);
    void* args[] = {&chargeBuffer.getDevicePointer(), &posq.getDevicePointer(), &atomIndexDevice.getDevicePointer(), &numAtoms};
    executeKernel(setChargesKernel, args, numAtoms);
}

bool HipContext::requestPosqCharges() {
    bool allow = !hasAssignedPosqCharges;
    hasAssignedPosqCharges = true;
    return allow;
}

void HipContext::addEnergyParameterDerivative(const string& param) {
    // See if this parameter has already been registered.

    for (int i = 0; i < energyParamDerivNames.size(); i++)
        if (param == energyParamDerivNames[i])
            return;
    energyParamDerivNames.push_back(param);
}

void HipContext::flushQueue() {
    hipStreamSynchronize(getCurrentStream());
}

vector<int> HipContext::getDevicePrecedence() {
    int numDevices;
    hipDeviceProp_t thisDevice;
    string errorMessage = "Error initializing Context";
    vector<pair<int, int> > devices;

    CHECK_RESULT(hipGetDeviceCount(&numDevices));
    for (int i = 0; i < numDevices; i++) {
        CHECK_RESULT(hipGetDeviceProperties(&thisDevice, i));
        int clock, multiprocessors, speed;
        clock = thisDevice.clockRate;
        multiprocessors = thisDevice.multiProcessorCount;
        speed = clock*multiprocessors;
        devices.push_back(std::make_pair(speed, -i));
    }

    // sort first by speed (higher is better), and finally device index (lower is better)
    std::sort(devices.begin(), devices.end());
    std::reverse(devices.begin(), devices.end());

    vector<int> precedence;
    for (int i = 0; i < static_cast<int>(devices.size()); i++) {
        precedence.push_back(-devices[i].second);
    }

    return precedence;
}
