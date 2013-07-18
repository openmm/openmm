/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2012 Stanford University and the Authors.      *
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
#include "OpenCLContext.h"
#include "OpenCLArray.h"
#include "OpenCLBondedUtilities.h"
#include "OpenCLForceInfo.h"
#include "OpenCLIntegrationUtilities.h"
#include "OpenCLKernelSources.h"
#include "OpenCLNonbondedUtilities.h"
#include "hilbert.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/VirtualSite.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <typeinfo>

using namespace OpenMM;
using namespace std;

#ifndef CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV
  #define CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV 0x4000
#endif
#ifndef CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV
  #define CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV 0x4001
#endif

const int OpenCLContext::ThreadBlockSize = 64;
const int OpenCLContext::TileSize = 32;

static void CL_CALLBACK errorCallback(const char* errinfo, const void* private_info, size_t cb, void* user_data) {
    string skip = "OpenCL Build Warning : Compiler build log:";
    if (strncmp(errinfo, skip.c_str(), skip.length()) == 0)
        return; // OS X Lion insists on calling this for every build warning, even though they aren't errors.
    std::cerr << "OpenCL internal error: " << errinfo << std::endl;
}

OpenCLContext::OpenCLContext(const System& system, int platformIndex, int deviceIndex, const string& precision, OpenCLPlatform::PlatformData& platformData) :
        system(system), time(0.0), platformData(platformData), stepCount(0), computeForceCount(0), stepsSinceReorder(99999), atomsWereReordered(false), posq(NULL),
        posqCorrection(NULL), velm(NULL), forceBuffers(NULL), longForceBuffer(NULL), energyBuffer(NULL), atomIndexDevice(NULL), integration(NULL),
        expression(NULL), bonded(NULL), nonbonded(NULL), thread(NULL) {
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
        throw OpenMMException("Illegal value for OpenCLPrecision: "+precision);
    try {
        contextIndex = platformData.contexts.size();
        std::vector<cl::Platform> platforms;
        cl::Platform::get(&platforms);
        if (platformIndex < 0 || platformIndex >= (int) platforms.size())
            throw OpenMMException("Illegal value for OpenCL platform index");
        string platformVendor = platforms[platformIndex].getInfo<CL_PLATFORM_VENDOR>();
        vector<cl::Device> devices;
        platforms[platformIndex].getDevices(CL_DEVICE_TYPE_ALL, &devices);
        const int minThreadBlockSize = 32;
        if (deviceIndex < 0 || deviceIndex >= (int) devices.size()) {
            // Try to figure out which device is the fastest.

            int bestSpeed = -1;
            for (int i = 0; i < (int) devices.size(); i++) {
                if (platformVendor == "Apple" && devices[i].getInfo<CL_DEVICE_VENDOR>() == "AMD")
                    continue; // Don't use AMD GPUs on OS X due to serious bugs.
                int maxSize = devices[i].getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>()[0];
                int processingElementsPerComputeUnit = 8;
                if (devices[i].getInfo<CL_DEVICE_TYPE>() != CL_DEVICE_TYPE_GPU) {
                    processingElementsPerComputeUnit = 1;
                }
                else if (devices[i].getInfo<CL_DEVICE_EXTENSIONS>().find("cl_nv_device_attribute_query") != string::npos) {
                    cl_uint computeCapabilityMajor;
                    clGetDeviceInfo(devices[i](), CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV, sizeof(cl_uint), &computeCapabilityMajor, NULL);
                    processingElementsPerComputeUnit = (computeCapabilityMajor < 2 ? 8 : 32);
                }
                else if (devices[i].getInfo<CL_DEVICE_EXTENSIONS>().find("cl_amd_device_attribute_query") != string::npos) {
                    // This attribute does not ensure that all queries are supported by the runtime (it may be an older runtime,
                    // or the CPU device) so still have to check for errors.
                    try {
#ifdef CL_DEVICE_SIMD_WIDTH_AMD
                        processingElementsPerComputeUnit =
                            // AMD GPUs either have a single VLIW SIMD or multiple scalar SIMDs.
                            // The SIMD width is the number of threads the SIMD executes per cycle.
                            // This will be less than the wavefront width since it takes several
                            // cycles to execute the full wavefront.
                            // The SIMD instruction width is the VLIW instruction width (or 1 for scalar),
                            // this is the number of ALUs that can be executing per instruction per thread. 
                            devices[i].getInfo<CL_DEVICE_SIMD_PER_COMPUTE_UNIT_AMD>() *
                            devices[i].getInfo<CL_DEVICE_SIMD_WIDTH_AMD>() *
                            devices[i].getInfo<CL_DEVICE_SIMD_INSTRUCTION_WIDTH_AMD>();
                        // Just in case any of the queries return 0.
                        if (processingElementsPerComputeUnit <= 0)
                            processingElementsPerComputeUnit = 1;
#endif
                    }
                    catch (cl::Error err) {
                        // Runtime does not support the queries so use default.
                    }
                }
                int speed = devices[i].getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>()*processingElementsPerComputeUnit*devices[i].getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>();
                if (maxSize >= minThreadBlockSize && speed > bestSpeed) {
                    deviceIndex = i;
                    bestSpeed = speed;
                }
            }
        }
        if (deviceIndex == -1)
            throw OpenMMException("No compatible OpenCL device is available");
        device = devices[deviceIndex];
        this->deviceIndex = deviceIndex;
        if (device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>() < minThreadBlockSize)
            throw OpenMMException("The specified OpenCL device is not compatible with OpenMM");
        compilationDefines["WORK_GROUP_SIZE"] = intToString(ThreadBlockSize);
        if (platformVendor.size() >= 5 && platformVendor.substr(0, 5) == "Intel")
            defaultOptimizationOptions = "";
        else
            defaultOptimizationOptions = "-cl-fast-relaxed-math";
        supports64BitGlobalAtomics = (device.getInfo<CL_DEVICE_EXTENSIONS>().find("cl_khr_int64_base_atomics") != string::npos);
        supportsDoublePrecision = (device.getInfo<CL_DEVICE_EXTENSIONS>().find("cl_khr_fp64") != string::npos);
        if ((useDoublePrecision || useMixedPrecision) && !supportsDoublePrecision)
            throw OpenMMException("This device does not support double precision");
        string vendor = device.getInfo<CL_DEVICE_VENDOR>();
        int numThreadBlocksPerComputeUnit = 6;
        if (vendor.size() >= 6 && vendor.substr(0, 6) == "NVIDIA") {
            compilationDefines["WARPS_ARE_ATOMIC"] = "";
            simdWidth = 32;
            if (device.getInfo<CL_DEVICE_EXTENSIONS>().find("cl_nv_device_attribute_query") != string::npos) {
                // Compute level 1.2 and later Nvidia GPUs support 64 bit atomics, even though they don't list the
                // proper extension as supported.  We only use them on compute level 2.0 or later, since they're very
                // slow on earlier GPUs.

                cl_uint computeCapabilityMajor;
                clGetDeviceInfo(device(), CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV, sizeof(cl_uint), &computeCapabilityMajor, NULL);
                if (computeCapabilityMajor > 1)
                    supports64BitGlobalAtomics = true;
            }
        }
        else if (vendor.size() >= 28 && vendor.substr(0, 28) == "Advanced Micro Devices, Inc.") {
            if (device.getInfo<CL_DEVICE_TYPE>() != CL_DEVICE_TYPE_GPU) {
                /// \todo Is 6 a good value for the OpenCL CPU device?
                // numThreadBlocksPerComputeUnit = ?;
                simdWidth = 1;
            }
            else {
                bool amdPostSdk2_4 = false;
                // Default to 1 which will use the default kernels.
                simdWidth = 1;
                if (device.getInfo<CL_DEVICE_EXTENSIONS>().find("cl_amd_device_attribute_query") != string::npos) {
                    // This attribute does not ensure that all queries are supported by the runtime so still have to
                    // check for errors.
                    try {
#ifdef CL_DEVICE_SIMD_PER_COMPUTE_UNIT_AMD
                        // Must catch cl:Error as will fail if runtime does not support queries.

                        cl_uint simdPerComputeUnit = device.getInfo<CL_DEVICE_SIMD_PER_COMPUTE_UNIT_AMD>();
                        simdWidth = device.getInfo<CL_DEVICE_WAVEFRONT_WIDTH_AMD>();

                        // If the GPU has multiple SIMDs per compute unit then it is uses the scalar instruction
                        // set instead of the VLIW instruction set. It therefore needs more thread blocks per
                        // compute unit to hide memory latency.
                        if (simdPerComputeUnit > 1)
                            numThreadBlocksPerComputeUnit = 4 * simdPerComputeUnit;

                        // If the queries are supported then must be newer than SDK 2.4.
                        amdPostSdk2_4 = true;
#endif
                    }
                    catch (cl::Error err) {
                        // Runtime does not support the query so is unlikely to be the newer scalar GPU.
                        // Stay with the default simdWidth and numThreadBlocksPerComputeUnit.
                    }
                }
                // AMD APP SDK 2.4 has a performance problem with atomics. Enable the work around. This is fixed after SDK 2.4.
                if (!amdPostSdk2_4)
                    compilationDefines["AMD_ATOMIC_WORK_AROUND"] = "";
            }
        }
        else
            simdWidth = 1;
        if (platformVendor == "Apple" && vendor == "AMD")
            compilationDefines["MAC_AMD_WORKAROUND"] = "";
        if (supports64BitGlobalAtomics)
            compilationDefines["SUPPORTS_64_BIT_ATOMICS"] = "";
        if (supportsDoublePrecision)
            compilationDefines["SUPPORTS_DOUBLE_PRECISION"] = "";
        if (simdWidth >= 32)
            compilationDefines["SYNC_WARPS"] = "";
        else
            compilationDefines["SYNC_WARPS"] = "barrier(CLK_LOCAL_MEM_FENCE)";
        vector<cl::Device> contextDevices;
        contextDevices.push_back(device);
        cl_context_properties cprops[] = {CL_CONTEXT_PLATFORM, (cl_context_properties) platforms[platformIndex](), 0};
        context = cl::Context(contextDevices, cprops, errorCallback);
        queue = cl::CommandQueue(context, device);
        numAtoms = system.getNumParticles();
        paddedNumAtoms = TileSize*((numAtoms+TileSize-1)/TileSize);
        numAtomBlocks = (paddedNumAtoms+(TileSize-1))/TileSize;
        numThreadBlocks = numThreadBlocksPerComputeUnit*device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
        bonded = new OpenCLBondedUtilities(*this);
        nonbonded = new OpenCLNonbondedUtilities(*this);
        if (useDoublePrecision) {
            posq = OpenCLArray::create<mm_double4>(*this, paddedNumAtoms, "posq");
            velm = OpenCLArray::create<mm_double4>(*this, paddedNumAtoms, "velm");
            compilationDefines["USE_DOUBLE_PRECISION"] = "1";
            compilationDefines["convert_real4"] = "convert_double4";
            compilationDefines["convert_mixed4"] = "convert_double4";
        }
        else if (useMixedPrecision) {
            posq = OpenCLArray::create<mm_float4>(*this, paddedNumAtoms, "posq");
            posqCorrection = OpenCLArray::create<mm_float4>(*this, paddedNumAtoms, "posq");
            velm = OpenCLArray::create<mm_double4>(*this, paddedNumAtoms, "velm");
            compilationDefines["USE_MIXED_PRECISION"] = "1";
            compilationDefines["convert_real4"] = "convert_float4";
            compilationDefines["convert_mixed4"] = "convert_double4";
        }
        else {
            posq = OpenCLArray::create<mm_float4>(*this, paddedNumAtoms, "posq");
            velm = OpenCLArray::create<mm_float4>(*this, paddedNumAtoms, "velm");
            compilationDefines["convert_real4"] = "convert_float4";
            compilationDefines["convert_mixed4"] = "convert_float4";
        }
        posCellOffsets.resize(paddedNumAtoms, mm_int4(0, 0, 0, 0));
    }
    catch (cl::Error err) {
        std::stringstream str;
        str<<"Error initializing context: "<<err.what()<<" ("<<err.err()<<")";
        throw OpenMMException(str.str());
    }

    // Create utility kernels that are used in multiple places.

    cl::Program utilities = createProgram(OpenCLKernelSources::utilities);
    clearBufferKernel = cl::Kernel(utilities, "clearBuffer");
    clearTwoBuffersKernel = cl::Kernel(utilities, "clearTwoBuffers");
    clearThreeBuffersKernel = cl::Kernel(utilities, "clearThreeBuffers");
    clearFourBuffersKernel = cl::Kernel(utilities, "clearFourBuffers");
    clearFiveBuffersKernel = cl::Kernel(utilities, "clearFiveBuffers");
    clearSixBuffersKernel = cl::Kernel(utilities, "clearSixBuffers");
    reduceReal4Kernel = cl::Kernel(utilities, "reduceReal4Buffer");
    if (supports64BitGlobalAtomics)
        reduceForcesKernel = cl::Kernel(utilities, "reduceForces");

    // Decide whether native_sqrt(), native_rsqrt(), and native_recip() are sufficiently accurate to use.

    if (!useDoublePrecision) {
        cl::Kernel accuracyKernel(utilities, "determineNativeAccuracy");
        OpenCLArray valuesArray(*this, 20, sizeof(mm_float8), "values");
        vector<mm_float8> values(valuesArray.getSize());
        float nextValue = 1e-4f;
        for (int i = 0; i < (int) values.size(); ++i) {
            values[i].s0 = nextValue;
            nextValue *= (float) M_PI;
        }
        valuesArray.upload(values);
        accuracyKernel.setArg<cl::Buffer>(0, valuesArray.getDeviceBuffer());
        accuracyKernel.setArg<cl_int>(1, values.size());
        executeKernel(accuracyKernel, values.size());
        valuesArray.download(values);
        double maxSqrtError = 0.0, maxRsqrtError = 0.0, maxRecipError = 0.0, maxExpError = 0.0, maxLogError = 0.0;
        for (int i = 0; i < (int) values.size(); ++i) {
            double v = values[i].s0;
            double correctSqrt = sqrt(v);
            maxSqrtError = max(maxSqrtError, fabs(correctSqrt-values[i].s1)/correctSqrt);
            maxRsqrtError = max(maxRsqrtError, fabs(1.0/correctSqrt-values[i].s2)*correctSqrt);
            maxRecipError = max(maxRecipError, fabs(1.0/v-values[i].s3)/values[i].s3);
            maxExpError = max(maxExpError, fabs(exp(v)-values[i].s4)/values[i].s4);
            maxLogError = max(maxLogError, fabs(log(v)-values[i].s5)/values[i].s5);
        }
        compilationDefines["SQRT"] = (maxSqrtError < 1e-6) ? "native_sqrt" : "sqrt";
        compilationDefines["RSQRT"] = (maxRsqrtError < 1e-6) ? "native_rsqrt" : "rsqrt";
        compilationDefines["RECIP"] = (maxRecipError < 1e-6) ? "native_recip" : "1.0f/";
        compilationDefines["EXP"] = (maxExpError < 1e-6) ? "native_exp" : "exp";
        compilationDefines["LOG"] = (maxLogError < 1e-6) ? "native_log" : "log";
    }
    else {
        compilationDefines["SQRT"] = "sqrt";
        compilationDefines["RSQRT"] = "rsqrt";
        compilationDefines["RECIP"] = "1.0/";
        compilationDefines["EXP"] = "exp";
        compilationDefines["LOG"] = "log";
    }
    
    // Create the work thread used for parallelization when running on multiple devices.
    
    thread = new WorkThread();
    
    // Create utilities objects.
    
    integration = new OpenCLIntegrationUtilities(*this, system);
    expression = new OpenCLExpressionUtilities(*this);
}

OpenCLContext::~OpenCLContext() {
    for (int i = 0; i < (int) forces.size(); i++)
        delete forces[i];
    for (int i = 0; i < (int) reorderListeners.size(); i++)
        delete reorderListeners[i];
    for (int i = 0; i < (int) preComputations.size(); i++)
        delete preComputations[i];
    for (int i = 0; i < (int) postComputations.size(); i++)
        delete postComputations[i];
    if (pinnedBuffer != NULL)
        delete pinnedBuffer;
    if (posq != NULL)
        delete posq;
    if (posqCorrection != NULL)
        delete posqCorrection;
    if (velm != NULL)
        delete velm;
    if (force != NULL)
        delete force;
    if (forceBuffers != NULL)
        delete forceBuffers;
    if (longForceBuffer != NULL)
        delete longForceBuffer;
    if (energyBuffer != NULL)
        delete energyBuffer;
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
}

void OpenCLContext::initialize() {
    bonded->initialize(system);
    numForceBuffers = platformData.contexts.size();
    numForceBuffers = std::max(numForceBuffers, bonded->getNumForceBuffers());
    for (int i = 0; i < (int) forces.size(); i++)
        numForceBuffers = std::max(numForceBuffers, forces[i]->getRequiredForceBuffers());
    if (useDoublePrecision) {
        forceBuffers = OpenCLArray::create<mm_double4>(*this, paddedNumAtoms*numForceBuffers, "forceBuffers");
        force = OpenCLArray::create<mm_double4>(*this, &forceBuffers->getDeviceBuffer(), paddedNumAtoms, "force");
        energyBuffer = OpenCLArray::create<cl_double>(*this, max(numThreadBlocks*ThreadBlockSize, nonbonded->getNumEnergyBuffers()), "energyBuffer");
    }
    else {
        forceBuffers = OpenCLArray::create<mm_float4>(*this, paddedNumAtoms*numForceBuffers, "forceBuffers");
        force = OpenCLArray::create<mm_float4>(*this, &forceBuffers->getDeviceBuffer(), paddedNumAtoms, "force");
        energyBuffer = OpenCLArray::create<cl_float>(*this, max(numThreadBlocks*ThreadBlockSize, nonbonded->getNumEnergyBuffers()), "energyBuffer");
    }
    if (supports64BitGlobalAtomics) {
        longForceBuffer = OpenCLArray::create<cl_long>(*this, 3*paddedNumAtoms, "longForceBuffer");
        reduceForcesKernel.setArg<cl::Buffer>(0, longForceBuffer->getDeviceBuffer());
        reduceForcesKernel.setArg<cl::Buffer>(1, forceBuffers->getDeviceBuffer());
        reduceForcesKernel.setArg<cl_int>(2, paddedNumAtoms);
        reduceForcesKernel.setArg<cl_int>(3, numForceBuffers);
        addAutoclearBuffer(*longForceBuffer);
    }
    addAutoclearBuffer(*forceBuffers);
    addAutoclearBuffer(*energyBuffer);
    int bufferBytes = max(velm->getSize()*velm->getElementSize(), energyBuffer->getSize()*energyBuffer->getElementSize());
    pinnedBuffer = new cl::Buffer(context, CL_MEM_ALLOC_HOST_PTR, bufferBytes);
    pinnedMemory = queue.enqueueMapBuffer(*pinnedBuffer, CL_TRUE, CL_MAP_READ | CL_MAP_WRITE, 0, bufferBytes);
    for (int i = 0; i < numAtoms; i++) {
        double mass = system.getParticleMass(i);
        if (useDoublePrecision || useMixedPrecision)
            ((mm_double4*) pinnedMemory)[i] = mm_double4(0.0, 0.0, 0.0, mass == 0.0 ? 0.0 : 1.0/mass);
        else
            ((mm_float4*) pinnedMemory)[i] = mm_float4(0.0f, 0.0f, 0.0f, mass == 0.0 ? 0.0f : (cl_float) (1.0/mass));
    }
    velm->upload(pinnedMemory);
    atomIndexDevice = OpenCLArray::create<cl_int>(*this, paddedNumAtoms, "atomIndexDevice");
    atomIndex.resize(paddedNumAtoms);
    for (int i = 0; i < paddedNumAtoms; ++i)
        atomIndex[i] = i;
    atomIndexDevice->upload(atomIndex);
    findMoleculeGroups();
    nonbonded->initialize(system);
}

void OpenCLContext::addForce(OpenCLForceInfo* force) {
    forces.push_back(force);
}

string OpenCLContext::replaceStrings(const string& input, const std::map<std::string, std::string>& replacements) const {
    string result = input;
    for (map<string, string>::const_iterator iter = replacements.begin(); iter != replacements.end(); iter++) {
        int index = -1;
        do {
            index = result.find(iter->first);
            if (index != result.npos)
                result.replace(index, iter->first.size(), iter->second);
        } while (index != result.npos);
    }
    return result;
}

cl::Program OpenCLContext::createProgram(const string source, const char* optimizationFlags) {
    return createProgram(source, map<string, string>(), optimizationFlags);
}

cl::Program OpenCLContext::createProgram(const string source, const map<string, string>& defines, const char* optimizationFlags) {
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
    if (supportsDoublePrecision)
        src << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n";
    if (useDoublePrecision) {
        src << "typedef double real;\n";
        src << "typedef double2 real2;\n";
        src << "typedef double4 real4;\n";
    }
    else {
        src << "typedef float real;\n";
        src << "typedef float2 real2;\n";
        src << "typedef float4 real4;\n";
    }
    if (useDoublePrecision || useMixedPrecision) {
        src << "typedef double mixed;\n";
        src << "typedef double2 mixed2;\n";
        src << "typedef double4 mixed4;\n";
    }
    else {
        src << "typedef float mixed;\n";
        src << "typedef float2 mixed2;\n";
        src << "typedef float4 mixed4;\n";
    }
    for (map<string, string>::const_iterator iter = defines.begin(); iter != defines.end(); ++iter) {
        src << "#define " << iter->first;
        if (!iter->second.empty())
            src << " " << iter->second;
        src << endl;
    }
    if (!defines.empty())
        src << endl;
    src << source << endl;
    // Get length before using c_str() to avoid length() call invalidating the c_str() value.
    string src_string = src.str();
    ::size_t src_length = src_string.length();
    cl::Program::Sources sources(1, make_pair(src_string.c_str(), src_length));
    cl::Program program(context, sources);
    try {
        program.build(vector<cl::Device>(1, device), options.c_str());
    } catch (cl::Error err) {
        throw OpenMMException("Error compiling kernel: "+program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device));
    }
    return program;
}

string OpenCLContext::doubleToString(double value) {
    stringstream s;
    s.precision(useDoublePrecision ? 16 : 8);
    s << scientific << value;
    if (!useDoublePrecision)
        s << "f";
    return s.str();
}

string OpenCLContext::intToString(int value) {
    stringstream s;
    s << value;
    return s.str();
}

void OpenCLContext::executeKernel(cl::Kernel& kernel, int workUnits, int blockSize) {
    if (blockSize == -1)
        blockSize = ThreadBlockSize;
    int size = std::min((workUnits+blockSize-1)/blockSize, numThreadBlocks)*blockSize;
    try {
        queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(size), cl::NDRange(blockSize));
    }
    catch (cl::Error err) {
        stringstream str;
        str<<"Error invoking kernel "<<kernel.getInfo<CL_KERNEL_FUNCTION_NAME>()<<": "<<err.what()<<" ("<<err.err()<<")";
        throw OpenMMException(str.str());
    }
}

void OpenCLContext::clearBuffer(OpenCLArray& array) {
    clearBuffer(array.getDeviceBuffer(), array.getSize()*array.getElementSize());
}

void OpenCLContext::clearBuffer(cl::Memory& memory, int size) {
    int words = size/4;
    clearBufferKernel.setArg<cl::Memory>(0, memory);
    clearBufferKernel.setArg<cl_int>(1, words);
    executeKernel(clearBufferKernel, words, 128);
}

void OpenCLContext::addAutoclearBuffer(OpenCLArray& array) {
    addAutoclearBuffer(array.getDeviceBuffer(), array.getSize()*array.getElementSize());
}

void OpenCLContext::addAutoclearBuffer(cl::Memory& memory, int size) {
    autoclearBuffers.push_back(&memory);
    autoclearBufferSizes.push_back(size/4);
}

void OpenCLContext::clearAutoclearBuffers() {
    int base = 0;
    int total = autoclearBufferSizes.size();
    while (total-base >= 6) {
        clearSixBuffersKernel.setArg<cl::Memory>(0, *autoclearBuffers[base]);
        clearSixBuffersKernel.setArg<cl_int>(1, autoclearBufferSizes[base]);
        clearSixBuffersKernel.setArg<cl::Memory>(2, *autoclearBuffers[base+1]);
        clearSixBuffersKernel.setArg<cl_int>(3, autoclearBufferSizes[base+1]);
        clearSixBuffersKernel.setArg<cl::Memory>(4, *autoclearBuffers[base+2]);
        clearSixBuffersKernel.setArg<cl_int>(5, autoclearBufferSizes[base+2]);
        clearSixBuffersKernel.setArg<cl::Memory>(6, *autoclearBuffers[base+3]);
        clearSixBuffersKernel.setArg<cl_int>(7, autoclearBufferSizes[base+3]);
        clearSixBuffersKernel.setArg<cl::Memory>(8, *autoclearBuffers[base+4]);
        clearSixBuffersKernel.setArg<cl_int>(9, autoclearBufferSizes[base+4]);
        clearSixBuffersKernel.setArg<cl::Memory>(10, *autoclearBuffers[base+5]);
        clearSixBuffersKernel.setArg<cl_int>(11, autoclearBufferSizes[base+5]);
        executeKernel(clearSixBuffersKernel, max(max(max(max(max(autoclearBufferSizes[base], autoclearBufferSizes[base+1]), autoclearBufferSizes[base+2]), autoclearBufferSizes[base+3]), autoclearBufferSizes[base+4]), autoclearBufferSizes[base+5]), 128);
        base += 6;
    }
    if (total-base == 5) {
        clearFiveBuffersKernel.setArg<cl::Memory>(0, *autoclearBuffers[base]);
        clearFiveBuffersKernel.setArg<cl_int>(1, autoclearBufferSizes[base]);
        clearFiveBuffersKernel.setArg<cl::Memory>(2, *autoclearBuffers[base+1]);
        clearFiveBuffersKernel.setArg<cl_int>(3, autoclearBufferSizes[base+1]);
        clearFiveBuffersKernel.setArg<cl::Memory>(4, *autoclearBuffers[base+2]);
        clearFiveBuffersKernel.setArg<cl_int>(5, autoclearBufferSizes[base+2]);
        clearFiveBuffersKernel.setArg<cl::Memory>(6, *autoclearBuffers[base+3]);
        clearFiveBuffersKernel.setArg<cl_int>(7, autoclearBufferSizes[base+3]);
        clearFiveBuffersKernel.setArg<cl::Memory>(8, *autoclearBuffers[base+4]);
        clearFiveBuffersKernel.setArg<cl_int>(9, autoclearBufferSizes[base+4]);
        executeKernel(clearFiveBuffersKernel, max(max(max(max(autoclearBufferSizes[base], autoclearBufferSizes[base+1]), autoclearBufferSizes[base+2]), autoclearBufferSizes[base+3]), autoclearBufferSizes[base+4]), 128);
    }
    else if (total-base == 4) {
        clearFourBuffersKernel.setArg<cl::Memory>(0, *autoclearBuffers[base]);
        clearFourBuffersKernel.setArg<cl_int>(1, autoclearBufferSizes[base]);
        clearFourBuffersKernel.setArg<cl::Memory>(2, *autoclearBuffers[base+1]);
        clearFourBuffersKernel.setArg<cl_int>(3, autoclearBufferSizes[base+1]);
        clearFourBuffersKernel.setArg<cl::Memory>(4, *autoclearBuffers[base+2]);
        clearFourBuffersKernel.setArg<cl_int>(5, autoclearBufferSizes[base+2]);
        clearFourBuffersKernel.setArg<cl::Memory>(6, *autoclearBuffers[base+3]);
        clearFourBuffersKernel.setArg<cl_int>(7, autoclearBufferSizes[base+3]);
        executeKernel(clearFourBuffersKernel, max(max(max(autoclearBufferSizes[base], autoclearBufferSizes[base+1]), autoclearBufferSizes[base+2]), autoclearBufferSizes[base+3]), 128);
    }
    else if (total-base == 3) {
        clearThreeBuffersKernel.setArg<cl::Memory>(0, *autoclearBuffers[base]);
        clearThreeBuffersKernel.setArg<cl_int>(1, autoclearBufferSizes[base]);
        clearThreeBuffersKernel.setArg<cl::Memory>(2, *autoclearBuffers[base+1]);
        clearThreeBuffersKernel.setArg<cl_int>(3, autoclearBufferSizes[base+1]);
        clearThreeBuffersKernel.setArg<cl::Memory>(4, *autoclearBuffers[base+2]);
        clearThreeBuffersKernel.setArg<cl_int>(5, autoclearBufferSizes[base+2]);
        executeKernel(clearThreeBuffersKernel, max(max(autoclearBufferSizes[base], autoclearBufferSizes[base+1]), autoclearBufferSizes[base+2]), 128);
    }
    else if (total-base == 2) {
        clearTwoBuffersKernel.setArg<cl::Memory>(0, *autoclearBuffers[base]);
        clearTwoBuffersKernel.setArg<cl_int>(1, autoclearBufferSizes[base]);
        clearTwoBuffersKernel.setArg<cl::Memory>(2, *autoclearBuffers[base+1]);
        clearTwoBuffersKernel.setArg<cl_int>(3, autoclearBufferSizes[base+1]);
        executeKernel(clearTwoBuffersKernel, max(autoclearBufferSizes[base], autoclearBufferSizes[base+1]), 128);
    }
    else if (total-base == 1) {
        clearBuffer(*autoclearBuffers[base], autoclearBufferSizes[base]);
    }
}

void OpenCLContext::reduceForces() {
    if (supports64BitGlobalAtomics)
        executeKernel(reduceForcesKernel, paddedNumAtoms, 128);
    else
        reduceBuffer(*forceBuffers, numForceBuffers);
}

void OpenCLContext::reduceBuffer(OpenCLArray& array, int numBuffers) {
    int bufferSize = array.getSize()/numBuffers;
    reduceReal4Kernel.setArg<cl::Buffer>(0, array.getDeviceBuffer());
    reduceReal4Kernel.setArg<cl_int>(1, bufferSize);
    reduceReal4Kernel.setArg<cl_int>(2, numBuffers);
    executeKernel(reduceReal4Kernel, bufferSize, 128);
}

void OpenCLContext::tagAtomsInMolecule(int atom, int molecule, vector<int>& atomMolecule, vector<vector<int> >& atomBonds) {
    // Recursively tag atoms as belonging to a particular molecule.

    atomMolecule[atom] = molecule;
    for (int i = 0; i < (int) atomBonds[atom].size(); i++)
        if (atomMolecule[atomBonds[atom][i]] == -1)
            tagAtomsInMolecule(atomBonds[atom][i], molecule, atomMolecule, atomBonds);
}

/**
 * This class ensures that atom reordering doesn't break virtual sites.
 */
class OpenCLContext::VirtualSiteInfo : public OpenCLForceInfo {
public:
    VirtualSiteInfo(const System& system) : OpenCLForceInfo(0) {
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


void OpenCLContext::findMoleculeGroups() {
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

        // Now tag atoms by which molecule they belong to.

        vector<int> atomMolecule(numAtoms, -1);
        int numMolecules = 0;
        for (int i = 0; i < numAtoms; i++)
            if (atomMolecule[i] == -1)
                tagAtomsInMolecule(i, numMolecules++, atomMolecule, atomBonds);
        vector<vector<int> > atomIndices(numMolecules);
        for (int i = 0; i < numAtoms; i++)
            atomIndices[atomMolecule[i]].push_back(i);

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

void OpenCLContext::invalidateMolecules() {
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
    
    vector<mm_int4> newCellOffsets(numAtoms);
    if (useDoublePrecision) {
        vector<mm_double4> oldPosq(paddedNumAtoms);
        vector<mm_double4> newPosq(paddedNumAtoms, mm_double4(0,0,0,0));
        vector<mm_double4> oldVelm(paddedNumAtoms);
        vector<mm_double4> newVelm(paddedNumAtoms, mm_double4(0,0,0,0));
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
        vector<mm_float4> oldPosq(paddedNumAtoms);
        vector<mm_float4> newPosq(paddedNumAtoms, mm_float4(0,0,0,0));
        vector<mm_float4> oldPosqCorrection(paddedNumAtoms);
        vector<mm_float4> newPosqCorrection(paddedNumAtoms, mm_float4(0,0,0,0));
        vector<mm_double4> oldVelm(paddedNumAtoms);
        vector<mm_double4> newVelm(paddedNumAtoms, mm_double4(0,0,0,0));
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
        vector<mm_float4> oldPosq(paddedNumAtoms);
        vector<mm_float4> newPosq(paddedNumAtoms, mm_float4(0,0,0,0));
        vector<mm_float4> oldVelm(paddedNumAtoms);
        vector<mm_float4> newVelm(paddedNumAtoms, mm_float4(0,0,0,0));
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

void OpenCLContext::reorderAtoms() {
    atomsWereReordered = false;
    if (numAtoms == 0 || nonbonded == NULL || !nonbonded->getUseCutoff() || stepsSinceReorder < 100) {
        stepsSinceReorder++;
        return;
    }
    atomsWereReordered = true;
    stepsSinceReorder = 0;
    if (useDoublePrecision)
        reorderAtomsImpl<cl_double, mm_double4, cl_double, mm_double4>();
    else if (useMixedPrecision)
        reorderAtomsImpl<cl_float, mm_float4, cl_double, mm_double4>();
    else
        reorderAtomsImpl<cl_float, mm_float4, cl_float, mm_float4>();
    nonbonded->updateNeighborListSize();
}

template <class Real, class Real4, class Mixed, class Mixed4>
void OpenCLContext::reorderAtomsImpl() {

    // Find the range of positions and the number of bins along each axis.

    vector<Real4> oldPosq(paddedNumAtoms);
    vector<Real4> oldPosqCorrection(paddedNumAtoms);
    vector<Mixed4> oldVelm(paddedNumAtoms);
    posq->download(oldPosq);
    velm->download(oldVelm);
    if (useMixedPrecision)
        posqCorrection->download(oldPosqCorrection);
    Real minx = oldPosq[0].x, maxx = oldPosq[0].x;
    Real miny = oldPosq[0].y, maxy = oldPosq[0].y;
    Real minz = oldPosq[0].z, maxz = oldPosq[0].z;
    if (nonbonded->getUsePeriodic()) {
        minx = miny = minz = 0.0;
        maxx = periodicBoxSizeDouble.x;
        maxy = periodicBoxSizeDouble.y;
        maxz = periodicBoxSizeDouble.z;
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
    vector<Real4> newPosq(paddedNumAtoms, Real4(0,0,0,0));
    vector<Real4> newPosqCorrection(paddedNumAtoms, Real4(0,0,0,0));
    vector<Mixed4> newVelm(paddedNumAtoms, Mixed4(0,0,0,0));
    vector<mm_int4> newCellOffsets(numAtoms);
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
        }
        if (nonbonded->getUsePeriodic()) {
            // Move each molecule position into the same box.

            for (int i = 0; i < numMolecules; i++) {
                int xcell = (int) floor(molPos[i].x*invPeriodicBoxSizeDouble.x);
                int ycell = (int) floor(molPos[i].y*invPeriodicBoxSizeDouble.y);
                int zcell = (int) floor(molPos[i].z*invPeriodicBoxSizeDouble.z);
                Real dx = xcell*periodicBoxSizeDouble.x;
                Real dy = ycell*periodicBoxSizeDouble.y;
                Real dz = zcell*periodicBoxSizeDouble.z;
                if (dx != 0.0f || dy != 0.0f || dz != 0.0f) {
                    molPos[i].x -= dx;
                    molPos[i].y -= dy;
                    molPos[i].z -= dz;
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
            binWidth = (Real) (0.2*nonbonded->getCutoffDistance());
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

void OpenCLContext::addReorderListener(ReorderListener* listener) {
    reorderListeners.push_back(listener);
}

void OpenCLContext::addPreComputation(ForcePreComputation* computation) {
    preComputations.push_back(computation);
}

void OpenCLContext::addPostComputation(ForcePostComputation* computation) {
    postComputations.push_back(computation);
}

struct OpenCLContext::WorkThread::ThreadData {
    ThreadData(std::queue<OpenCLContext::WorkTask*>& tasks, bool& waiting,  bool& finished,
            pthread_mutex_t& queueLock, pthread_cond_t& waitForTaskCondition, pthread_cond_t& queueEmptyCondition) :
        tasks(tasks), waiting(waiting), finished(finished), queueLock(queueLock),
        waitForTaskCondition(waitForTaskCondition), queueEmptyCondition(queueEmptyCondition) {
    }
    std::queue<OpenCLContext::WorkTask*>& tasks;
    bool& waiting;
    bool& finished;
    pthread_mutex_t& queueLock;
    pthread_cond_t& waitForTaskCondition;
    pthread_cond_t& queueEmptyCondition;
};

static void* threadBody(void* args) {
    OpenCLContext::WorkThread::ThreadData& data = *reinterpret_cast<OpenCLContext::WorkThread::ThreadData*>(args);
    while (!data.finished || data.tasks.size() > 0) {
        pthread_mutex_lock(&data.queueLock);
        while (data.tasks.empty() && !data.finished) {
            data.waiting = true;
            pthread_cond_signal(&data.queueEmptyCondition);
            pthread_cond_wait(&data.waitForTaskCondition, &data.queueLock);
        }
        OpenCLContext::WorkTask* task = NULL;
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

OpenCLContext::WorkThread::WorkThread() : waiting(true), finished(false) {
    pthread_mutex_init(&queueLock, NULL);
    pthread_cond_init(&waitForTaskCondition, NULL);
    pthread_cond_init(&queueEmptyCondition, NULL);
    ThreadData* data = new ThreadData(tasks, waiting, finished, queueLock, waitForTaskCondition, queueEmptyCondition);
    pthread_create(&thread, NULL, threadBody, data);
}

OpenCLContext::WorkThread::~WorkThread() {
    pthread_mutex_lock(&queueLock);
    finished = true;
    pthread_cond_broadcast(&waitForTaskCondition);
    pthread_mutex_unlock(&queueLock);
    pthread_join(thread, NULL);
    pthread_mutex_destroy(&queueLock);
    pthread_cond_destroy(&waitForTaskCondition);
    pthread_cond_destroy(&queueEmptyCondition);
}

void OpenCLContext::WorkThread::addTask(OpenCLContext::WorkTask* task) {
    pthread_mutex_lock(&queueLock);
    tasks.push(task);
    waiting = false;
    pthread_cond_signal(&waitForTaskCondition);
    pthread_mutex_unlock(&queueLock);
}

bool OpenCLContext::WorkThread::isWaiting() {
    return waiting;
}

bool OpenCLContext::WorkThread::isFinished() {
    return finished;
}

void OpenCLContext::WorkThread::flush() {
    pthread_mutex_lock(&queueLock);
    while (!waiting)
       pthread_cond_wait(&queueEmptyCondition, &queueLock);
    pthread_mutex_unlock(&queueLock);
}
