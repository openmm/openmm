/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2024 Stanford University and the Authors.      *
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
#include "OpenCLEvent.h"
#include "OpenCLForceInfo.h"
#include "OpenCLIntegrationUtilities.h"
#include "OpenCLKernelSources.h"
#include "OpenCLNonbondedUtilities.h"
#include "OpenCLProgram.h"
#include "openmm/common/ComputeArray.h"
#include "openmm/MonteCarloFlexibleBarostat.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/VirtualSite.h"
#include "openmm/internal/ContextImpl.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <typeinfo>

using namespace OpenMM;
using namespace std;

// Uncomment the following line to enable profiling of all kernel launches.  The results are written
// to stdout in the JSON format used by https://ui.perfetto.dev.
//#define ENABLE_PROFILING

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

static bool isSupported(cl::Platform platform) {
    string vendor = platform.getInfo<CL_PLATFORM_VENDOR>();
    return (vendor.find("NVIDIA") == 0 ||
            vendor.find("Advanced Micro Devices") == 0 ||
            vendor.find("Apple") == 0 ||
            vendor.find("Intel") == 0);
}

OpenCLContext::OpenCLContext(const System& system, int platformIndex, int deviceIndex, const string& precision, OpenCLPlatform::PlatformData& platformData, OpenCLContext* originalContext) :
        ComputeContext(system), platformData(platformData), numForceBuffers(0), hasAssignedPosqCharges(false), profileStartTime(0),
        integration(NULL), expression(NULL), bonded(NULL), nonbonded(NULL), pinnedBuffer(NULL) {
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
    try {
        contextIndex = platformData.contexts.size();
        std::vector<cl::Platform> platforms;
        cl::Platform::get(&platforms);
        if (platformIndex < -1 || platformIndex >= (int) platforms.size())
            throw OpenMMException("Illegal value for OpenCLPlatformIndex: "+intToString(platformIndex));
        if (platforms.size() > 1 && platformIndex == -1 && deviceIndex != -1)
            throw OpenMMException("Specified DeviceIndex but not OpenCLPlatformIndex.  When multiple platforms are available, a platform index is needed to specify a device.");
        const int minThreadBlockSize = 32;

        int bestSpeed = -1;
        int bestDevice = -1;
        int bestPlatform = -1;
        bool bestSupported = false;
        for (int j = 0; j < platforms.size(); j++) {
            // If they supplied a valid platformIndex, we only look through that platform
            if (j != platformIndex && platformIndex != -1)
                continue;

            // Always prefer a supported platform over an unsupported one.
            bool supported = isSupported(platforms[j]);
            if (!supported && bestSupported)
                continue;
            string platformVendor = platforms[j].getInfo<CL_PLATFORM_VENDOR>();
            vector<cl::Device> devices;
            try {
                platforms[j].getDevices(CL_DEVICE_TYPE_ALL, &devices);
            }
            catch (...) {
                // There are no devices available for this platform.
                continue;
            }
            if (deviceIndex < -1 || deviceIndex >= (int) devices.size())
                throw OpenMMException("Illegal value for DeviceIndex: "+intToString(deviceIndex));

            for (int i = 0; i < (int) devices.size(); i++) {
                // If they supplied a valid deviceIndex, we only look through that one
                if (i != deviceIndex && deviceIndex != -1)
                    continue;
                if (platformVendor == "Apple" && (devices[i].getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_CPU))
                    continue; // The CPU device on OS X won't work correctly.
                if (useMixedPrecision || useDoublePrecision) {
                    bool supportsDouble = (devices[i].getInfo<CL_DEVICE_EXTENSIONS>().find("cl_khr_fp64") != string::npos);
                    if (!supportsDouble)
                        continue; // This device does not support double precision.
                }
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
                    }
                    catch (cl::Error err) {
                        // Runtime does not support the queries so use default.
                    }
                }
                int speed = devices[i].getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>()*processingElementsPerComputeUnit*devices[i].getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>();
                if (maxSize >= minThreadBlockSize && (speed > bestSpeed || (supported && !bestSupported))) {
                    bestDevice = i;
                    bestSpeed = speed;
                    bestPlatform = j;
                    bestSupported = supported;
                }
            }
        }

        if (bestPlatform == -1)
            throw OpenMMException("No compatible OpenCL platform is available");

        if (bestDevice == -1)
            throw OpenMMException("No compatible OpenCL device is available");
        
        if (!bestSupported)
            cout << "WARNING: Using an unsupported OpenCL implementation.  Results may be incorrect." << endl;

        vector<cl::Device> devices;
        platforms[bestPlatform].getDevices(CL_DEVICE_TYPE_ALL, &devices);
        string platformVendor = platforms[bestPlatform].getInfo<CL_PLATFORM_VENDOR>();
        device = devices[bestDevice];

        this->deviceIndex = bestDevice;
        this->platformIndex = bestPlatform;
        if (device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>() < minThreadBlockSize)
            throw OpenMMException("The specified OpenCL device is not compatible with OpenMM");
        compilationDefines["WORK_GROUP_SIZE"] = intToString(ThreadBlockSize);
        if (platformVendor.size() >= 5 && platformVendor.substr(0, 5) == "Intel")
            defaultOptimizationOptions = "";
        else
            defaultOptimizationOptions = "-cl-mad-enable -cl-no-signed-zeros";
        supports64BitGlobalAtomics = (device.getInfo<CL_DEVICE_EXTENSIONS>().find("cl_khr_int64_base_atomics") != string::npos);
        supportsDoublePrecision = (device.getInfo<CL_DEVICE_EXTENSIONS>().find("cl_khr_fp64") != string::npos);
        if ((useDoublePrecision || useMixedPrecision) && !supportsDoublePrecision)
            throw OpenMMException("This device does not support double precision");
        string vendor = device.getInfo<CL_DEVICE_VENDOR>();
        int numThreadBlocksPerComputeUnit = 6;
      
        if (vendor.size() >= 5 && vendor.substr(0, 5) == "Apple") {
            simdWidth = 32;

            // 768 threads per GPU core.
            numThreadBlocksPerComputeUnit = 12;
        }
        else if (vendor.size() >= 6 && vendor.substr(0, 6) == "NVIDIA") {
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
                if (computeCapabilityMajor == 5) {
                    // Workaround for a bug in Maxwell on CUDA 6.x.

                    string platformVersion = platforms[bestPlatform].getInfo<CL_PLATFORM_VERSION>();
                    if (platformVersion.find("CUDA 6") != string::npos)
                        supports64BitGlobalAtomics = false;
                }
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
                        // Must catch cl:Error as will fail if runtime does not support queries.

                        cl_uint simdPerComputeUnit = device.getInfo<CL_DEVICE_SIMD_PER_COMPUTE_UNIT_AMD>();
                        simdWidth = device.getInfo<CL_DEVICE_WAVEFRONT_WIDTH_AMD>();

                        // If the GPU has multiple SIMDs per compute unit then it is uses the scalar instruction
                        // set instead of the VLIW instruction set. It therefore needs more thread blocks per
                        // compute unit to hide memory latency.
                        if (simdPerComputeUnit > 1) {
                            if (simdWidth == 32)
                                numThreadBlocksPerComputeUnit = 6*simdPerComputeUnit; // Navi seems to like more thread blocks than older GPUs
                            else
                                numThreadBlocksPerComputeUnit = 4*simdPerComputeUnit;
                        }

                        // If the queries are supported then must be newer than SDK 2.4.
                        amdPostSdk2_4 = true;
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
        if (supports64BitGlobalAtomics)
            compilationDefines["SUPPORTS_64_BIT_ATOMICS"] = "";
        if (supportsDoublePrecision)
            compilationDefines["SUPPORTS_DOUBLE_PRECISION"] = "";
        if (simdWidth >= 32)
            compilationDefines["SYNC_WARPS"] = "mem_fence(CLK_LOCAL_MEM_FENCE)";
        else
            compilationDefines["SYNC_WARPS"] = "barrier(CLK_LOCAL_MEM_FENCE)";
        vector<cl::Device> contextDevices;
        contextDevices.push_back(device);
        cl_context_properties cprops[] = {CL_CONTEXT_PLATFORM, (cl_context_properties) platforms[bestPlatform](), 0};
        if (originalContext == NULL) {
            context = cl::Context(contextDevices, cprops, errorCallback);
#ifdef ENABLE_PROFILING
            defaultQueue = cl::CommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE);
            printf("[ ");
#else
            defaultQueue = cl::CommandQueue(context, device);
#endif
        }
        else {
            context = originalContext->context;
            defaultQueue = originalContext->defaultQueue;
        }
        currentQueue = defaultQueue;
        numAtoms = system.getNumParticles();
        paddedNumAtoms = TileSize*((numAtoms+TileSize-1)/TileSize);
        numAtomBlocks = (paddedNumAtoms+(TileSize-1))/TileSize;
        numThreadBlocks = numThreadBlocksPerComputeUnit*device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
        if (useDoublePrecision) {
            posq.initialize<mm_double4>(*this, paddedNumAtoms, "posq");
            velm.initialize<mm_double4>(*this, paddedNumAtoms, "velm");
            compilationDefines["USE_DOUBLE_PRECISION"] = "1";
            compilationDefines["convert_real4"] = "convert_double4";
            compilationDefines["make_real2"] = "make_double2";
            compilationDefines["make_real3"] = "make_double3";
            compilationDefines["make_real4"] = "make_double4";
            compilationDefines["convert_mixed4"] = "convert_double4";
            compilationDefines["make_mixed2"] = "make_double2";
            compilationDefines["make_mixed3"] = "make_double3";
            compilationDefines["make_mixed4"] = "make_double4";
        }
        else if (useMixedPrecision) {
            posq.initialize<mm_float4>(*this, paddedNumAtoms, "posq");
            posqCorrection.initialize<mm_float4>(*this, paddedNumAtoms, "posq");
            velm.initialize<mm_double4>(*this, paddedNumAtoms, "velm");
            compilationDefines["USE_MIXED_PRECISION"] = "1";
            compilationDefines["convert_real4"] = "convert_float4";
            compilationDefines["make_real2"] = "make_float2";
            compilationDefines["make_real3"] = "make_float3";
            compilationDefines["make_real4"] = "make_float4";
            compilationDefines["convert_mixed4"] = "convert_double4";
            compilationDefines["make_mixed2"] = "make_double2";
            compilationDefines["make_mixed3"] = "make_double3";
            compilationDefines["make_mixed4"] = "make_double4";
        }
        else {
            posq.initialize<mm_float4>(*this, paddedNumAtoms, "posq");
            velm.initialize<mm_float4>(*this, paddedNumAtoms, "velm");
            compilationDefines["convert_real4"] = "convert_float4";
            compilationDefines["make_real2"] = "make_float2";
            compilationDefines["make_real3"] = "make_float3";
            compilationDefines["make_real4"] = "make_float4";
            compilationDefines["convert_mixed4"] = "convert_float4";
            compilationDefines["make_mixed2"] = "make_float2";
            compilationDefines["make_mixed3"] = "make_float3";
            compilationDefines["make_mixed4"] = "make_float4";
        }
        longForceBuffer.initialize<cl_long>(*this, 3*paddedNumAtoms, "longForceBuffer");
        posCellOffsets.resize(paddedNumAtoms, mm_int4(0, 0, 0, 0));
        atomIndexDevice.initialize<cl_int>(*this, paddedNumAtoms, "atomIndexDevice");
        atomIndex.resize(paddedNumAtoms);
        for (int i = 0; i < paddedNumAtoms; ++i)
            atomIndex[i] = i;
        atomIndexDevice.upload(atomIndex);
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
    reduceForcesKernel = cl::Kernel(utilities, "reduceForces");
    reduceEnergyKernel = cl::Kernel(utilities, "reduceEnergy");
    setChargesKernel = cl::Kernel(utilities, "setCharges");

    // Decide whether native_sqrt(), native_rsqrt(), and native_recip() are sufficiently accurate to use.

    if (!useDoublePrecision) {
        cl::Kernel accuracyKernel(utilities, "determineNativeAccuracy");
        OpenCLArray valuesArray(*this, 20, sizeof(mm_float8), "values");
        vector<mm_float8> values(valuesArray.getSize());
        float nextValue = 1e-4f;
        for (auto& val : values) {
            val.s0 = nextValue;
            nextValue *= (float) M_PI;
        }
        valuesArray.upload(values);
        accuracyKernel.setArg<cl::Buffer>(0, valuesArray.getDeviceBuffer());
        accuracyKernel.setArg<cl_int>(1, values.size());
        executeKernel(accuracyKernel, values.size());
        valuesArray.download(values);
        double maxSqrtError = 0.0, maxRsqrtError = 0.0, maxRecipError = 0.0, maxExpError = 0.0, maxLogError = 0.0;
        for (auto& val : values) {
            double v = val.s0;
            double correctSqrt = sqrt(v);
            maxSqrtError = max(maxSqrtError, fabs(correctSqrt-val.s1)/correctSqrt);
            maxRsqrtError = max(maxRsqrtError, fabs(1.0/correctSqrt-val.s2)*correctSqrt);
            maxRecipError = max(maxRecipError, fabs(1.0/v-val.s3)/val.s3);
            maxExpError = max(maxExpError, fabs(exp(v)-val.s4)/val.s4);
            maxLogError = max(maxLogError, fabs(log(v)-val.s5)/val.s5);
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
    compilationDefines["POW"] = "pow";
    compilationDefines["COS"] = "cos";
    compilationDefines["SIN"] = "sin";
    compilationDefines["TAN"] = "tan";
    compilationDefines["ACOS"] = "acos";
    compilationDefines["ASIN"] = "asin";
    compilationDefines["ATAN"] = "atan";
    compilationDefines["ERF"] = "erf";
    compilationDefines["ERFC"] = "erfc";

    // Set defines for applying periodic boundary conditions.

    Vec3 boxVectors[3];
    system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    boxIsTriclinic = (boxVectors[0][1] != 0.0 || boxVectors[0][2] != 0.0 ||
                      boxVectors[1][0] != 0.0 || boxVectors[1][2] != 0.0 ||
                      boxVectors[2][0] != 0.0 || boxVectors[2][1] != 0.0);
    for (int i = 0; i < system.getNumForces(); i++)
        if (dynamic_cast<const MonteCarloFlexibleBarostat*>(&system.getForce(i)) != NULL)
            boxIsTriclinic = true;
    if (boxIsTriclinic) {
        compilationDefines["APPLY_PERIODIC_TO_DELTA(delta)"] =
            "{"
            "real scale3 = floor(delta.z*invPeriodicBoxSize.z+0.5f); \\\n"
            "delta.xyz -= scale3*periodicBoxVecZ.xyz; \\\n"
            "real scale2 = floor(delta.y*invPeriodicBoxSize.y+0.5f); \\\n"
            "delta.xy -= scale2*periodicBoxVecY.xy; \\\n"
            "real scale1 = floor(delta.x*invPeriodicBoxSize.x+0.5f); \\\n"
            "delta.x -= scale1*periodicBoxVecX.x;}";
        compilationDefines["APPLY_PERIODIC_TO_POS(pos)"] =
            "{"
            "real scale3 = floor(pos.z*invPeriodicBoxSize.z); \\\n"
            "pos.xyz -= scale3*periodicBoxVecZ.xyz; \\\n"
            "real scale2 = floor(pos.y*invPeriodicBoxSize.y); \\\n"
            "pos.xy -= scale2*periodicBoxVecY.xy; \\\n"
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
            "delta.xyz -= floor(delta.xyz*invPeriodicBoxSize.xyz+0.5f)*periodicBoxSize.xyz;";
        compilationDefines["APPLY_PERIODIC_TO_POS(pos)"] =
            "pos.xyz -= floor(pos.xyz*invPeriodicBoxSize.xyz)*periodicBoxSize.xyz;";
        compilationDefines["APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, center)"] =
            "{"
            "pos.x -= floor((pos.x-center.x)*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x; \\\n"
            "pos.y -= floor((pos.y-center.y)*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y; \\\n"
            "pos.z -= floor((pos.z-center.z)*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;}";
    }

    // Create utilities objects.

    bonded = new OpenCLBondedUtilities(*this);
    nonbonded = new OpenCLNonbondedUtilities(*this);
    integration = new OpenCLIntegrationUtilities(*this, system);
    expression = new OpenCLExpressionUtilities(*this);
}

OpenCLContext::~OpenCLContext() {
    for (auto force : forces)
        delete force;
    for (auto listener : reorderListeners)
        delete listener;
    for (auto computation : preComputations)
        delete computation;
    for (auto computation : postComputations)
        delete computation;
    if (pinnedBuffer != NULL)
        delete pinnedBuffer;
    if (integration != NULL)
        delete integration;
    if (expression != NULL)
        delete expression;
    if (bonded != NULL)
        delete bonded;
    if (nonbonded != NULL)
        delete nonbonded;
#ifdef ENABLE_PROFILING
    printProfilingEvents();
    printf(" ]\n");
#endif
}

void OpenCLContext::initialize() {
    bonded->initialize(system);
    numForceBuffers = std::max(numForceBuffers, (int) platformData.contexts.size());
    int energyBufferSize = max(numThreadBlocks*ThreadBlockSize, nonbonded->getNumEnergyBuffers());
    int numComputeUnits = device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
    if (useDoublePrecision) {
        forceBuffers.initialize<mm_double4>(*this, paddedNumAtoms*numForceBuffers, "forceBuffers");
        force.initialize<mm_double4>(*this, &forceBuffers.getDeviceBuffer(), paddedNumAtoms, "force");
        energyBuffer.initialize<cl_double>(*this, energyBufferSize, "energyBuffer");
        energySum.initialize<cl_double>(*this, numComputeUnits, "energySum");
    }
    else if (useMixedPrecision) {
        forceBuffers.initialize<mm_float4>(*this, paddedNumAtoms*numForceBuffers, "forceBuffers");
        force.initialize<mm_float4>(*this, &forceBuffers.getDeviceBuffer(), paddedNumAtoms, "force");
        energyBuffer.initialize<cl_double>(*this, energyBufferSize, "energyBuffer");
        energySum.initialize<cl_double>(*this, numComputeUnits, "energySum");
    }
    else {
        forceBuffers.initialize<mm_float4>(*this, paddedNumAtoms*numForceBuffers, "forceBuffers");
        force.initialize<mm_float4>(*this, &forceBuffers.getDeviceBuffer(), paddedNumAtoms, "force");
        energyBuffer.initialize<cl_float>(*this, energyBufferSize, "energyBuffer");
        energySum.initialize<cl_float>(*this, numComputeUnits, "energySum");
    }
    reduceForcesKernel.setArg<cl::Buffer>(0, longForceBuffer.getDeviceBuffer());
    reduceForcesKernel.setArg<cl::Buffer>(1, forceBuffers.getDeviceBuffer());
    reduceForcesKernel.setArg<cl_int>(2, paddedNumAtoms);
    reduceForcesKernel.setArg<cl_int>(3, numForceBuffers);
    addAutoclearBuffer(longForceBuffer);
    addAutoclearBuffer(forceBuffers);
    addAutoclearBuffer(energyBuffer);
    int numEnergyParamDerivs = energyParamDerivNames.size();
    if (numEnergyParamDerivs > 0) {
        if (useDoublePrecision || useMixedPrecision)
            energyParamDerivBuffer.initialize<cl_double>(*this, numEnergyParamDerivs*energyBufferSize, "energyParamDerivBuffer");
        else
            energyParamDerivBuffer.initialize<cl_float>(*this, numEnergyParamDerivs*energyBufferSize, "energyParamDerivBuffer");
        addAutoclearBuffer(energyParamDerivBuffer);
    }
    int bufferBytes = max(max((int) velm.getSize()*velm.getElementSize(),
            energyBufferSize*energyBuffer.getElementSize()),
            (int) longForceBuffer.getSize()*longForceBuffer.getElementSize());
    pinnedBuffer = new cl::Buffer(context, CL_MEM_ALLOC_HOST_PTR, bufferBytes);
    pinnedMemory = currentQueue.enqueueMapBuffer(*pinnedBuffer, CL_TRUE, CL_MAP_READ | CL_MAP_WRITE, 0, bufferBytes);
    for (int i = 0; i < numAtoms; i++) {
        double mass = system.getParticleMass(i);
        if (useDoublePrecision || useMixedPrecision)
            ((mm_double4*) pinnedMemory)[i] = mm_double4(0.0, 0.0, 0.0, mass == 0.0 ? 0.0 : 1.0/mass);
        else
            ((mm_float4*) pinnedMemory)[i] = mm_float4(0.0f, 0.0f, 0.0f, mass == 0.0 ? 0.0f : (cl_float) (1.0/mass));
    }
    velm.upload(pinnedMemory);
    findMoleculeGroups();
    nonbonded->initialize(system);
}

void OpenCLContext::initializeContexts() {
    getPlatformData().initializeContexts(system);
}

void OpenCLContext::addForce(ComputeForceInfo* force) {
    ComputeContext::addForce(force);
    OpenCLForceInfo* clinfo = dynamic_cast<OpenCLForceInfo*>(force);
    if (clinfo != NULL)
        requestForceBuffers(clinfo->getRequiredForceBuffers());
}

void OpenCLContext::requestForceBuffers(int minBuffers) {
    numForceBuffers = std::max(numForceBuffers, minBuffers);
}

cl::Program OpenCLContext::createProgram(const string source, const char* optimizationFlags) {
    return createProgram(source, map<string, string>(), optimizationFlags);
}

cl::Program OpenCLContext::createProgram(const string source, const map<string, string>& defines, const char* optimizationFlags) {
    string options = (optimizationFlags == NULL ? defaultOptimizationOptions : string(optimizationFlags));
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
    src << OpenCLKernelSources::common << endl;
    for (auto& pair : defines) {
        src << "#define " << pair.first;
        if (!pair.second.empty())
            src << " " << pair.second;
        src << endl;
    }
    if (!defines.empty())
        src << endl;
    src << source << endl;
    cl::Program::Sources sources({src.str()});
    cl::Program program(context, sources);
    try {
        program.build(vector<cl::Device>(1, device), options.c_str());
    } catch (cl::Error err) {
        throw OpenMMException("Error compiling kernel: "+program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device));
    }
    return program;
}

vector<ComputeContext*> OpenCLContext::getAllContexts() {
    vector<ComputeContext*> result;
    for (OpenCLContext* c : platformData.contexts)
        result.push_back(c);
    return result;
}

double& OpenCLContext::getEnergyWorkspace() {
    return platformData.contextEnergy[contextIndex];
}

cl::CommandQueue& OpenCLContext::getQueue() {
    return currentQueue;
}

void OpenCLContext::setQueue(cl::CommandQueue& queue) {
    currentQueue = queue;
}

void OpenCLContext::restoreDefaultQueue() {
    currentQueue = defaultQueue;
}

OpenCLArray* OpenCLContext::createArray() {
    return new OpenCLArray();
}

ComputeEvent OpenCLContext::createEvent() {
    return shared_ptr<ComputeEventImpl>(new OpenCLEvent(*this));
}

ComputeProgram OpenCLContext::compileProgram(const std::string source, const std::map<std::string, std::string>& defines) {
    cl::Program program = createProgram(source, defines);
    return shared_ptr<ComputeProgramImpl>(new OpenCLProgram(*this, program));
}

OpenCLArray& OpenCLContext::unwrap(ArrayInterface& array) const {
    OpenCLArray* clarray;
    ComputeArray* wrapper = dynamic_cast<ComputeArray*>(&array);
    if (wrapper != NULL)
        clarray = dynamic_cast<OpenCLArray*>(&wrapper->getArray());
    else
        clarray = dynamic_cast<OpenCLArray*>(&array);
    if (clarray == NULL)
        throw OpenMMException("Array argument is not an OpenCLArray");
    return *clarray;
}

void OpenCLContext::executeKernel(cl::Kernel& kernel, int workUnits, int blockSize) {
    if (blockSize == -1)
        blockSize = ThreadBlockSize;
    int size = std::min((workUnits+blockSize-1)/blockSize, numThreadBlocks)*blockSize;
    try {
#ifdef ENABLE_PROFILING
    cl::Event event;
    currentQueue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(size), cl::NDRange(blockSize), NULL, &event);
    profilingEvents.push_back(event);
    profilingKernelNames.push_back(kernel.getInfo<CL_KERNEL_FUNCTION_NAME>());
    if (profilingEvents.size() >= 500)
        printProfilingEvents();
#else
        currentQueue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(size), cl::NDRange(blockSize));
#endif
    }
    catch (cl::Error err) {
        stringstream str;
        str<<"Error invoking kernel "<<kernel.getInfo<CL_KERNEL_FUNCTION_NAME>()<<": "<<err.what()<<" ("<<err.err()<<")";
        throw OpenMMException(str.str());
    }
}

void OpenCLContext::printProfilingEvents() {
    for (int i = 0; i < profilingEvents.size(); i++) {
        cl::Event event = profilingEvents[i];
        event.wait();
        cl_ulong start, end;
        event.getProfilingInfo(CL_PROFILING_COMMAND_START, &start);
        event.getProfilingInfo(CL_PROFILING_COMMAND_END, &end);
        if (profileStartTime == 0)
            profileStartTime = start;
        else
            printf(",\n");
        printf("{ \"pid\":1, \"tid\":1, \"ts\":%.6g, \"dur\":%g, \"ph\":\"X\", \"name\":\"%s\" }",
                0.001*(start-profileStartTime), 0.001*(end-start), profilingKernelNames[i].c_str());
    }
    profilingEvents.clear();
    profilingKernelNames.clear();
}

int OpenCLContext::computeThreadBlockSize(double memory) const {
    int maxShared = device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>();
    // On some implementations, more local memory gets used than we calculate by
    // adding up the sizes of the fields.  To be safe, include a factor of 0.5.
    int max = (int) (0.5*maxShared/memory);
    if (max < 64)
        return 32;
    int threads = 64;
    while (threads+64 < max)
        threads += 64;
    return threads;
}

void OpenCLContext::clearBuffer(ArrayInterface& array) {
    clearBuffer(unwrap(array).getDeviceBuffer(), array.getSize()*array.getElementSize());
}

void OpenCLContext::clearBuffer(cl::Memory& memory, int size) {
    int words = size/4;
    clearBufferKernel.setArg<cl::Memory>(0, memory);
    clearBufferKernel.setArg<cl_int>(1, words);
    executeKernel(clearBufferKernel, words, 128);
}

void OpenCLContext::addAutoclearBuffer(ArrayInterface& array) {
    addAutoclearBuffer(unwrap(array).getDeviceBuffer(), array.getSize()*array.getElementSize());
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
        clearBuffer(*autoclearBuffers[base], autoclearBufferSizes[base]*4);
    }
}

void OpenCLContext::reduceForces() {
    executeKernel(reduceForcesKernel, paddedNumAtoms, 128);
}

void OpenCLContext::reduceBuffer(OpenCLArray& array, OpenCLArray& longBuffer, int numBuffers) {
    int bufferSize = array.getSize()/numBuffers;
    reduceReal4Kernel.setArg<cl::Buffer>(0, array.getDeviceBuffer());
    reduceReal4Kernel.setArg<cl::Buffer>(1, longBuffer.getDeviceBuffer());
    reduceReal4Kernel.setArg<cl_int>(2, bufferSize);
    reduceReal4Kernel.setArg<cl_int>(3, numBuffers);
    executeKernel(reduceReal4Kernel, bufferSize, 128);
}

double OpenCLContext::reduceEnergy() {
    int workGroupSize  = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
    if (workGroupSize > 512)
        workGroupSize = 512;
    reduceEnergyKernel.setArg<cl::Buffer>(0, energyBuffer.getDeviceBuffer());
    reduceEnergyKernel.setArg<cl::Buffer>(1, energySum.getDeviceBuffer());
    reduceEnergyKernel.setArg<cl_int>(2, energyBuffer.getSize());
    reduceEnergyKernel.setArg<cl_int>(3, workGroupSize);
    reduceEnergyKernel.setArg(4, workGroupSize*energyBuffer.getElementSize(), NULL);
    executeKernel(reduceEnergyKernel, workGroupSize*energySum.getSize(), workGroupSize);
    energySum.download(pinnedMemory);
    double result = 0;
    if (getUseDoublePrecision() || getUseMixedPrecision()) {
        for (int i = 0; i < energySum.getSize(); i++)
            result += ((double*) pinnedMemory)[i];
    }
    else {
        for (int i = 0; i < energySum.getSize(); i++)
            result += ((float*) pinnedMemory)[i];
    }
    return result;
}

void OpenCLContext::setCharges(const vector<double>& charges) {
    if (!chargeBuffer.isInitialized())
        chargeBuffer.initialize(*this, numAtoms, useDoublePrecision ? sizeof(double) : sizeof(float), "chargeBuffer");
    vector<double> c(numAtoms);
    for (int i = 0; i < numAtoms; i++)
        c[i] = charges[i];
    chargeBuffer.upload(c, true);
    setChargesKernel.setArg<cl::Buffer>(0, chargeBuffer.getDeviceBuffer());
    setChargesKernel.setArg<cl::Buffer>(1, posq.getDeviceBuffer());
    setChargesKernel.setArg<cl::Buffer>(2, atomIndexDevice.getDeviceBuffer());
    setChargesKernel.setArg<cl_int>(3, numAtoms);
    executeKernel(setChargesKernel, numAtoms);
}

bool OpenCLContext::requestPosqCharges() {
    bool allow = !hasAssignedPosqCharges;
    hasAssignedPosqCharges = true;
    return allow;
}

void OpenCLContext::addEnergyParameterDerivative(const string& param) {
    // See if this parameter has already been registered.
    
    for (int i = 0; i < energyParamDerivNames.size(); i++)
        if (param == energyParamDerivNames[i])
            return;
    energyParamDerivNames.push_back(param);
}

void OpenCLContext::flushQueue() {
    getQueue().flush();
}
