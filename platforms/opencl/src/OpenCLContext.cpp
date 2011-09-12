/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
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
#include "OpenCLForceInfo.h"
#include "OpenCLIntegrationUtilities.h"
#include "OpenCLKernelSources.h"
#include "OpenCLNonbondedUtilities.h"
#include "hilbert.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

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
    std::cerr << "OpenCL internal error: " << errinfo << std::endl;
}

OpenCLContext::OpenCLContext(int numParticles, int deviceIndex, OpenCLPlatform::PlatformData& platformData) :
        time(0.0), platformData(platformData), stepCount(0), computeForceCount(0), posq(NULL), velm(NULL),
        forceBuffers(NULL), longForceBuffer(NULL), energyBuffer(NULL), atomIndex(NULL), integration(NULL),
        nonbonded(NULL), thread(NULL) {
    try {
        contextIndex = platformData.contexts.size();
        std::vector<cl::Platform> platforms;
        cl::Platform::get(&platforms);
        cl_context_properties cprops[] = {CL_CONTEXT_PLATFORM, (cl_context_properties) platforms[0](), 0};
        context = cl::Context(CL_DEVICE_TYPE_ALL, cprops, errorCallback);
        vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();
        const int minThreadBlockSize = 32;
        if (deviceIndex < 0 || deviceIndex >= (int) devices.size()) {
            // Try to figure out which device is the fastest.

            int bestSpeed = -1;
            for (int i = 0; i < (int) devices.size(); i++) {
                int maxSize = devices[i].getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>()[0];
                int processingElementsPerComputeUnit = (devices[i].getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_GPU ? 8 : 1);
                if (devices[i].getInfo<CL_DEVICE_EXTENSIONS>().find("cl_nv_device_attribute_query") != string::npos) {
                    cl_uint computeCapabilityMajor;
                    clGetDeviceInfo(devices[i](), CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV, sizeof(cl_uint), &computeCapabilityMajor, NULL);
                    processingElementsPerComputeUnit = (computeCapabilityMajor < 2 ? 8 : 32);
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
        compilationDefines["WORK_GROUP_SIZE"] = OpenCLExpressionUtilities::intToString(ThreadBlockSize);
        defaultOptimizationOptions = "-cl-fast-relaxed-math";
        supports64BitGlobalAtomics = (device.getInfo<CL_DEVICE_EXTENSIONS>().find("cl_khr_int64_base_atomics") != string::npos);
        string vendor = device.getInfo<CL_DEVICE_VENDOR>();
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
            // AMD APP SDK 2.4 has a performance problem with atomics. Enable the work around.
            compilationDefines["AMD_ATOMIC_WORK_AROUND"] = "";
            // AMD has both 32 and 64 width SIMDs. To determine need to create a kernel to query.
            // For now default to 1 which will use the default kernels.
            simdWidth = 1;
        }
        else
            simdWidth = 1;
        if (supports64BitGlobalAtomics)
            compilationDefines["SUPPORTS_64_BIT_ATOMICS"] = "";
        queue = cl::CommandQueue(context, device);
        numAtoms = numParticles;
        paddedNumAtoms = TileSize*((numParticles+TileSize-1)/TileSize);
        numAtomBlocks = (paddedNumAtoms+(TileSize-1))/TileSize;
        numThreadBlocks = 6*device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
        nonbonded = new OpenCLNonbondedUtilities(*this);
        posq = new OpenCLArray<mm_float4>(*this, paddedNumAtoms, "posq", true);
        velm = new OpenCLArray<mm_float4>(*this, paddedNumAtoms, "velm", true);
        posCellOffsets.resize(paddedNumAtoms, mm_int4(0, 0, 0, 0));
    }
    catch (cl::Error err) {
        std::stringstream str;
        str<<"Error initializing context: "<<err.what()<<" ("<<err.err()<<")";
        throw OpenMMException(str.str());
    }

    // Create utility kernels that are used in multiple places.

    utilities = createProgram(OpenCLKernelSources::utilities);
    clearBufferKernel = cl::Kernel(utilities, "clearBuffer");
    clearTwoBuffersKernel = cl::Kernel(utilities, "clearTwoBuffers");
    clearThreeBuffersKernel = cl::Kernel(utilities, "clearThreeBuffers");
    clearFourBuffersKernel = cl::Kernel(utilities, "clearFourBuffers");
    clearFiveBuffersKernel = cl::Kernel(utilities, "clearFiveBuffers");
    clearSixBuffersKernel = cl::Kernel(utilities, "clearSixBuffers");
    reduceFloat4Kernel = cl::Kernel(utilities, "reduceFloat4Buffer");
    reduceForcesKernel = cl::Kernel(utilities, "reduceForces");

    // Decide whether native_sqrt(), native_rsqrt(), and native_recip() are sufficiently accurate to use.

    cl::Kernel accuracyKernel(utilities, "determineNativeAccuracy");
    OpenCLArray<mm_float8> values(*this, 20, "values", true);
    float nextValue = 1e-4f;
    for (int i = 0; i < values.getSize(); ++i) {
        values[i].s0 = nextValue;
        nextValue *= (float) M_PI;
    }
    values.upload();
    accuracyKernel.setArg<cl::Buffer>(0, values.getDeviceBuffer());
    accuracyKernel.setArg<cl_int>(1, values.getSize());
    executeKernel(accuracyKernel, values.getSize());
    values.download();
    double maxSqrtError = 0.0, maxRsqrtError = 0.0, maxRecipError = 0.0, maxExpError = 0.0, maxLogError = 0.0;
    for (int i = 0; i < values.getSize(); ++i) {
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
    
    // Create the work thread used for parallelization when running on multiple devices.
    
    thread = new WorkThread();
}

OpenCLContext::~OpenCLContext() {
    for (int i = 0; i < (int) forces.size(); i++)
        delete forces[i];
    if (posq != NULL)
        delete posq;
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
    if (atomIndex != NULL)
        delete atomIndex;
    if (integration != NULL)
        delete integration;
    if (nonbonded != NULL)
        delete nonbonded;
    if (thread != NULL)
        delete thread;
}

void OpenCLContext::initialize(const System& system) {
    for (int i = 0; i < numAtoms; i++)
        (*velm)[i].w = (float) (1.0/system.getParticleMass(i));
    velm->upload();
    numForceBuffers = platformData.contexts.size();
    for (int i = 0; i < (int) forces.size(); i++)
        numForceBuffers = std::max(numForceBuffers, forces[i]->getRequiredForceBuffers());
    forceBuffers = new OpenCLArray<mm_float4>(*this, paddedNumAtoms*numForceBuffers, "forceBuffers", false);
    if (supports64BitGlobalAtomics) {
        longForceBuffer = new OpenCLArray<cl_long>(*this, 3*paddedNumAtoms, "longForceBuffer", false);
        reduceForcesKernel.setArg<cl::Buffer>(0, longForceBuffer->getDeviceBuffer());
        reduceForcesKernel.setArg<cl::Buffer>(1, forceBuffers->getDeviceBuffer());
        reduceForcesKernel.setArg<cl_int>(2, paddedNumAtoms);
        reduceForcesKernel.setArg<cl_int>(3, numForceBuffers);
        addAutoclearBuffer(longForceBuffer->getDeviceBuffer(), longForceBuffer->getSize()*2);
    }
    addAutoclearBuffer(forceBuffers->getDeviceBuffer(), forceBuffers->getSize()*4);
    force = new OpenCLArray<mm_float4>(*this, &forceBuffers->getDeviceBuffer(), paddedNumAtoms, "force", true);
    energyBuffer = new OpenCLArray<cl_float>(*this, max(numThreadBlocks*ThreadBlockSize, nonbonded->getNumEnergyBuffers()), "energyBuffer", true);
    addAutoclearBuffer(energyBuffer->getDeviceBuffer(), energyBuffer->getSize());
    atomIndex = new OpenCLArray<cl_int>(*this, paddedNumAtoms, "atomIndex", true);
    for (int i = 0; i < paddedNumAtoms; ++i)
        (*atomIndex)[i] = i;
    atomIndex->upload();
    findMoleculeGroups(system);
    integration = new OpenCLIntegrationUtilities(*this, system);
    nonbonded->initialize(system);
}

void OpenCLContext::addForce(OpenCLForceInfo* force) {
    forces.push_back(force);
}

string OpenCLContext::loadSourceFromFile(const string& filename) const {
    ifstream file((Platform::getDefaultPluginsDirectory()+"/opencl/"+filename).c_str());
    if (!file.is_open())
        throw OpenMMException("Unable to load kernel: "+filename);
    string kernel;
    string line;
    while (!file.eof()) {
        getline(file, line);
        kernel += line;
        kernel += '\n';
    }
    file.close();
    return kernel;
}

string OpenCLContext::loadSourceFromFile(const string& filename, const std::map<std::string, std::string>& replacements) const {
    return replaceStrings(loadSourceFromFile(filename), replacements);
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
    string options = (optimizationFlags == NULL ? defaultOptimizationOptions : optimizationFlags);
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

void OpenCLContext::clearBuffer(OpenCLArray<float>& array) {
    clearBuffer(array.getDeviceBuffer(), array.getSize());
}

void OpenCLContext::clearBuffer(OpenCLArray<mm_float4>& array) {
    clearBuffer(array.getDeviceBuffer(), array.getSize()*4);
}

void OpenCLContext::clearBuffer(cl::Memory& memory, int size) {
    clearBufferKernel.setArg<cl::Memory>(0, memory);
    clearBufferKernel.setArg<cl_int>(1, size);
    executeKernel(clearBufferKernel, size, 128);
}

void OpenCLContext::addAutoclearBuffer(cl::Memory& memory, int size) {
    autoclearBuffers.push_back(&memory);
    autoclearBufferSizes.push_back(size);
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

void OpenCLContext::reduceBuffer(OpenCLArray<mm_float4>& array, int numBuffers) {
    int bufferSize = array.getSize()/numBuffers;
    reduceFloat4Kernel.setArg<cl::Buffer>(0, array.getDeviceBuffer());
    reduceFloat4Kernel.setArg<cl_int>(1, bufferSize);
    reduceFloat4Kernel.setArg<cl_int>(2, numBuffers);
    executeKernel(reduceFloat4Kernel, bufferSize, 128);
}

void OpenCLContext::tagAtomsInMolecule(int atom, int molecule, vector<int>& atomMolecule, vector<vector<int> >& atomBonds) {
    // Recursively tag atoms as belonging to a particular molecule.

    atomMolecule[atom] = molecule;
    for (int i = 0; i < (int) atomBonds[atom].size(); i++)
        if (atomMolecule[atomBonds[atom][i]] == -1)
            tagAtomsInMolecule(atomBonds[atom][i], molecule, atomMolecule, atomBonds);
}

struct OpenCLContext::Molecule {
    vector<int> atoms;
    vector<int> constraints;
    vector<vector<int> > groups;
};

void OpenCLContext::findMoleculeGroups(const System& system) {
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

    vector<Molecule> molecules(numMolecules);
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

    // Sort them into groups of identical molecules.

    vector<Molecule> uniqueMolecules;
    vector<vector<int> > moleculeInstances;
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
                moleculeInstances[j].push_back(mol.atoms[0]);
                isNew = false;
            }
        }
        if (isNew) {
            uniqueMolecules.push_back(mol);
            moleculeInstances.push_back(vector<int>());
            moleculeInstances[moleculeInstances.size()-1].push_back(mol.atoms[0]);
        }
    }
    moleculeGroups.resize(moleculeInstances.size());
    for (int i = 0; i < (int) moleculeInstances.size(); i++)
    {
        moleculeGroups[i].instances = moleculeInstances[i];
        vector<int>& atoms = uniqueMolecules[i].atoms;
        moleculeGroups[i].atoms.resize(atoms.size());
        for (int j = 0; j < (int) atoms.size(); j++)
            moleculeGroups[i].atoms[j] = atoms[j]-atoms[0];
    }
}

void OpenCLContext::reorderAtoms() {
    if (numAtoms == 0 || nonbonded == NULL || !nonbonded->getUseCutoff())
        return;

    // Find the range of positions and the number of bins along each axis.

    posq->download();
    velm->download();
    float minx = posq->get(0).x, maxx = posq->get(0).x;
    float miny = posq->get(0).y, maxy = posq->get(0).y;
    float minz = posq->get(0).z, maxz = posq->get(0).z;
    if (nonbonded->getUsePeriodic()) {
        minx = miny = minz = 0.0;
        maxx = periodicBoxSize.x;
        maxy = periodicBoxSize.y;
        maxz = periodicBoxSize.z;
    }
    else {
        for (int i = 1; i < numAtoms; i++) {
            minx = min(minx, posq->get(i).x);
            maxx = max(maxx, posq->get(i).x);
            miny = min(miny, posq->get(i).y);
            maxy = max(maxy, posq->get(i).y);
            minz = min(minz, posq->get(i).z);
            maxz = max(maxz, posq->get(i).z);
        }
    }

    // Loop over each group of identical molecules and reorder them.

    vector<int> originalIndex(numAtoms);
    vector<mm_float4> newPosq(numAtoms);
    vector<mm_float4> newVelm(numAtoms);
    vector<mm_int4> newCellOffsets(numAtoms);
    for (int group = 0; group < (int) moleculeGroups.size(); group++) {
        // Find the center of each molecule.

        MoleculeGroup& mol = moleculeGroups[group];
        int numMolecules = mol.instances.size();
        vector<int>& atoms = mol.atoms;
        vector<mm_float4> molPos(numMolecules);
        for (int i = 0; i < numMolecules; i++) {
            molPos[i].x = 0.0f;
            molPos[i].y = 0.0f;
            molPos[i].z = 0.0f;
            for (int j = 0; j < (int)atoms.size(); j++) {
                int atom = atoms[j]+mol.instances[i];
                molPos[i].x += posq->get(atom).x;
                molPos[i].y += posq->get(atom).y;
                molPos[i].z += posq->get(atom).z;
            }
            molPos[i].x /= atoms.size();
            molPos[i].y /= atoms.size();
            molPos[i].z /= atoms.size();
        }
        if (nonbonded->getUsePeriodic()) {
            // Move each molecule position into the same box.

            for (int i = 0; i < numMolecules; i++) {
                int xcell = (int) floor(molPos[i].x/periodicBoxSize.x);
                int ycell = (int) floor(molPos[i].y/periodicBoxSize.y);
                int zcell = (int) floor(molPos[i].z/periodicBoxSize.z);
                float dx = xcell*periodicBoxSize.x;
                float dy = ycell*periodicBoxSize.y;
                float dz = zcell*periodicBoxSize.z;
                if (dx != 0.0f || dy != 0.0f || dz != 0.0f) {
                    molPos[i].x -= dx;
                    molPos[i].y -= dy;
                    molPos[i].z -= dz;
                    for (int j = 0; j < (int) atoms.size(); j++) {
                        int atom = atoms[j]+mol.instances[i];
                        mm_float4 p = posq->get(atom);
                        p.x -= dx;
                        p.y -= dy;
                        p.z -= dz;
                        posq->set(atom, p);
                        posCellOffsets[atom].x -= xcell;
                        posCellOffsets[atom].y -= ycell;
                        posCellOffsets[atom].z -= zcell;
                    }
                }
            }
        }

        // Select a bin for each molecule, then sort them by bin.

        bool useHilbert = (numMolecules > 5000 || atoms.size() > 8); // For small systems, a simple zigzag curve works better than a Hilbert curve.
        float binWidth;
        if (useHilbert)
            binWidth = (float)(max(max(maxx-minx, maxy-miny), maxz-minz)/255.0);
        else
            binWidth = (float)(0.2*nonbonded->getCutoffDistance());
        int xbins = 1 + (int) ((maxx-minx)/binWidth);
        int ybins = 1 + (int) ((maxy-miny)/binWidth);
        vector<pair<int, int> > molBins(numMolecules);
        bitmask_t coords[3];
        for (int i = 0; i < numMolecules; i++) {
            int x = (int) ((molPos[i].x-minx)/binWidth);
            int y = (int) ((molPos[i].y-miny)/binWidth);
            int z = (int) ((molPos[i].z-minz)/binWidth);
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
                int oldIndex = mol.instances[molBins[i].second]+atoms[j];
                int newIndex = mol.instances[i]+atoms[j];
                originalIndex[newIndex] = atomIndex->get(oldIndex);
                newPosq[newIndex] = posq->get(oldIndex);
                newVelm[newIndex] = velm->get(oldIndex);
                newCellOffsets[newIndex] = posCellOffsets[oldIndex];
            }
        }
    }

    // Update the streams.

    for (int i = 0; i < numAtoms; i++) {
        posq->set(i, newPosq[i]);
        velm->set(i, newVelm[i]);
        atomIndex->set(i, originalIndex[i]);
        posCellOffsets[i] = newCellOffsets[i];
    }
    posq->upload();
    velm->upload();
    atomIndex->upload();
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
