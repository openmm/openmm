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
#include "CudaContext.h"
#include "CudaArray.h"
//#include "CudaBondedUtilities.h"
#include "CudaExpressionUtilities.h"
#include "CudaForceInfo.h"
//#include "CudaIntegrationUtilities.h"
#include "CudaKernelSources.h"
//#include "CudaNonbondedUtilities.h"
#include "hilbert.h"
#include "openmm/OpenMMException.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/VirtualSite.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <typeinfo>


#define CHECK_RESULT(result) CHECK_RESULT2(result, errorMessage);
#define CHECK_RESULT2(result, prefix) \
    if (result != CUDA_SUCCESS) { \
        std::stringstream m; \
        m<<prefix<<": "<<result<<" ("<<__FILE__<<": "<<__LINE__<<")"; \
        throw OpenMMException(m.str());\
    }

using namespace OpenMM;
using namespace std;

const int CudaContext::ThreadBlockSize = 64;
const int CudaContext::TileSize = 32;
bool CudaContext::hasInitializedCuda = false;

CudaContext::CudaContext(const System& system, int deviceIndex, bool useBlockingSync, const string& precision, const string& compiler,
        const string& tempDir, CudaPlatform::PlatformData& platformData) : system(system), compiler(compiler),
        time(0.0), platformData(platformData), stepCount(0), computeForceCount(0), contextIsValid(false), atomsWereReordered(false), posq(NULL),
        velm(NULL), /*forceBuffers(NULL), longForceBuffer(NULL), energyBuffer(NULL), atomIndex(NULL), integration(NULL),
        bonded(NULL), nonbonded(NULL),*/ thread(NULL) {
    if (!hasInitializedCuda) {
        CHECK_RESULT2(cuInit(0), "Error initializing CUDA");
        hasInitializedCuda = true;
    }
    if (precision == "single") {
        useDoublePrecision = false;
        accumulateInDouble = false;
    }
    else if (precision == "mixed") {
        useDoublePrecision = false;
        accumulateInDouble = true;
    }
    else if (precision == "double") {
        useDoublePrecision = true;
        accumulateInDouble = true;
    }
    else
        throw OpenMMException("Illegal value for CudaPrecision: "+precision);
#ifdef WIN32
    this->tempDir = tempDir+"\";
#else
    this->tempDir = tempDir+"/";
#endif
    contextIndex = platformData.contexts.size();
    int numDevices;
    string errorMessage = "Error initializing Context";
    CHECK_RESULT(cuDeviceGetCount(&numDevices));
    if (deviceIndex < 0 || deviceIndex >= numDevices) {
        // Try to figure out which device is the fastest.

        int bestSpeed = -1;
        int bestCompute = -1;
        for (int i = 0; i < numDevices; i++) {
            CHECK_RESULT(cuDeviceGet(&device, i));
            int major, minor, clock, multiprocessors;
            CHECK_RESULT(cuDeviceComputeCapability(&major, &minor, device));
            if (major == 1 && minor < 2)
                continue; // 1.0 and 1.1 are not supported
            CHECK_RESULT(cuDeviceGetAttribute(&clock, CU_DEVICE_ATTRIBUTE_CLOCK_RATE, device));
            CHECK_RESULT(cuDeviceGetAttribute(&multiprocessors, CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT, device));
            int speed = clock*multiprocessors;
            if (major > bestCompute || (major == bestCompute && speed > bestSpeed)) {
                deviceIndex = i;
                bestSpeed = speed;
                bestCompute = major;
            }
        }
    }
    if (deviceIndex == -1)
        throw OpenMMException("No compatible CUDA device is available");
    CHECK_RESULT(cuDeviceGet(&device, deviceIndex));
    this->deviceIndex = deviceIndex;
    int major, minor;
    CHECK_RESULT(cuDeviceComputeCapability(&major, &minor, device));
    gpuArchitecture = CudaExpressionUtilities::intToString(major)+CudaExpressionUtilities::intToString(minor);
    compilationDefines["WORK_GROUP_SIZE"] = CudaExpressionUtilities::intToString(ThreadBlockSize);
    defaultOptimizationOptions = "--use_fast_math";
    int numThreadBlocksPerComputeUnit = 6;
    CHECK_RESULT(cuCtxCreate(&context, 0, device));
    contextIsValid = true;
    numAtoms = system.getNumParticles();
    paddedNumAtoms = TileSize*((numAtoms+TileSize-1)/TileSize);
    numAtomBlocks = (paddedNumAtoms+(TileSize-1))/TileSize;
    int multiprocessors;
    CHECK_RESULT(cuDeviceGetAttribute(&multiprocessors, CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT, device));
    numThreadBlocks = numThreadBlocksPerComputeUnit*multiprocessors;
//    bonded = new CudaBondedUtilities(*this);
//    nonbonded = new CudaNonbondedUtilities(*this);
    posq = CudaArray::create<float4>(paddedNumAtoms, "posq");
    velm = CudaArray::create<float4>(paddedNumAtoms, "velm");
    posCellOffsets.resize(paddedNumAtoms, make_int4(0, 0, 0, 0));

    // Create utility kernels that are used in multiple places.

    CUmodule utilities = createModule(CudaKernelSources::vectorOps+CudaKernelSources::utilities);
    cuModuleGetFunction(&clearBufferKernel, utilities, "clearBuffer");
    cuModuleGetFunction(&clearTwoBuffersKernel, utilities, "clearTwoBuffers");
    cuModuleGetFunction(&clearThreeBuffersKernel, utilities, "clearThreeBuffers");
    cuModuleGetFunction(&clearFourBuffersKernel, utilities, "clearFourBuffers");
    cuModuleGetFunction(&clearFiveBuffersKernel, utilities, "clearFiveBuffers");
    cuModuleGetFunction(&clearSixBuffersKernel, utilities, "clearSixBuffers");
    cuModuleGetFunction(&reduceFloat4Kernel, utilities, "reduceFloat4Buffer");
    cuModuleGetFunction(&reduceForcesKernel, utilities, "reduceForces");

    // Set defines based on the requested precision.

    compilationDefines["SQRT"] = useDoublePrecision ? "sqrt" : "sqrtf";
    compilationDefines["RSQRT"] = useDoublePrecision ? "rsqrt" : "rsqrtf";
    compilationDefines["RECIP"] = useDoublePrecision ? "1.0/" : "1.0f/";
    compilationDefines["EXP"] = useDoublePrecision ? "exp" : "expf";
    compilationDefines["LOG"] = useDoublePrecision ? "log" : "logf";
    
    // Create the work thread used for parallelization when running on multiple devices.
    
    thread = new WorkThread();
//    
//    // Create the integration utilities object.
//    
//    integration = new CudaIntegrationUtilities(*this, system);
}

CudaContext::~CudaContext() {
    for (int i = 0; i < (int) forces.size(); i++)
        delete forces[i];
    for (int i = 0; i < (int) reorderListeners.size(); i++)
        delete reorderListeners[i];
    if (posq != NULL)
        delete posq;
    if (velm != NULL)
        delete velm;
//    if (force != NULL)
//        delete force;
//    if (forceBuffers != NULL)
//        delete forceBuffers;
//    if (longForceBuffer != NULL)
//        delete longForceBuffer;
//    if (energyBuffer != NULL)
//        delete energyBuffer;
//    if (atomIndex != NULL)
//        delete atomIndex;
//    if (integration != NULL)
//        delete integration;
//    if (bonded != NULL)
//        delete bonded;
//    if (nonbonded != NULL)
//        delete nonbonded;
    if (thread != NULL)
        delete thread;
    string errorMessage = "Error deleting Context";
    if (contextIsValid)
        CHECK_RESULT(cuCtxDestroy(context));
}

//void CudaContext::initialize() {
//    for (int i = 0; i < numAtoms; i++) {
//        double mass = system.getParticleMass(i);
//        (*velm)[i].w = (float) (mass == 0.0 ? 0.0 : 1.0/mass);
//    }
//    velm->upload();
//    bonded->initialize(system);
//    numForceBuffers = platformData.contexts.size();
//    numForceBuffers = std::max(numForceBuffers, bonded->getNumForceBuffers());
//    for (int i = 0; i < (int) forces.size(); i++)
//        numForceBuffers = std::max(numForceBuffers, forces[i]->getRequiredForceBuffers());
//    forceBuffers = new CudaArray<mm_float4>(*this, paddedNumAtoms*numForceBuffers, "forceBuffers", false);
//    if (supports64BitGlobalAtomics) {
//        longForceBuffer = new CudaArray<cl_long>(*this, 3*paddedNumAtoms, "longForceBuffer", false);
//        reduceForcesKernel.setArg<cl::Buffer>(0, longForceBuffer->getDeviceBuffer());
//        reduceForcesKernel.setArg<cl::Buffer>(1, forceBuffers->getDeviceBuffer());
//        reduceForcesKernel.setArg<cl_int>(2, paddedNumAtoms);
//        reduceForcesKernel.setArg<cl_int>(3, numForceBuffers);
//        addAutoclearBuffer(longForceBuffer->getDeviceBuffer(), longForceBuffer->getSize()*2);
//    }
//    addAutoclearBuffer(forceBuffers->getDeviceBuffer(), forceBuffers->getSize()*4);
//    force = new CudaArray<mm_float4>(*this, &forceBuffers->getDeviceBuffer(), paddedNumAtoms, "force", true);
//    energyBuffer = new CudaArray<cl_float>(*this, max(numThreadBlocks*ThreadBlockSize, nonbonded->getNumEnergyBuffers()), "energyBuffer", true);
//    addAutoclearBuffer(energyBuffer->getDeviceBuffer(), energyBuffer->getSize());
//    atomIndex = new CudaArray<cl_int>(*this, paddedNumAtoms, "atomIndex", true);
//    for (int i = 0; i < paddedNumAtoms; ++i)
//        (*atomIndex)[i] = i;
//    atomIndex->upload();
//    findMoleculeGroups();
//    moleculesInvalid = false;
//    nonbonded->initialize(system);
//}

void CudaContext::addForce(CudaForceInfo* force) {
    forces.push_back(force);
}

string CudaContext::replaceStrings(const string& input, const std::map<std::string, std::string>& replacements) const {
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

CUmodule CudaContext::createModule(const string source, const char* optimizationFlags) {
    return createModule(source, map<string, string>(), optimizationFlags);
}

CUmodule CudaContext::createModule(const string source, const map<string, string>& defines, const char* optimizationFlags) {
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
    for (map<string, string>::const_iterator iter = defines.begin(); iter != defines.end(); ++iter) {
        src << "#define " << iter->first;
        if (!iter->second.empty())
            src << " " << iter->second;
        src << endl;
    }
    if (!defines.empty())
        src << endl;
    src << source << endl;
    
    // Write out the source to a temporary file.
    
    stringstream tempFileName;
    tempFileName << "openmmTempKernel" << this; // Include a pointer to this context as part of the filename to avoid collisions.
    string inputFile = (tempDir+tempFileName.str()+".cu");
    string outputFile = (tempDir+tempFileName.str()+".ptx");
    string logFile = (tempDir+tempFileName.str()+".log");
    ofstream out(inputFile.c_str());
    out << src.str();
    out.close();
#ifdef WIN32
#else
    string command = "\""+compiler+"\" --ptx -arch=compute_"+gpuArchitecture+" -o \""+outputFile+"\" "+options+" \""+inputFile+"\" 2> \""+logFile+"\"";
    int res = std::system(command.c_str());
#endif
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
        CUmodule module;
        CUresult result = cuModuleLoad(&module, outputFile.c_str());
        if (result != CUDA_SUCCESS) {
            std::stringstream m;
            m<<"Error loading CUDA module: "<<result;
            throw OpenMMException(m.str());
        }
        remove(inputFile.c_str());
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
//    
//    // Get length before using c_str() to avoid length() call invalidating the c_str() value.
//    string src_string = src.str();
//    ::size_t src_length = src_string.length();
//    cl::Program::Sources sources(1, make_pair(src_string.c_str(), src_length));
//    cl::Program program(context, sources);
//    try {
//        program.build(vector<cl::Device>(1, device), options.c_str());
//    } catch (cl::Error err) {
//        throw OpenMMException("Error compiling kernel: "+program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device));
//    }
}
//
//void CudaContext::executeKernel(cl::Kernel& kernel, int workUnits, int blockSize) {
//    if (blockSize == -1)
//        blockSize = ThreadBlockSize;
//    int size = std::min((workUnits+blockSize-1)/blockSize, numThreadBlocks)*blockSize;
//    try {
//        queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(size), cl::NDRange(blockSize));
//    }
//    catch (cl::Error err) {
//        stringstream str;
//        str<<"Error invoking kernel "<<kernel.getInfo<CL_KERNEL_FUNCTION_NAME>()<<": "<<err.what()<<" ("<<err.err()<<")";
//        throw OpenMMException(str.str());
//    }
//}
//
//void CudaContext::clearBuffer(CudaArray<float>& array) {
//    clearBuffer(array.getDeviceBuffer(), array.getSize());
//}
//
//void CudaContext::clearBuffer(CudaArray<mm_float4>& array) {
//    clearBuffer(array.getDeviceBuffer(), array.getSize()*4);
//}
//
//void CudaContext::clearBuffer(cl::Memory& memory, int size) {
//    clearBufferKernel.setArg<cl::Memory>(0, memory);
//    clearBufferKernel.setArg<cl_int>(1, size);
//    executeKernel(clearBufferKernel, size, 128);
//}
//
//void CudaContext::addAutoclearBuffer(cl::Memory& memory, int size) {
//    autoclearBuffers.push_back(&memory);
//    autoclearBufferSizes.push_back(size);
//}
//
//void CudaContext::clearAutoclearBuffers() {
//    int base = 0;
//    int total = autoclearBufferSizes.size();
//    while (total-base >= 6) {
//        clearSixBuffersKernel.setArg<cl::Memory>(0, *autoclearBuffers[base]);
//        clearSixBuffersKernel.setArg<cl_int>(1, autoclearBufferSizes[base]);
//        clearSixBuffersKernel.setArg<cl::Memory>(2, *autoclearBuffers[base+1]);
//        clearSixBuffersKernel.setArg<cl_int>(3, autoclearBufferSizes[base+1]);
//        clearSixBuffersKernel.setArg<cl::Memory>(4, *autoclearBuffers[base+2]);
//        clearSixBuffersKernel.setArg<cl_int>(5, autoclearBufferSizes[base+2]);
//        clearSixBuffersKernel.setArg<cl::Memory>(6, *autoclearBuffers[base+3]);
//        clearSixBuffersKernel.setArg<cl_int>(7, autoclearBufferSizes[base+3]);
//        clearSixBuffersKernel.setArg<cl::Memory>(8, *autoclearBuffers[base+4]);
//        clearSixBuffersKernel.setArg<cl_int>(9, autoclearBufferSizes[base+4]);
//        clearSixBuffersKernel.setArg<cl::Memory>(10, *autoclearBuffers[base+5]);
//        clearSixBuffersKernel.setArg<cl_int>(11, autoclearBufferSizes[base+5]);
//        executeKernel(clearSixBuffersKernel, max(max(max(max(max(autoclearBufferSizes[base], autoclearBufferSizes[base+1]), autoclearBufferSizes[base+2]), autoclearBufferSizes[base+3]), autoclearBufferSizes[base+4]), autoclearBufferSizes[base+5]), 128);
//        base += 6;
//    }
//    if (total-base == 5) {
//        clearFiveBuffersKernel.setArg<cl::Memory>(0, *autoclearBuffers[base]);
//        clearFiveBuffersKernel.setArg<cl_int>(1, autoclearBufferSizes[base]);
//        clearFiveBuffersKernel.setArg<cl::Memory>(2, *autoclearBuffers[base+1]);
//        clearFiveBuffersKernel.setArg<cl_int>(3, autoclearBufferSizes[base+1]);
//        clearFiveBuffersKernel.setArg<cl::Memory>(4, *autoclearBuffers[base+2]);
//        clearFiveBuffersKernel.setArg<cl_int>(5, autoclearBufferSizes[base+2]);
//        clearFiveBuffersKernel.setArg<cl::Memory>(6, *autoclearBuffers[base+3]);
//        clearFiveBuffersKernel.setArg<cl_int>(7, autoclearBufferSizes[base+3]);
//        clearFiveBuffersKernel.setArg<cl::Memory>(8, *autoclearBuffers[base+4]);
//        clearFiveBuffersKernel.setArg<cl_int>(9, autoclearBufferSizes[base+4]);
//        executeKernel(clearFiveBuffersKernel, max(max(max(max(autoclearBufferSizes[base], autoclearBufferSizes[base+1]), autoclearBufferSizes[base+2]), autoclearBufferSizes[base+3]), autoclearBufferSizes[base+4]), 128);
//    }
//    else if (total-base == 4) {
//        clearFourBuffersKernel.setArg<cl::Memory>(0, *autoclearBuffers[base]);
//        clearFourBuffersKernel.setArg<cl_int>(1, autoclearBufferSizes[base]);
//        clearFourBuffersKernel.setArg<cl::Memory>(2, *autoclearBuffers[base+1]);
//        clearFourBuffersKernel.setArg<cl_int>(3, autoclearBufferSizes[base+1]);
//        clearFourBuffersKernel.setArg<cl::Memory>(4, *autoclearBuffers[base+2]);
//        clearFourBuffersKernel.setArg<cl_int>(5, autoclearBufferSizes[base+2]);
//        clearFourBuffersKernel.setArg<cl::Memory>(6, *autoclearBuffers[base+3]);
//        clearFourBuffersKernel.setArg<cl_int>(7, autoclearBufferSizes[base+3]);
//        executeKernel(clearFourBuffersKernel, max(max(max(autoclearBufferSizes[base], autoclearBufferSizes[base+1]), autoclearBufferSizes[base+2]), autoclearBufferSizes[base+3]), 128);
//    }
//    else if (total-base == 3) {
//        clearThreeBuffersKernel.setArg<cl::Memory>(0, *autoclearBuffers[base]);
//        clearThreeBuffersKernel.setArg<cl_int>(1, autoclearBufferSizes[base]);
//        clearThreeBuffersKernel.setArg<cl::Memory>(2, *autoclearBuffers[base+1]);
//        clearThreeBuffersKernel.setArg<cl_int>(3, autoclearBufferSizes[base+1]);
//        clearThreeBuffersKernel.setArg<cl::Memory>(4, *autoclearBuffers[base+2]);
//        clearThreeBuffersKernel.setArg<cl_int>(5, autoclearBufferSizes[base+2]);
//        executeKernel(clearThreeBuffersKernel, max(max(autoclearBufferSizes[base], autoclearBufferSizes[base+1]), autoclearBufferSizes[base+2]), 128);
//    }
//    else if (total-base == 2) {
//        clearTwoBuffersKernel.setArg<cl::Memory>(0, *autoclearBuffers[base]);
//        clearTwoBuffersKernel.setArg<cl_int>(1, autoclearBufferSizes[base]);
//        clearTwoBuffersKernel.setArg<cl::Memory>(2, *autoclearBuffers[base+1]);
//        clearTwoBuffersKernel.setArg<cl_int>(3, autoclearBufferSizes[base+1]);
//        executeKernel(clearTwoBuffersKernel, max(autoclearBufferSizes[base], autoclearBufferSizes[base+1]), 128);
//    }
//    else if (total-base == 1) {
//        clearBuffer(*autoclearBuffers[base], autoclearBufferSizes[base]);
//    }
//}
//
//void CudaContext::reduceForces() {
//    if (supports64BitGlobalAtomics)
//        executeKernel(reduceForcesKernel, paddedNumAtoms, 128);
//    else
//        reduceBuffer(*forceBuffers, numForceBuffers);
//}
//
//void CudaContext::reduceBuffer(CudaArray<mm_float4>& array, int numBuffers) {
//    int bufferSize = array.getSize()/numBuffers;
//    reduceFloat4Kernel.setArg<cl::Buffer>(0, array.getDeviceBuffer());
//    reduceFloat4Kernel.setArg<cl_int>(1, bufferSize);
//    reduceFloat4Kernel.setArg<cl_int>(2, numBuffers);
//    executeKernel(reduceFloat4Kernel, bufferSize, 128);
//}
//
//void CudaContext::tagAtomsInMolecule(int atom, int molecule, vector<int>& atomMolecule, vector<vector<int> >& atomBonds) {
//    // Recursively tag atoms as belonging to a particular molecule.
//
//    atomMolecule[atom] = molecule;
//    for (int i = 0; i < (int) atomBonds[atom].size(); i++)
//        if (atomMolecule[atomBonds[atom][i]] == -1)
//            tagAtomsInMolecule(atomBonds[atom][i], molecule, atomMolecule, atomBonds);
//}
//
///**
// * This class ensures that atom reordering doesn't break virtual sites.
// */
//class CudaContext::VirtualSiteInfo : public CudaForceInfo {
//public:
//    VirtualSiteInfo(const System& system) : CudaForceInfo(0) {
//        for (int i = 0; i < system.getNumParticles(); i++) {
//            if (system.isVirtualSite(i)) {
//                siteTypes.push_back(&typeid(system.getVirtualSite(i)));
//                vector<int> particles;
//                particles.push_back(i);
//                for (int j = 0; j < system.getVirtualSite(i).getNumParticles(); j++)
//                    particles.push_back(system.getVirtualSite(i).getParticle(j));
//                siteParticles.push_back(particles);
//                vector<double> weights;
//                if (dynamic_cast<const TwoParticleAverageSite*>(&system.getVirtualSite(i)) != NULL) {
//                    // A two particle average.
//
//                    const TwoParticleAverageSite& site = dynamic_cast<const TwoParticleAverageSite&>(system.getVirtualSite(i));
//                    weights.push_back(site.getWeight(0));
//                    weights.push_back(site.getWeight(1));
//                }
//                else if (dynamic_cast<const ThreeParticleAverageSite*>(&system.getVirtualSite(i)) != NULL) {
//                    // A three particle average.
//
//                    const ThreeParticleAverageSite& site = dynamic_cast<const ThreeParticleAverageSite&>(system.getVirtualSite(i));
//                    weights.push_back(site.getWeight(0));
//                    weights.push_back(site.getWeight(1));
//                    weights.push_back(site.getWeight(2));
//                }
//                else if (dynamic_cast<const OutOfPlaneSite*>(&system.getVirtualSite(i)) != NULL) {
//                    // An out of plane site.
//
//                    const OutOfPlaneSite& site = dynamic_cast<const OutOfPlaneSite&>(system.getVirtualSite(i));
//                    weights.push_back(site.getWeight12());
//                    weights.push_back(site.getWeight13());
//                    weights.push_back(site.getWeightCross());
//                }
//                siteWeights.push_back(weights);
//            }
//        }
//    }
//    int getNumParticleGroups() {
//        return siteTypes.size();
//    }
//    void getParticlesInGroup(int index, std::vector<int>& particles) {
//        particles = siteParticles[index];
//    }
//    bool areGroupsIdentical(int group1, int group2) {
//        if (siteTypes[group1] != siteTypes[group2])
//            return false;
//        int numParticles = siteWeights[group1].size();
//        if (siteWeights[group2].size() != numParticles)
//            return false;
//        for (int i = 0; i < numParticles; i++)
//            if (siteWeights[group1][i] != siteWeights[group2][i])
//                return false;
//        return true;
//    }
//private:
//    vector<const type_info*> siteTypes;
//    vector<vector<int> > siteParticles;
//    vector<vector<double> > siteWeights;
//};
//
//
//void CudaContext::findMoleculeGroups() {
//    // The first time this is called, we need to identify all the molecules in the system.
//    
//    if (moleculeGroups.size() == 0) {
//        // Add a ForceInfo that makes sure reordering doesn't break virtual sites.
//
//        addForce(new VirtualSiteInfo(system));
//
//        // First make a list of every other atom to which each atom is connect by a constraint or force group.
//
//        vector<vector<int> > atomBonds(system.getNumParticles());
//        for (int i = 0; i < system.getNumConstraints(); i++) {
//            int particle1, particle2;
//            double distance;
//            system.getConstraintParameters(i, particle1, particle2, distance);
//            atomBonds[particle1].push_back(particle2);
//            atomBonds[particle2].push_back(particle1);
//        }
//        for (int i = 0; i < (int) forces.size(); i++) {
//            for (int j = 0; j < forces[i]->getNumParticleGroups(); j++) {
//                vector<int> particles;
//                forces[i]->getParticlesInGroup(j, particles);
//                for (int k = 0; k < (int) particles.size(); k++)
//                    for (int m = 0; m < (int) particles.size(); m++)
//                        if (k != m)
//                            atomBonds[particles[k]].push_back(particles[m]);
//            }
//        }
//
//        // Now tag atoms by which molecule they belong to.
//
//        vector<int> atomMolecule(numAtoms, -1);
//        int numMolecules = 0;
//        for (int i = 0; i < numAtoms; i++)
//            if (atomMolecule[i] == -1)
//                tagAtomsInMolecule(i, numMolecules++, atomMolecule, atomBonds);
//        vector<vector<int> > atomIndices(numMolecules);
//        for (int i = 0; i < numAtoms; i++)
//            atomIndices[atomMolecule[i]].push_back(i);
//
//        // Construct a description of each molecule.
//
//        molecules.resize(numMolecules);
//        for (int i = 0; i < numMolecules; i++) {
//            molecules[i].atoms = atomIndices[i];
//            molecules[i].groups.resize(forces.size());
//        }
//        for (int i = 0; i < system.getNumConstraints(); i++) {
//            int particle1, particle2;
//            double distance;
//            system.getConstraintParameters(i, particle1, particle2, distance);
//            molecules[atomMolecule[particle1]].constraints.push_back(i);
//        }
//        for (int i = 0; i < (int) forces.size(); i++)
//            for (int j = 0; j < forces[i]->getNumParticleGroups(); j++) {
//                vector<int> particles;
//                forces[i]->getParticlesInGroup(j, particles);
//                molecules[atomMolecule[particles[0]]].groups[i].push_back(j);
//            }
//    }
//
//    // Sort them into groups of identical molecules.
//
//    vector<Molecule> uniqueMolecules;
//    vector<vector<int> > moleculeInstances;
//    vector<vector<int> > moleculeOffsets;
//    for (int molIndex = 0; molIndex < (int) molecules.size(); molIndex++) {
//        Molecule& mol = molecules[molIndex];
//
//        // See if it is identical to another molecule.
//
//        bool isNew = true;
//        for (int j = 0; j < (int) uniqueMolecules.size() && isNew; j++) {
//            Molecule& mol2 = uniqueMolecules[j];
//            bool identical = (mol.atoms.size() == mol2.atoms.size() && mol.constraints.size() == mol2.constraints.size());
//
//            // See if the atoms are identical.
//
//            int atomOffset = mol2.atoms[0]-mol.atoms[0];
//            for (int i = 0; i < (int) mol.atoms.size() && identical; i++) {
//                if (mol.atoms[i] != mol2.atoms[i]-atomOffset || system.getParticleMass(mol.atoms[i]) != system.getParticleMass(mol2.atoms[i]))
//                    identical = false;
//                for (int k = 0; k < (int) forces.size(); k++)
//                    if (!forces[k]->areParticlesIdentical(mol.atoms[i], mol2.atoms[i]))
//                        identical = false;
//            }
//            
//            // See if the constraints are identical.
//
//            for (int i = 0; i < (int) mol.constraints.size() && identical; i++) {
//                int c1particle1, c1particle2, c2particle1, c2particle2;
//                double distance1, distance2;
//                system.getConstraintParameters(mol.constraints[i], c1particle1, c1particle2, distance1);
//                system.getConstraintParameters(mol2.constraints[i], c2particle1, c2particle2, distance2);
//                if (c1particle1 != c2particle1-atomOffset || c1particle2 != c2particle2-atomOffset || distance1 != distance2)
//                    identical = false;
//            }
//
//            // See if the force groups are identical.
//
//            for (int i = 0; i < (int) forces.size() && identical; i++) {
//                if (mol.groups[i].size() != mol2.groups[i].size())
//                    identical = false;
//                for (int k = 0; k < (int) mol.groups[i].size() && identical; k++)
//                    if (!forces[i]->areGroupsIdentical(mol.groups[i][k], mol2.groups[i][k]))
//                        identical = false;
//            }
//            if (identical) {
//                moleculeInstances[j].push_back(molIndex);
//                moleculeOffsets[j].push_back(mol.atoms[0]);
//                isNew = false;
//            }
//        }
//        if (isNew) {
//            uniqueMolecules.push_back(mol);
//            moleculeInstances.push_back(vector<int>());
//            moleculeInstances[moleculeInstances.size()-1].push_back(molIndex);
//            moleculeOffsets.push_back(vector<int>());
//            moleculeOffsets[moleculeOffsets.size()-1].push_back(mol.atoms[0]);
//        }
//    }
//    moleculeGroups.resize(moleculeInstances.size());
//    for (int i = 0; i < (int) moleculeInstances.size(); i++)
//    {
//        moleculeGroups[i].instances = moleculeInstances[i];
//        moleculeGroups[i].offsets = moleculeOffsets[i];
//        vector<int>& atoms = uniqueMolecules[i].atoms;
//        moleculeGroups[i].atoms.resize(atoms.size());
//        for (int j = 0; j < (int) atoms.size(); j++)
//            moleculeGroups[i].atoms[j] = atoms[j]-atoms[0];
//    }
//}
//
//void CudaContext::invalidateMolecules() {
//    moleculesInvalid = true;
//}
//
//
//void OpenCLContext::validateMolecules() {
//    moleculesInvalid = false;
//    if (numAtoms == 0 || nonbonded == NULL || !nonbonded->getUseCutoff())
//        return;
//    bool valid = true;
//    for (int group = 0; valid && group < (int) moleculeGroups.size(); group++) {
//        MoleculeGroup& mol = moleculeGroups[group];
//        vector<int>& instances = mol.instances;
//        vector<int>& offsets = mol.offsets;
//        vector<int>& atoms = mol.atoms;
//        int numMolecules = instances.size();
//        Molecule& m1 = molecules[instances[0]];
//        int offset1 = offsets[0];
//        for (int j = 1; valid && j < numMolecules; j++) {
//            // See if the atoms are identical.
//
//            Molecule& m2 = molecules[instances[j]];
//            int offset2 = offsets[j];
//            for (int i = 0; i < (int) atoms.size() && valid; i++) {
//                for (int k = 0; k < (int) forces.size(); k++)
//                    if (!forces[k]->areParticlesIdentical(atoms[i]+offset1, atoms[i]+offset2))
//                        valid = false;
//            }
//
//            // See if the force groups are identical.
//
//            for (int i = 0; i < (int) forces.size() && valid; i++) {
//                for (int k = 0; k < (int) m1.groups[i].size() && valid; k++)
//                    if (!forces[i]->areGroupsIdentical(m1.groups[i][k], m2.groups[i][k]))
//                        valid = false;
//            }
//        }
//    }
//    if (valid)
//        return;
//    
//    // The list of which molecules are identical is no longer valid.  We need to restore the
//    // atoms to their original order, rebuild the list of identical molecules, and sort them
//    // again.
//    
//    vector<mm_float4> newPosq(numAtoms);
//    vector<mm_float4> newVelm(numAtoms);
//    vector<mm_int4> newCellOffsets(numAtoms);
//    posq->download();
//    velm->download();
//    for (int i = 0; i < numAtoms; i++) {
//        int index = atomIndex->get(i);
//        newPosq[index] = posq->get(i);
//        newVelm[index] = velm->get(i);
//        newCellOffsets[index] = posCellOffsets[i];
//    }
//    for (int i = 0; i < numAtoms; i++) {
//        posq->set(i, newPosq[i]);
//        velm->set(i, newVelm[i]);
//        atomIndex->set(i, i);
//        posCellOffsets[i] = newCellOffsets[i];
//    }
//    posq->upload();
//    velm->upload();
//    atomIndex->upload();
//    findMoleculeGroups();
//    for (int i = 0; i < (int) reorderListeners.size(); i++)
//        reorderListeners[i]->execute();
//}
//
//void OpenCLContext::reorderAtoms(bool enforcePeriodic) {
//    if (numAtoms == 0 || nonbonded == NULL || !nonbonded->getUseCutoff())
//        return;
//    if (moleculesInvalid)
//        validateMolecules();
//    atomsWereReordered = true;
//
//    // Find the range of positions and the number of bins along each axis.
//
//    posq->download();
//    velm->download();
//    float minx = posq->get(0).x, maxx = posq->get(0).x;
//    float miny = posq->get(0).y, maxy = posq->get(0).y;
//    float minz = posq->get(0).z, maxz = posq->get(0).z;
//    if (nonbonded->getUsePeriodic()) {
//        minx = miny = minz = 0.0;
//        maxx = periodicBoxSize.x;
//        maxy = periodicBoxSize.y;
//        maxz = periodicBoxSize.z;
//    }
//    else {
//        for (int i = 1; i < numAtoms; i++) {
//            const mm_float4& pos = posq->get(i);
//            minx = min(minx, pos.x);
//            maxx = max(maxx, pos.x);
//            miny = min(miny, pos.y);
//            maxy = max(maxy, pos.y);
//            minz = min(minz, pos.z);
//            maxz = max(maxz, pos.z);
//        }
//    }
//
//    // Loop over each group of identical molecules and reorder them.
//
//    vector<int> originalIndex(numAtoms);
//    vector<mm_float4> newPosq(numAtoms);
//    vector<mm_float4> newVelm(numAtoms);
//    vector<mm_int4> newCellOffsets(numAtoms);
//    for (int group = 0; group < (int) moleculeGroups.size(); group++) {
//        // Find the center of each molecule.
//
//        MoleculeGroup& mol = moleculeGroups[group];
//        int numMolecules = mol.offsets.size();
//        vector<int>& atoms = mol.atoms;
//        vector<mm_float4> molPos(numMolecules);
//        float invNumAtoms = 1.0f/atoms.size();
//        for (int i = 0; i < numMolecules; i++) {
//            molPos[i].x = 0.0f;
//            molPos[i].y = 0.0f;
//            molPos[i].z = 0.0f;
//            for (int j = 0; j < (int)atoms.size(); j++) {
//                int atom = atoms[j]+mol.offsets[i];
//                const mm_float4& pos = posq->get(atom);
//                molPos[i].x += pos.x;
//                molPos[i].y += pos.y;
//                molPos[i].z += pos.z;
//            }
//            molPos[i].x *= invNumAtoms;
//            molPos[i].y *= invNumAtoms;
//            molPos[i].z *= invNumAtoms;
//        }
//        if (nonbonded->getUsePeriodic()) {
//            // Move each molecule position into the same box.
//
//            for (int i = 0; i < numMolecules; i++) {
//                int xcell = (int) floor(molPos[i].x*invPeriodicBoxSize.x);
//                int ycell = (int) floor(molPos[i].y*invPeriodicBoxSize.y);
//                int zcell = (int) floor(molPos[i].z*invPeriodicBoxSize.z);
//                float dx = xcell*periodicBoxSize.x;
//                float dy = ycell*periodicBoxSize.y;
//                float dz = zcell*periodicBoxSize.z;
//                if (dx != 0.0f || dy != 0.0f || dz != 0.0f) {
//                    molPos[i].x -= dx;
//                    molPos[i].y -= dy;
//                    molPos[i].z -= dz;
//                    if (enforcePeriodic) {
//                        for (int j = 0; j < (int) atoms.size(); j++) {
//                            int atom = atoms[j]+mol.offsets[i];
//                            mm_float4 p = posq->get(atom);
//                            p.x -= dx;
//                            p.y -= dy;
//                            p.z -= dz;
//                            posq->set(atom, p);
//                            posCellOffsets[atom].x -= xcell;
//                            posCellOffsets[atom].y -= ycell;
//                            posCellOffsets[atom].z -= zcell;
//                        }
//                    }
//                }
//            }
//        }
//
//        // Select a bin for each molecule, then sort them by bin.
//
//        bool useHilbert = (numMolecules > 5000 || atoms.size() > 8); // For small systems, a simple zigzag curve works better than a Hilbert curve.
//        float binWidth;
//        if (useHilbert)
//            binWidth = (float)(max(max(maxx-minx, maxy-miny), maxz-minz)/255.0);
//        else
//            binWidth = (float)(0.2*nonbonded->getCutoffDistance());
//        float invBinWidth = 1.0f/binWidth;
//        int xbins = 1 + (int) ((maxx-minx)*invBinWidth);
//        int ybins = 1 + (int) ((maxy-miny)*invBinWidth);
//        vector<pair<int, int> > molBins(numMolecules);
//        bitmask_t coords[3];
//        for (int i = 0; i < numMolecules; i++) {
//            int x = (int) ((molPos[i].x-minx)*invBinWidth);
//            int y = (int) ((molPos[i].y-miny)*invBinWidth);
//            int z = (int) ((molPos[i].z-minz)*invBinWidth);
//            int bin;
//            if (useHilbert) {
//                coords[0] = x;
//                coords[1] = y;
//                coords[2] = z;
//                bin = (int) hilbert_c2i(3, 8, coords);
//            }
//            else {
//                int yodd = y&1;
//                int zodd = z&1;
//                bin = z*xbins*ybins;
//                bin += (zodd ? ybins-y : y)*xbins;
//                bin += (yodd ? xbins-x : x);
//            }
//            molBins[i] = pair<int, int>(bin, i);
//        }
//        sort(molBins.begin(), molBins.end());
//
//        // Reorder the atoms.
//
//        for (int i = 0; i < numMolecules; i++) {
//            for (int j = 0; j < (int)atoms.size(); j++) {
//                int oldIndex = mol.offsets[molBins[i].second]+atoms[j];
//                int newIndex = mol.offsets[i]+atoms[j];
//                originalIndex[newIndex] = atomIndex->get(oldIndex);
//                newPosq[newIndex] = posq->get(oldIndex);
//                newVelm[newIndex] = velm->get(oldIndex);
//                newCellOffsets[newIndex] = posCellOffsets[oldIndex];
//            }
//        }
//    }
//
//    // Update the streams.
//
//    for (int i = 0; i < numAtoms; i++) {
//        posq->set(i, newPosq[i]);
//        velm->set(i, newVelm[i]);
//        atomIndex->set(i, originalIndex[i]);
//        posCellOffsets[i] = newCellOffsets[i];
//    }
//    posq->upload();
//    velm->upload();
//    atomIndex->upload();
//    for (int i = 0; i < (int) reorderListeners.size(); i++)
//        reorderListeners[i]->execute();
//}

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
