/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2009-2026 Stanford University and the Authors.      *
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

#include <cmath>
#include "MetalContext.h"
#include "MetalEvent.h"
#include "MetalFFT3D.h"
#include "MetalQueue.h"
#include "MetalKernels.h"
#include "MetalKernelSources.h"
#include "MetalProgram.h"
#include "MetalSort.h"
#include "openmm/common/ComputeArray.h"
#include "openmm/common/ContextSelector.h"
#include "SHA1.h"
#include "openmm/MonteCarloFlexibleBarostat.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/VirtualSite.h"
#include "openmm/internal/ContextImpl.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <set>
#include <sstream>
#include <typeinfo>
#include <sys/stat.h>
#include <unistd.h>
#include <IOKit/IOKitLib.h>

using namespace OpenMM;
using namespace std;

const int MetalContext::ThreadBlockSize = 64;
const int MetalContext::TileSize = sizeof(tileflags)*8;

// Uncomment the following line to enable printf() calls inside kernels.  This affects performance, so it should only
// be done for debugging.
//#define ENABLE_PRINTF

MetalContext::MetalContext(const System& system, const string& precision, MetalPlatform::PlatformData& platformData,
        MetalContext* originalContext) : ComputeContext(system), platformData(platformData), integration(NULL),
        expression(NULL), bonded(NULL), nonbonded(NULL) {
#ifdef ENABLE_PRINTF
    setenv("MTL_LOG_LEVEL", "MTLLogLevelDebug", 0);
    setenv("MTL_LOG_TO_STDERR", "1", 0);
    setenv("MTL_LOG_BUFFER_SIZE", "100000", 0);
    compilationDefines["printf"] = "os_log_default.log";
#endif
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
    contextIndex = platformData.contexts.size();
    if (originalContext == NULL) {
        device = MTL::CreateSystemDefaultDevice();
        defaultQueue = shared_ptr<ComputeQueueImpl>(new MetalQueue(*device));
        isLinkedContext = false;
    }
    else {
        device = originalContext->getDevice().retain();
        defaultQueue = originalContext->defaultQueue;
        isLinkedContext = true;
    }

    currentQueue = defaultQueue;
    numAtoms = system.getNumParticles();
    paddedNumAtoms = TileSize*((numAtoms+TileSize-1)/TileSize);
    numAtomBlocks = (paddedNumAtoms+(TileSize-1))/TileSize;

    // Determine the number of cores in the GPU.  Why does Apple make this so hard?

    numGpuCores = -1;
    io_iterator_t iterator;
    if (IOServiceGetMatchingServices(kIOMainPortDefault, IOServiceMatching("AGXAccelerator"), &iterator) == KERN_SUCCESS) {
        io_object_t entry;
        while ((entry = IOIteratorNext(iterator)) != 0) {
            CFTypeRef value = IORegistryEntryCreateCFProperty(entry, CFSTR("gpu-core-count"), kCFAllocatorDefault, 0);
            if (value != NULL) {
                if (CFGetTypeID(value) == CFNumberGetTypeID())
                    CFNumberGetValue((CFNumberRef) value, kCFNumberIntType, &numGpuCores);
                CFRelease(value);
            }
            IOObjectRelease(entry);
        }
        IOObjectRelease(iterator);
    }
    if (numGpuCores == -1)
        throw OpenMMException("Unable to determine number of GPU cores");
    numThreadBlocks = 12*numGpuCores;

    // Decide whether the fast versions of math routines are sufficiently accurate to use.

    ComputeProgram program = compileProgram(MetalKernelSources::utilities);
    ComputeKernel accuracyKernel = program->createKernel("determineNativeAccuracy");
    int numValues = 20;
    MetalArray valuesArray(*this, 6*numValues, sizeof(float), "values");
    vector<float> values(valuesArray.getSize(), 0.0);
    float nextValue = 1e-4f;
    for (int i = 0; i < numValues; i++) {
        values[6*i] = nextValue;
        nextValue *= (float) M_PI;
    }
    valuesArray.upload(values);
    accuracyKernel->addArg(valuesArray);
    accuracyKernel->addArg(numValues);
    accuracyKernel->execute(numValues);
    valuesArray.download(values);
    double maxSqrtError = 0.0, maxRsqrtError = 0.0, maxRecipError = 0.0, maxExpError = 0.0, maxLogError = 0.0;
    for (int i = 0; i < numValues; i++) {
        double v = values[6*i];
        double correctSqrt = sqrt(v);
        maxSqrtError = max(maxSqrtError, fabs(correctSqrt-values[6*i+1])/correctSqrt);
        maxRsqrtError = max(maxRsqrtError, fabs(1.0/correctSqrt-values[6*i+2])*correctSqrt);
        maxRecipError = max(maxRecipError, fabs(1.0/v-values[6*i+3])/values[6*i+3]);
        maxExpError = max(maxExpError, fabs(exp(v)-values[6*i+4])/values[6*i+4]);
        maxLogError = max(maxLogError, fabs(log(v)-values[6*i+5])/values[6*i+5]);
    }
    compilationDefines["SQRT"] = (maxSqrtError < 1e-6) ? "fast::sqrt" : "sqrt";
    compilationDefines["RSQRT"] = (maxRsqrtError < 1e-6) ? "fast::rsqrt" : "rsqrt";
    compilationDefines["RECIP(v)"] = (maxRecipError < 1e-6) ? "fast::divide(1.0, v)" : "(1.0/(v))";
    compilationDefines["EXP"] = (maxExpError < 1e-6) ? "fast::exp" : "exp";
    compilationDefines["LOG"] = (maxLogError < 1e-6) ? "fast::log" : "log";

    // Set defines based on the requested precision.

    compilationDefines["POW"] = "pow";
    compilationDefines["COS"] = "cos";
    compilationDefines["SIN"] = "sin";
    compilationDefines["TAN"] = "tan";
    compilationDefines["ACOS"] = "acos";
    compilationDefines["ASIN"] = "asin";
    compilationDefines["ATAN"] = "atan";
    compilationDefines["ERF"] = "erf";
    compilationDefines["ERFC"] = "erfc";
    compilationDefines["FMA"] = "fma";
    compilationDefines["FABS"] = "fabs";
    compilationDefines["make_real2"] = "make_float2";
    compilationDefines["make_real3"] = "make_float3";
    compilationDefines["make_real4"] = "make_float4";
    compilationDefines["make_mixed2"] = "make_float2";
    compilationDefines["make_mixed3"] = "make_float3";
    compilationDefines["make_mixed4"] = "make_float4";

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
    initializeKernels();

    // Create utilities objects.

    bonded = new BondedUtilities(*this);
    nonbonded = new MetalNonbondedUtilities(*this);
    integration = new MetalIntegrationUtilities(*this, system);
    expression = new ExpressionUtilities(*this);
    clearBuffer(posq);
}

MetalContext::~MetalContext() {
    for (auto force : forces)
        delete force;
    for (auto listener : reorderListeners)
        delete listener;
    for (auto computation : preComputations)
        delete computation;
    for (auto computation : postComputations)
        delete computation;
    if (integration != NULL)
        delete integration;
    if (expression != NULL)
        delete expression;
    if (bonded != NULL)
        delete bonded;
    if (nonbonded != NULL)
        delete nonbonded;
    device->release();
}

void MetalContext::initialize() {
    int numEnergyBuffers = max(numThreadBlocks*ThreadBlockSize, nonbonded->getNumEnergyBuffers());
    energyBuffer.initialize<float>(*this, numEnergyBuffers, "energyBuffer");
    energySum.initialize<float>(*this, numGpuCores, "energySum");
    int pinnedBufferSize = max(paddedNumAtoms*6, numEnergyBuffers);
    pinnedBuffer.resize(pinnedBufferSize*sizeof(float));
    for (int i = 0; i < numAtoms; i++) {
        double mass = system.getParticleMass(i);
        if (useMixedPrecision)
            ((mm_double4*) pinnedBuffer.data())[i] = mm_double4(0.0, 0.0, 0.0, mass == 0.0 ? 0.0 : 1.0/mass);
        else
            ((mm_float4*) pinnedBuffer.data())[i] = mm_float4(0.0f, 0.0f, 0.0f, mass == 0.0 ? 0.0f : (float) (1.0/mass));
    }
    velm.upload(pinnedBuffer.data());
    bonded->initialize(system);
    addAutoclearBuffer(longForceBuffer);
    addAutoclearBuffer(energyBuffer);
    int numEnergyParamDerivs = energyParamDerivNames.size();
    if (numEnergyParamDerivs > 0) {
        if (useMixedPrecision)
            energyParamDerivBuffer.initialize<double>(*this, numEnergyParamDerivs*numEnergyBuffers, "energyParamDerivBuffer");
        else
            energyParamDerivBuffer.initialize<float>(*this, numEnergyParamDerivs*numEnergyBuffers, "energyParamDerivBuffer");
        addAutoclearBuffer(energyParamDerivBuffer);
    }
    findMoleculeGroups();
    nonbonded->initialize(system);
}

void MetalContext::initializeContexts() {
    getPlatformData().initializeContexts(system);
}

FFT3D MetalContext::createFFT(int xsize, int ysize, int zsize, bool realToComplex) {
    return FFT3D(new MetalFFT3D(*this, xsize, ysize, zsize, realToComplex));
}

vector<ComputeContext*> MetalContext::getAllContexts() {
    vector<ComputeContext*> result;
    for (MetalContext* c : platformData.contexts)
        result.push_back(c);
    return result;
}

double& MetalContext::getEnergyWorkspace() {
    return platformData.contextEnergy[contextIndex];
}

ComputeQueue MetalContext::createQueue() {
    return shared_ptr<ComputeQueueImpl>(new MetalQueue(*device));
}

MetalArray* MetalContext::createArray() {
    return new MetalArray();
}

ComputeEvent MetalContext::createEvent() {
    return shared_ptr<ComputeEventImpl>(new MetalEvent(*this));
}

ComputeSort MetalContext::createSort(ComputeSortImpl::SortTrait* trait, unsigned int length, bool uniform) {
    return shared_ptr<ComputeSortImpl>(new MetalSort(*this, trait, length, uniform));
}

static string rewriteKernelArgs(const string& source) {
    static regex kernelMatcher("KERNEL void[\\s\\S]*?\\(([\\s\\S]*?)\\)");
    static regex argMatcher("\\#.*|[^,\\s#][^\\,#]*[^,\\s#]*");
    static regex wordMatcher("[^\\s]+");
    stringstream result;
    int pos = 0;

    // Loop over all kernels declared in the file.

    sregex_iterator nextKernel(source.begin(), source.end(), kernelMatcher);
    sregex_iterator end;
    while (nextKernel != end) {
        smatch kernelMatch = *nextKernel;
        result << source.substr(pos, kernelMatch.position(1)-pos);

        // Loop over arguments to the kernel.

        string argList = kernelMatch[1].str();
        pos = kernelMatch.position(1)+argList.size();
        sregex_iterator nextArg(argList.begin(), argList.end(), argMatcher);
        bool addComma = false;
        while (nextArg != end) {
            smatch argMatch = *nextArg;
            string arg = argMatch.str();
            if (arg.rfind("GLOBAL", 0) == 0 || arg.rfind("LOCAL", 0) == 0) {
                if (addComma)
                    result << ", ";
                addComma = true;
                result << arg;
            }
            else if (arg[0] == '#') {
                result << "\n" << arg << "\n";
            }
            else {
                // This is a primitive value.  We need to transform it into reference to constant memory.

                if (addComma)
                    result << ", ";
                addComma = true;
                result << "constant";
                sregex_iterator nextWord(arg.begin(), arg.end(), wordMatcher);
                bool addAmpersand = true;
                while (nextWord != end) {
                    smatch wordMatch = *nextWord;
                    string word = wordMatch.str();
                    if (word != "const") {
                        result << " " << word;
                        if (addAmpersand && word != "unsigned") {
                            result << "&";
                            addAmpersand = false;
                        }
                    }
                    ++nextWord;
                }
            }
            ++nextArg;
        }
        ++nextKernel;
    }
    result << source.substr(pos, source.size()-pos);
    return result.str();
}

ComputeProgram MetalContext::compileProgram(const string source, const map<string, string>& defines) {
    stringstream src;
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
    src << "typedef float real;\n";
    src << "typedef float2 real2;\n";
    src << "typedef float3 real3;\n";
    src << "typedef float4 real4;\n";
    src << "typedef float mixed;\n";
    src << "typedef float2 mixed2;\n";
    src << "typedef float3 mixed3;\n";
    src << "typedef float4 mixed4;\n";
    src << "typedef unsigned int tileflags;\n";
    src << MetalKernelSources::common << endl;
    for (auto& pair : defines) {
        src << "#define " << pair.first;
        if (!pair.second.empty())
            src << " " << pair.second;
        src << endl;
    }
    if (!defines.empty())
        src << endl;
    src << rewriteKernelArgs(source) << endl;

    // Compile the program.

    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();
    MTL::CompileOptions* options = MTL::CompileOptions::alloc()->init();
    options->setLanguageVersion(MTL::LanguageVersion3_2);
    options->setMathMode(MTL::MathModeSafe);
    options->setMathFloatingPointFunctions(MTL::MathFloatingPointFunctionsPrecise);
#ifdef ENABLE_PRINTF
    options->setEnableLogging(true);
#endif
    NS::Error* error = nullptr;
    MTL::Library* library = device->newLibrary(NS::String::string(src.str().c_str(), NS::UTF8StringEncoding), options, &error);
    options->release();
    string errorString = (error == nullptr ? "" : error->localizedDescription()->utf8String());
    pool->release();
    if (library == nullptr)
        throw OpenMMException("Error compiling program: "+errorString);
    return shared_ptr<ComputeProgramImpl>(new MetalProgram(*this, library));
}

MetalArray& MetalContext::unwrap(ArrayInterface& array) const {
    MetalArray* metalArray;
    ComputeArray* wrapper = dynamic_cast<ComputeArray*>(&array);
    if (wrapper != NULL)
        metalArray = dynamic_cast<MetalArray*>(&wrapper->getArray());
    else
        metalArray = dynamic_cast<MetalArray*>(&array);
    if (metalArray == NULL)
        throw OpenMMException("Array argument is not an MetalArray");
    return *metalArray;
}

int MetalContext::computeThreadBlockSize(double memory) const {
    int maxShared = 32768;
    int max = (int) (maxShared/memory);
    if (max < 64)
        return 32;
    int threads = 64;
    while (threads+64 < max)
        threads += 64;
    return threads;
}

double MetalContext::reduceEnergy() {
    int workGroupSize  = 512;
    reduceEnergyKernel->setArg(0, energyBuffer);
    reduceEnergyKernel->setArg(1, energySum);
    reduceEnergyKernel->setArg(2, energyBuffer.getSize());
    reduceEnergyKernel->setArg(3, workGroupSize);
    reduceEnergyKernel->execute(workGroupSize*energySum.getSize(), workGroupSize);
    energySum.download(pinnedBuffer.data());
    double result = 0;
    if (getUseMixedPrecision()) {
        for (int i = 0; i < energySum.getSize(); i++)
            result += ((double*) pinnedBuffer.data())[i];
    }
    else {
        for (int i = 0; i < energySum.getSize(); i++)
            result += ((float*) pinnedBuffer.data())[i];
    }
    return result;
}

void MetalContext::flushQueue() {
    MetalQueue* queue = dynamic_cast<MetalQueue*>(getCurrentQueue().get());
    queue->flush();
}
