/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2018 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "CpuPlatform.h"
#include "CpuKernelFactory.h"
#include "CpuKernels.h"
#include "CpuSETTLE.h"
#include "ReferenceConstraints.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/hardware.h"
#include "openmm/internal/vectorize.h"
#include <algorithm>
#include <sstream>
#include <stdlib.h>

using namespace OpenMM;
using namespace std;

#ifdef OPENMM_CPU_BUILDING_STATIC_LIBRARY
extern "C" void registerCpuPlatform() {
    if (CpuPlatform::isProcessorSupported())
        Platform::registerPlatform(new CpuPlatform());
}
#else
extern "C" OPENMM_EXPORT_CPU void registerPlatforms() {
    // Only register this platform if the CPU supports SSE 4.1.

    if (CpuPlatform::isProcessorSupported())
        Platform::registerPlatform(new CpuPlatform());
}
#endif

map<const ContextImpl*, CpuPlatform::PlatformData*> CpuPlatform::contextData;

CpuPlatform::CpuPlatform() {
    deprecatedPropertyReplacements["CpuThreads"] = CpuThreads();
    CpuKernelFactory* factory = new CpuKernelFactory();
    registerKernelFactory(CalcForcesAndEnergyKernel::Name(), factory);
    registerKernelFactory(CalcHarmonicAngleForceKernel::Name(), factory);
    registerKernelFactory(CalcPeriodicTorsionForceKernel::Name(), factory);
    registerKernelFactory(CalcRBTorsionForceKernel::Name(), factory);
    registerKernelFactory(CalcNonbondedForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomNonbondedForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomManyParticleForceKernel::Name(), factory);
    registerKernelFactory(CalcGBSAOBCForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomGBForceKernel::Name(), factory);
    registerKernelFactory(CalcGayBerneForceKernel::Name(), factory);
    registerKernelFactory(IntegrateLangevinStepKernel::Name(), factory);
    platformProperties.push_back(CpuThreads());
    platformProperties.push_back(CpuDeterministicForces());
    int threads = getNumProcessors();
    char* threadsEnv = getenv("OPENMM_CPU_THREADS");
    if (threadsEnv != NULL)
        stringstream(threadsEnv) >> threads;
    stringstream defaultThreads;
    defaultThreads << threads;
    setPropertyDefaultValue(CpuThreads(), defaultThreads.str());
    setPropertyDefaultValue(CpuDeterministicForces(), "false");
}

const string& CpuPlatform::getPropertyValue(const Context& context, const string& property) const {
    const ContextImpl& impl = getContextImpl(context);
    const PlatformData& data = getPlatformData(impl);
    string propertyName = property;
    if (deprecatedPropertyReplacements.find(property) != deprecatedPropertyReplacements.end())
        propertyName = deprecatedPropertyReplacements.find(property)->second;
    map<string, string>::const_iterator value = data.propertyValues.find(propertyName);
    if (value != data.propertyValues.end())
        return value->second;
    return ReferencePlatform::getPropertyValue(context, property);
}

double CpuPlatform::getSpeed() const {
    return 10;
}

bool CpuPlatform::supportsDoublePrecision() const {
    return false;
}

bool CpuPlatform::isProcessorSupported() {
    return isVec4Supported();
}

void CpuPlatform::contextCreated(ContextImpl& context, const map<string, string>& properties) const {
    ReferencePlatform::contextCreated(context, properties);
    const string& threadsPropValue = (properties.find(CpuThreads()) == properties.end() ?
            getPropertyDefaultValue(CpuThreads()) : properties.find(CpuThreads())->second);
    string deterministicForcesValue = (properties.find(CpuDeterministicForces()) == properties.end() ?
            getPropertyDefaultValue(CpuDeterministicForces()) : properties.find(CpuDeterministicForces())->second);
    int numThreads;
    stringstream(threadsPropValue) >> numThreads;
    transform(deterministicForcesValue.begin(), deterministicForcesValue.end(), deterministicForcesValue.begin(), ::tolower);
    bool deterministicForces = (deterministicForcesValue == "true");
    PlatformData* data = new PlatformData(context.getSystem().getNumParticles(), numThreads, deterministicForces);
    contextData[&context] = data;
    ReferenceConstraints& constraints = *(ReferenceConstraints*) reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData())->constraints;
    if (constraints.settle != NULL) {
        CpuSETTLE* parallelSettle = new CpuSETTLE(context.getSystem(), *(ReferenceSETTLEAlgorithm*) constraints.settle, data->threads);
        delete constraints.settle;
        constraints.settle = parallelSettle;
    }
}

void CpuPlatform::contextDestroyed(ContextImpl& context) const {
    PlatformData* data = contextData[&context];
    delete data;
    contextData.erase(&context);
    ReferencePlatform::PlatformData* refPlatformData = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    delete refPlatformData;
}

CpuPlatform::PlatformData& CpuPlatform::getPlatformData(ContextImpl& context) {
    return *contextData[&context];
}

const CpuPlatform::PlatformData& CpuPlatform::getPlatformData(const ContextImpl& context) {
    return *contextData[&context];
}

CpuPlatform::PlatformData::PlatformData(int numParticles, int numThreads, bool deterministicForces) : posq(4*numParticles), threads(numThreads),
        deterministicForces(deterministicForces), neighborList(NULL), cutoff(0.0), paddedCutoff(0.0), anyExclusions(false), currentPosqIndex(-1), nextPosqIndex(0) {
    numThreads = threads.getNumThreads();
    threadForce.resize(numThreads);
    for (int i = 0; i < numThreads; i++)
        threadForce[i].resize(4*numParticles);
    isPeriodic = false;
    stringstream threadsProperty;
    threadsProperty << numThreads;
    propertyValues[CpuThreads()] = threadsProperty.str();
    propertyValues[CpuDeterministicForces()] = deterministicForces ? "true" : "false";
}

CpuPlatform::PlatformData::~PlatformData() {
    if (neighborList != NULL)
        delete neighborList;
}

bool isVec8Supported();

void CpuPlatform::PlatformData::requestNeighborList(double cutoffDistance, double padding, bool useExclusions, const vector<set<int> >& exclusionList) {
    if (neighborList == NULL)
        neighborList = new CpuNeighborList(isVec8Supported() ? 8 : 4);
    if (cutoffDistance > cutoff)
        cutoff = cutoffDistance;
    if (cutoffDistance+padding > paddedCutoff)
        paddedCutoff = cutoffDistance+padding;
    if (useExclusions) {
        if (anyExclusions && exclusions != exclusionList)
            throw OpenMMException("All Forces must have identical exclusions");
        else {
            exclusions = exclusionList;
            anyExclusions = true;
        }
    }
    else if (!anyExclusions)
        exclusions = exclusionList;
}

int CpuPlatform::PlatformData::requestPosqIndex() {
    return nextPosqIndex++;
}