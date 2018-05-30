#ifndef OPENMM_CPUPLATFORM_H_
#define OPENMM_CPUPLATFORM_H_

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

#include "AlignedArray.h"
#include "CpuRandom.h"
#include "CpuNeighborList.h"
#include "ReferencePlatform.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/ThreadPool.h"
#include "windowsExportCpu.h"
#include <map>

namespace OpenMM {
    
/**
 * This Platform subclass uses CPU implementations of the OpenMM kernels.
 */

class OPENMM_EXPORT_CPU CpuPlatform : public ReferencePlatform {
public:
    class PlatformData;
    CpuPlatform();
    const std::string& getName() const {
        static const std::string name = "CPU";
        return name;
    }
    double getSpeed() const;
    const std::string& getPropertyValue(const Context& context, const std::string& property) const;
    bool supportsDoublePrecision() const;
    static bool isProcessorSupported();
    void contextCreated(ContextImpl& context, const std::map<std::string, std::string>& properties) const;
    void contextDestroyed(ContextImpl& context) const;
    /**
     * This is the name of the parameter for selecting the number of threads to use.
     */
    static const std::string& CpuThreads() {
        static const std::string key = "Threads";
        return key;
    }
    /**
     * This is the name of the parameter for requesting that force computations be deterministic.  Setting
     * this to "true" DOES NOT GUARANTEE that the forces will actually be fully deterministic, but it does
     * try to reduce the variation in them at the cost of a small loss in performance.
     */
    static const std::string& CpuDeterministicForces() {
        static const std::string key = "DeterministicForces";
        return key;
    }
    /**
     * We cannot use the standard mechanism for platform data, because that is already used by the superclass.
     * Instead, we maintain a table of ContextImpls to PlatformDatas.
     */
    static PlatformData& getPlatformData(ContextImpl& context);
    static const PlatformData& getPlatformData(const ContextImpl& context);
private:
    static std::map<const ContextImpl*, PlatformData*> contextData;
};

class CpuPlatform::PlatformData {
public:
    PlatformData(int numParticles, int numThreads, bool deterministicForces);
    ~PlatformData();
    void requestNeighborList(double cutoffDistance, double padding, bool useExclusions, const std::vector<std::set<int> >& exclusionList);
    int requestPosqIndex();
    AlignedArray<float> posq;
    std::vector<AlignedArray<float> > threadForce;
    ThreadPool threads;
    bool isPeriodic;
    CpuRandom random;
    std::map<std::string, std::string> propertyValues;
    CpuNeighborList* neighborList;
    double cutoff, paddedCutoff;
    bool anyExclusions, deterministicForces;
    int currentPosqIndex, nextPosqIndex;
    std::vector<std::set<int> > exclusions;
};

} // namespace OpenMM

#endif /*OPENMM_CPUPLATFORM_H_*/
