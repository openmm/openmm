/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011-2015 Stanford University and the Authors.      *
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

#include "CudaParallelKernels.h"
#include "CudaKernelSources.h"

using namespace OpenMM;
using namespace std;


#define CHECK_RESULT(result, prefix) \
if (result != CUDA_SUCCESS) { \
    std::stringstream m; \
    m<<prefix<<": "<<cu.getErrorString(result)<<" ("<<result<<")"<<" at "<<__FILE__<<":"<<__LINE__; \
    throw OpenMMException(m.str());\
}

/**
 * Get the current clock time, measured in microseconds.
 */
#ifdef _MSC_VER
    #include <Windows.h>
    static long long getTime() {
        FILETIME ft;
        GetSystemTimeAsFileTime(&ft); // 100-nanoseconds since 1-1-1601
        ULARGE_INTEGER result;
        result.LowPart = ft.dwLowDateTime;
        result.HighPart = ft.dwHighDateTime;
        return result.QuadPart/10;
    }
#else
    #include <sys/time.h> 
    static long long getTime() {
        struct timeval tod;
        gettimeofday(&tod, 0);
        return 1000000*tod.tv_sec+tod.tv_usec;
    }
#endif

class CudaParallelCalcForcesAndEnergyKernel::BeginComputationTask : public CudaContext::WorkTask {
public:
    BeginComputationTask(ContextImpl& context, CudaContext& cu, CudaCalcForcesAndEnergyKernel& kernel,
            bool includeForce, bool includeEnergy, int groups, void* pinnedMemory, CUevent event, int2& interactionCount) : context(context), cu(cu), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), groups(groups), pinnedMemory(pinnedMemory), event(event), interactionCount(interactionCount) {
    }
    void execute() {
        // Copy coordinates over to this device and execute the kernel.

        cu.setAsCurrent();
        if (cu.getContextIndex() > 0) {
            cuStreamWaitEvent(cu.getCurrentStream(), event, 0);
            if (!cu.getPlatformData().peerAccessSupported)
                cu.getPosq().upload(pinnedMemory, false);
        }
        kernel.beginComputation(context, includeForce, includeEnergy, groups);
        if (cu.getNonbondedUtilities().getUsePeriodic())
            cu.getNonbondedUtilities().getInteractionCount().download(&interactionCount, false);
    }
private:
    ContextImpl& context;
    CudaContext& cu;
    CudaCalcForcesAndEnergyKernel& kernel;
    bool includeForce, includeEnergy;
    int groups;
    void* pinnedMemory;
    CUevent event;
    int2& interactionCount;
};

class CudaParallelCalcForcesAndEnergyKernel::FinishComputationTask : public CudaContext::WorkTask {
public:
    FinishComputationTask(ContextImpl& context, CudaContext& cu, CudaCalcForcesAndEnergyKernel& kernel,
            bool includeForce, bool includeEnergy, int groups, double& energy, long long& completionTime, long long* pinnedMemory, CudaArray& contextForces, bool& valid, int2& interactionCount) :
            context(context), cu(cu), kernel(kernel), includeForce(includeForce), includeEnergy(includeEnergy), groups(groups), energy(energy),
            completionTime(completionTime), pinnedMemory(pinnedMemory), contextForces(contextForces), valid(valid), interactionCount(interactionCount) {
    }
    void execute() {
        // Execute the kernel, then download forces.
        
        energy += kernel.finishComputation(context, includeForce, includeEnergy, groups, valid);
        if (cu.getComputeForceCount() < 200) {
            // Record timing information for load balancing.  Since this takes time, only do it at the start of the simulation.

            CHECK_RESULT(cuCtxSynchronize(), "Error synchronizing CUDA context");
            completionTime = getTime();
        }
        if (includeForce) {
            if (cu.getContextIndex() > 0) {
                int numAtoms = cu.getPaddedNumAtoms();
                if (cu.getPlatformData().peerAccessSupported) {
                    int numBytes = numAtoms*3*sizeof(long long);
                    int offset = (cu.getContextIndex()-1)*numBytes;
                    CudaContext& context0 = *cu.getPlatformData().contexts[0];
                    CHECK_RESULT(cuMemcpy(contextForces.getDevicePointer()+offset, cu.getForce().getDevicePointer(), numBytes), "Error copying forces");
                }
                else
                    cu.getForce().download(&pinnedMemory[(cu.getContextIndex()-1)*numAtoms*3]);
            }
        }
        if (cu.getNonbondedUtilities().getUsePeriodic() && (interactionCount.x > cu.getNonbondedUtilities().getInteractingTiles().getSize() ||
                interactionCount.y > cu.getNonbondedUtilities().getSinglePairs().getSize())) {
            valid = false;
            cu.getNonbondedUtilities().updateNeighborListSize();
        }
    }
private:
    ContextImpl& context;
    CudaContext& cu;
    CudaCalcForcesAndEnergyKernel& kernel;
    bool includeForce, includeEnergy;
    int groups;
    double& energy;
    long long& completionTime;
    long long* pinnedMemory;
    CudaArray& contextForces;
    bool& valid;
    int2& interactionCount;
};

CudaParallelCalcForcesAndEnergyKernel::CudaParallelCalcForcesAndEnergyKernel(string name, const Platform& platform, CudaPlatform::PlatformData& data) :
        CalcForcesAndEnergyKernel(name, platform), data(data), completionTimes(data.contexts.size()), contextNonbondedFractions(data.contexts.size()),
        interactionCounts(NULL), contextForces(NULL), pinnedPositionBuffer(NULL), pinnedForceBuffer(NULL) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CudaCalcForcesAndEnergyKernel(name, platform, *data.contexts[i])));
}

CudaParallelCalcForcesAndEnergyKernel::~CudaParallelCalcForcesAndEnergyKernel() {
    data.contexts[0]->setAsCurrent();
    if (contextForces != NULL)
        delete contextForces;
    if (pinnedPositionBuffer != NULL)
        cuMemFreeHost(pinnedPositionBuffer);
    if (pinnedForceBuffer != NULL)
        cuMemFreeHost(pinnedForceBuffer);
    cuEventDestroy(event);
    cuStreamDestroy(peerCopyStream);
    if (interactionCounts != NULL)
        cuMemFreeHost(interactionCounts);
}

void CudaParallelCalcForcesAndEnergyKernel::initialize(const System& system) {
    CudaContext& cu = *data.contexts[0];
    cu.setAsCurrent();
    CUmodule module = cu.createModule(CudaKernelSources::parallel);
    sumKernel = cu.getKernel(module, "sumForces");
    int numContexts = data.contexts.size();
    for (int i = 0; i < numContexts; i++)
        getKernel(i).initialize(system);
    for (int i = 0; i < numContexts; i++)
        contextNonbondedFractions[i] = 1/(double) numContexts;
    CHECK_RESULT(cuEventCreate(&event, 0), "Error creating event");
    CHECK_RESULT(cuStreamCreate(&peerCopyStream, CU_STREAM_NON_BLOCKING), "Error creating stream");
    CHECK_RESULT(cuMemHostAlloc((void**) &interactionCounts, numContexts*sizeof(int2), 0), "Error creating interaction counts buffer");
}

void CudaParallelCalcForcesAndEnergyKernel::beginComputation(ContextImpl& context, bool includeForce, bool includeEnergy, int groups) {
    CudaContext& cu = *data.contexts[0];
    cu.setAsCurrent();
    if (contextForces == NULL) {
        contextForces = CudaArray::create<long long>(cu, 3*(data.contexts.size()-1)*cu.getPaddedNumAtoms(), "contextForces");
        CHECK_RESULT(cuMemHostAlloc((void**) &pinnedForceBuffer, 3*(data.contexts.size()-1)*cu.getPaddedNumAtoms()*sizeof(long long), CU_MEMHOSTALLOC_PORTABLE), "Error allocating pinned memory");
        CHECK_RESULT(cuMemHostAlloc(&pinnedPositionBuffer, cu.getPaddedNumAtoms()*(cu.getUseDoublePrecision() ? sizeof(double4) : sizeof(float4)), CU_MEMHOSTALLOC_PORTABLE), "Error allocating pinned memory");
    }

    // Copy coordinates over to each device and execute the kernel.
    
    if (!cu.getPlatformData().peerAccessSupported) {
        cu.getPosq().download(pinnedPositionBuffer, false);
        cuEventRecord(event, cu.getCurrentStream());
    }
    else {
        int numBytes = cu.getPosq().getSize()*cu.getPosq().getElementSize();
        cuEventRecord(event, cu.getCurrentStream());
        cuStreamWaitEvent(peerCopyStream, event, 0);
        for (int i = 1; i < (int) data.contexts.size(); i++)
            CHECK_RESULT(cuMemcpyAsync(data.contexts[i]->getPosq().getDevicePointer(), cu.getPosq().getDevicePointer(), numBytes, peerCopyStream), "Error copying positions");
        cuEventRecord(event, peerCopyStream);
    }
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        data.contextEnergy[i] = 0.0;
        CudaContext& cu = *data.contexts[i];
        CudaContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new BeginComputationTask(context, cu, getKernel(i), includeForce, includeEnergy, groups, pinnedPositionBuffer, event, interactionCounts[i]));
    }
}

double CudaParallelCalcForcesAndEnergyKernel::finishComputation(ContextImpl& context, bool includeForce, bool includeEnergy, int groups, bool& valid) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        CudaContext& cu = *data.contexts[i];
        CudaContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new FinishComputationTask(context, cu, getKernel(i), includeForce, includeEnergy, groups, data.contextEnergy[i], completionTimes[i], pinnedForceBuffer, *contextForces, valid, interactionCounts[i]));
    }
    data.syncContexts();
    double energy = 0.0;
    for (int i = 0; i < (int) data.contextEnergy.size(); i++)
        energy += data.contextEnergy[i];
    if (includeForce && valid) {
        // Sum the forces from all devices.
        
        CudaContext& cu = *data.contexts[0];
        if (!cu.getPlatformData().peerAccessSupported)
            contextForces->upload(pinnedForceBuffer, false);
        int bufferSize = 3*cu.getPaddedNumAtoms();
        int numBuffers = data.contexts.size()-1;
        void* args[] = {&cu.getForce().getDevicePointer(), &contextForces->getDevicePointer(), &bufferSize, &numBuffers};
        cu.executeKernel(sumKernel, args, bufferSize);
        
        // Balance work between the contexts by transferring a little nonbonded work from the context that
        // finished last to the one that finished first.
        
        if (cu.getComputeForceCount() < 200) {
            int firstIndex = 0, lastIndex = 0;
            for (int i = 0; i < (int) completionTimes.size(); i++) {
                if (completionTimes[i] < completionTimes[firstIndex])
                    firstIndex = i;
                if (completionTimes[i] > completionTimes[lastIndex])
                    lastIndex = i;
            }
            double fractionToTransfer = min(0.01, contextNonbondedFractions[lastIndex]);
            contextNonbondedFractions[firstIndex] += fractionToTransfer;
            contextNonbondedFractions[lastIndex] -= fractionToTransfer;
            double startFraction = 0.0;
            for (int i = 0; i < (int) contextNonbondedFractions.size(); i++) {
                double endFraction = startFraction+contextNonbondedFractions[i];
                if (i == contextNonbondedFractions.size()-1)
                    endFraction = 1.0; // Avoid roundoff error
                data.contexts[i]->getNonbondedUtilities().setAtomBlockRange(startFraction, endFraction);
                startFraction = endFraction;
            }
	}
    }
    return energy;
}

class CudaParallelCalcHarmonicBondForceKernel::Task : public CudaContext::WorkTask {
public:
    Task(ContextImpl& context, CudaCalcHarmonicBondForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CudaCalcHarmonicBondForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

CudaParallelCalcHarmonicBondForceKernel::CudaParallelCalcHarmonicBondForceKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data, const System& system) :
        CalcHarmonicBondForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CudaCalcHarmonicBondForceKernel(name, platform, *data.contexts[i], system)));
}

void CudaParallelCalcHarmonicBondForceKernel::initialize(const System& system, const HarmonicBondForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CudaParallelCalcHarmonicBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        CudaContext& cu = *data.contexts[i];
        CudaContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void CudaParallelCalcHarmonicBondForceKernel::copyParametersToContext(ContextImpl& context, const HarmonicBondForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force);
}

class CudaParallelCalcCustomBondForceKernel::Task : public CudaContext::WorkTask {
public:
    Task(ContextImpl& context, CudaCalcCustomBondForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CudaCalcCustomBondForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

CudaParallelCalcCustomBondForceKernel::CudaParallelCalcCustomBondForceKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data, const System& system) :
        CalcCustomBondForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CudaCalcCustomBondForceKernel(name, platform, *data.contexts[i], system)));
}

void CudaParallelCalcCustomBondForceKernel::initialize(const System& system, const CustomBondForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CudaParallelCalcCustomBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        CudaContext& cu = *data.contexts[i];
        CudaContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void CudaParallelCalcCustomBondForceKernel::copyParametersToContext(ContextImpl& context, const CustomBondForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force);
}

class CudaParallelCalcHarmonicAngleForceKernel::Task : public CudaContext::WorkTask {
public:
    Task(ContextImpl& context, CudaCalcHarmonicAngleForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CudaCalcHarmonicAngleForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

CudaParallelCalcHarmonicAngleForceKernel::CudaParallelCalcHarmonicAngleForceKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data, const System& system) :
        CalcHarmonicAngleForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CudaCalcHarmonicAngleForceKernel(name, platform, *data.contexts[i], system)));
}

void CudaParallelCalcHarmonicAngleForceKernel::initialize(const System& system, const HarmonicAngleForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CudaParallelCalcHarmonicAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        CudaContext& cu = *data.contexts[i];
        CudaContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void CudaParallelCalcHarmonicAngleForceKernel::copyParametersToContext(ContextImpl& context, const HarmonicAngleForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force);
}

class CudaParallelCalcCustomAngleForceKernel::Task : public CudaContext::WorkTask {
public:
    Task(ContextImpl& context, CudaCalcCustomAngleForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CudaCalcCustomAngleForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

CudaParallelCalcCustomAngleForceKernel::CudaParallelCalcCustomAngleForceKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data, const System& system) :
        CalcCustomAngleForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CudaCalcCustomAngleForceKernel(name, platform, *data.contexts[i], system)));
}

void CudaParallelCalcCustomAngleForceKernel::initialize(const System& system, const CustomAngleForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CudaParallelCalcCustomAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        CudaContext& cu = *data.contexts[i];
        CudaContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void CudaParallelCalcCustomAngleForceKernel::copyParametersToContext(ContextImpl& context, const CustomAngleForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force);
}

class CudaParallelCalcPeriodicTorsionForceKernel::Task : public CudaContext::WorkTask {
public:
    Task(ContextImpl& context, CudaCalcPeriodicTorsionForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CudaCalcPeriodicTorsionForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

CudaParallelCalcPeriodicTorsionForceKernel::CudaParallelCalcPeriodicTorsionForceKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data, const System& system) :
        CalcPeriodicTorsionForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CudaCalcPeriodicTorsionForceKernel(name, platform, *data.contexts[i], system)));
}

void CudaParallelCalcPeriodicTorsionForceKernel::initialize(const System& system, const PeriodicTorsionForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CudaParallelCalcPeriodicTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        CudaContext& cu = *data.contexts[i];
        CudaContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void CudaParallelCalcPeriodicTorsionForceKernel::copyParametersToContext(ContextImpl& context, const PeriodicTorsionForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force);
}

class CudaParallelCalcRBTorsionForceKernel::Task : public CudaContext::WorkTask {
public:
    Task(ContextImpl& context, CudaCalcRBTorsionForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CudaCalcRBTorsionForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

CudaParallelCalcRBTorsionForceKernel::CudaParallelCalcRBTorsionForceKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data, const System& system) :
        CalcRBTorsionForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CudaCalcRBTorsionForceKernel(name, platform, *data.contexts[i], system)));
}

void CudaParallelCalcRBTorsionForceKernel::initialize(const System& system, const RBTorsionForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CudaParallelCalcRBTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        CudaContext& cu = *data.contexts[i];
        CudaContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void CudaParallelCalcRBTorsionForceKernel::copyParametersToContext(ContextImpl& context, const RBTorsionForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force);
}

class CudaParallelCalcCMAPTorsionForceKernel::Task : public CudaContext::WorkTask {
public:
    Task(ContextImpl& context, CudaCalcCMAPTorsionForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CudaCalcCMAPTorsionForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

CudaParallelCalcCMAPTorsionForceKernel::CudaParallelCalcCMAPTorsionForceKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data, const System& system) :
        CalcCMAPTorsionForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CudaCalcCMAPTorsionForceKernel(name, platform, *data.contexts[i], system)));
}

void CudaParallelCalcCMAPTorsionForceKernel::initialize(const System& system, const CMAPTorsionForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CudaParallelCalcCMAPTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        CudaContext& cu = *data.contexts[i];
        CudaContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void CudaParallelCalcCMAPTorsionForceKernel::copyParametersToContext(ContextImpl& context, const CMAPTorsionForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force);
}

class CudaParallelCalcCustomTorsionForceKernel::Task : public CudaContext::WorkTask {
public:
    Task(ContextImpl& context, CudaCalcCustomTorsionForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CudaCalcCustomTorsionForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

CudaParallelCalcCustomTorsionForceKernel::CudaParallelCalcCustomTorsionForceKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data, const System& system) :
        CalcCustomTorsionForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CudaCalcCustomTorsionForceKernel(name, platform, *data.contexts[i], system)));
}

void CudaParallelCalcCustomTorsionForceKernel::initialize(const System& system, const CustomTorsionForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CudaParallelCalcCustomTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        CudaContext& cu = *data.contexts[i];
        CudaContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void CudaParallelCalcCustomTorsionForceKernel::copyParametersToContext(ContextImpl& context, const CustomTorsionForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force);
}

class CudaParallelCalcNonbondedForceKernel::Task : public CudaContext::WorkTask {
public:
    Task(ContextImpl& context, CudaCalcNonbondedForceKernel& kernel, bool includeForce,
            bool includeEnergy, bool includeDirect, bool includeReciprocal, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), includeDirect(includeDirect), includeReciprocal(includeReciprocal), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy, includeDirect, includeReciprocal);
    }
private:
    ContextImpl& context;
    CudaCalcNonbondedForceKernel& kernel;
    bool includeForce, includeEnergy, includeDirect, includeReciprocal;
    double& energy;
};

CudaParallelCalcNonbondedForceKernel::CudaParallelCalcNonbondedForceKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data, const System& system) :
        CalcNonbondedForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CudaCalcNonbondedForceKernel(name, platform, *data.contexts[i], system)));
}

void CudaParallelCalcNonbondedForceKernel::initialize(const System& system, const NonbondedForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CudaParallelCalcNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy, bool includeDirect, bool includeReciprocal) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        CudaContext& cu = *data.contexts[i];
        CudaContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, includeDirect, includeReciprocal, data.contextEnergy[i]));
    }
    return 0.0;
}

void CudaParallelCalcNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const NonbondedForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force);
}

void CudaParallelCalcNonbondedForceKernel::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    dynamic_cast<const CudaCalcNonbondedForceKernel&>(kernels[0].getImpl()).getPMEParameters(alpha, nx, ny, nz);
}

class CudaParallelCalcCustomNonbondedForceKernel::Task : public CudaContext::WorkTask {
public:
    Task(ContextImpl& context, CudaCalcCustomNonbondedForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CudaCalcCustomNonbondedForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

CudaParallelCalcCustomNonbondedForceKernel::CudaParallelCalcCustomNonbondedForceKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data, const System& system) :
        CalcCustomNonbondedForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CudaCalcCustomNonbondedForceKernel(name, platform, *data.contexts[i], system)));
}

void CudaParallelCalcCustomNonbondedForceKernel::initialize(const System& system, const CustomNonbondedForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CudaParallelCalcCustomNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        CudaContext& cu = *data.contexts[i];
        CudaContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void CudaParallelCalcCustomNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const CustomNonbondedForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force);
}

class CudaParallelCalcCustomExternalForceKernel::Task : public CudaContext::WorkTask {
public:
    Task(ContextImpl& context, CudaCalcCustomExternalForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CudaCalcCustomExternalForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

CudaParallelCalcCustomExternalForceKernel::CudaParallelCalcCustomExternalForceKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data, const System& system) :
        CalcCustomExternalForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CudaCalcCustomExternalForceKernel(name, platform, *data.contexts[i], system)));
}

void CudaParallelCalcCustomExternalForceKernel::initialize(const System& system, const CustomExternalForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CudaParallelCalcCustomExternalForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        CudaContext& cu = *data.contexts[i];
        CudaContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void CudaParallelCalcCustomExternalForceKernel::copyParametersToContext(ContextImpl& context, const CustomExternalForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force);
}

class CudaParallelCalcCustomHbondForceKernel::Task : public CudaContext::WorkTask {
public:
    Task(ContextImpl& context, CudaCalcCustomHbondForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CudaCalcCustomHbondForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

CudaParallelCalcCustomHbondForceKernel::CudaParallelCalcCustomHbondForceKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data, const System& system) :
        CalcCustomHbondForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CudaCalcCustomHbondForceKernel(name, platform, *data.contexts[i], system)));
}

void CudaParallelCalcCustomHbondForceKernel::initialize(const System& system, const CustomHbondForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CudaParallelCalcCustomHbondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        CudaContext& cu = *data.contexts[i];
        CudaContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void CudaParallelCalcCustomHbondForceKernel::copyParametersToContext(ContextImpl& context, const CustomHbondForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force);
}

class CudaParallelCalcCustomCompoundBondForceKernel::Task : public CudaContext::WorkTask {
public:
    Task(ContextImpl& context, CudaCalcCustomCompoundBondForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CudaCalcCustomCompoundBondForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

CudaParallelCalcCustomCompoundBondForceKernel::CudaParallelCalcCustomCompoundBondForceKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data, const System& system) :
        CalcCustomCompoundBondForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CudaCalcCustomCompoundBondForceKernel(name, platform, *data.contexts[i], system)));
}

void CudaParallelCalcCustomCompoundBondForceKernel::initialize(const System& system, const CustomCompoundBondForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CudaParallelCalcCustomCompoundBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        CudaContext& cu = *data.contexts[i];
        CudaContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void CudaParallelCalcCustomCompoundBondForceKernel::copyParametersToContext(ContextImpl& context, const CustomCompoundBondForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force);
}
