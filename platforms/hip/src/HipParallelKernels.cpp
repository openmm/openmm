/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011-2024 Stanford University and the Authors.      *
 * Portions copyright (c) 2020-2021 Advanced Micro Devices, Inc.              *
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

#include "HipParallelKernels.h"
#include "HipKernelSources.h"
#include "openmm/common/ContextSelector.h"

using namespace OpenMM;
using namespace std;


#define CHECK_RESULT(result, prefix) \
if (result != hipSuccess) { \
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

class HipParallelCalcForcesAndEnergyKernel::BeginComputationTask : public HipContext::WorkTask {
public:
    BeginComputationTask(ContextImpl& context, HipContext& cu, HipCalcForcesAndEnergyKernel& kernel,
            bool includeForce, bool includeEnergy, int groups, void* pinnedMemory, hipEvent_t event) : context(context), cu(cu), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), groups(groups), pinnedMemory(pinnedMemory), event(event) {
    }
    void execute() {
        // Copy coordinates over to this device and execute the kernel.

        ContextSelector selector(cu);
        if (cu.getContextIndex() > 0) {
            hipStreamWaitEvent(cu.getCurrentStream(), event, 0);
            if (!cu.getPlatformData().peerAccessSupported)
                cu.getPosq().upload(pinnedMemory, false);
        }
        kernel.beginComputation(context, includeForce, includeEnergy, groups);
    }
private:
    ContextImpl& context;
    HipContext& cu;
    HipCalcForcesAndEnergyKernel& kernel;
    bool includeForce, includeEnergy;
    int groups;
    void* pinnedMemory;
    hipEvent_t event;
};

class HipParallelCalcForcesAndEnergyKernel::FinishComputationTask : public HipContext::WorkTask {
public:
    FinishComputationTask(ContextImpl& context, HipContext& cu, HipCalcForcesAndEnergyKernel& kernel,
            bool includeForce, bool includeEnergy, int groups, double& energy, long long& completionTime, long long* pinnedMemory, HipArray& contextForces,
            bool& valid, hipStream_t stream, hipEvent_t event, hipEvent_t localEvent) :
            context(context), cu(cu), kernel(kernel), includeForce(includeForce), includeEnergy(includeEnergy), groups(groups), energy(energy),
            completionTime(completionTime), pinnedMemory(pinnedMemory), contextForces(contextForces), valid(valid),
            stream(stream), event(event), localEvent(localEvent) {
    }
    void execute() {
        // Execute the kernel, then download forces.

        ContextSelector selector(cu);
        energy += kernel.finishComputation(context, includeForce, includeEnergy, groups, valid);
        if (cu.getComputeForceCount() < 200) {
            // Record timing information for load balancing.  Since this takes time, only do it at the start of the simulation.

            CHECK_RESULT(hipStreamSynchronize(cu.getCurrentStream()), "Error synchronizing HIP context");
            completionTime = getTime();
        }
        if (includeForce) {
            if (cu.getContextIndex() > 0) {
                hipEventRecord(localEvent, cu.getCurrentStream());
                hipStreamWaitEvent(stream, localEvent, 0);
                int numAtoms = cu.getPaddedNumAtoms();
                if (cu.getPlatformData().peerAccessSupported) {
                    int numBytes = numAtoms*3*sizeof(long long);
                    int offset = (cu.getContextIndex()-1)*numBytes;
                    CHECK_RESULT(hipMemcpyAsync(static_cast<char*>(contextForces.getDevicePointer())+offset,
                                           cu.getForce().getDevicePointer(), numBytes, hipMemcpyDeviceToDevice, stream), "Error copying forces");
                    hipEventRecord(event, stream);
                }
                else
                    cu.getForce().download(&pinnedMemory[(cu.getContextIndex()-1)*numAtoms*3]);
            }
        }
    }
private:
    ContextImpl& context;
    HipContext& cu;
    HipCalcForcesAndEnergyKernel& kernel;
    bool includeForce, includeEnergy;
    int groups;
    double& energy;
    long long& completionTime;
    long long* pinnedMemory;
    HipArray& contextForces;
    bool& valid;
    hipStream_t stream;
    hipEvent_t event;
    hipEvent_t localEvent;
};

HipParallelCalcForcesAndEnergyKernel::HipParallelCalcForcesAndEnergyKernel(string name, const Platform& platform, HipPlatform::PlatformData& data) :
        CalcForcesAndEnergyKernel(name, platform), data(data), completionTimes(data.contexts.size()), contextNonbondedFractions(data.contexts.size()),
        pinnedPositionBuffer(NULL), pinnedForceBuffer(NULL) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new HipCalcForcesAndEnergyKernel(name, platform, *data.contexts[i])));
}

HipParallelCalcForcesAndEnergyKernel::~HipParallelCalcForcesAndEnergyKernel() {
    ContextSelector selector(*data.contexts[0]);
    if (pinnedPositionBuffer != NULL)
        hipHostFree(pinnedPositionBuffer);
    if (pinnedForceBuffer != NULL)
        hipHostFree(pinnedForceBuffer);
    hipEventDestroy(event);
    for (int i = 0; i < peerCopyEvent.size(); i++)
        hipEventDestroy(peerCopyEvent[i]);
    for (int i = 0; i < peerCopyEventLocal.size(); i++)
        hipEventDestroy(peerCopyEventLocal[i]);
    for (int i = 0; i < peerCopyStream.size(); i++)
        hipStreamDestroy(peerCopyStream[i]);
}

void HipParallelCalcForcesAndEnergyKernel::initialize(const System& system) {
    HipContext& cu = *data.contexts[0];
    ContextSelector selector(cu);
    hipModule_t module = cu.createModule(HipKernelSources::parallel);
    sumKernel = cu.getKernel(module, "sumForces");
    int numContexts = data.contexts.size();
    for (int i = 0; i < numContexts; i++)
        getKernel(i).initialize(system);
    for (int i = 0; i < numContexts; i++)
        contextNonbondedFractions[i] = 1/(double) numContexts;
    CHECK_RESULT(hipEventCreateWithFlags(&event, cu.getEventFlags()), "Error creating event");
    peerCopyEvent.resize(numContexts);
    peerCopyEventLocal.resize(numContexts);
    peerCopyStream.resize(numContexts);
    for (int i = 0; i < numContexts; i++) {
        HipContext& cuLocal = *data.contexts[i];
        ContextSelector selectorLocal(cuLocal);
        CHECK_RESULT(hipEventCreateWithFlags(&peerCopyEvent[i], cu.getEventFlags()), "Error creating event");
        CHECK_RESULT(hipStreamCreateWithFlags(&peerCopyStream[i], hipStreamNonBlocking), "Error creating stream");
        CHECK_RESULT(hipEventCreateWithFlags(&peerCopyEventLocal[i], cu.getEventFlags()), "Error creating event");
    }
}

void HipParallelCalcForcesAndEnergyKernel::beginComputation(ContextImpl& context, bool includeForce, bool includeEnergy, int groups) {
    HipContext& cu = *data.contexts[0];
    ContextSelector selector(cu);
    if (!contextForces.isInitialized()) {
        contextForces.initialize<long long>(cu, 3*(data.contexts.size()-1)*cu.getPaddedNumAtoms(), "contextForces");
        if (!cu.getPlatformData().peerAccessSupported) {
            CHECK_RESULT(hipHostMalloc((void**) &pinnedForceBuffer, 3*(data.contexts.size()-1)*cu.getPaddedNumAtoms()*sizeof(long long), hipHostMallocPortable), "Error allocating pinned memory");
            CHECK_RESULT(hipHostMalloc(&pinnedPositionBuffer, cu.getPaddedNumAtoms()*(cu.getUseDoublePrecision() ? sizeof(double4) : sizeof(float4)), hipHostMallocPortable), "Error allocating pinned memory");
        }
    }

    // Copy coordinates over to each device and execute the kernel.

    if (!cu.getPlatformData().peerAccessSupported) {
        cu.getPosq().download(pinnedPositionBuffer, false);
        hipEventRecord(event, cu.getCurrentStream());
    }
    else {
        int numBytes = cu.getPosq().getSize()*cu.getPosq().getElementSize();
        hipEventRecord(event, cu.getCurrentStream());
        for (int i = 1; i < (int) data.contexts.size(); i++) {
            hipStreamWaitEvent(peerCopyStream[i], event, 0);
            CHECK_RESULT(hipMemcpyAsync(
                data.contexts[i]->getPosq().getDevicePointer(),
                cu.getPosq().getDevicePointer(), numBytes,
                hipMemcpyDeviceToDevice, peerCopyStream[i]), "Error copying positions");
            hipEventRecord(peerCopyEvent[i], peerCopyStream[i]);
        }
    }
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        data.contextEnergy[i] = 0.0;
        HipContext& cu = *data.contexts[i];
        ComputeContext::WorkThread& thread = cu.getWorkThread();
        hipEvent_t waitEvent = (cu.getPlatformData().peerAccessSupported ? peerCopyEvent[i] : event);
        thread.addTask(new BeginComputationTask(context, cu, getKernel(i), includeForce, includeEnergy, groups, pinnedPositionBuffer, waitEvent));
    }
    data.syncContexts();
}

double HipParallelCalcForcesAndEnergyKernel::finishComputation(ContextImpl& context, bool includeForce, bool includeEnergy, int groups, bool& valid) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        HipContext& cu = *data.contexts[i];
        ComputeContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new FinishComputationTask(context, cu, getKernel(i), includeForce, includeEnergy, groups, data.contextEnergy[i], completionTimes[i],
                pinnedForceBuffer, contextForces, valid, peerCopyStream[i], peerCopyEvent[i], peerCopyEventLocal[i]));
    }
    data.syncContexts();
    HipContext& cu = *data.contexts[0];
    ContextSelector selector(cu);
    if (cu.getPlatformData().peerAccessSupported)
        for (int i = 1; i < data.contexts.size(); i++)
            hipStreamWaitEvent(cu.getCurrentStream(), peerCopyEvent[i], 0);
    double energy = 0.0;
    for (int i = 0; i < (int) data.contextEnergy.size(); i++)
        energy += data.contextEnergy[i];
    if (includeForce && valid) {
        // Sum the forces from all devices.

        if (!cu.getPlatformData().peerAccessSupported)
            contextForces.upload(pinnedForceBuffer, false);
        int bufferSize = 3*cu.getPaddedNumAtoms();
        int numBuffers = data.contexts.size()-1;
        void* args[] = {&cu.getForce().getDevicePointer(), &contextForces.getDevicePointer(), &bufferSize, &numBuffers};
        cu.executeKernel(sumKernel, args, bufferSize);

        // Balance work between the contexts by transferring a little nonbonded work from the context that
        // finished last to the one that finished first.

        if (cu.getComputeForceCount() < 200) {
            int firstIndex = 0, lastIndex = 0;
            const double eps = 0.001;
            for (int i = 0; i < (int) completionTimes.size(); i++) {
                if (completionTimes[i] < completionTimes[firstIndex])
                    firstIndex = i;
                if (contextNonbondedFractions[lastIndex] < eps || completionTimes[i] > completionTimes[lastIndex])
                    lastIndex = i;
            }
            double fractionToTransfer = min(cu.getComputeForceCount() < 100 ? 0.01 : 0.001, contextNonbondedFractions[lastIndex]);
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

class HipParallelCalcHarmonicBondForceKernel::Task : public HipContext::WorkTask {
public:
    Task(ContextImpl& context, CommonCalcHarmonicBondForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CommonCalcHarmonicBondForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

HipParallelCalcHarmonicBondForceKernel::HipParallelCalcHarmonicBondForceKernel(std::string name, const Platform& platform, HipPlatform::PlatformData& data, const System& system) :
        CalcHarmonicBondForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CommonCalcHarmonicBondForceKernel(name, platform, *data.contexts[i], system)));
}

void HipParallelCalcHarmonicBondForceKernel::initialize(const System& system, const HarmonicBondForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double HipParallelCalcHarmonicBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        HipContext& cu = *data.contexts[i];
        ComputeContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void HipParallelCalcHarmonicBondForceKernel::copyParametersToContext(ContextImpl& context, const HarmonicBondForce& force, int firstBond, int lastBond) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force, firstBond, lastBond);
}

class HipParallelCalcCustomBondForceKernel::Task : public HipContext::WorkTask {
public:
    Task(ContextImpl& context, CommonCalcCustomBondForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CommonCalcCustomBondForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

HipParallelCalcCustomBondForceKernel::HipParallelCalcCustomBondForceKernel(std::string name, const Platform& platform, HipPlatform::PlatformData& data, const System& system) :
        CalcCustomBondForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CommonCalcCustomBondForceKernel(name, platform, *data.contexts[i], system)));
}

void HipParallelCalcCustomBondForceKernel::initialize(const System& system, const CustomBondForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double HipParallelCalcCustomBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        HipContext& cu = *data.contexts[i];
        ComputeContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void HipParallelCalcCustomBondForceKernel::copyParametersToContext(ContextImpl& context, const CustomBondForce& force, int firstBond, int lastBond) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force, firstBond, lastBond);
}

class HipParallelCalcHarmonicAngleForceKernel::Task : public HipContext::WorkTask {
public:
    Task(ContextImpl& context, CommonCalcHarmonicAngleForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CommonCalcHarmonicAngleForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

HipParallelCalcHarmonicAngleForceKernel::HipParallelCalcHarmonicAngleForceKernel(std::string name, const Platform& platform, HipPlatform::PlatformData& data, const System& system) :
        CalcHarmonicAngleForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CommonCalcHarmonicAngleForceKernel(name, platform, *data.contexts[i], system)));
}

void HipParallelCalcHarmonicAngleForceKernel::initialize(const System& system, const HarmonicAngleForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double HipParallelCalcHarmonicAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        HipContext& cu = *data.contexts[i];
        ComputeContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void HipParallelCalcHarmonicAngleForceKernel::copyParametersToContext(ContextImpl& context, const HarmonicAngleForce& force, int firstAngle, int lastAngle) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force, firstAngle, lastAngle);
}

class HipParallelCalcCustomAngleForceKernel::Task : public HipContext::WorkTask {
public:
    Task(ContextImpl& context, CommonCalcCustomAngleForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CommonCalcCustomAngleForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

HipParallelCalcCustomAngleForceKernel::HipParallelCalcCustomAngleForceKernel(std::string name, const Platform& platform, HipPlatform::PlatformData& data, const System& system) :
        CalcCustomAngleForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CommonCalcCustomAngleForceKernel(name, platform, *data.contexts[i], system)));
}

void HipParallelCalcCustomAngleForceKernel::initialize(const System& system, const CustomAngleForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double HipParallelCalcCustomAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        HipContext& cu = *data.contexts[i];
        ComputeContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void HipParallelCalcCustomAngleForceKernel::copyParametersToContext(ContextImpl& context, const CustomAngleForce& force, int firstAngle, int lastAngle) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force, firstAngle, lastAngle);
}

class HipParallelCalcPeriodicTorsionForceKernel::Task : public HipContext::WorkTask {
public:
    Task(ContextImpl& context, CommonCalcPeriodicTorsionForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CommonCalcPeriodicTorsionForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

HipParallelCalcPeriodicTorsionForceKernel::HipParallelCalcPeriodicTorsionForceKernel(std::string name, const Platform& platform, HipPlatform::PlatformData& data, const System& system) :
        CalcPeriodicTorsionForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CommonCalcPeriodicTorsionForceKernel(name, platform, *data.contexts[i], system)));
}

void HipParallelCalcPeriodicTorsionForceKernel::initialize(const System& system, const PeriodicTorsionForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double HipParallelCalcPeriodicTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        HipContext& cu = *data.contexts[i];
        ComputeContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void HipParallelCalcPeriodicTorsionForceKernel::copyParametersToContext(ContextImpl& context, const PeriodicTorsionForce& force, int firstTorsion, int lastTorsion) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force, firstTorsion, lastTorsion);
}

class HipParallelCalcRBTorsionForceKernel::Task : public HipContext::WorkTask {
public:
    Task(ContextImpl& context, CommonCalcRBTorsionForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CommonCalcRBTorsionForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

HipParallelCalcRBTorsionForceKernel::HipParallelCalcRBTorsionForceKernel(std::string name, const Platform& platform, HipPlatform::PlatformData& data, const System& system) :
        CalcRBTorsionForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CommonCalcRBTorsionForceKernel(name, platform, *data.contexts[i], system)));
}

void HipParallelCalcRBTorsionForceKernel::initialize(const System& system, const RBTorsionForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double HipParallelCalcRBTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        HipContext& cu = *data.contexts[i];
        ComputeContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void HipParallelCalcRBTorsionForceKernel::copyParametersToContext(ContextImpl& context, const RBTorsionForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force);
}

class HipParallelCalcCMAPTorsionForceKernel::Task : public HipContext::WorkTask {
public:
    Task(ContextImpl& context, CommonCalcCMAPTorsionForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CommonCalcCMAPTorsionForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

HipParallelCalcCMAPTorsionForceKernel::HipParallelCalcCMAPTorsionForceKernel(std::string name, const Platform& platform, HipPlatform::PlatformData& data, const System& system) :
        CalcCMAPTorsionForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CommonCalcCMAPTorsionForceKernel(name, platform, *data.contexts[i], system)));
}

void HipParallelCalcCMAPTorsionForceKernel::initialize(const System& system, const CMAPTorsionForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double HipParallelCalcCMAPTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        HipContext& cu = *data.contexts[i];
        ComputeContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void HipParallelCalcCMAPTorsionForceKernel::copyParametersToContext(ContextImpl& context, const CMAPTorsionForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force);
}

class HipParallelCalcCustomTorsionForceKernel::Task : public HipContext::WorkTask {
public:
    Task(ContextImpl& context, CommonCalcCustomTorsionForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CommonCalcCustomTorsionForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

HipParallelCalcCustomTorsionForceKernel::HipParallelCalcCustomTorsionForceKernel(std::string name, const Platform& platform, HipPlatform::PlatformData& data, const System& system) :
        CalcCustomTorsionForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CommonCalcCustomTorsionForceKernel(name, platform, *data.contexts[i], system)));
}

void HipParallelCalcCustomTorsionForceKernel::initialize(const System& system, const CustomTorsionForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double HipParallelCalcCustomTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        HipContext& cu = *data.contexts[i];
        ComputeContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void HipParallelCalcCustomTorsionForceKernel::copyParametersToContext(ContextImpl& context, const CustomTorsionForce& force, int firstTorsion, int lastTorsion) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force, firstTorsion, lastTorsion);
}

class HipParallelCalcNonbondedForceKernel::Task : public HipContext::WorkTask {
public:
    Task(ContextImpl& context, HipCalcNonbondedForceKernel& kernel, bool includeForce,
            bool includeEnergy, bool includeDirect, bool includeReciprocal, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), includeDirect(includeDirect), includeReciprocal(includeReciprocal), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy, includeDirect, includeReciprocal);
    }
private:
    ContextImpl& context;
    HipCalcNonbondedForceKernel& kernel;
    bool includeForce, includeEnergy, includeDirect, includeReciprocal;
    double& energy;
};

HipParallelCalcNonbondedForceKernel::HipParallelCalcNonbondedForceKernel(std::string name, const Platform& platform, HipPlatform::PlatformData& data, const System& system) :
        CalcNonbondedForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new HipCalcNonbondedForceKernel(name, platform, *data.contexts[i], system)));
}

void HipParallelCalcNonbondedForceKernel::initialize(const System& system, const NonbondedForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double HipParallelCalcNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy, bool includeDirect, bool includeReciprocal) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        HipContext& cu = *data.contexts[i];
        ComputeContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, includeDirect, includeReciprocal, data.contextEnergy[i]));
    }
    return 0.0;
}

void HipParallelCalcNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const NonbondedForce& force, int firstParticle, int lastParticle, int firstException, int lastException) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force, firstParticle, lastParticle, firstException, lastException);
}

void HipParallelCalcNonbondedForceKernel::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    dynamic_cast<const HipCalcNonbondedForceKernel&>(kernels[0].getImpl()).getPMEParameters(alpha, nx, ny, nz);
}

void HipParallelCalcNonbondedForceKernel::getLJPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    dynamic_cast<const HipCalcNonbondedForceKernel&>(kernels[0].getImpl()).getLJPMEParameters(alpha, nx, ny, nz);
}

class HipParallelCalcCustomNonbondedForceKernel::Task : public HipContext::WorkTask {
public:
    Task(ContextImpl& context, CommonCalcCustomNonbondedForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CommonCalcCustomNonbondedForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

HipParallelCalcCustomNonbondedForceKernel::HipParallelCalcCustomNonbondedForceKernel(std::string name, const Platform& platform, HipPlatform::PlatformData& data, const System& system) :
        CalcCustomNonbondedForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CommonCalcCustomNonbondedForceKernel(name, platform, *data.contexts[i], system)));
}

void HipParallelCalcCustomNonbondedForceKernel::initialize(const System& system, const CustomNonbondedForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double HipParallelCalcCustomNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        HipContext& cu = *data.contexts[i];
        ComputeContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void HipParallelCalcCustomNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const CustomNonbondedForce& force, int firstParticle, int lastParticle) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force, firstParticle, lastParticle);
}

class HipParallelCalcCustomExternalForceKernel::Task : public HipContext::WorkTask {
public:
    Task(ContextImpl& context, CommonCalcCustomExternalForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CommonCalcCustomExternalForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

HipParallelCalcCustomExternalForceKernel::HipParallelCalcCustomExternalForceKernel(std::string name, const Platform& platform, HipPlatform::PlatformData& data, const System& system) :
        CalcCustomExternalForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CommonCalcCustomExternalForceKernel(name, platform, *data.contexts[i], system)));
}

void HipParallelCalcCustomExternalForceKernel::initialize(const System& system, const CustomExternalForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double HipParallelCalcCustomExternalForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        HipContext& cu = *data.contexts[i];
        ComputeContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void HipParallelCalcCustomExternalForceKernel::copyParametersToContext(ContextImpl& context, const CustomExternalForce& force, int firstParticle, int lastParticle) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force, firstParticle, lastParticle);
}

class HipParallelCalcCustomHbondForceKernel::Task : public HipContext::WorkTask {
public:
    Task(ContextImpl& context, CommonCalcCustomHbondForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CommonCalcCustomHbondForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

HipParallelCalcCustomHbondForceKernel::HipParallelCalcCustomHbondForceKernel(std::string name, const Platform& platform, HipPlatform::PlatformData& data, const System& system) :
        CalcCustomHbondForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CommonCalcCustomHbondForceKernel(name, platform, *data.contexts[i], system)));
}

void HipParallelCalcCustomHbondForceKernel::initialize(const System& system, const CustomHbondForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double HipParallelCalcCustomHbondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        HipContext& cu = *data.contexts[i];
        ComputeContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void HipParallelCalcCustomHbondForceKernel::copyParametersToContext(ContextImpl& context, const CustomHbondForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force);
}

class HipParallelCalcCustomCompoundBondForceKernel::Task : public HipContext::WorkTask {
public:
    Task(ContextImpl& context, CommonCalcCustomCompoundBondForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CommonCalcCustomCompoundBondForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

HipParallelCalcCustomCompoundBondForceKernel::HipParallelCalcCustomCompoundBondForceKernel(std::string name, const Platform& platform, HipPlatform::PlatformData& data, const System& system) :
        CalcCustomCompoundBondForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CommonCalcCustomCompoundBondForceKernel(name, platform, *data.contexts[i], system)));
}

void HipParallelCalcCustomCompoundBondForceKernel::initialize(const System& system, const CustomCompoundBondForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double HipParallelCalcCustomCompoundBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        HipContext& cu = *data.contexts[i];
        ComputeContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void HipParallelCalcCustomCompoundBondForceKernel::copyParametersToContext(ContextImpl& context, const CustomCompoundBondForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force);
}
