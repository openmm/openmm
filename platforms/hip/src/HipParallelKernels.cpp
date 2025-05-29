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
#include "openmm/internal/timer.h"

using namespace OpenMM;
using namespace std;


#define CHECK_RESULT(result, prefix) \
if (result != hipSuccess) { \
    std::stringstream m; \
    m<<prefix<<": "<<cu.getErrorString(result)<<" ("<<result<<")"<<" at "<<__FILE__<<":"<<__LINE__; \
    throw OpenMMException(m.str());\
}

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
            bool includeForce, bool includeEnergy, int groups, double& energy, double& completionTime, long long* pinnedMemory, HipArray& contextForces,
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
            completionTime = getCurrentTime();
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
    double& completionTime;
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
