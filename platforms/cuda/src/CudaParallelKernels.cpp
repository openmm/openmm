/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011-2024 Stanford University and the Authors.      *
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
#include "openmm/common/ContextSelector.h"
#include "openmm/internal/timer.h"

using namespace OpenMM;
using namespace std;


#define CHECK_RESULT(result, prefix) \
if (result != CUDA_SUCCESS) { \
    std::stringstream m; \
    m<<prefix<<": "<<cu.getErrorString(result)<<" ("<<result<<")"<<" at "<<__FILE__<<":"<<__LINE__; \
    throw OpenMMException(m.str());\
}

class CudaParallelCalcForcesAndEnergyKernel::BeginComputationTask : public CudaContext::WorkTask {
public:
    BeginComputationTask(ContextImpl& context, CudaContext& cu, CudaCalcForcesAndEnergyKernel& kernel,
            bool includeForce, bool includeEnergy, int groups, void* pinnedMemory, CUevent event, int2& interactionCount) : context(context), cu(cu), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), groups(groups), pinnedMemory(pinnedMemory), event(event), interactionCount(interactionCount) {
    }
    void execute() {
        // Copy coordinates over to this device and execute the kernel.

        ContextSelector selector(cu);
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
            bool includeForce, bool includeEnergy, int groups, double& energy, double& completionTime, long long* pinnedMemory, CudaArray& contextForces,
            bool& valid, int2& interactionCount, CUstream stream, CUevent event, CUevent localEvent, bool loadBalance) :
            context(context), cu(cu), kernel(kernel), includeForce(includeForce), includeEnergy(includeEnergy), groups(groups), energy(energy),
            completionTime(completionTime), pinnedMemory(pinnedMemory), contextForces(contextForces), valid(valid), interactionCount(interactionCount),
            stream(stream), event(event), localEvent(localEvent), loadBalance(loadBalance) {
    }
    void execute() {
        // Execute the kernel, then download forces.
        
        ContextSelector selector(cu);
        energy += kernel.finishComputation(context, includeForce, includeEnergy, groups, valid);
        if (loadBalance) {
            // Record timing information for load balancing.  Since this takes time, only do it at the start of the simulation.

            CHECK_RESULT(cuCtxSynchronize(), "Error synchronizing CUDA context");
            completionTime = getCurrentTime();
        }
        if (includeForce) {
            if (cu.getContextIndex() > 0) {
                cuEventRecord(localEvent, cu.getCurrentStream());
                cuStreamWaitEvent(stream, localEvent, 0);
                int numAtoms = cu.getPaddedNumAtoms();
                if (cu.getPlatformData().peerAccessSupported) {
                    int numBytes = numAtoms*3*sizeof(long long);
                    int offset = (cu.getContextIndex()-1)*numBytes;
                    CudaContext& context0 = *cu.getPlatformData().contexts[0];
                    CHECK_RESULT(cuMemcpyAsync(contextForces.getDevicePointer()+offset, cu.getForce().getDevicePointer(), numBytes, stream), "Error copying forces");
                    cuEventRecord(event, stream);
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
    bool includeForce, includeEnergy, loadBalance;
    int groups;
    double& energy;
    double& completionTime;
    long long* pinnedMemory;
    CudaArray& contextForces;
    bool& valid;
    int2& interactionCount;
    CUstream stream;
    CUevent event;
    CUevent localEvent;
};

CudaParallelCalcForcesAndEnergyKernel::CudaParallelCalcForcesAndEnergyKernel(string name, const Platform& platform, CudaPlatform::PlatformData& data) :
        CalcForcesAndEnergyKernel(name, platform), data(data), completionTimes(data.contexts.size()), contextNonbondedFractions(data.contexts.size()),
        interactionCounts(NULL), pinnedPositionBuffer(NULL), pinnedForceBuffer(NULL) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CudaCalcForcesAndEnergyKernel(name, platform, *data.contexts[i])));
}

CudaParallelCalcForcesAndEnergyKernel::~CudaParallelCalcForcesAndEnergyKernel() {
    ContextSelector selector(*data.contexts[0]);
    if (pinnedPositionBuffer != NULL)
        cuMemFreeHost(pinnedPositionBuffer);
    if (pinnedForceBuffer != NULL)
        cuMemFreeHost(pinnedForceBuffer);
    cuEventDestroy(event);
    for (int i = 0; i < peerCopyEvent.size(); i++)
        cuEventDestroy(peerCopyEvent[i]);
    for (int i = 0; i < peerCopyEventLocal.size(); i++)
        cuEventDestroy(peerCopyEventLocal[i]);
    for (int i = 0; i < peerCopyStream.size(); i++)
        cuStreamDestroy(peerCopyStream[i]);
    if (interactionCounts != NULL)
        cuMemFreeHost(interactionCounts);
}

void CudaParallelCalcForcesAndEnergyKernel::initialize(const System& system) {
    CudaContext& cu = *data.contexts[0];
    ContextSelector selector(cu);
    CUmodule module = cu.createModule(CudaKernelSources::parallel);
    sumKernel = cu.getKernel(module, "sumForces");
    int numContexts = data.contexts.size();
    for (int i = 0; i < numContexts; i++)
        getKernel(i).initialize(system);
    for (int i = 0; i < contextNonbondedFractions.size(); i++) {
        double x0 = i/(double) contextNonbondedFractions.size();
        double x1 = (i+1)/(double) contextNonbondedFractions.size();
        contextNonbondedFractions[i] = x1*x1 - x0*x0;
    }
    CHECK_RESULT(cuEventCreate(&event, cu.getEventFlags()), "Error creating event");
    peerCopyEvent.resize(numContexts);
    peerCopyEventLocal.resize(numContexts);
    peerCopyStream.resize(numContexts);
    for (int i = 0; i < numContexts; i++) {
        CHECK_RESULT(cuEventCreate(&peerCopyEvent[i], cu.getEventFlags()), "Error creating event");
        CHECK_RESULT(cuStreamCreate(&peerCopyStream[i], CU_STREAM_NON_BLOCKING), "Error creating stream");
    }
    for (int i = 0; i < numContexts; i++) {
        CudaContext& cuLocal = *data.contexts[i];
        ContextSelector selectorLocal(cuLocal);
        CHECK_RESULT(cuEventCreate(&peerCopyEventLocal[i], cu.getEventFlags()), "Error creating event");
    }
    CHECK_RESULT(cuMemHostAlloc((void**) &interactionCounts, numContexts*sizeof(int2), 0), "Error creating interaction counts buffer");
}

void CudaParallelCalcForcesAndEnergyKernel::beginComputation(ContextImpl& context, bool includeForce, bool includeEnergy, int groups) {
    CudaContext& cu = *data.contexts[0];
    ContextSelector selector(cu);
    if (!contextForces.isInitialized()) {
        contextForces.initialize<long long>(cu, 3*(data.contexts.size()-1)*cu.getPaddedNumAtoms(), "contextForces");
        CHECK_RESULT(cuMemHostAlloc((void**) &pinnedForceBuffer, 3*(data.contexts.size()-1)*cu.getPaddedNumAtoms()*sizeof(long long), CU_MEMHOSTALLOC_PORTABLE), "Error allocating pinned memory");
        CHECK_RESULT(cuMemHostAlloc(&pinnedPositionBuffer, cu.getPaddedNumAtoms()*(cu.getUseDoublePrecision() ? sizeof(double4) : sizeof(float4)), CU_MEMHOSTALLOC_PORTABLE), "Error allocating pinned memory");
    }
    loadBalance = (cu.getComputeForceCount() < 200 || cu.getComputeForceCount()%30 == 0);

    // Copy coordinates over to each device and execute the kernel.
    
    if (!cu.getPlatformData().peerAccessSupported) {
        cu.getPosq().download(pinnedPositionBuffer, false);
        cuEventRecord(event, cu.getCurrentStream());
    }
    else {
        int numBytes = cu.getPosq().getSize()*cu.getPosq().getElementSize();
        cuEventRecord(event, cu.getCurrentStream());
        for (int i = 1; i < (int) data.contexts.size(); i++) {
            cuStreamWaitEvent(peerCopyStream[i], event, 0);
            CHECK_RESULT(cuMemcpyAsync(data.contexts[i]->getPosq().getDevicePointer(), cu.getPosq().getDevicePointer(), numBytes, peerCopyStream[i]), "Error copying positions");
            cuEventRecord(peerCopyEvent[i], peerCopyStream[i]);
        }
    }
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        data.contextEnergy[i] = 0.0;
        CudaContext& cu = *data.contexts[i];
        ComputeContext::WorkThread& thread = cu.getWorkThread();
        CUevent waitEvent = (cu.getPlatformData().peerAccessSupported ? peerCopyEvent[i] : event);
        thread.addTask(new BeginComputationTask(context, cu, getKernel(i), includeForce, includeEnergy, groups, pinnedPositionBuffer, waitEvent, interactionCounts[i]));
    }
    data.syncContexts();
}

double CudaParallelCalcForcesAndEnergyKernel::finishComputation(ContextImpl& context, bool includeForce, bool includeEnergy, int groups, bool& valid) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        CudaContext& cu = *data.contexts[i];
        ComputeContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new FinishComputationTask(context, cu, getKernel(i), includeForce, includeEnergy, groups, data.contextEnergy[i], completionTimes[i],
                pinnedForceBuffer, contextForces, valid, interactionCounts[i], peerCopyStream[i], peerCopyEvent[i], peerCopyEventLocal[i], loadBalance));
    }
    data.syncContexts();
    CudaContext& cu = *data.contexts[0];
    ContextSelector selector(cu);
    if (cu.getPlatformData().peerAccessSupported)
        for (int i = 1; i < data.contexts.size(); i++)
            cuStreamWaitEvent(cu.getCurrentStream(), peerCopyEvent[i], 0);
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
        
        if (loadBalance) {
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
        ComputeContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, includeDirect, includeReciprocal, data.contextEnergy[i]));
    }
    return 0.0;
}

void CudaParallelCalcNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const NonbondedForce& force, int firstParticle, int lastParticle, int firstException, int lastException) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force, firstParticle, lastParticle, firstException, lastException);
}

void CudaParallelCalcNonbondedForceKernel::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    dynamic_cast<const CudaCalcNonbondedForceKernel&>(kernels[0].getImpl()).getPMEParameters(alpha, nx, ny, nz);
}

void CudaParallelCalcNonbondedForceKernel::getLJPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    dynamic_cast<const CudaCalcNonbondedForceKernel&>(kernels[0].getImpl()).getLJPMEParameters(alpha, nx, ny, nz);
}
