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

#include "OpenCLParallelKernels.h"
#include "openmm/internal/timer.h"

using namespace OpenMM;
using namespace std;

class OpenCLParallelCalcForcesAndEnergyKernel::BeginComputationTask : public OpenCLContext::WorkTask {
public:
    BeginComputationTask(ContextImpl& context, OpenCLContext& cl, OpenCLCalcForcesAndEnergyKernel& kernel,
            bool includeForce, bool includeEnergy, int groups, void* pinnedMemory, int& numTiles) : context(context), cl(cl), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), groups(groups), pinnedMemory(pinnedMemory), numTiles(numTiles) {
    }
    void execute() {
        // Copy coordinates over to this device and execute the kernel.

        if (cl.getContextIndex() > 0)
            cl.getQueue().enqueueWriteBuffer(cl.getPosq().getDeviceBuffer(), CL_FALSE, 0, cl.getPaddedNumAtoms()*cl.getPosq().getElementSize(), pinnedMemory);
        kernel.beginComputation(context, includeForce, includeEnergy, groups);
        if (cl.getNonbondedUtilities().getUsePeriodic())
            cl.getNonbondedUtilities().getInteractionCount().download(&numTiles, false);
    }
private:
    ContextImpl& context;
    OpenCLContext& cl;
    OpenCLCalcForcesAndEnergyKernel& kernel;
    bool includeForce, includeEnergy;
    int groups;
    void* pinnedMemory;
    int& numTiles;
};

class OpenCLParallelCalcForcesAndEnergyKernel::FinishComputationTask : public OpenCLContext::WorkTask {
public:
    FinishComputationTask(ContextImpl& context, OpenCLContext& cl, OpenCLCalcForcesAndEnergyKernel& kernel,
            bool includeForce, bool includeEnergy, int groups, double& energy, double& completionTime, void* pinnedMemory, bool& valid, int& numTiles) :
            context(context), cl(cl), kernel(kernel), includeForce(includeForce), includeEnergy(includeEnergy), groups(groups), energy(energy),
            completionTime(completionTime), pinnedMemory(pinnedMemory), valid(valid), numTiles(numTiles) {
    }
    void execute() {
        // Execute the kernel, then download forces.
        
        energy += kernel.finishComputation(context, includeForce, includeEnergy, groups, valid);
        if (includeForce) {
            if (cl.getContextIndex() > 0) {
                int numAtoms = cl.getPaddedNumAtoms();
                void* dest = (cl.getUseDoublePrecision() ? (void*) &((mm_double4*) pinnedMemory)[(cl.getContextIndex()-1)*numAtoms] : (void*) &((mm_float4*) pinnedMemory)[(cl.getContextIndex()-1)*numAtoms]);
                cl.getQueue().enqueueReadBuffer(cl.getForce().getDeviceBuffer(), CL_TRUE, 0,
                        numAtoms*cl.getForce().getElementSize(), dest);
            }
            else
                cl.getQueue().finish();
        }
        completionTime = getCurrentTime();
        if (cl.getNonbondedUtilities().getUsePeriodic() && numTiles > cl.getNonbondedUtilities().getInteractingTiles().getSize()) {
            valid = false;
            cl.getNonbondedUtilities().updateNeighborListSize();
        }
    }
private:
    ContextImpl& context;
    OpenCLContext& cl;
    OpenCLCalcForcesAndEnergyKernel& kernel;
    bool includeForce, includeEnergy;
    int groups;
    double& energy;
    double& completionTime;
    void* pinnedMemory;
    bool& valid;
    int& numTiles;
};

OpenCLParallelCalcForcesAndEnergyKernel::OpenCLParallelCalcForcesAndEnergyKernel(string name, const Platform& platform, OpenCLPlatform::PlatformData& data) :
        CalcForcesAndEnergyKernel(name, platform), data(data), completionTimes(data.contexts.size()), contextNonbondedFractions(data.contexts.size()),
        tileCounts(data.contexts.size()), pinnedPositionBuffer(NULL), pinnedPositionMemory(NULL), pinnedForceBuffer(NULL), pinnedForceMemory(NULL) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new OpenCLCalcForcesAndEnergyKernel(name, platform, *data.contexts[i])));
}

OpenCLParallelCalcForcesAndEnergyKernel::~OpenCLParallelCalcForcesAndEnergyKernel() {
    if (pinnedPositionBuffer != NULL)
        delete pinnedPositionBuffer;
    if (pinnedForceBuffer != NULL)
        delete pinnedForceBuffer;
}

void OpenCLParallelCalcForcesAndEnergyKernel::initialize(const System& system) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system);
    for (int i = 0; i < contextNonbondedFractions.size(); i++) {
        double x0 = i/(double) contextNonbondedFractions.size();
        double x1 = (i+1)/(double) contextNonbondedFractions.size();
        contextNonbondedFractions[i] = x1*x1 - x0*x0;
    }
}

void OpenCLParallelCalcForcesAndEnergyKernel::beginComputation(ContextImpl& context, bool includeForce, bool includeEnergy, int groups) {
    OpenCLContext& cl0 = *data.contexts[0];
    int elementSize = (cl0.getUseDoublePrecision() ? sizeof(mm_double4) : sizeof(mm_float4));
    if (!contextForces.isInitialized()) {
        contextForces.initialize<mm_float4>(cl0, &cl0.getForceBuffers().getDeviceBuffer(),
                data.contexts.size()*cl0.getPaddedNumAtoms(), "contextForces");
        int bufferBytes = (data.contexts.size()-1)*cl0.getPaddedNumAtoms()*elementSize;
        pinnedPositionBuffer = new cl::Buffer(cl0.getContext(), CL_MEM_ALLOC_HOST_PTR, bufferBytes);
        pinnedPositionMemory = cl0.getQueue().enqueueMapBuffer(*pinnedPositionBuffer, CL_TRUE, CL_MAP_READ | CL_MAP_WRITE, 0, bufferBytes);
        pinnedForceBuffer = new cl::Buffer(cl0.getContext(), CL_MEM_ALLOC_HOST_PTR, bufferBytes);
        pinnedForceMemory = cl0.getQueue().enqueueMapBuffer(*pinnedForceBuffer, CL_TRUE, CL_MAP_READ | CL_MAP_WRITE, 0, bufferBytes);
    }

    // Copy coordinates over to each device and execute the kernel.
    
    cl0.getQueue().enqueueReadBuffer(cl0.getPosq().getDeviceBuffer(), CL_TRUE, 0, cl0.getPaddedNumAtoms()*elementSize, pinnedPositionMemory);
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        data.contextEnergy[i] = 0.0;
        OpenCLContext& cl = *data.contexts[i];
        ComputeContext::WorkThread& thread = cl.getWorkThread();
        thread.addTask(new BeginComputationTask(context, cl, getKernel(i), includeForce, includeEnergy, groups, pinnedPositionMemory, tileCounts[i]));
    }
}

double OpenCLParallelCalcForcesAndEnergyKernel::finishComputation(ContextImpl& context, bool includeForce, bool includeEnergy, int groups, bool& valid) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        OpenCLContext& cl = *data.contexts[i];
        ComputeContext::WorkThread& thread = cl.getWorkThread();
        thread.addTask(new FinishComputationTask(context, cl, getKernel(i), includeForce, includeEnergy, groups, data.contextEnergy[i], completionTimes[i], pinnedForceMemory, valid, tileCounts[i]));
    }
    data.syncContexts();
    double energy = 0.0;
    for (int i = 0; i < (int) data.contextEnergy.size(); i++)
        energy += data.contextEnergy[i];
    if (includeForce && valid) {
        // Sum the forces from all devices.
        
        OpenCLContext& cl = *data.contexts[0];
        int numAtoms = cl.getPaddedNumAtoms();
        int elementSize = (cl.getUseDoublePrecision() ? sizeof(mm_double4) : sizeof(mm_float4));
        cl.getQueue().enqueueWriteBuffer(contextForces.getDeviceBuffer(), CL_FALSE, numAtoms*elementSize,
                numAtoms*(data.contexts.size()-1)*elementSize, pinnedForceMemory);
        cl.reduceBuffer(contextForces, cl.getLongForceBuffer(), data.contexts.size());
        
        // Balance work between the contexts by transferring a little nonbonded work from the context that
        // finished last to the one that finished first.
        
        if (cl.getComputeForceCount() < 200 || cl.getComputeForceCount()%30 == 0) {
            int firstIndex = 0, lastIndex = 0;
            for (int i = 0; i < (int) completionTimes.size(); i++) {
                if (completionTimes[i] < completionTimes[firstIndex])
                    firstIndex = i;
                if (completionTimes[i] > completionTimes[lastIndex])
                    lastIndex = i;
            }
            double fractionToTransfer = min(0.001, contextNonbondedFractions[lastIndex]);
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

class OpenCLParallelCalcNonbondedForceKernel::Task : public OpenCLContext::WorkTask {
public:
    Task(ContextImpl& context, OpenCLCalcNonbondedForceKernel& kernel, bool includeForce,
            bool includeEnergy, bool includeDirect, bool includeReciprocal, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), includeDirect(includeDirect), includeReciprocal(includeReciprocal), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy, includeDirect, includeReciprocal);
    }
private:
    ContextImpl& context;
    OpenCLCalcNonbondedForceKernel& kernel;
    bool includeForce, includeEnergy, includeDirect, includeReciprocal;
    double& energy;
};

OpenCLParallelCalcNonbondedForceKernel::OpenCLParallelCalcNonbondedForceKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data, const System& system) :
        CalcNonbondedForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new OpenCLCalcNonbondedForceKernel(name, platform, *data.contexts[i], system)));
}

void OpenCLParallelCalcNonbondedForceKernel::initialize(const System& system, const NonbondedForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double OpenCLParallelCalcNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy, bool includeDirect, bool includeReciprocal) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        OpenCLContext& cl = *data.contexts[i];
        ComputeContext::WorkThread& thread = cl.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, includeDirect, includeReciprocal, data.contextEnergy[i]));
    }
    return 0.0;
}

void OpenCLParallelCalcNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const NonbondedForce& force, int firstParticle, int lastParticle, int firstException, int lastException) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force, firstParticle, lastParticle, firstException, lastException);
}

void OpenCLParallelCalcNonbondedForceKernel::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    dynamic_cast<const OpenCLCalcNonbondedForceKernel&>(kernels[0].getImpl()).getPMEParameters(alpha, nx, ny, nz);
}

void OpenCLParallelCalcNonbondedForceKernel::getLJPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    dynamic_cast<const OpenCLCalcNonbondedForceKernel&>(kernels[0].getImpl()).getLJPMEParameters(alpha, nx, ny, nz);
}
