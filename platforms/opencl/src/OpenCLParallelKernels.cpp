/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011 Stanford University and the Authors.           *
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

using namespace OpenMM;
using namespace std;

/**
 * Get the current clock time, measured in microseconds.
 */
#ifdef _MSC_VER
    #include <Windows.h>
    static long long getTime() {
        FILETIME ft;
        GetSystemTimeAsFileTime(&ft);	 // 100-nanoseconds since 1-1-1601
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

class OpenCLParallelCalcForcesAndEnergyKernel::BeginComputationTask : public OpenCLContext::WorkTask {
public:
    BeginComputationTask(ContextImpl& context, OpenCLContext& cl, OpenCLCalcForcesAndEnergyKernel& kernel,
            bool includeForce, bool includeEnergy) : context(context), cl(cl), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy) {
    }
    void execute() {
        // Copy coordinates over to this device and execute the kernel.

        if (cl.getContextIndex() > 0)
            cl.getPosq().upload(cl.getPlatformData().contexts[0]->getPosq().getHostBuffer(), false);
        kernel.beginComputation(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    OpenCLContext& cl;
    OpenCLCalcForcesAndEnergyKernel& kernel;
    bool includeForce, includeEnergy;
};

class OpenCLParallelCalcForcesAndEnergyKernel::FinishComputationTask : public OpenCLContext::WorkTask {
public:
    FinishComputationTask(ContextImpl& context, OpenCLContext& cl, OpenCLCalcForcesAndEnergyKernel& kernel,
            bool includeForce, bool includeEnergy, double& energy, long long& completionTime, OpenCLArray<mm_float4>& contextForces) :
            context(context), cl(cl), kernel(kernel), includeForce(includeForce), includeEnergy(includeEnergy), energy(energy),
            completionTime(completionTime), contextForces(contextForces) {
    }
    void execute() {
        // Execute the kernel, then download forces.
        
        energy += kernel.finishComputation(context, includeForce, includeEnergy);
        if (includeForce) {
            if (cl.getContextIndex() > 0)
                cl.getForce().download(&contextForces[cl.getContextIndex()*cl.getPaddedNumAtoms()]);
            else
                cl.getQueue().finish();
        }
        completionTime = getTime();
    }
private:
    ContextImpl& context;
    OpenCLContext& cl;
    OpenCLCalcForcesAndEnergyKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
    long long& completionTime;
    OpenCLArray<mm_float4>& contextForces;
};

OpenCLParallelCalcForcesAndEnergyKernel::OpenCLParallelCalcForcesAndEnergyKernel(string name, const Platform& platform, OpenCLPlatform::PlatformData& data) :
        CalcForcesAndEnergyKernel(name, platform), data(data), completionTimes(data.contexts.size()), contextTiles(data.contexts.size()), contextForces(NULL) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new OpenCLCalcForcesAndEnergyKernel(name, platform, *data.contexts[i])));
}

OpenCLParallelCalcForcesAndEnergyKernel::~OpenCLParallelCalcForcesAndEnergyKernel() {
    if (contextForces != NULL)
        delete contextForces;
}

void OpenCLParallelCalcForcesAndEnergyKernel::initialize(const System& system) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system);
}

void OpenCLParallelCalcForcesAndEnergyKernel::beginComputation(ContextImpl& context, bool includeForce, bool includeEnergy) {
    // Copy coordinates over to each device and execute the kernel.
    
    data.contexts[0]->getPosq().download();
    if (contextForces == NULL) {
        OpenCLContext& cl = *data.contexts[0];
        contextForces = new OpenCLArray<mm_float4>(cl, &cl.getForceBuffers().getDeviceBuffer(),
                data.contexts.size()*cl.getPaddedNumAtoms(), "contextForces", true);
    }
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        data.contextEnergy[i] = 0.0;
        OpenCLContext& cl = *data.contexts[i];
        OpenCLContext::WorkThread& thread = cl.getWorkThread();
        thread.addTask(new BeginComputationTask(context, cl, getKernel(i), includeForce, includeEnergy));
    }
}

double OpenCLParallelCalcForcesAndEnergyKernel::finishComputation(ContextImpl& context, bool includeForce, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        OpenCLContext& cl = *data.contexts[i];
        OpenCLContext::WorkThread& thread = cl.getWorkThread();
        thread.addTask(new FinishComputationTask(context, cl, getKernel(i), includeForce, includeEnergy, data.contextEnergy[i], completionTimes[i], *contextForces));
    }
    data.syncContexts();
    double energy = 0.0;
    for (int i = 0; i < (int) data.contextEnergy.size(); i++)
        energy += data.contextEnergy[i];
    if (includeForce) {
        // Sum the forces from all devices.
        
        OpenCLContext& cl = *data.contexts[0];
        int numAtoms = cl.getPaddedNumAtoms();
        cl.getQueue().enqueueWriteBuffer(contextForces->getDeviceBuffer(), CL_FALSE, numAtoms*sizeof(mm_float4),
                numAtoms*(data.contexts.size()-1)*sizeof(mm_float4), &(*contextForces)[numAtoms]);
        cl.reduceBuffer(*contextForces, data.contexts.size());
        
        // Balance work between the contexts by transferring a few nonbonded tiles from the context that
        // finished last to the one that finished first.
        
        int firstIndex = 0, lastIndex = 0;
        int totalTiles = 0;
        for (int i = 0; i < (int) completionTimes.size(); i++) {
            if (completionTimes[i] < completionTimes[firstIndex])
                firstIndex = i;
            if (completionTimes[i] > completionTimes[lastIndex])
                lastIndex = i;
            contextTiles[i] = data.contexts[i]->getNonbondedUtilities().getNumTiles();
            totalTiles += contextTiles[i];
        }
        int tilesToTransfer = totalTiles/1000;
        if (tilesToTransfer < 1)
            tilesToTransfer = 1;
        if (tilesToTransfer > contextTiles[lastIndex])
            tilesToTransfer = contextTiles[lastIndex];
        contextTiles[firstIndex] += tilesToTransfer;
        contextTiles[lastIndex] -= tilesToTransfer;
        int startIndex = 0;
        for (int i = 0; i < (int) contextTiles.size(); i++) {
            data.contexts[i]->getNonbondedUtilities().setTileRange(startIndex, contextTiles[i]);
            startIndex += contextTiles[i];
        }
    }
    return energy;
}

class OpenCLParallelCalcHarmonicBondForceKernel::Task : public OpenCLContext::WorkTask {
public:
    Task(ContextImpl& context, OpenCLCalcHarmonicBondForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    OpenCLCalcHarmonicBondForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

OpenCLParallelCalcHarmonicBondForceKernel::OpenCLParallelCalcHarmonicBondForceKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data, System& system) :
        CalcHarmonicBondForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new OpenCLCalcHarmonicBondForceKernel(name, platform, *data.contexts[i], system)));
}

void OpenCLParallelCalcHarmonicBondForceKernel::initialize(const System& system, const HarmonicBondForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double OpenCLParallelCalcHarmonicBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        OpenCLContext& cl = *data.contexts[i];
        OpenCLContext::WorkThread& thread = cl.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

class OpenCLParallelCalcCustomBondForceKernel::Task : public OpenCLContext::WorkTask {
public:
    Task(ContextImpl& context, OpenCLCalcCustomBondForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    OpenCLCalcCustomBondForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

OpenCLParallelCalcCustomBondForceKernel::OpenCLParallelCalcCustomBondForceKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data, System& system) :
        CalcCustomBondForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new OpenCLCalcCustomBondForceKernel(name, platform, *data.contexts[i], system)));
}

void OpenCLParallelCalcCustomBondForceKernel::initialize(const System& system, const CustomBondForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double OpenCLParallelCalcCustomBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        OpenCLContext& cl = *data.contexts[i];
        OpenCLContext::WorkThread& thread = cl.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

class OpenCLParallelCalcHarmonicAngleForceKernel::Task : public OpenCLContext::WorkTask {
public:
    Task(ContextImpl& context, OpenCLCalcHarmonicAngleForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    OpenCLCalcHarmonicAngleForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

OpenCLParallelCalcHarmonicAngleForceKernel::OpenCLParallelCalcHarmonicAngleForceKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data, System& system) :
        CalcHarmonicAngleForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new OpenCLCalcHarmonicAngleForceKernel(name, platform, *data.contexts[i], system)));
}

void OpenCLParallelCalcHarmonicAngleForceKernel::initialize(const System& system, const HarmonicAngleForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double OpenCLParallelCalcHarmonicAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        OpenCLContext& cl = *data.contexts[i];
        OpenCLContext::WorkThread& thread = cl.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

class OpenCLParallelCalcCustomAngleForceKernel::Task : public OpenCLContext::WorkTask {
public:
    Task(ContextImpl& context, OpenCLCalcCustomAngleForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    OpenCLCalcCustomAngleForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

OpenCLParallelCalcCustomAngleForceKernel::OpenCLParallelCalcCustomAngleForceKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data, System& system) :
        CalcCustomAngleForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new OpenCLCalcCustomAngleForceKernel(name, platform, *data.contexts[i], system)));
}

void OpenCLParallelCalcCustomAngleForceKernel::initialize(const System& system, const CustomAngleForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double OpenCLParallelCalcCustomAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        OpenCLContext& cl = *data.contexts[i];
        OpenCLContext::WorkThread& thread = cl.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

class OpenCLParallelCalcPeriodicTorsionForceKernel::Task : public OpenCLContext::WorkTask {
public:
    Task(ContextImpl& context, OpenCLCalcPeriodicTorsionForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    OpenCLCalcPeriodicTorsionForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

OpenCLParallelCalcPeriodicTorsionForceKernel::OpenCLParallelCalcPeriodicTorsionForceKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data, System& system) :
        CalcPeriodicTorsionForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new OpenCLCalcPeriodicTorsionForceKernel(name, platform, *data.contexts[i], system)));
}

void OpenCLParallelCalcPeriodicTorsionForceKernel::initialize(const System& system, const PeriodicTorsionForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double OpenCLParallelCalcPeriodicTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        OpenCLContext& cl = *data.contexts[i];
        OpenCLContext::WorkThread& thread = cl.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

class OpenCLParallelCalcRBTorsionForceKernel::Task : public OpenCLContext::WorkTask {
public:
    Task(ContextImpl& context, OpenCLCalcRBTorsionForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    OpenCLCalcRBTorsionForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

OpenCLParallelCalcRBTorsionForceKernel::OpenCLParallelCalcRBTorsionForceKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data, System& system) :
        CalcRBTorsionForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new OpenCLCalcRBTorsionForceKernel(name, platform, *data.contexts[i], system)));
}

void OpenCLParallelCalcRBTorsionForceKernel::initialize(const System& system, const RBTorsionForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double OpenCLParallelCalcRBTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        OpenCLContext& cl = *data.contexts[i];
        OpenCLContext::WorkThread& thread = cl.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

class OpenCLParallelCalcCMAPTorsionForceKernel::Task : public OpenCLContext::WorkTask {
public:
    Task(ContextImpl& context, OpenCLCalcCMAPTorsionForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    OpenCLCalcCMAPTorsionForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

OpenCLParallelCalcCMAPTorsionForceKernel::OpenCLParallelCalcCMAPTorsionForceKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data, System& system) :
        CalcCMAPTorsionForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new OpenCLCalcCMAPTorsionForceKernel(name, platform, *data.contexts[i], system)));
}

void OpenCLParallelCalcCMAPTorsionForceKernel::initialize(const System& system, const CMAPTorsionForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double OpenCLParallelCalcCMAPTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        OpenCLContext& cl = *data.contexts[i];
        OpenCLContext::WorkThread& thread = cl.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

class OpenCLParallelCalcCustomTorsionForceKernel::Task : public OpenCLContext::WorkTask {
public:
    Task(ContextImpl& context, OpenCLCalcCustomTorsionForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    OpenCLCalcCustomTorsionForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

OpenCLParallelCalcCustomTorsionForceKernel::OpenCLParallelCalcCustomTorsionForceKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data, System& system) :
        CalcCustomTorsionForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new OpenCLCalcCustomTorsionForceKernel(name, platform, *data.contexts[i], system)));
}

void OpenCLParallelCalcCustomTorsionForceKernel::initialize(const System& system, const CustomTorsionForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double OpenCLParallelCalcCustomTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        OpenCLContext& cl = *data.contexts[i];
        OpenCLContext::WorkThread& thread = cl.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

class OpenCLParallelCalcNonbondedForceKernel::Task : public OpenCLContext::WorkTask {
public:
    Task(ContextImpl& context, OpenCLCalcNonbondedForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    OpenCLCalcNonbondedForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

OpenCLParallelCalcNonbondedForceKernel::OpenCLParallelCalcNonbondedForceKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data, System& system) :
        CalcNonbondedForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new OpenCLCalcNonbondedForceKernel(name, platform, *data.contexts[i], system)));
}

void OpenCLParallelCalcNonbondedForceKernel::initialize(const System& system, const NonbondedForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double OpenCLParallelCalcNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        OpenCLContext& cl = *data.contexts[i];
        OpenCLContext::WorkThread& thread = cl.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

class OpenCLParallelCalcCustomNonbondedForceKernel::Task : public OpenCLContext::WorkTask {
public:
    Task(ContextImpl& context, OpenCLCalcCustomNonbondedForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    OpenCLCalcCustomNonbondedForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

OpenCLParallelCalcCustomNonbondedForceKernel::OpenCLParallelCalcCustomNonbondedForceKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data, System& system) :
        CalcCustomNonbondedForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new OpenCLCalcCustomNonbondedForceKernel(name, platform, *data.contexts[i], system)));
}

void OpenCLParallelCalcCustomNonbondedForceKernel::initialize(const System& system, const CustomNonbondedForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double OpenCLParallelCalcCustomNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        OpenCLContext& cl = *data.contexts[i];
        OpenCLContext::WorkThread& thread = cl.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

class OpenCLParallelCalcCustomExternalForceKernel::Task : public OpenCLContext::WorkTask {
public:
    Task(ContextImpl& context, OpenCLCalcCustomExternalForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    OpenCLCalcCustomExternalForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

OpenCLParallelCalcCustomExternalForceKernel::OpenCLParallelCalcCustomExternalForceKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data, System& system) :
        CalcCustomExternalForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new OpenCLCalcCustomExternalForceKernel(name, platform, *data.contexts[i], system)));
}

void OpenCLParallelCalcCustomExternalForceKernel::initialize(const System& system, const CustomExternalForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double OpenCLParallelCalcCustomExternalForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        OpenCLContext& cl = *data.contexts[i];
        OpenCLContext::WorkThread& thread = cl.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

class OpenCLParallelCalcCustomHbondForceKernel::Task : public OpenCLContext::WorkTask {
public:
    Task(ContextImpl& context, OpenCLCalcCustomHbondForceKernel& kernel, bool includeForce,
            bool includeEnergy, double& energy) : context(context), kernel(kernel),
            includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    OpenCLCalcCustomHbondForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

OpenCLParallelCalcCustomHbondForceKernel::OpenCLParallelCalcCustomHbondForceKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data, System& system) :
        CalcCustomHbondForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new OpenCLCalcCustomHbondForceKernel(name, platform, *data.contexts[i], system)));
}

void OpenCLParallelCalcCustomHbondForceKernel::initialize(const System& system, const CustomHbondForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double OpenCLParallelCalcCustomHbondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        OpenCLContext& cl = *data.contexts[i];
        OpenCLContext::WorkThread& thread = cl.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}
