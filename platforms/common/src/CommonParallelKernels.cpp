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

#include "openmm/common/CommonParallelKernels.h"

using namespace OpenMM;
using namespace std;

class CommonParallelCalcHarmonicBondForceKernel::Task : public ComputeContext::WorkTask {
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

CommonParallelCalcHarmonicBondForceKernel::CommonParallelCalcHarmonicBondForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) :
        CalcHarmonicBondForceKernel(name, platform), cc(cc) {
    for (ComputeContext* context : cc.getAllContexts())
        kernels.push_back(Kernel(new CommonCalcHarmonicBondForceKernel(name, platform, *context, system)));
}

void CommonParallelCalcHarmonicBondForceKernel::initialize(const System& system, const HarmonicBondForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CommonParallelCalcHarmonicBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < cc.getNumContexts(); i++) {
        ComputeContext::WorkThread& thread = cc.getAllContexts()[i]->getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, cc.getEnergyWorkspace()));
    }
    return 0.0;
}

void CommonParallelCalcHarmonicBondForceKernel::copyParametersToContext(ContextImpl& context, const HarmonicBondForce& force, int firstBond, int lastBond) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force, firstBond, lastBond);
}

class CommonParallelCalcCustomBondForceKernel::Task : public ComputeContext::WorkTask {
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

CommonParallelCalcCustomBondForceKernel::CommonParallelCalcCustomBondForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) :
        CalcCustomBondForceKernel(name, platform), cc(cc) {
    for (ComputeContext* context : cc.getAllContexts())
        kernels.push_back(Kernel(new CommonCalcCustomBondForceKernel(name, platform, *context, system)));
}

void CommonParallelCalcCustomBondForceKernel::initialize(const System& system, const CustomBondForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CommonParallelCalcCustomBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < cc.getNumContexts(); i++) {
        ComputeContext::WorkThread& thread = cc.getAllContexts()[i]->getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, cc.getEnergyWorkspace()));
    }
    return 0.0;
}

void CommonParallelCalcCustomBondForceKernel::copyParametersToContext(ContextImpl& context, const CustomBondForce& force, int firstBond, int lastBond) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force, firstBond, lastBond);
}

class CommonParallelCalcHarmonicAngleForceKernel::Task : public ComputeContext::WorkTask {
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

CommonParallelCalcHarmonicAngleForceKernel::CommonParallelCalcHarmonicAngleForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) :
        CalcHarmonicAngleForceKernel(name, platform), cc(cc) {
    for (ComputeContext* context : cc.getAllContexts())
        kernels.push_back(Kernel(new CommonCalcHarmonicAngleForceKernel(name, platform, *context, system)));
}

void CommonParallelCalcHarmonicAngleForceKernel::initialize(const System& system, const HarmonicAngleForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CommonParallelCalcHarmonicAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < cc.getNumContexts(); i++) {
        ComputeContext::WorkThread& thread = cc.getAllContexts()[i]->getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, cc.getEnergyWorkspace()));
    }
    return 0.0;
}

void CommonParallelCalcHarmonicAngleForceKernel::copyParametersToContext(ContextImpl& context, const HarmonicAngleForce& force, int firstAngle, int lastAngle) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force, firstAngle, lastAngle);
}

class CommonParallelCalcCustomAngleForceKernel::Task : public ComputeContext::WorkTask {
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

CommonParallelCalcCustomAngleForceKernel::CommonParallelCalcCustomAngleForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) :
        CalcCustomAngleForceKernel(name, platform), cc(cc) {
    for (ComputeContext* context : cc.getAllContexts())
        kernels.push_back(Kernel(new CommonCalcCustomAngleForceKernel(name, platform, *context, system)));
}

void CommonParallelCalcCustomAngleForceKernel::initialize(const System& system, const CustomAngleForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CommonParallelCalcCustomAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < cc.getNumContexts(); i++) {
        ComputeContext::WorkThread& thread = cc.getAllContexts()[i]->getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, cc.getEnergyWorkspace()));
    }
    return 0.0;
}

void CommonParallelCalcCustomAngleForceKernel::copyParametersToContext(ContextImpl& context, const CustomAngleForce& force, int firstAngle, int lastAngle) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force, firstAngle, lastAngle);
}

class CommonParallelCalcPeriodicTorsionForceKernel::Task : public ComputeContext::WorkTask {
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

CommonParallelCalcPeriodicTorsionForceKernel::CommonParallelCalcPeriodicTorsionForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) :
        CalcPeriodicTorsionForceKernel(name, platform), cc(cc) {
    for (ComputeContext* context : cc.getAllContexts())
        kernels.push_back(Kernel(new CommonCalcPeriodicTorsionForceKernel(name, platform, *context, system)));
}

void CommonParallelCalcPeriodicTorsionForceKernel::initialize(const System& system, const PeriodicTorsionForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CommonParallelCalcPeriodicTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < cc.getNumContexts(); i++) {
        ComputeContext::WorkThread& thread = cc.getAllContexts()[i]->getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, cc.getEnergyWorkspace()));
    }
    return 0.0;
}

void CommonParallelCalcPeriodicTorsionForceKernel::copyParametersToContext(ContextImpl& context, const PeriodicTorsionForce& force, int firstTorsion, int lastTorsion) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force, firstTorsion, lastTorsion);
}

class CommonParallelCalcRBTorsionForceKernel::Task : public ComputeContext::WorkTask {
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

CommonParallelCalcRBTorsionForceKernel::CommonParallelCalcRBTorsionForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) :
        CalcRBTorsionForceKernel(name, platform), cc(cc) {
    for (ComputeContext* context : cc.getAllContexts())
        kernels.push_back(Kernel(new CommonCalcRBTorsionForceKernel(name, platform, *context, system)));
}

void CommonParallelCalcRBTorsionForceKernel::initialize(const System& system, const RBTorsionForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CommonParallelCalcRBTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < cc.getNumContexts(); i++) {
        ComputeContext::WorkThread& thread = cc.getAllContexts()[i]->getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, cc.getEnergyWorkspace()));
    }
    return 0.0;
}

void CommonParallelCalcRBTorsionForceKernel::copyParametersToContext(ContextImpl& context, const RBTorsionForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force);
}

class CommonParallelCalcCMAPTorsionForceKernel::Task : public ComputeContext::WorkTask {
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

CommonParallelCalcCMAPTorsionForceKernel::CommonParallelCalcCMAPTorsionForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) :
        CalcCMAPTorsionForceKernel(name, platform), cc(cc) {
    for (ComputeContext* context : cc.getAllContexts())
        kernels.push_back(Kernel(new CommonCalcCMAPTorsionForceKernel(name, platform, *context, system)));
}

void CommonParallelCalcCMAPTorsionForceKernel::initialize(const System& system, const CMAPTorsionForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CommonParallelCalcCMAPTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < cc.getNumContexts(); i++) {
        ComputeContext::WorkThread& thread = cc.getAllContexts()[i]->getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, cc.getEnergyWorkspace()));
    }
    return 0.0;
}

void CommonParallelCalcCMAPTorsionForceKernel::copyParametersToContext(ContextImpl& context, const CMAPTorsionForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force);
}

class CommonParallelCalcCustomTorsionForceKernel::Task : public ComputeContext::WorkTask {
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

CommonParallelCalcCustomTorsionForceKernel::CommonParallelCalcCustomTorsionForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) :
        CalcCustomTorsionForceKernel(name, platform), cc(cc) {
    for (ComputeContext* context : cc.getAllContexts())
        kernels.push_back(Kernel(new CommonCalcCustomTorsionForceKernel(name, platform, *context, system)));
}

void CommonParallelCalcCustomTorsionForceKernel::initialize(const System& system, const CustomTorsionForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CommonParallelCalcCustomTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < cc.getNumContexts(); i++) {
        ComputeContext::WorkThread& thread = cc.getAllContexts()[i]->getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, cc.getEnergyWorkspace()));
    }
    return 0.0;
}

void CommonParallelCalcCustomTorsionForceKernel::copyParametersToContext(ContextImpl& context, const CustomTorsionForce& force, int firstTorsion, int lastTorsion) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force, firstTorsion, lastTorsion);
}

class CommonParallelCalcCustomNonbondedForceKernel::Task : public ComputeContext::WorkTask {
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

CommonParallelCalcCustomNonbondedForceKernel::CommonParallelCalcCustomNonbondedForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) :
        CalcCustomNonbondedForceKernel(name, platform), cc(cc) {
    for (ComputeContext* context : cc.getAllContexts())
        kernels.push_back(Kernel(new CommonCalcCustomNonbondedForceKernel(name, platform, *context, system)));
}

void CommonParallelCalcCustomNonbondedForceKernel::initialize(const System& system, const CustomNonbondedForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CommonParallelCalcCustomNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < cc.getNumContexts(); i++) {
        ComputeContext::WorkThread& thread = cc.getAllContexts()[i]->getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, cc.getEnergyWorkspace()));
    }
    return 0.0;
}

void CommonParallelCalcCustomNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const CustomNonbondedForce& force, int firstParticle, int lastParticle) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force, firstParticle, lastParticle);
}

class CommonParallelCalcCustomExternalForceKernel::Task : public ComputeContext::WorkTask {
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

CommonParallelCalcCustomExternalForceKernel::CommonParallelCalcCustomExternalForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) :
        CalcCustomExternalForceKernel(name, platform), cc(cc) {
    for (ComputeContext* context : cc.getAllContexts())
        kernels.push_back(Kernel(new CommonCalcCustomExternalForceKernel(name, platform, *context, system)));
}

void CommonParallelCalcCustomExternalForceKernel::initialize(const System& system, const CustomExternalForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CommonParallelCalcCustomExternalForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < cc.getNumContexts(); i++) {
        ComputeContext::WorkThread& thread = cc.getAllContexts()[i]->getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, cc.getEnergyWorkspace()));
    }
    return 0.0;
}

void CommonParallelCalcCustomExternalForceKernel::copyParametersToContext(ContextImpl& context, const CustomExternalForce& force, int firstParticle, int lastParticle) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force, firstParticle, lastParticle);
}

class CommonParallelCalcCustomHbondForceKernel::Task : public ComputeContext::WorkTask {
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

CommonParallelCalcCustomHbondForceKernel::CommonParallelCalcCustomHbondForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) :
        CalcCustomHbondForceKernel(name, platform), cc(cc) {
    for (ComputeContext* context : cc.getAllContexts())
        kernels.push_back(Kernel(new CommonCalcCustomHbondForceKernel(name, platform, *context, system)));
}

void CommonParallelCalcCustomHbondForceKernel::initialize(const System& system, const CustomHbondForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CommonParallelCalcCustomHbondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < cc.getNumContexts(); i++) {
        ComputeContext::WorkThread& thread = cc.getAllContexts()[i]->getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, cc.getEnergyWorkspace()));
    }
    return 0.0;
}

void CommonParallelCalcCustomHbondForceKernel::copyParametersToContext(ContextImpl& context, const CustomHbondForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force);
}

class CommonParallelCalcCustomCompoundBondForceKernel::Task : public ComputeContext::WorkTask {
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

CommonParallelCalcCustomCompoundBondForceKernel::CommonParallelCalcCustomCompoundBondForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) :
        CalcCustomCompoundBondForceKernel(name, platform), cc(cc) {
    for (ComputeContext* context : cc.getAllContexts())
        kernels.push_back(Kernel(new CommonCalcCustomCompoundBondForceKernel(name, platform, *context, system)));
}

void CommonParallelCalcCustomCompoundBondForceKernel::initialize(const System& system, const CustomCompoundBondForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CommonParallelCalcCustomCompoundBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < cc.getNumContexts(); i++) {
        ComputeContext::WorkThread& thread = cc.getAllContexts()[i]->getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, cc.getEnergyWorkspace()));
    }
    return 0.0;
}

void CommonParallelCalcCustomCompoundBondForceKernel::copyParametersToContext(ContextImpl& context, const CustomCompoundBondForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force);
}
