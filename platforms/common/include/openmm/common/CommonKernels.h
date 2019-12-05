#ifndef OPENMM_COMMONKERNELS_H_
#define OPENMM_COMMONKERNELS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2019 Stanford University and the Authors.      *
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

#include "openmm/common/ComputeArray.h"
#include "openmm/common/ComputeContext.h"
#include "openmm/common/ComputeParameterSet.h"
#include "openmm/Platform.h"
#include "openmm/kernels.h"

namespace OpenMM {


/**
 * This kernel is invoked by HarmonicBondForce to calculate the forces acting on the system and the energy of the system.
 */
class CommonCalcHarmonicBondForceKernel : public CalcHarmonicBondForceKernel {
public:
    CommonCalcHarmonicBondForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) : CalcHarmonicBondForceKernel(name, platform),
            hasInitializedKernel(false), cc(cc), system(system) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the HarmonicBondForce this kernel will be used for
     */
    void initialize(const System& system, const HarmonicBondForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the HarmonicBondForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const HarmonicBondForce& force);
private:
    class ForceInfo;
    int numBonds;
    bool hasInitializedKernel;
    ComputeContext& cc;
    ForceInfo* info;
    const System& system;
    ComputeArray params;
};

/**
 * This kernel is invoked by CustomBondForce to calculate the forces acting on the system and the energy of the system.
 */
class CommonCalcCustomBondForceKernel : public CalcCustomBondForceKernel {
public:
    CommonCalcCustomBondForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) : CalcCustomBondForceKernel(name, platform),
            hasInitializedKernel(false), cc(cc), system(system), params(NULL) {
    }
    ~CommonCalcCustomBondForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomBondForce this kernel will be used for
     */
    void initialize(const System& system, const CustomBondForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomBondForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const CustomBondForce& force);
private:
    class ForceInfo;
    int numBonds;
    bool hasInitializedKernel;
    ComputeContext& cc;
    ForceInfo* info;
    const System& system;
    ComputeParameterSet* params;
    ComputeArray globals;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
};

/**
 * This kernel is invoked by HarmonicAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class CommonCalcHarmonicAngleForceKernel : public CalcHarmonicAngleForceKernel {
public:
    CommonCalcHarmonicAngleForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) : CalcHarmonicAngleForceKernel(name, platform),
            hasInitializedKernel(false), cc(cc), system(system) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the HarmonicAngleForce this kernel will be used for
     */
    void initialize(const System& system, const HarmonicAngleForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the HarmonicAngleForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const HarmonicAngleForce& force);
private:
    class ForceInfo;
    int numAngles;
    bool hasInitializedKernel;
    ComputeContext& cc;
    ForceInfo* info;
    const System& system;
    ComputeArray params;
};

/**
 * This kernel is invoked by CustomAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class CommonCalcCustomAngleForceKernel : public CalcCustomAngleForceKernel {
public:
    CommonCalcCustomAngleForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) : CalcCustomAngleForceKernel(name, platform),
            hasInitializedKernel(false), cc(cc), system(system), params(NULL) {
    }
    ~CommonCalcCustomAngleForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomAngleForce this kernel will be used for
     */
    void initialize(const System& system, const CustomAngleForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomAngleForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const CustomAngleForce& force);
private:
    class ForceInfo;
    int numAngles;
    bool hasInitializedKernel;
    ComputeContext& cc;
    ForceInfo* info;
    const System& system;
    ComputeParameterSet* params;
    ComputeArray globals;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
};

/**
 * This kernel is invoked by PeriodicTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CommonCalcPeriodicTorsionForceKernel : public CalcPeriodicTorsionForceKernel {
public:
    CommonCalcPeriodicTorsionForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) : CalcPeriodicTorsionForceKernel(name, platform),
            hasInitializedKernel(false), cc(cc), system(system) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the PeriodicTorsionForce this kernel will be used for
     */
    void initialize(const System& system, const PeriodicTorsionForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the PeriodicTorsionForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const PeriodicTorsionForce& force);
private:
    class ForceInfo;
    int numTorsions;
    bool hasInitializedKernel;
    ComputeContext& cc;
    ForceInfo* info;
    const System& system;
    ComputeArray params;
};

/**
 * This kernel is invoked by RBTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CommonCalcRBTorsionForceKernel : public CalcRBTorsionForceKernel {
public:
    CommonCalcRBTorsionForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) : CalcRBTorsionForceKernel(name, platform),
            hasInitializedKernel(false), cc(cc), system(system) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the RBTorsionForce this kernel will be used for
     */
    void initialize(const System& system, const RBTorsionForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the RBTorsionForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const RBTorsionForce& force);
private:
    class ForceInfo;
    int numTorsions;
    bool hasInitializedKernel;
    ComputeContext& cc;
    ForceInfo* info;
    const System& system;
    ComputeArray params1;
    ComputeArray params2;
};

/**
 * This kernel is invoked by CustomTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CommonCalcCustomTorsionForceKernel : public CalcCustomTorsionForceKernel {
public:
    CommonCalcCustomTorsionForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) : CalcCustomTorsionForceKernel(name, platform),
            hasInitializedKernel(false), cc(cc), system(system), params(NULL) {
    }
    ~CommonCalcCustomTorsionForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomTorsionForce this kernel will be used for
     */
    void initialize(const System& system, const CustomTorsionForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomTorsionForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const CustomTorsionForce& force);
private:
    class ForceInfo;
    int numTorsions;
    bool hasInitializedKernel;
    ComputeContext& cc;
    ForceInfo* info;
    const System& system;
    ComputeParameterSet* params;
    ComputeArray globals;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
};

/**
 * This kernel is invoked by CMAPTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CommonCalcCMAPTorsionForceKernel : public CalcCMAPTorsionForceKernel {
public:
    CommonCalcCMAPTorsionForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) : CalcCMAPTorsionForceKernel(name, platform),
            hasInitializedKernel(false), cc(cc), system(system) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CMAPTorsionForce this kernel will be used for
     */
    void initialize(const System& system, const CMAPTorsionForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CMAPTorsionForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const CMAPTorsionForce& force);
private:
    class ForceInfo;
    int numTorsions;
    bool hasInitializedKernel;
    ComputeContext& cc;
    ForceInfo* info;
    const System& system;
    std::vector<mm_int2> mapPositionsVec;
    ComputeArray coefficients;
    ComputeArray mapPositions;
    ComputeArray torsionMaps;
};

/**
 * This kernel is invoked by CustomExternalForce to calculate the forces acting on the system and the energy of the system.
 */
class CommonCalcCustomExternalForceKernel : public CalcCustomExternalForceKernel {
public:
    CommonCalcCustomExternalForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) : CalcCustomExternalForceKernel(name, platform),
            hasInitializedKernel(false), cc(cc), system(system), params(NULL) {
    }
    ~CommonCalcCustomExternalForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomExternalForce this kernel will be used for
     */
    void initialize(const System& system, const CustomExternalForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomExternalForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const CustomExternalForce& force);
private:
    class ForceInfo;
    int numParticles;
    bool hasInitializedKernel;
    ComputeContext& cc;
    ForceInfo* info;
    const System& system;
    ComputeParameterSet* params;
    ComputeArray globals;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
};

/**
 * This kernel is invoked by CustomCompoundBondForce to calculate the forces acting on the system.
 */
class CommonCalcCustomCompoundBondForceKernel : public CalcCustomCompoundBondForceKernel {
public:
    CommonCalcCustomCompoundBondForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) : CalcCustomCompoundBondForceKernel(name, platform),
            cc(cc), params(NULL), system(system) {
    }
    ~CommonCalcCustomCompoundBondForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomCompoundBondForce this kernel will be used for
     */
    void initialize(const System& system, const CustomCompoundBondForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomCompoundBondForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const CustomCompoundBondForce& force);

private:
    class ForceInfo;
    int numBonds;
    ComputeContext& cc;
    ForceInfo* info;
    ComputeParameterSet* params;
    ComputeArray globals;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
    std::vector<ComputeArray> tabulatedFunctions;
    const System& system;
};

/**
 * This kernel is invoked by CustomCentroidBondForce to calculate the forces acting on the system.
 */
class CommonCalcCustomCentroidBondForceKernel : public CalcCustomCentroidBondForceKernel {
public:
    CommonCalcCustomCentroidBondForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) : CalcCustomCentroidBondForceKernel(name, platform),
            cc(cc), params(NULL), system(system) {
    }
    ~CommonCalcCustomCentroidBondForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomCentroidBondForce this kernel will be used for
     */
    void initialize(const System& system, const CustomCentroidBondForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomCentroidBondForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const CustomCentroidBondForce& force);

private:
    class ForceInfo;
    int numGroups, numBonds;
    bool needEnergyParamDerivs;
    ComputeContext& cc;
    ForceInfo* info;
    ComputeParameterSet* params;
    ComputeArray globals, groupParticles, groupWeights, groupOffsets;
    ComputeArray groupForces, bondGroups, centerPositions;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
    std::vector<ComputeArray> tabulatedFunctions;
    std::vector<void*> groupForcesArgs;
    ComputeKernel computeCentersKernel, groupForcesKernel, applyForcesKernel;
    const System& system;
};

/**
 * This kernel is invoked by CustomManyParticleForce to calculate the forces acting on the system.
 */
class CommonCalcCustomManyParticleForceKernel : public CalcCustomManyParticleForceKernel {
public:
    CommonCalcCustomManyParticleForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) : CalcCustomManyParticleForceKernel(name, platform),
            hasInitializedKernel(false), cc(cc), params(NULL), system(system) {
    }
    ~CommonCalcCustomManyParticleForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomManyParticleForce this kernel will be used for
     */
    void initialize(const System& system, const CustomManyParticleForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomManyParticleForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const CustomManyParticleForce& force);

private:
    class ForceInfo;
    ComputeContext& cc;
    ForceInfo* info;
    bool hasInitializedKernel;
    NonbondedMethod nonbondedMethod;
    int maxNeighborPairs, forceWorkgroupSize, findNeighborsWorkgroupSize;
    ComputeParameterSet* params;
    ComputeArray globals, particleTypes,  orderIndex, particleOrder;
    ComputeArray exclusions, exclusionStartIndex, blockCenter, blockBoundingBox;
    ComputeArray neighborPairs, numNeighborPairs, neighborStartIndex, numNeighborsForAtom, neighbors;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
    std::vector<ComputeArray> tabulatedFunctions;
    const System& system;
    ComputeKernel forceKernel, blockBoundsKernel, neighborsKernel, startIndicesKernel, copyPairsKernel;
    ComputeEvent event;
};

/**
 * This kernel is invoked to remove center of mass motion from the system.
 */
class CommonRemoveCMMotionKernel : public RemoveCMMotionKernel {
public:
    CommonRemoveCMMotionKernel(std::string name, const Platform& platform, ComputeContext& cc) : RemoveCMMotionKernel(name, platform), cc(cc) {
    }
    /**
     * Initialize the kernel, setting up the particle masses.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CMMotionRemover this kernel will be used for
     */
    void initialize(const System& system, const CMMotionRemover& force);
    /**
     * Execute the kernel.
     *
     * @param context    the context in which to execute this kernel
     */
    void execute(ContextImpl& context);
private:
    ComputeContext& cc;
    int frequency;
    ComputeArray cmMomentum;
    ComputeKernel kernel1, kernel2;
};

} // namespace OpenMM

#endif /*OPENMM_COMMONKERNELS_H_*/
