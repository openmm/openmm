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
 * Portions copyright (c) 2008-2025 Stanford University and the Authors.      *
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
#include "openmm/internal/CompiledExpressionSet.h"
#include "openmm/internal/ContextImpl.h"
#include "lepton/CompiledExpression.h"
#include "lepton/CustomFunction.h"
#include "lepton/ExpressionProgram.h"

namespace OpenMM {


/**
 * This kernel provides methods for setting and retrieving various state data: time, positions,
 * velocities, and forces.
 */
class CommonUpdateStateDataKernel : public UpdateStateDataKernel {
public:
    CommonUpdateStateDataKernel(std::string name, const Platform& platform, ComputeContext& cc) : UpdateStateDataKernel(name, platform), cc(cc) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     */
    void initialize(const System& system);
    /**
     * Get the current time (in picoseconds).
     *
     * @param context    the context in which to execute this kernel
     */
    double getTime(const ContextImpl& context) const;
    /**
     * Set the current time (in picoseconds).
     *
     * @param context    the context in which to execute this kernel
     */
    void setTime(ContextImpl& context, double time);
    /**
     * Get the current step count
     *
     * @param context    the context in which to execute this kernel
     */
    long long getStepCount(const ContextImpl& context) const;
    /**
     * Set the current step count
     *
     * @param context    the context in which to execute this kernel
     */
    void setStepCount(const ContextImpl& context, long long count);
    /**
     * Get the positions of all particles.
     *
     * @param positions  on exit, this contains the particle positions
     */
    void getPositions(ContextImpl& context, std::vector<Vec3>& positions);
    /**
     * Set the positions of all particles.
     *
     * @param positions  a vector containg the particle positions
     */
    void setPositions(ContextImpl& context, const std::vector<Vec3>& positions);
    /**
     * Get the velocities of all particles.
     *
     * @param velocities  on exit, this contains the particle velocities
     */
    void getVelocities(ContextImpl& context, std::vector<Vec3>& velocities);
    /**
     * Set the velocities of all particles.
     *
     * @param velocities  a vector containg the particle velocities
     */
    void setVelocities(ContextImpl& context, const std::vector<Vec3>& velocities);
    /**
     * Compute velocities, shifted in time to account for a leapfrog integrator.  The shift
     * is based on the most recently computed forces.
     * 
     * @param context     the context in which to execute this kernel
     * @param timeShift   the amount by which to shift the velocities in time
     * @param velocities  the shifted velocities are returned in this
     */
    void computeShiftedVelocities(ContextImpl& context, double timeShift, std::vector<Vec3>& velocities);
    /**
     * Get the current forces on all particles.
     *
     * @param forces  on exit, this contains the forces
     */
    void getForces(ContextImpl& context, std::vector<Vec3>& forces);
    /**
     * Get the current derivatives of the energy with respect to context parameters.
     *
     * @param derivs  on exit, this contains the derivatives
     */
    void getEnergyParameterDerivatives(ContextImpl& context, std::map<std::string, double>& derivs);
    /**
     * Get the current periodic box vectors.
     *
     * @param a      on exit, this contains the vector defining the first edge of the periodic box
     * @param b      on exit, this contains the vector defining the second edge of the periodic box
     * @param c      on exit, this contains the vector defining the third edge of the periodic box
     */
    void getPeriodicBoxVectors(ContextImpl& context, Vec3& a, Vec3& b, Vec3& c) const;
    /**
     * Set the current periodic box vectors.
     *
     * @param a      the vector defining the first edge of the periodic box
     * @param b      the vector defining the second edge of the periodic box
     * @param c      the vector defining the third edge of the periodic box
     */
    void setPeriodicBoxVectors(ContextImpl& context, const Vec3& a, const Vec3& b, const Vec3& c);
    /**
     * Create a checkpoint recording the current state of the Context.
     * 
     * @param stream    an output stream the checkpoint data should be written to
     */
    void createCheckpoint(ContextImpl& context, std::ostream& stream);
    /**
     * Load a checkpoint that was written by createCheckpoint().
     * 
     * @param stream    an input stream the checkpoint data should be read from
     */
    void loadCheckpoint(ContextImpl& context, std::istream& stream);
private:
    ComputeContext& cc;
};

/**
 * This kernel modifies the positions of particles to enforce distance constraints.
 */
class CommonApplyConstraintsKernel : public ApplyConstraintsKernel {
public:
    CommonApplyConstraintsKernel(std::string name, const Platform& platform, ComputeContext& cc) : ApplyConstraintsKernel(name, platform),
            cc(cc), hasInitializedKernel(false) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     */
    void initialize(const System& system);
    /**
     * Update particle positions to enforce constraints.
     *
     * @param context    the context in which to execute this kernel
     * @param tol        the distance tolerance within which constraints must be satisfied.
     */
    void apply(ContextImpl& context, double tol);
    /**
     * Update particle velocities to enforce constraints.
     *
     * @param context    the context in which to execute this kernel
     * @param tol        the velocity tolerance within which constraints must be satisfied.
     */
    void applyToVelocities(ContextImpl& context, double tol);
private:
    ComputeContext& cc;
    bool hasInitializedKernel;
    ComputeKernel applyDeltasKernel;
};

/**
 * This kernel recomputes the positions of virtual sites.
 */
class CommonVirtualSitesKernel : public VirtualSitesKernel {
public:
    CommonVirtualSitesKernel(std::string name, const Platform& platform, ComputeContext& cc) : VirtualSitesKernel(name, platform), cc(cc) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     */
    void initialize(const System& system);
    /**
     * Compute the virtual site locations.
     *
     * @param context    the context in which to execute this kernel
     */
    void computePositions(ContextImpl& context);
private:
    ComputeContext& cc;
};

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
     * @param firstBond  the index of the first bond whose parameters might have changed
     * @param lastBond   the index of the last bond whose parameters might have changed
     */
    void copyParametersToContext(ContextImpl& context, const HarmonicBondForce& force, int firstBond, int lastBond);
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
     * @param firstBond  the index of the first bond whose parameters might have changed
     * @param lastBond   the index of the last bond whose parameters might have changed
     */
    void copyParametersToContext(ContextImpl& context, const CustomBondForce& force, int firstBond, int lastBond);
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
     * @param firstAngle the index of the first bond whose parameters might have changed
     * @param lastAngle  the index of the last bond whose parameters might have changed
     */
    void copyParametersToContext(ContextImpl& context, const HarmonicAngleForce& force, int firstAngle, int lastAngle);
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
     * @param firstAngle the index of the first bond whose parameters might have changed
     * @param lastAngle  the index of the last bond whose parameters might have changed
     */
    void copyParametersToContext(ContextImpl& context, const CustomAngleForce& force, int firstAngle, int lastAngle);
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
     * @param context      the context to copy parameters to
     * @param force        the PeriodicTorsionForce to copy the parameters from
     * @param firstTorsion the index of the first torsion whose parameters might have changed
     * @param lastTorsion  the index of the last torsion whose parameters might have changed
     */
    void copyParametersToContext(ContextImpl& context, const PeriodicTorsionForce& force, int firstTorsion, int lastTorsion);
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
     * @param context      the context to copy parameters to
     * @param force        the CustomTorsionForce to copy the parameters from
     * @param firstTorsion the index of the first torsion whose parameters might have changed
     * @param lastTorsion  the index of the last torsion whose parameters might have changed
     */
    void copyParametersToContext(ContextImpl& context, const CustomTorsionForce& force, int firstTorsion, int lastTorsion);
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
     * @param context        the context to copy parameters to
     * @param force          the CustomExternalForce to copy the parameters from
     * @param firstParticle  the index of the first particle whose parameters might have changed
     * @param lastParticle   the index of the last particle whose parameters might have changed
     */
    void copyParametersToContext(ContextImpl& context, const CustomExternalForce& force, int firstParticle, int lastParticle);
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
    std::vector<ComputeArray> tabulatedFunctionArrays;
    std::map<std::string, int> tabulatedFunctionUpdateCount;
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
    std::vector<ComputeArray> tabulatedFunctionArrays;
    std::map<std::string, int> tabulatedFunctionUpdateCount;
    std::vector<void*> groupForcesArgs;
    ComputeKernel computeCentersKernel, groupForcesKernel, applyForcesKernel;
    const System& system;
};

/**
 * This kernel is invoked by GBSAOBCForce to calculate the forces acting on the system.
 */
class CommonCalcGBSAOBCForceKernel : public CalcGBSAOBCForceKernel {
public:
    CommonCalcGBSAOBCForceKernel(std::string name, const Platform& platform, ComputeContext& cc) : CalcGBSAOBCForceKernel(name, platform), cc(cc),
            hasCreatedKernels(false) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the GBSAOBCForce this kernel will be used for
     */
    void initialize(const System& system, const GBSAOBCForce& force);
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
     * @param force      the GBSAOBCForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const GBSAOBCForce& force);
private:
    class ForceInfo;
    double prefactor, surfaceAreaFactor, cutoff;
    bool hasCreatedKernels;
    int maxTiles;
    ComputeContext& cc;
    ForceInfo* info;
    ComputeArray params, charges, bornSum, bornRadii, bornForce, obcChain;
    ComputeKernel computeBornSumKernel, reduceBornSumKernel, force1Kernel, reduceBornForceKernel;
};

/**
 * This kernel is invoked by GayBerneForce to calculate the forces acting on the system.
 */
class CommonCalcGayBerneForceKernel : public CalcGayBerneForceKernel {
public:
    CommonCalcGayBerneForceKernel(std::string name, const Platform& platform, ComputeContext& cc) : CalcGayBerneForceKernel(name, platform), cc(cc),
            hasInitializedKernels(false) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the GayBerneForce this kernel will be used for
     */
    void initialize(const System& system, const GayBerneForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the GayBerneForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const GayBerneForce& force);
private:
    class ForceInfo;
    class ReorderListener;
    void sortAtoms();
    ComputeContext& cc;
    ForceInfo* info;
    bool hasInitializedKernels;
    int numRealParticles, maxNeighborBlocks;
    GayBerneForce::NonbondedMethod nonbondedMethod;
    ComputeArray sortedParticles, axisParticleIndices, sigParams, epsParams;
    ComputeArray scale, exceptionParticles, exceptionParams;
    ComputeArray aMatrix, bMatrix, gMatrix;
    ComputeArray exclusions, exclusionStartIndex, blockCenter, blockBoundingBox;
    ComputeArray neighbors, neighborIndex, neighborBlockCount;
    ComputeArray sortedPos, torque;
    std::vector<bool> isRealParticle;
    std::vector<std::pair<int, int> > exceptionAtoms;
    std::vector<std::pair<int, int> > excludedPairs;
    ComputeKernel framesKernel, blockBoundsKernel, neighborsKernel, forceKernel, torqueKernel;
    ComputeEvent event;
};

/**
 * This kernel is invoked by CustomCVForce to calculate the forces acting on the system and the energy of the system.
 */
class CommonCalcCustomCVForceKernel : public CalcCustomCVForceKernel {
public:
    CommonCalcCustomCVForceKernel(std::string name, const Platform& platform, ComputeContext& cc) : CalcCustomCVForceKernel(name, platform),
            cc(cc), hasInitializedListeners(false) {
    }
    ~CommonCalcCustomCVForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomCVForce this kernel will be used for
     * @param innerContext   the context created by the CustomCVForce for computing collective variables
     */
    void initialize(const System& system, const CustomCVForce& force, ContextImpl& innerContext);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param innerContext   the context created by the CustomCVForce for computing collective variables
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, ContextImpl& innerContext, bool includeForces, bool includeEnergy);
    /**
     * Copy state information to the inner context.
     *
     * @param context        the context in which to execute this kernel
     * @param innerContext   the context created by the CustomCVForce for computing collective variables
     */
    void copyState(ContextImpl& context, ContextImpl& innerContext);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomCVForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const CustomCVForce& force);
    /**
     * Get the ComputeContext corresponding to the inner Context.
     */
    virtual ComputeContext& getInnerComputeContext(ContextImpl& innerContext) = 0;
private:
    class ForceInfo;
    class ReorderListener;
    class TabulatedFunctionWrapper;
    ComputeContext& cc;
    bool hasInitializedListeners;
    Lepton::CompiledExpression energyExpression;
    std::vector<std::string> variableNames, paramDerivNames, globalParameterNames;
    std::vector<Lepton::CompiledExpression> variableDerivExpressions;
    std::vector<Lepton::CompiledExpression> paramDerivExpressions;
    std::vector<ComputeArray> cvForces;
    std::vector<double> globalValues, cvValues;
    std::vector<Lepton::CustomFunction*> tabulatedFunctions;
    ComputeArray invAtomOrder;
    ComputeArray innerInvAtomOrder;
    ComputeKernel copyStateKernel, copyForcesKernel, addForcesKernel;
};

/**
 * This kernel is invoked by VerletIntegrator to take one time step.
 */
class CommonIntegrateVerletStepKernel : public IntegrateVerletStepKernel {
public:
    CommonIntegrateVerletStepKernel(std::string name, const Platform& platform, ComputeContext& cc) : IntegrateVerletStepKernel(name, platform), cc(cc),
            hasInitializedKernels(false) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the VerletIntegrator this kernel will be used for
     */
    void initialize(const System& system, const VerletIntegrator& integrator);
    /**
     * Execute the kernel.
     *
     * @param context    the context in which to execute this kernel
     * @param integrator the VerletIntegrator this kernel is being used for
     */
    void execute(ContextImpl& context, const VerletIntegrator& integrator);
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the VerletIntegrator this kernel is being used for
     */
    double computeKineticEnergy(ContextImpl& context, const VerletIntegrator& integrator);
private:
    ComputeContext& cc;
    bool hasInitializedKernels;
    ComputeKernel kernel1, kernel2;
};

/**
 * This kernel is invoked by LangevinMiddleIntegrator to take one time step.
 */
class CommonIntegrateLangevinMiddleStepKernel : public IntegrateLangevinMiddleStepKernel {
public:
    CommonIntegrateLangevinMiddleStepKernel(std::string name, const Platform& platform, ComputeContext& cc) : IntegrateLangevinMiddleStepKernel(name, platform), cc(cc),
            hasInitializedKernels(false) {
    }
    /**
     * Initialize the kernel, setting up the particle masses.
     * 
     * @param system     the System this kernel will be applied to
     * @param integrator the LangevinMiddleIntegrator this kernel will be used for
     */
    void initialize(const System& system, const LangevinMiddleIntegrator& integrator);
    /**
     * Execute the kernel.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the LangevinMiddleIntegrator this kernel is being used for
     */
    void execute(ContextImpl& context, const LangevinMiddleIntegrator& integrator);
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the LangevinMiddleIntegrator this kernel is being used for
     */
    double computeKineticEnergy(ContextImpl& context, const LangevinMiddleIntegrator& integrator);
private:
    ComputeContext& cc;
    double prevTemp, prevFriction, prevStepSize;
    bool hasInitializedKernels;
    ComputeArray params, oldDelta;
    ComputeKernel kernel1, kernel2, kernel3;
};

/**
 * This kernel is invoked by BrownianIntegrator to take one time step.
 */
class CommonIntegrateBrownianStepKernel : public IntegrateBrownianStepKernel {
public:
    CommonIntegrateBrownianStepKernel(std::string name, const Platform& platform, ComputeContext& cc) : IntegrateBrownianStepKernel(name, platform), cc(cc),
            hasInitializedKernels(false), prevTemp(-1), prevFriction(-1), prevStepSize(-1) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the BrownianIntegrator this kernel will be used for
     */
    void initialize(const System& system, const BrownianIntegrator& integrator);
    /**
     * Execute the kernel.
     *
     * @param context    the context in which to execute this kernel
     * @param integrator the BrownianIntegrator this kernel is being used for
     */
    void execute(ContextImpl& context, const BrownianIntegrator& integrator);
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the BrownianIntegrator this kernel is being used for
     */
    double computeKineticEnergy(ContextImpl& context, const BrownianIntegrator& integrator);
private:
    ComputeContext& cc;
    double prevTemp, prevFriction, prevStepSize;
    bool hasInitializedKernels;
    ComputeKernel kernel1, kernel2;
};

/**
 * This kernel is invoked by VerletIntegrator to take one time step.
 */
class CommonIntegrateVariableVerletStepKernel : public IntegrateVariableVerletStepKernel {
public:
    CommonIntegrateVariableVerletStepKernel(std::string name, const Platform& platform, ComputeContext& cc) : IntegrateVariableVerletStepKernel(name, platform), cc(cc),
            hasInitializedKernels(false) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the VariableVerletIntegrator this kernel will be used for
     */
    void initialize(const System& system, const VariableVerletIntegrator& integrator);
    /**
     * Execute the kernel.
     *
     * @param context    the context in which to execute this kernel
     * @param integrator the VariableVerletIntegrator this kernel is being used for
     * @param maxTime    the maximum time beyond which the simulation should not be advanced
     * @return the size of the step that was taken
     */
    double execute(ContextImpl& context, const VariableVerletIntegrator& integrator, double maxTime);
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the VariableVerletIntegrator this kernel is being used for
     */
    double computeKineticEnergy(ContextImpl& context, const VariableVerletIntegrator& integrator);
private:
    ComputeContext& cc;
    bool hasInitializedKernels;
    int blockSize;
    ComputeKernel kernel1, kernel2, selectSizeKernel;
};

/**
 * This kernel is invoked by VariableLangevinIntegrator to take one time step.
 */
class CommonIntegrateVariableLangevinStepKernel : public IntegrateVariableLangevinStepKernel {
public:
    CommonIntegrateVariableLangevinStepKernel(std::string name, const Platform& platform, ComputeContext& cc) : IntegrateVariableLangevinStepKernel(name, platform), cc(cc),
            hasInitializedKernels(false) {
    }
    /**
     * Initialize the kernel, setting up the particle masses.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the VariableLangevinIntegrator this kernel will be used for
     */
    void initialize(const System& system, const VariableLangevinIntegrator& integrator);
    /**
     * Execute the kernel.
     *
     * @param context    the context in which to execute this kernel
     * @param integrator the VariableLangevinIntegrator this kernel is being used for
     * @param maxTime    the maximum time beyond which the simulation should not be advanced
     * @return the size of the step that was taken
     */
    double execute(ContextImpl& context, const VariableLangevinIntegrator& integrator, double maxTime);
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the VariableLangevinIntegrator this kernel is being used for
     */
    double computeKineticEnergy(ContextImpl& context, const VariableLangevinIntegrator& integrator);
private:
    ComputeContext& cc;
    bool hasInitializedKernels;
    int blockSize;
    ComputeArray params, oldDelta;
    ComputeKernel kernel1, kernel2, kernel3, selectSizeKernel;
    double prevTemp, prevFriction, prevErrorTol;
};

/**
 * This kernel is invoked by DPDIntegrator to take one time step.
 */
class CommonIntegrateDPDStepKernel : public IntegrateDPDStepKernel {
public:
    CommonIntegrateDPDStepKernel(std::string name, const Platform& platform, ComputeContext& cc) : IntegrateDPDStepKernel(name, platform), cc(cc),
            hasInitializedKernels(false) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param integrator the DPDIntegrator this kernel will be used for
     */
    void initialize(const System& system, const DPDIntegrator& integrator);
    /**
     * Execute the kernel.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the DPDIntegrator this kernel is being used for
     */
    void execute(ContextImpl& context, const DPDIntegrator& integrator);
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the DPDIntegrator this kernel is being used for
     */
    double computeKineticEnergy(ContextImpl& context, const DPDIntegrator& integrator);
private:
    class ReorderListener;
    ComputeContext& cc;
    bool hasInitializedKernels;
    int numTypes, randomSeed, blockSize;
    double maxCutoff;
    ComputeArray particleType, pairParams, oldDelta, velDelta, tileCounter;
    ComputeKernel kernel1, kernel2, kernel3, kernel4;
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

/**
 * This kernel is invoked by RMSDForce to calculate the forces acting on the system and the energy of the system.
 */
class CommonCalcRMSDForceKernel : public CalcRMSDForceKernel {
public:
    CommonCalcRMSDForceKernel(std::string name, const Platform& platform, ComputeContext& cc) : CalcRMSDForceKernel(name, platform), cc(cc) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the RMSDForce this kernel will be used for
     */
    void initialize(const System& system, const RMSDForce& force);
    /**
     * Record the reference positions and particle indices.
     */
    void recordParameters(const RMSDForce& force);
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
     * This is the internal implementation of execute(), templatized on whether we're
     * using single or double precision.
     */
    template <class REAL>
    double executeImpl(ContextImpl& context);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the RMSDForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const RMSDForce& force);
private:
    class ForceInfo;
    ComputeContext& cc;
    ForceInfo* info;
    int blockSize;
    double sumNormRef;
    ComputeArray referencePos, particles, buffer;
    ComputeKernel kernel1, kernel2;
};

/**
 * This kernel is invoked by AndersenThermostat at the start of each time step to adjust the particle velocities.
 */
class CommonApplyAndersenThermostatKernel : public ApplyAndersenThermostatKernel {
public:
    CommonApplyAndersenThermostatKernel(std::string name, const Platform& platform, ComputeContext& cc) : ApplyAndersenThermostatKernel(name, platform), cc(cc) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param thermostat the AndersenThermostat this kernel will be used for
     */
    void initialize(const System& system, const AndersenThermostat& thermostat);
    /**
     * Execute the kernel.
     *
     * @param context    the context in which to execute this kernel
     */
    void execute(ContextImpl& context);
private:
    ComputeContext& cc;
    int randomSeed;
    ComputeArray atomGroups;
    ComputeKernel kernel;
};

/**
 * This kernel is invoked by MonteCarloBarostat to adjust the periodic box volume
 */
class CommonApplyMonteCarloBarostatKernel : public ApplyMonteCarloBarostatKernel {
public:
    CommonApplyMonteCarloBarostatKernel(std::string name, const Platform& platform, ComputeContext& cc) : ApplyMonteCarloBarostatKernel(name, platform), cc(cc),
            hasInitializedKernels(false) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system          the System this kernel will be applied to
     * @param barostat        the MonteCarloBarostat this kernel will be used for
     * @param rigidMolecules  whether molecules should be kept rigid while scaling coordinates
     * @param components      the number of box components the barostat operates one (1 for isotropic scaling,
     *                        3 for anisotropic, 6 for both lengths and angles)
     */
    void initialize(const System& system, const Force& barostat, int components, bool rigidMolecules=true);
    /**
     * Save the coordinates before attempting a Monte Carlo step.  This allows us to restore them
     * if the step is rejected.
     *
     * @param context    the context in which to execute this kernel
     */
    void saveCoordinates(ContextImpl& context);
    /**
     * Attempt a Monte Carlo step, scaling particle positions (or cluster centers) by a specified value.
     * This version scales the x, y, and z positions independently.
     * This is called BEFORE the periodic box size is modified.  It should begin by translating each particle
     * or cluster into the first periodic box, so that coordinates will still be correct after the box size
     * is changed.
     *
     * @param context    the context in which to execute this kernel
     * @param scaleX     the scale factor by which to multiply particle x-coordinate
     * @param scaleY     the scale factor by which to multiply particle y-coordinate
     * @param scaleZ     the scale factor by which to multiply particle z-coordinate
     */
    void scaleCoordinates(ContextImpl& context, double scaleX, double scaleY, double scaleZ);
    /**
     * Reject the most recent Monte Carlo step, restoring the particle positions to where they were when
     * saveCoordinates() was last called.
     *
     * @param context    the context in which to execute this kernel
     */
    void restoreCoordinates(ContextImpl& context);
    /**
     * Compute the kinetic energy of the system.  If initialize() was called with rigidMolecules=true, this
     * should include only the translational center of mass motion of molecules.  Otherwise it should include
     * the total kinetic energy of all particles.  This is used when computing instantaneous pressure.
     * 
     * @param context    the context in which to execute this kernel
     * @param ke         a vector to store the kinetic energy components into.  On output, its length will
     *                   equal the number of components passed to initialize().
     */
    void computeKineticEnergy(ContextImpl& context, std::vector<double>& ke);
private:
    ComputeContext& cc;
    bool hasInitializedKernels, rigidMolecules, atomsWereReordered;
    int numMolecules, components;
    ComputeArray savedPositions, savedFloatForces, savedLongForces, savedVelocities;
    ComputeArray moleculeAtoms, moleculeStartIndex;
    std::vector<ComputeArray> energyBuffers;
    ComputeKernel kernel, kineticEnergyKernel;
    std::vector<int> lastAtomOrder;
    std::vector<mm_int4> lastPosCellOffsets;
};

/**
 * This kernel is invoked by ATMForce to calculate the forces acting on the system and the energy of the system.
 */
class CommonCalcATMForceKernel : public CalcATMForceKernel {
public:
    CommonCalcATMForceKernel(std::string name, const Platform& platform, ComputeContext& cc): CalcATMForceKernel(name, platform), hasInitializedKernel(false), cc(cc) {
    }

    ~CommonCalcATMForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the ATMForce this kernel will be used for
     */
    void initialize(const System& system, const ATMForce& force);
    /**
     * Scale the forces from the inner contexts and apply them to the main context.
     *
     * @param context        the context in which to execute this kernel
     * @param innerContext1  the first inner context
     * @param innerContext2  the second inner context
     * @param dEdu0          the derivative of the final energy with respect to the first inner context's energy
     * @param dEdu1          the derivative of the final energy with respect to the second inner context's energy
     * @param energyParamDerivs  derivatives of the final energy with respect to global parameters
     */
    void applyForces(ContextImpl& context, ContextImpl& innerContext0, ContextImpl& innerContext1,
                     double dEdu0, double dEdu1, const std::map<std::string, double>& energyParamDerivs);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the ATMForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const ATMForce& force);
    /**
     * Copy state information to the inner contexts.
     *
     * @param context        the context in which to execute this kernel
     * @param innerContext1  the first context created by the ATMForce for computing displaced energy
     * @param innerContext2  the second context created by the ATMForce for computing displaced energy
     */
    void copyState(ContextImpl& context, ContextImpl& innerContext1, ContextImpl& innerContext2);
    /**
     * Get the ComputeContext corresponding to the inner Context.
     */
    virtual ComputeContext& getInnerComputeContext(ContextImpl& innerContext) = 0;
    
private:
    class ReorderListener;
    
    void initKernels(ContextImpl& context, ContextImpl& innerContext0, ContextImpl& innerContext1);
    
    bool hasInitializedKernel;
    ComputeContext& cc;
    ComputeArray displ1;
    ComputeArray displ0;
    ComputeArray invAtomOrder, inner0InvAtomOrder, inner1InvAtomOrder;
    ComputeKernel copyStateKernel;
    ComputeKernel hybridForceKernel;

    int numParticles;
};

/**
 * This kernel is invoked by CustomCPPForce to calculate the forces acting on the system and the energy of the system.
 */
class CommonCalcCustomCPPForceKernel : public CalcCustomCPPForceKernel {
public:
    CommonCalcCustomCPPForceKernel(std::string name, const Platform& platform, OpenMM::ContextImpl& contextImpl, ComputeContext& cc) :
            CalcCustomCPPForceKernel(name, platform), contextImpl(contextImpl), cc(cc), force(NULL) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomCPPForceImpl this kernel will be used for
     */
    void initialize(const System& system, CustomCPPForceImpl& force);
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
     * The is called by the pre-computation to start the calculation running.
     */
    void beginComputation(bool includeForces, bool includeEnergy, int groups);
    /**
     * This is called by the worker thread to do the computation.
     */
    void executeOnWorkerThread(bool includeForces);
    /**
     * This is called by the post-computation to add the forces to the main array.
     */
    double addForces(bool includeForces, bool includeEnergy, int groups);
private:
    class ExecuteTask;
    class StartCalculationPreComputation;
    class AddForcesPostComputation;
    OpenMM::ContextImpl& contextImpl;
    ComputeContext& cc;
    CustomCPPForceImpl* force;
    ComputeArray forcesArray;
    ComputeKernel addForcesKernel;
    std::vector<Vec3> positionsVec, forcesVec;
    std::vector<float> floatForces;
    int forceGroupFlag;
    double energy;
};

} // namespace OpenMM

#endif /*OPENMM_COMMONKERNELS_H_*/
