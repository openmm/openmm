#ifndef OPENMM_KERNELS_H_
#define OPENMM_KERNELS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2018 Stanford University and the Authors.      *
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

#include "openmm/AndersenThermostat.h"
#include "openmm/BrownianIntegrator.h"
#include "openmm/CMAPTorsionForce.h"
#include "openmm/CMMotionRemover.h"
#include "openmm/CustomAngleForce.h"
#include "openmm/CustomBondForce.h"
#include "openmm/CustomCentroidBondForce.h"
#include "openmm/CustomCVForce.h"
#include "openmm/CustomCompoundBondForce.h"
#include "openmm/CustomExternalForce.h"
#include "openmm/CustomGBForce.h"
#include "openmm/CustomHbondForce.h"
#include "openmm/CustomIntegrator.h"
#include "openmm/CustomNonbondedForce.h"
#include "openmm/CustomManyParticleForce.h"
#include "openmm/CustomTorsionForce.h"
#include "openmm/GayBerneForce.h"
#include "openmm/GBSAOBCForce.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/KernelImpl.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/MonteCarloBarostat.h"
#include "openmm/PeriodicTorsionForce.h"
#include "openmm/RBTorsionForce.h"
#include "openmm/RMSDForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/VariableLangevinIntegrator.h"
#include "openmm/VariableVerletIntegrator.h"
#include "openmm/VerletIntegrator.h"
#include <iosfwd>
#include <set>
#include <string>
#include <vector>

namespace OpenMM {

/**
 * This kernel is invoked at the beginning and end of force and energy computations.  It gives the
 * Platform a chance to clear buffers and do other initialization at the beginning, and to do any
 * necessary work at the end to determine the final results.
 */
class CalcForcesAndEnergyKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcForcesAndEnergyKernel";
    }
    CalcForcesAndEnergyKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     */
    virtual void initialize(const System& system) = 0;
    /**
     * This is called at the beginning of each force/energy computation, before calcForcesAndEnergy() has been called on
     * any ForceImpl.
     * 
     * @param context       the context in which to execute this kernel
     * @param includeForce  true if forces should be computed
     * @param includeEnergy true if potential energy should be computed
     * @param groups        a set of bit flags for which force groups to include
     */
    virtual void beginComputation(ContextImpl& context, bool includeForce, bool includeEnergy, int groups) = 0;
    /**
     * This is called at the end of each force/energy computation, after calcForcesAndEnergy() has been called on
     * every ForceImpl.
     *
     * @param context       the context in which to execute this kernel
     * @param includeForce  true if forces should be computed
     * @param includeEnergy true if potential energy should be computed
     * @param groups        a set of bit flags for which force groups to include
     * @param valid         the method may set this to false to indicate the results are invalid and the force/energy
     *                      calculation should be repeated
     * @return the potential energy of the system.  This value is added to all values returned by ForceImpls'
     * calcForcesAndEnergy() methods.  That is, each force kernel may <i>either</i> return its contribution to the
     * energy directly, <i>or</i> add it to an internal buffer so that it will be included here.
     */
    virtual double finishComputation(ContextImpl& context, bool includeForce, bool includeEnergy, int groups, bool& valid) = 0;
};

/**
 * This kernel provides methods for setting and retrieving various state data: time, positions,
 * velocities, and forces.
 */
class UpdateStateDataKernel : public KernelImpl {
public:
    static std::string Name() {
        return "UpdateTime";
    }
    UpdateStateDataKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     */
    virtual void initialize(const System& system) = 0;
    /**
     * Get the current time (in picoseconds).
     *
     * @param context    the context in which to execute this kernel
     */
    virtual double getTime(const ContextImpl& context) const = 0;
    /**
     * Set the current time (in picoseconds).
     *
     * @param context    the context in which to execute this kernel
     * @param time       the time
     */
    virtual void setTime(ContextImpl& context, double time) = 0;
    /**
     * Get the positions of all particles.
     *
     * @param positions  on exit, this contains the particle positions
     */
    virtual void getPositions(ContextImpl& context, std::vector<Vec3>& positions) = 0;
    /**
     * Set the positions of all particles.
     *
     * @param positions  a vector containing the particle positions
     */
    virtual void setPositions(ContextImpl& context, const std::vector<Vec3>& positions) = 0;
    /**
     * Get the velocities of all particles.
     *
     * @param velocities  on exit, this contains the particle velocities
     */
    virtual void getVelocities(ContextImpl& context, std::vector<Vec3>& velocities) = 0;
    /**
     * Set the velocities of all particles.
     *
     * @param velocities  a vector containing the particle velocities
     */
    virtual void setVelocities(ContextImpl& context, const std::vector<Vec3>& velocities) = 0;
    /**
     * Get the current forces on all particles.
     *
     * @param forces  on exit, this contains the forces
     */
    virtual void getForces(ContextImpl& context, std::vector<Vec3>& forces) = 0;
    /**
     * Get the current derivatives of the energy with respect to context parameters.
     *
     * @param derivs  on exit, this contains the derivatives
     */
    virtual void getEnergyParameterDerivatives(ContextImpl& context, std::map<std::string, double>& derivs) = 0;
    /**
     * Get the current periodic box vectors.
     *
     * @param a      on exit, this contains the vector defining the first edge of the periodic box
     * @param b      on exit, this contains the vector defining the second edge of the periodic box
     * @param c      on exit, this contains the vector defining the third edge of the periodic box
     */
    virtual void getPeriodicBoxVectors(ContextImpl& context, Vec3& a, Vec3& b, Vec3& c) const = 0;
    /**
     * Set the current periodic box vectors.
     *
     * @param a      the vector defining the first edge of the periodic box
     * @param b      the vector defining the second edge of the periodic box
     * @param c      the vector defining the third edge of the periodic box
     */
    virtual void setPeriodicBoxVectors(ContextImpl& context, const Vec3& a, const Vec3& b, const Vec3& c) = 0;
    /**
     * Create a checkpoint recording the current state of the Context.
     * 
     * @param stream    an output stream the checkpoint data should be written to
     */
    virtual void createCheckpoint(ContextImpl& context, std::ostream& stream) = 0;
    /**
     * Load a checkpoint that was written by createCheckpoint().
     * 
     * @param stream    an input stream the checkpoint data should be read from
     */
    virtual void loadCheckpoint(ContextImpl& context, std::istream& stream) = 0;
};

/**
 * This kernel modifies the positions of particles to enforce distance constraints.
 */
class ApplyConstraintsKernel : public KernelImpl {
public:
    static std::string Name() {
        return "ApplyConstraints";
    }
    ApplyConstraintsKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     */
    virtual void initialize(const System& system) = 0;
    /**
     * Update particle positions to enforce constraints.
     *
     * @param context    the context in which to execute this kernel
     * @param tol        the distance tolerance within which constraints must be satisfied.
     */
    virtual void apply(ContextImpl& context, double tol) = 0;
    /**
     * Update particle velocities to enforce constraints.
     *
     * @param context    the context in which to execute this kernel
     * @param tol        the velocity tolerance within which constraints must be satisfied.
     */
    virtual void applyToVelocities(ContextImpl& context, double tol) = 0;
};

/**
 * This kernel recomputes the positions of virtual sites.
 */
class VirtualSitesKernel : public KernelImpl {
public:
    static std::string Name() {
        return "VirtualSites";
    }
    VirtualSitesKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     */
    virtual void initialize(const System& system) = 0;
    /**
     * Compute the virtual site locations.
     *
     * @param context    the context in which to execute this kernel
     */
    virtual void computePositions(ContextImpl& context) = 0;
};

/**
 * This kernel is invoked by HarmonicBondForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcHarmonicBondForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcHarmonicBondForce";
    }
    CalcHarmonicBondForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the HarmonicBondForce this kernel will be used for
     */
    virtual void initialize(const System& system, const HarmonicBondForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     * 
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the HarmonicBondForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const HarmonicBondForce& force) = 0;
};

/**
 * This kernel is invoked by CustomBondForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcCustomBondForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcCustomBondForce";
    }
    CalcCustomBondForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomBondForce this kernel will be used for
     */
    virtual void initialize(const System& system, const CustomBondForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomBondForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const CustomBondForce& force) = 0;
};

/**
 * This kernel is invoked by HarmonicAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcHarmonicAngleForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcHarmonicAngleForce";
    }
    CalcHarmonicAngleForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the HarmonicAngleForce this kernel will be used for
     */
    virtual void initialize(const System& system, const HarmonicAngleForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the HarmonicAngleForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const HarmonicAngleForce& force) = 0;
};

/**
 * This kernel is invoked by CustomAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcCustomAngleForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcCustomAngleForce";
    }
    CalcCustomAngleForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomAngleForce this kernel will be used for
     */
    virtual void initialize(const System& system, const CustomAngleForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomAngleForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const CustomAngleForce& force) = 0;
};

/**
 * This kernel is invoked by PeriodicTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcPeriodicTorsionForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcPeriodicTorsionForce";
    }
    CalcPeriodicTorsionForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the PeriodicTorsionForce this kernel will be used for
     */
    virtual void initialize(const System& system, const PeriodicTorsionForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the PeriodicTorsionForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const PeriodicTorsionForce& force) = 0;
};

/**
 * This kernel is invoked by RBTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcRBTorsionForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcRBTorsionForce";
    }
    CalcRBTorsionForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the RBTorsionForce this kernel will be used for
     */
    virtual void initialize(const System& system, const RBTorsionForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the RBTorsionForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const RBTorsionForce& force) = 0;
};

/**
 * This kernel is invoked by CMAPTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcCMAPTorsionForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcCMAPTorsionForce";
    }
    CalcCMAPTorsionForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CMAPTorsionForce this kernel will be used for
     */
    virtual void initialize(const System& system, const CMAPTorsionForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CMAPTorsionForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const CMAPTorsionForce& force) = 0;
};

/**
 * This kernel is invoked by CustomTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcCustomTorsionForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcCustomTorsionForce";
    }
    CalcCustomTorsionForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomTorsionForce this kernel will be used for
     */
    virtual void initialize(const System& system, const CustomTorsionForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomTorsionForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const CustomTorsionForce& force) = 0;
};

/**
 * This kernel is invoked by NonbondedForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcNonbondedForceKernel : public KernelImpl {
public:
    enum NonbondedMethod {
        NoCutoff = 0,
        CutoffNonPeriodic = 1,
        CutoffPeriodic = 2,
        Ewald = 3,
        PME = 4,
        LJPME = 5
    };
    static std::string Name() {
        return "CalcNonbondedForce";
    }
    CalcNonbondedForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the NonbondedForce this kernel will be used for
     */
    virtual void initialize(const System& system, const NonbondedForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @param includeDirect  true if direct space interactions should be included
     * @param includeReciprocal  true if reciprocal space interactions should be included
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy, bool includeDirect, bool includeReciprocal) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the NonbondedForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const NonbondedForce& force) = 0;
    /**
     * Get the parameters being used for PME.
     *
     * @param alpha   the separation parameter
     * @param nx      the number of grid points along the X axis
     * @param ny      the number of grid points along the Y axis
     * @param nz      the number of grid points along the Z axis
     */
    virtual void getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const = 0;
    /**
     * Get the parameters being used for the dispersion terms in LJPME.
     *
     * @param alpha   the separation parameter
     * @param nx      the number of grid points along the X axis
     * @param ny      the number of grid points along the Y axis
     * @param nz      the number of grid points along the Z axis
     */
    virtual void getLJPMEParameters(double& alpha, int& nx, int& ny, int& nz) const = 0;
};

/**
 * This kernel is invoked by CustomNonbondedForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcCustomNonbondedForceKernel : public KernelImpl {
public:
    enum NonbondedMethod {
        NoCutoff = 0,
        CutoffNonPeriodic = 1,
        CutoffPeriodic = 2
    };
    static std::string Name() {
        return "CalcCustomNonbondedForce";
    }
    CalcCustomNonbondedForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomNonbondedForce this kernel will be used for
     */
    virtual void initialize(const System& system, const CustomNonbondedForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomNonbondedForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const CustomNonbondedForce& force) = 0;
};

/**
 * This kernel is invoked by GBSAOBCForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcGBSAOBCForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcGBSAOBCForce";
    }
    CalcGBSAOBCForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the GBSAOBCForce this kernel will be used for
     */
    virtual void initialize(const System& system, const GBSAOBCForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the GBSAOBCForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const GBSAOBCForce& force) = 0;
};

/**
 * This kernel is invoked by CustomGBForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcCustomGBForceKernel : public KernelImpl {
public:
    enum NonbondedMethod {
        NoCutoff = 0,
        CutoffNonPeriodic = 1,
        CutoffPeriodic = 2
    };
    static std::string Name() {
        return "CalcCustomGBForce";
    }
    CalcCustomGBForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomGBForce this kernel will be used for
     */
    virtual void initialize(const System& system, const CustomGBForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomGBForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const CustomGBForce& force) = 0;
};

/**
 * This kernel is invoked by CustomExternalForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcCustomExternalForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcCustomExternalForce";
    }
    CalcCustomExternalForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomExternalForce this kernel will be used for
     */
    virtual void initialize(const System& system, const CustomExternalForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomExternalForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const CustomExternalForce& force) = 0;
};

/**
 * This kernel is invoked by CustomHbondForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcCustomHbondForceKernel : public KernelImpl {
public:
    enum NonbondedMethod {
        NoCutoff = 0,
        CutoffNonPeriodic = 1,
        CutoffPeriodic = 2
    };
    static std::string Name() {
        return "CalcCustomHbondForce";
    }
    CalcCustomHbondForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomHbondForce this kernel will be used for
     */
    virtual void initialize(const System& system, const CustomHbondForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomHbondForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const CustomHbondForce& force) = 0;
};

/**
 * This kernel is invoked by CustomCentroidBondForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcCustomCentroidBondForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcCustomCentroidBondForce";
    }
    CalcCustomCentroidBondForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomCentroidBondForce this kernel will be used for
     */
    virtual void initialize(const System& system, const CustomCentroidBondForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomCentroidBondForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const CustomCentroidBondForce& force) = 0;
};

/**
 * This kernel is invoked by CustomCompoundBondForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcCustomCompoundBondForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcCustomCompoundBondForce";
    }
    CalcCustomCompoundBondForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomCompoundBondForce this kernel will be used for
     */
    virtual void initialize(const System& system, const CustomCompoundBondForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomCompoundBondForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const CustomCompoundBondForce& force) = 0;
};

/**
 * This kernel is invoked by CustomManyParticleForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcCustomManyParticleForceKernel : public KernelImpl {
public:
    enum NonbondedMethod {
        NoCutoff = 0,
        CutoffNonPeriodic = 1,
        CutoffPeriodic = 2
    };
    static std::string Name() {
        return "CalcCustomManyParticleForce";
    }
    CalcCustomManyParticleForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomManyParticleForce this kernel will be used for
     */
    virtual void initialize(const System& system, const CustomManyParticleForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomManyParticleForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const CustomManyParticleForce& force) = 0;
};

/**
 * This kernel is invoked by GayBerneForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcGayBerneForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcGayBerneForce";
    }
    CalcGayBerneForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the GayBerneForce this kernel will be used for
     */
    virtual void initialize(const System& system, const GayBerneForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the GayBerneForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const GayBerneForce& force) = 0;
};

/**
 * This kernel is invoked by CustomCVForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcCustomCVForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcCustomCVForce";
    }
    CalcCustomCVForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomCVForce this kernel will be used for
     * @param innerContext   the context created by the CustomCVForce for computing collective variables
     */
    virtual void initialize(const System& system, const CustomCVForce& force, ContextImpl& innerContext) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param innerContext   the context created by the CustomCVForce for computing collective variables
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, ContextImpl& innerContext, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy state information to the inner context.
     *
     * @param context        the context in which to execute this kernel
     * @param innerContext   the context created by the CustomCVForce for computing collective variables
     */
    virtual void copyState(ContextImpl& context, ContextImpl& innerContext) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomCVForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const CustomCVForce& force) = 0;
};

/**
 * This kernel is invoked by RMSDForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcRMSDForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcRMSDForce";
    }
    CalcRMSDForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the RMSDForce this kernel will be used for
     */
    virtual void initialize(const System& system, const RMSDForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the RMSDForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const RMSDForce& force) = 0;
};

/**
 * This kernel is invoked by VerletIntegrator to take one time step.
 */
class IntegrateVerletStepKernel : public KernelImpl {
public:
    static std::string Name() {
        return "IntegrateVerletStep";
    }
    IntegrateVerletStepKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param integrator the VerletIntegrator this kernel will be used for
     */
    virtual void initialize(const System& system, const VerletIntegrator& integrator) = 0;
    /**
     * Execute the kernel.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the VerletIntegrator this kernel is being used for
     */
    virtual void execute(ContextImpl& context, const VerletIntegrator& integrator) = 0;
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the VerletIntegrator this kernel is being used for
     */
    virtual double computeKineticEnergy(ContextImpl& context, const VerletIntegrator& integrator) = 0;
};

/**
 * This kernel is invoked by LangevinIntegrator to take one time step.
 */
class IntegrateLangevinStepKernel : public KernelImpl {
public:
    static std::string Name() {
        return "IntegrateLangevinStep";
    }
    IntegrateLangevinStepKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param integrator the LangevinIntegrator this kernel will be used for
     */
    virtual void initialize(const System& system, const LangevinIntegrator& integrator) = 0;
    /**
     * Execute the kernel.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the LangevinIntegrator this kernel is being used for
     */
    virtual void execute(ContextImpl& context, const LangevinIntegrator& integrator) = 0;
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the LangevinIntegrator this kernel is being used for
     */
    virtual double computeKineticEnergy(ContextImpl& context, const LangevinIntegrator& integrator) = 0;
};

/**
 * This kernel is invoked by BrownianIntegrator to take one time step.
 */
class IntegrateBrownianStepKernel : public KernelImpl {
public:
    static std::string Name() {
        return "IntegrateBrownianStep";
    }
    IntegrateBrownianStepKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param integrator the BrownianIntegrator this kernel will be used for
     */
    virtual void initialize(const System& system, const BrownianIntegrator& integrator) = 0;
    /**
     * Execute the kernel.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the BrownianIntegrator this kernel is being used for
     */
    virtual void execute(ContextImpl& context, const BrownianIntegrator& integrator) = 0;
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the BrownianIntegrator this kernel is being used for
     */
    virtual double computeKineticEnergy(ContextImpl& context, const BrownianIntegrator& integrator) = 0;
};

/**
 * This kernel is invoked by VariableLangevinIntegrator to take one time step.
 */
class IntegrateVariableLangevinStepKernel : public KernelImpl {
public:
    static std::string Name() {
        return "IntegrateVariableLangevinStep";
    }
    IntegrateVariableLangevinStepKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the VariableLangevinIntegrator this kernel will be used for
     */
    virtual void initialize(const System& system, const VariableLangevinIntegrator& integrator) = 0;
    /**
     * Execute the kernel.
     *
     * @param context    the context in which to execute this kernel
     * @param integrator the VariableLangevinIntegrator this kernel is being used for
     * @param maxTime    the maximum time beyond which the simulation should not be advanced
     * @return the size of the step that was taken
     */
    virtual double execute(ContextImpl& context, const VariableLangevinIntegrator& integrator, double maxTime) = 0;
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the VariableLangevinIntegrator this kernel is being used for
     */
    virtual double computeKineticEnergy(ContextImpl& context, const VariableLangevinIntegrator& integrator) = 0;
};

/**
 * This kernel is invoked by VariableVerletIntegrator to take one time step.
 */
class IntegrateVariableVerletStepKernel : public KernelImpl {
public:
    static std::string Name() {
        return "IntegrateVariableVerletStep";
    }
    IntegrateVariableVerletStepKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the VariableVerletIntegrator this kernel will be used for
     */
    virtual void initialize(const System& system, const VariableVerletIntegrator& integrator) = 0;
    /**
     * Execute the kernel.
     *
     * @param context    the context in which to execute this kernel
     * @param integrator the VariableVerletIntegrator this kernel is being used for
     * @param maxTime    the maximum time beyond which the simulation should not be advanced
     * @return the size of the step that was taken
     */
    virtual double execute(ContextImpl& context, const VariableVerletIntegrator& integrator, double maxTime) = 0;
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the VariableVerletIntegrator this kernel is being used for
     */
    virtual double computeKineticEnergy(ContextImpl& context, const VariableVerletIntegrator& integrator) = 0;
};

/**
 * This kernel is invoked by CustomIntegrator to take one time step.
 */
class IntegrateCustomStepKernel : public KernelImpl {
public:
    static std::string Name() {
        return "IntegrateCustomStep";
    }
    IntegrateCustomStepKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param integrator the CustomIntegrator this kernel will be used for
     */
    virtual void initialize(const System& system, const CustomIntegrator& integrator) = 0;
    /**
     * Execute the kernel.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the CustomIntegrator this kernel is being used for
     * @param forcesAreValid if the context has been modified since the last time step, this will be
     *                       false to show that cached forces are invalid and must be recalculated.
     *                       On exit, this should specify whether the cached forces are valid at the
     *                       end of the step.
     */
    virtual void execute(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid) = 0;
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the CustomIntegrator this kernel is being used for
     * @param forcesAreValid if the context has been modified since the last time step, this will be
     *                       false to show that cached forces are invalid and must be recalculated.
     *                       On exit, this should specify whether the cached forces are valid at the
     *                       end of the step.
     */
    virtual double computeKineticEnergy(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid) = 0;
    /**
     * Get the values of all global variables.
     *
     * @param context   the context in which to execute this kernel
     * @param values    on exit, this contains the values
     */
    virtual void getGlobalVariables(ContextImpl& context, std::vector<double>& values) const = 0;
    /**
     * Set the values of all global variables.
     *
     * @param context   the context in which to execute this kernel
     * @param values    a vector containing the values
     */
    virtual void setGlobalVariables(ContextImpl& context, const std::vector<double>& values) = 0;
    /**
     * Get the values of a per-DOF variable.
     *
     * @param context   the context in which to execute this kernel
     * @param variable  the index of the variable to get
     * @param values    on exit, this contains the values
     */
    virtual void getPerDofVariable(ContextImpl& context, int variable, std::vector<Vec3>& values) const = 0;
    /**
     * Set the values of a per-DOF variable.
     *
     * @param context   the context in which to execute this kernel
     * @param variable  the index of the variable to get
     * @param values    a vector containing the values
     */
    virtual void setPerDofVariable(ContextImpl& context, int variable, const std::vector<Vec3>& values) = 0;
};

/**
 * This kernel is invoked by AndersenThermostat at the start of each time step to adjust the particle velocities.
 */
class ApplyAndersenThermostatKernel : public KernelImpl {
public:
    static std::string Name() {
        return "ApplyAndersenThermostat";
    }
    ApplyAndersenThermostatKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param thermostat the AndersenThermostat this kernel will be used for
     */
    virtual void initialize(const System& system, const AndersenThermostat& thermostat) = 0;
    /**
     * Execute the kernel.
     * 
     * @param context    the context in which to execute this kernel
     */
    virtual void execute(ContextImpl& context) = 0;
};

/**
 * This kernel is invoked by MonteCarloBarostat to adjust the periodic box volume
 */
class ApplyMonteCarloBarostatKernel : public KernelImpl {
public:
    static std::string Name() {
        return "ApplyMonteCarloBarostat";
    }
    ApplyMonteCarloBarostatKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param barostat   the MonteCarloBarostat this kernel will be used for
     */
    virtual void initialize(const System& system, const Force& barostat) = 0;
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
    virtual void scaleCoordinates(ContextImpl& context, double scaleX, double scaleY, double scaleZ) = 0;
    /**
     * Reject the most recent Monte Carlo step, restoring the particle positions to where they were before
     * scaleCoordinates() was last called.
     *
     * @param context    the context in which to execute this kernel
     */
    virtual void restoreCoordinates(ContextImpl& context) = 0;
};

/**
 * This kernel is invoked to remove center of mass motion from the system.
 */
class RemoveCMMotionKernel : public KernelImpl {
public:
    static std::string Name() {
        return "RemoveCMMotion";
    }
    RemoveCMMotionKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the CMMotionRemover this kernel will be used for
     */
    virtual void initialize(const System& system, const CMMotionRemover& force) = 0;
    /**
     * Execute the kernel.
     * 
     * @param context    the context in which to execute this kernel
     */
    virtual void execute(ContextImpl& context) = 0;
};

/**
 * This kernel performs the reciprocal space calculation for PME.  In most cases, this
 * calculation is done directly by CalcNonbondedForceKernel so this kernel is unneeded.
 * In some cases it may want to outsource the work to a different kernel.  In particular,
 * GPU based platforms sometimes use a CPU based implementation provided by a separate
 * plugin.
 */
class CalcPmeReciprocalForceKernel : public KernelImpl {
public:
    class IO;
    static std::string Name() {
        return "CalcPmeReciprocalForce";
    }
    CalcPmeReciprocalForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param gridx        the x size of the PME grid
     * @param gridy        the y size of the PME grid
     * @param gridz        the z size of the PME grid
     * @param numParticles the number of particles in the system
     * @param alpha        the Ewald blending parameter
     * @param deterministic whether it should attempt to make the resulting forces deterministic
     */
    virtual void initialize(int gridx, int gridy, int gridz, int numParticles, double alpha, bool deterministic) = 0;
    /**
     * Begin computing the force and energy.
     *
     * @param io                  an object that coordinates data transfer
     * @param periodicBoxVectors  the vectors defining the periodic box (measured in nm)
     * @param includeEnergy       true if potential energy should be computed
     */
    virtual void beginComputation(IO& io, const Vec3* periodicBoxVectors, bool includeEnergy) = 0;
    /**
     * Finish computing the force and energy.
     * 
     * @param io   an object that coordinates data transfer
     * @return the potential energy due to the PME reciprocal space interactions
     */
    virtual double finishComputation(IO& io) = 0;
    /**
     * Get the parameters being used for PME.
     * 
     * @param alpha   the separation parameter
     * @param nx      the number of grid points along the X axis
     * @param ny      the number of grid points along the Y axis
     * @param nz      the number of grid points along the Z axis
     */
    virtual void getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const = 0;
};

/**
 * Any class that uses CalcPmeReciprocalForceKernel should create an implementation of this
 * class, then pass it to the kernel to manage communication with it.
 */
class CalcPmeReciprocalForceKernel::IO {
public:
    /**
     * Get a pointer to the atom charges and positions.  This array should contain four
     * elements for each atom: x, y, z, and q in that order.
     */
    virtual float* getPosq() = 0;
    /**
     * Record the forces calculated by the kernel.
     * 
     * @param force    an array containing four elements for each atom.  The first three
     *                 are the x, y, and z components of the force, while the fourth element
     *                 should be ignored.
     */
    virtual void setForce(float* force) = 0;
};


/**
 * This kernel performs the dispersion reciprocal space calculation for LJPME.  In most cases, this
 * calculation is done directly by CalcNonbondedForceKernel so this kernel is unneeded.
 * In some cases it may want to outsource the work to a different kernel.  In particular,
 * GPU based platforms sometimes use a CPU based implementation provided by a separate
 * plugin.
 */
class CalcDispersionPmeReciprocalForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcDispersionPmeReciprocalForce";
    }
    CalcDispersionPmeReciprocalForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param gridx        the x size of the PME grid
     * @param gridy        the y size of the PME grid
     * @param gridz        the z size of the PME grid
     * @param numParticles the number of particles in the system
     * @param alpha        the Ewald blending parameter
     * @param deterministic whether it should attempt to make the resulting forces deterministic
     */
    virtual void initialize(int gridx, int gridy, int gridz, int numParticles, double alpha, bool deterministic) = 0;
    /**
     * Begin computing the force and energy.
     *
     * @param io                  an object that coordinates data transfer
     * @param periodicBoxVectors  the vectors defining the periodic box (measured in nm)
     * @param includeEnergy       true if potential energy should be computed
     */
    virtual void beginComputation(CalcPmeReciprocalForceKernel::IO& io, const Vec3* periodicBoxVectors, bool includeEnergy) = 0;
    /**
     * Finish computing the force and energy.
     * 
     * @param io   an object that coordinates data transfer
     * @return the potential energy due to the PME reciprocal space interactions
     */
    virtual double finishComputation(CalcPmeReciprocalForceKernel::IO& io) = 0;
    /**
     * Get the parameters being used for PME.
     * 
     * @param alpha   the separation parameter
     * @param nx      the number of grid points along the X axis
     * @param ny      the number of grid points along the Y axis
     * @param nz      the number of grid points along the Z axis
     */
    virtual void getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const = 0;
};

} // namespace OpenMM

#endif /*OPENMM_KERNELS_H_*/
