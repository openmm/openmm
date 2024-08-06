#ifndef OPENMM_REFERENCEKERNELS_H_
#define OPENMM_REFERENCEKERNELS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2024 Stanford University and the Authors.      *
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

#include "ReferencePlatform.h"
#include "openmm/kernels.h"
#include "openmm/internal/CustomCPPForceImpl.h"
#include "openmm/internal/CustomNonbondedForceImpl.h"
#include "SimTKOpenMMRealType.h"
#include "ReferenceNeighborList.h"
#include "lepton/CompiledExpression.h"
#include "lepton/CustomFunction.h"
#include <array>
#include <utility>

namespace OpenMM {

class ReferenceObc;
class ReferenceAndersenThermostat;
class ReferenceLangevinMiddleDynamics;
class ReferenceCustomBondIxn;
class ReferenceCustomAngleIxn;
class ReferenceCustomTorsionIxn;
class ReferenceCustomExternalIxn;
class ReferenceCustomCentroidBondIxn;
class ReferenceCustomCompoundBondIxn;
class ReferenceCustomCVForce;
class ReferenceCustomHbondIxn;
class ReferenceCustomManyParticleIxn;
class ReferenceGayBerneForce;
class ReferenceBrownianDynamics;
class ReferenceConstraintAlgorithm;
class ReferenceNoseHooverChain;
class ReferenceMonteCarloBarostat;
class ReferenceNoseHooverDynamics;
class ReferenceVariableStochasticDynamics;
class ReferenceVariableVerletDynamics;
class ReferenceVerletDynamics;
class ReferenceCustomDynamics;

/**
 * This kernel is invoked at the beginning and end of force and energy computations.  It gives the
 * Platform a chance to clear buffers and do other initialization at the beginning, and to do any
 * necessary work at the end to determine the final results.
 */
class ReferenceCalcForcesAndEnergyKernel : public CalcForcesAndEnergyKernel {
public:
    ReferenceCalcForcesAndEnergyKernel(std::string name, const Platform& platform) : CalcForcesAndEnergyKernel(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     */
    void initialize(const System& system);
    /**
     * This is called at the beginning of each force/energy computation, before calcForcesAndEnergy() has been called on
     * any ForceImpl.
     *
     * @param context       the context in which to execute this kernel
     * @param includeForce  true if forces should be computed
     * @param includeEnergy true if potential energy should be computed
     * @param groups        a set of bit flags for which force groups to include
     */
    void beginComputation(ContextImpl& context, bool includeForce, bool includeEnergy, int groups);
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
    double finishComputation(ContextImpl& context, bool includeForce, bool includeEnergy, int groups, bool& valid);
private:
    std::vector<Vec3> savedForces;
};

/**
 * This kernel provides methods for setting and retrieving various state data: time, positions,
 * velocities, and forces.
 */
class ReferenceUpdateStateDataKernel : public UpdateStateDataKernel {
public:
    ReferenceUpdateStateDataKernel(std::string name, const Platform& platform, ReferencePlatform::PlatformData& data) : UpdateStateDataKernel(name, platform), data(data) {
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
    ReferencePlatform::PlatformData& data;
    std::vector<double> masses;
};

/**
 * This kernel modifies the positions of particles to enforce distance constraints.
 */
class ReferenceApplyConstraintsKernel : public ApplyConstraintsKernel {
public:
    ReferenceApplyConstraintsKernel(std::string name, const Platform& platform, ReferencePlatform::PlatformData& data) :
            ApplyConstraintsKernel(name, platform), data(data) {
    }
    ~ReferenceApplyConstraintsKernel();
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
    ReferencePlatform::PlatformData& data;
    std::vector<double> masses;
    std::vector<double> inverseMasses;
};

/**
 * This kernel recomputes the positions of virtual sites.
 */
class ReferenceVirtualSitesKernel : public VirtualSitesKernel {
public:
    ReferenceVirtualSitesKernel(std::string name, const Platform& platform) : VirtualSitesKernel(name, platform) {
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
};

/**
 * This kernel is invoked by HarmonicBondForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcHarmonicBondForceKernel : public CalcHarmonicBondForceKernel {
public:
    ReferenceCalcHarmonicBondForceKernel(std::string name, const Platform& platform) : CalcHarmonicBondForceKernel(name, platform) {
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
    int numBonds;
    std::vector<std::vector<int> >bondIndexArray;
    std::vector<std::vector<double> >bondParamArray;
    bool usePeriodic;
};

/**
 * This kernel is invoked by CustomBondForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcCustomBondForceKernel : public CalcCustomBondForceKernel {
public:
    ReferenceCalcCustomBondForceKernel(std::string name, const Platform& platform) : CalcCustomBondForceKernel(name, platform), ixn(NULL) {
    }
    ~ReferenceCalcCustomBondForceKernel();
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
     * @param firstBond  the index of the first bond whose parameters might have changed
     * @param lastBond   the index of the last bond whose parameters might have changed
     */
    void copyParametersToContext(ContextImpl& context, const CustomBondForce& force, int firstBond, int lastBond);
private:
    int numBonds;
    ReferenceCustomBondIxn* ixn;
    std::vector<std::vector<int> >bondIndexArray;
    std::vector<std::vector<double> >bondParamArray;
    Lepton::CompiledExpression energyExpression, forceExpression;
    std::vector<Lepton::CompiledExpression> energyParamDerivExpressions;
    std::vector<std::string> parameterNames, globalParameterNames, energyParamDerivNames;
    bool usePeriodic;
};

/**
 * This kernel is invoked by HarmonicAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcHarmonicAngleForceKernel : public CalcHarmonicAngleForceKernel {
public:
    ReferenceCalcHarmonicAngleForceKernel(std::string name, const Platform& platform) : CalcHarmonicAngleForceKernel(name, platform) {
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
     * @param firstAngle the index of the first bond whose parameters might have changed
     * @param lastAngle  the index of the last bond whose parameters might have changed
     */
    void copyParametersToContext(ContextImpl& context, const HarmonicAngleForce& force, int firstAngle, int lastAngle);
private:
    int numAngles;
    std::vector<std::vector<int> >angleIndexArray;
    std::vector<std::vector<double> >angleParamArray;
    bool usePeriodic;
};

/**
 * This kernel is invoked by CustomAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcCustomAngleForceKernel : public CalcCustomAngleForceKernel {
public:
    ReferenceCalcCustomAngleForceKernel(std::string name, const Platform& platform) : CalcCustomAngleForceKernel(name, platform), ixn(NULL) {
    }
    ~ReferenceCalcCustomAngleForceKernel();
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
     * @param firstAngle the index of the first bond whose parameters might have changed
     * @param lastAngle  the index of the last bond whose parameters might have changed
     */
    void copyParametersToContext(ContextImpl& context, const CustomAngleForce& force, int firstAngle, int lastAngle);
private:
    int numAngles;
    ReferenceCustomAngleIxn* ixn;
    std::vector<std::vector<int> >angleIndexArray;
    std::vector<std::vector<double> >angleParamArray;
    Lepton::CompiledExpression energyExpression, forceExpression;
    std::vector<Lepton::CompiledExpression> energyParamDerivExpressions;
    std::vector<std::string> parameterNames, globalParameterNames, energyParamDerivNames;
    bool usePeriodic;
};

/**
 * This kernel is invoked by PeriodicTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcPeriodicTorsionForceKernel : public CalcPeriodicTorsionForceKernel {
public:
    ReferenceCalcPeriodicTorsionForceKernel(std::string name, const Platform& platform) : CalcPeriodicTorsionForceKernel(name, platform) {
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
    int numTorsions;
    std::vector<std::vector<int> >torsionIndexArray;
    std::vector<std::vector<double> >torsionParamArray;
    bool usePeriodic;
};

/**
 * This kernel is invoked by RBTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcRBTorsionForceKernel : public CalcRBTorsionForceKernel {
public:
    ReferenceCalcRBTorsionForceKernel(std::string name, const Platform& platform) : CalcRBTorsionForceKernel(name, platform) {
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
    int numTorsions;
    std::vector<std::vector<int> >torsionIndexArray;
    std::vector<std::vector<double> >torsionParamArray;
    bool usePeriodic;
};

/**
 * This kernel is invoked by CMAPTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcCMAPTorsionForceKernel : public CalcCMAPTorsionForceKernel {
public:
    ReferenceCalcCMAPTorsionForceKernel(std::string name, const Platform& platform) : CalcCMAPTorsionForceKernel(name, platform) {
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
    std::vector<std::vector<std::vector<double> > > coeff;
    std::vector<int> torsionMaps;
    std::vector<std::vector<int> > torsionIndices;
    bool usePeriodic;
};

/**
 * This kernel is invoked by CustomTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcCustomTorsionForceKernel : public CalcCustomTorsionForceKernel {
public:
    ReferenceCalcCustomTorsionForceKernel(std::string name, const Platform& platform) : CalcCustomTorsionForceKernel(name, platform), ixn(NULL) {
    }
    ~ReferenceCalcCustomTorsionForceKernel();
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
    int numTorsions;
    ReferenceCustomTorsionIxn* ixn;
    std::vector<std::vector<int> >torsionIndexArray;
    std::vector<std::vector<double> >torsionParamArray;
    Lepton::CompiledExpression energyExpression, forceExpression;
    std::vector<Lepton::CompiledExpression> energyParamDerivExpressions;
    std::vector<std::string> parameterNames, globalParameterNames, energyParamDerivNames;
    bool usePeriodic;
};

/**
 * This kernel is invoked by NonbondedForce to calculate the forces acting on the system.
 */
class ReferenceCalcNonbondedForceKernel : public CalcNonbondedForceKernel {
public:
    ReferenceCalcNonbondedForceKernel(std::string name, const Platform& platform) : CalcNonbondedForceKernel(name, platform) {
    }
    ~ReferenceCalcNonbondedForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the NonbondedForce this kernel will be used for
     */
    void initialize(const System& system, const NonbondedForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @param includeReciprocal  true if reciprocal space interactions should be included
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy, bool includeDirect, bool includeReciprocal);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context        the context to copy parameters to
     * @param force          the NonbondedForce to copy the parameters from
     * @param firstParticle  the index of the first particle whose parameters might have changed
     * @param lastParticle   the index of the last particle whose parameters might have changed
     * @param firstException the index of the first exception whose parameters might have changed
     * @param lastException  the index of the last exception whose parameters might have changed
     */
    void copyParametersToContext(ContextImpl& context, const NonbondedForce& force, int firstParticle, int lastParticle, int firstException, int lastException);
    /**
     * Get the parameters being used for PME.
     * 
     * @param alpha   the separation parameter
     * @param nx      the number of grid points along the X axis
     * @param ny      the number of grid points along the Y axis
     * @param nz      the number of grid points along the Z axis
     */
    void getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const;
    /**
     * Get the dispersion parameters being used for the dispersion term in LJPME.
     *
     * @param alpha   the separation parameter
     * @param nx      the number of grid points along the X axis
     * @param ny      the number of grid points along the Y axis
     * @param nz      the number of grid points along the Z axis
     */
    void getLJPMEParameters(double& alpha, int& nx, int& ny, int& nz) const;
private:
    void computeParameters(ContextImpl& context);
    int numParticles, num14;
    std::vector<std::vector<int> >bonded14IndexArray;
    std::vector<std::vector<double> > particleParamArray, bonded14ParamArray;
    std::vector<std::array<double, 3> > baseParticleParams, baseExceptionParams;
    std::map<std::pair<std::string, int>, std::array<double, 3> > particleParamOffsets, exceptionParamOffsets;
    double nonbondedCutoff, switchingDistance, rfDielectric, ewaldAlpha, ewaldDispersionAlpha, dispersionCoefficient;
    int kmax[3], gridSize[3], dispersionGridSize[3];
    bool useSwitchingFunction, exceptionsArePeriodic;
    std::vector<std::set<int> > exclusions;
    NonbondedMethod nonbondedMethod;
    NeighborList* neighborList;
};

/**
 * This kernel is invoked by CustomNonbondedForce to calculate the forces acting on the system.
 */
class ReferenceCalcCustomNonbondedForceKernel : public CalcCustomNonbondedForceKernel {
public:
    ReferenceCalcCustomNonbondedForceKernel(std::string name, const Platform& platform) : CalcCustomNonbondedForceKernel(name, platform), forceCopy(NULL) {
    }
    ~ReferenceCalcCustomNonbondedForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomNonbondedForce this kernel will be used for
     */
    void initialize(const System& system, const CustomNonbondedForce& force);
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
     * @param force      the CustomNonbondedForce to copy the parameters from
     * @param firstParticle  the index of the first particle whose parameters might have changed
     * @param lastParticle   the index of the last particle whose parameters might have changed
     */
    void copyParametersToContext(ContextImpl& context, const CustomNonbondedForce& force, int firstParticle, int lastParticle);
private:
    void createExpressions(const CustomNonbondedForce& force);
    int numParticles;
    std::vector<std::vector<double> > particleParamArray;
    double nonbondedCutoff, switchingDistance, periodicBoxSize[3], longRangeCoefficient;
    bool useSwitchingFunction, hasInitializedLongRangeCorrection;
    CustomNonbondedForce* forceCopy;
    CustomNonbondedForceImpl::LongRangeCorrectionData longRangeCorrectionData;
    std::map<std::string, double> globalParamValues;
    std::vector<std::set<int> > exclusions;
    Lepton::CompiledExpression energyExpression, forceExpression;
    std::vector<Lepton::CompiledExpression> computedValueExpressions, energyParamDerivExpressions;
    std::vector<std::string> parameterNames, globalParameterNames, computedValueNames, energyParamDerivNames;
    std::vector<std::pair<std::set<int>, std::set<int> > > interactionGroups;
    std::vector<double> longRangeCoefficientDerivs;
    std::map<std::string, int> tabulatedFunctionUpdateCount;
    NonbondedMethod nonbondedMethod;
    NeighborList* neighborList;
};

/**
 * This kernel is invoked by GBSAOBCForce to calculate the forces acting on the system.
 */
class ReferenceCalcGBSAOBCForceKernel : public CalcGBSAOBCForceKernel {
public:
    ReferenceCalcGBSAOBCForceKernel(std::string name, const Platform& platform) : CalcGBSAOBCForceKernel(name, platform) {
    }
    ~ReferenceCalcGBSAOBCForceKernel();
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
    ReferenceObc* obc;
    std::vector<double> charges;
    bool isPeriodic;
};

/**
 * This kernel is invoked by CustomGBForce to calculate the forces acting on the system.
 */
class ReferenceCalcCustomGBForceKernel : public CalcCustomGBForceKernel {
public:
    ReferenceCalcCustomGBForceKernel(std::string name, const Platform& platform) : CalcCustomGBForceKernel(name, platform) {
    }
    ~ReferenceCalcCustomGBForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomGBForce this kernel will be used for
     */
    void initialize(const System& system, const CustomGBForce& force);
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
     * @param force      the CustomGBForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const CustomGBForce& force);
private:
    void createExpressions(const CustomGBForce& force);
    int numParticles;
    bool isPeriodic;
    std::vector<std::vector<double> > particleParamArray;
    double nonbondedCutoff;
    std::vector<std::set<int> > exclusions;
    std::vector<std::string> particleParameterNames, globalParameterNames, energyParamDerivNames, valueNames;
    std::vector<Lepton::CompiledExpression> valueExpressions;
    std::vector<std::vector<Lepton::CompiledExpression> > valueDerivExpressions;
    std::vector<std::vector<Lepton::CompiledExpression> > valueGradientExpressions;
    std::vector<std::vector<Lepton::CompiledExpression> > valueParamDerivExpressions;
    std::vector<OpenMM::CustomGBForce::ComputationType> valueTypes;
    std::vector<Lepton::CompiledExpression> energyExpressions;
    std::vector<std::vector<Lepton::CompiledExpression> > energyDerivExpressions;
    std::vector<std::vector<Lepton::CompiledExpression> > energyGradientExpressions;
    std::vector<std::vector<Lepton::CompiledExpression> > energyParamDerivExpressions;
    std::vector<OpenMM::CustomGBForce::ComputationType> energyTypes;
    std::map<std::string, int> tabulatedFunctionUpdateCount;
    NonbondedMethod nonbondedMethod;
    NeighborList* neighborList;
};

/**
 * This kernel is invoked by CustomExternalForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcCustomExternalForceKernel : public CalcCustomExternalForceKernel {
public:
    ReferenceCalcCustomExternalForceKernel(std::string name, const Platform& platform) : CalcCustomExternalForceKernel(name, platform), ixn(NULL) {
    }
    ~ReferenceCalcCustomExternalForceKernel();
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
    int numParticles;
    ReferenceCustomExternalIxn* ixn;
    std::vector<int> particles;
    std::vector<std::vector<double> > particleParamArray;
    Lepton::CompiledExpression energyExpression, forceExpressionX, forceExpressionY, forceExpressionZ;
    std::vector<std::string> parameterNames, globalParameterNames;
    Vec3* boxVectors;
};

/**
 * This kernel is invoked by CustomHbondForce to calculate the forces acting on the system.
 */
class ReferenceCalcCustomHbondForceKernel : public CalcCustomHbondForceKernel {
public:
    ReferenceCalcCustomHbondForceKernel(std::string name, const Platform& platform) : CalcCustomHbondForceKernel(name, platform), ixn(NULL) {
    }
    ~ReferenceCalcCustomHbondForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomHbondForce this kernel will be used for
     */
    void initialize(const System& system, const CustomHbondForce& force);
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
     * @param force      the CustomHbondForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const CustomHbondForce& force);
private:
    void createInteraction(const CustomHbondForce& force);
    int numDonors, numAcceptors, numParticles;
    bool isPeriodic;
    std::vector<std::vector<int> > donorParticles, acceptorParticles;
    std::vector<std::vector<double> > donorParamArray, acceptorParamArray;
    double nonbondedCutoff;
    ReferenceCustomHbondIxn* ixn;
    std::vector<std::set<int> > exclusions;
    std::vector<std::string> globalParameterNames;
    std::map<std::string, int> tabulatedFunctionUpdateCount;
};

/**
 * This kernel is invoked by CustomCentroidBondForce to calculate the forces acting on the system.
 */
class ReferenceCalcCustomCentroidBondForceKernel : public CalcCustomCentroidBondForceKernel {
public:
    ReferenceCalcCustomCentroidBondForceKernel(std::string name, const Platform& platform) : CalcCustomCentroidBondForceKernel(name, platform), ixn(NULL) {
    }
    ~ReferenceCalcCustomCentroidBondForceKernel();
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
    void createInteraction(const CustomCentroidBondForce& force);
    int numBonds, numParticles;
    std::vector<std::vector<int> > bondGroups;
    std::vector<std::vector<int> > groupAtoms;
    std::vector<std::vector<double> > normalizedWeights;
    std::vector<std::vector<double> > bondParamArray;
    ReferenceCustomCentroidBondIxn* ixn;
    std::vector<std::string> globalParameterNames, energyParamDerivNames;
    std::map<std::string, int> tabulatedFunctionUpdateCount;
    bool usePeriodic;
    Vec3* boxVectors;
};

/**
 * This kernel is invoked by CustomCompoundBondForce to calculate the forces acting on the system.
 */
class ReferenceCalcCustomCompoundBondForceKernel : public CalcCustomCompoundBondForceKernel {
public:
    ReferenceCalcCustomCompoundBondForceKernel(std::string name, const Platform& platform) : CalcCustomCompoundBondForceKernel(name, platform), ixn(NULL) {
    }
    ~ReferenceCalcCustomCompoundBondForceKernel();
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
    void createInteraction(const CustomCompoundBondForce& force);
    int numBonds;
    std::vector<std::vector<int> > bondParticles;
    std::vector<std::vector<double> > bondParamArray;
    ReferenceCustomCompoundBondIxn* ixn;
    std::vector<std::string> globalParameterNames, energyParamDerivNames;
    std::map<std::string, int> tabulatedFunctionUpdateCount;
    bool usePeriodic;
    Vec3* boxVectors;
};

/**
 * This kernel is invoked by CustomManyParticleForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcCustomManyParticleForceKernel : public CalcCustomManyParticleForceKernel {
public:
    ReferenceCalcCustomManyParticleForceKernel(std::string name, const Platform& platform) : CalcCustomManyParticleForceKernel(name, platform), ixn(NULL) {
    }
    ~ReferenceCalcCustomManyParticleForceKernel();
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
    int numParticles;
    double cutoffDistance;
    std::vector<std::vector<double> > particleParamArray;
    ReferenceCustomManyParticleIxn* ixn;
    std::vector<std::string> globalParameterNames;
    std::map<std::string, int> tabulatedFunctionUpdateCount;
    NonbondedMethod nonbondedMethod;
};

/**
 * This kernel is invoked by GayBerneForce to calculate the forces acting on the system.
 */
class ReferenceCalcGayBerneForceKernel : public CalcGayBerneForceKernel {
public:
    ReferenceCalcGayBerneForceKernel(std::string name, const Platform& platform) : CalcGayBerneForceKernel(name, platform), ixn(NULL) {
    }
    ~ReferenceCalcGayBerneForceKernel();
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
    ReferenceGayBerneForce* ixn;
};

/**
 * This kernel is invoked by CustomCVForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcCustomCVForceKernel : public CalcCustomCVForceKernel {
public:
    ReferenceCalcCustomCVForceKernel(std::string name, const Platform& platform) : CalcCustomCVForceKernel(name, platform), ixn(NULL) {
    }
    ~ReferenceCalcCustomCVForceKernel();
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
private:
    ReferenceCustomCVForce* ixn;
    std::vector<std::string> globalParameterNames, energyParamDerivNames;
};

/**
 * This kernel is invoked by RMSDForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcRMSDForceKernel : public CalcRMSDForceKernel {
public:
    ReferenceCalcRMSDForceKernel(std::string name, const Platform& platform) : CalcRMSDForceKernel(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the RMSDForce this kernel will be used for
     */
    void initialize(const System& system, const RMSDForce& force);
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
     * @param force      the RMSDForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const RMSDForce& force);
private:
    std::vector<Vec3> referencePos;
    std::vector<int> particles;
};

/**
 * This kernel is invoked by VerletIntegrator to take one time step.
 */
class ReferenceIntegrateVerletStepKernel : public IntegrateVerletStepKernel {
public:
    ReferenceIntegrateVerletStepKernel(std::string name, const Platform& platform, ReferencePlatform::PlatformData& data) : IntegrateVerletStepKernel(name, platform),
        data(data), dynamics(0) {
    }
    ~ReferenceIntegrateVerletStepKernel();
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
    ReferencePlatform::PlatformData& data;
    ReferenceVerletDynamics* dynamics;
    std::vector<double> masses;
    double prevStepSize;
}
;
/**
 * This kernel is invoked by NoseHooverIntegrator to take one time step.
 */
class ReferenceIntegrateNoseHooverStepKernel : public IntegrateNoseHooverStepKernel {
public:
    ReferenceIntegrateNoseHooverStepKernel(std::string name, const Platform& platform, ReferencePlatform::PlatformData& data) : IntegrateNoseHooverStepKernel(name, platform),
        data(data), dynamics(0) {
    }
    ~ReferenceIntegrateNoseHooverStepKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param integrator the NoseHooverIntegrator this kernel will be used for
     */
    void initialize(const System& system, const NoseHooverIntegrator& integrator);
    /**
     * Execute the kernel.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the VerletIntegrator this kernel is being used for
     * @param forcesAreValid a reference to the parent integrator's boolean for keeping
     *                       track of the validity of the current forces.
     */
    void execute(ContextImpl& context, const NoseHooverIntegrator& integrator, bool &forcesAreValid);
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the NoseHooverIntegrator this kernel is being used for
     */
    double computeKineticEnergy(ContextImpl& context, const NoseHooverIntegrator& integrator);
    /**
     * Execute the kernel that propagates the Nose Hoover chain and determines the velocity scale factor.
     * 
     * @param context  the context in which to execute this kernel
     * @param noseHooverChain the object describing the chain to be propagated.
     * @param kineticEnergy the {center of mass, relative} kineticEnergies of the particles being thermostated by this chain.
     * @param timeStep the time step used by the integrator.
     * @return the velocity scale factor to apply to the particles associated with this heat bath.
     */
    std::pair<double, double> propagateChain(ContextImpl& context, const NoseHooverChain &noseHooverChain, std::pair<double, double> kineticEnergy, double timeStep);
    /**
     * Execute the kernal that computes the total (kinetic + potential) heat bath energy.
     *
     * @param context the context in which to execute this kernel
     * @param noseHooverChain the chain whose energy is to be determined.
     * @return the total heat bath energy.
     */
    double computeHeatBathEnergy(ContextImpl& context, const NoseHooverChain &noseHooverChain);
    /**
     * Execute the kernel that computes the kinetic energy for a subset of atoms,
     * or the relative kinetic energy of Drude particles with respect to their parent atoms
     *
     * @param context the context in which to execute this kernel
     * @param noseHooverChain the chain whose energy is to be determined.
     * @param downloadValue whether the computed value should be downloaded and returned.
     */
    std::pair<double, double> computeMaskedKineticEnergy(ContextImpl& context, const NoseHooverChain &noseHooverChain, bool downloadValue);
    /**
     * Execute the kernel that scales the velocities of particles associated with a nose hoover chain
     *
     * @param context the context in which to execute this kernel
     * @param noseHooverChain the chain whose energy is to be determined.
     * @param scaleFactor the multiplicative factor by which {absolute, relative} velocities are scaled.
     */
    void scaleVelocities(ContextImpl& context, const NoseHooverChain &noseHooverChain, std::pair<double, double> scaleFactor);
    /**
     * Write the chain states to a checkpoint.
     */
    void createCheckpoint(ContextImpl& context, std::ostream& stream) const;
    /**
     * Load the chain states from a checkpoint.
     */
    void loadCheckpoint(ContextImpl& context, std::istream& stream);
    /**
     * Get the internal states of all chains.
     * 
     * @param context       the context for which to get the states
     * @param positions     element [i][j] contains the position of bead j for chain i
     * @param velocities    element [i][j] contains the velocity of bead j for chain i
     */
    void getChainStates(ContextImpl& context, std::vector<std::vector<double> >& positions, std::vector<std::vector<double> >& velocities) const;
    /**
     * Set the internal states of all chains.
     * 
     * @param context       the context for which to get the states
     * @param positions     element [i][j] contains the position of bead j for chain i
     * @param velocities    element [i][j] contains the velocity of bead j for chain i
     */
    void setChainStates(ContextImpl& context, const std::vector<std::vector<double> >& positions, const std::vector<std::vector<double> >& velocities);
private:
    ReferencePlatform::PlatformData& data;
    ReferenceNoseHooverChain* chainPropagator;
    ReferenceNoseHooverDynamics* dynamics;
    std::vector<double> masses;
    std::vector<std::vector<double> > chainPositions;
    std::vector<std::vector<double> > chainVelocities;
    double prevStepSize;
};

/**
 * This kernel is invoked by LangevinMiddleIntegrator to take one time step.
 */
class ReferenceIntegrateLangevinMiddleStepKernel : public IntegrateLangevinMiddleStepKernel {
public:
    ReferenceIntegrateLangevinMiddleStepKernel(std::string name, const Platform& platform, ReferencePlatform::PlatformData& data) : IntegrateLangevinMiddleStepKernel(name, platform),
        data(data), dynamics(0) {
    }
    ~ReferenceIntegrateLangevinMiddleStepKernel();
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
    ReferencePlatform::PlatformData& data;
    ReferenceLangevinMiddleDynamics* dynamics;
    std::vector<double> masses;
    double prevTemp, prevFriction, prevStepSize;
};

/**
 * This kernel is invoked by BrownianIntegrator to take one time step.
 */
class ReferenceIntegrateBrownianStepKernel : public IntegrateBrownianStepKernel {
public:
    ReferenceIntegrateBrownianStepKernel(std::string name, const Platform& platform, ReferencePlatform::PlatformData& data) : IntegrateBrownianStepKernel(name, platform),
        data(data), dynamics(0) {
    }
    ~ReferenceIntegrateBrownianStepKernel();
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
    ReferencePlatform::PlatformData& data;
    ReferenceBrownianDynamics* dynamics;
    std::vector<double> masses;
    double prevTemp, prevFriction, prevStepSize;
};

/**
 * This kernel is invoked by VariableLangevinIntegrator to take one time step.
 */
class ReferenceIntegrateVariableLangevinStepKernel : public IntegrateVariableLangevinStepKernel {
public:
    ReferenceIntegrateVariableLangevinStepKernel(std::string name, const Platform& platform, ReferencePlatform::PlatformData& data) : IntegrateVariableLangevinStepKernel(name, platform),
        data(data), dynamics(0) {
    }
    ~ReferenceIntegrateVariableLangevinStepKernel();
    /**
     * Initialize the kernel.
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
    ReferencePlatform::PlatformData& data;
    ReferenceVariableStochasticDynamics* dynamics;
    std::vector<double> masses;
    double prevTemp, prevFriction, prevErrorTol;
};

/**
 * This kernel is invoked by VariableVerletIntegrator to take one time step.
 */
class ReferenceIntegrateVariableVerletStepKernel : public IntegrateVariableVerletStepKernel {
public:
    ReferenceIntegrateVariableVerletStepKernel(std::string name, const Platform& platform, ReferencePlatform::PlatformData& data) : IntegrateVariableVerletStepKernel(name, platform),
        data(data), dynamics(0) {
    }
    ~ReferenceIntegrateVariableVerletStepKernel();
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
    ReferencePlatform::PlatformData& data;
    ReferenceVariableVerletDynamics* dynamics;
    std::vector<double> masses;
    double prevErrorTol;
};

/**
 * This kernel is invoked by CustomIntegrator to take one time step.
 */
class ReferenceIntegrateCustomStepKernel : public IntegrateCustomStepKernel {
public:
    ReferenceIntegrateCustomStepKernel(std::string name, const Platform& platform, ReferencePlatform::PlatformData& data) : IntegrateCustomStepKernel(name, platform),
        data(data), dynamics(0) {
    }
    ~ReferenceIntegrateCustomStepKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param integrator the CustomIntegrator this kernel will be used for
     */
    void initialize(const System& system, const CustomIntegrator& integrator);
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
    void execute(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid);
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
    double computeKineticEnergy(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid);
    /**
     * Get the values of all global variables.
     *
     * @param context   the context in which to execute this kernel
     * @param values    on exit, this contains the values
     */
    void getGlobalVariables(ContextImpl& context, std::vector<double>& values) const;
    /**
     * Set the values of all global variables.
     *
     * @param context   the context in which to execute this kernel
     * @param values    a vector containing the values
     */
    void setGlobalVariables(ContextImpl& context, const std::vector<double>& values);
    /**
     * Get the values of a per-DOF variable.
     *
     * @param context   the context in which to execute this kernel
     * @param variable  the index of the variable to get
     * @param values    on exit, this contains the values
     */
    void getPerDofVariable(ContextImpl& context, int variable, std::vector<Vec3>& values) const;
    /**
     * Set the values of a per-DOF variable.
     *
     * @param context   the context in which to execute this kernel
     * @param variable  the index of the variable to get
     * @param values    a vector containing the values
     */
    void setPerDofVariable(ContextImpl& context, int variable, const std::vector<Vec3>& values);
private:
    ReferencePlatform::PlatformData& data;
    ReferenceCustomDynamics* dynamics;
    std::vector<double> masses, globalValues;
    std::vector<std::vector<OpenMM::Vec3> > perDofValues; 
};

/**
 * This kernel is invoked by AndersenThermostat at the start of each time step to adjust the particle velocities.
 */
class ReferenceApplyAndersenThermostatKernel : public ApplyAndersenThermostatKernel {
public:
    ReferenceApplyAndersenThermostatKernel(std::string name, const Platform& platform) : ApplyAndersenThermostatKernel(name, platform), thermostat(0) {
    }
    ~ReferenceApplyAndersenThermostatKernel();
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
    ReferenceAndersenThermostat* thermostat;
    std::vector<std::vector<int> > particleGroups;
    std::vector<double> masses;
};

/**
 * This kernel is invoked by MonteCarloBarostat to adjust the periodic box volume
 */
class ReferenceApplyMonteCarloBarostatKernel : public ApplyMonteCarloBarostatKernel {
public:
    ReferenceApplyMonteCarloBarostatKernel(std::string name, const Platform& platform) : ApplyMonteCarloBarostatKernel(name, platform), barostat(NULL) {
    }
    ~ReferenceApplyMonteCarloBarostatKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param barostat   the MonteCarloBarostat this kernel will be used for
     * @param rigidMolecules  whether molecules should be kept rigid while scaling coordinates
     */
    void initialize(const System& system, const Force& barostat, bool rigidMolecules=true);
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
private:
    bool rigidMolecules;
    ReferenceMonteCarloBarostat* barostat;
};

/**
 * This kernel is invoked to remove center of mass motion from the system.
 */
class ReferenceRemoveCMMotionKernel : public RemoveCMMotionKernel {
public:
    ReferenceRemoveCMMotionKernel(std::string name, const Platform& platform, ReferencePlatform::PlatformData& data) : RemoveCMMotionKernel(name, platform), data(data) {
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
    ReferencePlatform::PlatformData& data;
    std::vector<double> masses;
    int frequency;
};

/**
 * This kernel is invoked by ATMForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcATMForceKernel : public CalcATMForceKernel {
public:
    ReferenceCalcATMForceKernel(std::string name, const Platform& platform) : CalcATMForceKernel(name, platform) {
    }
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
    void copyState(ContextImpl& context, ContextImpl& innerContext0, ContextImpl& innerContext1);
private:
    int numParticles;
    std::vector<Vec3> displ1;
    std::vector<Vec3> displ0;
};

/**
 * This kernel is invoked by CustomCPPForceImpl to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcCustomCPPForceKernel : public CalcCustomCPPForceKernel {
public:
    ReferenceCalcCustomCPPForceKernel(std::string name, const Platform& platform) : CalcCustomCPPForceKernel(name, platform), force(NULL) {
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
private:
    CustomCPPForceImpl* force;
    std::vector<Vec3> forces;
};

} // namespace OpenMM

#endif /*OPENMM_REFERENCEKERNELS_H_*/
