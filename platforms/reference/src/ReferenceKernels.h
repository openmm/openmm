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
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "kernels.h"
#include "SimTKUtilities/SimTKOpenMMRealType.h"

class CpuObc;
class ReferenceBrownianDynamics;
class ReferenceStochasticDynamics;
class ReferenceShakeAlgorithm;
class ReferenceVerletDynamics;

namespace OpenMM {

/**
 * This kernel is invoked by StandardMMForceField to calculate the forces acting on the system.
 */
class ReferenceCalcStandardMMForceFieldKernel : public CalcStandardMMForceFieldKernel {
public:
    ReferenceCalcStandardMMForceFieldKernel(std::string name, const Platform& platform) : CalcStandardMMForceFieldKernel(name, platform) {
    }
    ~ReferenceCalcStandardMMForceFieldKernel();
    /**
     * Initialize the kernel, setting up the values of all the force field parameters.
     * 
     * @param bondIndices               the two atoms connected by each bond term
     * @param bondParameters            the force parameters (length, k) for each bond term
     * @param angleIndices              the three atoms connected by each angle term
     * @param angleParameters           the force parameters (angle, k) for each angle term
     * @param periodicTorsionIndices    the four atoms connected by each periodic torsion term
     * @param periodicTorsionParameters the force parameters (k, phase, periodicity) for each periodic torsion term
     * @param rbTorsionIndices          the four atoms connected by each Ryckaert-Bellemans torsion term
     * @param rbTorsionParameters       the coefficients (in order of increasing powers) for each Ryckaert-Bellemans torsion term
     * @param bonded14Indices           each element contains the indices of two atoms whose nonbonded interactions should be reduced since
     *                                  they form a bonded 1-4 pair
     * @param lj14Scale                 the factor by which van der Waals interactions should be reduced for bonded 1-4 pairs
     * @param coulomb14Scale            the factor by which Coulomb interactions should be reduced for bonded 1-4 pairs
     * @param exclusions                the i'th element lists the indices of all atoms with which the i'th atom should not interact through
     *                                  nonbonded forces.  Bonded 1-4 pairs are also included in this list, since they should be omitted from
     *                                  the standard nonbonded calculation.
     * @param nonbondedParameters       the nonbonded force parameters (charge, sigma, epsilon) for each atom
     */
    void initialize(const std::vector<std::vector<int> >& bondIndices, const std::vector<std::vector<double> >& bondParameters,
            const std::vector<std::vector<int> >& angleIndices, const std::vector<std::vector<double> >& angleParameters,
            const std::vector<std::vector<int> >& periodicTorsionIndices, const std::vector<std::vector<double> >& periodicTorsionParameters,
            const std::vector<std::vector<int> >& rbTorsionIndices, const std::vector<std::vector<double> >& rbTorsionParameters,
            const std::vector<std::vector<int> >& bonded14Indices, double lj14Scale, double coulomb14Scale,
            const std::vector<std::set<int> >& exclusions, const std::vector<std::vector<double> >& nonbondedParameters);
    /**
     * Execute the kernel to calculate the forces.
     * 
     * @param positions   a Stream of type Double3 containing the position (x, y, z) of each atom
     * @param forces      a Stream of type Double3 containing the force (x, y, z) on each atom.  On entry, this contains the forces that
     *                    have been calculated so far.  The kernel should add its own forces to the values already in the stream.
     */
    void executeForces(const Stream& positions, Stream& forces);
    /**
     * Execute the kernel to calculate the energy.
     * 
     * @param positions   a Stream of type Double3 containing the position (x, y, z) of each atom
     * @return the potential energy due to the StandardMMForceField
     */
    double executeEnergy(const Stream& positions);
private:
    int numAtoms, numBonds, numAngles, numPeriodicTorsions, numRBTorsions, num14;
    int **bondIndexArray, **angleIndexArray, **periodicTorsionIndexArray, **rbTorsionIndexArray, **exclusionArray, **bonded14IndexArray;
    RealOpenMM **bondParamArray, **angleParamArray, **periodicTorsionParamArray, **rbTorsionParamArray, **atomParamArray, **bonded14ParamArray;
};

/**
 * This kernel is invoked by GBSAOBCForceField to calculate the forces acting on the system.
 */
class ReferenceCalcGBSAOBCForceFieldKernel : public CalcGBSAOBCForceFieldKernel {
public:
    ReferenceCalcGBSAOBCForceFieldKernel(std::string name, const Platform& platform) : CalcGBSAOBCForceFieldKernel(name, platform) {
    }
    ~ReferenceCalcGBSAOBCForceFieldKernel();
    /**
     * Initialize the kernel, setting up the values of all the force field parameters.
     * 
     * @param atomParameters      the force parameters (charge, atomic radius, scaling factor) for each atom
     * @param solventDielectric   the dielectric constant of the solvent
     * @param soluteDielectric    the dielectric constant of the solute
     */
    void initialize(const std::vector<std::vector<double> >& atomParameters, double solventDielectric, double soluteDielectric);
    /**
     * Execute the kernel to calculate the forces.
     * 
     * @param positions   a Stream of type Double3 containing the position (x, y, z) of each atom
     * @param forces      a Stream of type Double3 containing the force (x, y, z) on each atom.  On entry, this contains the forces that
     *                    have been calculated so far.  The kernel should add its own forces to the values already in the stream.
     */
    void executeForces(const Stream& positions, Stream& forces);
    /**
     * Execute the kernel to calculate the energy.
     * 
     * @param positions   a Stream of type Double3 containing the position (x, y, z) of each atom
     * @return the potential energy due to the GBSAOBCForceField
     */
    double executeEnergy(const Stream& positions);
private:
    CpuObc* obc;
    std::vector<RealOpenMM> charges;
};

/**
 * This kernel is invoked by VerletIntegrator to take one time step.
 */
class ReferenceIntegrateVerletStepKernel : public IntegrateVerletStepKernel {
public:
    ReferenceIntegrateVerletStepKernel(std::string name, const Platform& platform) : IntegrateVerletStepKernel(name, platform),
        dynamics(0), shake(0), masses(0), shakeParameters(0), constraintIndices(0) {
    }
    ~ReferenceIntegrateVerletStepKernel();
    /**
     * Initialize the kernel, setting up all parameters related to integrator.
     * 
     * @param masses             the mass of each atom
     * @param constraintIndices  each element contains the indices of two atoms whose distance should be constrained
     * @param constraintLengths  the required distance between each pair of constrained atoms
     */
    void initialize(const std::vector<double>& masses, const std::vector<std::vector<int> >& constraintIndices,
            const std::vector<double>& constraintLengths);
    /**
     * Execute the kernel.
     * 
     * @param positions          a Stream of type Double3 containing the position (x, y, z) of each atom
     * @param velocities         a Stream of type Double3 containing the velocity (x, y, z) of each atom
     * @param forces             a Stream of type Double3 containing the force (x, y, z) on each atom
     * @param stepSize           the integration step size
     */
    void execute(Stream& positions, Stream& velocities, const Stream& forces, double stepSize);
private:
    ReferenceVerletDynamics* dynamics;
    ReferenceShakeAlgorithm* shake;
    RealOpenMM* masses;
    RealOpenMM** shakeParameters;
    int** constraintIndices;
    int numConstraints;
    double prevStepSize;
};

/**
 * This kernel is invoked by LangevinIntegrator to take one time step.
 */
class ReferenceIntegrateLangevinStepKernel : public IntegrateLangevinStepKernel {
public:
    ReferenceIntegrateLangevinStepKernel(std::string name, const Platform& platform) : IntegrateLangevinStepKernel(name, platform),
        dynamics(0), shake(0), masses(0), shakeParameters(0), constraintIndices(0) {
    }
    ~ReferenceIntegrateLangevinStepKernel();
    /**
     * Initialize the kernel, setting up all parameters related to integrator.
     * 
     * @param masses             the mass of each atom
     * @param constraintIndices  each element contains the indices of two atoms whose distance should be constrained
     * @param constraintLengths  the required distance between each pair of constrained atoms
     */
    void initialize(const std::vector<double>& masses, const std::vector<std::vector<int> >& constraintIndices,
            const std::vector<double>& constraintLengths);
    /**
     * Execute the kernel.
     * 
     * @param positions          a Stream of type Double3 containing the position (x, y, z) of each atom
     * @param velocities         a Stream of type Double3 containing the velocity (x, y, z) of each atom
     * @param forces             a Stream of type Double3 containing the force (x, y, z) on each atom
     * @param temperature        the temperature of the heat bath
     * @param friction           the friction coefficient coupling the system to the heat bath
     * @param stepSize           the integration step size
     */
    void execute(Stream& positions, Stream& velocities, const Stream& forces, double temperature, double friction, double stepSize);
private:
    ReferenceStochasticDynamics* dynamics;
    ReferenceShakeAlgorithm* shake;
    RealOpenMM* masses;
    RealOpenMM** shakeParameters;
    int** constraintIndices;
    int numConstraints;
    double prevTemp, prevFriction, prevStepSize;
};

/**
 * This kernel is invoked by BrownianIntegrator to take one time step.
 */
class ReferenceIntegrateBrownianStepKernel : public IntegrateBrownianStepKernel {
public:
    ReferenceIntegrateBrownianStepKernel(std::string name, const Platform& platform) : IntegrateBrownianStepKernel(name, platform),
        dynamics(0), shake(0), masses(0), shakeParameters(0), constraintIndices(0) {
    }
    ~ReferenceIntegrateBrownianStepKernel();
    /**
     * Initialize the kernel, setting up all parameters related to integrator.
     * 
     * @param masses             the mass of each atom
     * @param constraintIndices  each element contains the indices of two atoms whose distance should be constrained
     * @param constraintLengths  the required distance between each pair of constrained atoms
     */
    void initialize(const std::vector<double>& masses, const std::vector<std::vector<int> >& constraintIndices,
            const std::vector<double>& constraintLengths);
    /**
     * Execute the kernel.
     * 
     * @param positions          a Stream of type Double3 containing the position (x, y, z) of each atom
     * @param velocities         a Stream of type Double3 containing the velocity (x, y, z) of each atom
     * @param forces             a Stream of type Double3 containing the force (x, y, z) on each atom
     * @param temperature        the temperature of the heat bath
     * @param friction           the friction coefficient coupling the system to the heat bath
     * @param stepSize           the integration step size
     */
    void execute(Stream& positions, Stream& velocities, const Stream& forces, double temperature, double friction, double stepSize);
private:
    ReferenceBrownianDynamics* dynamics;
    ReferenceShakeAlgorithm* shake;
    RealOpenMM* masses;
    RealOpenMM** shakeParameters;
    int** constraintIndices;
    int numConstraints;
    double prevTemp, prevFriction, prevStepSize;
};

/**
 * This kernel is invoked by AndersenThermostat at the start of each time step to adjust the atom velocities.
 */
class ReferenceApplyAndersenThermostatKernel : public ApplyAndersenThermostatKernel {
public:
    ReferenceApplyAndersenThermostatKernel(std::string name, const Platform& platform) : ApplyAndersenThermostatKernel(name, platform) {
    }
    /**
     * Initialize the kernel, setting up the values of unchanging parameters.
     * 
     * @param masses the mass of each atom
     */
    void initialize(const std::vector<double>& masses);
    /**
     * Execute the kernel.
     * 
     * @param velocities         a Stream of type Double3 containing the velocity (x, y, z) of each atom
     * @param temperature        the temperature of the heat bath
     * @param collisionFrequency the frequency at which atom collide with particles in the heat bath
     * @param stepSize           the integration step size
     */
    void execute(Stream& velocities, double temperature, double collisionFrequency, double stepSize);
};

/**
 * This kernel is invoked to calculate the kinetic energy of the system.
 */
class ReferenceCalcKineticEnergyKernel : public CalcKineticEnergyKernel {
public:
    ReferenceCalcKineticEnergyKernel(std::string name, const Platform& platform) : CalcKineticEnergyKernel(name, platform) {
    }
    /**
     * Initialize the kernel, setting up the atomic masses.
     * 
     * @param masses the mass of each atom
     */
    void initialize(const std::vector<double>& masses);
    /**
     * Execute the kernel.
     * 
     * @param velocities a Stream of type Double3 containing the velocity (x, y, z) of each atom
     * @return the kinetic energy of the system
     */
    double execute(const Stream& velocities);
private:
    std::vector<double> masses;
};

} // namespace OpenMM

#endif /*OPENMM_REFERENCEKERNELS_H_*/
