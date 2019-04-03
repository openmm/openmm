#ifndef AMOEBA_OPENMM_REFERENCE_KERNELS_H_
#define AMOEBA_OPENMM_REFERENCE_KERNELS_H_

/* -------------------------------------------------------------------------- *
 *                              OpenMMAmoeba                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2018 Stanford University and the Authors.      *
 * Authors:                                                                   *
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

#include "openmm/System.h"
#include "openmm/amoebaKernels.h"
#include "openmm/AmoebaMultipoleForce.h"
#include "openmm/HippoNonbondedForce.h"
#include "AmoebaReferenceMultipoleForce.h"
#include "AmoebaReferenceHippoNonbondedForce.h"
#include "ReferenceNeighborList.h"
#include "SimTKOpenMMRealType.h"

namespace OpenMM {

/**
 * This kernel is invoked by AmoebaBondForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcAmoebaBondForceKernel : public CalcAmoebaBondForceKernel {
public:
    ReferenceCalcAmoebaBondForceKernel(std::string name, 
                                               const Platform& platform,
                                               const System& system);
    ~ReferenceCalcAmoebaBondForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the AmoebaBondForce this kernel will be used for
     */
    void initialize(const System& system, const AmoebaBondForce& force);
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
     * @param force      the AmoebaBondForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const AmoebaBondForce& force);
private:
    int numBonds;
    std::vector<int>   particle1;
    std::vector<int>   particle2;
    std::vector<double> length;
    std::vector<double> kQuadratic;
    double globalBondCubic;
    double globalBondQuartic;
    const System& system;
    bool usePeriodic;
};

/**
 * This kernel is invoked by AmoebaAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcAmoebaAngleForceKernel : public CalcAmoebaAngleForceKernel {
public:
    ReferenceCalcAmoebaAngleForceKernel(std::string name, const Platform& platform, const System& system);
    ~ReferenceCalcAmoebaAngleForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the AmoebaAngleForce this kernel will be used for
     */
    void initialize(const System& system, const AmoebaAngleForce& force);
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
     * @param force      the AmoebaAngleForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const AmoebaAngleForce& force);
private:
    int numAngles;
    std::vector<int>   particle1;
    std::vector<int>   particle2;
    std::vector<int>   particle3;
    std::vector<double> angle;
    std::vector<double> kQuadratic;
    double globalAngleCubic;
    double globalAngleQuartic;
    double globalAnglePentic;
    double globalAngleSextic;
    const System& system;
    bool usePeriodic;
};

/**
 * This kernel is invoked by AmoebaInPlaneAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcAmoebaInPlaneAngleForceKernel : public CalcAmoebaInPlaneAngleForceKernel {
public:
    ReferenceCalcAmoebaInPlaneAngleForceKernel(std::string name, const Platform& platform, const System& system);
    ~ReferenceCalcAmoebaInPlaneAngleForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the AmoebaInPlaneAngleForce this kernel will be used for
     */
    void initialize(const System& system, const AmoebaInPlaneAngleForce& force);
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
     * @param force      the AmoebaInPlaneAngleForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const AmoebaInPlaneAngleForce& force);
private:
    int numAngles;
    std::vector<int>   particle1;
    std::vector<int>   particle2;
    std::vector<int>   particle3;
    std::vector<int>   particle4;
    std::vector<double> angle;
    std::vector<double> kQuadratic;
    double globalInPlaneAngleCubic;
    double globalInPlaneAngleQuartic;
    double globalInPlaneAnglePentic;
    double globalInPlaneAngleSextic;
    const System& system;
    bool usePeriodic;
};

/**
 * This kernel is invoked by AmoebaPiTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcAmoebaPiTorsionForceKernel : public CalcAmoebaPiTorsionForceKernel {
public:
    ReferenceCalcAmoebaPiTorsionForceKernel(std::string name, const Platform& platform, const System& system);
    ~ReferenceCalcAmoebaPiTorsionForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the AmoebaPiTorsionForce this kernel will be used for
     */
    void initialize(const System& system, const AmoebaPiTorsionForce& force);
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
     * @param force      the AmoebaPiTorsionForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const AmoebaPiTorsionForce& force);
private:
    int numPiTorsions;
    std::vector<int>   particle1;
    std::vector<int>   particle2;
    std::vector<int>   particle3;
    std::vector<int>   particle4;
    std::vector<int>   particle5;
    std::vector<int>   particle6;
    std::vector<double> kTorsion;
    const System& system;
    bool usePeriodic;
};

/**
 * This kernel is invoked by AmoebaStretchBendForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcAmoebaStretchBendForceKernel : public CalcAmoebaStretchBendForceKernel {
public:
    ReferenceCalcAmoebaStretchBendForceKernel(std::string name, const Platform& platform, const System& system);
    ~ReferenceCalcAmoebaStretchBendForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the AmoebaStretchBendForce this kernel will be used for
     */
    void initialize(const System& system, const AmoebaStretchBendForce& force);
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
     * @param force      the AmoebaStretchBendForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const AmoebaStretchBendForce& force);
private:
    int numStretchBends;
    std::vector<int>   particle1;
    std::vector<int>   particle2;
    std::vector<int>   particle3;
    std::vector<double> lengthABParameters;
    std::vector<double> lengthCBParameters;
    std::vector<double> angleParameters;
    std::vector<double> k1Parameters;
    std::vector<double> k2Parameters;
    const System& system;
    bool usePeriodic;
};

/**
 * This kernel is invoked by AmoebaOutOfPlaneBendForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcAmoebaOutOfPlaneBendForceKernel : public CalcAmoebaOutOfPlaneBendForceKernel {
public:
    ReferenceCalcAmoebaOutOfPlaneBendForceKernel(std::string name, const Platform& platform, const System& system);
    ~ReferenceCalcAmoebaOutOfPlaneBendForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the AmoebaOutOfPlaneBendForce this kernel will be used for
     */
    void initialize(const System& system, const AmoebaOutOfPlaneBendForce& force);
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
     * @param force      the AmoebaOutOfPlaneBendForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const AmoebaOutOfPlaneBendForce& force);
private:
    int numOutOfPlaneBends;
    std::vector<int>   particle1;
    std::vector<int>   particle2;
    std::vector<int>   particle3;
    std::vector<int>   particle4;
    std::vector<double> kParameters;
    double globalOutOfPlaneBendAngleCubic;
    double globalOutOfPlaneBendAngleQuartic;
    double globalOutOfPlaneBendAnglePentic;
    double globalOutOfPlaneBendAngleSextic;
    const System& system;
    bool usePeriodic;
};

/**
 * This kernel is invoked by AmoebaTorsionTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcAmoebaTorsionTorsionForceKernel : public CalcAmoebaTorsionTorsionForceKernel {
public:
    ReferenceCalcAmoebaTorsionTorsionForceKernel(std::string name, const Platform& platform, const System& system);
    ~ReferenceCalcAmoebaTorsionTorsionForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the AmoebaTorsionTorsionForce this kernel will be used for
     */
    void initialize(const System& system, const AmoebaTorsionTorsionForce& force);
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
    int numTorsionTorsions;
    std::vector<int>   particle1;
    std::vector<int>   particle2;
    std::vector<int>   particle3;
    std::vector<int>   particle4;
    std::vector<int>   particle5;
    std::vector<int>   chiralCheckAtom;
    std::vector<int>   gridIndices;

    int numTorsionTorsionGrids;
    std::vector< std::vector< std::vector< std::vector<double> > > > torsionTorsionGrids;

    const System& system;
    bool usePeriodic;
};

/**
 * This kernel is invoked by AmoebaMultipoleForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcAmoebaMultipoleForceKernel : public CalcAmoebaMultipoleForceKernel {
public:
    ReferenceCalcAmoebaMultipoleForceKernel(std::string name, const Platform& platform, const System& system);
    ~ReferenceCalcAmoebaMultipoleForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the AmoebaMultipoleForce this kernel will be used for
     */
    void initialize(const System& system, const AmoebaMultipoleForce& force);
    /**
     * Setup for AmoebaReferenceMultipoleForce instance. 
     *
     * @param context        the current context
     *
     * @return pointer to initialized instance of AmoebaReferenceMultipoleForce
     */
    AmoebaReferenceMultipoleForce* setupAmoebaReferenceMultipoleForce(ContextImpl& context);
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
     * Get the induced dipole moments of all particles.
     * 
     * @param context    the Context for which to get the induced dipoles
     * @param dipoles    the induced dipole moment of particle i is stored into the i'th element
     */
    void getInducedDipoles(ContextImpl& context, std::vector<Vec3>& dipoles);
    /**
     * Get the fixed dipole moments of all particles in the global reference frame.
     * 
     * @param context    the Context for which to get the fixed dipoles
     * @param dipoles    the fixed dipole moment of particle i is stored into the i'th element
     */
    void getLabFramePermanentDipoles(ContextImpl& context, std::vector<Vec3>& dipoles);
    /**
     * Get the total dipole moments of all particles in the global reference frame.
     * 
     * @param context    the Context for which to get the fixed dipoles
     * @param dipoles    the fixed dipole moment of particle i is stored into the i'th element
     */
    void getTotalDipoles(ContextImpl& context, std::vector<Vec3>& dipoles);
    /** 
     * Calculate the electrostatic potential given vector of grid coordinates.
     *
     * @param context                      context
     * @param inputGrid                    input grid coordinates
     * @param outputElectrostaticPotential output potential 
     */
    void getElectrostaticPotential(ContextImpl& context, const std::vector< Vec3 >& inputGrid,
                                   std::vector< double >& outputElectrostaticPotential);

    /**
     * Get the system multipole moments.
     *
     * @param context                context 
     * @param outputMultipoleMoments vector of multipole moments:
                                     (charge,
                                      dipole_x, dipole_y, dipole_z,
                                      quadrupole_xx, quadrupole_xy, quadrupole_xz,
                                      quadrupole_yx, quadrupole_yy, quadrupole_yz,
                                      quadrupole_zx, quadrupole_zy, quadrupole_zz)
     */
    void getSystemMultipoleMoments(ContextImpl& context, std::vector< double >& outputMultipoleMoments);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the AmoebaMultipoleForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const AmoebaMultipoleForce& force);
    /**
     * Get the parameters being used for PME.
     * 
     * @param alpha   the separation parameter
     * @param nx      the number of grid points along the X axis
     * @param ny      the number of grid points along the Y axis
     * @param nz      the number of grid points along the Z axis
     */
    void getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const;

private:

    int numMultipoles;
    AmoebaMultipoleForce::NonbondedMethod nonbondedMethod;
    AmoebaMultipoleForce::PolarizationType polarizationType;
    std::vector<double> charges;
    std::vector<double> dipoles;
    std::vector<double> quadrupoles;
    std::vector<double> tholes;
    std::vector<double> dampingFactors;
    std::vector<double> polarity;
    std::vector<int>   axisTypes;
    std::vector<int>   multipoleAtomZs;
    std::vector<int>   multipoleAtomXs;
    std::vector<int>   multipoleAtomYs;
    std::vector< std::vector< std::vector<int> > > multipoleAtomCovalentInfo;

    int mutualInducedMaxIterations;
    double mutualInducedTargetEpsilon;
    std::vector<double> extrapolationCoefficients;

    bool usePme;
    double alphaEwald;
    double cutoffDistance;
    std::vector<int> pmeGridDimension;

    const System& system;
};

/**
 * This kernel is invoked to calculate the vdw forces acting on the system and the energy of the system.
 */
class ReferenceCalcAmoebaVdwForceKernel : public CalcAmoebaVdwForceKernel {
public:
    ReferenceCalcAmoebaVdwForceKernel(std::string name, const Platform& platform, const System& system);
    ~ReferenceCalcAmoebaVdwForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the AmoebaVdwForce this kernel will be used for
     */
    void initialize(const System& system, const AmoebaVdwForce& force);
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
     * @param force      the AmoebaVdwForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const AmoebaVdwForce& force);
private:
    int numParticles;
    int useCutoff;
    int usePBC;
    double cutoff;
    double dispersionCoefficient;
    std::vector<int> indexIVs;
    std::vector< std::set<int> > allExclusions;
    std::vector<double> sigmas;
    std::vector<double> epsilons;
    std::vector<double> reductions;
    std::string sigmaCombiningRule;
    std::string epsilonCombiningRule;
    const System& system;
    NeighborList* neighborList;
};

/**
 * This kernel is invoked to calculate the WCA dispersion forces acting on the system and the energy of the system.
 */
class ReferenceCalcAmoebaWcaDispersionForceKernel : public CalcAmoebaWcaDispersionForceKernel {
public:
    ReferenceCalcAmoebaWcaDispersionForceKernel(std::string name, const Platform& platform, const System& system);
    ~ReferenceCalcAmoebaWcaDispersionForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the AmoebaMultipoleForce this kernel will be used for
     */
    void initialize(const System& system, const AmoebaWcaDispersionForce& force);
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
     * @param force      the AmoebaWcaDispersionForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const AmoebaWcaDispersionForce& force);
private:

    int numParticles;
    std::vector<double> radii;
    std::vector<double> epsilons;
    double epso; 
    double epsh; 
    double rmino; 
    double rminh; 
    double awater; 
    double shctd; 
    double dispoff;
    double slevy;
    double totalMaximumDispersionEnergy;
    const System& system;
};

/**
 * This kernel is invoked to calculate the Gerneralized Kirkwood forces acting on the system and the energy of the system.
 */
class ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel : public CalcAmoebaGeneralizedKirkwoodForceKernel {
public:
    ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel(std::string name, const Platform& platform, const System& system);
    ~ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the AmoebaMultipoleForce this kernel will be used for
     */
    void initialize(const System& system, const AmoebaGeneralizedKirkwoodForce& force);
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
     *  Get the 'include cavity term' flag.
     *
     *  @return includeCavityTerm
     */
    int getIncludeCavityTerm() const;

    /**
     *  Get the number of particles.
     *
     *  @return number of particles
     */
    int getNumParticles() const;

    /**
     *  Get Direct Polarization flag.
     *
     *  @return directPolarization
     *
     */
    int getDirectPolarization() const;

    /**
     *  Get the solute dielectric.
     *
     *  @return soluteDielectric
     *
     */
    double getSoluteDielectric() const;

    /**
     *  Get the solvent dielectric.
     *
     *  @return solventDielectric
     *
     */
    double getSolventDielectric() const;

    /**
     *  Get the dielectric offset.
     *
     *  @return dielectricOffset
     *
     */
    double getDielectricOffset() const;

    /**
     *  Get the probe radius.
     *
     *  @return probeRadius
     *
     */
    double getProbeRadius() const;

    /**
     *  Get the surface area factor.
     *
     *  @return surfaceAreaFactor
     *
     */
    double getSurfaceAreaFactor() const;

    /**
     *  Get the vector of particle radii.
     *
     *  @param atomicRadii vector of atomic radii
     *
     */
    void getAtomicRadii(std::vector<double>& atomicRadii) const;

    /**
     *  Get the vector of scale factors.
     *
     *  @param scaleFactors vector of scale factors
     *
     */
    void getScaleFactors(std::vector<double>& scaleFactors) const;

    /**
     *  Get the vector of charges.
     *
     *  @param charges vector of charges
     *
     */
    void getCharges(std::vector<double>& charges) const;

    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the AmoebaGeneralizedKirkwoodForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const AmoebaGeneralizedKirkwoodForce& force);

private:

    int numParticles;
    std::vector<double> atomicRadii;
    std::vector<double> scaleFactors;
    std::vector<double> charges;
    double soluteDielectric;
    double solventDielectric;
    double dielectricOffset;
    double probeRadius;
    double surfaceAreaFactor;
    int includeCavityTerm;
    int directPolarization;
    const System& system;
};

/**
 * This kernel is invoked by HippoNonbondedForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcHippoNonbondedForceKernel : public CalcHippoNonbondedForceKernel {
public:
    ReferenceCalcHippoNonbondedForceKernel(std::string name, const Platform& platform, const System& system);
    ~ReferenceCalcHippoNonbondedForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the HippoNonbondedForce this kernel will be used for
     */
    void initialize(const System& system, const HippoNonbondedForce& force);
    /**
     * Setup for AmoebaReferenceHippoNonbondedForce instance. 
     *
     * @param context        the current context
     *
     * @return pointer to initialized instance of AmoebaReferenceHippoNonbondedForce
     */
    void setupAmoebaReferenceHippoNonbondedForce(ContextImpl& context);
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
     * Get the induced dipole moments of all particles.
     * 
     * @param context    the Context for which to get the induced dipoles
     * @param dipoles    the induced dipole moment of particle i is stored into the i'th element
     */
    void getInducedDipoles(ContextImpl& context, std::vector<Vec3>& dipoles);
    /**
     * Get the fixed dipole moments of all particles in the global reference frame.
     * 
     * @param context    the Context for which to get the fixed dipoles
     * @param dipoles    the fixed dipole moment of particle i is stored into the i'th element
     */
    void getLabFramePermanentDipoles(ContextImpl& context, std::vector<Vec3>& dipoles);
    /**
     * Get the total dipole moments of all particles in the global reference frame.
     * 
     * @param context    the Context for which to get the fixed dipoles
     * @param dipoles    the fixed dipole moment of particle i is stored into the i'th element
     */
    void getTotalDipoles(ContextImpl& context, std::vector<Vec3>& dipoles);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the HippoNonbondedForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const HippoNonbondedForce& force);
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
     * Get the parameters being used for dispersion PME.
     * 
     * @param alpha   the separation parameter
     * @param nx      the number of grid points along the X axis
     * @param ny      the number of grid points along the Y axis
     * @param nz      the number of grid points along the Z axis
     */
    void getDPMEParameters(double& alpha, int& nx, int& ny, int& nz) const;

private:

    AmoebaReferenceHippoNonbondedForce* ixn;
    int numParticles;
};

} // namespace OpenMM

#endif /*AMOEBA_OPENMM_REFERENCE_KERNELS_H*/
