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
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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
#include "AmoebaReferenceMultipoleForce.h"
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
    std::vector<RealOpenMM> length;
    std::vector<RealOpenMM> kQuadratic;
    RealOpenMM globalBondCubic;
    RealOpenMM globalBondQuartic;
    const System& system;
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
    std::vector<RealOpenMM> angle;
    std::vector<RealOpenMM> kQuadratic;
    RealOpenMM globalAngleCubic;
    RealOpenMM globalAngleQuartic;
    RealOpenMM globalAnglePentic;
    RealOpenMM globalAngleSextic;
    const System& system;
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
    std::vector<RealOpenMM> angle;
    std::vector<RealOpenMM> kQuadratic;
    RealOpenMM globalInPlaneAngleCubic;
    RealOpenMM globalInPlaneAngleQuartic;
    RealOpenMM globalInPlaneAnglePentic;
    RealOpenMM globalInPlaneAngleSextic;
    const System& system;
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
    std::vector<RealOpenMM> kTorsion;
    const System& system;
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
    std::vector<RealOpenMM> lengthABParameters;
    std::vector<RealOpenMM> lengthCBParameters;
    std::vector<RealOpenMM> angleParameters;
    std::vector<RealOpenMM> kParameters;
    const System& system;
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
    std::vector<RealOpenMM> kParameters;
    RealOpenMM globalOutOfPlaneBendAngleCubic;
    RealOpenMM globalOutOfPlaneBendAngleQuartic;
    RealOpenMM globalOutOfPlaneBendAnglePentic;
    RealOpenMM globalOutOfPlaneBendAngleSextic;
    const System& system;
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
    std::vector< std::vector< std::vector< std::vector<RealOpenMM> > > > torsionTorsionGrids;

    const System& system;
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
    AmoebaReferenceMultipoleForce* setupAmoebaReferenceMultipoleForce(ContextImpl& context );
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
     * Calculate the electrostatic potential given vector of grid coordinates.
     *
     * @param context                      context
     * @param inputGrid                    input grid coordinates
     * @param outputElectrostaticPotential output potential 
     */
    void getElectrostaticPotential(ContextImpl& context, const std::vector< Vec3 >& inputGrid,
                                   std::vector< double >& outputElectrostaticPotential );

    /**
     * Get the system multipole moments.
     *
     * @param context                context 
     * @param outputMultipoleMonents vector of multipole moments:
                                     (charge,
                                      dipole_x, dipole_y, dipole_z,
                                      quadrupole_xx, quadrupole_xy, quadrupole_xz,
                                      quadrupole_yx, quadrupole_yy, quadrupole_yz,
                                      quadrupole_zx, quadrupole_zy, quadrupole_zz )
     */
    void getSystemMultipoleMoments(ContextImpl& context, std::vector< double >& outputMultipoleMoments);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the AmoebaMultipoleForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const AmoebaMultipoleForce& force);

private:

    int numMultipoles;
    AmoebaMultipoleForce::NonbondedMethod nonbondedMethod;
    AmoebaMultipoleForce::PolarizationType polarizationType;
    std::vector<RealOpenMM> charges;
    std::vector<RealOpenMM> dipoles;
    std::vector<RealOpenMM> quadrupoles;
    std::vector<RealOpenMM> tholes;
    std::vector<RealOpenMM> dampingFactors;
    std::vector<RealOpenMM> polarity;
    std::vector<int>   axisTypes;
    std::vector<int>   multipoleAtomZs;
    std::vector<int>   multipoleAtomXs;
    std::vector<int>   multipoleAtomYs;
    std::vector< std::vector< std::vector<int> > > multipoleAtomCovalentInfo;

    int mutualInducedMaxIterations;
    RealOpenMM mutualInducedTargetEpsilon;

    bool usePme;
    RealOpenMM alphaEwald;
    RealOpenMM cutoffDistance;
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
    std::vector<RealOpenMM> sigmas;
    std::vector<RealOpenMM> epsilons;
    std::vector<RealOpenMM> reductions;
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
    std::vector<RealOpenMM> radii;
    std::vector<RealOpenMM> epsilons;
    RealOpenMM epso; 
    RealOpenMM epsh; 
    RealOpenMM rmino; 
    RealOpenMM rminh; 
    RealOpenMM awater; 
    RealOpenMM shctd; 
    RealOpenMM dispoff;
    RealOpenMM slevy;
    RealOpenMM totalMaximumDispersionEnergy;
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
    int getIncludeCavityTerm( void ) const;

    /**
     *  Get the number of particles.
     *
     *  @return number of particles
     */
    int getNumParticles( void ) const;

    /**
     *  Get Direct Polarization flag.
     *
     *  @return directPolarization
     *
     */
    int getDirectPolarization( void ) const;

    /**
     *  Get the solute dielectric.
     *
     *  @return soluteDielectric
     *
     */
    RealOpenMM getSoluteDielectric( void ) const;

    /**
     *  Get the solvent dielectric.
     *
     *  @return solventDielectric
     *
     */
    RealOpenMM getSolventDielectric( void ) const;

    /**
     *  Get the dielectric offset.
     *
     *  @return dielectricOffset
     *
     */
    RealOpenMM getDielectricOffset( void ) const;

    /**
     *  Get the probe radius.
     *
     *  @return probeRadius
     *
     */
    RealOpenMM getProbeRadius( void ) const;

    /**
     *  Get the surface area factor.
     *
     *  @return surfaceAreaFactor
     *
     */
    RealOpenMM getSurfaceAreaFactor( void ) const;

    /**
     *  Get the vector of particle radii.
     *
     *  @param atomicRadii vector of atomic radii
     *
     */
    void getAtomicRadii( std::vector<RealOpenMM>& atomicRadii ) const;

    /**
     *  Get the vector of scale factors.
     *
     *  @param scaleFactors vector of scale factors
     *
     */
    void getScaleFactors( std::vector<RealOpenMM>& scaleFactors ) const;

    /**
     *  Get the vector of charges.
     *
     *  @param charges vector of charges
     *
     */
    void getCharges( std::vector<RealOpenMM>& charges ) const;

    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the AmoebaGeneralizedKirkwoodForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const AmoebaGeneralizedKirkwoodForce& force);

private:

    int numParticles;
    std::vector<RealOpenMM> atomicRadii;
    std::vector<RealOpenMM> scaleFactors;
    std::vector<RealOpenMM> charges;
    RealOpenMM soluteDielectric;
    RealOpenMM solventDielectric;
    RealOpenMM dielectricOffset;
    RealOpenMM probeRadius;
    RealOpenMM surfaceAreaFactor;
    int includeCavityTerm;
    int directPolarization;
    const System& system;
};

} // namespace OpenMM

#endif /*AMOEBA_OPENMM_REFERENCE_KERNELS_H*/
