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
#include "SimTKUtilities/SimTKOpenMMRealType.h"

namespace OpenMM {

/**
 * This kernel is invoked by AmoebaHarmonicBondForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcAmoebaHarmonicBondForceKernel : public CalcAmoebaHarmonicBondForceKernel {
public:
    ReferenceCalcAmoebaHarmonicBondForceKernel(std::string name, 
                                               const Platform& platform,
                                               System& system);
    ~ReferenceCalcAmoebaHarmonicBondForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the AmoebaHarmonicBondForce this kernel will be used for
     */
    void initialize(const System& system, const AmoebaHarmonicBondForce& force);
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
    int numBonds;
    std::vector<int>   particle1;
    std::vector<int>   particle2;
    std::vector<RealOpenMM> length;
    std::vector<RealOpenMM> kQuadratic;
    RealOpenMM globalHarmonicBondCubic;
    RealOpenMM globalHarmonicBondQuartic;
    System& system;
};

/**
 * This kernel is invoked by AmoebaUreyBradleyForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcAmoebaUreyBradleyForceKernel : public CalcAmoebaUreyBradleyForceKernel {
public:
    ReferenceCalcAmoebaUreyBradleyForceKernel(std::string name, 
                                               const Platform& platform,
                                               System& system);
    ~ReferenceCalcAmoebaUreyBradleyForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the AmoebaUreyBradleyForce this kernel will be used for
     */
    void initialize(const System& system, const AmoebaUreyBradleyForce& force);
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
    int numIxns;
    std::vector<int>   particle1;
    std::vector<int>   particle2;
    std::vector<RealOpenMM> length;
    std::vector<RealOpenMM> kQuadratic;
    RealOpenMM globalUreyBradleyCubic;
    RealOpenMM globalUreyBradleyQuartic;
    System& system;
};

/**
 * This kernel is invoked by AmoebaHarmonicAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcAmoebaHarmonicAngleForceKernel : public CalcAmoebaHarmonicAngleForceKernel {
public:
    ReferenceCalcAmoebaHarmonicAngleForceKernel(std::string name, const Platform& platform, System& system);
    ~ReferenceCalcAmoebaHarmonicAngleForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the AmoebaHarmonicAngleForce this kernel will be used for
     */
    void initialize(const System& system, const AmoebaHarmonicAngleForce& force);
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
    int numAngles;
    std::vector<int>   particle1;
    std::vector<int>   particle2;
    std::vector<int>   particle3;
    std::vector<RealOpenMM> angle;
    std::vector<RealOpenMM> kQuadratic;
    RealOpenMM globalHarmonicAngleCubic;
    RealOpenMM globalHarmonicAngleQuartic;
    RealOpenMM globalHarmonicAnglePentic;
    RealOpenMM globalHarmonicAngleSextic;
    System& system;
};

/**
 * This kernel is invoked by AmoebaHarmonicInPlaneAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcAmoebaHarmonicInPlaneAngleForceKernel : public CalcAmoebaHarmonicInPlaneAngleForceKernel {
public:
    ReferenceCalcAmoebaHarmonicInPlaneAngleForceKernel(std::string name, const Platform& platform, System& system);
    ~ReferenceCalcAmoebaHarmonicInPlaneAngleForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the AmoebaHarmonicInPlaneAngleForce this kernel will be used for
     */
    void initialize(const System& system, const AmoebaHarmonicInPlaneAngleForce& force);
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
    int numAngles;
    std::vector<int>   particle1;
    std::vector<int>   particle2;
    std::vector<int>   particle3;
    std::vector<int>   particle4;
    std::vector<RealOpenMM> angle;
    std::vector<RealOpenMM> kQuadratic;
    RealOpenMM globalHarmonicInPlaneAngleCubic;
    RealOpenMM globalHarmonicInPlaneAngleQuartic;
    RealOpenMM globalHarmonicInPlaneAnglePentic;
    RealOpenMM globalHarmonicInPlaneAngleSextic;
    System& system;
};

/**
 * This kernel is invoked by AmoebaTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcAmoebaTorsionForceKernel : public CalcAmoebaTorsionForceKernel {
public:
    ReferenceCalcAmoebaTorsionForceKernel(std::string name, const Platform& platform, System& system);
    ~ReferenceCalcAmoebaTorsionForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the AmoebaTorsionForce this kernel will be used for
     */
    void initialize(const System& system, const AmoebaTorsionForce& force);
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
    int numTorsions;
    std::vector<int>   particle1;
    std::vector<int>   particle2;
    std::vector<int>   particle3;
    std::vector<int>   particle4;
    std::vector< std::vector<RealOpenMM> > torsionParameters1;
    std::vector< std::vector<RealOpenMM> > torsionParameters2;
    std::vector< std::vector<RealOpenMM> > torsionParameters3;
    System& system;
};

/**
 * This kernel is invoked by AmoebaPiTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcAmoebaPiTorsionForceKernel : public CalcAmoebaPiTorsionForceKernel {
public:
    ReferenceCalcAmoebaPiTorsionForceKernel(std::string name, const Platform& platform, System& system);
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
private:
    int numPiTorsions;
    std::vector<int>   particle1;
    std::vector<int>   particle2;
    std::vector<int>   particle3;
    std::vector<int>   particle4;
    std::vector<int>   particle5;
    std::vector<int>   particle6;
    std::vector<RealOpenMM> kTorsion;
    System& system;
};

/**
 * This kernel is invoked by AmoebaStretchBendForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcAmoebaStretchBendForceKernel : public CalcAmoebaStretchBendForceKernel {
public:
    ReferenceCalcAmoebaStretchBendForceKernel(std::string name, const Platform& platform, System& system);
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
private:
    int numStretchBends;
    std::vector<int>   particle1;
    std::vector<int>   particle2;
    std::vector<int>   particle3;
    std::vector<RealOpenMM> lengthABParameters;
    std::vector<RealOpenMM> lengthCBParameters;
    std::vector<RealOpenMM> angleParameters;
    std::vector<RealOpenMM> kParameters;
    System& system;
};

/**
 * This kernel is invoked by AmoebaOutOfPlaneBendForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcAmoebaOutOfPlaneBendForceKernel : public CalcAmoebaOutOfPlaneBendForceKernel {
public:
    ReferenceCalcAmoebaOutOfPlaneBendForceKernel(std::string name, const Platform& platform, System& system);
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
    System& system;
};

/**
 * This kernel is invoked by AmoebaTorsionTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcAmoebaTorsionTorsionForceKernel : public CalcAmoebaTorsionTorsionForceKernel {
public:
    ReferenceCalcAmoebaTorsionTorsionForceKernel(std::string name, const Platform& platform, System& system);
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

    System& system;
};

/**
 * This kernel is invoked by AmoebaMultipoleForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcAmoebaMultipoleForceKernel : public CalcAmoebaMultipoleForceKernel {
public:
    ReferenceCalcAmoebaMultipoleForceKernel(std::string name, const Platform& platform, System& system);
    ~ReferenceCalcAmoebaMultipoleForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the AmoebaMultipoleForce this kernel will be used for
     */
    void initialize(const System& system, const AmoebaMultipoleForce& force);
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
     * Execute the kernel to calculate the electrostatic potential
     *
     * @param context        the context in which to execute this kernel
     * @param inputGrid      input grid coordinates
     * @param outputElectrostaticPotential output potential 
     */
    void getElectrostaticPotential(ContextImpl& context, const std::vector< Vec3 >& inputGrid,
                                   std::vector< double >& outputElectrostaticPotential );

    /**
     * Get the system multipole moments
     *
     * @param origin       origin
     * @param context      context
     * @param outputMultipoleMonents (charge,
                                      dipole_x, dipole_y, dipole_z,
                                      quadrupole_xx, quadrupole_xy, quadrupole_xz,
                                      quadrupole_yx, quadrupole_yy, quadrupole_yz,
                                      quadrupole_zx, quadrupole_zy, quadrupole_zz )
     */
    void getSystemMultipoleMoments(ContextImpl& context, const Vec3& origin, std::vector< double >& outputMultipoleMonents);

private:
    int numMultipoles;
    int polarizationType;
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

    //int iterativeMethod;
    int nonbondedMethod;
    int mutualInducedMaxIterations;
    RealOpenMM mutualInducedTargetEpsilon;

    System& system;
};

// /**
//  * This kernel is invoked by AmoebaMultipoleForce to calculate the forces acting on the system and the energy of the system.
//  */
// class ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel : public CalcAmoebaGeneralizedKirkwoodForceKernel {
// public:
//     ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel(std::string name, const Platform& platform, System& system);
//     ~ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel();
//     /**
//      * Initialize the kernel.
//      * 
//      * @param system     the System this kernel will be applied to
//      * @param force      the AmoebaMultipoleForce this kernel will be used for
//      */
//     void initialize(const System& system, const AmoebaGeneralizedKirkwoodForce& force);
//     /**
//      * Execute the kernel to calculate the forces and/or energy.
//      *
//      * @param context        the context in which to execute this kernel
//      * @param includeForces  true if forces should be calculated
//      * @param includeEnergy  true if the energy should be calculated
//      * @return the potential energy due to the force
//      */
//     double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
// private:
//     System& system;
// };

/**
 * This kernel is invoked to calculate the vdw forces acting on the system and the energy of the system.
 */
class ReferenceCalcAmoebaVdwForceKernel : public CalcAmoebaVdwForceKernel {
public:
    ReferenceCalcAmoebaVdwForceKernel(std::string name, const Platform& platform, System& system);
    ~ReferenceCalcAmoebaVdwForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the AmoebaMultipoleForce this kernel will be used for
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
private:
    int numParticles;
    std::vector<int> indexIVs;
    std::vector<int> indexClasses;
    std::vector< std::vector<int> > allExclusions;
    std::vector<RealOpenMM> sigmas;
    std::vector<RealOpenMM> epsilons;
    std::vector<RealOpenMM> reductions;
    std::string sigmaCombiningRule;
    std::string epsilonCombiningRule;
    System& system;
};

/**
 * This kernel is invoked to calculate the WCA dispersion forces acting on the system and the energy of the system.
 */
class ReferenceCalcAmoebaWcaDispersionForceKernel : public CalcAmoebaWcaDispersionForceKernel {
public:
    ReferenceCalcAmoebaWcaDispersionForceKernel(std::string name, const Platform& platform, System& system);
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
    System& system;
};

} // namespace OpenMM

#endif /*AMOEBA_OPENMM_REFERENCE_KERNELS_H*/
