#ifndef AMOEBA_OPENMM_CUDAKERNELS_H_
#define AMOEBA_OPENMM_CUDAKERNELS_H_

/* -------------------------------------------------------------------------- *
 *                              OpenMMAmoeba                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2018 Stanford University and the Authors.      *
 * Authors: Mark Friedrichs, Peter Eastman                                    *
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

#include "openmm/amoebaKernels.h"
#include "openmm/kernels.h"
#include "openmm/System.h"
#include "CudaArray.h"
#include "CudaContext.h"
#include "CudaNonbondedUtilities.h"
#include "CudaSort.h"
#include <cufft.h>

namespace OpenMM {

class CudaCalcAmoebaGeneralizedKirkwoodForceKernel;

/**
 * This kernel is invoked by AmoebaBondForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcAmoebaBondForceKernel : public CalcAmoebaBondForceKernel {
public:
    CudaCalcAmoebaBondForceKernel(std::string name, 
                                          const Platform& platform,
                                          CudaContext& cu,
                                          const System& system);
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
    class ForceInfo;
    int numBonds;
    CudaContext& cu;
    const System& system;
    CudaArray params;
};

/**
 * This kernel is invoked by AmoebaAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcAmoebaAngleForceKernel : public CalcAmoebaAngleForceKernel {
public:
    CudaCalcAmoebaAngleForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
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
    class ForceInfo;
    int numAngles;
    CudaContext& cu;
    const System& system;
    CudaArray params;
};

/**
 * This kernel is invoked by AmoebaInPlaneAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcAmoebaInPlaneAngleForceKernel : public CalcAmoebaInPlaneAngleForceKernel {
public:
    CudaCalcAmoebaInPlaneAngleForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
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
    class ForceInfo;
    int numAngles;
    CudaContext& cu;
    const System& system;
    CudaArray params;
};

/**
 * This kernel is invoked by AmoebaPiTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcAmoebaPiTorsionForceKernel : public CalcAmoebaPiTorsionForceKernel {
public:
    CudaCalcAmoebaPiTorsionForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
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
    class ForceInfo;
    int numPiTorsions;
    CudaContext& cu;
    const System& system;
    CudaArray params;
};

/**
 * This kernel is invoked by AmoebaStretchBendForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcAmoebaStretchBendForceKernel : public CalcAmoebaStretchBendForceKernel {
public:
    CudaCalcAmoebaStretchBendForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
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
    class ForceInfo;
    int numStretchBends;
    CudaContext& cu;
    const System& system;
    CudaArray params1; // Equilibrium values
    CudaArray params2; // force constants
};

/**
 * This kernel is invoked by AmoebaOutOfPlaneBendForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcAmoebaOutOfPlaneBendForceKernel : public CalcAmoebaOutOfPlaneBendForceKernel {
public:
    CudaCalcAmoebaOutOfPlaneBendForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
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
    class ForceInfo;
    int numOutOfPlaneBends;
    CudaContext& cu;
    const System& system;
    CudaArray params;
};

/**
 * This kernel is invoked by AmoebaTorsionTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcAmoebaTorsionTorsionForceKernel : public CalcAmoebaTorsionTorsionForceKernel {
public:
    CudaCalcAmoebaTorsionTorsionForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
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
    class ForceInfo;
    int numTorsionTorsions;
    int numTorsionTorsionGrids;
    CudaContext& cu;
    const System& system;
    CudaArray gridValues;
    CudaArray gridParams;
    CudaArray torsionParams;
};

/**
 * This kernel is invoked by AmoebaMultipoleForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcAmoebaMultipoleForceKernel : public CalcAmoebaMultipoleForceKernel {
public:
    CudaCalcAmoebaMultipoleForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
    ~CudaCalcAmoebaMultipoleForceKernel();
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
     * Get the LabFrame dipole moments of all particles.
     * 
     * @param context    the Context for which to get the induced dipoles
     * @param dipoles    the induced dipole moment of particle i is stored into the i'th element
     */
    void getLabFramePermanentDipoles(ContextImpl& context, std::vector<Vec3>& dipoles);
    /**
     * Get the induced dipole moments of all particles.
     * 
     * @param context    the Context for which to get the induced dipoles
     * @param dipoles    the induced dipole moment of particle i is stored into the i'th element
     */
    void getInducedDipoles(ContextImpl& context, std::vector<Vec3>& dipoles);
    /**
     * Get the total dipole moments of all particles.
     * 
     * @param context    the Context for which to get the induced dipoles
     * @param dipoles    the induced dipole moment of particle i is stored into the i'th element
     */
    void getTotalDipoles(ContextImpl& context, std::vector<Vec3>& dipoles);
    /**
     * Execute the kernel to calculate the electrostatic potential
     *
     * @param context        the context in which to execute this kernel
     * @param inputGrid      input grid coordinates
     * @param outputElectrostaticPotential output potential 
     */
    void getElectrostaticPotential(ContextImpl& context, const std::vector< Vec3 >& inputGrid,
                                   std::vector< double >& outputElectrostaticPotential);

   /** 
     * Get the system multipole moments
     *
     * @param context      context
     * @param outputMultipoleMoments (charge,
     *                                dipole_x, dipole_y, dipole_z,
     *                                quadrupole_xx, quadrupole_xy, quadrupole_xz,
     *                                quadrupole_yx, quadrupole_yy, quadrupole_yz,
     *                                quadrupole_zx, quadrupole_zy, quadrupole_zz)
     */
    void getSystemMultipoleMoments(ContextImpl& context, std::vector<double>& outputMultipoleMoments);
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
    class ForceInfo;
    void initializeScaleFactors();
    void computeInducedField(void** recipBoxVectorPointer);
    bool iterateDipolesByDIIS(int iteration);
    void computeExtrapolatedDipoles(void** recipBoxVectorPointer);
    void ensureMultipolesValid(ContextImpl& context);
    template <class T, class T4, class M4> void computeSystemMultipoleMoments(ContextImpl& context, std::vector<double>& outputMultipoleMoments);
    int numMultipoles, maxInducedIterations, maxExtrapolationOrder;
    int fixedFieldThreads, inducedFieldThreads, electrostaticsThreads;
    int gridSizeX, gridSizeY, gridSizeZ;
    double alpha, inducedEpsilon;
    bool usePME, hasQuadrupoles, hasInitializedScaleFactors, hasInitializedFFT, multipolesAreValid, hasCreatedEvent;
    AmoebaMultipoleForce::PolarizationType polarizationType;
    CudaContext& cu;
    const System& system;
    std::vector<int3> covalentFlagValues;
    std::vector<int2> polarizationFlagValues;
    CudaArray multipoleParticles;
    CudaArray molecularDipoles;
    CudaArray molecularQuadrupoles;
    CudaArray labFrameDipoles;
    CudaArray labFrameQuadrupoles;
    CudaArray sphericalDipoles;
    CudaArray sphericalQuadrupoles;
    CudaArray fracDipoles;
    CudaArray fracQuadrupoles;
    CudaArray field;
    CudaArray fieldPolar;
    CudaArray inducedField;
    CudaArray inducedFieldPolar;
    CudaArray torque;
    CudaArray dampingAndThole;
    CudaArray inducedDipole;
    CudaArray inducedDipolePolar;
    CudaArray inducedDipoleErrors;
    CudaArray prevDipoles;
    CudaArray prevDipolesPolar;
    CudaArray prevDipolesGk;
    CudaArray prevDipolesGkPolar;
    CudaArray prevErrors;
    CudaArray diisMatrix;
    CudaArray diisCoefficients;
    CudaArray extrapolatedDipole;
    CudaArray extrapolatedDipolePolar;
    CudaArray extrapolatedDipoleGk;
    CudaArray extrapolatedDipoleGkPolar;
    CudaArray inducedDipoleFieldGradient;
    CudaArray inducedDipoleFieldGradientPolar;
    CudaArray inducedDipoleFieldGradientGk;
    CudaArray inducedDipoleFieldGradientGkPolar;
    CudaArray extrapolatedDipoleFieldGradient;
    CudaArray extrapolatedDipoleFieldGradientPolar;
    CudaArray extrapolatedDipoleFieldGradientGk;
    CudaArray extrapolatedDipoleFieldGradientGkPolar;
    CudaArray polarizability;
    CudaArray covalentFlags;
    CudaArray polarizationGroupFlags;
    CudaArray pmeGrid;
    CudaArray pmeBsplineModuliX;
    CudaArray pmeBsplineModuliY;
    CudaArray pmeBsplineModuliZ;
    CudaArray pmePhi;
    CudaArray pmePhid;
    CudaArray pmePhip;
    CudaArray pmePhidp;
    CudaArray pmeCphi;
    CudaArray lastPositions;
    cufftHandle fft;
    CUfunction computeMomentsKernel, recordInducedDipolesKernel, computeFixedFieldKernel, computeInducedFieldKernel, updateInducedFieldKernel, electrostaticsKernel, mapTorqueKernel;
    CUfunction pmeSpreadFixedMultipolesKernel, pmeSpreadInducedDipolesKernel, pmeFinishSpreadChargeKernel, pmeConvolutionKernel;
    CUfunction pmeFixedPotentialKernel, pmeInducedPotentialKernel, pmeFixedForceKernel, pmeInducedForceKernel, pmeRecordInducedFieldDipolesKernel, computePotentialKernel;
    CUfunction recordDIISDipolesKernel, buildMatrixKernel, solveMatrixKernel;
    CUfunction initExtrapolatedKernel, iterateExtrapolatedKernel, computeExtrapolatedKernel, addExtrapolatedGradientKernel;
    CUfunction pmeTransformMultipolesKernel, pmeTransformPotentialKernel;
    CUevent syncEvent;
    CudaCalcAmoebaGeneralizedKirkwoodForceKernel* gkKernel;
    static const int PmeOrder = 5;
    static const int MaxPrevDIISDipoles = 20;
};

/**
 * This kernel is invoked by AmoebaMultipoleForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcAmoebaGeneralizedKirkwoodForceKernel : public CalcAmoebaGeneralizedKirkwoodForceKernel {
public:
    CudaCalcAmoebaGeneralizedKirkwoodForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
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
     * Perform the computation of Born radii.
     */
    void computeBornRadii();
    /**
     * Perform the final parts of the force/energy computation.
     */
    void finishComputation(CudaArray& torque, CudaArray& labFrameDipoles, CudaArray& labFrameQuadrupoles, CudaArray& inducedDipole, CudaArray& inducedDipolePolar, CudaArray& dampingAndThole, CudaArray& covalentFlags, CudaArray& polarizationGroupFlags);
    CudaArray& getBornRadii() {
        return bornRadii;
    }
    CudaArray& getField() {
        return field;
    }
    CudaArray& getInducedField() {
        return inducedField;
    }
    CudaArray& getInducedFieldPolar() {
        return inducedFieldPolar;
    }
    CudaArray& getInducedDipoles() {
        return inducedDipoleS;
    }
    CudaArray& getInducedDipolesPolar() {
        return inducedDipolePolarS;
    }
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the AmoebaGeneralizedKirkwoodForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const AmoebaGeneralizedKirkwoodForce& force);
private:
    class ForceInfo;
    CudaContext& cu;
    const System& system;
    bool includeSurfaceArea, hasInitializedKernels;
    int computeBornSumThreads, gkForceThreads, chainRuleThreads, ediffThreads;
    AmoebaMultipoleForce::PolarizationType polarizationType;
    std::map<std::string, std::string> defines;
    CudaArray params;
    CudaArray bornSum;
    CudaArray bornRadii;
    CudaArray bornForce;
    CudaArray field;
    CudaArray inducedField;
    CudaArray inducedFieldPolar;
    CudaArray inducedDipoleS;
    CudaArray inducedDipolePolarS;
    CUfunction computeBornSumKernel, reduceBornSumKernel, surfaceAreaKernel, gkForceKernel, chainRuleKernel, ediffKernel;
};

/**
 * This kernel is invoked to calculate the vdw forces acting on the system and the energy of the system.
 */
class CudaCalcAmoebaVdwForceKernel : public CalcAmoebaVdwForceKernel {
public:
    CudaCalcAmoebaVdwForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
    ~CudaCalcAmoebaVdwForceKernel();
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
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the AmoebaVdwForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const AmoebaVdwForce& force);
private:
    class ForceInfo;
    CudaContext& cu;
    const System& system;
    bool hasInitializedNonbonded;
    double dispersionCoefficient;
    CudaArray sigmaEpsilon;
    CudaArray bondReductionAtoms;
    CudaArray bondReductionFactors;
    CudaArray tempPosq;
    CudaArray tempForces;
    CudaNonbondedUtilities* nonbonded;
    CUfunction prepareKernel, spreadKernel;
};

/**
 * This kernel is invoked to calculate the WCA dispersion forces acting on the system and the energy of the system.
 */
class CudaCalcAmoebaWcaDispersionForceKernel : public CalcAmoebaWcaDispersionForceKernel {
public:
    CudaCalcAmoebaWcaDispersionForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
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
    class ForceInfo;
    CudaContext& cu;
    const System& system;
    double totalMaximumDispersionEnergy;
    CudaArray radiusEpsilon;
    CUfunction forceKernel;
};

/**
 * This kernel is invoked by HippoNonbondedForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcHippoNonbondedForceKernel : public CalcHippoNonbondedForceKernel {
public:
    CudaCalcHippoNonbondedForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
    ~CudaCalcHippoNonbondedForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the HippoNonbondedForce this kernel will be used for
     */
    void initialize(const System& system, const HippoNonbondedForce& force);
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
     * Calculate the electrostatic potential given vector of grid coordinates.
     *
     * @param context                      context
     * @param inputGrid                    input grid coordinates
     * @param outputElectrostaticPotential output potential 
     */
    void getElectrostaticPotential(ContextImpl& context, const std::vector< Vec3 >& inputGrid,
                                   std::vector< double >& outputElectrostaticPotential);
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
    class ForceInfo;
    class TorquePostComputation;
    class SortTrait : public CudaSort::SortTrait {
        int getDataSize() const {return 8;}
        int getKeySize() const {return 4;}
        const char* getDataType() const {return "int2";}
        const char* getKeyType() const {return "int";}
        const char* getMinKey() const {return "(-2147483647-1)";}
        const char* getMaxKey() const {return "2147483647";}
        const char* getMaxValue() const {return "make_int2(2147483647, 2147483647)";}
        const char* getSortKey() const {return "value.y";}
    };
    void computeInducedField(void** recipBoxVectorPointer, int optOrder);
    void computeExtrapolatedDipoles(void** recipBoxVectorPointer);
    void ensureMultipolesValid(ContextImpl& context);
    void addTorquesToForces();
    void createFieldKernel(const std::string& interactionSrc, std::vector<CudaArray*> params, CudaArray& fieldBuffer,
        CUfunction& kernel, std::vector<void*>& args, CUfunction& exceptionKernel, std::vector<void*>& exceptionArgs,
        CudaArray& exceptionScale);
    int numParticles, maxExtrapolationOrder, maxTiles;
    int gridSizeX, gridSizeY, gridSizeZ;
    int dispersionGridSizeX, dispersionGridSizeY, dispersionGridSizeZ;
    double pmeAlpha, dpmeAlpha, cutoff;
    bool usePME, hasInitializedKernels, hasInitializedFFT, multipolesAreValid;
    std::vector<double> extrapolationCoefficients;
    CudaContext& cu;
    const System& system;
    CudaArray multipoleParticles;
    CudaArray coreCharge, valenceCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability;
    CudaArray localDipoles, labDipoles, fracDipoles;
    CudaArray localQuadrupoles, labQuadrupoles[5], fracQuadrupoles;
    CudaArray field;
    CudaArray inducedField;
    CudaArray torque;
    CudaArray inducedDipole;
    CudaArray extrapolatedDipole, extrapolatedPhi;
    CudaArray pmeGrid1, pmeGrid2;
    CudaArray pmeAtomGridIndex;
    CudaArray pmeBsplineModuliX, pmeBsplineModuliY, pmeBsplineModuliZ;
    CudaArray dpmeBsplineModuliX, dpmeBsplineModuliY, dpmeBsplineModuliZ;
    CudaArray pmePhi, pmePhidp, pmeCphi;
    CudaArray lastPositions;
    CudaArray exceptionScales[6];
    CudaArray exceptionAtoms;
    CudaSort* sort;
    cufftHandle fftForward, fftBackward, dfftForward, dfftBackward;
    CUfunction computeMomentsKernel, fixedFieldKernel, fixedFieldExceptionKernel, mutualFieldKernel, mutualFieldExceptionKernel, computeExceptionsKernel;
    CUfunction recordInducedDipolesKernel, mapTorqueKernel;
    CUfunction pmeSpreadFixedMultipolesKernel, pmeSpreadInducedDipolesKernel, pmeFinishSpreadChargeKernel, pmeConvolutionKernel;
    CUfunction pmeFixedPotentialKernel, pmeInducedPotentialKernel, pmeFixedForceKernel, pmeInducedForceKernel, pmeRecordInducedFieldDipolesKernel;
    CUfunction pmeSelfEnergyKernel;
    CUfunction dpmeGridIndexKernel, dpmeSpreadChargeKernel, dpmeFinishSpreadChargeKernel, dpmeEvalEnergyKernel, dpmeConvolutionKernel, dpmeInterpolateForceKernel;
    CUfunction initExtrapolatedKernel, iterateExtrapolatedKernel, computeExtrapolatedKernel, polarizationEnergyKernel;
    CUfunction pmeTransformMultipolesKernel, pmeTransformPotentialKernel;
    std::vector<void*> fixedFieldArgs, fixedFieldExceptionArgs, mutualFieldArgs, mutualFieldExceptionArgs, computeExceptionsArgs;
    static const int PmeOrder = 5;
};

} // namespace OpenMM

#endif /*AMOEBA_OPENMM_CUDAKERNELS_H*/
