#ifndef AMOEBA_OPENMM_COMMONKERNELS_H_
#define AMOEBA_OPENMM_COMMONKERNELS_H_

/* -------------------------------------------------------------------------- *
 *                              OpenMMAmoeba                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2021 Stanford University and the Authors.      *
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
#include "openmm/common/ComputeContext.h"
#include "openmm/common/ComputeArray.h"
#include "openmm/common/NonbondedUtilities.h"

namespace OpenMM {

class CommonCalcAmoebaGeneralizedKirkwoodForceKernel;

/**
 * This kernel is invoked by AmoebaTorsionTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CommonCalcAmoebaTorsionTorsionForceKernel : public CalcAmoebaTorsionTorsionForceKernel {
public:
    CommonCalcAmoebaTorsionTorsionForceKernel(const std::string& name, const Platform& platform, ComputeContext& cc, const System& system);
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
    ComputeContext& cc;
    const System& system;
    ComputeArray gridValues;
    ComputeArray gridParams;
    ComputeArray torsionParams;
};

/**
 * This kernel is invoked by AmoebaMultipoleForce to calculate the forces acting on the system and the energy of the system.
 */
class CommonCalcAmoebaMultipoleForceKernel : public CalcAmoebaMultipoleForceKernel {
public:
    CommonCalcAmoebaMultipoleForceKernel(const std::string& name, const Platform& platform, ComputeContext& cc, const System& system);
    ~CommonCalcAmoebaMultipoleForceKernel();
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
    /**
     * Compute the FFT.
     */
    virtual void computeFFT(bool forward) = 0;
    /**
     * Get whether charge spreading should be done in fixed point.
     */
    virtual bool useFixedPointChargeSpreading() const = 0;
protected:
    class ForceInfo;
    void initializeScaleFactors();
    void computeInducedField();
    bool iterateDipolesByDIIS(int iteration);
    void computeExtrapolatedDipoles();
    void ensureMultipolesValid(ContextImpl& context);
    template <class T, class T4, class M4> void computeSystemMultipoleMoments(ContextImpl& context, std::vector<double>& outputMultipoleMoments);
    int numMultipoles, maxInducedIterations, maxExtrapolationOrder;
    int fixedFieldThreads, inducedFieldThreads, electrostaticsThreads;
    int gridSizeX, gridSizeY, gridSizeZ;
    double pmeAlpha, inducedEpsilon;
    bool usePME, hasQuadrupoles, hasInitializedScaleFactors, multipolesAreValid, hasCreatedEvent;
    AmoebaMultipoleForce::PolarizationType polarizationType;
    ComputeContext& cc;
    const System& system;
    std::vector<mm_int4> covalentFlagValues;
    std::vector<mm_int2> polarizationFlagValues;
    ComputeArray multipoleParticles;
    ComputeArray localDipoles;
    ComputeArray localQuadrupoles;
    ComputeArray labDipoles;
    ComputeArray labQuadrupoles;
    ComputeArray sphericalDipoles;
    ComputeArray sphericalQuadrupoles;
    ComputeArray fracDipoles;
    ComputeArray fracQuadrupoles;
    ComputeArray field;
    ComputeArray fieldPolar;
    ComputeArray inducedField;
    ComputeArray inducedFieldPolar;
    ComputeArray torque;
    ComputeArray dampingAndThole;
    ComputeArray inducedDipole;
    ComputeArray inducedDipolePolar;
    ComputeArray inducedDipoleErrors;
    ComputeArray prevDipoles;
    ComputeArray prevDipolesPolar;
    ComputeArray prevDipolesGk;
    ComputeArray prevDipolesGkPolar;
    ComputeArray prevErrors;
    ComputeArray diisMatrix;
    ComputeArray diisCoefficients;
    ComputeArray extrapolatedDipole;
    ComputeArray extrapolatedDipolePolar;
    ComputeArray extrapolatedDipoleGk;
    ComputeArray extrapolatedDipoleGkPolar;
    ComputeArray inducedDipoleFieldGradient;
    ComputeArray inducedDipoleFieldGradientPolar;
    ComputeArray inducedDipoleFieldGradientGk;
    ComputeArray inducedDipoleFieldGradientGkPolar;
    ComputeArray extrapolatedDipoleFieldGradient;
    ComputeArray extrapolatedDipoleFieldGradientPolar;
    ComputeArray extrapolatedDipoleFieldGradientGk;
    ComputeArray extrapolatedDipoleFieldGradientGkPolar;
    ComputeArray polarizability;
    ComputeArray covalentFlags;
    ComputeArray polarizationGroupFlags;
    ComputeArray pmeGrid1;
    ComputeArray pmeGrid2;
    ComputeArray pmeGridLong;
    ComputeArray pmeBsplineModuliX;
    ComputeArray pmeBsplineModuliY;
    ComputeArray pmeBsplineModuliZ;
    ComputeArray pmePhi;
    ComputeArray pmePhid;
    ComputeArray pmePhip;
    ComputeArray pmePhidp;
    ComputeArray pmeCphi;
    ComputeArray lastPositions;
    ComputeKernel computeMomentsKernel, recordInducedDipolesKernel, mapTorqueKernel, computePotentialKernel, electrostaticsKernel;
    ComputeKernel computeFixedFieldKernel, computeInducedFieldKernel, updateInducedFieldKernel;
    ComputeKernel recordDIISDipolesKernel, buildMatrixKernel, solveMatrixKernel;
    ComputeKernel initExtrapolatedKernel, iterateExtrapolatedKernel, computeExtrapolatedKernel, addExtrapolatedGradientKernel;
    ComputeKernel pmeSpreadFixedMultipolesKernel, pmeSpreadInducedDipolesKernel, pmeFinishSpreadChargeKernel, pmeConvolutionKernel;
    ComputeKernel pmeFixedPotentialKernel, pmeInducedPotentialKernel, pmeFixedForceKernel, pmeInducedForceKernel, pmeRecordInducedFieldDipolesKernel;
    ComputeKernel pmeTransformMultipolesKernel, pmeTransformPotentialKernel;
    ComputeEvent syncEvent;
    CommonCalcAmoebaGeneralizedKirkwoodForceKernel* gkKernel;
    static const int PmeOrder = 5;
    static const int MaxPrevDIISDipoles = 20;
};

/**
 * This kernel is invoked by AmoebaMultipoleForce to calculate the forces acting on the system and the energy of the system.
 */
class CommonCalcAmoebaGeneralizedKirkwoodForceKernel : public CalcAmoebaGeneralizedKirkwoodForceKernel {
public:
    CommonCalcAmoebaGeneralizedKirkwoodForceKernel(const std::string& name, const Platform& platform, ComputeContext& cc, const System& system);
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
    void computeBornRadii(ComputeArray& torque, ComputeArray& labFrameDipoles, ComputeArray& labFrameQuadrupoles, ComputeArray& inducedDipole, ComputeArray& inducedDipolePolar, ComputeArray& dampingAndThole, ComputeArray& covalentFlags, ComputeArray& polarizationGroupFlags);
    /**
     * Perform the final parts of the force/energy computation.
     */
    void finishComputation();
    ComputeArray& getBornRadii() {
        return bornRadii;
    }
    ComputeArray& getField() {
        return field;
    }
    ComputeArray& getInducedField() {
        return inducedField;
    }
    ComputeArray& getInducedFieldPolar() {
        return inducedFieldPolar;
    }
    ComputeArray& getInducedDipoles() {
        return inducedDipoleS;
    }
    ComputeArray& getInducedDipolesPolar() {
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
    ComputeContext& cc;
    const System& system;
    bool includeSurfaceArea, hasInitializedKernels;
    int computeBornSumThreads, gkForceThreads, chainRuleThreads, ediffThreads;
    AmoebaMultipoleForce::PolarizationType polarizationType;
    std::map<std::string, std::string> defines;
    ComputeArray params;
    ComputeArray bornSum;
    ComputeArray bornRadii;
    ComputeArray bornForce;
    ComputeArray field;
    ComputeArray inducedField;
    ComputeArray inducedFieldPolar;
    ComputeArray inducedDipoleS;
    ComputeArray inducedDipolePolarS;
    ComputeKernel computeBornSumKernel, reduceBornSumKernel, surfaceAreaKernel, gkForceKernel, chainRuleKernel, ediffKernel;
};

/**
 * This kernel is invoked to calculate the vdw forces acting on the system and the energy of the system.
 */
class CommonCalcAmoebaVdwForceKernel : public CalcAmoebaVdwForceKernel {
public:
    CommonCalcAmoebaVdwForceKernel(const std::string& name, const Platform& platform, ComputeContext& cc, const System& system);
    ~CommonCalcAmoebaVdwForceKernel();
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
    class ForceInfo;
    ComputeContext& cc;
    const System& system;
    bool hasInitializedNonbonded;

    // True if the AmoebaVdwForce AlchemicalMethod is not None.
    bool hasAlchemical;
    // Device memory for the alchemical state.
    ComputeArray vdwLambda;
    // Only update device memory when lambda changes.
    float currentVdwLambda;
    // Per particle alchemical flag.
    ComputeArray isAlchemical;

    double dispersionCoefficient;
    ComputeArray sigmaEpsilon, atomType;
    ComputeArray bondReductionAtoms;
    ComputeArray bondReductionFactors;
    ComputeArray tempPosq;
    ComputeArray tempForces;
    NonbondedUtilities* nonbonded;
    ComputeKernel prepareKernel, spreadKernel;
};

/**
 * This kernel is invoked to calculate the WCA dispersion forces acting on the system and the energy of the system.
 */
class CommonCalcAmoebaWcaDispersionForceKernel : public CalcAmoebaWcaDispersionForceKernel {
public:
    CommonCalcAmoebaWcaDispersionForceKernel(const std::string& name, const Platform& platform, ComputeContext& cc, const System& system);
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
    ComputeContext& cc;
    const System& system;
    double totalMaximumDispersionEnergy;
    int forceThreadBlockSize;
    ComputeArray radiusEpsilon;
    ComputeKernel forceKernel;
};

/**
 * This kernel is invoked by HippoNonbondedForce to calculate the forces acting on the system and the energy of the system.
 */
class CommonCalcHippoNonbondedForceKernel : public CalcHippoNonbondedForceKernel {
public:
    CommonCalcHippoNonbondedForceKernel(const std::string& name, const Platform& platform, ComputeContext& cc, const System& system);
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the HippoNonbondedForce this kernel will be used for
     */
    void initialize(const System& system, const HippoNonbondedForce& force);
    /**
     * Compute the FFT.
     */
    virtual void computeFFT(bool forward, bool dispersion) = 0;
    /**
     * Get whether charge spreading should be done in fixed point.
     */
    virtual bool useFixedPointChargeSpreading() const = 0;
    /**
     * Sort the atom grid indices.
     */
    virtual void sortGridIndex() = 0;
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
protected:
    class ForceInfo;
    class TorquePostComputation;
    void computeInducedField(int optOrder);
    void computeExtrapolatedDipoles();
    void ensureMultipolesValid(ContextImpl& context);
    void addTorquesToForces();
    void createFieldKernel(const std::string& interactionSrc, std::vector<ComputeArray*> params, ComputeArray& fieldBuffer,
        ComputeKernel& kernel, ComputeKernel& exceptionKernel, ComputeArray& exceptionScale);
    int numParticles, maxExtrapolationOrder, maxTiles, fieldThreadBlockSize;
    int gridSizeX, gridSizeY, gridSizeZ;
    int dispersionGridSizeX, dispersionGridSizeY, dispersionGridSizeZ;
    double pmeAlpha, dpmeAlpha, cutoff;
    bool usePME, hasInitializedKernels, multipolesAreValid;
    std::vector<double> extrapolationCoefficients;
    ComputeContext& cc;
    const System& system;
    ComputeArray multipoleParticles;
    ComputeArray coreCharge, valenceCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability;
    ComputeArray localDipoles, labDipoles, fracDipoles;
    ComputeArray localQuadrupoles, labQuadrupoles[5], fracQuadrupoles;
    ComputeArray field;
    ComputeArray inducedField;
    ComputeArray torque;
    ComputeArray inducedDipole;
    ComputeArray extrapolatedDipole, extrapolatedPhi;
    ComputeArray pmeGrid1, pmeGrid2, pmeGridLong;
    ComputeArray pmeAtomGridIndex;
    ComputeArray pmeBsplineModuliX, pmeBsplineModuliY, pmeBsplineModuliZ;
    ComputeArray dpmeBsplineModuliX, dpmeBsplineModuliY, dpmeBsplineModuliZ;
    ComputeArray pmePhi, pmePhidp, pmeCphi;
    ComputeArray lastPositions;
    ComputeArray exceptionScales[6];
    ComputeArray exceptionAtoms;
    ComputeKernel computeMomentsKernel, recordInducedDipolesKernel, mapTorqueKernel;
    ComputeKernel fixedFieldKernel, fixedFieldExceptionKernel, mutualFieldKernel, mutualFieldExceptionKernel, computeExceptionsKernel;
    ComputeKernel pmeSpreadFixedMultipolesKernel, pmeSpreadInducedDipolesKernel, pmeFinishSpreadChargeKernel, pmeConvolutionKernel;
    ComputeKernel pmeFixedPotentialKernel, pmeInducedPotentialKernel, pmeFixedForceKernel, pmeInducedForceKernel, pmeRecordInducedFieldDipolesKernel;
    ComputeKernel pmeSelfEnergyKernel, pmeTransformMultipolesKernel, pmeTransformPotentialKernel;
    ComputeKernel dpmeGridIndexKernel, dpmeSpreadChargeKernel, dpmeFinishSpreadChargeKernel, dpmeEvalEnergyKernel, dpmeConvolutionKernel, dpmeInterpolateForceKernel;
    ComputeKernel initExtrapolatedKernel, iterateExtrapolatedKernel, computeExtrapolatedKernel, polarizationEnergyKernel;
    static const int PmeOrder = 5;
};

} // namespace OpenMM

#endif /*AMOEBA_OPENMM_COMMONKERNELS_H_*/