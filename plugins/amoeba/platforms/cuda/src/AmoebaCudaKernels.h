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
 * Portions copyright (c) 2008-2020 Stanford University and the Authors.      *
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
#include "AmoebaCommonKernels.h"
#include <cufft.h>

namespace OpenMM {

class CudaCalcAmoebaGeneralizedKirkwoodForceKernel;

/**
 * This kernel is invoked by AmoebaMultipoleForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcAmoebaMultipoleForceKernel : public CommonCalcAmoebaMultipoleForceKernel {
public:
    CudaCalcAmoebaMultipoleForceKernel(const std::string& name, const Platform& platform, CudaContext& cu, const System& system) :
            CommonCalcAmoebaMultipoleForceKernel(name, platform, cu, system), hasInitializedFFT(false) {
    }
    ~CudaCalcAmoebaMultipoleForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the AmoebaMultipoleForce this kernel will be used for
     */
    void initialize(const System& system, const AmoebaMultipoleForce& force);
    /**
     * Compute the FFT in the forward direction.
     */
    void computeForwardFFT();
    /**
     * Compute the FFT in the inverse direction.
     */
    void computeInverseFFT();
private:
    bool hasInitializedFFT;
    cufftHandle fft;
};

/**
 * This kernel is invoked by HippoNonbondedForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcHippoNonbondedForceKernel : public CalcHippoNonbondedForceKernel {
public:
    CudaCalcHippoNonbondedForceKernel(const std::string& name, const Platform& platform, CudaContext& cu, const System& system);
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
