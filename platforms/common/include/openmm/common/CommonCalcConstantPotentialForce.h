#ifndef OPENMM_COMMONCALCCONSTANTPOTENTIALFORCE_H_
#define OPENMM_COMMONCALCCONSTANTPOTENTIALFORCE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2008-2025 Stanford University and the Authors.      *
 * Authors: Evan Pretti                                                       *
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

#include "openmm/kernels.h"
#include "openmm/common/ComputeArray.h"
#include "openmm/common/ComputeContext.h"
#include "openmm/common/ComputeEvent.h"
#include "openmm/common/ComputeQueue.h"
#include "openmm/common/ComputeSort.h"
#include "openmm/common/FFT3D.h"
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace OpenMM {

class CommonCalcConstantPotentialForceKernel;

/**
 * A generic charge solver for the constant potential method.
 */
class CommonConstantPotentialSolver {
public:
    /**
     * Creates a CommonConstantPotentialSolver.
     * 
     * @param cc                     the compute context (should be selected)
     * @param numParticles           the number of particles
     * @param numElectrodeParticles  the number of electrode (fluctuating-charge) particles (should be positive)
     * @param paddedProblemSize      the padded number of electrode particles
     */
    CommonConstantPotentialSolver(ComputeContext& cc, int numParticles, int numElectrodeParticles, int paddedProblemSize);
    virtual ~CommonConstantPotentialSolver();
    /**
     * Sets up solver kernels.  Should be called by derived classes if
     * overridden.
     * 
     * @param kernel  main constant potential kernel
     */
    virtual void compileKernels(CommonCalcConstantPotentialForceKernel& kernel);
    /**
     * Mark precomputed data stored by the solver as invalid due to a change in
     * electrode parameters.
     */
    void invalidate();
    /**
     * Mark any existing solution stored by the solver as invalid due to a
     * change in system parameters.
     */
    void discardSavedSolution();
    /**
     * Solves for charges if needed.
     * 
     * @param kernel  main constant potential kernel
     */
    void solve(CommonCalcConstantPotentialForceKernel& kernel);
    /**
     * Solves for charges.
     * 
     * @param kernel  main constant potential kernel
     */
    virtual void solveImpl(CommonCalcConstantPotentialForceKernel& kernel) = 0;
    /**
     * Retrieves a list of arrays of initial guess charges to be reordered.
     * 
     * @param arrays  Arrays to be reordered.
     */
    virtual void getGuessChargeArrays(std::vector<ComputeArray*>& arrays);

protected:
    int numParticles, numElectrodeParticles, paddedProblemSize;
    bool valid, hasSavedSolution;
    Vec3 savedBoxVectors[3];
    ComputeArray savedPositions;
    ComputeArray savedCharges;
    ComputeArray checkSavedPositionsKernelResult;
    ComputeKernel checkSavedPositionsKernel;
};

/**
 * A constant potential solver using direct inversion of the Coulomb matrix.
 * Suitable only when electrode particle positions are fixed.
 */
class CommonConstantPotentialMatrixSolver : public CommonConstantPotentialSolver {
private:
    Vec3 boxVectors[3];
    ComputeArray electrodePosData;
    ComputeArray capacitance;
    ComputeArray constraintVector;
    ComputeArray checkSavedElectrodePositionsKernelResult;
    ComputeKernel checkSavedElectrodePositionsKernel;
    ComputeKernel saveElectrodePositionsKernel;
    ComputeKernel solveKernel;

public:
    /**
     * Creates a CommonConstantPotentialMatrixSolver.
     * 
     * @param cc                     the compute context (should be selected)
     * @param numParticles           the number of particles
     * @param numElectrodeParticles  the number of electrode (fluctuating-charge) particles
     * @param paddedProblemSize      the padded number of electrode particles
     */
    CommonConstantPotentialMatrixSolver(ComputeContext& cc, int numParticles, int numElectrodeParticles, int paddedProblemSize);
    /**
     * Sets up solver kernels.
     * 
     * @param kernel  main constant potential kernel
     */
    void compileKernels(CommonCalcConstantPotentialForceKernel& kernel);
    /**
     * Solves for charges.
     * 
     * @param kernel  main constant potential kernel
     */
    void solveImpl(CommonCalcConstantPotentialForceKernel& kernel);

private:
    /**
     * Solves for charges.
     * 
     * @param kernel  main constant potential kernel
     */
    void ensureValid(CommonCalcConstantPotentialForceKernel& kernel);
};

/**
 * A constant potential solver using direct inversion of the Coulomb matrix.
 * Suitable only when electrode particle positions are fixed.
 */
class CommonConstantPotentialCGSolver : public CommonConstantPotentialSolver {
private:
    Vec3 boxVectors[3];
    bool precondRequested, precondActivated;
    int threadBlockCount, threadBlockSize;
    ComputeArray precondVector;
    ComputeArray q;
    ComputeArray grad;
    ComputeArray projGrad;
    ComputeArray precGrad;
    ComputeArray qStep;
    ComputeArray gradStep;
    ComputeArray grad0;
    ComputeArray qLast;
    ComputeArray blockSums1;
    ComputeArray blockSums2;
    ComputeArray convergedResult;
    ComputeKernel solveInitializeStep1Kernel;
    ComputeKernel solveInitializeStep2Kernel;
    ComputeKernel solveInitializeStep3Kernel;
    ComputeKernel solveLoopStep1Kernel;
    ComputeKernel solveLoopStep2Kernel;
    ComputeKernel solveLoopStep3Kernel;
    ComputeKernel solveLoopStep4Kernel;
    ComputeKernel solveLoopStep5Kernel;
    ComputeEvent convergedDownloadStartEvent;
    ComputeEvent convergedDownloadFinishEvent;
    ComputeQueue convergedDownloadQueue;

public:
    /**
     * Creates a CommonConstantPotentialCGSolver.
     * 
     * @param cc                     the compute context (should be selected)
     * @param numParticles           the number of particles
     * @param numElectrodeParticles  the number of electrode (fluctuating-charge) particles
     * @param paddedProblemSize      the padded number of electrode particles
     * @param precond                whether or not to use a preconditioner
     */
    CommonConstantPotentialCGSolver(ComputeContext& cc, int numParticles, int numElectrodeParticles, int paddedProblemSize, bool precond);
    /**
     * Sets up solver kernels.
     * 
     * @param kernel  main constant potential kernel
     */
    virtual void compileKernels(CommonCalcConstantPotentialForceKernel& kernel);
    /**
     * Solves for charges.
     * 
     * @param kernel  main constant potential kernel
     */
    void solveImpl(CommonCalcConstantPotentialForceKernel& kernel);
    /**
     * Retrieves a list of arrays of initial guess charges to be reordered.
     * 
     * @param arrays  Arrays to be reordered.
     */
    virtual void getGuessChargeArrays(std::vector<ComputeArray*>& arrays);
private:
    /**
     * Solves for charges.
     * 
     * @param kernel  main constant potential kernel
     */
    void ensureValid(CommonCalcConstantPotentialForceKernel& kernel);
};

/**
 * This kernel is invoked by ConstantPotentialForce to calculate the forces acting on the system.
 */
class CommonCalcConstantPotentialForceKernel : public CalcConstantPotentialForceKernel {
    friend class CommonConstantPotentialSolver;
    friend class CommonConstantPotentialMatrixSolver;
    friend class CommonConstantPotentialCGSolver;

public:
    CommonCalcConstantPotentialForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) : CalcConstantPotentialForceKernel(name, platform),
            cc(cc), info(NULL), solver(NULL), hasInitializedKernel(false) {
    }
    ~CommonCalcConstantPotentialForceKernel();
    /**
     * Initialize the kernel.  Subclasses should call this from their initialize() method.
     *
     * @param system       the System this kernel will be applied to
     * @param force        the ConstantPotentialForce this kernel will be used for
     * @param deviceIsCpu  whether the device this calculation is running on is a CPU
     * @param useFixedPointChargeSpreading  whether PME charge spreading should be done in fixed point or floating point
     */
    void commonInitialize(const System& system, const ConstantPotentialForce& force, bool deviceIsCpu, bool useFixedPointChargeSpreading);
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
     * @param force          the ConstantPotentialForce to copy the parameters from
     * @param firstParticle  the index of the first particle whose parameters might have changed
     * @param lastParticle   the index of the last particle whose parameters might have changed
     * @param firstException the index of the first exception whose parameters might have changed
     * @param lastException  the index of the last exception whose parameters might have changed
     * @param firstElectrode the index of the first electrode whose parameters might have changed
     * @param lastElectrode  the index of the last electrode whose parameters might have changed
     */
    void copyParametersToContext(ContextImpl& context, const ConstantPotentialForce& force, int firstParticle, int lastParticle, int firstException, int lastException, int firstElectrode, int lastElectrode);
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
     * Get the charges on all particles.
     *
     * @param context       the context to copy parameters to
     * @param[out] charges  a vector to populate with particle charges
     */
    void getCharges(ContextImpl& context, std::vector<double>& charges);
private:
    void ensureInitialized(ContextImpl& context);
    double doEnergyForces(bool includeForces, bool includeEnergy);
    void initDoDerivatives();
    void doDerivatives(bool init = true);
    void pmeSetup();
    void pmeCompileKernels();
    void initPmeExecute();
    void pmeExecute(bool includeEnergy, bool includeForces, bool includeChargeDerivatives, bool init = true);
    void setKernelInputs(bool includeEnergy, bool includeForces);
    void ensureValidNeighborList();
    class SortTrait : public ComputeSortImpl::SortTrait {
        int getDataSize() const {return 8;}
        int getKeySize() const {return 4;}
        const char* getDataType() const {return "int2";}
        const char* getKeyType() const {return "int";}
        const char* getMinKey() const {return "(-2147483647-1)";}
        const char* getMaxKey() const {return "2147483647";}
        const char* getMaxValue() const {return "make_int2(2147483647, 2147483647)";}
        const char* getSortKey() const {return "value.y";}
    };
    class ForceInfo;
    class ReorderListener;
    class InvalidatePostComputation;
    ComputeContext& cc;
    ForceInfo* info;
    CommonConstantPotentialSolver* solver;
    bool hasInitializedKernel, hasElectrodes, mustUpdateNonElectrodeCharges, mustUpdateElectrodeCharges, pmeShouldSort;
    int numParticles, numElectrodeParticles, numElectrodes, chunkSize, chunkCount, paddedProblemSize;
    ComputeArray charges;
    ComputeArray nonElectrodeCharges;
    ComputeArray electrodeCharges;
    ComputeArray chargeDerivatives;
    ComputeArray chargeDerivativesFixed;
    ComputeArray totalChargeBuffer;
    ComputeArray sysToElec;
    ComputeArray elecToSys;
    ComputeArray sysElec;
    ComputeArray elecElec;
    ComputeArray electrodeParams;
    ComputeArray exclusionScales;
    ComputeArray exceptionScales;
    ComputeArray pmeGrid1;
    ComputeArray pmeGrid2;
    ComputeArray pmeBsplineModuliX;
    ComputeArray pmeBsplineModuliY;
    ComputeArray pmeBsplineModuliZ;
    ComputeArray pmeAtomGridIndex;
    ComputeArray posCellOffsets;
    ComputeSort sort;
    FFT3D fft;
    ComputeKernel updateNonElectrodeChargesKernel;
    ComputeKernel updateElectrodeChargesKernel;
    ComputeKernel getTotalChargeKernel;
    ComputeKernel evaluateSelfEnergyForcesKernel;
    ComputeKernel evaluateDirectDerivativesKernel;
    ComputeKernel finishDerivativesKernel;
    ComputeKernel pmeGridIndexKernel;
    ComputeKernel pmeSpreadChargeKernel;
    ComputeKernel pmeFinishSpreadChargeKernel;
    ComputeKernel pmeConvolutionKernel;
    ComputeKernel pmeEvalEnergyKernel;
    ComputeKernel pmeInterpolateForceKernel;
    ComputeKernel pmeInterpolateChargeDerivativesKernel;
    std::vector<mm_int4> hostPosCellOffsets;
    std::vector<double> setCharges, hostNonElectrodeCharges, hostElectrodeCharges;
    std::vector<std::pair<int, int> > exclusions;
    std::vector<int> exceptions;
    std::vector<int> hostSysToElec, hostElecToSys, hostSysElec, hostElecElec;
    std::vector<mm_double4> hostElectrodeParams;
    std::map<int, int> exceptionIndex;
    std::map<std::string, std::string> pmeDefines;
    ConstantPotentialForce::ConstantPotentialMethod method;
    Vec3 boxVectors[3], externalField;
    double cutoff, ewaldAlpha, chargeTarget, cgErrorTol;
    int maxThreadBlockSize, gridSizeX, gridSizeY, gridSizeZ;
    bool deviceIsCpu, useFixedPointChargeSpreading, usePosqCharges, useChargeConstraint;
    int forceGroup;
    static const int PmeOrder = 5;
    static const double SELF_ALPHA_SCALE, SELF_ETA_SCALE, SELF_TF_SCALE, PLASMA_SCALE;
};

} // namespace OpenMM

#endif /* OPENMM_COMMONCALCCONSTANTPOTENTIALFORCE_H_ */
