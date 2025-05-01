#ifndef OPENMM_COMMONCALCNONBONDEDFORCEKERNEL_H_
#define OPENMM_COMMONCALCNONBONDEDFORCEKERNEL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2025 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
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

/**
 * This kernel is invoked by NonbondedForce to calculate the forces acting on the system.
 */
class CommonCalcNonbondedForceKernel : public CalcNonbondedForceKernel {
public:
    CommonCalcNonbondedForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) : CalcNonbondedForceKernel(name, platform),
            hasInitializedKernel(false), cc(cc), pmeio(NULL) {
    }
    ~CommonCalcNonbondedForceKernel();
    /**
     * Initialize the kernel.  Subclasses should call this from their initialize() method.
     *
     * @param system       the System this kernel will be applied to
     * @param force        the NonbondedForce this kernel will be used for
     * @param usePmeQueue  whether to perform PME on a separate queue
     * @param deviceIsCpu  whether the device this calculation is running on is a CPU
     * @param useFixedPointChargeSpreading  whether PME charge spreading should be done in fixed point or floating point
     * @param useCpuPme    whether to perform the PME reciprocal space calculation on the CPU
     */
    void commonInitialize(const System& system, const NonbondedForce& force, bool usePmeQueue, bool deviceIsCpu, bool useFixedPointChargeSpreading, bool useCpuPme);
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
     * Get the parameters being used for the dispersion term in LJPME.
     *
     * @param alpha   the separation parameter
     * @param nx      the number of grid points along the X axis
     * @param ny      the number of grid points along the Y axis
     * @param nz      the number of grid points along the Z axis
     */
    void getLJPMEParameters(double& alpha, int& nx, int& ny, int& nz) const;
private:
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
    class PmeIO;
    class PmePreComputation;
    class PmePostComputation;
    class SyncQueuePreComputation;
    class SyncQueuePostComputation;
    ComputeContext& cc;
    ForceInfo* info;
    bool hasInitializedKernel;
    ComputeArray charges;
    ComputeArray sigmaEpsilon;
    ComputeArray exceptionParams;
    ComputeArray exclusionAtoms;
    ComputeArray exclusionParams;
    ComputeArray baseParticleParams;
    ComputeArray baseExceptionParams;
    ComputeArray particleParamOffsets;
    ComputeArray exceptionParamOffsets;
    ComputeArray particleOffsetIndices;
    ComputeArray exceptionOffsetIndices;
    ComputeArray globalParams;
    ComputeArray cosSinSums;
    ComputeArray pmeGrid1;
    ComputeArray pmeGrid2;
    ComputeArray pmeBsplineModuliX;
    ComputeArray pmeBsplineModuliY;
    ComputeArray pmeBsplineModuliZ;
    ComputeArray pmeDispersionBsplineModuliX;
    ComputeArray pmeDispersionBsplineModuliY;
    ComputeArray pmeDispersionBsplineModuliZ;
    ComputeArray pmeAtomGridIndex;
    ComputeArray pmeEnergyBuffer;
    ComputeArray chargeBuffer;
    ComputeSort sort;
    ComputeQueue pmeQueue;
    ComputeEvent pmeSyncEvent, paramsSyncEvent;
    FFT3D fft, dispersionFft;
    Kernel cpuPme;
    PmeIO* pmeio;
    SyncQueuePostComputation* syncQueue;
    ComputeKernel computeParamsKernel, computeExclusionParamsKernel, computePlasmaCorrectionKernel;
    ComputeKernel ewaldSumsKernel, ewaldForcesKernel;
    ComputeKernel pmeGridIndexKernel, pmeDispersionGridIndexKernel;
    ComputeKernel pmeSpreadChargeKernel, pmeDispersionSpreadChargeKernel;
    ComputeKernel pmeFinishSpreadChargeKernel, pmeDispersionFinishSpreadChargeKernel;
    ComputeKernel pmeConvolutionKernel, pmeDispersionConvolutionKernel;
    ComputeKernel pmeEvalEnergyKernel, pmeDispersionEvalEnergyKernel;
    ComputeKernel pmeInterpolateForceKernel, pmeDispersionInterpolateForceKernel;
    std::map<std::string, std::string> pmeDefines;
    std::vector<std::pair<int, int> > exceptionAtoms;
    std::vector<std::string> paramNames;
    std::vector<double> paramValues;
    std::map<int, int> exceptionIndex;
    double ewaldSelfEnergy, dispersionCoefficient, alpha, dispersionAlpha, totalCharge;
    int gridSizeX, gridSizeY, gridSizeZ;
    int dispersionGridSizeX, dispersionGridSizeY, dispersionGridSizeZ;
    bool usePmeQueue, deviceIsCpu, useFixedPointChargeSpreading, useCpuPme;
    bool hasCoulomb, hasLJ, doLJPME, usePosqCharges, recomputeParams, hasOffsets;
    NonbondedMethod nonbondedMethod;
    static const int PmeOrder = 5;
};

} // namespace OpenMM

#endif /*OPENMM_COMMONCALCNONBONDEDFORCEKERNEL_H_*/
