#ifndef OPENMM_OPENCLKERNELS_H_
#define OPENMM_OPENCLKERNELS_H_

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

#include "OpenCLPlatform.h"
#include "OpenCLArray.h"
#include "OpenCLContext.h"
#include "OpenCLFFT3D.h"
#include "OpenCLSort.h"
#include "openmm/kernels.h"
#include "openmm/System.h"
#include "openmm/common/CommonKernels.h"

namespace OpenMM {

/**
 * This kernel is invoked at the beginning and end of force and energy computations.  It gives the
 * Platform a chance to clear buffers and do other initialization at the beginning, and to do any
 * necessary work at the end to determine the final results.
 */
class OpenCLCalcForcesAndEnergyKernel : public CalcForcesAndEnergyKernel {
public:
    OpenCLCalcForcesAndEnergyKernel(std::string name, const Platform& platform, OpenCLContext& cl) : CalcForcesAndEnergyKernel(name, platform), cl(cl) {
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
   OpenCLContext& cl;
};

/**
 * This kernel is invoked by NonbondedForce to calculate the forces acting on the system.
 */
class OpenCLCalcNonbondedForceKernel : public CalcNonbondedForceKernel {
public:
    OpenCLCalcNonbondedForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, const System& system) : CalcNonbondedForceKernel(name, platform),
            hasInitializedKernel(false), cl(cl), sort(NULL), fft(NULL), dispersionFft(NULL), pmeio(NULL), usePmeQueue(false) {
    }
    ~OpenCLCalcNonbondedForceKernel();
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
    class SortTrait : public OpenCLSort::SortTrait {
        int getDataSize() const {return 8;}
        int getKeySize() const {return 4;}
        const char* getDataType() const {return "int2";}
        const char* getKeyType() const {return "int";}
        const char* getMinKey() const {return "INT_MIN";}
        const char* getMaxKey() const {return "INT_MAX";}
        const char* getMaxValue() const {return "(int2) (INT_MAX, INT_MAX)";}
        const char* getSortKey() const {return "value.y";}
    };
    class ForceInfo;
    class PmeIO;
    class PmePreComputation;
    class PmePostComputation;
    class SyncQueuePreComputation;
    class SyncQueuePostComputation;
    OpenCLContext& cl;
    ForceInfo* info;
    bool hasInitializedKernel;
    OpenCLArray charges;
    OpenCLArray sigmaEpsilon;
    OpenCLArray exceptionParams;
    OpenCLArray exclusionAtoms;
    OpenCLArray exclusionParams;
    OpenCLArray baseParticleParams;
    OpenCLArray baseExceptionParams;
    OpenCLArray particleParamOffsets;
    OpenCLArray exceptionParamOffsets;
    OpenCLArray particleOffsetIndices;
    OpenCLArray exceptionOffsetIndices;
    OpenCLArray globalParams;
    OpenCLArray cosSinSums;
    OpenCLArray pmeGrid1;
    OpenCLArray pmeGrid2;
    OpenCLArray pmeBsplineModuliX;
    OpenCLArray pmeBsplineModuliY;
    OpenCLArray pmeBsplineModuliZ;
    OpenCLArray pmeDispersionBsplineModuliX;
    OpenCLArray pmeDispersionBsplineModuliY;
    OpenCLArray pmeDispersionBsplineModuliZ;
    OpenCLArray pmeBsplineTheta;
    OpenCLArray pmeAtomRange;
    OpenCLArray pmeAtomGridIndex;
    OpenCLArray pmeEnergyBuffer;
    OpenCLSort* sort;
    cl::CommandQueue pmeQueue;
    cl::Event pmeSyncEvent;
    OpenCLFFT3D* fft;
    OpenCLFFT3D* dispersionFft;
    Kernel cpuPme;
    PmeIO* pmeio;
    SyncQueuePostComputation* syncQueue;
    cl::Kernel computeParamsKernel, computeExclusionParamsKernel;
    cl::Kernel ewaldSumsKernel;
    cl::Kernel ewaldForcesKernel;
    cl::Kernel pmeAtomRangeKernel;
    cl::Kernel pmeDispersionAtomRangeKernel;
    cl::Kernel pmeZIndexKernel;
    cl::Kernel pmeDispersionZIndexKernel;
    cl::Kernel pmeGridIndexKernel;
    cl::Kernel pmeDispersionGridIndexKernel;
    cl::Kernel pmeSpreadChargeKernel;
    cl::Kernel pmeDispersionSpreadChargeKernel;
    cl::Kernel pmeFinishSpreadChargeKernel;
    cl::Kernel pmeDispersionFinishSpreadChargeKernel;
    cl::Kernel pmeConvolutionKernel;
    cl::Kernel pmeDispersionConvolutionKernel;
    cl::Kernel pmeEvalEnergyKernel;
    cl::Kernel pmeDispersionEvalEnergyKernel;
    cl::Kernel pmeInterpolateForceKernel;
    cl::Kernel pmeDispersionInterpolateForceKernel;
    std::map<std::string, std::string> pmeDefines;
    std::vector<std::pair<int, int> > exceptionAtoms;
    std::vector<std::string> paramNames;
    std::vector<double> paramValues;
    double ewaldSelfEnergy, dispersionCoefficient, alpha, dispersionAlpha;
    int gridSizeX, gridSizeY, gridSizeZ;
    int dispersionGridSizeX, dispersionGridSizeY, dispersionGridSizeZ;
    bool hasCoulomb, hasLJ, usePmeQueue, doLJPME, usePosqCharges, recomputeParams, hasOffsets;
    NonbondedMethod nonbondedMethod;
    static const int PmeOrder = 5;
};

/**
 * This kernel is invoked by CustomCVForce to calculate the forces acting on the system and the energy of the system.
 */
class OpenCLCalcCustomCVForceKernel : public CommonCalcCustomCVForceKernel {
public:
    OpenCLCalcCustomCVForceKernel(std::string name, const Platform& platform, ComputeContext& cc) : CommonCalcCustomCVForceKernel(name, platform, cc) {
    }
    ComputeContext& getInnerComputeContext(ContextImpl& innerContext) {
        return *reinterpret_cast<OpenCLPlatform::PlatformData*>(innerContext.getPlatformData())->contexts[0];
    }
};

/**
 * This kernel is invoked by ATMForce to calculate the forces acting on the system and the energy of the system.
 */
class OpenCLCalcATMForceKernel : public CommonCalcATMForceKernel {
public:
    OpenCLCalcATMForceKernel(std::string name, const Platform& platform, ComputeContext& cc) : CommonCalcATMForceKernel(name, platform, cc) {
    }
    ComputeContext& getInnerComputeContext(ContextImpl& innerContext) {
        return *reinterpret_cast<OpenCLPlatform::PlatformData*>(innerContext.getPlatformData())->contexts[0];
    }
};

  
} // namespace OpenMM

#endif /*OPENMM_OPENCLKERNELS_H_*/
