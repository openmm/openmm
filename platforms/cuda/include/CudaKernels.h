#ifndef OPENMM_CUDAKERNELS_H_
#define OPENMM_CUDAKERNELS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2021 Stanford University and the Authors.      *
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

#include "CudaPlatform.h"
#include "CudaArray.h"
#include "CudaContext.h"
#include "CudaFFT3D.h"
#include "CudaSort.h"
#include "openmm/kernels.h"
#include "openmm/System.h"
#include "openmm/common/CommonKernels.h"
#include <cufft.h>

namespace OpenMM {

/**
 * This abstract class defines an interface for code that can compile CUDA kernels.  This allows a plugin to take advantage of runtime compilation
 * when running on recent versions of CUDA.
 */
class CudaCompilerKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CudaCompilerKernel";
    }
    CudaCompilerKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Compile a kernel to PTX.
     *
     * @param source     the source code for the kernel
     * @param options    the flags to be passed to the compiler
     * @param cu         the CudaContext for which the kernel is being compiled
     */
    virtual std::string createModule(const std::string& source, const std::string& flags, CudaContext& cu) = 0;
    /**
     * Get the maximum architecture version the compiler supports.
     */
    virtual int getMaxSupportedArchitecture() const = 0;
};

/**
 * This kernel is invoked at the beginning and end of force and energy computations.  It gives the
 * Platform a chance to clear buffers and do other initialization at the beginning, and to do any
 * necessary work at the end to determine the final results.
 */
class CudaCalcForcesAndEnergyKernel : public CalcForcesAndEnergyKernel {
public:
    CudaCalcForcesAndEnergyKernel(std::string name, const Platform& platform, CudaContext& cu) : CalcForcesAndEnergyKernel(name, platform), cu(cu) {
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
   CudaContext& cu;
};

/**
 * This kernel provides methods for setting and retrieving various state data: time, positions,
 * velocities, and forces.
 */
class CudaUpdateStateDataKernel : public UpdateStateDataKernel {
public:
    CudaUpdateStateDataKernel(std::string name, const Platform& platform, CudaContext& cu) : UpdateStateDataKernel(name, platform), cu(cu) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     */
    void initialize(const System& system);
    /**
     * Get the current time (in picoseconds).
     *
     * @param context    the context in which to execute this kernel
     */
    double getTime(const ContextImpl& context) const;
    /**
     * Set the current time (in picoseconds).
     *
     * @param context    the context in which to execute this kernel
     */
    void setTime(ContextImpl& context, double time);
    /**
     * Get the positions of all particles.
     *
     * @param positions  on exit, this contains the particle positions
     */
    void getPositions(ContextImpl& context, std::vector<Vec3>& positions);
    /**
     * Set the positions of all particles.
     *
     * @param positions  a vector containg the particle positions
     */
    void setPositions(ContextImpl& context, const std::vector<Vec3>& positions);
    /**
     * Get the velocities of all particles.
     *
     * @param velocities  on exit, this contains the particle velocities
     */
    void getVelocities(ContextImpl& context, std::vector<Vec3>& velocities);
    /**
     * Set the velocities of all particles.
     *
     * @param velocities  a vector containg the particle velocities
     */
    void setVelocities(ContextImpl& context, const std::vector<Vec3>& velocities);
    /**
     * Get the current forces on all particles.
     *
     * @param forces  on exit, this contains the forces
     */
    void getForces(ContextImpl& context, std::vector<Vec3>& forces);
    /**
     * Get the current derivatives of the energy with respect to context parameters.
     *
     * @param derivs  on exit, this contains the derivatives
     */
    void getEnergyParameterDerivatives(ContextImpl& context, std::map<std::string, double>& derivs);
    /**
     * Get the current periodic box vectors.
     *
     * @param a      on exit, this contains the vector defining the first edge of the periodic box
     * @param b      on exit, this contains the vector defining the second edge of the periodic box
     * @param c      on exit, this contains the vector defining the third edge of the periodic box
     */
    void getPeriodicBoxVectors(ContextImpl& context, Vec3& a, Vec3& b, Vec3& c) const;
    /**
     * Set the current periodic box vectors.
     *
     * @param a      the vector defining the first edge of the periodic box
     * @param b      the vector defining the second edge of the periodic box
     * @param c      the vector defining the third edge of the periodic box
     */
    void setPeriodicBoxVectors(ContextImpl& context, const Vec3& a, const Vec3& b, const Vec3& c);
    /**
     * Create a checkpoint recording the current state of the Context.
     * 
     * @param stream    an output stream the checkpoint data should be written to
     */
    void createCheckpoint(ContextImpl& context, std::ostream& stream);
    /**
     * Load a checkpoint that was written by createCheckpoint().
     * 
     * @param stream    an input stream the checkpoint data should be read from
     */
    void loadCheckpoint(ContextImpl& context, std::istream& stream);
private:
    CudaContext& cu;
};

/**
 * This kernel is invoked by NonbondedForce to calculate the forces acting on the system.
 */
class CudaCalcNonbondedForceKernel : public CalcNonbondedForceKernel {
public:
    CudaCalcNonbondedForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcNonbondedForceKernel(name, platform),
            cu(cu), hasInitializedFFT(false), sort(NULL), dispersionFft(NULL), fft(NULL), pmeio(NULL), usePmeStream(false) {
    }
    ~CudaCalcNonbondedForceKernel();
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
     * @param context    the context to copy parameters to
     * @param force      the NonbondedForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const NonbondedForce& force);
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
     * Get the dispersion parameters being used for the dispersion term in LJPME.
     * 
     * @param alpha   the separation parameter
     * @param nx      the number of grid points along the X axis
     * @param ny      the number of grid points along the Y axis
     * @param nz      the number of grid points along the Z axis
     */
    void getLJPMEParameters(double& alpha, int& nx, int& ny, int& nz) const;
private:
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
    class ForceInfo;
    class PmeIO;
    class PmePreComputation;
    class PmePostComputation;
    class SyncStreamPreComputation;
    class SyncStreamPostComputation;
    CudaContext& cu;
    ForceInfo* info;
    bool hasInitializedFFT;
    CudaArray charges;
    CudaArray sigmaEpsilon;
    CudaArray exceptionParams;
    CudaArray exclusionAtoms;
    CudaArray exclusionParams;
    CudaArray baseParticleParams;
    CudaArray baseExceptionParams;
    CudaArray particleParamOffsets;
    CudaArray exceptionParamOffsets;
    CudaArray particleOffsetIndices;
    CudaArray exceptionOffsetIndices;
    CudaArray globalParams;
    CudaArray cosSinSums;
    CudaArray pmeGrid1;
    CudaArray pmeGrid2;
    CudaArray pmeBsplineModuliX;
    CudaArray pmeBsplineModuliY;
    CudaArray pmeBsplineModuliZ;
    CudaArray pmeDispersionBsplineModuliX;
    CudaArray pmeDispersionBsplineModuliY;
    CudaArray pmeDispersionBsplineModuliZ;
    CudaArray pmeAtomGridIndex;
    CudaArray pmeEnergyBuffer;
    CudaSort* sort;
    Kernel cpuPme;
    PmeIO* pmeio;
    CUstream pmeStream;
    CUevent pmeSyncEvent, paramsSyncEvent;
    CudaFFT3D* fft;
    cufftHandle fftForward;
    cufftHandle fftBackward;
    CudaFFT3D* dispersionFft;
    cufftHandle dispersionFftForward;
    cufftHandle dispersionFftBackward;
    CUfunction computeParamsKernel, computeExclusionParamsKernel;
    CUfunction ewaldSumsKernel;
    CUfunction ewaldForcesKernel;
    CUfunction pmeGridIndexKernel;
    CUfunction pmeDispersionGridIndexKernel;
    CUfunction pmeSpreadChargeKernel;
    CUfunction pmeDispersionSpreadChargeKernel;
    CUfunction pmeFinishSpreadChargeKernel;
    CUfunction pmeDispersionFinishSpreadChargeKernel;
    CUfunction pmeEvalEnergyKernel;
    CUfunction pmeEvalDispersionEnergyKernel;
    CUfunction pmeConvolutionKernel;
    CUfunction pmeDispersionConvolutionKernel;
    CUfunction pmeInterpolateForceKernel;
    CUfunction pmeInterpolateDispersionForceKernel;
    std::vector<std::pair<int, int> > exceptionAtoms;
    std::vector<std::string> paramNames;
    std::vector<double> paramValues;
    double ewaldSelfEnergy, dispersionCoefficient, alpha, dispersionAlpha;
    int interpolateForceThreads;
    int gridSizeX, gridSizeY, gridSizeZ;
    int dispersionGridSizeX, dispersionGridSizeY, dispersionGridSizeZ;
    bool hasCoulomb, hasLJ, usePmeStream, useCudaFFT, doLJPME, usePosqCharges, recomputeParams, hasOffsets;
    NonbondedMethod nonbondedMethod;
    static const int PmeOrder = 5;
};

/**
 * This kernel is invoked by CustomCVForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcCustomCVForceKernel : public CommonCalcCustomCVForceKernel {
public:
    CudaCalcCustomCVForceKernel(std::string name, const Platform& platform, ComputeContext& cc) : CommonCalcCustomCVForceKernel(name, platform, cc) {
    }
    ComputeContext& getInnerComputeContext(ContextImpl& innerContext) {
        return *reinterpret_cast<CudaPlatform::PlatformData*>(innerContext.getPlatformData())->contexts[0];
    }
};

} // namespace OpenMM

#endif /*OPENMM_CUDAKERNELS_H_*/

