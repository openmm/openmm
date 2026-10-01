#ifndef OPENMM_METALCONTEXT_H_
#define OPENMM_METALCONTEXT_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2009-2026 Stanford University and the Authors.      *
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

#include "Metal.hpp"
#include <map>
#include <string>
#include <utility>
#include "MetalArray.h"
#include "MetalIntegrationUtilities.h"
#include "MetalNonbondedUtilities.h"
#include "MetalPlatform.h"
#include "openmm/OpenMMException.h"
#include "openmm/common/BondedUtilities.h"
#include "openmm/common/ComputeContext.h"
#include "openmm/common/ExpressionUtilities.h"
#include "openmm/Kernel.h"

typedef unsigned int tileflags;

namespace OpenMM {

/**
 * This class contains the information associated with a Context by the Metal Platform.
 */

class MetalContext : public ComputeContext {
public:
    class WorkTask;
    class WorkThread;
    class ReorderListener;
    class ForcePreComputation;
    class ForcePostComputation;
    static const int ThreadBlockSize;
    static const int TileSize;
    MetalContext(const System& system, const std::string& precision, MetalPlatform::PlatformData& platformData, MetalContext* originalContext);
    ~MetalContext();
    /**
     * This is called to initialize internal data structures after all Forces in the system
     * have been initialized.
     */
    void initialize();
    /**
     * Get the MTL::Device associated with this object.
     */
    MTL::Device& getDevice() {
        return *device;
    }
    /**
     * Get the PlatformData object this context is part of.
     */
    MetalPlatform::PlatformData& getPlatformData() {
        return platformData;
    }
    /**
     * Get the number of contexts being used for the current simulation.
     * This is relevant when a simulation is parallelized across multiple devices.  In that case,
     * one MetalContext is created for each device.
     */
    int getNumContexts() const {
        return platformData.contexts.size();
    }
    /**
     * Get the index of this context in the list stored in the PlatformData.
     */
    int getContextIndex() const {
        return contextIndex;
    }
    /**
     * Get a list of all contexts being used for the current simulation.
     * This is relevant when a simulation is parallelized across multiple devices.  In that case,
     * one ComputeContext is created for each device.
     */
    std::vector<ComputeContext*> getAllContexts();
    /**
     * Get the ContextImpl is ComputeContext is associated with.
     */
    ContextImpl* getContextImpl() {
        return platformData.context;
    }
    /**
     * Get a workspace used for accumulating energy when a simulation is parallelized across
     * multiple devices.
     */
    double& getEnergyWorkspace();
    /**
     * Create a new ComputeQueue for use with this context.
     */
    ComputeQueue createQueue();
    /**
     * Construct an uninitialized array of the appropriate class for this platform.  The returned
     * value should be created on the heap with the "new" operator.
     */
    MetalArray* createArray();
    /**
     * Construct a ComputeEvent object of the appropriate class for this platform.
     */
    ComputeEvent createEvent();
    /**
     * Construct a ComputeSort object of the appropriate class for this platform.
     * 
     * @param trait      a SortTrait defining the type of data to sort.  It should have been allocated
     *                   on the heap with the "new" operator.  This object takes over ownership of it,
     *                   and deletes it when the ComputeSort is deleted.
     * @param length     the length of the arrays this object will be used to sort
     * @param uniform    whether the input data is expected to follow a uniform or nonuniform
     *                   distribution.  This argument is used only as a hint.  It allows parts
     *                   of the algorithm to be tuned for faster performance on the expected
     *                   distribution.
     */
    ComputeSort createSort(ComputeSortImpl::SortTrait* trait, unsigned int length, bool uniform=true);
    /**
     * Compile source code to create a ComputeProgram.
     *
     * @param source             the source code of the program
     * @param defines            a set of preprocessor definitions (name, value) to define when compiling the program
     */
    ComputeProgram compileProgram(const std::string source, const std::map<std::string, std::string>& defines=std::map<std::string, std::string>());
    /**
     * Convert an array to an MetalArray.  If the argument is already an MetalArray, this simply casts it.
     * If the argument is a ComputeArray that wraps a MetalArray, this returns the wrapped array.  For any
     * other argument, this throws an exception.
     */
    MetalArray& unwrap(ArrayInterface& array) const;
    /**
     * Get the array which contains the force on each atom (represented as three long longs in 64 bit fixed point).
     */
    ArrayInterface& getForce() {
        return longForceBuffer;
    }
    /**
     * The Metal platform does not use floating point force buffers, so this throws an exception.
     */
    ArrayInterface& getFloatForceBuffer() {
        throw OpenMMException("Metal platform does not use floating point force buffers");
    }
    /**
     * The Metal platform does not use floating point force buffers, so this throws an exception.
     */
    ArrayInterface& getForceBuffers() {
        throw OpenMMException("Metal platform does not use floating point force buffers");
    }
    /**
     * Because Apple GPUs use unified memory, there is no such thing as pinned memory.  This just returns a pointer
     * to a block of ordinary memory.  It is guaranteed to be at least as large as any of the arrays returned by methods
     * of this class.
     */
    void* getPinnedBuffer() {
        return pinnedBuffer.data();
    }
    /**
     * Get a shared ThreadPool that code can use to parallelize operations.
     * 
     * Because this object is freely available to all code, care is needed to avoid conflicts.  Only use it
     * from the main thread, and make sure all operations are complete before you invoke any other code that
     * might make use of it
     */
    ThreadPool& getThreadPool() {
        return getPlatformData().threads;
    }
    /**
     * Execute a kernel.
     *
     * @param kernel       the kernel to execute
     * @param arguments    an array of pointers to the kernel arguments
     * @param threads      the maximum number of threads that should be used
     * @param blockSize    the size of each thread block to use
     * @param sharedSize   the amount of dynamic shared memory to allocated for the kernel, in bytes
     */
//    void executeKernel(CUfunction kernel, void** arguments, int workUnits, int blockSize = -1, unsigned int sharedSize = 0);
    /**
     * Compute the largest thread block size that can be used for a kernel that requires a particular amount of
     * shared memory per thread.
     * 
     * @param memory        the number of bytes of shared memory per thread
     */
    int computeThreadBlockSize(double memory) const;
    /**
     * Sum the buffer containing energy.
     */
    double reduceEnergy();
    /**
     * Get the maximum number of threads in a thread block supported by this device.
     */
    int getMaxThreadBlockSize() const {
        return 1024;
    }
    /**
     * Get whether the device being used is a CPU.  In some cases, different algorithms
     * may be more efficient on CPUs and GPUs.
     */
    bool getIsCPU() const {
        return false;
    }
    /**
     * Get the SIMD width of the device being used.
     */
    int getSIMDWidth() const {
        return 32;
    }
    /**
     * Get whether the device being used supports 64 bit atomic operations on global memory.
     */
    bool getSupports64BitGlobalAtomics() const {
        return false;
    }
    /**
     * Get whether the device being used supports double precision math.
     */
    bool getSupportsDoublePrecision() const {
        return false;
    }
    /**
     * Get the number of cores in the GPU.
     */
    int getNumGPUCores() const {
        return numGpuCores;
    }
    /**
     * Get the vectors defining the periodic box.
     */
    void getPeriodicBoxVectors(Vec3& a, Vec3& b, Vec3& c) const {
        a = Vec3(periodicBoxVecX.x, periodicBoxVecX.y, periodicBoxVecX.z);
        b = Vec3(periodicBoxVecY.x, periodicBoxVecY.y, periodicBoxVecY.z);
        c = Vec3(periodicBoxVecZ.x, periodicBoxVecZ.y, periodicBoxVecZ.z);
    }
    /**
     * Set the vectors defining the periodic box.
     */
    void setPeriodicBoxVectors(const Vec3& a, const Vec3& b, const Vec3& c) {
        periodicBoxVecX = mm_double4(a[0], a[1], a[2], 0.0);
        periodicBoxVecY = mm_double4(b[0], b[1], b[2], 0.0);
        periodicBoxVecZ = mm_double4(c[0], c[1], c[2], 0.0);
        periodicBoxVecXFloat = mm_float4((float) a[0], (float) a[1], (float) a[2], 0.0f);
        periodicBoxVecYFloat = mm_float4((float) b[0], (float) b[1], (float) b[2], 0.0f);
        periodicBoxVecZFloat = mm_float4((float) c[0], (float) c[1], (float) c[2], 0.0f);
        periodicBoxSize = mm_double4(a[0], b[1], c[2], 0.0);
        invPeriodicBoxSize = mm_double4(1.0/a[0], 1.0/b[1], 1.0/c[2], 0.0);
        periodicBoxSizeFloat = mm_float4((float) a[0], (float) b[1], (float) c[2], 0.0f);
        invPeriodicBoxSizeFloat = mm_float4(1.0f/(float) a[0], 1.0f/(float) b[1], 1.0f/(float) c[2], 0.0f);
    }
    /**
     * Get the size of the periodic box.
     */
    mm_double4 getPeriodicBoxSize() const {
        return periodicBoxSize;
    }
    /**
     * Get the inverse of the size of the periodic box.
     */
    mm_double4 getInvPeriodicBoxSize() const {
        return invPeriodicBoxSize;
    }
    /**
     * Get the MetalIntegrationUtilities for this context.
     */
    MetalIntegrationUtilities& getIntegrationUtilities() {
        return *integration;
    }
    /**
     * Get the MetalExpressionUtilities for this context.
     */
    ExpressionUtilities& getExpressionUtilities() {
        return *expression;
    }
    /**
     * Get the MetalBondedUtilities for this context.
     */
    BondedUtilities& getBondedUtilities() {
        return *bonded;
    }
    /**
     * Get the MetalNonbondedUtilities for this context.
     */
    MetalNonbondedUtilities& getNonbondedUtilities() {
        return *nonbonded;
    }
    /**
     * Create a new NonbondedUtilities for use with this context.  This should be called
     * only in unusual situations, when a Force needs its own NonbondedUtilities object
     * separate from the standard one.  The caller is responsible for deleting the object
     * when it is no longer needed.
     */
    MetalNonbondedUtilities* createNonbondedUtilities() {
        return new MetalNonbondedUtilities(*this);
    }
    /**
     * Create an object for performing 3D FFTs.  The caller is responsible for deleting
     * the object when it is no longer needed.
     *
     * @param xsize   the first dimension of the data sets on which FFTs will be performed
     * @param ysize   the second dimension of the data sets on which FFTs will be performed
     * @param zsize   the third dimension of the data sets on which FFTs will be performed
     * @param realToComplex  if true, a real-to-complex transform will be done.  Otherwise, it is complex-to-complex.
     */
    FFT3D createFFT(int xsize, int ysize, int zsize, bool realToComplex=false);
    /**
     * This should be called by the Integrator from its own initialize() method.
     * It ensures all contexts are fully initialized.
     */
    void initializeContexts();
    /**
     * Wait until all work that has been queued (kernel executions, asynchronous data transfers, etc.)
     * has been submitted to the device.  This does not mean it has necessarily been completed.
     * Calling this periodically may improve the responsiveness of the computer's GUI, but at the
     * expense of reduced simulation performance.
     */
    void flushQueue();
private:
    MetalPlatform::PlatformData& platformData;
    int contextIndex, numGpuCores;
    bool isLinkedContext;
    mm_float4 periodicBoxVecXFloat, periodicBoxVecYFloat, periodicBoxVecZFloat, periodicBoxSizeFloat, invPeriodicBoxSizeFloat;
    mm_double4 periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ, periodicBoxSize, invPeriodicBoxSize;
    std::string defaultOptimizationOptions;
    std::map<std::string, std::string> compilationDefines;
    MTL::Device* device;
    std::vector<char> pinnedBuffer;
    MetalIntegrationUtilities* integration;
    ExpressionUtilities* expression;
    BondedUtilities* bonded;
    MetalNonbondedUtilities* nonbonded;
};

/**
 * This class exists only for backward compatibility.  Use ComputeContext::WorkTask instead.
 */
class MetalContext::WorkTask : public ComputeContext::WorkTask {
};

/**
 * This class exists only for backward compatibility.  Use ComputeContext::ReorderListener instead.
 */
class MetalContext::ReorderListener : public ComputeContext::ReorderListener {
};

/**
 * This class exists only for backward compatibility.  Use ComputeContext::ForcePreComputation instead.
 */
class MetalContext::ForcePreComputation : public ComputeContext::ForcePreComputation {
};

/**
 * This class exists only for backward compatibility.  Use ComputeContext::ForcePostComputation instead.
 */
class MetalContext::ForcePostComputation : public ComputeContext::ForcePostComputation {
};

} // namespace OpenMM

#endif /*OPENMM_METALCONTEXT_H_*/
