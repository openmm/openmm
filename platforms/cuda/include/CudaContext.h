#ifndef OPENMM_CUDACONTEXT_H_
#define OPENMM_CUDACONTEXT_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2013 Stanford University and the Authors.      *
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

#include <map>
#include <queue>
#include <string>
#define __CL_ENABLE_EXCEPTIONS
#ifdef _MSC_VER
    // Prevent Windows from defining macros that interfere with other code.
    #define NOMINMAX
#endif
#include <pthread.h>
#include <cuda.h>
#include <builtin_types.h>
#include <vector_functions.h>
#include "windowsExportCuda.h"
#include "CudaPlatform.h"

typedef unsigned int tileflags;

namespace OpenMM {

class CudaArray;
class CudaForceInfo;
class CudaExpressionUtilities;
class CudaIntegrationUtilities;
class CudaBondedUtilities;
class CudaNonbondedUtilities;
class System;

/**
 * This class contains the information associated with a Context by the CUDA Platform.  Each CudaContext is
 * specific to a particular device, and manages data structures and kernels for that device.  When running a simulation
 * in parallel on multiple devices, there is a separate CudaContext for each one.  The list of all contexts is
 * stored in the CudaPlatform::PlatformData.
 * <p>
 * In addition, a worker thread is created for each CudaContext.  This is used for parallel computations, so that
 * blocking calls to one device will not block other devices.  When only a single device is being used, the worker
 * thread is not used and calculations are performed on the main application thread.
 */

class OPENMM_EXPORT_CUDA CudaContext {
public:
    class WorkTask;
    class WorkThread;
    class ReorderListener;
    class ForcePreComputation;
    class ForcePostComputation;
    static const int ThreadBlockSize;
    static const int TileSize;
    CudaContext(const System& system, int deviceIndex, bool useBlockingSync, const std::string& precision,
            const std::string& compiler, const std::string& tempDir, const std::string& hostCompiler, CudaPlatform::PlatformData& platformData);
    ~CudaContext();
    /**
     * This is called to initialize internal data structures after all Forces in the system
     * have been initialized.
     */
    void initialize();
    /**
     * Add a CudaForce to this context.
     */
    void addForce(CudaForceInfo* force);
    /**
     * Get the CUcontext associated with this object.
     */
    CUcontext getContext() {
        return context;
    }
    /**
     * Get whether the CUcontext associated with this object is currently a valid contex.
     */
    bool getContextIsValid() const {
        return contextIsValid;
    }
    /**
     * Set the CUcontext associated with this object to be the current context.  If the context is not
     * valid, this returns without doing anything.
     */
    void setAsCurrent();
    /**
     * Get the CUdevice associated with this object.
     */
    CUdevice getDevice() {
        return device;
    }
    /**
     * Get the compute capability of the device associated with this object.
     */
    double getComputeCapability() const {
        return computeCapability;
    }
    /**
     * Get the index of the CUdevice associated with this object.
     */
    int getDeviceIndex() const {
        return deviceIndex;
    }
    /**
     * Get the PlatformData object this context is part of.
     */
    CudaPlatform::PlatformData& getPlatformData() {
        return platformData;
    }
    /**
     * Get the index of this context in the list stored in the PlatformData.
     */
    int getContextIndex() const {
        return contextIndex;
    }
    /**
     * Get the stream currently being used for execution.
     */
    CUstream getCurrentStream();
    /**
     * Set the stream to use for execution.
     */
    void setCurrentStream(CUstream stream);
    /**
     * Reset the context to using the default stream for execution.
     */
    void restoreDefaultStream();
    /**
     * Get the array which contains the position (the xyz components) and charge (the w component) of each atom.
     */
    CudaArray& getPosq() {
        return *posq;
    }
    /**
     * Get the array which contains a correction to the position of each atom.  This only exists if getUseMixedPrecision() returns true.
     */
    CudaArray& getPosqCorrection() {
        return *posqCorrection;
    }
    /**
     * Get the array which contains the velocity (the xyz components) and inverse mass (the w component) of each atom.
     */
    CudaArray& getVelm() {
        return *velm;
    }
    /**
     * Get the array which contains the force on each atom (represented as three long longs in 64 bit fixed point).
     */
    CudaArray& getForce() {
        return *force;
    }
    /**
     * Get the array which contains the buffer in which energy is computed.
     */
    CudaArray& getEnergyBuffer() {
        return *energyBuffer;
    }
    /**
     * Get a pointer to a block of pinned memory that can be used for efficient transfers between host and device.
     * This is guaranteed to be at least as large as any of the arrays returned by methods of this class.
     */
    void* getPinnedBuffer() {
        return pinnedBuffer;
    }
    /**
     * Get the host-side vector which contains the index of each atom.
     */
    const std::vector<int>& getAtomIndex() const {
        return atomIndex;
    }
    /**
     * Get the array which contains the index of each atom.
     */
    CudaArray& getAtomIndexArray() {
        return *atomIndexDevice;
    }
    /**
     * Get the number of cells by which the positions are offset.
     */
    std::vector<int4>& getPosCellOffsets() {
        return posCellOffsets;
    }
    /**
     * Replace all occurrences of a list of substrings.
     *
     * @param input   a string to process
     * @param replacements a set of strings that should be replaced with new strings wherever they appear in the input string
     * @return a new string produced by performing the replacements
     */
    std::string replaceStrings(const std::string& input, const std::map<std::string, std::string>& replacements) const;
    /**
     * Create a CUDA module from source code.
     *
     * @param source             the source code of the module
     * @param optimizationFlags  the optimization flags to pass to the CUDA compiler.  If this is
     *                           omitted, a default set of options will be used
     */
    CUmodule createModule(const std::string source, const char* optimizationFlags = NULL);
    /**
     * Create a CUDA module from source code.
     *
     * @param source             the source code of the module
     * @param defines            a set of preprocessor definitions (name, value) to define when compiling the program
     * @param optimizationFlags  the optimization flags to pass to the CUDA compiler.  If this is
     *                           omitted, a default set of options will be used
     */
    CUmodule createModule(const std::string source, const std::map<std::string, std::string>& defines, const char* optimizationFlags = NULL);
    /**
     * Get a kernel from a CUDA module.
     *
     * @param module    the module to get the kernel from
     * @param name      the name of the kernel to get
     */
    CUfunction getKernel(CUmodule& module, const std::string& name);
    /**
     * Execute a kernel.
     *
     * @param kernel       the kernel to execute
     * @param arguments    an array of pointers to the kernel arguments
     * @param threads      the maximum number of threads that should be used
     * @param blockSize    the size of each thread block to use
     * @param sharedSize   the amount of dynamic shared memory to allocated for the kernel, in bytes
     */
    void executeKernel(CUfunction kernel, void** arguments, int workUnits, int blockSize = -1, unsigned int sharedSize = 0);
    /**
     * Compute the largest thread block size that can be used for a kernel that requires a particular amount of
     * shared memory per thread.
     * 
     * @param memory        the number of bytes of shared memory per thread
     * @param preferShared  whether the kernel is set to prefer shared memory over cache
     */
    int computeThreadBlockSize(double memory, bool preferShared=true) const;
    /**
     * Set all elements of an array to 0.
     */
    void clearBuffer(CudaArray& array);
    /**
     * Set all elements of an array to 0.
     *
     * @param memory     the memory to clear
     * @param size       the size of the buffer in bytes
     */
    void clearBuffer(CUdeviceptr memory, int size);
    /**
     * Register a buffer that should be automatically cleared (all elements set to 0) at the start of each force or energy computation.
     */
    void addAutoclearBuffer(CudaArray& array);
    /**
     * Register a buffer that should be automatically cleared (all elements set to 0) at the start of each force or energy computation.
     *
     * @param memory     the memory to clear
     * @param size       the size of the buffer in bytes
     */
    void addAutoclearBuffer(CUdeviceptr memory, int size);
    /**
     * Clear all buffers that have been registered with addAutoclearBuffer().
     */
    void clearAutoclearBuffers();
    /**
     * Get the current simulation time.
     */
    double getTime() {
        return time;
    }
    /**
     * Set the current simulation time.
     */
    void setTime(double t) {
        time = t;
    }
    /**
     * Get the number of integration steps that have been taken.
     */
    int getStepCount() {
        return stepCount;
    }
    /**
     * Set the number of integration steps that have been taken.
     */
    void setStepCount(int steps) {
        stepCount = steps;
    }
    /**
     * Get the number of times forces or energy has been computed.
     */
    int getComputeForceCount() {
        return computeForceCount;
    }
    /**
     * Set the number of times forces or energy has been computed.
     */
    void setComputeForceCount(int count) {
        computeForceCount = count;
    }
    /**
     * Get the number of time steps since the atoms were reordered.
     */
    int getStepsSinceReorder() const {
        return stepsSinceReorder;
    }
    /**
     * Set the number of time steps since the atoms were reordered.
     */
    void setStepsSinceReorder(int steps) {
        stepsSinceReorder = steps;
    }
    /**
     * Get the number of atoms.
     */
    int getNumAtoms() const {
        return numAtoms;
    }
    /**
     * Get the number of atoms, rounded up to a multiple of TileSize.  This is the actual size of
     * most arrays with one element per atom.
     */
    int getPaddedNumAtoms() const {
        return paddedNumAtoms;
    }
    /**
     * Get the number of blocks of TileSize atoms.
     */
    int getNumAtomBlocks() const {
        return numAtomBlocks;
    }
    /**
     * Get the standard number of thread blocks to use when executing kernels.
     */
    int getNumThreadBlocks() const {
        return numThreadBlocks;
    }
    /**
     * Get whether double precision is being used.
     */
    bool getUseDoublePrecision() {
        return useDoublePrecision;
    }
    /**
     * Get whether mixed precision is being used.
     */
    bool getUseMixedPrecision() {
        return useMixedPrecision;
    }
    /**
     * Convert a number to a string in a format suitable for including in a kernel.
     * This takes into account whether the context uses single or double precision.
     */
    std::string doubleToString(double value);
    /**
     * Convert a number to a string in a format suitable for including in a kernel.
     */
    std::string intToString(int value);
    /**
     * Convert a CUDA result code to the corresponding string description.
     */
    static std::string getErrorString(CUresult result);
    
    /**
     * Get the size of the periodic box.
     */
    double4 getPeriodicBoxSize() const {
        return periodicBoxSize;
    }
    /**
     * Set the size of the periodic box.
     */
    void setPeriodicBoxSize(double xsize, double ysize, double zsize) {
        periodicBoxSize = make_double4(xsize, ysize, zsize, 0.0);
        invPeriodicBoxSize = make_double4(1.0/xsize, 1.0/ysize, 1.0/zsize, 0.0);
        periodicBoxSizeFloat = make_float4((float) xsize, (float) ysize, (float) zsize, 0.0f);
        invPeriodicBoxSizeFloat = make_float4(1.0f/(float) xsize, 1.0f/(float) ysize, 1.0f/(float) zsize, 0.0f);
    }
    /**
     * Get the inverse of the size of the periodic box.
     */
    double4 getInvPeriodicBoxSize() const {
        return invPeriodicBoxSize;
    }
    /**
     * Get a pointer to the size of the periodic box, represented as either a float4 or double4 depending on
     * this context's precision.  This value is suitable for passing to kernels as an argument.
     */
    void* getPeriodicBoxSizePointer() {
        return (useDoublePrecision ? reinterpret_cast<void*>(&periodicBoxSize) : reinterpret_cast<void*>(&periodicBoxSizeFloat));
    }
    /**
     * Get a pointer to the inverse of the size of the periodic box, represented as either a float4 or double4 depending on
     * this context's precision.  This value is suitable for passing to kernels as an argument.
     */
    void* getInvPeriodicBoxSizePointer() {
        return (useDoublePrecision ? reinterpret_cast<void*>(&invPeriodicBoxSize) : reinterpret_cast<void*>(&invPeriodicBoxSizeFloat));
    }
    /**
     * Get the CudaIntegrationUtilities for this context.
     */
    CudaIntegrationUtilities& getIntegrationUtilities() {
        return *integration;
    }
    /**
     * Get the CudaExpressionUtilities for this context.
     */
    CudaExpressionUtilities& getExpressionUtilities() {
        return *expression;
    }
    /**
     * Get the CudaBondedUtilities for this context.
     */
    CudaBondedUtilities& getBondedUtilities() {
        return *bonded;
    }
    /**
     * Get the CudaNonbondedUtilities for this context.
     */
    CudaNonbondedUtilities& getNonbondedUtilities() {
        return *nonbonded;
    }
    /**
     * Get the thread used by this context for executing parallel computations.
     */
    WorkThread& getWorkThread() {
        return *thread;
    }
    /**
     * Get whether atoms were reordered during the most recent force/energy computation.
     */
    bool getAtomsWereReordered() const {
        return atomsWereReordered;
    }
    /**
     * Set whether atoms were reordered during the most recent force/energy computation.
     */
    void setAtomsWereReordered(bool wereReordered) {
        atomsWereReordered = wereReordered;
    }
    /**
     * Reorder the internal arrays of atoms to try to keep spatially contiguous atoms close
     * together in the arrays.
     */
    void reorderAtoms();
    /**
     * Add a listener that should be called whenever atoms get reordered.  The CudaContext
     * assumes ownership of the object, and deletes it when the context itself is deleted.
     */
    void addReorderListener(ReorderListener* listener);
    /**
     * Get the list of ReorderListeners.
     */
    std::vector<ReorderListener*>& getReorderListeners() {
        return reorderListeners;
    }
    /**
     * Add a pre-computation that should be called at the very start of force and energy evaluations.
     * The CudaContext assumes ownership of the object, and deletes it when the context itself is deleted.
     */
    void addPreComputation(ForcePreComputation* computation);
    /**
     * Get the list of ForcePreComputations.
     */
    std::vector<ForcePreComputation*>& getPreComputations() {
        return preComputations;
    }
    /**
     * Add a post-computation that should be called at the very end of force and energy evaluations.
     * The CudaContext assumes ownership of the object, and deletes it when the context itself is deleted.
     */
    void addPostComputation(ForcePostComputation* computation);
    /**
     * Get the list of ForcePostComputations.
     */
    std::vector<ForcePostComputation*>& getPostComputations() {
        return postComputations;
    }
    /**
     * Mark that the current molecule definitions (and hence the atom order) may be invalid.
     * This should be called whenever force field parameters change.  It will cause the definitions
     * and order to be revalidated.
     */
    void invalidateMolecules();
private:
    struct Molecule;
    struct MoleculeGroup;
    class VirtualSiteInfo;
    void findMoleculeGroups();
    /**
     * Ensure that all molecules marked as "identical" really are identical.  This should be
     * called whenever force field parameters change.  If necessary, it will rebuild the list
     * of molecules and resort the atoms.
     */
    void validateMolecules();
    /**
     * This is the internal implementation of reorderAtoms(), templatized by the numerical precision in use.
     */
    template <class Real, class Real4, class Mixed, class Mixed4>
    void reorderAtomsImpl();
    static bool hasInitializedCuda;
    const System& system;
    double time, computeCapability;
    CudaPlatform::PlatformData& platformData;
    int deviceIndex;
    int contextIndex;
    int stepCount;
    int computeForceCount;
    int stepsSinceReorder;
    int numAtoms;
    int paddedNumAtoms;
    int numAtomBlocks;
    int numThreadBlocks;
    bool useBlockingSync, useDoublePrecision, useMixedPrecision, contextIsValid, atomsWereReordered;
    std::string compiler, tempDir, cacheDir, gpuArchitecture;
    float4 periodicBoxSizeFloat, invPeriodicBoxSizeFloat;
    double4 periodicBoxSize, invPeriodicBoxSize;
    std::string defaultOptimizationOptions;
    std::map<std::string, std::string> compilationDefines;
    CUcontext context;
    CUdevice device;
    CUstream currentStream;
    CUfunction clearBufferKernel;
    CUfunction clearTwoBuffersKernel;
    CUfunction clearThreeBuffersKernel;
    CUfunction clearFourBuffersKernel;
    CUfunction clearFiveBuffersKernel;
    CUfunction clearSixBuffersKernel;
    std::vector<CudaForceInfo*> forces;
    std::vector<Molecule> molecules;
    std::vector<MoleculeGroup> moleculeGroups;
    std::vector<int4> posCellOffsets;
    void* pinnedBuffer;
    CudaArray* posq;
    CudaArray* posqCorrection;
    CudaArray* velm;
    CudaArray* force;
    CudaArray* energyBuffer;
    CudaArray* atomIndexDevice;
    std::vector<int> atomIndex;
    std::vector<CUdeviceptr> autoclearBuffers;
    std::vector<int> autoclearBufferSizes;
    std::vector<ReorderListener*> reorderListeners;
    std::vector<ForcePreComputation*> preComputations;
    std::vector<ForcePostComputation*> postComputations;
    CudaIntegrationUtilities* integration;
    CudaExpressionUtilities* expression;
    CudaBondedUtilities* bonded;
    CudaNonbondedUtilities* nonbonded;
    WorkThread* thread;
};

struct CudaContext::Molecule {
    std::vector<int> atoms;
    std::vector<int> constraints;
    std::vector<std::vector<int> > groups;
};

struct CudaContext::MoleculeGroup {
    std::vector<int> atoms;
    std::vector<int> instances;
    std::vector<int> offsets;
};

/**
 * This abstract class defines a task to be executed on the worker thread.
 */
class CudaContext::WorkTask {
public:
    virtual void execute() = 0;
    virtual ~WorkTask() {
    }
};

class CudaContext::WorkThread {
public:
    struct ThreadData;
    WorkThread();
    ~WorkThread();
    /**
     * Request that a task be executed on the worker thread.  The argument should have been allocated on the
     * heap with the "new" operator.  After its execute() method finishes, the object will be deleted automatically.
     */
    void addTask(CudaContext::WorkTask* task);
    /**
     * Get whether the worker thread is idle, waiting for a task to be added.
     */
    bool isWaiting();
    /**
     * Get whether the worker thread has exited.
     */
    bool isFinished();
    /**
     * Block until all tasks have finished executing and the worker thread is idle.
     */
    void flush();
private:
    std::queue<CudaContext::WorkTask*> tasks;
    bool waiting, finished;
    pthread_mutex_t queueLock;
    pthread_cond_t waitForTaskCondition, queueEmptyCondition;
    pthread_t thread;
};

/**
 * This abstract class defines a function to be executed whenever atoms get reordered.
 * Objects that need to know when reordering happens should create a ReorderListener
 * and register it by calling addReorderListener().
 */
class CudaContext::ReorderListener {
public:
    virtual void execute() = 0;
    virtual ~ReorderListener() {
    }
};

/**
 * This abstract class defines a function to be executed at the very beginning of force and
 * energy evaluation, before any other calculation has been done.  It is useful for operations
 * that need to be performed at a nonstandard point in the process.  After creating a
 * ForcePreComputation, register it by calling addForcePreComputation().
 */
class CudaContext::ForcePreComputation {
public:
    virtual ~ForcePreComputation() {
    }
    /**
     * @param includeForce  true if forces should be computed
     * @param includeEnergy true if potential energy should be computed
     * @param groups        a set of bit flags for which force groups to include
     */
    virtual void computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) = 0;
};

/**
 * This abstract class defines a function to be executed at the very end of force and
 * energy evaluation, after all other calculations have been done.  It is useful for operations
 * that need to be performed at a nonstandard point in the process.  After creating a
 * ForcePostComputation, register it by calling addForcePostComputation().
 */
class CudaContext::ForcePostComputation {
public:
    virtual ~ForcePostComputation() {
    }
    /**
     * @param includeForce  true if forces should be computed
     * @param includeEnergy true if potential energy should be computed
     * @param groups        a set of bit flags for which force groups to include
     * @return an optional contribution to add to the potential energy.
     */
    virtual double computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) = 0;
};

} // namespace OpenMM

#endif /*OPENMM_CUDACONTEXT_H_*/
