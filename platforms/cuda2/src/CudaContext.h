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
 * Portions copyright (c) 2009-2012 Stanford University and the Authors.      *
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
#include <pthread.h>
#define __CL_ENABLE_EXCEPTIONS
#ifdef _MSC_VER
    // Prevent Windows from defining macros that interfere with other code.
    #define NOMINMAX
#endif
#include <cuda.h>
#include <builtin_types.h>
#include <vector_functions.h>
#include "openmm/internal/windowsExport.h"
#include "CudaPlatform.h"

namespace OpenMM {

class CudaArray;
class CudaForceInfo;
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

class OPENMM_EXPORT CudaContext {
public:
    class WorkTask;
    class WorkThread;
    class ReorderListener;
    static const int ThreadBlockSize;
    static const int TileSize;
    CudaContext(const System& system, int deviceIndex, bool useBlockingSync, const std::string& precision,
            const std::string& compiler, const std::string& tempDir, CudaPlatform::PlatformData& platformData);
    ~CudaContext();
//    /**
//     * This is called to initialize internal data structures after all Forces in the system
//     * have been initialized.
//     */
//    void initialize();
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
     * Get the CUdevice associated with this object.
     */
    CUdevice getDevice() {
        return device;
    }
    /**
     * Get the index of the CUdevice associated with this object.
     */
    int getDeviceIndex() {
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
     * Get the array which contains the position (the xyz components) and charge (the w component) of each atom.
     */
    CudaArray& getPosq() {
        return *posq;
    }
    /**
     * Get the array which contains the velocity (the xyz components) and inverse mass (the w component) of each atom.
     */
    CudaArray& getVelm() {
        return *velm;
    }
//    /**
//     * Get the array which contains the force on each atom.
//     */
//    CudaArray<mm_float4>& getForce() {
//        return *force;
//    }
//    /**
//     * Get the array which contains the buffers in which forces are computed.
//     */
//    CudaArray<mm_float4>& getForceBuffers() {
//        return *forceBuffers;
//    }
//    /**
//     * Get the array which contains a contribution to each force represented as 64 bit fixed point.
//     */
//    CudaArray<cl_long>& getLongForceBuffer() {
//        return *longForceBuffer;
//    }
//    /**
//     * Get the array which contains the buffer in which energy is computed.
//     */
//    CudaArray<cl_float>& getEnergyBuffer() {
//        return *energyBuffer;
//    }
//    /**
//     * Get the array which contains the index of each atom.
//     */
//    CudaArray<cl_int>& getAtomIndex() {
//        return *atomIndex;
//    }
//    /**
//     * Get the number of cells by which the positions are offset.
//     */
//    std::vector<mm_int4>& getPosCellOffsets() {
//        return posCellOffsets;
//    }
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
//    /**
//     * Execute a kernel.
//     *
//     * @param kernel       the kernel to execute
//     * @param workUnits    the maximum number of work units that should be used
//     * @param blockSize    the size of each thread block to use
//     */
//    void executeKernel(cl::Kernel& kernel, int workUnits, int blockSize = -1);
//    /**
//     * Set all elements of an array to 0.
//     */
//    void clearBuffer(CudaArray<float>& array);
//    /**
//     * Set all elements of an array to 0.
//     */
//    void clearBuffer(CudaArray<mm_float4>& array);
//    /**
//     * Set all elements of an array to 0.
//     *
//     * @param memory     the Memory to clear
//     * @param size       the number of float elements in the buffer
//     */
//    void clearBuffer(cl::Memory& memory, int size);
//    /**
//     * Register a buffer that should be automatically cleared (all elements set to 0) at the start of each force or energy computation.
//     *
//     * @param memory     the Memory to clear
//     * @param size       the number of float elements in the buffer
//     */
//    void addAutoclearBuffer(cl::Memory& memory, int size);
//    /**
//     * Clear all buffers that have been registered with addAutoclearBuffer().
//     */
//    void clearAutoclearBuffers();
//    /**
//     * Given a collection of buffers packed into an array, sum them and store
//     * the sum in the first buffer.
//     *
//     * @param array       the array containing the buffers to reduce
//     * @param numBuffers  the number of buffers packed into the array
//     */
//    void reduceBuffer(CudaArray<mm_float4>& array, int numBuffers);
//    /**
//     * Sum the buffesr containing forces.
//     */
//    void reduceForces();
//    /**
//     * Get the current simulation time.
//     */
//    double getTime() {
//        return time;
//    }
//    /**
//     * Set the current simulation time.
//     */
//    void setTime(double t) {
//        time = t;
//    }
//    /**
//     * Get the number of integration steps that have been taken.
//     */
//    int getStepCount() {
//        return stepCount;
//    }
//    /**
//     * Set the number of integration steps that have been taken.
//     */
//    void setStepCount(int steps) {
//        stepCount = steps;
//    }
//    /**
//     * Get the number of times forces or energy has been computed.
//     */
//    int getComputeForceCount() {
//        return computeForceCount;
//    }
//    /**
//     * Set the number of times forces or energy has been computed.
//     */
//    void setComputeForceCount(int count) {
//        computeForceCount = count;
//    }
//    /**
//     * Get the number of atoms.
//     */
//    int getNumAtoms() const {
//        return numAtoms;
//    }
//    /**
//     * Get the number of atoms, rounded up to a multiple of TileSize.  This is the actual size of
//     * most arrays with one element per atom.
//     */
//    int getPaddedNumAtoms() const {
//        return paddedNumAtoms;
//    }
//    /**
//     * Get the number of blocks of TileSize atoms.
//     */
//    int getNumAtomBlocks() const {
//        return numAtomBlocks;
//    }
//    /**
//     * Get the standard number of thread blocks to use when executing kernels.
//     */
//    int getNumThreadBlocks() const {
//        return numThreadBlocks;
//    }
//    /**
//     * Get the number of force buffers.
//     */
//    int getNumForceBuffers() const {
//        return numForceBuffers;
//    }
//    /**
//     * Get the SIMD width of the device being used.
//     */
//    int getSIMDWidth() const {
//        return simdWidth;
//    }
//    /**
//     * Get whether the device being used supports 64 bit atomic operations on global memory.
//     */
//    bool getSupports64BitGlobalAtomics() {
//        return supports64BitGlobalAtomics;
//    }
//    /**
//     * Get whether the device being used supports double precision math.
//     */
//    bool getSupportsDoublePrecision() {
//        return supportsDoublePrecision;
//    }
//    /**
//     * Get the size of the periodic box.
//     */
//    mm_float4 getPeriodicBoxSize() const {
//        return periodicBoxSize;
//    }
//    /**
//     * Set the size of the periodic box.
//     */
//    void setPeriodicBoxSize(double xsize, double ysize, double zsize) {
//        periodicBoxSize = mm_float4((float) xsize, (float) ysize, (float) zsize, 0);
//        invPeriodicBoxSize = mm_float4((float) (1.0/xsize), (float) (1.0/ysize), (float) (1.0/zsize), 0);
//    }
//    /**
//     * Get the inverse of the size of the periodic box.
//     */
//    mm_float4 getInvPeriodicBoxSize() const {
//        return invPeriodicBoxSize;
//    }
//    /**
//     * Get the CudaIntegrationUtilities for this context.
//     */
//    CudaIntegrationUtilities& getIntegrationUtilities() {
//        return *integration;
//    }
//    /**
//     * Get the CudaBondedUtilities for this context.
//     */
//    CudaBondedUtilities& getBondedUtilities() {
//        return *bonded;
//    }
//    /**
//     * Get the CudaNonbondedUtilities for this context.
//     */
//    CudaNonbondedUtilities& getNonbondedUtilities() {
//        return *nonbonded;
//    }
//    /**
//     * Get the thread used by this context for executing parallel computations.
//     */
//    WorkThread& getWorkThread() {
//        return *thread;
//    }
//    /**
//     * Get whether atoms were reordered during the most recent force/energy computation.
//     */
//    bool getAtomsWereReordered() const {
//        return atomsWereReordered;
//    }
//    /**
//     * Set whether atoms were reordered during the most recent force/energy computation.
//     */
//    void setAtomsWereReordered(bool wereReordered) {
//        atomsWereReordered = wereReordered;
//    }
//    /**
//     * Reorder the internal arrays of atoms to try to keep spatially contiguous atoms close
//     * together in the arrays.
//     * 
//     * @param enforcePeriodic    if true, the atom positions may be altered to enforce periodic boundary conditions
//     */
//    void reorderAtoms(bool enforcePeriodic);
//    /**
//     * Add a listener that should be called whenever atoms get reordered.  The CudaContext
//     * assumes ownership of the object, and deletes it when the context itself is deleted.
//     */
//    void addReorderListener(ReorderListener* listener);
//    /**
//     * Get the list of ReorderListeners.
//     */
//    std::vector<ReorderListener*>& getReorderListeners() {
//        return reorderListeners;
//    }
//    /**
//     * Mark that the current molecule definitions (and hence the atom order) may be invalid.
//     * This should be called whenever force field parameters change.  It will cause the definitions
//     * and order to be revalidated the next to reorderAtoms() is called.
//     */
//    void invalidateMolecules();
//    /**
//     * Get whether the current molecule definitions are valid.
//     */
//    bool getMoleculesAreInvalid() {
//        return moleculesInvalid;
//    }
private:
    struct Molecule;
    struct MoleculeGroup;
    class VirtualSiteInfo;
//    void findMoleculeGroups();
//    static void tagAtomsInMolecule(int atom, int molecule, std::vector<int>& atomMolecule, std::vector<std::vector<int> >& atomBonds);
//    /**
//     * Ensure that all molecules marked as "identical" really are identical.  This should be
//     * called whenever force field parameters change.  If necessary, it will rebuild the list
//     * of molecules and resort the atoms.
//     */
//    void validateMolecules();
    static bool hasInitializedCuda;
    const System& system;
    double time;
    CudaPlatform::PlatformData& platformData;
    int deviceIndex;
    int contextIndex;
    int stepCount;
    int computeForceCount;
    int numAtoms;
    int paddedNumAtoms;
    int numAtomBlocks;
    int numThreadBlocks;
//    int numForceBuffers;
//    int simdWidth;
    bool useBlockingSync, useDoublePrecision, accumulateInDouble, contextIsValid, atomsWereReordered, moleculesInvalid;
    std::string compiler, tempDir, gpuArchitecture;
    float4 periodicBoxSize;
    float4 invPeriodicBoxSize;
    std::string defaultOptimizationOptions;
    std::map<std::string, std::string> compilationDefines;
    CUcontext context;
    CUdevice device;
    CUfunction clearBufferKernel;
    CUfunction clearTwoBuffersKernel;
    CUfunction clearThreeBuffersKernel;
    CUfunction clearFourBuffersKernel;
    CUfunction clearFiveBuffersKernel;
    CUfunction clearSixBuffersKernel;
    CUfunction reduceFloat4Kernel;
    CUfunction reduceForcesKernel;
    std::vector<CudaForceInfo*> forces;
    std::vector<Molecule> molecules;
    std::vector<MoleculeGroup> moleculeGroups;
    std::vector<int4> posCellOffsets;
    CudaArray* posq;
    CudaArray* velm;
//    CudaArray<mm_float4>* force;
//    CudaArray<mm_float4>* forceBuffers;
//    CudaArray<cl_long>* longForceBuffer;
//    CudaArray<cl_float>* energyBuffer;
//    CudaArray<cl_int>* atomIndex;
//    std::vector<cl::Memory*> autoclearBuffers;
//    std::vector<int> autoclearBufferSizes;
    std::vector<ReorderListener*> reorderListeners;
//    CudaIntegrationUtilities* integration;
//    CudaBondedUtilities* bonded;
//    CudaNonbondedUtilities* nonbonded;
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
 * Objects that need to know when reordering happens should create a reorderListener
 * and register it by calling addReorderListener().
 */
class CudaContext::ReorderListener {
public:
    virtual void execute() = 0;
    virtual ~ReorderListener() {
    }
};

} // namespace OpenMM

#endif /*OPENMM_CUDACONTEXT_H_*/
