#ifndef OPENMM_COMPUTECONTEXT_H_
#define OPENMM_COMPUTECONTEXT_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2019-2024 Stanford University and the Authors.      *
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

#ifdef _MSC_VER
    // Prevent Windows from defining macros that interfere with other code.
    #define NOMINMAX
#endif
#include "openmm/common/ArrayInterface.h"
#include "openmm/common/BondedUtilities.h"
#include "openmm/common/ComputeEvent.h"
#include "openmm/common/ComputeForceInfo.h"
#include "openmm/common/ComputeProgram.h"
#include "openmm/common/ComputeVectorTypes.h"
#include "openmm/common/IntegrationUtilities.h"
#include "openmm/common/NonbondedUtilities.h"
#include "openmm/Vec3.h"
#include <pthread.h>
#include <map>
#include <queue>
#include <string>
#include <vector>

namespace OpenMM {

class ExpressionUtilities;
class System;
class ThreadPool;

/**
 * This abstract class defines the interface by which platforms compile and execute
 * kernels.  It also manages the arrays use for storing standard information, like
 * positions and forces.
 */

class OPENMM_EXPORT_COMMON ComputeContext {
public:
    class WorkTask;
    class WorkThread;
    class ReorderListener;
    class ForcePreComputation;
    class ForcePostComputation;
    static const int ThreadBlockSize;
    static const int TileSize;
    ComputeContext(const System& system);
    virtual ~ComputeContext();
    /**
     * Add a ComputeForceInfo to this context.  Force kernels call this during initialization
     * to provide information about particular forces.
     */
    virtual void addForce(ComputeForceInfo* force);
    /**
     * Get all ComputeForceInfos that have been added to this context.
     */
    std::vector<ComputeForceInfo*>& getForceInfos() {
        return forces;
    }
    /**
     * Request that the context provide at least a particular number of force buffers.
     * This is only meaningful for devices that do not support 64 bit atomic operations.
     * On other devices, this will typically have no effect.  Force kernels should call
     * this during initialization.
     */
    virtual void requestForceBuffers(int minBuffers) {
    }
    /**
     * Set this as the current context for the calling thread.  This should be called before
     * doing any computation when you do not know what other code has just been executing on
     * the thread.  Platforms that rely on binding contexts to threads (such as CUDA) need to
     * implement this.
     * 
     * @deprecated It is recommended to use pushAsCurrent() and popAsCurrent() instead, or even better to create a ContextSelector.
     * This provides better interoperability with other libraries that use CUDA and create
     * their own contexts.
     */
    virtual void setAsCurrent() {
    }
    /**
     * Set this as the current context for the calling thread, maintaining any previous context
     * on a stack.  This should be called before doing any computation when you do not know what
     * other code has just been executing on the thread.  It must be paired with popAsCurrent()
     * when you are done to restore the previous context.  Alternatively, you can create a
     * ContextSelector object to automate this for a block of code.
     * 
     * Platforms that rely on binding contexts to threads (such as CUDA) need to implement this.
     */
    virtual void pushAsCurrent() {
    }
    /**
     * Restore a previous context that was replaced by pushAsCurrent().  Platforms that rely on binding
     * contexts to threads (such as CUDA) need to implement this.
     */
    virtual void popAsCurrent() {
    }
    /**
     * Get the number of contexts being used for the current simulation.
     * This is relevant when a simulation is parallelized across multiple devices.  In that case,
     * one ComputeContext is created for each device.
     */
    virtual int getNumContexts() const = 0;
    /**
     * Get the index of this context in the list of ones being used for the current simulation.
     * This is relevant when a simulation is parallelized across multiple devices.  In that case,
     * one ComputeContext is created for each device.
     */
    virtual int getContextIndex() const = 0;
    /**
     * Get a list of all contexts being used for the current simulation.
     * This is relevant when a simulation is parallelized across multiple devices.  In that case,
     * one ComputeContext is created for each device.
     */
    virtual std::vector<ComputeContext*> getAllContexts() = 0;
    /**
     * Get a workspace used for accumulating energy when a simulation is parallelized across
     * multiple devices.
     */
    virtual double& getEnergyWorkspace() = 0;
    /**
     * Construct an uninitialized array of the appropriate class for this platform.  The returned
     * value should be created on the heap with the "new" operator.
     */
    virtual ArrayInterface* createArray() = 0;
    /**
     * Construct a ComputeEvent object of the appropriate class for this platform.
     */
    virtual ComputeEvent createEvent() = 0;
    /**
     * Compile source code to create a ComputeProgram.
     *
     * @param source             the source code of the program
     * @param defines            a set of preprocessor definitions (name, value) to define when compiling the program
     */
    virtual ComputeProgram compileProgram(const std::string source, const std::map<std::string, std::string>& defines=std::map<std::string, std::string>()) = 0;
    /**
     * Compute the largest thread block size that can be used for a kernel that requires a particular amount of
     * shared memory per thread.
     * 
     * @param memory        the number of bytes of shared memory per thread
     */
    virtual int computeThreadBlockSize(double memory) const = 0;
    /**
     * Set all elements of an array to 0.
     */
    virtual void clearBuffer(ArrayInterface& array) = 0;
    /**
     * Register an array that should be automatically cleared (all elements set to 0) at the start of each force or energy computation.
     */
    virtual void addAutoclearBuffer(ArrayInterface& array) = 0;
    /**
     * Get whether the device being used is a CPU.  In some cases, different algorithms
     * may be more efficient on CPUs and GPUs.
     */
    virtual bool getIsCPU() const = 0;
    /**
     * Get the SIMD width of the device being used.
     */
    virtual int getSIMDWidth() const = 0;
    /**
     * Get whether the device being used supports 64 bit atomic operations on global memory.
     */
    virtual bool getSupports64BitGlobalAtomics() const = 0;
    /**
     * Get whether the device being used supports double precision math.
     */
    virtual bool getSupportsDoublePrecision() const = 0;
    /**
     * Get whether double precision is being used.
     */
    virtual bool getUseDoublePrecision() const = 0;
    /**
     * Get whether mixed precision is being used.
     */
    virtual bool getUseMixedPrecision() const = 0;
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
    long long getStepCount() {
        return stepCount;
    }
    /**
     * Set the number of integration steps that have been taken.
     */
    void setStepCount(long long steps) {
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
     * 
     * Calling this method might or might not actually change the atom order.  It uses
     * internal heuristics to decide when and how often to update the order.  If you
     * want to guarantee that reordering will definitely be done, call forceReorder() before
     * calling this.
     */
    void reorderAtoms();
    /**
     * Calling this method guarantees that the next call to reorderAtoms() will actually
     * perform reordering.
     */
    void forceReorder();
    /**
     * Add a listener that should be called whenever atoms get reordered.  The OpenCLContext
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
     * The OpenCLContext assumes ownership of the object, and deletes it when the context itself is deleted.
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
     * The OpenCLContext assumes ownership of the object, and deletes it when the context itself is deleted.
     */
    void addPostComputation(ForcePostComputation* computation);
    /**
     * Get the list of ForcePostComputations.
     */
    std::vector<ForcePostComputation*>& getPostComputations() {
        return postComputations;
    }
    /**
     * Get the flag that marks whether the current force evaluation is valid.
     */
    bool getForcesValid() const {
        return forcesValid;
    }
    /**
     * Get the flag that marks whether the current force evaluation is valid.
     */
    void setForcesValid(bool valid) {
        forcesValid = valid;
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
    virtual int getNumAtomBlocks() const = 0;
    /**
     * Get the standard number of thread blocks to use when executing kernels.
     */
    virtual int getNumThreadBlocks() const = 0;
    /**
     * Get the maximum number of threads in a thread block supported by this device.
     */
    virtual int getMaxThreadBlockSize() const = 0;
    /**
     * Get the array which contains the position (the xyz components) and charge (the w component) of each atom.
     */
    virtual ArrayInterface& getPosq() = 0;
    /**
     * Get the array which contains a correction to the position of each atom.  This only exists if getUseMixedPrecision() returns true.
     */
    virtual ArrayInterface& getPosqCorrection() = 0;
    /**
     * Get the array which contains the velocity (the xyz components) and inverse mass (the w component) of each atom.
     */
    virtual ArrayInterface& getVelm() = 0;
    /**
     * On devices that do not support 64 bit atomics, this returns an array containing buffers of type real4 in which
     * forces can be accumulated.  On platforms that do not use floating point force buffers, this will throw an exception.
     */
    virtual ArrayInterface& getForceBuffers() = 0;
    /**
     * Get the array which contains a contribution to each force represented as a real4.  On platforms that do not use
     * floating point force buffers, this will throw an exception.
     */
    virtual ArrayInterface& getFloatForceBuffer() = 0;
    /**
     * Get the array which contains a contribution to each force represented as 64 bit fixed point.
     */
    virtual ArrayInterface& getLongForceBuffer() = 0;
    /**
     * Get the array which contains the buffer in which energy is computed.
     */
    virtual ArrayInterface& getEnergyBuffer() = 0;
    /**
     * Get the array which contains the buffer in which derivatives of the energy with respect to parameters are computed.
     */
    virtual ArrayInterface& getEnergyParamDerivBuffer() = 0;
    /**
     * Get a pointer to a block of pinned memory that can be used for asynchronous transfers between host and device.
     * This is guaranteed to be at least as large as any of the arrays returned by methods of this class.
     * 
     * Because this buffer is freely available to all code, care is needed to avoid conflicts.  Only access this
     * buffer from the main thread, and make sure all transfers are complete before you invoke any other code that
     * might make use of it
     */
    virtual void* getPinnedBuffer() = 0;
    /**
     * Get a shared ThreadPool that code can use to parallelize operations.
     * 
     * Because this object is freely available to all code, care is needed to avoid conflicts.  Only use it
     * from the main thread, and make sure all operations are complete before you invoke any other code that
     * might make use of it
     */
    virtual ThreadPool& getThreadPool() = 0;
    /**
     * Get the host-side vector which contains the index of each atom.
     */
    const std::vector<int>& getAtomIndex() const {
        return atomIndex;
    }
    /**
     * Set the vector which contains the index of each atom.
     */
    void setAtomIndex(std::vector<int>& index);
    /**
     * Get the array which contains the index of each atom.
     */
    virtual ArrayInterface& getAtomIndexArray() = 0;
    /**
     * Get the number of cells by which the positions are offset.
     */
    std::vector<mm_int4>& getPosCellOffsets() {
        return posCellOffsets;
    }
    /**
     * Set the number of cells by which the positions are offset.
     */
    void setPosCellOffsets(std::vector<mm_int4>& offsets) {
        posCellOffsets = offsets;
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
     * Convert a number to a string in a format suitable for including in a kernel.
     * This takes into account whether the context uses single or double precision.
     * If mixedIsDouble is true, a double precision constant will also be produced
     * in mixed precision mode.
     */
    std::string doubleToString(double value, bool mixedIsDouble=false) const;
    /**
     * Convert a number to a string in a format suitable for including in a kernel.
     */
    std::string intToString(int value) const;
    /**
     * Get whether the periodic box is triclinic.
     */
    virtual bool getBoxIsTriclinic() const = 0;
    /**
     * Get the vectors defining the periodic box.
     */
    virtual void getPeriodicBoxVectors(Vec3& a, Vec3& b, Vec3& c) const = 0;
    /**
     * Set the vectors defining the periodic box.
     */
    virtual void setPeriodicBoxVectors(const Vec3& a, const Vec3& b, const Vec3& c) = 0; 
    /**
     * Get the IntegrationUtilities for this context.
     */
    virtual IntegrationUtilities& getIntegrationUtilities() = 0;
    /**
     * Get the ExpressionUtilities for this context.
     */
    virtual ExpressionUtilities& getExpressionUtilities() = 0;
    /**
     * Get the BondedUtilities for this context.
     */
    virtual BondedUtilities& getBondedUtilities() = 0;
    /**
     * Get the NonbondedUtilities for this context.
     */
    virtual NonbondedUtilities& getNonbondedUtilities() = 0;
    /**
     * Create a new NonbondedUtilities for use with this context.  This should be called
     * only in unusual situations, when a Force needs its own NonbondedUtilities object
     * separate from the standard one.  The caller is responsible for deleting the object
     * when it is no longer needed.
     */
    virtual NonbondedUtilities* createNonbondedUtilities() = 0;
    /**
     * Get the smallest legal size for a dimension of the grid.
     */
    virtual int findLegalFFTDimension(int minimum);
    /**
     * This should be called by the Integrator from its own initialize() method.
     * It ensures all contexts are fully initialized.
     */
    virtual void initializeContexts() = 0;
    /**
     * Get the thread used by this context for executing parallel computations.
     */
    WorkThread& getWorkThread() {
        return *thread;
    }
    /**
     * Get the names of all parameters with respect to which energy derivatives are computed.
     */
    virtual const std::vector<std::string>& getEnergyParamDerivNames() const = 0;
    /**
     * Get a workspace data structure used for accumulating the values of derivatives of the energy
     * with respect to parameters.
     */
    virtual std::map<std::string, double>& getEnergyParamDerivWorkspace() = 0;
    /**
     * Register that the derivative of potential energy with respect to a context parameter
     * will need to be calculated.  If this is called multiple times for a single parameter,
     * it is only added to the list once.
     * 
     * @param param    the name of the parameter to add
     */
    virtual void addEnergyParameterDerivative(const std::string& param) = 0;
    /**
     * Mark that the current molecule definitions (and hence the atom order) may be invalid.
     * This should be called whenever force field parameters change.  It will cause the definitions
     * and order to be revalidated.
     * 
     * If you know which force has changed, calling the alternate form that takes a ComputeForceInfo
     * is more efficient.
     */
    void invalidateMolecules();
    /**
     * Mark that the current molecule definitions from one particular force (and hence the atom order)
     * may be invalid.  This should be called whenever force field parameters change.  It will cause the
     * definitions and order to be revalidated.
     */
    bool invalidateMolecules(ComputeForceInfo* force, bool checkAtoms=true, bool checkGroups=true);
    /**
     * Make sure the current atom order is valid, based on the forces.  If not, perform reordering
     * to generate a new valid order.  This method is only needed in very unusual situations.
     */
    void validateAtomOrder();
    /**
     * Wait until all work that has been queued (kernel executions, asynchronous data transfers, etc.)
     * has been submitted to the device.  This does not mean it has necessarily been completed.
     * Calling this periodically may improve the responsiveness of the computer's GUI, but at the
     * expense of reduced simulation performance.
     */
    virtual void flushQueue() = 0;
protected:
    struct Molecule;
    struct MoleculeGroup;
    class VirtualSiteInfo;
    void findMoleculeGroups();
    void resetAtomOrder();
    /**
     * This is the internal implementation of reorderAtoms(), templatized by the numerical precision in use.
     */
    template <class Real, class Real4, class Mixed, class Mixed4>
    void reorderAtomsImpl();
    const System& system;
    double time;
    int numAtoms, paddedNumAtoms, computeForceCount, stepsSinceReorder;
    long long stepCount;
    bool forceNextReorder, atomsWereReordered, forcesValid;
    std::vector<ComputeForceInfo*> forces;
    std::vector<Molecule> molecules;
    std::vector<MoleculeGroup> moleculeGroups;
    std::vector<int> atomIndex;
    std::vector<mm_int4> posCellOffsets;
    std::vector<ReorderListener*> reorderListeners;
    std::vector<ForcePreComputation*> preComputations;
    std::vector<ForcePostComputation*> postComputations;
    WorkThread* thread;
};

struct ComputeContext::Molecule {
    std::vector<int> atoms;
    std::vector<int> constraints;
    std::vector<std::vector<int> > groups;
};

struct ComputeContext::MoleculeGroup {
    std::vector<int> atoms;
    std::vector<int> instances;
    std::vector<int> offsets;
};

/**
 * This abstract class defines a task to be executed on the worker thread.
 */
class OPENMM_EXPORT_COMMON ComputeContext::WorkTask {
public:
    virtual void execute() = 0;
    virtual ~WorkTask() {
    }
};

class OPENMM_EXPORT_COMMON ComputeContext::WorkThread {
public:
    struct ThreadData;
    WorkThread();
    ~WorkThread();
    /**
     * Request that a task be executed on the worker thread.  The argument should have been allocated on the
     * heap with the "new" operator.  After its execute() method finishes, the object will be deleted automatically.
     */
    void addTask(ComputeContext::WorkTask* task);
    /**
     * Get whether the worker thread is idle, waiting for a task to be added.
     */
    bool isWaiting();
    /**
     * Get whether the worker thread has exited.
     */
    bool isFinished();
    /**
     * Get whether the thread invoking this method is the worker thread.
     */
    bool isCurrentThread();
    /**
     * Block until all tasks have finished executing and the worker thread is idle.
     */
    void flush();
private:
    std::queue<ComputeContext::WorkTask*> tasks;
    bool waiting, finished, threwException;
    OpenMMException stashedException;
    pthread_mutex_t queueLock;
    pthread_cond_t waitForTaskCondition, queueEmptyCondition;
    pthread_t thread;
};

/**
 * This abstract class defines a function to be executed whenever atoms get reordered.
 * Objects that need to know when reordering happens should create a ReorderListener
 * and register it by calling addReorderListener().
 */
class OPENMM_EXPORT_COMMON ComputeContext::ReorderListener {
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
class OPENMM_EXPORT_COMMON ComputeContext::ForcePreComputation {
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
class OPENMM_EXPORT_COMMON ComputeContext::ForcePostComputation {
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

#endif /*OPENMM_COMPUTECONTEXT_H_*/
