#ifndef OPENMM_OPENCLCONTEXT_H_
#define OPENMM_OPENCLCONTEXT_H_

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
#include <pthread.h>
#define __CL_ENABLE_EXCEPTIONS
#ifdef _MSC_VER
    // Prevent Windows from defining macros that interfere with other code.
    #define NOMINMAX
#endif
#include <cl.hpp>
#include "windowsExportOpenCL.h"
#include "OpenCLPlatform.h"

namespace OpenMM {

class OpenCLArray;
class OpenCLForceInfo;
class OpenCLIntegrationUtilities;
class OpenCLExpressionUtilities;
class OpenCLBondedUtilities;
class OpenCLNonbondedUtilities;
class System;

/**
 * We can't use predefined vector types like cl_float4, since different OpenCL implementations currently define
 * them in incompatible ways.  Hopefully that will be fixed in the future.  In the mean time, we define our own
 * types to represent them on the host.
 */

struct mm_float2 {
    cl_float x, y;
    mm_float2() {
    }
    mm_float2(cl_float x, cl_float y) : x(x), y(y) {
    }
};
struct mm_float4 {
    cl_float x, y, z, w;
    mm_float4() {
    }
    mm_float4(cl_float x, cl_float y, cl_float z, cl_float w) : x(x), y(y), z(z), w(w) {
    }
};
struct mm_float8 {
    cl_float s0, s1, s2, s3, s4, s5, s6, s7;
    mm_float8() {
    }
    mm_float8(cl_float s0, cl_float s1, cl_float s2, cl_float s3, cl_float s4, cl_float s5, cl_float s6, cl_float s7) :
        s0(s0), s1(s1), s2(s2), s3(s3), s4(s4), s5(s5), s6(s6), s7(s7) {
    }
};
struct mm_float16 {
    cl_float s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15;
    mm_float16() {
    }
    mm_float16(cl_float s0, cl_float s1, cl_float s2, cl_float s3, cl_float s4, cl_float s5, cl_float s6, cl_float s7,
            cl_float s8, cl_float s9, cl_float s10, cl_float s11, cl_float s12, cl_float s13, cl_float s14, cl_float s15) :
        s0(s0), s1(s1), s2(s2), s3(s3), s4(s4), s5(s5), s6(s6), s7(s7),
        s8(s8), s9(s9), s10(s10), s11(s11), s12(s12), s13(s13), s14(s14), s15(15) {
    }
};
struct mm_double2 {
    cl_double x, y;
    mm_double2() {
    }
    mm_double2(cl_double x, cl_double y) : x(x), y(y) {
    }
};
struct mm_double4 {
    cl_double x, y, z, w;
    mm_double4() {
    }
    mm_double4(cl_double x, cl_double y, cl_double z, cl_double w) : x(x), y(y), z(z), w(w) {
    }
};
struct mm_ushort2 {
    cl_ushort x, y;
    mm_ushort2() {
    }
    mm_ushort2(cl_ushort x, cl_ushort y) : x(x), y(y) {
    }
};
struct mm_int2 {
    cl_int x, y;
    mm_int2() {
    }
    mm_int2(cl_int x, cl_int y) : x(x), y(y) {
    }
};
struct mm_int4 {
    cl_int x, y, z, w;
    mm_int4() {
    }
    mm_int4(cl_int x, cl_int y, cl_int z, cl_int w) : x(x), y(y), z(z), w(w) {
    }
};
struct mm_int8 {
    cl_int s0, s1, s2, s3, s4, s5, s6, s7;
    mm_int8() {
    }
    mm_int8(cl_int s0, cl_int s1, cl_int s2, cl_int s3, cl_int s4, cl_int s5, cl_int s6, cl_int s7) :
        s0(s0), s1(s1), s2(s2), s3(s3), s4(s4), s5(s5), s6(s6), s7(s7) {
    }
};
struct mm_int16 {
    cl_int s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15;
    mm_int16() {
    }
    mm_int16(cl_int s0, cl_int s1, cl_int s2, cl_int s3, cl_int s4, cl_int s5, cl_int s6, cl_int s7,
            cl_int s8, cl_int s9, cl_int s10, cl_int s11, cl_int s12, cl_int s13, cl_int s14, cl_int s15) :
        s0(s0), s1(s1), s2(s2), s3(s3), s4(s4), s5(s5), s6(s6), s7(s7),
        s8(s8), s9(s9), s10(s10), s11(s11), s12(s12), s13(s13), s14(s14), s15(15) {
    }
};

/**
 * This class contains the information associated with a Context by the OpenCL Platform.  Each OpenCLContext is
 * specific to a particular device, and manages data structures and kernels for that device.  When running a simulation
 * in parallel on multiple devices, there is a separate OpenCLContext for each one.  The list of all contexts is
 * stored in the OpenCLPlatform::PlatformData.
 * <p>
 * In addition, a worker thread is created for each OpenCLContext.  This is used for parallel computations, so that
 * blocking calls to one device will not block other devices.  When only a single device is being used, the worker
 * thread is not used and calculations are performed on the main application thread.
 */

class OPENMM_EXPORT_OPENCL OpenCLContext {
public:
    class WorkTask;
    class WorkThread;
    class ReorderListener;
    class ForcePreComputation;
    class ForcePostComputation;
    static const int ThreadBlockSize;
    static const int TileSize;
    OpenCLContext(const System& system, int platformIndex, int deviceIndex, const std::string& precision, OpenCLPlatform::PlatformData& platformData);
    ~OpenCLContext();
    /**
     * This is called to initialize internal data structures after all Forces in the system
     * have been initialized.
     */
    void initialize();
    /**
     * Add an OpenCLForce to this context.
     */
    void addForce(OpenCLForceInfo* force);
    /**
     * Get the cl::Context associated with this object.
     */
    cl::Context& getContext() {
        return context;
    }
    /**
     * Get the cl::Device associated with this object.
     */
    cl::Device& getDevice() {
        return device;
    }
    /**
     * Get the index of the cl::Device associated with this object.
     */
    int getDeviceIndex() {
        return deviceIndex;
    }
    /**
     * Get the index of the cl::Platform associated with this object.
     */
    int getPlatformIndex() {
        return platformIndex;
    }
    /**
     * Get the PlatformData object this context is part of.
     */
    OpenCLPlatform::PlatformData& getPlatformData() {
        return platformData;
    }
    /**
     * Get the index of this context in the list stored in the PlatformData.
     */
    int getContextIndex() const {
        return contextIndex;
    }
    /**
     * Get the cl::CommandQueue associated with this object.
     */
    cl::CommandQueue& getQueue() {
        return queue;
    }
    /**
     * Get the array which contains the position (the xyz components) and charge (the w component) of each atom.
     */
    OpenCLArray& getPosq() {
        return *posq;
    }
    /**
     * Get the array which contains a correction to the position of each atom.  This only exists if getUseMixedPrecision() returns true.
     */
    OpenCLArray& getPosqCorrection() {
        return *posqCorrection;
    }
    /**
     * Get the array which contains the velocity (the xyz components) and inverse mass (the w component) of each atom.
     */
    OpenCLArray& getVelm() {
        return *velm;
    }
    /**
     * Get the array which contains the force on each atom.
     */
    OpenCLArray& getForce() {
        return *force;
    }
    /**
     * Get the array which contains the buffers in which forces are computed.
     */
    OpenCLArray& getForceBuffers() {
        return *forceBuffers;
    }
    /**
     * Get the array which contains a contribution to each force represented as 64 bit fixed point.
     */
    OpenCLArray& getLongForceBuffer() {
        return *longForceBuffer;
    }
    /**
     * Get the array which contains the buffer in which energy is computed.
     */
    OpenCLArray& getEnergyBuffer() {
        return *energyBuffer;
    }
    /**
     * Get a pointer to a block of pinned memory that can be used for efficient transfers between host and device.
     * This is guaranteed to be at least as large as any of the arrays returned by methods of this class.
     */
    void* getPinnedBuffer() {
        return pinnedMemory;
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
    OpenCLArray& getAtomIndexArray() {
        return *atomIndexDevice;
    }
    /**
     * Get the number of cells by which the positions are offset.
     */
    std::vector<mm_int4>& getPosCellOffsets() {
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
     * Create an OpenCL Program from source code.
     *
     * @param source             the source code of the program
     * @param optimizationFlags  the optimization flags to pass to the OpenCL compiler.  If this is
     *                           omitted, a default set of options will be used
     */
    cl::Program createProgram(const std::string source, const char* optimizationFlags = NULL);
    /**
     * Create an OpenCL Program from source code.
     *
     * @param source             the source code of the program
     * @param defines            a set of preprocessor definitions (name, value) to define when compiling the program
     * @param optimizationFlags  the optimization flags to pass to the OpenCL compiler.  If this is
     *                           omitted, a default set of options will be used
     */
    cl::Program createProgram(const std::string source, const std::map<std::string, std::string>& defines, const char* optimizationFlags = NULL);
    /**
     * Execute a kernel.
     *
     * @param kernel       the kernel to execute
     * @param workUnits    the maximum number of work units that should be used
     * @param blockSize    the size of each thread block to use
     */
    void executeKernel(cl::Kernel& kernel, int workUnits, int blockSize = -1);
    /**
     * Set all elements of an array to 0.
     */
    void clearBuffer(OpenCLArray& array);
    /**
     * Set all elements of an array to 0.
     *
     * @param memory     the Memory to clear
     * @param size       the size of the buffer in bytes
     */
    void clearBuffer(cl::Memory& memory, int size);
    /**
     * Register a buffer that should be automatically cleared (all elements set to 0) at the start of each force or energy computation.
     */
    void addAutoclearBuffer(OpenCLArray& array);
    /**
     * Register a buffer that should be automatically cleared (all elements set to 0) at the start of each force or energy computation.
     *
     * @param memory     the Memory to clear
     * @param size       the size of the buffer in bytes
     */
    void addAutoclearBuffer(cl::Memory& memory, int size);
    /**
     * Clear all buffers that have been registered with addAutoclearBuffer().
     */
    void clearAutoclearBuffers();
    /**
     * Given a collection of floating point buffers packed into an array, sum them and store
     * the sum in the first buffer.
     *
     * @param array       the array containing the buffers to reduce
     * @param numBuffers  the number of buffers packed into the array
     */
    void reduceBuffer(OpenCLArray& array, int numBuffers);
    /**
     * Sum the buffesr containing forces.
     */
    void reduceForces();
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
     * Get the number of force buffers.
     */
    int getNumForceBuffers() const {
        return numForceBuffers;
    }
    /**
     * Get the SIMD width of the device being used.
     */
    int getSIMDWidth() const {
        return simdWidth;
    }
    /**
     * Get whether the device being used supports 64 bit atomic operations on global memory.
     */
    bool getSupports64BitGlobalAtomics() {
        return supports64BitGlobalAtomics;
    }
    /**
     * Get whether the device being used supports double precision math.
     */
    bool getSupportsDoublePrecision() {
        return supportsDoublePrecision;
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
     * Get the size of the periodic box.
     */
    mm_float4 getPeriodicBoxSize() const {
        return periodicBoxSize;
    }
    /**
     * Get the size of the periodic box.
     */
    mm_double4 getPeriodicBoxSizeDouble() const {
        return periodicBoxSizeDouble;
    }
    /**
     * Set the size of the periodic box.
     */
    void setPeriodicBoxSize(double xsize, double ysize, double zsize) {
        periodicBoxSize = mm_float4((float) xsize, (float) ysize, (float) zsize, 0);
        invPeriodicBoxSize = mm_float4((float) (1.0/xsize), (float) (1.0/ysize), (float) (1.0/zsize), 0);
        periodicBoxSizeDouble = mm_double4(xsize, ysize, zsize, 0);
        invPeriodicBoxSizeDouble = mm_double4(1.0/xsize, 1.0/ysize, 1.0/zsize, 0);
    }
    /**
     * Get the inverse of the size of the periodic box.
     */
    mm_float4 getInvPeriodicBoxSize() const {
        return invPeriodicBoxSize;
    }
    /**
     * Get the inverse of the size of the periodic box.
     */
    mm_double4 getInvPeriodicBoxSizeDouble() const {
        return invPeriodicBoxSizeDouble;
    }
    /**
     * Get the OpenCLIntegrationUtilities for this context.
     */
    OpenCLIntegrationUtilities& getIntegrationUtilities() {
        return *integration;
    }
    /**
     * Get the OpenCLExpressionUtilities for this context.
     */
    OpenCLExpressionUtilities& getExpressionUtilities() {
        return *expression;
    }
    /**
     * Get the OpenCLBondedUtilities for this context.
     */
    OpenCLBondedUtilities& getBondedUtilities() {
        return *bonded;
    }
    /**
     * Get the OpenCLNonbondedUtilities for this context.
     */
    OpenCLNonbondedUtilities& getNonbondedUtilities() {
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
    const System& system;
    double time;
    OpenCLPlatform::PlatformData& platformData;
    int deviceIndex;
    int platformIndex;
    int contextIndex;
    int stepCount;
    int computeForceCount;
    int stepsSinceReorder;
    int numAtoms;
    int paddedNumAtoms;
    int numAtomBlocks;
    int numThreadBlocks;
    int numForceBuffers;
    int simdWidth;
    bool supports64BitGlobalAtomics, supportsDoublePrecision, useDoublePrecision, useMixedPrecision, atomsWereReordered;
    mm_float4 periodicBoxSize, invPeriodicBoxSize;
    mm_double4 periodicBoxSizeDouble, invPeriodicBoxSizeDouble;
    std::string defaultOptimizationOptions;
    std::map<std::string, std::string> compilationDefines;
    cl::Context context;
    cl::Device device;
    cl::CommandQueue queue;
    cl::Kernel clearBufferKernel;
    cl::Kernel clearTwoBuffersKernel;
    cl::Kernel clearThreeBuffersKernel;
    cl::Kernel clearFourBuffersKernel;
    cl::Kernel clearFiveBuffersKernel;
    cl::Kernel clearSixBuffersKernel;
    cl::Kernel reduceReal4Kernel;
    cl::Kernel reduceForcesKernel;
    std::vector<OpenCLForceInfo*> forces;
    std::vector<Molecule> molecules;
    std::vector<MoleculeGroup> moleculeGroups;
    std::vector<mm_int4> posCellOffsets;
    cl::Buffer* pinnedBuffer;
    void* pinnedMemory;
    OpenCLArray* posq;
    OpenCLArray* posqCorrection;
    OpenCLArray* velm;
    OpenCLArray* force;
    OpenCLArray* forceBuffers;
    OpenCLArray* longForceBuffer;
    OpenCLArray* energyBuffer;
    OpenCLArray* atomIndexDevice;
    std::vector<int> atomIndex;
    std::vector<cl::Memory*> autoclearBuffers;
    std::vector<int> autoclearBufferSizes;
    std::vector<ReorderListener*> reorderListeners;
    std::vector<ForcePreComputation*> preComputations;
    std::vector<ForcePostComputation*> postComputations;
    OpenCLIntegrationUtilities* integration;
    OpenCLExpressionUtilities* expression;
    OpenCLBondedUtilities* bonded;
    OpenCLNonbondedUtilities* nonbonded;
    WorkThread* thread;
};

struct OpenCLContext::Molecule {
    std::vector<int> atoms;
    std::vector<int> constraints;
    std::vector<std::vector<int> > groups;
};

struct OpenCLContext::MoleculeGroup {
    std::vector<int> atoms;
    std::vector<int> instances;
    std::vector<int> offsets;
};

/**
 * This abstract class defines a task to be executed on the worker thread.
 */
class OpenCLContext::WorkTask {
public:
    virtual void execute() = 0;
    virtual ~WorkTask() {
    }
};

class OpenCLContext::WorkThread {
public:
    struct ThreadData;
    WorkThread();
    ~WorkThread();
    /**
     * Request that a task be executed on the worker thread.  The argument should have been allocated on the
     * heap with the "new" operator.  After its execute() method finishes, the object will be deleted automatically.
     */
    void addTask(OpenCLContext::WorkTask* task);
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
    std::queue<OpenCLContext::WorkTask*> tasks;
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
class OpenCLContext::ReorderListener {
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
class OpenCLContext::ForcePreComputation {
public:
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
class OpenCLContext::ForcePostComputation {
public:
    /**
     * @param includeForce  true if forces should be computed
     * @param includeEnergy true if potential energy should be computed
     * @param groups        a set of bit flags for which force groups to include
     * @return an optional contribution to add to the potential energy.
     */
    virtual double computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) = 0;
};

} // namespace OpenMM

#endif /*OPENMM_OPENCLCONTEXT_H_*/
