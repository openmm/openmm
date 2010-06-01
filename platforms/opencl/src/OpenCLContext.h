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
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
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
#include <string>
#define __CL_ENABLE_EXCEPTIONS
#ifdef _MSC_VER
    // Prevent Windows from defining macros that interfere with other code.
    #define NOMINMAX
#endif
#include <cl.hpp>
#include "openmm/internal/windowsExport.h"

namespace OpenMM {

template <class T>
class OpenCLArray;
class OpenCLForceInfo;
class OpenCLIntegrationUtilities;
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

/**
 * This class contains the information associated with a Context by the OpenCL Platform.
 */

class OPENMM_EXPORT OpenCLContext {
public:
    static const int ThreadBlockSize = 64;
    static const int TileSize = 32;
    OpenCLContext(int numParticles, int deviceIndex);
    ~OpenCLContext();
    /**
     * This is called to initialize internal data structures after all Forces in the system
     * have been initialized.
     */
    void initialize(const System& system);
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
     * Get the cl::CommandQueue associated with this object.
     */
    cl::CommandQueue& getQueue() {
        return queue;
    }
    /**
     * Get the array which contains the position and charge of each atom.
     */
    OpenCLArray<mm_float4>& getPosq() {
        return *posq;
    }
    /**
     * Get the array which contains the velocity and inverse mass of each atom.
     */
    OpenCLArray<mm_float4>& getVelm() {
        return *velm;
    }
    /**
     * Get the array which contains the force on each atom.
     */
    OpenCLArray<mm_float4>& getForce() {
        return *force;
    }
    /**
     * Get the array which contains the buffers in which forces are computed.
     */
    OpenCLArray<mm_float4>& getForceBuffers() {
        return *forceBuffers;
    }
    /**
     * Get the array which contains the buffer in which energy is computed.
     */
    OpenCLArray<cl_float>& getEnergyBuffer() {
        return *energyBuffer;
    }
    /**
     * Get the array which contains the index of each atom.
     */
    OpenCLArray<cl_int>& getAtomIndex() {
        return *atomIndex;
    }
    /**
     * Get the number of cells by which the positions are offset.
     */
    std::vector<mm_int4>& getPosCellOffsets() {
        return posCellOffsets;
    }
    /**
     * Load OpenCL source code from a file in the kernels directory.
     */
    std::string loadSourceFromFile(const std::string& filename) const;
    /**
     * Load OpenCL source code from a file in the kernels directory.
     *
     * @param filename     the file to load
     * @param replacements a set of strings that should be replaced with new strings wherever they appear in the
     */
    std::string loadSourceFromFile(const std::string& filename, const std::map<std::string, std::string>& replacements) const;
    /**
     * Replace all occurance of a list of substrings.
     *
     * @param input   a string to process
     * @param replacements a set of strings that should be replaced with new strings wherever they appear in the input string
     * @return a new string produced by performing the replacements
     */
    std::string replaceStrings(const std::string& input, const std::map<std::string, std::string>& replacements) const;
    /**
     * Create an OpenCL Program from source code.
     */
    cl::Program createProgram(const std::string source);
    /**
     * Create an OpenCL Program from source code.
     *
     * @param defines    a set of preprocessor definitions (name, value) to define when compiling the program
     */
    cl::Program createProgram(const std::string source, const std::map<std::string, std::string>& defines);
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
    void clearBuffer(OpenCLArray<float>& array);
    /**
     * Set all elements of an array to 0.
     */
    void clearBuffer(OpenCLArray<mm_float4>& array);
    /**
     * Set all elements of an array to 0.
     *
     * @param memory     the Memory to clear
     * @param size       the number of float elements in the buffer
     */
    void clearBuffer(cl::Memory& memory, int size);
    /**
     * Register a buffer that should be automatically cleared (all elements set to 0) at the start of each force or energy computation.
     *
     * @param memory     the Memory to clear
     * @param size       the number of float elements in the buffer
     */
    void addAutoclearBuffer(cl::Memory& memory, int size);
    /**
     * Clear all buffers that have been registered with addAutoclearBuffer().
     */
    void clearAutoclearBuffers();
    /**
     * Given a collection of buffers packed into an array, sum them and store
     * the sum in the first buffer.
     *
     * @param array       the array containing the buffers to reduce
     * @param numBuffers  the number of buffers packed into the array
     */
    void reduceBuffer(OpenCLArray<mm_float4>& array, int numBuffers);
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
     * Get the size of the periodic box.
     */
    mm_float4 getPeriodicBoxSize() const {
        return periodicBoxSize;
    }
    /**
     * Set the size of the periodic box.
     */
    void setPeriodicBoxSize(double xsize, double ysize, double zsize) {
        periodicBoxSize = mm_float4((float) xsize, (float) ysize, (float) zsize, 0);
        invPeriodicBoxSize = mm_float4((float) (1.0/xsize), (float) (1.0/ysize), (float) (1.0/zsize), 0);
    }
    /**
     * Get the inverse of the size of the periodic box.
     */
    mm_float4 getInvPeriodicBoxSize() const {
        return invPeriodicBoxSize;
    }
    /**
     * Get the OpenCLIntegrationUtilities for this context.
     */
    OpenCLIntegrationUtilities& getIntegrationUtilities() {
        return *integration;
    }
    /**
     * Get the OpenCLNonbondedUtilities for this context.
     */
    OpenCLNonbondedUtilities& getNonbondedUtilities() {
        return *nonbonded;
    }
    /**
     * Reorder the internal arrays of atoms to try to keep spatially contiguous atoms close
     * together in the arrays.
     */
    void reorderAtoms();
private:
    struct Molecule;
    struct MoleculeGroup;
    void findMoleculeGroups(const System& system);
    static void tagAtomsInMolecule(int atom, int molecule, std::vector<int>& atomMolecule, std::vector<std::vector<int> >& atomBonds);
    double time;
    int deviceIndex;
    int stepCount;
    int computeForceCount;
    int numAtoms;
    int paddedNumAtoms;
    int numAtomBlocks;
    int numThreadBlocks;
    int numForceBuffers;
    int simdWidth;
    mm_float4 periodicBoxSize;
    mm_float4 invPeriodicBoxSize;
    std::string compilationOptions;
    cl::Context context;
    cl::Device device;
    cl::CommandQueue queue;
    cl::Program utilities;
    cl::Kernel clearBufferKernel;
    cl::Kernel clearTwoBuffersKernel;
    cl::Kernel clearThreeBuffersKernel;
    cl::Kernel clearFourBuffersKernel;
    cl::Kernel reduceFloat4Kernel;
    std::vector<OpenCLForceInfo*> forces;
    std::vector<MoleculeGroup> moleculeGroups;
    std::vector<mm_int4> posCellOffsets;
    OpenCLArray<mm_float4>* posq;
    OpenCLArray<mm_float4>* velm;
    OpenCLArray<mm_float4>* force;
    OpenCLArray<mm_float4>* forceBuffers;
    OpenCLArray<cl_float>* energyBuffer;
    OpenCLArray<cl_int>* atomIndex;
    std::vector<cl::Memory*> autoclearBuffers;
    std::vector<int> autoclearBufferSizes;
    OpenCLIntegrationUtilities* integration;
    OpenCLNonbondedUtilities* nonbonded;
};

struct OpenCLContext::MoleculeGroup {
    std::vector<int> atoms;
    std::vector<int> instances;
};

} // namespace OpenMM

#endif /*OPENMM_OPENCLCONTEXT_H_*/
