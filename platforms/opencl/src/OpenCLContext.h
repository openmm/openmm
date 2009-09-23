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

#define __CL_ENABLE_EXCEPTIONS
#include <cl.hpp>

namespace OpenMM {

template <class T>
class OpenCLArray;

/**
 * We can't use cl_float4, since different OpenCL implementations currently define it in
 * incompatible ways.  Hopefully that will be fixed in the future.  In the mean time, we
 * define our own type to represent float4 on the host.
 */

typedef struct {
    cl_float x, y, z, w;
} mm_float4;

/**
 * This class contains the information associated with a Context by the OpenCL Platform.
 */

class OpenCLContext {
public:
    static const int ThreadBlockSize = 64;
    static const int TileSize = 32;
    OpenCLContext(int numParticles, int platformIndex, int deviceIndex);
    ~OpenCLContext();
    /**
     * Get the cl::Context associated with this object.
     */
    cl::Context& getContext() {
        return context;
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
     * Get the array which contains the velocity and massof each atom.
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
     * Get the array which contains the index of each atom.
     */
    OpenCLArray<cl_int>& getAtomIndex() {
        return *atomIndex;
    }
    /**
     * Load OpenCL source code from a file in the kernels directory.
     */
    std::string loadSourceFromFile(const std::string& filename) const;
    /**
     * Create an OpenCL Program from source code.
     */
    cl::Program createProgram(const std::string source);
    /**
     * Set all elements of an array to 0.
     */
    void clearBuffer(OpenCLArray<float>& array);
    /**
     * Set all elements of an array to 0.
     */
    void clearBuffer(OpenCLArray<mm_float4>& array);
    int numAtoms;
    int paddedNumAtoms;
    int numAtomBlocks;
    int numTiles;
    int numThreadBlocks;
    int numForceBuffers;
    bool forceBufferPerWarp;
private:
    cl::Context context;
    cl::Device device;
    cl::CommandQueue queue;
    cl::Program utilities;
    cl::Kernel clearBufferKernel;
    OpenCLArray<mm_float4>* posq;
    OpenCLArray<mm_float4>* velm;
    OpenCLArray<mm_float4>* force;
    OpenCLArray<mm_float4>* forceBuffers;
    OpenCLArray<cl_int>* atomIndex;
};

} // namespace OpenMM

#endif /*OPENMM_OPENCLCONTEXT_H_*/
