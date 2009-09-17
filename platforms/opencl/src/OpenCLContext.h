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
#include <CL/cl.hpp>

namespace OpenMM {

template <class T>
class OpenCLArray;

/**
 * This class contains the information associated with a Context by the OpenCL Platform.
 */

class OpenCLContext {
public:
    OpenCLContext(int numParticles, int platformIndex, int deviceIndex);
    ~OpenCLContext();
    /**
     * Get the cl::Context associated with this object.
     */
    cl::Context& getContext() {
        return *context;
    }
    /**
     * Get the cl::CommandQueue associated with this object.
     */
    cl::CommandQueue& getQueue() {
        return *queue;
    }
    /**
     * Get the array which contains the position and charge of each atom.
     */
    OpenCLArray<cl_float4>& getPosq() {
        return *posq;
    }
    /**
     * Get the array which contains the velocity and massof each atom.
     */
    OpenCLArray<cl_float4>& getVelm() {
        return *velm;
    }
    /**
     * Get the array which contains the force on each atom.
     */
    OpenCLArray<cl_float4>& getForce() {
        return *force;
    }
    /**
     * Get the array which contains the index of each atom.
     */
    OpenCLArray<cl_int>& getAtomIndex() {
        return *atomIndex;
    }
private:
    cl::Context* context;
    cl::CommandQueue* queue;
    OpenCLArray<cl_float4>* posq;
    OpenCLArray<cl_float4>* velm;
    OpenCLArray<cl_float4>* force;
    OpenCLArray<cl_int>* atomIndex;
};

} // namespace OpenMM

#endif /*OPENMM_OPENCLCONTEXT_H_*/
