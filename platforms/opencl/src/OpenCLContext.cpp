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

#include "OpenCLContext.h"
#include "OpenCLArray.h"

using namespace OpenMM;

OpenCLContext::OpenCLContext(int numParticles, int platformIndex, int deviceIndex) {
    // TODO Select the platform and device correctly
    context = new cl::Context(CL_DEVICE_TYPE_CPU);
    queue = new cl::CommandQueue(getContext(), getContext().getInfo<CL_CONTEXT_DEVICES>()[0]);
    posq = new OpenCLArray<cl_float4>(*this, numParticles, "posq", true);
    velm = new OpenCLArray<cl_float4>(*this, numParticles, "velm", true);
    force = new OpenCLArray<cl_float4>(*this, numParticles, "force", true);
    atomIndex = new OpenCLArray<cl_int>(*this, numParticles, "atomIndex", true);
}

OpenCLContext::~OpenCLContext() {
    delete context;
    delete queue;
    delete posq;
    delete velm;
    delete force;
    delete atomIndex;
}

