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
#include "openmm/Platform.h"
#include <fstream>
#include <iostream>

using namespace OpenMM;
using namespace std;

OpenCLContext::OpenCLContext(int numParticles, int platformIndex, int deviceIndex) {
    // TODO Select the platform and device correctly
    context = cl::Context(CL_DEVICE_TYPE_CPU);
    device = context.getInfo<CL_CONTEXT_DEVICES>()[0];
    queue = cl::CommandQueue(context, device);
    numAtoms = numParticles;
    paddedNumAtoms = TileSize*((numParticles+TileSize-1)/TileSize);
    numAtomBlocks = (paddedNumAtoms+(TileSize-1))/TileSize;
    numTiles = numAtomBlocks*(numAtomBlocks+1)/2;
    numThreadBlocks = 8*device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
    forceBufferPerWarp = true;
    numForceBuffers = numThreadBlocks*ThreadBlockSize/TileSize;
    if (numForceBuffers >= numAtomBlocks) {
        // For small systems, it is more efficient to have one force buffer per block of 32 atoms instead of one per warp.

        forceBufferPerWarp = false;
        numForceBuffers = numAtomBlocks;
    }
    posq = new OpenCLArray<mm_float4>(*this, paddedNumAtoms, "posq", true);
    velm = new OpenCLArray<mm_float4>(*this, paddedNumAtoms, "velm", true);
    forceBuffers = new OpenCLArray<mm_float4>(*this, paddedNumAtoms*numForceBuffers, "forceBuffers", false);
    force = new OpenCLArray<mm_float4>(*this, &forceBuffers->getDeviceBuffer(), paddedNumAtoms, "force", true);
    atomIndex = new OpenCLArray<cl_int>(*this, paddedNumAtoms, "atomIndex", true);
    for (int i = 0; i < paddedNumAtoms; ++i)
        atomIndex->set(i, i);
    atomIndex->upload();

    // Create utility kernels that are used in multiple places.

    utilities = createProgram(loadSourceFromFile("utilities.cl"));
    clearBufferKernel = cl::Kernel(utilities, "clearBuffer");
}

OpenCLContext::~OpenCLContext() {
    delete posq;
    delete velm;
    delete force;
    delete atomIndex;
}

string OpenCLContext::loadSourceFromFile(const string& filename) const {
    ifstream file((Platform::getDefaultPluginsDirectory()+"/opencl/"+filename).c_str());
    if (!file.is_open())
        throw OpenMMException("Unable to load kernel: "+filename);
    string kernel;
    string line;
    while (!file.eof()) {
        getline(file, line);
        kernel += line;
        kernel += '\n';
    }
    file.close();
    return kernel;
}

cl::Program OpenCLContext::createProgram(const std::string source) {
    cl::Program::Sources sources(1, make_pair(source.c_str(), source.size()));
    cl::Program program(context, sources);
    try {
        program.build(vector<cl::Device>(1, device));
    } catch (cl::Error err) {
        throw OpenMMException("Error compiling kernel: "+program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device));
    }
    return program;
}

void OpenCLContext::clearBuffer(OpenCLArray<float>& array) {
    clearBufferKernel.setArg<cl::Buffer>(0, array.getDeviceBuffer());
    clearBufferKernel.setArg<cl_int>(1, array.getSize());
    queue.enqueueNDRangeKernel(clearBufferKernel, cl::NullRange, cl::NDRange(numThreadBlocks*ThreadBlockSize), cl::NDRange(ThreadBlockSize));
}

void OpenCLContext::clearBuffer(OpenCLArray<mm_float4>& array) {
    clearBufferKernel.setArg<cl::Buffer>(0, array.getDeviceBuffer());
    clearBufferKernel.setArg<cl_int>(1, array.getSize()*4);
    queue.enqueueNDRangeKernel(clearBufferKernel, cl::NullRange, cl::NDRange(numThreadBlocks*ThreadBlockSize), cl::NDRange(ThreadBlockSize));
}