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
#include "OpenCLForceInfo.h"
#include "OpenCLIntegrationUtilities.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include <fstream>
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

OpenCLContext::OpenCLContext(int numParticles, int deviceIndex) : time(0.0), stepCount(0) {
    context = cl::Context(CL_DEVICE_TYPE_ALL);
    vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();
    const int minThreadBlockSize = 32;
    if (deviceIndex < 0 || deviceIndex >= devices.size()) {
        // Try to figure out which device is the fastest.

        int bestSpeed = 0;
        for (int i = 0; i < devices.size(); i++) {
            int maxSize = devices[i].getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>()[0];
            int speed = devices[i].getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>()*devices[i].getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>();
            if (maxSize >= minThreadBlockSize && speed > bestSpeed)
                deviceIndex = i;
        }
    }
    if (deviceIndex == -1)
        throw OpenMMException("No compatible OpenCL device is available");
    device = devices[deviceIndex];
    if (device.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>()[0] < minThreadBlockSize)
        throw OpenMMException("The specified OpenCL device is not compatible with OpenMM");
    queue = cl::CommandQueue(context, device);
    numAtoms = numParticles;
    paddedNumAtoms = TileSize*((numParticles+TileSize-1)/TileSize);
    numAtomBlocks = (paddedNumAtoms+(TileSize-1))/TileSize;
    numTiles = numAtomBlocks*(numAtomBlocks+1)/2;
    numThreadBlocks = device.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>()[0]/ThreadBlockSize;

    // Create utility kernels that are used in multiple places.

    utilities = createProgram(loadSourceFromFile("utilities.cl"));
    clearBufferKernel = cl::Kernel(utilities, "clearBuffer");
    reduceFloat4Kernel = cl::Kernel(utilities, "reduceFloat4Buffer");
}

OpenCLContext::~OpenCLContext() {
    for (int i = 0; i < (int) forces.size(); i++)
        delete forces[i];
    delete posq;
    delete velm;
    delete force;
    delete forceBuffers;
    delete energyBuffer;
    delete atomIndex;
    delete integration;
}

void OpenCLContext::initialize(const System& system) {
//    forceBufferPerWarp = true;
//    numForceBuffers = numThreadBlocks*ThreadBlockSize/TileSize;
//    if (numForceBuffers >= numAtomBlocks) {
//        // For small systems, it is more efficient to have one force buffer per block of 32 atoms instead of one per warp.
//
//        forceBufferPerWarp = false;
//        numForceBuffers = numAtomBlocks;
//    }
    posq = new OpenCLArray<mm_float4>(*this, paddedNumAtoms, "posq", true);
    velm = new OpenCLArray<mm_float4>(*this, paddedNumAtoms, "velm", true);
    for (int i = 0; i < numAtoms; i++)
        (*velm)[i].w = (float) (1.0/system.getParticleMass(i));
    velm->upload();
    numForceBuffers = 1;
    for (int i = 0; i < (int) forces.size(); i++)
        numForceBuffers = std::max(numForceBuffers, forces[i]->getRequiredForceBuffers());
    forceBuffers = new OpenCLArray<mm_float4>(*this, paddedNumAtoms*numForceBuffers, "forceBuffers", false);
    force = new OpenCLArray<mm_float4>(*this, &forceBuffers->getDeviceBuffer(), paddedNumAtoms, "force", true);
    energyBuffer = new OpenCLArray<cl_float>(*this, numThreadBlocks*ThreadBlockSize, "energyBuffer", true);
    atomIndex = new OpenCLArray<cl_int>(*this, paddedNumAtoms, "atomIndex", true);
    for (int i = 0; i < paddedNumAtoms; ++i)
        (*atomIndex)[i] = i;
    atomIndex->upload();
    integration = new OpenCLIntegrationUtilities(*this, system);
}

void OpenCLContext::addForce(OpenCLForceInfo* force) {
    forces.push_back(force);
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

void OpenCLContext::executeKernel(cl::Kernel& kernel, int workUnits) {
    int size = std::min((workUnits+ThreadBlockSize-1)/ThreadBlockSize, numThreadBlocks)*ThreadBlockSize;
    try {
        queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(size), cl::NDRange(ThreadBlockSize));
    }
    catch (cl::Error err) {
        stringstream str;
        str<<"Error invoking kernel "<<kernel.getInfo<CL_KERNEL_FUNCTION_NAME>()<<": "<<err.what()<<" ("<<err.err()<<")";
        throw OpenMMException(str.str());
    }
}

void OpenCLContext::clearBuffer(OpenCLArray<float>& array) {
    clearBufferKernel.setArg<cl::Buffer>(0, array.getDeviceBuffer());
    clearBufferKernel.setArg<cl_int>(1, array.getSize());
    executeKernel(clearBufferKernel, array.getSize());
}

void OpenCLContext::clearBuffer(OpenCLArray<mm_float4>& array) {
    clearBufferKernel.setArg<cl::Buffer>(0, array.getDeviceBuffer());
    clearBufferKernel.setArg<cl_int>(1, array.getSize()*4);
    executeKernel(clearBufferKernel, array.getSize()*4);
}

void OpenCLContext::reduceBuffer(OpenCLArray<mm_float4>& array, int numBuffers) {
    int bufferSize = array.getSize()/numBuffers;
    reduceFloat4Kernel.setArg<cl::Buffer>(0, array.getDeviceBuffer());
    reduceFloat4Kernel.setArg<cl_int>(1, bufferSize);
    reduceFloat4Kernel.setArg<cl_int>(2, numBuffers);
    executeKernel(reduceFloat4Kernel, bufferSize);
}
