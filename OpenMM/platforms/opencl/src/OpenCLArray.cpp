/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
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

#include "OpenCLArray.h"
#include <iostream>
#include <sstream>
#include <vector>

using namespace OpenMM;

OpenCLArray::OpenCLArray(OpenCLContext& context, int size, int elementSize, const std::string& name, cl_int flags) :
        context(context), size(size), elementSize(elementSize), name(name), ownsBuffer(true) {
    try {
        buffer = new cl::Buffer(context.getContext(), flags, size*elementSize);
    }
    catch (cl::Error err) {
        std::stringstream str;
        str<<"Error creating array "<<name<<": "<<err.what()<<" ("<<err.err()<<")";
        throw OpenMMException(str.str());
    }
}

OpenCLArray::OpenCLArray(OpenCLContext& context, cl::Buffer* buffer, int size, int elementSize, const std::string& name) :
        context(context), buffer(buffer), size(size), elementSize(elementSize), name(name), ownsBuffer(false) {
}

OpenCLArray::~OpenCLArray() {
    if (ownsBuffer)
        delete buffer;
}

void OpenCLArray::upload(const void* data, bool blocking) {
    try {
        context.getQueue().enqueueWriteBuffer(*buffer, blocking ? CL_TRUE : CL_FALSE, 0, size*elementSize, data);
    }
    catch (cl::Error err) {
        std::stringstream str;
        str<<"Error uploading array "<<name<<": "<<err.what()<<" ("<<err.err()<<")";
        throw OpenMMException(str.str());
    }
}

void OpenCLArray::download(void* data, bool blocking) const {
    try {
        context.getQueue().enqueueReadBuffer(*buffer, blocking ? CL_TRUE : CL_FALSE, 0, size*elementSize, data);
    }
    catch (cl::Error err) {
        std::stringstream str;
        str<<"Error downloading array "<<name<<": "<<err.what()<<" ("<<err.err()<<")";
        throw OpenMMException(str.str());
    }
}

void OpenCLArray::copyTo(OpenCLArray& dest) const {
    if (dest.getSize() != size || dest.getElementSize() != elementSize)
        throw OpenMMException("Error copying array "+name+" to "+dest.getName()+": The destination array does not match the size of the array");
    try {
        context.getQueue().enqueueCopyBuffer(*buffer, dest.getDeviceBuffer(), 0, 0, size*elementSize);
    }
    catch (cl::Error err) {
        std::stringstream str;
        str<<"Error copying array "<<name<<" to "<<dest.getName()<<": "<<err.what()<<" ("<<err.err()<<")";
        throw OpenMMException(str.str());
    }
}
