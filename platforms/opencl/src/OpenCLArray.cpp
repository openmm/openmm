/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2012-2018 Stanford University and the Authors.      *
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
#include "OpenCLContext.h"
#include <iostream>
#include <sstream>
#include <vector>

using namespace OpenMM;

OpenCLArray::OpenCLArray() : buffer(NULL), ownsBuffer(false) {
}

OpenCLArray::OpenCLArray(OpenCLContext& context, int size, int elementSize, const std::string& name, cl_int flags) : buffer(NULL) {
    initialize(context, size, elementSize, name, flags);
}

OpenCLArray::OpenCLArray(OpenCLContext& context, cl::Buffer* buffer, int size, int elementSize, const std::string& name) : buffer(NULL) {
    initialize(context, buffer, size, elementSize, name);
}

OpenCLArray::~OpenCLArray() {
    if (buffer != NULL && ownsBuffer)
        delete buffer;
}

void OpenCLArray::initialize(OpenCLContext& context, int size, int elementSize, const std::string& name, cl_int flags) {
    if (buffer != NULL)
        throw OpenMMException("OpenCLArray has already been initialized");
    this->context = &context;
    this->size = size;
    this->elementSize = elementSize;
    this->name = name;
    this->flags = flags;
    ownsBuffer = true;
    try {
        buffer = new cl::Buffer(context.getContext(), flags, size*elementSize);
    }
    catch (cl::Error err) {
        std::stringstream str;
        str<<"Error creating array "<<name<<": "<<err.what()<<" ("<<err.err()<<")";
        throw OpenMMException(str.str());
    }
}

void OpenCLArray::initialize(OpenCLContext& context, cl::Buffer* buffer, int size, int elementSize, const std::string& name) {
    if (this->buffer != NULL)
        throw OpenMMException("OpenCLArray has already been initialized");
    this->context = &context;
    this->buffer = buffer;
    this->size = size;
    this->elementSize = elementSize;
    this->name = name;
    ownsBuffer = false;
}

void OpenCLArray::resize(int size) {
    if (buffer == NULL)
        throw OpenMMException("OpenCLArray has not been initialized");
    if (!ownsBuffer)
        throw OpenMMException("Cannot resize an array that does not own its storage");
    delete buffer;
    buffer = NULL;
    initialize(*context, size, elementSize, name, flags);
}

void OpenCLArray::upload(const void* data, bool blocking) {
    if (buffer == NULL)
        throw OpenMMException("OpenCLArray has not been initialized");
    try {
        context->getQueue().enqueueWriteBuffer(*buffer, blocking ? CL_TRUE : CL_FALSE, 0, size*elementSize, data);
    }
    catch (cl::Error err) {
        std::stringstream str;
        str<<"Error uploading array "<<name<<": "<<err.what()<<" ("<<err.err()<<")";
        throw OpenMMException(str.str());
    }
}

void OpenCLArray::download(void* data, bool blocking) const {
    if (buffer == NULL)
        throw OpenMMException("OpenCLArray has not been initialized");
    try {
        context->getQueue().enqueueReadBuffer(*buffer, blocking ? CL_TRUE : CL_FALSE, 0, size*elementSize, data);
    }
    catch (cl::Error err) {
        std::stringstream str;
        str<<"Error downloading array "<<name<<": "<<err.what()<<" ("<<err.err()<<")";
        throw OpenMMException(str.str());
    }
}

void OpenCLArray::copyTo(OpenCLArray& dest) const {
    if (buffer == NULL)
        throw OpenMMException("OpenCLArray has not been initialized");
    if (dest.getSize() != size || dest.getElementSize() != elementSize)
        throw OpenMMException("Error copying array "+name+" to "+dest.getName()+": The destination array does not match the size of the array");
    try {
        context->getQueue().enqueueCopyBuffer(*buffer, dest.getDeviceBuffer(), 0, 0, size*elementSize);
    }
    catch (cl::Error err) {
        std::stringstream str;
        str<<"Error copying array "<<name<<" to "<<dest.getName()<<": "<<err.what()<<" ("<<err.err()<<")";
        throw OpenMMException(str.str());
    }
}
