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

#include "CudaArray.h"
#include "CudaContext.h"
#include <iostream>
#include <sstream>
#include <vector>

using namespace OpenMM;

CudaArray::CudaArray() : pointer(0), ownsMemory(false) {
}

CudaArray::CudaArray(CudaContext& context, int size, int elementSize, const std::string& name) : pointer(0) {
    initialize(context, size, elementSize, name);
}

CudaArray::~CudaArray() {
    if (pointer != 0 && ownsMemory && context->getContextIsValid()) {
        context->setAsCurrent();
        CUresult result = cuMemFree(pointer);
        if (result != CUDA_SUCCESS) {
            std::stringstream str;
            str<<"Error deleting array "<<name<<": "<<CudaContext::getErrorString(result)<<" ("<<result<<")";
            throw OpenMMException(str.str());
        }
    }
}

void CudaArray::initialize(CudaContext& context, int size, int elementSize, const std::string& name) {
    if (this->pointer != 0)
        throw OpenMMException("CudaArray has already been initialized");
    this->context = &context;
    this->size = size;
    this->elementSize = elementSize;
    this->name = name;
    ownsMemory = true;
    CUresult result = cuMemAlloc(&pointer, size*elementSize);
    if (result != CUDA_SUCCESS) {
        std::stringstream str;
        str<<"Error creating array "<<name<<": "<<CudaContext::getErrorString(result)<<" ("<<result<<")";
        throw OpenMMException(str.str());
    }
}

void CudaArray::resize(int size) {
    if (pointer == 0)
        throw OpenMMException("CudaArray has not been initialized");
    if (!ownsMemory)
        throw OpenMMException("Cannot resize an array that does not own its storage");
    CUresult result = cuMemFree(pointer);
    if (result != CUDA_SUCCESS) {
        std::stringstream str;
        str<<"Error deleting array "<<name<<": "<<CudaContext::getErrorString(result)<<" ("<<result<<")";
        throw OpenMMException(str.str());
    }
    pointer = 0;
    initialize(*context, size, elementSize, name);
}

void CudaArray::upload(const void* data, bool blocking) {
    if (pointer == 0)
        throw OpenMMException("CudaArray has not been initialized");
    CUresult result;
    if (blocking)
        result = cuMemcpyHtoD(pointer, data, size*elementSize);
    else
        result = cuMemcpyHtoDAsync(pointer, data, size*elementSize, context->getCurrentStream());
    if (result != CUDA_SUCCESS) {
        std::stringstream str;
        str<<"Error uploading array "<<name<<": "<<CudaContext::getErrorString(result)<<" ("<<result<<")";
        throw OpenMMException(str.str());
    }
}

void CudaArray::download(void* data, bool blocking) const {
    if (pointer == 0)
        throw OpenMMException("CudaArray has not been initialized");
    CUresult result;
    if (blocking)
        result = cuMemcpyDtoH(data, pointer, size*elementSize);
    else
        result = cuMemcpyDtoHAsync(data, pointer, size*elementSize, context->getCurrentStream());
    if (result != CUDA_SUCCESS) {
        std::stringstream str;
        str<<"Error downloading array "<<name<<": "<<CudaContext::getErrorString(result)<<" ("<<result<<")";
        throw OpenMMException(str.str());
    }
}

void CudaArray::copyTo(CudaArray& dest) const {
    if (pointer == 0)
        throw OpenMMException("CudaArray has not been initialized");
    if (dest.getSize() != size || dest.getElementSize() != elementSize)
        throw OpenMMException("Error copying array "+name+" to "+dest.getName()+": The destination array does not match the size of the array");
    CUresult result = cuMemcpyDtoDAsync(dest.getDevicePointer(), pointer, size*elementSize, context->getCurrentStream());
    if (result != CUDA_SUCCESS) {
        std::stringstream str;
        str<<"Error copying array "<<name<<" to "<<dest.getName()<<": "<<CudaContext::getErrorString(result)<<" ("<<result<<")";
        throw OpenMMException(str.str());
    }
}
