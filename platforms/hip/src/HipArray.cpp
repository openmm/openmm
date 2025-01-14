/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2012-2022 Stanford University and the Authors.      *
 * Portions copyright (c) 2020-2023 Advanced Micro Devices, Inc.              *
 * Authors: Peter Eastman, Nicholas Curtis                                    *
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

#include "HipArray.h"
#include "HipContext.h"
#include "openmm/common/ContextSelector.h"
#include <iostream>
#include <sstream>
#include <vector>

using namespace OpenMM;

HipArray::HipArray() : pointer(0), ownsMemory(false) {
}

HipArray::HipArray(HipContext& context, size_t size, int elementSize, const std::string& name) : pointer(0) {
    initialize(context, size, elementSize, name);
}

HipArray::~HipArray() {
    if (pointer != 0 && ownsMemory && context->getContextIsValid()) {
        ContextSelector selector(*context);
        hipError_t result = hipFree(pointer);
        if (result != hipSuccess) {
            std::stringstream str;
            str<<"Error deleting array "<<name<<": "<<HipContext::getErrorString(result)<<" ("<<result<<")";
            throw OpenMMException(str.str());
        }
    }
}

void HipArray::initialize(ComputeContext& context, size_t size, int elementSize, const std::string& name) {
    if (this->pointer != 0)
        throw OpenMMException("HipArray has already been initialized");
    this->context = &dynamic_cast<HipContext&>(context);
    this->size = size;
    this->elementSize = elementSize;
    this->name = name;
    ownsMemory = true;
    ContextSelector selector(*this->context);
    hipError_t result = hipMalloc(&pointer, size*elementSize);
    if (result != hipSuccess) {
        std::stringstream str;
        str<<"Error creating array "<<name<<": "<<HipContext::getErrorString(result)<<" ("<<result<<")";
        throw OpenMMException(str.str());
    }
}

void HipArray::resize(size_t size) {
    if (pointer == 0)
        throw OpenMMException("HipArray has not been initialized");
    if (!ownsMemory)
        throw OpenMMException("Cannot resize an array that does not own its storage");
    ContextSelector selector(*context);
    hipError_t result = hipFree(pointer);
    if (result != hipSuccess) {
        std::stringstream str;
        str<<"Error deleting array "<<name<<": "<<HipContext::getErrorString(result)<<" ("<<result<<")";
        throw OpenMMException(str.str());
    }
    pointer = 0;
    initialize(*context, size, elementSize, name);
}

ComputeContext& HipArray::getContext() {
    return *context;
}

void HipArray::uploadSubArray(const void* data, int offset, int elements, bool blocking) {
    if (pointer == 0)
        throw OpenMMException("HipArray has not been initialized");
    if (offset < 0 || offset+elements > getSize())
        throw OpenMMException("uploadSubArray: data exceeds range of array");
    hipError_t result;
    result = hipMemcpyAsync(reinterpret_cast<char*>(pointer)+offset*elementSize, const_cast<void*>(data), elements*elementSize, hipMemcpyHostToDevice, context->getCurrentStream());
    if (blocking && result == hipSuccess)
        result = hipStreamSynchronize(context->getCurrentStream());
    if (result != hipSuccess) {
        std::stringstream str;
        str<<"Error uploading array "<<name<<": "<<HipContext::getErrorString(result)<<" ("<<result<<")";
        throw OpenMMException(str.str());
    }
}

void HipArray::download(void* data, bool blocking) const {
    if (pointer == 0)
        throw OpenMMException("HipArray has not been initialized");
    hipError_t result;
    result = hipMemcpyAsync(data, pointer, size*elementSize, hipMemcpyDeviceToHost, context->getCurrentStream());
    if (blocking && result == hipSuccess)
        result = hipStreamSynchronize(context->getCurrentStream());
    if (result != hipSuccess) {
        std::stringstream str;
        str<<"Error downloading array "<<name<<": "<<HipContext::getErrorString(result)<<" ("<<result<<")";
        throw OpenMMException(str.str());
    }
}

void HipArray::copyTo(ArrayInterface& dest) const {
    if (pointer == 0)
        throw OpenMMException("HipArray has not been initialized");
    if (dest.getSize() != size || dest.getElementSize() != elementSize)
        throw OpenMMException("Error copying array "+name+" to "+dest.getName()+": The destination array does not match the size of the array");
    HipArray& cuDest = context->unwrap(dest);
    hipError_t result = hipMemcpyAsync(cuDest.getDevicePointer(), pointer, size*elementSize, hipMemcpyDeviceToDevice, context->getCurrentStream());
    if (result != hipSuccess) {
        std::stringstream str;
        str<<"Error copying array "<<name<<" to "<<dest.getName()<<": "<<HipContext::getErrorString(result)<<" ("<<result<<")";
        throw OpenMMException(str.str());
    }
}
