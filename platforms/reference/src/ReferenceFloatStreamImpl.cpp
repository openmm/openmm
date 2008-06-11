/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "ReferenceFloatStreamImpl.h"

using namespace OpenMM;

ReferenceFloatStreamImpl::ReferenceFloatStreamImpl(std::string name, int size, Stream::DataType type, const Platform& platform) : StreamImpl(name, size, type, platform) {
    switch (type) {
    case Stream::Float:
    case Stream::Float2:
    case Stream::Float3:
    case Stream::Float4:
        baseType = Stream::Float;
        break;
    case Stream::Double:
    case Stream::Double2:
    case Stream::Double3:
    case Stream::Double4:
        baseType = Stream::Double;
        break;
    }
    switch (type) {
    case Stream::Float:
    case Stream::Double:
        width = 1;
        break;
    case Stream::Float2:
    case Stream::Double2:
        width = 2;
        break;
    case Stream::Float3:
    case Stream::Double3:
        width = 3;
        break;
    case Stream::Float4:
    case Stream::Double4:
        width = 4;
        break;
    }
    data = new RealOpenMM*[size];
    for (int i = 0; i < size; ++i)
        data[i] = new RealOpenMM[width];
}

ReferenceFloatStreamImpl::~ReferenceFloatStreamImpl() {
    delete [] data;
}

void ReferenceFloatStreamImpl::loadFromArray(const void* array) {
    if (baseType == Stream::Float) {
        float* arrayData = (float*) array;
        for (int i = 0; i < getSize(); ++i)
            for (int j = 0; j < width; ++j)
                data[i][j] = arrayData[i*width+j];
    }
    else {
        double* arrayData = (double*) array;
        for (int i = 0; i < getSize(); ++i)
            for (int j = 0; j < width; ++j)
                data[i][j] = arrayData[i*width+j];
    }
}

void ReferenceFloatStreamImpl::saveToArray(void* array) {
    if (baseType == Stream::Float) {
        float* arrayData = (float*) array;
        for (int i = 0; i < getSize(); ++i)
            for (int j = 0; j < width; ++j)
                arrayData[i*width+j] = data[i][j];
    }
    else {
        double* arrayData = (double*) array;
        for (int i = 0; i < getSize(); ++i)
            for (int j = 0; j < width; ++j)
                arrayData[i*width+j] = data[i][j];
    }
}

void ReferenceFloatStreamImpl::fillWithValue(void* value) {
    if (baseType == Stream::Float) {
        float valueData = *((float*) value);
        for (int i = 0; i < getSize(); ++i)
            for (int j = 0; j < width; ++j)
                data[i][j] = valueData;
    }
    else {
        double valueData = *((double*) value);
        for (int i = 0; i < getSize(); ++i)
            for (int j = 0; j < width; ++j)
                data[i][j] = valueData;
    }

}

const RealOpenMM* const * ReferenceFloatStreamImpl::getData() const {
    return data;
}

RealOpenMM** ReferenceFloatStreamImpl::getData() {
    return data;
}

