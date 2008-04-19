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

#include "ReferenceIntStreamImpl.h"

using namespace OpenMM;

ReferenceIntStreamImpl::ReferenceIntStreamImpl(std::string name, int size, Stream::DataType type, const Platform& platform) : StreamImpl(name, size, type, platform) {
    switch (type) {
    case Stream::Integer:
        width = 1;
        break;
    case Stream::Integer2:
        width = 2;
        break;
    case Stream::Integer3:
        width = 3;
        break;
    case Stream::Integer4:
        width = 4;
        break;
    }
    data = new int*[size];
    for (int i = 0; i < size; ++i)
        data[i] = new int[width];
}

ReferenceIntStreamImpl::~ReferenceIntStreamImpl() {
    delete data;
}

void ReferenceIntStreamImpl::loadFromArray(const void* array) {
    int* arrayData = (int*) array;
    for (int i = 0; i < getSize(); ++i)
        for (int j = 0; j < width; ++j)
            data[i][j] = arrayData[i*width+j];
}

void ReferenceIntStreamImpl::saveToArray(void* array) {
    int* arrayData = (int*) array;
    for (int i = 0; i < getSize(); ++i)
        for (int j = 0; j < width; ++j)
            arrayData[i*width+j] = data[i][j];
}

void ReferenceIntStreamImpl::fillWithValue(void* value) {
    int valueData = *((int*) value);
    for (int i = 0; i < getSize(); ++i)
        for (int j = 0; j < width; ++j)
            data[i][j] = valueData;
}

int** ReferenceIntStreamImpl::getData() {
    return data;
}

