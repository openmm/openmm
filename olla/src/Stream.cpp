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

#include "Stream.h"
#include "StreamImpl.h"

using namespace OpenMM;
using namespace std;

Stream::Stream() : impl(0) {
}

Stream::Stream(StreamImpl* impl) : impl(impl) {
}

Stream::Stream(const Stream& copy) : impl(copy.impl) {
    impl->referenceCount++;
}

Stream::~Stream() {
    if (impl) {
        impl->referenceCount--;
        if (impl->referenceCount == 0)
            delete impl;
    }
}

Stream& Stream::operator=(const Stream& copy) {
    impl = copy.impl;
    impl->referenceCount++;
    return *this;
}

string Stream::Stream::getName() const {
    return impl->getName();
}

int Stream::getSize() const {
    return impl->getSize();
}

Stream::DataType Stream::getDataType() const {
    return impl->getDataType();
}

void Stream::loadFromArray(const void* array) {
    impl->loadFromArray(array);
}

void Stream::saveToArray(void* array) {
    impl->saveToArray(array);
}

void Stream::fillWithValue(void* value) {
    impl->fillWithValue(value);
}

StreamImpl& Stream::getImpl() {
    return *impl;
}
