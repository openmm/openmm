/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Mark Friedrichs                                    *
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

#include "OpenMMException.h"
#include "BrookStreamFactory.h"
#include "BrookFloatStreamImpl.h"
#include "BrookIntStreamImpl.h"

using namespace OpenMM;

StreamImpl* BrookStreamFactory::createStreamImpl(std::string name, int size, Stream::DataType type, int streamWidth, const Platform& platform, OpenMMContextImpl& context) const {
    switch (type) {
    case Stream::Float:
    case Stream::Float2:
    case Stream::Float3:
    case Stream::Float4:
    case Stream::Double:
    case Stream::Double2:
    case Stream::Double3:
    case Stream::Double4:
        return new BrookFloatStreamImpl(name, size, type, streamWidth, platform);
    case Stream::Integer:
    case Stream::Integer2:
    case Stream::Integer3:
    case Stream::Integer4:
        return new BrookIntStreamImpl(name, size, type, streamWidth, platform);
    }
    throw OpenMMException("Tried to create a Stream with an illegal DataType.");
}
