#ifndef OPENMM_STREAMFACTORY_H_
#define OPENMM_STREAMFACTORY_H_

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

#include "StreamImpl.h"
#include "openmm/internal/windowsExport.h"

namespace OpenMM {

/**
 * A StreamFactory is an object that can create StreamImpls.  This is an abstract class.
 * Each Platform maintains a list of StreamFactory objects that it uses to create
 * StreamImpls as needed.
 */

class OPENMM_EXPORT StreamFactory {
public:
    /**
     * Create a StreamImpl.
     * 
     * @param name the name of the stream to create
     * @param size the number of elements in the stream
     * @param type the data type of each element in the stream
     * @param context the context the kernel will belong to
     */
    virtual StreamImpl* createStreamImpl(std::string name, int size, Stream::DataType type, const Platform& platform, ContextImpl& context) const = 0;
    virtual ~StreamFactory() {
    }
};

} // namespace OpenMM

#endif /*OPENMM_STREAMFACTORY_H_*/
