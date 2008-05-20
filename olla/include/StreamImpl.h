#ifndef OPENMM_STREAMIMPL_H_
#define OPENMM_STREAMIMPL_H_

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

#include "Platform.h"
#include "Stream.h"
#include <string>
#include <cassert>

namespace OpenMM {

/**
 * A StreamImpl defines the internal implementation of a Stream object.
 */

class StreamImpl {
public:
    /**
     * Create a StreamImpl.
     * 
     * @param name     the name of the stream to create
     * @param size     the number of elements in the stream
     * @param type     the data type of each element in the stream
     * @param platform the Platform that created this kernel
     */
    StreamImpl(std::string name, int size, Stream::DataType type, const Platform& platform);
    virtual ~StreamImpl() {
        assert(referenceCount == 0);
    }
    /**
     * Get the name of this stream.
     */
    std::string getName() const;
    /**
     * Get the number of elements in this stream.
     */
    int getSize() const;
    /**
     * Get the data type of each element in the stream.
     */
    Stream::DataType getDataType() const;
    /**
     * Get the Platform that created this KernelImpl.
     */
    const Platform& getPlatform();
    /**
     * Copy the contents of an array into this stream.
     * 
     * @param  array a pointer to the start of the array.  The array is assumed to have the same length as this stream,
     * and to contain elements of the correct data type for this stream.  If the stream has a compound data type, all
     * the values should be packed into a single array: all the values for the first element, followed by all the values
     * for the next element, etc.
     */
    virtual void loadFromArray(const void* array) = 0;
    /**
     * Copy the contents of this stream into an array.
     * 
     * @param  array a pointer to the start of the array.  The array is assumed to have the same length as this stream,
     * and to contain elements of the correct data type for this stream.  If the stream has a compound data type, all
     * the values should be packed into a single array: all the values for the first element, followed by all the values
     * for the next element, etc.
     */
    virtual void saveToArray(void* array) = 0;
    /**
     * Set every element of this stream to the same value.
     * 
     * @param a pointer to the value.  It is assumed to be of the correct data type for this stream.
     */
    virtual void fillWithValue(void* value) = 0;
private:
    friend class Stream;
    std::string name;
    int size;
    Stream::DataType type;
    const Platform* platform;
    int referenceCount;
};

} // namespace OpenMM

#endif /*OPENMM_STREAMIMPL_H_*/
