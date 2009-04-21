#ifndef OPENMM_STREAM_H_
#define OPENMM_STREAM_H_

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

#include <string>
#include "openmm/internal/windowsExport.h"

namespace OpenMM {

class OPENMM_EXPORT StreamImpl;

/**
 * A Stream encapsulates a particular implementation of a stream of data elements.  A Stream has three
 * fundamental properties:
 * 
 * <ul>
 * <li>A name.  Although a Stream's name has no impact on the API for working with it, a particular platform
 * may choose to implement different Streams in different ways.  The name allows it to distinguish them.</li>
 * <li>A size.  This is the number of data elements contained in the stream.</li>
 * <li>A data type.  There are three basic data types: 32 bit floating point (Float), 64 bit floating point
 * (Double), and 32 bit signed integers (Integer).  In addition, there are compound data types which allow
 * a single stream element to contain multiple values.  For example, the data type Float3 indicates that
 * each stream element is an array of three 32 bit floats.</li>
 * </ul>
 * Stream objects are created by Platforms:
 * 
 * <pre>
 * Stream stream = platform.createStream(streamName, size, dataType);
 * </pre>
 * 
 * A Stream is an opaque reference to the stream data.  You cannot access stream elements directly.  Instead,
 * you must query or modify the stream data by calling methods such as loadFromArray() and saveToArray().
 * These are potentially expensive operations, since they may involve transferring data between main memory
 * and video memory, so they should be called sparing.
 * 
 * It is important to remember that the data type dictates only the API for working with the Stream, not how
 * it is represented internally.  For example, even if a Stream object has type Double, a particular Platform
 * may choose to implement it internally with single precision values.
 */

class OPENMM_EXPORT Stream {
public:
    Stream();
    Stream(const Stream& copy);
    ~Stream();
    Stream& operator=(const Stream& copy);
    /**
     * This is an enumeration of the allowed data types for a Stream.
     */
    enum DataType {Float, Float2, Float3, Float4, Double, Double2, Double3, Double4, Integer, Integer2, Integer3, Integer4};
    /**
     * Get the name of this Stream.
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
     * Copy the contents of an array into this stream.
     * 
     * @param  array a pointer to the start of the array.  The array is assumed to have the same length as this stream,
     * and to contain elements of the correct data type for this stream.  If the stream has a compound data type, all
     * the values should be packed into a single array: all the values for the first element, followed by all the values
     * for the next element, etc.
     */
    void loadFromArray(const void* array);
    /**
     * Copy the contents of this stream into an array.
     * 
     * @param  array a pointer to the start of the array.  The array is assumed to have the same length as this stream,
     * and to contain elements of the correct data type for this stream.  If the stream has a compound data type, all
     * the values should be packed into a single array: all the values for the first element, followed by all the values
     * for the next element, etc.
     */
    void saveToArray(void* array) const;
    /**
     * Set every element of this stream to the same value.
     * 
     * @param a pointer to the value.  It is assumed to be of the correct data type for this stream.
     */
    void fillWithValue(void* value);
    /**
     * Get the object which implements this Stream.
     */
    const StreamImpl& getImpl() const;
    /**
     * Get the object which implements this Stream.
     */
    StreamImpl& getImpl();
private:
    friend class Platform;
    /**
     * Create a StreamImpl.
     * 
     * @param name the name of the stream to create
     */
    Stream(StreamImpl* impl);
    StreamImpl* impl;
};

} // namespace OpenMM

#endif /*OPENMM_STREAM_H_*/
