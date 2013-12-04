#ifndef OPENMM_ALIGNEDARRAY_H_
#define OPENMM_ALIGNEDARRAY_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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

namespace OpenMM {

/**
 * This class represents an array in memory whose starting point is guaranteed to
 * be aligned with a 16 byte boundary.  This can improve the performance of vectorized
 * code, since loads and stores are more efficient.
 */
template <class T>
class AlignedArray {
public:
    /**
     * Default constructor, to allow AlignedArrays to be used inside collections.
     */
    AlignedArray() : dataSize(0), baseData(0), data(0) {
    }
    /**
     * Create an Aligned array that contains a specified number of elements.
     */
    AlignedArray(int size) {
        allocate(size);
    }
    ~AlignedArray() {
        if (baseData != 0)
            delete[] baseData;
    }
    /**
     * Get the number of elements in the array.
     */
    int size() const {
        return dataSize;
    }
    /**
     * Change the size of the array.  This may cause all contents to be lost.
     */
    void resize(int size) {
        if (dataSize == size)
            return;
        if (baseData != 0)
            delete[] baseData;
        allocate(size);
    }
    /**
     * Get a reference to an element of the array.
     */
    T& operator[](int i) {
        return data[i];
    }
    /**
     * Get a const reference to an element of the array.
     */
    const T& operator[](int i) const {
        return data[i];
    }
private:
    void allocate(int size) {
        dataSize = size;
        baseData = new char[size*sizeof(T)+16];
        char* offsetData = baseData+15;
        offsetData -= (long long)offsetData&0xF;
        data = (T*) offsetData;
    }
    int dataSize;
    char* baseData;
    T* data;
};

} // namespace OpenMM

#endif /*OPENMM_ALIGNEDARRAY_H_*/

