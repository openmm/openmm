#ifndef OPENMM_VECTORIZE_AVX2_H_
#define OPENMM_VECTORIZE_AVX2_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2014 Stanford University and the Authors.      *
 * Authors: Daniel Towner                                                     *
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

#include "vectorizeAvx.h"
#include <immintrin.h>

// This file defines classes and functions to simplify vectorizing code with AVX.

bool isAvx2Supported() {

    // Provide an alternative implementation of CPUID to support AVX2. On older
    // non-Windows OSes the hardware.h support for CPUID doesn't set the CX register
    // properly and gives the wrong answer when detecting AVX2 and beyond. On Windows
    // the cpuid seems to work as expected so can be used.
#if !(defined(_WIN32) || defined(WIN32))
    auto cpuid = [](int output[4], int functionnumber) {
        int a, b, c, d;
        __asm("cpuid" : "=a"(a),"=b"(b),"=c"(c),"=d"(d) : "a"(functionnumber), "c"(0) : );
        output[0] = a;
        output[1] = b;
        output[2] = c;
        output[3] = d;
    };
#endif

    int cpuInfo[4];
    cpuid(cpuInfo, 0);

    if (cpuInfo[0] >= 7) {
        cpuInfo[2] = 0;
        cpuid(cpuInfo, 7);
        return ((cpuInfo[1] & ((int) 1 << 5)) != 0);
    }

    return false;
}

/** 
 * Derive from fvec8 so that default implementations of everything are provided,
 * but can be overriden with AVX2-specific variants where possible.
 */
class fvecAvx2 : public fvec8 {
public:

    fvecAvx2() = default;
    fvecAvx2(fvec8 v) : fvec8(v) {}
    fvecAvx2(float v) : fvec8(v) {}
    fvecAvx2(float v1, float v2, float v3, float v4, float v5, float v6, float v7, float v8) : fvec8(v8, v7, v6, v5, v4, v3, v2, v1) {}
    fvecAvx2(__m256 v) : fvec8(v) {}
    fvecAvx2(const float* v) : fvec8(v) {}

    /** Create a vector by gathering individual indexes of data from a table. Element i of the vector will
     * be loaded from table[idx[i]].
     * @param table The table from which to do a lookup.
     * @param indexes The indexes to gather.
     */
    fvecAvx2(const float* table, const int idx[8])
        : fvec8(_mm256_i32gather_ps(table, _mm256_loadu_si256((const __m256i*)idx), 4)) {}

    static fvecAvx2 expandBitsToMask(int bitmask);
};

inline fvecAvx2 fvecAvx2::expandBitsToMask(int bitmask) {
    // Put a copy of all bits into each vector element and then shift so that the
    // appropriate sub-bit becomes the MSB. For masking purposes, only the MSB matters and
    // the other bits can be completely arbitrary.
    const auto msb = _mm256_sllv_epi32(_mm256_set1_epi8(bitmask),
                                       _mm256_setr_epi32(7, 6, 5, 4, 3, 2, 1, 0));

    return _mm256_castsi256_ps(msb);
}

#endif /*OPENMM_VECTORIZE_AVX2_H_*/
