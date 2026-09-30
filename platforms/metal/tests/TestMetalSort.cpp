/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from platforms/opencl/tests/TestOpenCLSort.cpp.
 * Original OpenCL Platform tests:
 * Portions copyright (c) 2008-2021 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 *
 * Metal Platform tests:
 * Portions copyright (c) 2026 Chun-Chi Hung.
 * Authors: Chun-Chi Hung
 *
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

#include "MetalArray.h"
#include "MetalContext.h"
#include "MetalSort.h"
#include "openmm/internal/AssertionUtilities.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using namespace OpenMM;
using namespace std;

/** @brief Sort ordinary floats with the unchanged OpenCL bucket/bitonic algorithms. */
class FloatSortTrait : public ComputeSortImpl::SortTrait {
public:
    int getDataSize() const override { return sizeof(float); }
    int getKeySize() const override { return sizeof(float); }
    const char* getDataType() const override { return "float"; }
    const char* getKeyType() const override { return "float"; }
    const char* getMinKey() const override { return "-MAXFLOAT"; }
    const char* getMaxKey() const override { return "MAXFLOAT"; }
    const char* getMaxValue() const override { return "MAXFLOAT"; }
    const char* getSortKey() const override { return "value"; }
};

/** @brief Exercise the unsigned-key format used by the nonbonded block sorter. */
class UnsignedSortTrait : public FloatSortTrait {
public:
    const char* getDataType() const override { return "unsigned int"; }
    const char* getKeyType() const override { return "unsigned int"; }
    const char* getMinKey() const override { return "0"; }
    const char* getMaxKey() const override { return "0xFFFFFFFFu"; }
    const char* getMaxValue() const override { return "0xFFFFFFFFu"; }
};

/** @brief Verify exact ordering and value preservation against a host-side oracle. */
void verifySorting(MetalContext& context, vector<float> values, bool uniform) {
    MetalArray data(context, values.size(), sizeof(float), "sortData");
    data.upload(values);
    MetalSort sorter(context, new FloatSortTrait(), values.size(), uniform);
    sorter.sort(data);
    vector<float> result;
    data.download(result);
    sort(values.begin(), values.end());
    ASSERT_EQUAL(values.size(), result.size());
    for (size_t i = 0; i < result.size(); i++)
        ASSERT_EQUAL(values[i], result[i]);
}

/** @brief Check unsigned comparison and preservation of all 32 key bits. */
void verifyUnsignedSorting(MetalContext& context, vector<uint32_t> values, bool uniform) {
    MetalArray data(context, values.size(), sizeof(uint32_t), "unsignedSortData");
    data.upload(values);
    MetalSort sorter(context, new UnsignedSortTrait(), values.size(), uniform);
    sorter.sort(data);
    vector<uint32_t> result;
    data.download(result);
    sort(values.begin(), values.end());
    ASSERT_EQUAL(values.size(), result.size());
    for (size_t i = 0; i < result.size(); i++)
        ASSERT_EQUAL(values[i], result[i]);
}

int main() {
    try {
        System system;
        system.addParticle(1.0);
        MetalContext context(system);
        mt19937 generator(0);
        uniform_real_distribution<float> distribution(0.001f, 1.0f);
        for (int length : {0, 1, 2, 31, 32, 33, 500, 1024, 1025, 3000, 3001, 10000}) {
            vector<float> values(length);
            for (float& value : values)
                value = log(distribution(generator));
            verifySorting(context, values, true);
            verifySorting(context, values, false);
        }
        for (int length : {31, 32, 33, 10000}) {
            vector<float> equal(length, 3.5f);
            verifySorting(context, equal, true);
            verifySorting(context, equal, false);
            vector<float> duplicates(length);
            for (int i = 0; i < length; i++)
                duplicates[i] = (i%17)-8;
            verifySorting(context, duplicates, true);
            verifySorting(context, duplicates, false);
        }
        for (int length : {3, 31, 32, 33, 1025, 3000, 3001, 10000}) {
            vector<uint32_t> values(length);
            for (uint32_t& value : values)
                value = generator();
            values[0] = 0;
            values[1] = 0xFFFFFFFFu;
            values[2] = 0x80000000u;
            verifyUnsignedSorting(context, values, true);
            verifyUnsignedSorting(context, values, false);
        }
    }
    catch (const exception& error) {
        cerr << "exception: " << error.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
