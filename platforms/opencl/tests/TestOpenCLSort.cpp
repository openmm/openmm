/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
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

/**
 * This tests the OpenCL implementation of sorting.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "OpenCLArray.h"
#include "OpenCLContext.h"
#include "OpenCLSort.h"
#include "sfmt/SFMT.h"
#include "openmm/System.h"
#include <iostream>
#include <cmath>
#include <set>

using namespace OpenMM;
using namespace std;

static OpenCLPlatform platform;

class SortTrait : public OpenCLSort::SortTrait {
    int getDataSize() const {return 4;}
    int getKeySize() const {return 4;}
    const char* getDataType() const {return "float";}
    const char* getKeyType() const {return "float";}
    const char* getMinKey() const {return "-MAXFLOAT";}
    const char* getMaxKey() const {return "MAXFLOAT";}
    const char* getMaxValue() const {return "MAXFLOAT";}
    const char* getSortKey() const {return "value";}
};

void verifySorting(vector<float> array) {
    // Sort the array.

    System system;
    system.addParticle(0.0);
    OpenCLPlatform::PlatformData platformData(system, "", "", platform.getPropertyDefaultValue("OpenCLPrecision"), "false", "false", 1, NULL);
    OpenCLContext& context = *platformData.contexts[0];
    context.initialize();
    OpenCLArray data(context, array.size(), sizeof(float), "sortData");
    data.upload(array);
    OpenCLSort sort(context, new SortTrait(), array.size());
    sort.sort(data);
    vector<float> sorted;
    data.download(sorted);

    // Verify that it is in sorted order.

    for (int i = 1; i < (int) sorted.size(); i++)
        ASSERT(sorted[i-1] <= sorted[i]);

    // Make sure the sorted array contains the same values as the original one.

    multiset<float> elements1(array.begin(), array.end());
    multiset<float> elements2(sorted.begin(), sorted.end());
    ASSERT(elements1 == elements2);
}

void testUniformValues() {
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    vector<float> array(10000);
    for (int i = 0; i < (int) array.size(); i++)
        array[i] = (float) genrand_real2(sfmt);
    verifySorting(array);
}

void testLogValues() {
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    vector<float> array(10000);
    for (int i = 0; i < (int) array.size(); i++)
        array[i] = (float) log(genrand_real2(sfmt));
    verifySorting(array);
}

void testShortList() {
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
 
    vector<float> array(500);
    for (int i = 0; i < (int) array.size(); i++)
        array[i] = (float) log(genrand_real2(sfmt));
    verifySorting(array);
}

int main(int argc, char* argv[]) {
    try {
        if (argc > 1)
            platform.setPropertyDefaultValue("OpenCLPrecision", string(argv[1]));
        testUniformValues();
        testLogValues();
        testShortList();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}


