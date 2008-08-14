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

/**
 * This tests the reference implementation of random number generation.
 */

#include "../../../tests/AssertionUtilities.h"
#include "../src/SimTKUtilities/SimTKOpenMMUtilities.h"
#include <iostream>

using namespace OpenMM;
using namespace std;

void testGaussian() {
    mt_init(0);
    const int numValues = 10000000;
    double mean = 0.0;
    double var = 0.0;
    double skew = 0.0;
    double kurtosis = 0.0;
    unsigned long jran = 12399103;
    for (int i = 0; i < numValues; i++) {
        double value = SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
        mean += value;
        var += value*value;
        skew += value*value*value;
        kurtosis += value*value*value*value;
    }
    mean /= numValues;
    var /= numValues;
    skew /= numValues;
    kurtosis /= numValues;
    double c2 = var-mean*mean;
    double c3 = skew-3*var*mean+2*mean*mean*mean;
    double c4 = kurtosis-4*skew*mean-3*var*var+12*var*mean*mean-6*mean*mean*mean*mean;
    ASSERT_EQUAL_TOL(0.0, mean, 0.01);
    ASSERT_EQUAL_TOL(1.0, c2, 0.01);
    ASSERT_EQUAL_TOL(0.0, c3, 0.01);
    ASSERT_EQUAL_TOL(0.0, c4, 0.01);
}

int main() {
    try {
        testGaussian();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
