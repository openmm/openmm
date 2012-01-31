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
 * This tests the Cuda implementation of random number generation.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "../src/kernels/gputypes.h"
#include "../src/kernels/cudaKernels.h"
#include <iostream>

using namespace OpenMM;
using namespace std;

static const float KILO                     =    1e3;                      // Thousand
static const float BOLTZMANN                =    1.380658e-23f;            // (J/K)    
static const float AVOGADRO                 =    6.0221367e23f;            // ()        
static const float RGAS                     =    BOLTZMANN * AVOGADRO;     // (J/(mol K))
static const float BOLTZ                    =    (RGAS / KILO);            // (kJ/(mol K)) 

void testGaussian() {
    _gpuContext* gpu = (_gpuContext*) gpuInit(5000, 0);
    gpuSetConstants(gpu);
    kGenerateRandoms(gpu);
    const int numValues = 4*gpu->psRandom4->_length;
    gpu->psRandom4->Download();
    float* data = reinterpret_cast<float*>(gpu->psRandom4->_pSysData);
    double mean = 0.0;
    double var = 0.0;
    double skew = 0.0;
    double kurtosis = 0.0;
    for (int i = 0; i < numValues; i++) {
        double value = data[i];
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
    ASSERT_EQUAL_TOL(0.0, mean, 3.0/sqrt((double)numValues));
    ASSERT_EQUAL_TOL(1.0, c2, 3.0/pow(numValues, 1.0/3.0));
    ASSERT_EQUAL_TOL(0.0, c3, 3.0/pow(numValues, 1.0/4.0));
    ASSERT_EQUAL_TOL(0.0, c4, 3.0/pow(numValues, 1.0/4.0));
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
