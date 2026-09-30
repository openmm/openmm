/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Metal Platform code: Portions copyright (c) 2026 Chun-Chi Hung.             *
 * Authors: Chun-Chi Hung                                                     *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING    *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OF THE SOFTWARE.*
 * -------------------------------------------------------------------------- */

#include "MetalContext.h"
#include "MetalKernel.h"
#include "MetalReductionOptimizations.h"
#include "CommonKernelSources.h"
#include "MetalOpenCLKernelSources.h"
#include "openmm/common/ComputeArray.h"
#include "openmm/internal/AssertionUtilities.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <numeric>

using namespace OpenMM;
using namespace std;

/** @brief Verify independent source selection and immunity of unrelated same-name helpers. */
void testSelection() {
    using Settings = MetalReductionOptimizations::Settings;
    Settings none;
    for (const string* source : {&CommonKernelSources::rg, &CommonKernelSources::rmsd,
            &CommonKernelSources::orientationRestraintForce, &CommonKernelSources::customCentroidBond,
            &CommonKernelSources::lcpo, &CommonKernelSources::customManyParticle, &MetalOpenCLKernelSources::sort})
        ASSERT_EQUAL(*source, MetalReductionOptimizations::apply(*source, none));
    Settings all;
    const string unrelated = "DEVICE real reduceValue(real value, LOCAL_ARG volatile real* temp) { return value; }";
    ASSERT_EQUAL(unrelated, MetalReductionOptimizations::apply(unrelated, all));
}

int main(int argc, char* argv[]) {
    try {
        testSelection();
        if (argc == 2 && string(argv[1]) == "--source-only") {
            cout << "Metal reduction/scan source-selection tests passed" << endl;
            return 0;
        }
    }
    catch (const exception& error) {
        cerr << "exception: " << error.what() << endl;
        return 1;
    }
    cout << "Metal reduction/scan optimization tests passed" << endl;
    return 0;
}
