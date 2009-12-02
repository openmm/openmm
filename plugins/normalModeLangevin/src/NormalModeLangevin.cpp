/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Chris Sweet                                                       *
 * Contributors: Christopher Bruns                                            *
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

#include <iostream>
#include "OpenMM.h"
#include "ReferenceNMLKernelFactory.h"

using namespace OpenMM;
using namespace std;
 
extern "C" void OPENMM_EXPORT registerPlatforms() {
    cout << "calling NML registerPlatforms()..." << endl;
}

extern "C" void OPENMM_EXPORT registerKernelFactories() {
    cout << "Initializing Normal Mode Langevin OpenMM plugin..." << endl;
    for (int p = 0; p < Platform::getNumPlatforms(); ++p) {
        cout << "Platform " << p << " name = " << Platform::getPlatform(p).getName() << endl;
    }

    // Only register cuda kernels if cuda platform is found
    cout << "NML looking for Cuda plugin..." << endl;
    try {
        Platform& platform = Platform::getPlatformByName("Cuda");
        cout << "NML found Cuda platform..." << endl;
        // platform.registerKernelFactory("CudaNMLKernelFactory", new CudaNMLKernelFactory());
    } catch (const std::exception& exc) { // non fatal
        cout << "NML Cuda platform not found. " << exc.what() << endl;
    }

    cout << "NML looking for Reference plugin..." << endl;
    try {
        Platform& platform = Platform::getPlatformByName("Reference");
        cout << "NML found Reference platform..." << endl;
        ReferenceNMLKernelFactory* factory = new ReferenceNMLKernelFactory();
        platform.registerKernelFactory("IntegrateNMLStep", factory);
    } catch (const std::exception& exc) { // non fatal
        cout << "NML Reference platform not found. " << exc.what() << endl;
    }
}

