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

#include "CpuPlatform.h"
#include "CpuKernelFactory.h"
#include "CpuKernels.h"
#include "openmm/internal/hardware.h"

using namespace OpenMM;

extern "C" OPENMM_EXPORT_CPU void registerPlatforms() {
    // Only register this platform if the CPU supports SSE 4.1.

    int cpuInfo[4];
    cpuid(cpuInfo, 0);
    if (cpuInfo[0] >= 1) {
        cpuid(cpuInfo, 1);
        if ((cpuInfo[2] & ((int) 1 << 19)) != 0)
            Platform::registerPlatform(new CpuPlatform());
    }
}

CpuPlatform::CpuPlatform() {
    CpuKernelFactory* factory = new CpuKernelFactory();
    registerKernelFactory(CalcNonbondedForceKernel::Name(), factory);
}

double CpuPlatform::getSpeed() const {
    return 10;
}

bool CpuPlatform::supportsDoublePrecision() const {
    return false;
}
