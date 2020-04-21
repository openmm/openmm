
/* Portions copyright (c) 2006-2015 Stanford University and Simbios.
 * Contributors: Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "CpuNonbondedForceFvec.h"

#include "openmm/internal/vectorize8.h"

#include <iostream>

using OpenMM::CpuNonbondedForce;

CpuNonbondedForce* createCpuNonbondedForceSse();

CpuNonbondedForce* createCpuNonbondedForceAvx() {
    return new OpenMM::CpuNonbondedForceFvec<fvec8>();
}

class CpuProperties
{
public:
    CpuProperties();
    int numBlockElements = 0;
    OpenMM::CpuNonbondedForce* cpuForce = 0;
};

static const CpuProperties cpuProperties;

CpuProperties::CpuProperties()
{
    bool supportsSse = false;
    bool supportsAvx = false;
    bool supportsAvx2 = false;

    int cpuInfo[4];
    cpuid(cpuInfo, 0);
    if (cpuInfo[0] >= 1) {
        cpuid(cpuInfo, 1);

        supportsSse = ((cpuInfo[2] & ((int) 1 << 19)) != 0);

#ifdef __AVX__
        supportsAvx = ((cpuInfo[2] & ((int) 1 << 28)) != 0);
#endif
    }

    if (supportsAvx) {
        cpuForce = createCpuNonbondedForceAvx();
        numBlockElements = 8;
    }
    else if (supportsSse)
    {
        cpuForce = createCpuNonbondedForceSse();
        numBlockElements = 4;
    }

    std::cout << "SSE: " << supportsSse << std::endl;
    std::cout << "AVX: " << supportsAvx << std::endl;
    std::cout << "AVX2: " << supportsAvx2 << std::endl;

}

OpenMM::CpuNonbondedForce* createCpuNonbondedForceVec() {
    return cpuProperties.cpuForce;
}

int getVecBlockSize() {
    return cpuProperties.numBlockElements;
}
