
/* Portions copyright (c) 2006-2015 Stanford University and Simbios.
 * Contributors: Daniel Towner
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
#include "openmm/OpenMMException.h"

#ifdef __AVX2__

#include "openmm/internal/vectorizeAvx2.h"
OpenMM::CpuNonbondedForce* createCpuNonbondedForceAvx2() {
    return new OpenMM::CpuNonbondedForceFvec<fvecAvx2>();
}

#else

bool isAvx2Supported() {
    return false;
}

OpenMM::CpuNonbondedForce* createCpuNonbondedForceAvx2() {
   throw OpenMM::OpenMMException("Internal error: OpenMM was compiled without AVX2 support");
}
#endif
