/* Portions copyright (c) 2013 Stanford University and Simbios.
 * Authors: Peter Eastman
 * Contributors: 
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

#include "CpuRandom.h"
#include "openmm/internal/OSRngSeed.h"
#include "openmm/OpenMMException.h"
#include <cmath>

using namespace std;
using namespace OpenMM;

CpuRandom::CpuRandom() : hasInitialized(false) {
}

CpuRandom::~CpuRandom() {
    for (auto random : threadRandom)
        delete random;
}

void CpuRandom::initialize(int seed, int numThreads) {
    if (hasInitialized) {
        if (seed == randomSeed)
            return; // Already initialized with the same seed.
        throw OpenMMException("Random number generator initialized twice with different seeds");
    }
    randomSeed = seed;
    hasInitialized = true;
    threadRandom.resize(numThreads);
    nextGaussian.resize(numThreads);
    nextGaussianIsValid.resize(numThreads, false);

    /* Use a quick and dirty RNG to pick seeds for the real random number generator.
     * A random seed of 0 means pick a unique seed
     */

    unsigned int r = (unsigned int) seed;
    if (r == 0)
        r = (unsigned int) osrngseed();
    for (int i = 0; i < numThreads; i++) {
        r = (1664525*r + 1013904223) & 0xFFFFFFFF;
        threadRandom[i] = new OpenMM_SFMT::SFMT();
        init_gen_rand(r, *threadRandom[i]);
    }
}

float CpuRandom::getGaussianRandom(int threadIndex) {
    if (nextGaussianIsValid[threadIndex]) {
        nextGaussianIsValid[threadIndex] = false;
        return nextGaussian[threadIndex];
    }
    
    // Use the polar form of the Box-Muller transformation to generate two Gaussian random numbers.
    
    float x, y, r2;
    do {
        x = 2.0f*(float) genrand_real2(*threadRandom[threadIndex])-1.0f;
        y = 2.0f*(float) genrand_real2(*threadRandom[threadIndex])-1.0f;
        r2 = x*x + y*y;
    } while (r2 >= 1.0f || r2 == 0.0f);
    float multiplier = sqrtf((-2.0f*logf(r2))/r2);
    nextGaussian[threadIndex] = y*multiplier;
    nextGaussianIsValid[threadIndex] = true;
    return x*multiplier;
}

float CpuRandom::getUniformRandom(int threadIndex) {
    return genrand_real2(*threadRandom[threadIndex]);
}
