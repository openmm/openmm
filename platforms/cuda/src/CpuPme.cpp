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

#include "CpuPme.h"
#include "../src/SimTKUtilities/SimTKOpenMMRealType.h"
#include <smmintrin.h>

using namespace OpenMM;
using namespace std;

static const int PME_ORDER = 5;

static float extract_float(__m128 v, unsigned int element) {
    float f[4];
    _mm_store_ps(f, v);
    return f[element];
}

CpuPme::CpuPme(int gridx, int gridy, int gridz, int numParticles, double alpha) :
        gridx(gridx), gridy(gridy), gridz(gridz), numParticles(numParticles), alpha(alpha) {
    realGrid.resize(gridx*gridy*gridz);
}

CpuPme::~CpuPme() {
}

void spreadCharge(float* posq, float* grid, int gridx, int gridy, int gridz, int numParticles, Vec3 periodicBoxSize) {
    float temp[16];
    __m128 boxSize = _mm_set_ps(0, (float) periodicBoxSize[2], (float) periodicBoxSize[1], (float) periodicBoxSize[0]);
    __m128 invBoxSize = _mm_set_ps(0, (float) (1/periodicBoxSize[2]), (float) (1/periodicBoxSize[1]), (float) (1/periodicBoxSize[0]));
    __m128 gridSize = _mm_set_ps(0, gridz, gridy, gridx);
    __m128 gridSizeInt = _mm_set_epi32(0, gridz, gridy, gridx);
    __m128 one  = _mm_set1_ps(1);
    __m128 scale = _mm_set1_ps(1.0f/(PME_ORDER-1));
    const float epsilonFactor = sqrt(ONE_4PI_EPS0);
    for (int i = 0; i < numParticles; i++) {
        // Find the position relative to the nearest grid point.
        
        __m128 pos = _mm_load_ps(&posq[4*i]);
        __m128 posInBox = _mm_sub_ps(pos, _mm_mul_ps(boxSize, _mm_floor_ps(_mm_mul_ps(pos, invBoxSize))));
        __m128 t = _mm_mul_ps(_mm_mul_ps(posInBox, invBoxSize), gridSize);
        __m128 ti = _mm_cvttps_epi32(t);
        __m128 dr = _mm_sub_ps(t, _mm_cvtepi32_ps(ti));
        __m128 gridIndex = _mm_sub_epi32(ti, _mm_and_si128(gridSizeInt, _mm_cmpeq_epi32(ti, gridSizeInt)));
        
        // Compute the B-spline coefficients.
        
        __m128 data[PME_ORDER];
        data[PME_ORDER-1] = _mm_setzero_ps();
        data[1] = dr;
        data[0] = _mm_sub_ps(one, dr);
        for (int j = 3; j < PME_ORDER; j++) {
            __m128 div = _mm_set1_ps(1.0f/(j-1));
            data[j-1] = _mm_mul_ps(_mm_mul_ps(div, dr), data[j-2]);
            for (int k = 1; k < j-1; k++)
                data[j-k-1] = _mm_mul_ps(div, _mm_add_ps(_mm_mul_ps(_mm_add_ps(dr, _mm_set1_ps(k)), data[j-k-2]), _mm_mul_ps(_mm_sub_ps(_mm_set1_ps(j-k), dr), data[j-k-1])));
            data[0] = _mm_mul_ps(_mm_mul_ps(div, _mm_sub_ps(one, dr)), data[0]);
        }
        data[PME_ORDER-1] = _mm_mul_ps(_mm_mul_ps(scale, dr), data[PME_ORDER-2]);
        for (int j = 1; j < (PME_ORDER-1); j++)
            data[PME_ORDER-j-1] = _mm_mul_ps(scale, _mm_add_ps(_mm_mul_ps(_mm_add_ps(dr, _mm_set1_ps(j)), data[PME_ORDER-j-2]), _mm_mul_ps(_mm_sub_ps(_mm_set1_ps(PME_ORDER-j), dr), data[PME_ORDER-j-1])));
        data[0] = _mm_mul_ps(_mm_mul_ps(scale, _mm_sub_ps(one, dr)), data[0]);
        
        // Spread the charges.
        
        int gridIndexX = _mm_extract_epi32(gridIndex, 0);
        int gridIndexY = _mm_extract_epi32(gridIndex, 1);
        int gridIndexZ = _mm_extract_epi32(gridIndex, 2);
        float charge = epsilonFactor*posq[4*i+3];
        __m128 zdata0to3 = _mm_set_epi32(_mm_extract_ps(data[3], 2), _mm_extract_ps(data[2], 2), _mm_extract_ps(data[1], 2), _mm_extract_ps(data[0], 2));
        float zdata4 = extract_float(data[4], 2);
        for (int ix = 0; ix < PME_ORDER; ix++) {
            int xbase = gridIndexX+ix;
            xbase -= (xbase >= gridx ? gridx : 0);
            xbase = xbase*gridy*gridz;
            float xdata = extract_float(data[ix], 0);

            for (int iy = 0; iy < PME_ORDER; iy++) {
                int ybase = gridIndexY+iy;
                ybase -= (ybase >= gridy ? gridy : 0);
                ybase = xbase + ybase*gridz;
                float multiplier = charge*xdata*extract_float(data[iy], 1);

                __m128 add0to3 = _mm_mul_ps(zdata0to3, _mm_set1_ps(multiplier));
                if (gridIndexZ+4 < gridz)
                    _mm_storeu_ps(&grid[ybase+gridIndexZ], add0to3);
                else {
                    _mm_store_ps(temp, add0to3);
                    int zindex = gridIndexZ;
                    grid[ybase+zindex] += temp[0];
                    zindex++;
                    zindex -= (zindex >= gridz ? gridz : 0);
                    grid[ybase+zindex] += temp[1];
                    zindex++;
                    zindex -= (zindex >= gridz ? gridz : 0);
                    grid[ybase+zindex] += temp[2];
                    zindex++;
                    zindex -= (zindex >= gridz ? gridz : 0);
                    grid[ybase+zindex] += temp[3];
                }
                int zindex = gridIndexZ+4;
                zindex -= (zindex >= gridz ? gridz : 0);
                grid[ybase+zindex] += multiplier*zdata4;
            }
        }
    }
}

#include <sys/time.h>
double diff(struct timeval t1, struct timeval t2) {
    return t2.tv_usec-t1.tv_usec+1e6*(t2.tv_sec-t1.tv_sec);
}

double CpuPme::computeForceAndEnergy(float* posq, float* force, Vec3 periodicBoxSize) {
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);
    spreadCharge(posq, &realGrid[0], gridx, gridy, gridz, numParticles, periodicBoxSize);
    gettimeofday(&t2, NULL);
    printf("time %g\n", diff(t1, t2));
    return 0;
}
