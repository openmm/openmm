/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2018 Stanford University and the Authors.      *
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

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "CpuPmeKernels.h"
#include "SimTKOpenMMRealType.h"
#include "openmm/internal/hardware.h"
#include "openmm/internal/vectorize.h"
#include "openmm/OpenMMException.h"
#include <cmath>
#include <algorithm>
#include <cstring>
#include <sstream>
#include <cstdlib>

using namespace OpenMM;
using namespace std;

static const int PME_ORDER = 5;

bool CpuCalcDispersionPmeReciprocalForceKernel::hasInitializedThreads = false;
int CpuCalcDispersionPmeReciprocalForceKernel::numThreads = 0;

static void spreadCharge(float* posq, float* grid, int gridx, int gridy, int gridz, int numParticles, Vec3* periodicBoxVectors, Vec3* recipBoxVectors,
        atomic<int>& atomicCounter, const float epsilonFactor, int threadIndex, int numThreads, bool deterministic) {
    float temp[4];
    fvec4 boxSize((float) periodicBoxVectors[0][0], (float) periodicBoxVectors[1][1], (float) periodicBoxVectors[2][2], 0);
    fvec4 invBoxSize((float) recipBoxVectors[0][0], (float) recipBoxVectors[1][1], (float) recipBoxVectors[2][2], 0);
    fvec4 recipBoxVec0((float) recipBoxVectors[0][0], (float) recipBoxVectors[0][1], (float) recipBoxVectors[0][2], 0);
    fvec4 recipBoxVec1((float) recipBoxVectors[1][0], (float) recipBoxVectors[1][1], (float) recipBoxVectors[1][2], 0);
    fvec4 recipBoxVec2((float) recipBoxVectors[2][0], (float) recipBoxVectors[2][1], (float) recipBoxVectors[2][2], 0);
    fvec4 gridSize(gridx, gridy, gridz, 0);
    ivec4 gridSizeInt(gridx, gridy, gridz, 0);
    fvec4 one(1);
    fvec4 scale(1.0f/(PME_ORDER-1));
    float posInBox[4] = {0,0,0,0};
    memset(grid, 0, sizeof(float)*gridx*gridy*gridz);

    int i = threadIndex;
    while (true) {
        if (!deterministic)
            i = atomicCounter++;
        if (i >= numParticles)
            break;

        // Find the position relative to the nearest grid point.

        fvec4 pos(&posq[4*i]);
        (pos-boxSize*floor(pos*invBoxSize)).store(posInBox);
        fvec4 t = posInBox[0]*recipBoxVec0 + posInBox[1]*recipBoxVec1 + posInBox[2]*recipBoxVec2;
        t = (t-floor(t))*gridSize;
        ivec4 ti = t;
        fvec4 dr = t-ti;
        ivec4 gridIndex = ti-(gridSizeInt&ti==gridSizeInt);
        
        // Compute the B-spline coefficients.

        fvec4 data[PME_ORDER];
        data[PME_ORDER-1] = 0.0f;
        data[1] = dr;
        data[0] = one-dr;
        for (int j = 3; j < PME_ORDER; j++) {
            fvec4 div(1.0f/(j-1));
            data[j-1] = div*dr*data[j-2];
            for (int k = 1; k < j-1; k++)
                data[j-k-1] = div*((dr+k)*data[j-k-2]+(fvec4(j-k)-dr)*data[j-k-1]);
            data[0] = div*(one-dr)*data[0];
        }
        data[PME_ORDER-1] = scale*dr*data[PME_ORDER-2];
        for (int j = 1; j < (PME_ORDER-1); j++)
            data[PME_ORDER-j-1] = scale*((dr+j)*data[PME_ORDER-j-2]+(fvec4(PME_ORDER-j)-dr)*data[PME_ORDER-j-1]);
        data[0] = scale*(one-dr)*data[0];
        
        // Spread the charges.
        
        int gridIndexX = gridIndex[0];
        int gridIndexY = gridIndex[1];
        int gridIndexZ = gridIndex[2];
        if (gridIndexX < 0)
            return; // This happens when a simulation blows up and coordinates become NaN.
        int zindex[PME_ORDER];
        for (int j = 0; j < PME_ORDER; j++) {
            zindex[j] = gridIndexZ+j;
            zindex[j] -= (zindex[j] >= gridz ? gridz : 0);
        }
        float charge = epsilonFactor*posq[4*i+3];
        fvec4 zdata0to3(data[0][2], data[1][2], data[2][2], data[3][2]);
        float zdata4 = data[4][2];
        if (gridIndexZ+4 < gridz) {
            for (int ix = 0; ix < PME_ORDER; ix++) {
                int xbase = gridIndexX+ix;
                xbase -= (xbase >= gridx ? gridx : 0);
                xbase = xbase*gridy*gridz;
                float xdata = charge*data[ix][0];
                for (int iy = 0; iy < PME_ORDER; iy++) {
                    int ybase = gridIndexY+iy;
                    ybase -= (ybase >= gridy ? gridy : 0);
                    ybase = xbase + ybase*gridz;
                    float multiplier = xdata*data[iy][1];
                    fvec4 add0to3 = zdata0to3*multiplier;
                    (fvec4(&grid[ybase+gridIndexZ])+add0to3).store(&grid[ybase+gridIndexZ]);
                    grid[ybase+zindex[4]] += multiplier*zdata4;
                }
            }
        }
        else {
            for (int ix = 0; ix < PME_ORDER; ix++) {
                int xbase = gridIndexX+ix;
                xbase -= (xbase >= gridx ? gridx : 0);
                xbase = xbase*gridy*gridz;
                float xdata = charge*data[ix][0];
                for (int iy = 0; iy < PME_ORDER; iy++) {
                    int ybase = gridIndexY+iy;
                    ybase -= (ybase >= gridy ? gridy : 0);
                    ybase = xbase + ybase*gridz;
                    float multiplier = xdata*data[iy][1];
                    fvec4 add0to3 = zdata0to3*multiplier;
                    add0to3.store(temp);
                    grid[ybase+zindex[0]] += temp[0];
                    grid[ybase+zindex[1]] += temp[1];
                    grid[ybase+zindex[2]] += temp[2];
                    grid[ybase+zindex[3]] += temp[3];
                    grid[ybase+zindex[4]] += multiplier*zdata4;
                }
            }
        }
        if (deterministic)
            i += numThreads;
    }
}

#define FAST_ERFC 1
static void computeReciprocalDispersionEterm(int start, int end, int gridx, int gridy, int gridz, vector<float>& recipEterm, double alpha, vector<float>* bsplineModuli, Vec3* periodicBoxVectors, Vec3* recipBoxVectors) {
    const unsigned int zsize = gridz/2+1;
    const unsigned int yzsize = gridy*zsize;
    const float scaleFactor = (float)  -2.0f*M_PI*sqrtf(M_PI) / (6.0*periodicBoxVectors[0][0]*periodicBoxVectors[1][1]*periodicBoxVectors[2][2]);

    float bfac = M_PI / alpha;
    float fac1 = 2.0f*M_PI*M_PI*M_PI*sqrtf(M_PI);
    float fac2 = alpha*alpha*alpha;
    float fac3 = -2.0f*alpha*M_PI*M_PI;
    float b, m, m3, expterm, erfcterm, t;

    for (int kx = start; kx < end; kx++) {
        int mx = (kx < (gridx+1)/2) ? kx : kx-gridx;
        float mhx = mx*(float)recipBoxVectors[0][0];
        float bx = bsplineModuli[0][kx];
        for (int ky = 0; ky < gridy; ky++) {
            int my = (ky < (gridy+1)/2) ? ky : ky-gridy;
            float mhy = mx*(float)recipBoxVectors[1][0] + my*(float)recipBoxVectors[1][1];
            float mhx2y2 = mhx*mhx + mhy*mhy;
            float bxby = bx*bsplineModuli[1][ky];
            for (int kz = 0; kz < zsize; kz++) {
                int index = kx*yzsize + ky*zsize + kz;
                int mz = (kz < (gridz+1)/2) ? kz : kz-gridz;
                float mhz = mx*(float)recipBoxVectors[2][0] + my*(float)recipBoxVectors[2][1] + mz*(float)recipBoxVectors[2][2];
                float bz = bsplineModuli[2][kz];
                float m2 = mhx2y2 + mhz*mhz;
                float denom = scaleFactor/(bxby*bz);

                m = sqrtf(m2);
                m3 = m*m2;
                b = bfac*m;
                expterm = exp(-b*b);

#if FAST_ERFC
                // This approximation for erfc is from Abramowitz and Stegun (1964) p. 299.  They cite the following as
                // the original source: C. Hastings, Jr., Approximations for Digital Computers (1955).  It has a maximum
                // error of 1.5e-7.  Stolen by ACS from the CUDA platform's AMOEBA plugin.
                t = 1.0f/(1.0f+0.3275911f*b);
                erfcterm = (0.254829592f+(-0.284496736f+(1.421413741f+(-1.453152027f+1.061405429f*t)*t)*t)*t)*t*expterm;
#else
                erfcterm = erfc(b);
#endif
                recipEterm[index] = (fac1*erfcterm*m3 + expterm*(fac2 + fac3*m2)) * denom;
            }
        }
    }
}

static void computeReciprocalEterm(int start, int end, int gridx, int gridy, int gridz, vector<float>& recipEterm, double alpha, vector<float>* bsplineModuli, Vec3* periodicBoxVectors, Vec3* recipBoxVectors) {
    const unsigned int zsize = gridz/2+1;
    const unsigned int yzsize = gridy*zsize;
    const float scaleFactor = (float) (M_PI*periodicBoxVectors[0][0]*periodicBoxVectors[1][1]*periodicBoxVectors[2][2]);
    const float recipExpFactor = (float) (M_PI*M_PI/(alpha*alpha));

    int firstz = (start == 0 ? 1 : 0);
    for (int kx = start; kx < end; kx++) {
        int mx = (kx < (gridx+1)/2) ? kx : kx-gridx;
        float mhx = mx*(float)recipBoxVectors[0][0];
        float bx = scaleFactor*bsplineModuli[0][kx];
        for (int ky = 0; ky < gridy; ky++) {
            int my = (ky < (gridy+1)/2) ? ky : ky-gridy;
            float mhy = mx*(float)recipBoxVectors[1][0] + my*(float)recipBoxVectors[1][1];
            float mhx2y2 = mhx*mhx + mhy*mhy;
            float bxby = bx*bsplineModuli[1][ky];
            for (int kz = firstz; kz < zsize; kz++) {
                int index = kx*yzsize + ky*zsize + kz;
                int mz = (kz < (gridz+1)/2) ? kz : kz-gridz;
                float mhz = mx*(float)recipBoxVectors[2][0] + my*(float)recipBoxVectors[2][1] + mz*(float)recipBoxVectors[2][2];
                float bz = bsplineModuli[2][kz];
                float m2 = mhx2y2 + mhz*mhz;
                float denom = m2*bxby*bz;
                recipEterm[index] = exp(-recipExpFactor*m2)/denom;
            }
            firstz = 0;
        }
    }
}

static double reciprocalEnergy(int start, int end, fftwf_complex* grid, vector<float>& recipEterm, int gridx, int gridy, int gridz, double alpha, vector<float>* bsplineModuli, Vec3* periodicBoxVectors, Vec3* recipBoxVectors) {
    const unsigned int zsizeHalf = gridz/2+1;
    const unsigned int yzsizeHalf = gridy*zsizeHalf;

    double energy = 0.0;

    int firstz = (start == 0 ? 1 : 0);
    for (int kx = start; kx < end; kx++) {
        for (int ky = 0; ky < gridy; ky++) {
            for (int kz = firstz; kz < gridz; kz++) {
                int kx1, ky1, kz1;
                if (kz >= gridz/2+1) {
                    kx1 = (kx == 0 ? kx : gridx-kx);
                    ky1 = (ky == 0 ? ky : gridy-ky);
                    kz1 = gridz-kz;
                }
                else {
                    kx1 = kx;
                    ky1 = ky;
                    kz1 = kz;
                }
                int index = kx1*yzsizeHalf + ky1*zsizeHalf + kz1;
                float gridReal = grid[index][0];
                float gridImag = grid[index][1];
                energy += recipEterm[index]*(gridReal*gridReal+gridImag*gridImag);
            }
            firstz = 0;
        }
    }
    return 0.5*energy;
}


static double reciprocalDispersionEnergy(int start, int end, fftwf_complex* grid, const vector<float>& recipEterm, int gridx, int gridy, int gridz, double alpha, vector<float>* bsplineModuli, Vec3* periodicBoxVectors, Vec3* recipBoxVectors) {
    const unsigned int zsizeHalf = gridz/2+1;
    const unsigned int yzsizeHalf = gridy*zsizeHalf;

    double energy = 0.0;

    for (int kx = start; kx < end; kx++) {
        for (int ky = 0; ky < gridy; ky++) {
            for (int kz = 0; kz < gridz; kz++) {
                int kx1, ky1, kz1;
                if (kz >= gridz/2+1) {
                    kx1 = (kx == 0 ? kx : gridx-kx);
                    ky1 = (ky == 0 ? ky : gridy-ky);
                    kz1 = gridz-kz;
                }
                else {
                    kx1 = kx;
                    ky1 = ky;
                    kz1 = kz;
                }
                int index = kx1*yzsizeHalf + ky1*zsizeHalf + kz1;
                float gridReal = grid[index][0];
                float gridImag = grid[index][1];
                energy += recipEterm[index]*(gridReal*gridReal+gridImag*gridImag);
            }
        }
    }
    return 0.5f*energy;
}


static void reciprocalConvolution(int start, int end, fftwf_complex* grid, vector<float>& recipEterm) {
    for (int index = start; index < end; index++) {
        float eterm = recipEterm[index];
        grid[index][0] *= eterm;
        grid[index][1] *= eterm;
    }
}

static void interpolateForces(float* posq, float* force, float* grid, int gridx, int gridy, int gridz, int numParticles, Vec3* periodicBoxVectors, Vec3* recipBoxVectors, atomic<int>& atomicCounter, const float epsilonFactor) {
    fvec4 boxSize((float) periodicBoxVectors[0][0], (float) periodicBoxVectors[1][1], (float) periodicBoxVectors[2][2], 0);
    fvec4 invBoxSize((float) recipBoxVectors[0][0], (float) recipBoxVectors[1][1], (float) recipBoxVectors[2][2], 0);
    fvec4 recipBoxVec0((float) recipBoxVectors[0][0], (float) recipBoxVectors[0][1], (float) recipBoxVectors[0][2], 0);
    fvec4 recipBoxVec1((float) recipBoxVectors[1][0], (float) recipBoxVectors[1][1], (float) recipBoxVectors[1][2], 0);
    fvec4 recipBoxVec2((float) recipBoxVectors[2][0], (float) recipBoxVectors[2][1], (float) recipBoxVectors[2][2], 0);
    fvec4 gridSize(gridx, gridy, gridz, 0);
    ivec4 gridSizeInt(gridx, gridy, gridz, 0);
    fvec4 one(1);
    fvec4 scale(1.0f/(PME_ORDER-1));
    while (true) {
        int i = atomicCounter++;
        if (i >= numParticles)
            break;

        // Find the position relative to the nearest grid point.
        
        fvec4 pos(&posq[4*i]);
        float posInBox[4];
        (pos-boxSize*floor(pos*invBoxSize)).store(posInBox);
        fvec4 t = posInBox[0]*recipBoxVec0 + posInBox[1]*recipBoxVec1 + posInBox[2]*recipBoxVec2;
        t = (t-floor(t))*gridSize;
        ivec4 ti = t;
        fvec4 dr = t-ti;
        ivec4 gridIndex = ti-(gridSizeInt&ti==gridSizeInt);
        
        // Compute the B-spline coefficients.
        
        fvec4 data[PME_ORDER];
        fvec4 ddata[PME_ORDER];
        data[PME_ORDER-1] = 0.0f;
        data[1] = dr;
        data[0] = one-dr;
        for (int j = 3; j < PME_ORDER; j++) {
            fvec4 div(1.0f/(j-1));
            data[j-1] = div*dr*data[j-2];
            for (int k = 1; k < j-1; k++)
                data[j-k-1] = div*((dr+k)*data[j-k-2]+(fvec4(j-k)-dr)*data[j-k-1]);
            data[0] = div*(one-dr)*data[0];
        }
        ddata[0] = -data[0];
        for (int j = 1; j < PME_ORDER; j++)
            ddata[j] = data[j-1]-data[j];
        data[PME_ORDER-1] = scale*dr*data[PME_ORDER-2];
        for (int j = 1; j < (PME_ORDER-1); j++)
            data[PME_ORDER-j-1] = scale*((dr+j)*data[PME_ORDER-j-2]+(fvec4(PME_ORDER-j)-dr)*data[PME_ORDER-j-1]);
        data[0] = scale*(one-dr)*data[0];
                
        // Compute the force on this atom.
        
        int gridIndexX = gridIndex[0];
        int gridIndexY = gridIndex[1];
        int gridIndexZ = gridIndex[2];
        if (gridIndexX < 0)
            return; // This happens when a simulation blows up and coordinates become NaN.
        int zindex[PME_ORDER];
        for (int j = 0; j < PME_ORDER; j++) {
            zindex[j] = gridIndexZ+j;
            zindex[j] -= (zindex[j] >= gridz ? gridz : 0);
        }
        fvec4 zdata[PME_ORDER];
        for (int j = 0; j < PME_ORDER; j++)
            zdata[j] = fvec4(data[j][2], data[j][2], ddata[j][2], 0);
        fvec4 f = 0.0f;
        for (int ix = 0; ix < PME_ORDER; ix++) {
            int xbase = gridIndexX+ix;
            xbase -= (xbase >= gridx ? gridx : 0);
            xbase = xbase*gridy*gridz;
            float dx = data[ix][0];
            float ddx = ddata[ix][0];
            fvec4 xdata(ddx, dx, dx, 0);

            for (int iy = 0; iy < PME_ORDER; iy++) {
                int ybase = gridIndexY+iy;
                ybase -= (ybase >= gridy ? gridy : 0);
                ybase = xbase + ybase*gridz;
                float dy = data[iy][1];
                float ddy = ddata[iy][1];
                fvec4 xydata = xdata*fvec4(dy, ddy, dy, 0);

                for (int iz = 0; iz < PME_ORDER; iz++) {
                    fvec4 gridValue(grid[ybase+zindex[iz]]);
                    f = f+xydata*zdata[iz]*gridValue;
                }
            }
        }
        f *= -epsilonFactor*posq[4*i+3];
        float fc[4];
        f.store(fc);
        force[4*i+0] = fc[0]*gridx*(float)recipBoxVectors[0][0];
        force[4*i+1] = fc[0]*gridx*(float)recipBoxVectors[1][0]+fc[1]*gridy*(float)recipBoxVectors[1][1];
        force[4*i+2] = fc[0]*gridx*(float)recipBoxVectors[2][0]+fc[1]*gridy*(float)recipBoxVectors[2][1]+fc[2]*gridz*(float)recipBoxVectors[2][2];
    }
}

static void* threadBody(void* args) {
    CpuCalcPmeReciprocalForceKernel& owner = *reinterpret_cast<CpuCalcPmeReciprocalForceKernel*>(args);
    owner.runMainThread();
    return 0;
}

void CpuCalcPmeReciprocalForceKernel::initialize(int xsize, int ysize, int zsize, int numParticles, double alpha, bool deterministic) {
    if (!hasInitializedThreads) {
        numThreads = getNumProcessors();
        char* threadsEnv = getenv("OPENMM_CPU_THREADS");
        if (threadsEnv != NULL)
            stringstream(threadsEnv) >> numThreads;
        fftwf_init_threads();
        hasInitializedThreads = true;
    }
    threadEnergy.resize(numThreads);
    gridx = findFFTDimension(xsize, false);
    gridy = findFFTDimension(ysize, false);
    gridz = findFFTDimension(zsize, true);
    this->numParticles = numParticles;
    this->alpha = alpha;
    this->deterministic = deterministic;
    force.resize(4*numParticles);
    recipEterm.resize(gridx*gridy*gridz);
    
    // Initialize threads.
    
    isFinished = false;
    pthread_cond_init(&startCondition, NULL);
    pthread_cond_init(&endCondition, NULL);
    pthread_mutex_init(&lock, NULL);
    pthread_create(&mainThread, NULL, threadBody, this);
    
    // Wait until the main thread is up and running.
    
    pthread_mutex_lock(&lock);
    while (!isFinished)
        pthread_cond_wait(&endCondition, &lock);
    pthread_mutex_unlock(&lock);
    
    // Initialize FFTW.
    
    for (int i = 0; i < numThreads; i++)
        tempGrid.push_back((float*) fftwf_malloc(sizeof(float)*(gridx*gridy*gridz+3)));
    realGrid = tempGrid[0];
    complexGrid = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*gridx*gridy*(gridz/2+1));
    fftwf_plan_with_nthreads(numThreads);
    forwardFFT = fftwf_plan_dft_r2c_3d(gridx, gridy, gridz, realGrid, complexGrid, FFTW_MEASURE);
    backwardFFT = fftwf_plan_dft_c2r_3d(gridx, gridy, gridz, complexGrid, realGrid, FFTW_MEASURE);
    hasCreatedPlan = true;
    
    // Initialize the b-spline moduli.

    int maxSize = std::max(std::max(gridx, gridy), gridz);
    vector<double> data(PME_ORDER);
    vector<double> ddata(PME_ORDER);
    vector<double> bsplinesData(maxSize);
    data[PME_ORDER-1] = 0.0;
    data[1] = 0.0;
    data[0] = 1.0;
    for (int i = 3; i < PME_ORDER; i++) {
        double div = 1.0/(i-1.0);
        data[i-1] = 0.0;
        for (int j = 1; j < (i-1); j++)
            data[i-j-1] = div*(j*data[i-j-2]+(i-j)*data[i-j-1]);
        data[0] = div*data[0];
    }

    // Differentiate.

    ddata[0] = -data[0];
    for (int i = 1; i < PME_ORDER; i++)
        ddata[i] = data[i-1]-data[i];
    double div = 1.0/(PME_ORDER-1);
    data[PME_ORDER-1] = 0.0;
    for (int i = 1; i < (PME_ORDER-1); i++)
        data[PME_ORDER-i-1] = div*(i*data[PME_ORDER-i-2]+(PME_ORDER-i)*data[PME_ORDER-i-1]);
    data[0] = div*data[0];
    for (int i = 0; i < maxSize; i++)
        bsplinesData[i] = 0.0;
    for (int i = 1; i <= PME_ORDER; i++)
        bsplinesData[i] = data[i-1];

    // Evaluate the actual bspline moduli for X/Y/Z.

    bsplineModuli[0].resize(gridx);
    bsplineModuli[1].resize(gridy);
    bsplineModuli[2].resize(gridz);
    for (int dim = 0; dim < 3; dim++) {
        int ndata = bsplineModuli[dim].size();
        vector<float>& moduli = bsplineModuli[dim];
        for (int i = 0; i < ndata; i++) {
            double sc = 0.0;
            double ss = 0.0;
            for (int j = 0; j < ndata; j++) {
                double arg = (2.0*M_PI*i*j)/ndata;
                sc += bsplinesData[j]*cos(arg);
                ss += bsplinesData[j]*sin(arg);
            }
            moduli[i] = (float) (sc*sc+ss*ss);
        }
        for (int i = 0; i < ndata; i++)
            if (moduli[i] < 1.0e-7f)
                moduli[i] = (moduli[i-1]+moduli[i+1])*0.5f;
    }
}

CpuCalcPmeReciprocalForceKernel::~CpuCalcPmeReciprocalForceKernel() {
    isDeleted = true;
    pthread_mutex_lock(&lock);
    pthread_cond_broadcast(&startCondition);
    pthread_mutex_unlock(&lock);
    pthread_join(mainThread, NULL);
    pthread_mutex_destroy(&lock);
    pthread_cond_destroy(&startCondition);
    pthread_cond_destroy(&endCondition);
    for (auto grid : tempGrid)
        fftwf_free(grid);
    if (complexGrid != NULL)
        fftwf_free(complexGrid);
    if (hasCreatedPlan) {
        fftwf_destroy_plan(forwardFFT);
        fftwf_destroy_plan(backwardFFT);
    }
}

void CpuCalcPmeReciprocalForceKernel::runMainThread() {
    // This is the main thread that coordinates all the other ones.

    pthread_mutex_lock(&lock);
    isFinished = true;
    pthread_cond_signal(&endCondition);
    ThreadPool threads(numThreads);
    while (true) {
        // Wait for the signal to start.

        pthread_cond_wait(&startCondition, &lock);
        if (isDeleted)
            break;
        posq = io->getPosq();
        atomicCounter = 0;
        threads.execute([&] (ThreadPool& threads, int threadIndex) { runWorkerThread(threads, threadIndex); }); // Signal threads to perform charge spreading.
        threads.waitForThreads();
        threads.resumeThreads(); // Signal threads to sum the charge grids.
        threads.waitForThreads();
        fftwf_execute_dft_r2c(forwardFFT, realGrid, complexGrid);
        if (lastBoxVectors[0] != periodicBoxVectors[0] || lastBoxVectors[1] != periodicBoxVectors[1] || lastBoxVectors[2] != periodicBoxVectors[2]) {
            threads.resumeThreads(); // Signal threads to compute the reciprocal scale factors.
            threads.waitForThreads();
        }
        if (includeEnergy) {
            threads.resumeThreads(); // Signal threads to compute energy.
            threads.waitForThreads();
            for (auto e : threadEnergy)
                energy += e;
        }
        threads.resumeThreads(); // Signal threads to perform reciprocal convolution.
        threads.waitForThreads();
        fftwf_execute_dft_c2r(backwardFFT, complexGrid, realGrid);
        atomicCounter = 0;
        threads.resumeThreads(); // Signal threads to interpolate forces.
        threads.waitForThreads();
        isFinished = true;
        lastBoxVectors[0] = periodicBoxVectors[0];
        lastBoxVectors[1] = periodicBoxVectors[1];
        lastBoxVectors[2] = periodicBoxVectors[2];
        pthread_cond_signal(&endCondition);
    }
    pthread_mutex_unlock(&lock);
}

void CpuCalcPmeReciprocalForceKernel::runWorkerThread(ThreadPool& threads, int index) {
    int gridxStart = (index*gridx)/numThreads;
    int gridxEnd = ((index+1)*gridx)/numThreads;
    int gridSize = (gridx*gridy*gridz+3)/4;
    int gridStart = 4*((index*gridSize)/numThreads);
    int gridEnd = 4*(((index+1)*gridSize)/numThreads);
    int complexSize = gridx*gridy*(gridz/2+1);
    int complexStart = std::max(1, ((index*complexSize)/numThreads));
    int complexEnd = (((index+1)*complexSize)/numThreads);
    const float epsilonFactor = sqrt(ONE_4PI_EPS0);
    spreadCharge(posq, tempGrid[index], gridx, gridy, gridz, numParticles, periodicBoxVectors, recipBoxVectors, atomicCounter, epsilonFactor, index, numThreads, deterministic);
    threads.syncThreads();
    int numGrids = tempGrid.size();
    for (int i = gridStart; i < gridEnd; i += 4) {
        fvec4 sum(&realGrid[i]);
        for (int j = 1; j < numGrids; j++)
            sum += fvec4(&tempGrid[j][i]);
        sum.store(&realGrid[i]);
    }
    threads.syncThreads();
    if (lastBoxVectors[0] != periodicBoxVectors[0] || lastBoxVectors[1] != periodicBoxVectors[1] || lastBoxVectors[2] != periodicBoxVectors[2]) {
        computeReciprocalEterm(gridxStart, gridxEnd, gridx, gridy, gridz, recipEterm, alpha, bsplineModuli, periodicBoxVectors, recipBoxVectors);
        threads.syncThreads();
    }
    if (includeEnergy) {
        threadEnergy[index] = reciprocalEnergy(gridxStart, gridxEnd, complexGrid, recipEterm, gridx, gridy, gridz, alpha, bsplineModuli, periodicBoxVectors, recipBoxVectors);
        threads.syncThreads();
    }
    reciprocalConvolution(complexStart, complexEnd, complexGrid, recipEterm);
    threads.syncThreads();
    interpolateForces(posq, &force[0], realGrid, gridx, gridy, gridz, numParticles, periodicBoxVectors, recipBoxVectors, atomicCounter, epsilonFactor);
}

void CpuCalcPmeReciprocalForceKernel::beginComputation(IO& io, const Vec3* periodicBoxVectors, bool includeEnergy) {
    this->io = &io;
    this->periodicBoxVectors[0] = periodicBoxVectors[0];
    this->periodicBoxVectors[1] = periodicBoxVectors[1];
    this->periodicBoxVectors[2] = periodicBoxVectors[2];
    this->includeEnergy = includeEnergy;
    energy = 0.0;

    // Invert the box vectors.

    double determinant = periodicBoxVectors[0][0]*periodicBoxVectors[1][1]*periodicBoxVectors[2][2];
    double scale = 1.0/determinant;
    recipBoxVectors[0] = Vec3(periodicBoxVectors[1][1]*periodicBoxVectors[2][2], 0, 0)*scale;
    recipBoxVectors[1] = Vec3(-periodicBoxVectors[1][0]*periodicBoxVectors[2][2], periodicBoxVectors[0][0]*periodicBoxVectors[2][2], 0)*scale;
    recipBoxVectors[2] = Vec3(periodicBoxVectors[1][0]*periodicBoxVectors[2][1]-periodicBoxVectors[1][1]*periodicBoxVectors[2][0], -periodicBoxVectors[0][0]*periodicBoxVectors[2][1], periodicBoxVectors[0][0]*periodicBoxVectors[1][1])*scale;

    // Do the calculation.

    pthread_mutex_lock(&lock);
    isFinished = false;
    pthread_cond_signal(&startCondition);
    pthread_mutex_unlock(&lock);
}

double CpuCalcPmeReciprocalForceKernel::finishComputation(IO& io) {
    pthread_mutex_lock(&lock);
    while (!isFinished) {
        pthread_cond_wait(&endCondition, &lock);
    }
    pthread_mutex_unlock(&lock);
    io.setForce(&force[0]);
    return energy;
}

bool CpuCalcPmeReciprocalForceKernel::isProcessorSupported() {
    return isVec4Supported();
}

void CpuCalcPmeReciprocalForceKernel::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    alpha = this->alpha;
    nx = gridx;
    ny = gridy;
    nz = gridz;
}

int CpuCalcPmeReciprocalForceKernel::findFFTDimension(int minimum, bool isZ) {
    if (minimum < 1)
        return 1;
    while (true) {
        // Attempt to factor the current value.

        if (isZ && minimum%2 == 1) {
            // Force the last dimension to be even, since this produces better performance in FFTW.

            minimum++;
            continue;
        }
        int unfactored = minimum;
        for (int factor = 2; factor < 8; factor++) {
            while (unfactored > 1 && unfactored%factor == 0)
                unfactored /= factor;
        }
        if (unfactored == 1 || unfactored == 11 || unfactored == 13)
            return minimum;
        minimum++;
    }
}

/*
 * Everything below here is just a clone of the above, but to handle the dispersion term
 * instead of electrostatics.
 */

bool CpuCalcPmeReciprocalForceKernel::hasInitializedThreads = false;
int CpuCalcPmeReciprocalForceKernel::numThreads = 0;


class CpuCalcDispersionPmeReciprocalForceKernel::ComputeTask : public ThreadPool::Task {
public:
    ComputeTask(CpuCalcDispersionPmeReciprocalForceKernel& owner) : owner(owner) {
    }
    void execute(ThreadPool& threads, int threadIndex) {
        owner.runWorkerThread(threads, threadIndex);
    }
    CpuCalcDispersionPmeReciprocalForceKernel& owner;
};

static void* dispersionThreadBody(void* args) {
    CpuCalcDispersionPmeReciprocalForceKernel& owner = *reinterpret_cast<CpuCalcDispersionPmeReciprocalForceKernel*>(args);
    owner.runMainThread();
    return 0;
}

void CpuCalcDispersionPmeReciprocalForceKernel::initialize(int xsize, int ysize, int zsize, int numParticles, double alpha, bool deterministic) {
    if (!hasInitializedThreads) {
        numThreads = getNumProcessors();
        char* threadsEnv = getenv("OPENMM_CPU_THREADS");
        if (threadsEnv != NULL)
            stringstream(threadsEnv) >> numThreads;
        fftwf_init_threads();
        hasInitializedThreads = true;
    }
    threadEnergy.resize(numThreads);
    gridx = findFFTDimension(xsize, false);
    gridy = findFFTDimension(ysize, false);
    gridz = findFFTDimension(zsize, true);
    this->numParticles = numParticles;
    this->alpha = alpha;
    this->deterministic = deterministic;
    force.resize(4*numParticles);
    recipEterm.resize(gridx*gridy*gridz);
    
    // Initialize threads.
    
    isFinished = false;
    pthread_cond_init(&startCondition, NULL);
    pthread_cond_init(&endCondition, NULL);
    pthread_mutex_init(&lock, NULL);
    pthread_create(&mainThread, NULL, dispersionThreadBody, this);
    
    // Wait until the main thread is up and running.
    
    pthread_mutex_lock(&lock);
    while (!isFinished)
        pthread_cond_wait(&endCondition, &lock);
    pthread_mutex_unlock(&lock);
    
    // Initialize FFTW.
    
    for (int i = 0; i < numThreads; i++)
        tempGrid.push_back((float*) fftwf_malloc(sizeof(float)*(gridx*gridy*gridz+3)));
    realGrid = tempGrid[0];
    complexGrid = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*gridx*gridy*(gridz/2+1));
    fftwf_plan_with_nthreads(numThreads);
    forwardFFT = fftwf_plan_dft_r2c_3d(gridx, gridy, gridz, realGrid, complexGrid, FFTW_MEASURE);
    backwardFFT = fftwf_plan_dft_c2r_3d(gridx, gridy, gridz, complexGrid, realGrid, FFTW_MEASURE);
    hasCreatedPlan = true;
    
    // Initialize the b-spline moduli.

    int maxSize = std::max(std::max(gridx, gridy), gridz);
    vector<double> data(PME_ORDER);
    vector<double> ddata(PME_ORDER);
    vector<double> bsplinesData(maxSize);
    data[PME_ORDER-1] = 0.0;
    data[1] = 0.0;
    data[0] = 1.0;
    for (int i = 3; i < PME_ORDER; i++) {
        double div = 1.0/(i-1.0);
        data[i-1] = 0.0;
        for (int j = 1; j < (i-1); j++)
            data[i-j-1] = div*(j*data[i-j-2]+(i-j)*data[i-j-1]);
        data[0] = div*data[0];
    }

    // Differentiate.

    ddata[0] = -data[0];
    for (int i = 1; i < PME_ORDER; i++)
        ddata[i] = data[i-1]-data[i];
    double div = 1.0/(PME_ORDER-1);
    data[PME_ORDER-1] = 0.0;
    for (int i = 1; i < (PME_ORDER-1); i++)
        data[PME_ORDER-i-1] = div*(i*data[PME_ORDER-i-2]+(PME_ORDER-i)*data[PME_ORDER-i-1]);
    data[0] = div*data[0];
    for (int i = 0; i < maxSize; i++)
        bsplinesData[i] = 0.0;
    for (int i = 1; i <= PME_ORDER; i++)
        bsplinesData[i] = data[i-1];

    // Evaluate the actual bspline moduli for X/Y/Z.

    bsplineModuli[0].resize(gridx);
    bsplineModuli[1].resize(gridy);
    bsplineModuli[2].resize(gridz);
    for (int dim = 0; dim < 3; dim++) {
        int ndata = bsplineModuli[dim].size();
        vector<float>& moduli = bsplineModuli[dim];
        for (int i = 0; i < ndata; i++) {
            double sc = 0.0;
            double ss = 0.0;
            for (int j = 0; j < ndata; j++) {
                double arg = (2.0*M_PI*i*j)/ndata;
                sc += bsplinesData[j]*cos(arg);
                ss += bsplinesData[j]*sin(arg);
            }
            moduli[i] = (float) (sc*sc+ss*ss);
        }
        for (int i = 0; i < ndata; i++)
            if (moduli[i] < 1.0e-7f)
                moduli[i] = (moduli[i-1]+moduli[i+1])*0.5f;
    }
}

CpuCalcDispersionPmeReciprocalForceKernel::~CpuCalcDispersionPmeReciprocalForceKernel() {
    isDeleted = true;
    pthread_mutex_lock(&lock);
    pthread_cond_broadcast(&startCondition);
    pthread_mutex_unlock(&lock);
    pthread_join(mainThread, NULL);
    pthread_mutex_destroy(&lock);
    pthread_cond_destroy(&startCondition);
    pthread_cond_destroy(&endCondition);
    for (auto grid : tempGrid)
        fftwf_free(grid);
    if (complexGrid != NULL)
        fftwf_free(complexGrid);
    if (hasCreatedPlan) {
        fftwf_destroy_plan(forwardFFT);
        fftwf_destroy_plan(backwardFFT);
    }
}

void CpuCalcDispersionPmeReciprocalForceKernel::runMainThread() {
    // This is the main thread that coordinates all the other ones.

    pthread_mutex_lock(&lock);
    isFinished = true;
    pthread_cond_signal(&endCondition);
    ThreadPool threads(numThreads);
    while (true) {
        // Wait for the signal to start.

        pthread_cond_wait(&startCondition, &lock);
        if (isDeleted)
            break;
        posq = io->getPosq();
        ComputeTask task(*this);
        atomicCounter = 0;
        threads.execute(task); // Signal threads to perform charge spreading.
        threads.waitForThreads();
        threads.resumeThreads(); // Signal threads to sum the charge grids.
        threads.waitForThreads();
        fftwf_execute_dft_r2c(forwardFFT, realGrid, complexGrid);
        if (lastBoxVectors[0] != periodicBoxVectors[0] || lastBoxVectors[1] != periodicBoxVectors[1] || lastBoxVectors[2] != periodicBoxVectors[2]) {
            threads.resumeThreads(); // Signal threads to compute the reciprocal scale factors.
            threads.waitForThreads();
        }
        if (includeEnergy) {
            threads.resumeThreads(); // Signal threads to compute energy.
            threads.waitForThreads();
            for (auto e : threadEnergy)
                energy += e;
        }
        threads.resumeThreads(); // Signal threads to perform reciprocal convolution.
        threads.waitForThreads();
        fftwf_execute_dft_c2r(backwardFFT, complexGrid, realGrid);
        atomicCounter = 0;
        threads.resumeThreads(); // Signal threads to interpolate forces.
        threads.waitForThreads();
        isFinished = true;
        lastBoxVectors[0] = periodicBoxVectors[0];
        lastBoxVectors[1] = periodicBoxVectors[1];
        lastBoxVectors[2] = periodicBoxVectors[2];
        pthread_cond_signal(&endCondition);
    }
    pthread_mutex_unlock(&lock);
}

void CpuCalcDispersionPmeReciprocalForceKernel::runWorkerThread(ThreadPool& threads, int index) {
    int gridxStart = (index*gridx)/numThreads;
    int gridxEnd = ((index+1)*gridx)/numThreads;
    int gridSize = (gridx*gridy*gridz+3)/4;
    int gridStart = 4*((index*gridSize)/numThreads);
    int gridEnd = 4*(((index+1)*gridSize)/numThreads);
    int complexSize = gridx*gridy*(gridz/2+1);
    int complexStart = std::max(1, ((index*complexSize)/numThreads));
    int complexEnd = (((index+1)*complexSize)/numThreads);
    const float epsilonFactor = 1.0f;
    spreadCharge(posq, tempGrid[index], gridx, gridy, gridz, numParticles, periodicBoxVectors, recipBoxVectors, atomicCounter, epsilonFactor, index, numThreads, deterministic);
    threads.syncThreads();
    int numGrids = tempGrid.size();
    for (int i = gridStart; i < gridEnd; i += 4) {
        fvec4 sum(&realGrid[i]);
        for (int j = 1; j < numGrids; j++)
            sum += fvec4(&tempGrid[j][i]);
        sum.store(&realGrid[i]);
    }
    threads.syncThreads();
    if (lastBoxVectors[0] != periodicBoxVectors[0] || lastBoxVectors[1] != periodicBoxVectors[1] || lastBoxVectors[2] != periodicBoxVectors[2]) {
        computeReciprocalDispersionEterm(gridxStart, gridxEnd, gridx, gridy, gridz, recipEterm, alpha, bsplineModuli, periodicBoxVectors, recipBoxVectors);
        threads.syncThreads();
    }
    if (includeEnergy) {
        threadEnergy[index] = reciprocalDispersionEnergy(gridxStart, gridxEnd, complexGrid, recipEterm, gridx, gridy, gridz, alpha, bsplineModuli, periodicBoxVectors, recipBoxVectors);
        threads.syncThreads();
    }
    // For dispersion, we include the {0,0,0} term, so the start point needs to be redefined
    complexStart = (index*complexSize)/numThreads;
    reciprocalConvolution(complexStart, complexEnd, complexGrid, recipEterm);
    threads.syncThreads();
    interpolateForces(posq, &force[0], realGrid, gridx, gridy, gridz, numParticles, periodicBoxVectors, recipBoxVectors, atomicCounter, epsilonFactor);
}

void CpuCalcDispersionPmeReciprocalForceKernel::beginComputation(CalcPmeReciprocalForceKernel::IO& io, const Vec3* periodicBoxVectors, bool includeEnergy) {
    this->io = &io;
    this->periodicBoxVectors[0] = periodicBoxVectors[0];
    this->periodicBoxVectors[1] = periodicBoxVectors[1];
    this->periodicBoxVectors[2] = periodicBoxVectors[2];
    this->includeEnergy = includeEnergy;
    energy = 0.0;

    // Invert the box vectors.

    double determinant = periodicBoxVectors[0][0]*periodicBoxVectors[1][1]*periodicBoxVectors[2][2];
    double scale = 1.0/determinant;
    recipBoxVectors[0] = Vec3(periodicBoxVectors[1][1]*periodicBoxVectors[2][2], 0, 0)*scale;
    recipBoxVectors[1] = Vec3(-periodicBoxVectors[1][0]*periodicBoxVectors[2][2], periodicBoxVectors[0][0]*periodicBoxVectors[2][2], 0)*scale;
    recipBoxVectors[2] = Vec3(periodicBoxVectors[1][0]*periodicBoxVectors[2][1]-periodicBoxVectors[1][1]*periodicBoxVectors[2][0], -periodicBoxVectors[0][0]*periodicBoxVectors[2][1], periodicBoxVectors[0][0]*periodicBoxVectors[1][1])*scale;

    // Do the calculation.

    pthread_mutex_lock(&lock);
    isFinished = false;
    pthread_cond_signal(&startCondition);
    pthread_mutex_unlock(&lock);
}

double CpuCalcDispersionPmeReciprocalForceKernel::finishComputation(CalcPmeReciprocalForceKernel::IO& io) {
    pthread_mutex_lock(&lock);
    while (!isFinished) {
        pthread_cond_wait(&endCondition, &lock);
    }
    pthread_mutex_unlock(&lock);
    io.setForce(&force[0]);
    return energy;
}

bool CpuCalcDispersionPmeReciprocalForceKernel::isProcessorSupported() {
    return isVec4Supported();
}

void CpuCalcDispersionPmeReciprocalForceKernel::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    alpha = this->alpha;
    nx = gridx;
    ny = gridy;
    nz = gridz;
}

int CpuCalcDispersionPmeReciprocalForceKernel::findFFTDimension(int minimum, bool isZ) {
    if (minimum < 1)
        return 1;
    while (true) {
        // Attempt to factor the current value.

        if (isZ && minimum%2 == 1) {
            // Force the last dimension to be even, since this produces better performance in FFTW.

            minimum++;
            continue;
        }
        int unfactored = minimum;
        for (int factor = 2; factor < 8; factor++) {
            while (unfactored > 1 && unfactored%factor == 0)
                unfactored /= factor;
        }
        if (unfactored == 1 || unfactored == 11 || unfactored == 13)
            return minimum;
        minimum++;
    }
}
