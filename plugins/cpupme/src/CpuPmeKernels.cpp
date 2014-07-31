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

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "CpuPmeKernels.h"
#include "SimTKOpenMMRealType.h"
#include "openmm/internal/hardware.h"
#include "openmm/internal/vectorize.h"
#include <cmath>
#include <cstring>

using namespace OpenMM;
using namespace std;

static const int PME_ORDER = 5;

bool CpuCalcPmeReciprocalForceKernel::hasInitializedThreads = false;
int CpuCalcPmeReciprocalForceKernel::numThreads = 0;

static void spreadCharge(int start, int end, float* posq, float* grid, int gridx, int gridy, int gridz, int numParticles, Vec3 periodicBoxSize) {
    float temp[4];
    fvec4 boxSize((float) periodicBoxSize[0], (float) periodicBoxSize[1], (float) periodicBoxSize[2], 0);
    fvec4 invBoxSize((float) (1/periodicBoxSize[0]), (float) (1/periodicBoxSize[1]), (float) (1/periodicBoxSize[2]), 0);
    fvec4 gridSize(gridx, gridy, gridz, 0);
    ivec4 gridSizeInt(gridx, gridy, gridz, 0);
    fvec4 one(1);
    fvec4 scale(1.0f/(PME_ORDER-1));
    const float epsilonFactor = sqrt(ONE_4PI_EPS0);
    memset(grid, 0, sizeof(float)*gridx*gridy*gridz);
    for (int i = start; i < end; i++) {
        // Find the position relative to the nearest grid point.
        
        fvec4 pos(&posq[4*i]);
        fvec4 posFloor = floor(pos*invBoxSize);
        fvec4 posInBox = pos-boxSize*posFloor;
        fvec4 t = posInBox*invBoxSize*gridSize;
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
    }
}

static void computeReciprocalEterm(int start, int end, int gridx, int gridy, int gridz, vector<float>& recipEterm, double alpha, vector<float>* bsplineModuli, Vec3 periodicBoxSize) {
    const unsigned int zsize = gridz/2+1;
    const unsigned int yzsize = gridy*zsize;
    const float scaleFactor = (float) (M_PI*periodicBoxSize[0]*periodicBoxSize[1]*periodicBoxSize[2]);
    const float recipExpFactor = (float) (M_PI*M_PI/(alpha*alpha));
    const float invPeriodicBoxSizeX = (float) (1.0/periodicBoxSize[0]);
    const float invPeriodicBoxSizeY = (float) (1.0/periodicBoxSize[1]);
    const float invPeriodicBoxSizeZ = (float) (1.0/periodicBoxSize[2]);

    int firstz = (start == 0 ? 1 : 0);
    for (int kx = start; kx < end; kx++) {
        int mx = (kx < (gridx+1)/2) ? kx : kx-gridx;
        float mhx = mx*invPeriodicBoxSizeX;
        float bx = scaleFactor*bsplineModuli[0][kx];
        for (int ky = 0; ky < gridy; ky++) {
            int my = (ky < (gridy+1)/2) ? ky : ky-gridy;
            float mhy = my*invPeriodicBoxSizeY;
            float mhx2y2 = mhx*mhx + mhy*mhy;
            float bxby = bx*bsplineModuli[1][ky];
            for (int kz = firstz; kz < zsize; kz++) {
                int index = kx*yzsize + ky*zsize + kz;
                int mz = (kz < (gridz+1)/2) ? kz : kz-gridz;
                float mhz = mz*invPeriodicBoxSizeZ;
                float bz = bsplineModuli[2][kz];
                float m2 = mhx2y2 + mhz*mhz;
                float denom = m2*bxby*bz;
                recipEterm[index] = exp(-recipExpFactor*m2)/denom;
            }
            firstz = 0;
        }
    }
}

static float reciprocalEnergy(int start, int end, fftwf_complex* grid, int gridx, int gridy, int gridz, double alpha, vector<float>* bsplineModuli, Vec3 periodicBoxSize) {
    const unsigned int zsizeHalf = gridz/2+1;
    const unsigned int yzsizeHalf = gridy*zsizeHalf;
    const float scaleFactor = (float) (M_PI*periodicBoxSize[0]*periodicBoxSize[1]*periodicBoxSize[2]);
    const float recipExpFactor = (float) (M_PI*M_PI/(alpha*alpha));
    const float invPeriodicBoxSizeX = (float) (1.0/periodicBoxSize[0]);
    const float invPeriodicBoxSizeY = (float) (1.0/periodicBoxSize[1]);
    const float invPeriodicBoxSizeZ = (float) (1.0/periodicBoxSize[2]);
    float energy = 0.0f;

    int firstz = (start == 0 ? 1 : 0);
    for (int kx = start; kx < end; kx++) {
        int mx = (kx < (gridx+1)/2) ? kx : kx-gridx;
        float mhx = mx*invPeriodicBoxSizeX;
        float bx = scaleFactor*bsplineModuli[0][kx];
        for (int ky = 0; ky < gridy; ky++) {
            int my = (ky < (gridy+1)/2) ? ky : ky-gridy;
            float mhy = my*invPeriodicBoxSizeY;
            float mhx2y2 = mhx*mhx + mhy*mhy;
            float bxby = bx*bsplineModuli[1][ky];
            for (int kz = firstz; kz < gridz; kz++) {
                int mz = (kz < (gridz+1)/2) ? kz : kz-gridz;
                float mhz = mz*invPeriodicBoxSizeZ;
                float bz = bsplineModuli[2][kz];
                float m2 = mhx2y2 + mhz*mhz;
                float denom = m2*bxby*bz;
                float eterm = exp(-recipExpFactor*m2)/denom;
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
                energy += eterm*(gridReal*gridReal+gridImag*gridImag);
            }
            firstz = 0;
        }
    }
    return 0.5f*energy;
}

static void reciprocalConvolution(int start, int end, fftwf_complex* grid, int gridx, int gridy, int gridz, vector<float>& recipEterm) {
    const unsigned int zsize = gridz/2+1;
    const unsigned int yzsize = gridy*zsize;

    int firstz = (start == 0 ? 1 : 0);
    for (int kx = start; kx < end; kx++) {
        for (int ky = 0; ky < gridy; ky++) {
            for (int kz = firstz; kz < zsize; kz++) {
                int index = kx*yzsize + ky*zsize + kz;
                float eterm = recipEterm[index];
                grid[index][0] *= eterm;
                grid[index][1] *= eterm;
            }
            firstz = 0;
        }
    }
}

static void interpolateForces(int start, int end, float* posq, float* force, float* grid, int gridx, int gridy, int gridz, int numParticles, Vec3 periodicBoxSize) {
    fvec4 boxSize((float) periodicBoxSize[0], (float) periodicBoxSize[1], (float) periodicBoxSize[2], 0);
    fvec4 invBoxSize((float) (1/periodicBoxSize[0]), (float) (1/periodicBoxSize[1]), (float) (1/periodicBoxSize[2]), 0);
    fvec4 gridSize(gridx, gridy, gridz, 0);
    ivec4 gridSizeInt(gridx, gridy, gridz, 0);
    fvec4 one(1);
    fvec4 scale(1.0f/(PME_ORDER-1));
    const float epsilonFactor = sqrt(ONE_4PI_EPS0);
    for (int i = start; i < end; i++) {
        // Find the position relative to the nearest grid point.
        
        fvec4 pos(&posq[4*i]);
        fvec4 posFloor = floor(pos*invBoxSize);
        fvec4 posInBox = pos-boxSize*posFloor;
        fvec4 t = posInBox*invBoxSize*gridSize;
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
        f = invBoxSize*gridSize*f*(-epsilonFactor*posq[4*i+3]);
        f.store(&force[4*i]);
    }
}

class CpuCalcPmeReciprocalForceKernel::ThreadData {
public:
    CpuCalcPmeReciprocalForceKernel& owner;
    int index;
    float* tempGrid;
    ThreadData(CpuCalcPmeReciprocalForceKernel& owner, int index) : owner(owner), index(index), tempGrid(NULL) {
    }
};

static void* threadBody(void* args) {
    CpuCalcPmeReciprocalForceKernel::ThreadData& data = *reinterpret_cast<CpuCalcPmeReciprocalForceKernel::ThreadData*>(args);
    data.owner.runThread(data.index);
    if (data.tempGrid != NULL)
        fftwf_free(data.tempGrid);
    delete &data;
    return 0;
}

void CpuCalcPmeReciprocalForceKernel::initialize(int xsize, int ysize, int zsize, int numParticles, double alpha) {
    if (!hasInitializedThreads) {
        numThreads = getNumProcessors();
        fftwf_init_threads();
        hasInitializedThreads = true;
    }
    gridx = findFFTDimension(xsize);
    gridy = findFFTDimension(ysize);
    gridz = findFFTDimension(zsize);
    this->numParticles = numParticles;
    this->alpha = alpha;
    force.resize(4*numParticles);
    recipEterm.resize(gridx*gridy*gridz);
    
    // Initialize threads.
    
    pthread_cond_init(&startCondition, NULL);
    pthread_cond_init(&endCondition, NULL);
    pthread_cond_init(&mainThreadStartCondition, NULL);
    pthread_cond_init(&mainThreadEndCondition, NULL);
    pthread_mutex_init(&lock, NULL);
    thread.resize(numThreads);
    for (int i = 0; i < numThreads; i++) {
        ThreadData* data = new ThreadData(*this, i);
        threadData.push_back(data);
        pthread_create(&thread[i], NULL, threadBody, data);
        data->tempGrid = (float*) fftwf_malloc(sizeof(float)*(gridx*gridy*gridz+3));
    }
    pthread_create(&mainThread, NULL, threadBody, new ThreadData(*this, -1));
    
    // Initialize FFTW.
    
    realGrid = threadData[0]->tempGrid;
    complexGrid = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*gridx*gridy*(gridz/2+1));
    fftwf_plan_with_nthreads(numThreads);
    forwardFFT = fftwf_plan_dft_r2c_3d(gridx, gridy, gridz, realGrid, complexGrid, FFTW_MEASURE);
    backwardFFT = fftwf_plan_dft_c2r_3d(gridx, gridy, gridz, complexGrid, realGrid, FFTW_MEASURE);
    hasCreatedPlan = true;
    
    // Initialize the b-spline moduli.

    int maxSize = max(max(gridx, gridy), gridz);
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
    pthread_cond_broadcast(&mainThreadStartCondition);
    pthread_mutex_unlock(&lock);
    for (int i = 0; i < (int) thread.size(); i++)
        pthread_join(thread[i], NULL);
    pthread_join(mainThread, NULL);
    pthread_mutex_destroy(&lock);
    pthread_cond_destroy(&startCondition);
    pthread_cond_destroy(&endCondition);
    pthread_cond_destroy(&mainThreadStartCondition);
    pthread_cond_destroy(&mainThreadEndCondition);
    if (complexGrid != NULL)
        fftwf_free(complexGrid);
    if (hasCreatedPlan) {
        fftwf_destroy_plan(forwardFFT);
        fftwf_destroy_plan(backwardFFT);
    }
}

void CpuCalcPmeReciprocalForceKernel::runThread(int index) {
    if (index == -1) {
        // This is the main thread that coordinates all the other ones.
        
        pthread_mutex_lock(&lock);
        while (true) {
            // Wait for the signal to start.
            
            pthread_cond_wait(&mainThreadStartCondition, &lock);
            if (isDeleted)
                break;
            posq = io->getPosq();
            advanceThreads(); // Signal threads to perform charge spreading.
            advanceThreads(); // Signal threads to sum the charge grids.
            fftwf_execute_dft_r2c(forwardFFT, realGrid, complexGrid);
            if (lastBoxSize != periodicBoxSize)
                advanceThreads(); // Signal threads to compute the reciprocal scale factors.
            if (includeEnergy)
                advanceThreads(); // Signal threads to compute energy.
            advanceThreads(); // Signal threads to perform reciprocal convolution.
            fftwf_execute_dft_c2r(backwardFFT, complexGrid, realGrid);
            advanceThreads(); // Signal threads to interpolate forces.
            isFinished = true;
            lastBoxSize = periodicBoxSize;
            pthread_cond_signal(&mainThreadEndCondition);
        }
        pthread_mutex_unlock(&lock);
    }
    else {
        // This is a worker thread.
        
        int particleStart = (index*numParticles)/numThreads;
        int particleEnd = ((index+1)*numParticles)/numThreads;
        int gridxStart = (index*gridx)/numThreads;
        int gridxEnd = ((index+1)*gridx)/numThreads;
        int gridSize = (gridx*gridy*gridz+3)/4;
        int gridStart = 4*((index*gridSize)/numThreads);
        int gridEnd = 4*(((index+1)*gridSize)/numThreads);
        while (true) {
            threadWait();
            if (isDeleted)
                break;
            spreadCharge(particleStart, particleEnd, posq, threadData[index]->tempGrid, gridx, gridy, gridz, numParticles, periodicBoxSize);
            threadWait();
            int numGrids = threadData.size();
            for (int i = gridStart; i < gridEnd; i += 4) {
                fvec4 sum(&realGrid[i]);
                for (int j = 1; j < numGrids; j++)
                    sum += fvec4(&threadData[j]->tempGrid[i]);
                sum.store(&realGrid[i]);
            }
            threadWait();
            if (lastBoxSize != periodicBoxSize) {
                computeReciprocalEterm(gridxStart, gridxEnd, gridx, gridy, gridz, recipEterm, alpha, bsplineModuli, periodicBoxSize);
                threadWait();
            }
            if (includeEnergy) {
                double threadEnergy = reciprocalEnergy(gridxStart, gridxEnd, complexGrid, gridx, gridy, gridz, alpha, bsplineModuli, periodicBoxSize);
                pthread_mutex_lock(&lock);
                energy += threadEnergy;
                pthread_mutex_unlock(&lock);
                threadWait();
            }
            reciprocalConvolution(gridxStart, gridxEnd, complexGrid, gridx, gridy, gridz, recipEterm);
            threadWait();
            interpolateForces(particleStart, particleEnd, posq, &force[0], realGrid, gridx, gridy, gridz, numParticles, periodicBoxSize);
        }
    }
}

void CpuCalcPmeReciprocalForceKernel::threadWait() {
    pthread_mutex_lock(&lock);
    waitCount++;
    pthread_cond_signal(&endCondition);
    pthread_cond_wait(&startCondition, &lock);
    pthread_mutex_unlock(&lock);
}

void CpuCalcPmeReciprocalForceKernel::advanceThreads() {
    waitCount = 0;
    pthread_cond_broadcast(&startCondition);
    while (waitCount < numThreads) {
        pthread_cond_wait(&endCondition, &lock);
    }
}

void CpuCalcPmeReciprocalForceKernel::beginComputation(IO& io, Vec3 periodicBoxSize, bool includeEnergy) {
    this->io = &io;
    this->periodicBoxSize = periodicBoxSize;
    this->includeEnergy = includeEnergy;
    energy = 0.0;
    pthread_mutex_lock(&lock);
    isFinished = false;
    pthread_cond_signal(&mainThreadStartCondition);
    pthread_mutex_unlock(&lock);
}

double CpuCalcPmeReciprocalForceKernel::finishComputation(IO& io) {
    pthread_mutex_lock(&lock);
    while (!isFinished) {
        pthread_cond_wait(&mainThreadEndCondition, &lock);
    }
    pthread_mutex_unlock(&lock);
    io.setForce(&force[0]);
    return energy;
}

bool CpuCalcPmeReciprocalForceKernel::isProcessorSupported() {
    return isVec4Supported();
}

int CpuCalcPmeReciprocalForceKernel::findFFTDimension(int minimum) {
    if (minimum < 1)
        return 1;
    while (true) {
        // Attempt to factor the current value.

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
