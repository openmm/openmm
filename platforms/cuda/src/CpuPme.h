#ifndef OPENMM_CPU_PME_H_
#define OPENMM_CPU_PME_H_

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

#include "windowsExportCuda.h"
#include "openmm/Vec3.h"
#include <complex>
#include <fftw3.h>
#include <pthread.h>
#include <vector>

namespace OpenMM {

/**
 */

class OPENMM_EXPORT_CUDA CpuPme {
public:
    class ThreadData;
    /**
     */
    CpuPme(int gridx, int gridy, int gridz, int numParticles, double alpha);
    ~CpuPme();
    double computeForceAndEnergy(float* posq, float* force, Vec3 periodicBoxSize, bool includeEnergy);
    void runThread(int index);
private:
    void threadWait();
    void advanceThreads();
    static bool hasInitializedThreads;
    static int numThreads;
    int gridx, gridy, gridz, numParticles;
    double alpha;
    bool hasCreatedPlan, isFinished;
    std::vector<float> bsplineModuli[3];
    float* realGrid;
    fftwf_complex* complexGrid;
    fftwf_plan forwardFFT, backwardFFT;
    int waitCount;
    pthread_cond_t startCondition, endCondition;
    pthread_mutex_t lock;
    std::vector<pthread_t> thread;
    std::vector<ThreadData*> threadData;
    // The following variables are used to store information about the calculation currently being performed.
    float energy;
    float* posq;
    float* force;
    Vec3 periodicBoxSize;
    bool includeEnergy;
};

} // namespace OpenMM

#endif /*OPENMM_CPU_PME_H_*/
