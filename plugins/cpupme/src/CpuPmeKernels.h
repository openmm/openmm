#ifndef OPENMM_CPU_PME_KERNELS_H_
#define OPENMM_CPU_PME_KERNELS_H_

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

#define NOMINMAX
#include "internal/windowsExportPme.h"
#include "openmm/kernels.h"
#include "openmm/Vec3.h"
#include <fftw3.h>
#include <pthread.h>
#include <vector>

namespace OpenMM {

/**
 * This is an optimized CPU implementation of CalcPmeReciprocalForceKernel.  It is both
 * vectorized (requiring SSE 4.1) and multithreaded.  It uses FFTW to perform the FFTs.
 */

class OPENMM_EXPORT_PME CpuCalcPmeReciprocalForceKernel : public CalcPmeReciprocalForceKernel {
public:
    class ThreadData;
    CpuCalcPmeReciprocalForceKernel(std::string name, const Platform& platform) : CalcPmeReciprocalForceKernel(name, platform),
            hasCreatedPlan(false), isDeleted(false), realGrid(NULL), complexGrid(NULL) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param gridx        the x size of the PME grid
     * @param gridy        the y size of the PME grid
     * @param gridz        the z size of the PME grid
     * @param numParticles the number of particles in the system
     * @param alpha        the Ewald blending parameter
     */
    void initialize(int xsize, int ysize, int zsize, int numParticles, double alpha);
    ~CpuCalcPmeReciprocalForceKernel();
    /**
     * Begin computing the force and energy.
     * 
     * @param io               an object that coordinates data transfer
     * @param periodicBoxSize  the size of the periodic box (measured in nm)
     * @param includeEnergy    true if potential energy should be computed
     */
    void beginComputation(IO& io, Vec3 periodicBoxSize, bool includeEnergy);
    /**
     * Finish computing the force and energy.
     * 
     * @param io   an object that coordinates data transfer
     * @return the potential energy due to the PME reciprocal space interactions
     */
    double finishComputation(IO& io);
    /**
     * This routine contains the code executed by each thread.
     */
    void runThread(int index);
    /**
     * Get whether the current CPU supports all features needed by this kernel.
     */
    static bool isProcessorSupported();
private:
    /**
     * This is called by the worker threads to wait until the master thread instructs them to advance.
     */
    void threadWait();
    /**
     * This is called by the master thread to instruct all the worker threads to advance.
     */
    void advanceThreads();
    /**
     * Select a size for one grid dimension that FFTW can handle efficiently.
     */
    int findFFTDimension(int minimum);
    static bool hasInitializedThreads;
    static int numThreads;
    int gridx, gridy, gridz, numParticles;
    double alpha;
    bool hasCreatedPlan, isFinished, isDeleted;
    std::vector<float> force;
    std::vector<float> bsplineModuli[3];
    std::vector<float> recipEterm;
    Vec3 lastBoxSize;
    float* realGrid;
    fftwf_complex* complexGrid;
    fftwf_plan forwardFFT, backwardFFT;
    int waitCount;
    pthread_cond_t startCondition, endCondition;
    pthread_cond_t mainThreadStartCondition, mainThreadEndCondition;
    pthread_mutex_t lock;
    pthread_t mainThread;
    std::vector<pthread_t> thread;
    std::vector<ThreadData*> threadData;
    // The following variables are used to store information about the calculation currently being performed.
    IO* io;
    float energy;
    float* posq;
    Vec3 periodicBoxSize;
    bool includeEnergy;
};

} // namespace OpenMM

#endif /*OPENMM_CPU_PME_KERNELS_H_*/
