#ifndef OPENMM_CUDAPLATFORM_H_
#define OPENMM_CUDAPLATFORM_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "openmm/Platform.h"
#include "CudaStreamFactory.h"

struct _gpuContext;

namespace OpenMM {
    
class KernelImpl;

/**
 * This Platform subclass uses CUDA implementations of the OpenMM kernels to run on NVidia GPUs.
 */

class OPENMM_EXPORT CudaPlatform : public Platform {
public:
    class PlatformData;
    CudaPlatform();
    std::string getName() const {
        return "Cuda";
    }
    double getSpeed() const {
        return 100;
    }
    bool supportsDoublePrecision() const;
    const StreamFactory& getDefaultStreamFactory() const;
    void contextCreated(OpenMMContextImpl& context) const;
    void contextDestroyed(OpenMMContextImpl& context) const;
private:
    CudaStreamFactory defaultStreamFactory;
};

class CudaPlatform::PlatformData {
public:
    PlatformData(_gpuContext* gpu) : gpu(gpu), removeCM(false), nonbondedMethod(0), hasBonds(false), hasAngles(false),
            hasPeriodicTorsions(false), hasRB(false), hasNonbonded(false), primaryKernel(NULL), stepCount(0) {
    }
    _gpuContext* gpu;
    KernelImpl* primaryKernel;
    bool removeCM;
    bool hasBonds, hasAngles, hasPeriodicTorsions, hasRB, hasNonbonded;
    int nonbondedMethod;
    int cmMotionFrequency;
    int stepCount;
};

} // namespace OpenMM

#endif /*OPENMM_CUDAPLATFORM_H_*/
