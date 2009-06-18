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
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
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
            hasPeriodicTorsions(false), hasRB(false), hasNonbonded(false), primaryKernel(NULL), stepCount(0), time(0.0) {
    }
    _gpuContext* gpu;
    KernelImpl* primaryKernel;
    bool removeCM;
    bool hasBonds, hasAngles, hasPeriodicTorsions, hasRB, hasNonbonded;
    int nonbondedMethod;
    int cmMotionFrequency;
    int stepCount;
    double time;
};

} // namespace OpenMM

#endif /*OPENMM_CUDAPLATFORM_H_*/
