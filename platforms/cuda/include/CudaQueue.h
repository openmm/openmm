#ifndef OPENMM_CUDAQUEUE_H_
#define OPENMM_CUDAQUEUE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
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

#include "openmm/common/ComputeQueue.h"
#include <cuda.h>

namespace OpenMM {

/**
 * This is the CUDA implementation of the ComputeQueue interface.  It wraps a CUstream.
 */

class CudaQueue : public ComputeQueueImpl {
public:
    /**
     * Create a CudaQueue that wraps an existing CUstream.
     */
    CudaQueue(CUstream stream);
    /**
     * Create a CudaQueue that create a new CUstream.
     */
    CudaQueue();
    ~CudaQueue();
    /**
     * Get the CUstream.
     */
    CUstream getStream() {
        return stream;
    }
private:
    CUstream stream;
    bool initialized;
};

} // namespace OpenMM

#endif /*OPENMM_CUDAQUEUE_H_*/
