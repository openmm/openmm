#ifndef OPENMM_METALQUEUE_H_
#define OPENMM_METALQUEUE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2026 Stanford University and the Authors.           *
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
#include "Metal.hpp"

namespace OpenMM {

/**
 * This is the Metal implementation of the ComputeQueue interface.  It wraps a Metal CommandQueue.
 */

class MetalQueue : public ComputeQueueImpl {
public:
    /**
     * Create a MetalQueue.
     */
    MetalQueue(MTL::Device& device);
    ~MetalQueue();
    /**
     * Get an encoder that can be used for launching kernels on this queue.
     */
    MTL::ComputeCommandEncoder& getEncoder();
    /**
     * Get a command buffer that can be used for launching kernels on this queue.
     */
    MTL::CommandBuffer& getCommandBuffer();
    /**
     * Flush the queue, ensuring that all work that has been queued has been submitted to the device, and optionally
     * unit it has completed.
     *
     * @param sync   If true, wait until all work has completed.  If false (the default), ensure the work has been
     *               submitted but do not wait for it to complete.
     */
    void flush(bool sync=false);
private:
    void ensureEncoderExists();
    MTL::CommandQueue* queue;
    MTL::CommandBuffer* commandBuffer;
    MTL::ComputeCommandEncoder* encoder;
};

} // namespace OpenMM

#endif /*OPENMM_METALQUEUE_H_*/
