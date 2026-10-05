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

#include "MetalQueue.h"
#include "MetalContext.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

MetalQueue::MetalQueue(MTL::Device& device) : queue(nullptr), commandBuffer(nullptr), encoder(nullptr) {
    queue = device.newCommandQueue();
}

MetalQueue::~MetalQueue() {
    flush();
    if (queue != nullptr)
        queue->release();
}

MTL::CommandQueue& MetalQueue::getQueue() {
    return *queue;
}

MTL::ComputeCommandEncoder& MetalQueue::getEncoder() {
    if (encoder == nullptr) {
        NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();
        if (commandBuffer == nullptr)
            commandBuffer = queue->commandBuffer()->retain();
        encoder = commandBuffer->computeCommandEncoder()->retain();
        pool->release();
    }
    return *encoder;
}

MTL::CommandBuffer& MetalQueue::getCommandBuffer() {
    if (commandBuffer == nullptr) {
        NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();
        commandBuffer = queue->commandBuffer()->retain();
        pool->release();
    }
    return *commandBuffer;
}

void MetalQueue::flush(bool sync) {
    if (encoder != nullptr) {
        encoder->endEncoding();
        commandBuffer->commit();
        if (sync)
            commandBuffer->waitUntilCompleted();
        encoder->release();
        commandBuffer->release();
        encoder = nullptr;
        commandBuffer = nullptr;
    }
}
