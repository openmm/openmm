/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 *                                                                            *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM CUDA Platform.                                      *
 * Source: platforms/cuda/src/CudaEvent.cpp                                   *
 *                                                                            *
 * Original CUDA Platform code:                                               *
 * Portions copyright (c) 2019-2025 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 *                                                                            *
 * Metal Platform code:                                                       *
 * Portions copyright (c) 2026 Chun-Chi Hung.                                 *
 * Authors: Chun-Chi Hung                                                     *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the               *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.       *
 * -------------------------------------------------------------------------- */

#include "MetalEvent.h"
#include "MetalContext.h"
#include "MetalQueue.h"
#include "openmm/OpenMMException.h"
#import <Metal/Metal.h>

using namespace OpenMM;

class MetalEvent::Impl {
public:
    id<MTLEvent> event;
    id<MTLCommandBuffer> marker;
    ComputeQueue queue;
};

MetalEvent::MetalEvent(MetalContext& context) : context(context), impl(new Impl()) {
}

MetalEvent::~MetalEvent() {
}

void MetalEvent::enqueue() {
    MetalQueue& queue = context.getCurrentMetalQueue();
    id<MTLCommandQueue> commandQueue = (__bridge id<MTLCommandQueue>) queue.getQueue();
    // Separate recordings must not signal each other when made on different queues.
    id<MTLEvent> event = [commandQueue.device newEvent];
    id<MTLCommandBuffer> marker = [commandQueue commandBuffer];
    if (event == nil || marker == nil)
        throw OpenMMException("Error recording Metal event");
    [marker encodeSignalEvent:event value:1];
    queue.submit((__bridge void*) marker);
    impl->event = event;
    impl->marker = marker;
    impl->queue = context.getCurrentQueue();
}

void MetalEvent::wait() {
    if (impl->marker != nil)
        static_cast<MetalQueue&>(*impl->queue).wait((__bridge void*) impl->marker);
}

void MetalEvent::queueWait(ComputeQueue queue) {
    MetalQueue* target = dynamic_cast<MetalQueue*>(queue.get());
    if (target == nullptr)
        throw OpenMMException("Metal event requires a Metal command queue");
    id<MTLCommandQueue> commandQueue = (__bridge id<MTLCommandQueue>) target->getQueue();
    if (commandQueue.device != (__bridge id<MTLDevice>) context.getDevice())
        throw OpenMMException("Metal event and command queue belong to different devices");
    if (impl->marker == nil)
        return;
    id<MTLCommandBuffer> buffer = [commandQueue commandBuffer];
    if (buffer == nil)
        throw OpenMMException("Error creating Metal event wait command buffer");
    [buffer encodeWaitForEvent:impl->event value:1];
    target->submit((__bridge void*) buffer);
}
