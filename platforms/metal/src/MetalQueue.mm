/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 *                                                                            *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM CUDA Platform.                                      *
 * Source: platforms/cuda/src/CudaQueue.cpp                                   *
 *                                                                            *
 * Original CUDA Platform code:                                               *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
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

#include "MetalQueue.h"
#include "openmm/OpenMMException.h"
#import <Metal/Metal.h>
#include <algorithm>
#include <deque>
#include <string>

using namespace OpenMM;
using namespace std;

class MetalQueue::Impl {
public:
    id<MTLCommandQueue> queue;
    deque<id<MTLCommandBuffer>> pending;
};

MetalQueue::MetalQueue(void* device) : impl(new Impl()) {
    if (device == nullptr)
        throw OpenMMException("Cannot create a Metal queue without a device");
    impl->queue = [(__bridge id<MTLDevice>) device newCommandQueue];
    if (impl->queue == nil)
        throw OpenMMException("Error creating Metal command queue");
}

MetalQueue::~MetalQueue() {
    // Finish outstanding device work; explicit waits report execution errors.
    for (id<MTLCommandBuffer> buffer : impl->pending)
        [buffer waitUntilCompleted];
}

void* MetalQueue::getQueue() const {
    return (__bridge void*) impl->queue;
}

void MetalQueue::submit(void* commandBuffer) {
    id<MTLCommandBuffer> buffer = (__bridge id<MTLCommandBuffer>) commandBuffer;
    if (buffer == nil || buffer.commandQueue != impl->queue)
        throw OpenMMException("Metal command buffer does not belong to this queue");
    if (buffer.status >= MTLCommandBufferStatusCommitted)
        throw OpenMMException("Metal command buffer has already been submitted");
    while (!impl->pending.empty() && impl->pending.front().status >= MTLCommandBufferStatusCompleted)
        wait((__bridge void*) impl->pending.front());
    impl->pending.push_back(buffer);
    [buffer commit];
}

void MetalQueue::finish() {
    if (!impl->pending.empty())
        wait((__bridge void*) impl->pending.back());
}

void MetalQueue::wait(void* commandBuffer) {
    id<MTLCommandBuffer> marker = (__bridge id<MTLCommandBuffer>) commandBuffer;
    if (marker == nil || marker.commandQueue != impl->queue || marker.status < MTLCommandBufferStatusCommitted)
        throw OpenMMException("Cannot wait for an unsubmitted or foreign Metal command buffer");
    auto markerPosition = find(impl->pending.begin(), impl->pending.end(), marker);
    string error;
    // A marker may already have been reaped by a later submission.
    size_t count = (markerPosition == impl->pending.end() ? 0 : markerPosition-impl->pending.begin()+1);
    for (size_t i = 0; i < count; i++) {
        id<MTLCommandBuffer> buffer = impl->pending.front();
        // Check each command so a successful marker cannot hide an earlier error.
        [buffer waitUntilCompleted];
        if (buffer.status == MTLCommandBufferStatusError && error.empty()) {
            const char* message = buffer.error.localizedDescription.UTF8String;
            error = (message == nullptr ? "Unknown GPU error" : message);
        }
        impl->pending.pop_front();
    }
    [marker waitUntilCompleted];
    if (marker.status == MTLCommandBufferStatusError && error.empty()) {
        const char* message = marker.error.localizedDescription.UTF8String;
        error = (message == nullptr ? "Unknown GPU error" : message);
    }
    if (!error.empty())
        throw OpenMMException("Error executing Metal command buffer: "+error);
}
