/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 *                                                                            *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
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
#include "openmm/internal/AssertionUtilities.h"
#import <Metal/Metal.h>
#include <functional>
#include <iostream>
#include <thread>

using namespace OpenMM;
using namespace std;

/** Record a blit and notify the queue that this operation has ended. */
id<MTLCommandBuffer> fill(MetalQueue& queue, id<MTLBuffer> buffer, unsigned char value) {
    auto queueLock = queue.lock();
    id<MTLCommandBuffer> command = (__bridge id<MTLCommandBuffer>) queue.getCommandBuffer();
    id<MTLBlitCommandEncoder> encoder = [command blitCommandEncoder];
    ASSERT(encoder != nil);
    [encoder fillBuffer:buffer range:NSMakeRange(0, buffer.length) value:value];
    [encoder endEncoding];
    queue.submit((__bridge void*) command);
    return command;
}

void checkBytes(id<MTLBuffer> buffer, unsigned char expected) {
    const unsigned char* bytes = static_cast<const unsigned char*>(buffer.contents);
    for (NSUInteger i = 0; i < buffer.length; i++)
        ASSERT_EQUAL(expected, bytes[i]);
}

void expectException(const function<void()>& operation) {
    try {
        operation();
    }
    catch (const OpenMMException&) {
        return;
    }
    throw OpenMMException("Expected a Metal queue exception");
}

void testSubmissionPolicy(id<MTLDevice> device) {
    MetalQueue queue((__bridge void*) device);
    id<MTLBuffer> buffer = [device newBufferWithLength:256 options:MTLResourceStorageModeShared];
    id<MTLCommandBuffer> first = fill(queue, buffer, 7);
    id<MTLCommandBuffer> second = fill(queue, buffer, 9);
    ASSERT(first != second);
    ASSERT(first.status >= MTLCommandBufferStatusCommitted);
    ASSERT(second.status >= MTLCommandBufferStatusCommitted);
    queue.finish();
    checkBytes(buffer, 9);
    queue.finish(); // An empty queue is valid.
}

void testBlockingBoundary(id<MTLDevice> device) {
    MetalQueue queue((__bridge void*) device);
    id<MTLBuffer> buffer = [device newBufferWithLength:256 options:MTLResourceStorageModeShared];
    id<MTLCommandBuffer> command = fill(queue, buffer, 11);
    queue.wait((__bridge void*) command);
    checkBytes(buffer, 11);
    // A completed and already-reaped marker can still be waited for.
    queue.wait((__bridge void*) command);
}

void testEventOrdering(id<MTLDevice> device) {
    MetalQueue producer((__bridge void*) device), consumer((__bridge void*) device);
    id<MTLBuffer> source = [device newBufferWithLength:256 options:MTLResourceStorageModePrivate];
    id<MTLBuffer> result = [device newBufferWithLength:256 options:MTLResourceStorageModeShared];
    id<MTLEvent> event = [device newEvent];
    id<MTLCommandBuffer> recorded = fill(producer, source, 37);

    // External markers must follow the preceding operations before signaling.
    id<MTLCommandQueue> producerQueue = (__bridge id<MTLCommandQueue>) producer.getQueue();
    id<MTLCommandBuffer> marker = [producerQueue commandBuffer];
    [marker encodeSignalEvent:event value:1];
    producer.submit((__bridge void*) marker);
    ASSERT(recorded.status >= MTLCommandBufferStatusCommitted);
    ASSERT(marker.status >= MTLCommandBufferStatusCommitted);

    // An external wait must follow previously submitted consumer commands too.
    id<MTLCommandBuffer> prior = fill(consumer, result, 0);
    id<MTLCommandQueue> consumerQueue = (__bridge id<MTLCommandQueue>) consumer.getQueue();
    id<MTLCommandBuffer> wait = [consumerQueue commandBuffer];
    [wait encodeWaitForEvent:event value:1];
    consumer.submit((__bridge void*) wait);
    ASSERT(prior.status >= MTLCommandBufferStatusCommitted);
    id<MTLCommandBuffer> copy = (__bridge id<MTLCommandBuffer>) consumer.getCommandBuffer();
    id<MTLBlitCommandEncoder> encoder = [copy blitCommandEncoder];
    [encoder copyFromBuffer:source sourceOffset:0 toBuffer:result destinationOffset:0 size:256];
    [encoder endEncoding];
    consumer.submit((__bridge void*) copy);
    consumer.finish();
    producer.finish();
    checkBytes(result, 37);
}

void testValidation(id<MTLDevice> device) {
    MetalQueue queue((__bridge void*) device), other((__bridge void*) device);
    id<MTLCommandQueue> nativeQueue = (__bridge id<MTLCommandQueue>) queue.getQueue();
    id<MTLCommandBuffer> external = [nativeQueue commandBuffer];
    queue.submit((__bridge void*) external);
    id<MTLCommandBuffer> unsubmitted = [nativeQueue commandBuffer];
    id<MTLBuffer> buffer = [device newBufferWithLength:256 options:MTLResourceStorageModeShared];
    fill(queue, buffer, 19);
    expectException([&] { queue.submit(nullptr); });
    expectException([&] { queue.wait(nullptr); });
    expectException([&] { queue.wait((__bridge void*) unsubmitted); });
    expectException([&] { other.submit((__bridge void*) unsubmitted); });
    expectException([&] { queue.submit(other.getCommandBuffer()); });
    expectException([&] { queue.submit((__bridge void*) external); });
    expectException([&] { other.wait((__bridge void*) external); });
    queue.finish();
    checkBytes(buffer, 19);
}

void testDestructor(id<MTLDevice> device) {
    id<MTLBuffer> buffer = [device newBufferWithLength:256 options:MTLResourceStorageModeShared];
    {
        MetalQueue queue((__bridge void*) device);
        fill(queue, buffer, 23);
    }
    checkBytes(buffer, 23);
}

/** Common's callback worker can upload while the main thread records kernels. */
void testConcurrentRecording(id<MTLDevice> device) {
    MetalQueue queue((__bridge void*) device);
    id<MTLBuffer> first = [device newBufferWithLength:256 options:MTLResourceStorageModeShared];
    id<MTLBuffer> second = [device newBufferWithLength:256 options:MTLResourceStorageModeShared];
    thread worker([&] {
        @autoreleasepool {
            for (int i = 0; i < 130; i++)
                fill(queue, first, i);
        }
    });
    for (int i = 0; i < 130; i++)
        fill(queue, second, 255-i);
    worker.join();
    queue.finish();
    checkBytes(first, 129);
    checkBytes(second, 126);
}

int main() {
    @autoreleasepool {
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        if (device == nil) {
            cout << "No Metal device is available" << endl;
            return 77;
        }
        try {
            testSubmissionPolicy(device);
            testBlockingBoundary(device);
            testEventOrdering(device);
            testValidation(device);
            testDestructor(device);
            testConcurrentRecording(device);
        }
        catch (const exception& error) {
            cerr << error.what() << endl;
            return 1;
        }
        cout << "Metal queue tests passed" << endl;
    }
    return 0;
}
