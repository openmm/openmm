/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 *                                                                            *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM CUDA Platform.                                      *
 * Source: platforms/cuda/include/CudaQueue.h                                 *
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

#ifndef OPENMM_METALQUEUE_H_
#define OPENMM_METALQUEUE_H_

#include "openmm/common/ComputeQueue.h"
#include <memory>
#include <mutex>

namespace OpenMM {

/**
 * @brief Metal command queue with submission tracking and GPU error reporting.
 *
 * Adapted from CudaQueue. Native handles are kept opaque so callers can include
 * this header from ordinary C++ code.
 * Recording operations hold lock() across encoding and submission. This permits
 * Common's worker-thread uploads to share a queue with the simulation thread.
 */
class MetalQueue : public ComputeQueueImpl {
public:
    /**
     * @brief Create and own a native command queue on the supplied device.
     * @param device Borrowed, non-null native @c id<MTLDevice> handle.
     * @throws OpenMMException If the device is null or queue creation fails.
     */
    explicit MetalQueue(void* device);
    /**
     * @brief Wait for tracked submissions before releasing native resources.
     * @note Destruction does not report GPU errors; use wait() or finish() explicitly.
     */
    ~MetalQueue();
    /**
     * @return A borrowed native @c id<MTLCommandQueue> handle; do not release it.
     * @note Submit its command buffers through submit() so completion and errors are tracked.
     */
    void* getQueue() const;
    /** @brief Serialize a complete encoding operation, including its submission.
     *  @note Recursive so an operation may call the synchronized queue methods.
     */
    std::unique_lock<std::recursive_mutex> lock();
    /**
     * @brief Get the command buffer for the next compute or blit operation.
     * @return A borrowed, uncommitted native @c id<MTLCommandBuffer> handle.
     * @note Hold lock(), end the encoder, and call submit() after each operation. With
     *       OPENMM_METAL_RECORD_AND_COMMIT=0, consecutive operations share a
     *       command buffer until a submission boundary or 64 operations.
     */
    void* getCommandBuffer();
    /**
     * @brief Finish recording an operation without waiting for GPU execution.
     * @param commandBuffer Borrowed native @c id<MTLCommandBuffer> from this queue,
     *                      not yet committed.
     * @throws OpenMMException If the buffer is null, foreign, or already committed,
     *         or an earlier completed submission reports an error.
     * @note OPENMM_METAL_RECORD_AND_COMMIT=1 commits each operation. When it is 0,
     *       buffers obtained from getCommandBuffer() are batched. Externally
     *       created buffers, including event markers and waits, always flush
     *       earlier recordings and commit immediately to preserve ordering.
     */
    void submit(void* commandBuffer);
    /**
     * @brief Commit recorded operations without waiting for GPU execution.
     * @throws OpenMMException If an earlier completed submission reports an error.
     * @note This is not ComputeContext::flushQueue(), which also waits for execution.
     */
    void flush();
    /**
     * @brief Submit recorded operations, wait for all submissions, and report errors.
     * @throws OpenMMException If GPU execution fails.
     * @note An empty queue requires no wait. Other queues are not synchronized.
     */
    void finish();
    /**
     * @brief Wait through a submitted marker, checking preceding tracked commands for errors.
     * @param commandBuffer Borrowed native @c id<MTLCommandBuffer> committed on
     *                      this queue, or the current recording buffer from
     *                      getCommandBuffer(), which is submitted before waiting.
     * @throws OpenMMException If the marker is null, foreign, or an uncommitted
     *         external buffer, or GPU execution fails. Tracked commands through
     *         the marker are drained before their first execution error is reported.
     * @note Later submissions are not waited for. Errors from already-reaped preceding
     *       commands are not retained; the marker's own status is always checked.
     *       Waiting for a recording buffer submits and waits for its entire batch.
     */
    void wait(void* commandBuffer);
private:
    /** Commit one buffer and track it until execution completes. */
    void commit(void* commandBuffer);
    class Impl;
    std::unique_ptr<Impl> impl;
};

} // namespace OpenMM

#endif /*OPENMM_METALQUEUE_H_*/
