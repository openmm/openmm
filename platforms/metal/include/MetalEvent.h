/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 *                                                                            *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM CUDA Platform.                                      *
 * Source: platforms/cuda/include/CudaEvent.h                                 *
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

#ifndef OPENMM_METALEVENT_H_
#define OPENMM_METALEVENT_H_

#include "openmm/common/ComputeEvent.h"
#include <memory>

namespace OpenMM {

class MetalContext;

/**
 * @brief Record a Metal queue position for host waits and cross-queue GPU dependencies.
 *
 * Adapted from CudaEvent. Each enqueue() uses a fresh native event so recording
 * again does not alter GPU waits already submitted for an earlier recording.
 * @warning The associated MetalContext must outlive this event. Host access must
 *          be serialized with recording and queue submission.
 */
class MetalEvent : public ComputeEventImpl {
public:
    /**
     * @brief Create an initially unrecorded event.
     * @param context Borrowed context whose current queue is used by enqueue().
     */
    explicit MetalEvent(MetalContext& context);
    /** @brief Release this event's native handles and recorded source-queue reference. */
    ~MetalEvent();
    /**
     * @brief Record a marker after preceding work on the context's current queue.
     * @throws OpenMMException If the current queue is invalid or recording/submission fails.
     * @note This does not wait for GPU completion. Subsequent wait() and queueWait()
     *       calls refer to this latest recording, even if the current queue changes.
     */
    void enqueue() override;
    /**
     * @brief Block the host until the recorded source queue completes through the marker.
     * @throws OpenMMException If the source queue reports an execution error.
     * @note An unrecorded event is a no-op; later source-queue commands are not waited for.
     */
    void wait() override;
    /**
     * @brief Submit a GPU dependency before subsequent work on a target queue.
     * @param queue Target MetalQueue on the same device as this event's context.
     * @throws OpenMMException If the queue is null, has a different backend/device,
     *         or wait-command creation/submission fails.
     * @note This does not block the host or change the context's selected queue.
     *       For an unrecorded event, the target is validated but no wait is submitted.
     */
    void queueWait(ComputeQueue queue) override;
private:
    class Impl;
    MetalContext& context;
    std::unique_ptr<Impl> impl;
};

} // namespace OpenMM

#endif /*OPENMM_METALEVENT_H_*/
