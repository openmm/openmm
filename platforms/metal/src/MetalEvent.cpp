/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2019-2026 Stanford University and the Authors.      *
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

#include "MetalEvent.h"
#include "MetalQueue.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

MetalEvent::MetalEvent(MetalContext& context) : context(context), event(nullptr), currentBuffer(nullptr), value(0) {
    event = context.getDevice().newSharedEvent();
    if (event == nullptr)
        throw OpenMMException("Error creating Metal event");
}

MetalEvent::~MetalEvent() {
    if (event != nullptr)
        event->release();
    if (currentBuffer != nullptr)
        currentBuffer->release();
}

void MetalEvent::enqueue() {
    MetalQueue& queue = dynamic_cast<MetalQueue&>(*context.getCurrentQueue());
    queue.flush();
    currentBuffer = &queue.getCommandBuffer();
    currentBuffer->retain();
    currentBuffer->encodeSignalEvent(event, ++value);
    queue.flush();
}

void MetalEvent::wait() {
    currentBuffer->waitUntilCompleted();
    currentBuffer->release();
    currentBuffer = nullptr;
}

void MetalEvent::queueWait(ComputeQueue queue) {
    MetalQueue& metalQueue = dynamic_cast<MetalQueue&>(*queue);
    metalQueue.flush();
    metalQueue.getCommandBuffer().encodeWait(event, value);
    currentBuffer->release();
    currentBuffer = nullptr;
}
