#ifndef OPENMM_METALEVENT_H_
#define OPENMM_METALEVENT_H_

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

#include "MetalContext.h"
#include "openmm/common/ComputeEvent.h"

namespace OpenMM {

/**
 * This is the Metal implementation of the ComputeEventImpl interface.
 */

class MetalEvent : public ComputeEventImpl {
public:
    MetalEvent(MetalContext& context);
    ~MetalEvent();
    /**
     * Place the event into the device's execution queue.
     */
    void enqueue();
    /**
     * Block until all operations started before the call to enqueue() have completed.
     */
    void wait();
    /**
     * Enqueue a barrier that causes a specified ComputeQueue to block until all
     * operations started before the call to enqueue() have completed.
     */
    void queueWait(ComputeQueue queue);
private:
    MetalContext& context;
    MTL::SharedEvent* event;
    MTL::CommandBuffer* currentBuffer;
    unsigned long value;
};

} // namespace OpenMM

#endif /*OPENMM_METALEVENT_H_*/
