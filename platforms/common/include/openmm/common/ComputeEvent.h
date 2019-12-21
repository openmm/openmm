#ifndef OPENMM_COMPUTEEVENT_H_
#define OPENMM_COMPUTEEVENT_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2019 Stanford University and the Authors.           *
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

#include <memory>

namespace OpenMM {

/**
 * This abstract class represents an event for synchronization between the host and
 * device.  It is created by calling createEvent() on a ComputeContext, which returns
 * an instance of a platform-specific subclass.  To use it, call enqueue() immediately
 * after starting an asynchronous operation, such as a kernel invocation or non-blocking
 * data transfer.  Then at a later point call wait().  This will cause the host to block
 * until all operations started before the call to enequeue() have completed.
 * 
 * Instead of referring to this class directly, it is best to use a ComputeEvent, which is
 * a typedef for a shared_ptr to a ComputeEventImpl.  This allows you to treat it as having
 * value semantics, and frees you from having to manage memory.  
 */

class OPENMM_EXPORT_COMMON ComputeEventImpl {
public:
    virtual ~ComputeEventImpl() {
    }
    /**
     * Place the event into the device's execution queue.
     */
    virtual void enqueue() = 0;
    /**
     * Block until all operations started before the call to enqueue() have completed.
     */
    virtual void wait() = 0;
};

typedef std::shared_ptr<ComputeEventImpl> ComputeEvent;

} // namespace OpenMM

#endif /*OPENMM_COMPUTEEVENT_H_*/
