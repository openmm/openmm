#ifndef OPENMM_CUDAEVENT_H_
#define OPENMM_CUDAEVENT_H_

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

#include "CudaContext.h"
#include "openmm/common/ComputeEvent.h"

namespace OpenMM {

/**
 * This is the CUDA implementation of the ComputeKernelImpl interface. 
 */

class CudaEvent : public ComputeEventImpl {
public:
    CudaEvent(CudaContext& context);
    ~CudaEvent();
    /**
     * Place the event into the device's execution queue.
     */
    void enqueue();
    /**
     * Block until all operations started before the call to enqueue() have completed.
     */
    void wait();
private:
    CudaContext& context;
    CUevent event;
    bool eventCreated;
};

} // namespace OpenMM

#endif /*OPENMM_CUDAEVENT_H_*/
