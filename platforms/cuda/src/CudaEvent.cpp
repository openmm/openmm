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

#include "CudaEvent.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

CudaEvent::CudaEvent(CudaContext& context) : context(context), eventCreated(false) {
    CUresult result = cuEventCreate(&event, CU_EVENT_DISABLE_TIMING);
    if (result != CUDA_SUCCESS)
        throw OpenMMException("Error creating CUDA event:"+CudaContext::getErrorString(result));
    eventCreated = true;
}

CudaEvent::~CudaEvent() {
    if (eventCreated)
        cuEventDestroy(event);
}

void CudaEvent::enqueue() {
    cuEventRecord(event, 0);
}

void CudaEvent::wait() {
    cuEventSynchronize(event);
}
