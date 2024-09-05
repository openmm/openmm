/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2019 Stanford University and the Authors.           *
 * Portions copyright (c) 2020-2023 Advanced Micro Devices, Inc.              *
 * Authors: Peter Eastman, Nicholas Curtis                                    *
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

#include "HipEvent.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

HipEvent::HipEvent(HipContext& context) : context(context), eventCreated(false) {
    hipError_t result = hipEventCreateWithFlags(&event, context.getEventFlags());
    if (result != hipSuccess)
        throw OpenMMException("Error creating HIP event:"+HipContext::getErrorString(result));
    eventCreated = true;
}

HipEvent::~HipEvent() {
    if (eventCreated)
        hipEventDestroy(event);
}

void HipEvent::enqueue() {
    hipEventRecord(event, context.getCurrentStream());
}

void HipEvent::wait() {
    hipEventSynchronize(event);
}
