/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
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

#include "HipQueue.h"
#include "HipContext.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

HipQueue::HipQueue(hipStream_t stream) : stream(stream), initialized(false) {
}

HipQueue::HipQueue() : initialized(false) {
    hipError_t result = hipStreamCreateWithFlags(&stream, hipStreamNonBlocking);
    if (result != hipSuccess)
        throw OpenMMException("Error creating HIP stream: "+HipContext::getErrorString(result));
    initialized = true;
}

HipQueue::~HipQueue() {
    if (initialized)
        hipStreamDestroy(stream);
}
