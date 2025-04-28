#ifndef OPENMM_COMPUTEQUEUE_H_
#define OPENMM_COMPUTEQUEUE_H_

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

#include "openmm/common/windowsExportCommon.h"
#include <memory>

namespace OpenMM {

/**
 * This abstract class represents a queue within which kernels can be executed.  Call
 * createQueue() on a ComputeContext to create an instance of a platform-specific
 * subclass.  You can then pass it to the ComputeContext's setQueue() method to cause
 * kernels to be launched on it.
 * 
 * Instead of referring to this class directly, it is best to use ComputeQueue, which is
 * a typedef for a shared_ptr to a ComputeQueueImpl.  This allows you to treat it as having
 * value semantics, and frees you from having to manage memory.  
 */

class OPENMM_EXPORT_COMMON ComputeQueueImpl {
public:
    virtual ~ComputeQueueImpl() {
    }
};

typedef std::shared_ptr<ComputeQueueImpl> ComputeQueue;

} // namespace OpenMM

#endif /*OPENMM_COMPUTEQUEUE_H_*/
