#ifndef OPENMM_COMPUTEPROGRAM_H_
#define OPENMM_COMPUTEPROGRAM_H_

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

#include "openmm/common/ComputeKernel.h"
#include <memory>

namespace OpenMM {

/**
 * This abstract class represents a compiled program that can be executed on a computing
 * device.  A ComputeProgramImpl is created by calling compileProgram() on a ComputeContext,
 * which returns an instance of a platform-specific subclass.  The source code for a
 * ComputeProgramImpl typically contains one or more kernels.  Call createKernel() to get
 * ComputeKernels for the kernels, which can then be executed.
 * 
 * Instead of referring to this class directly, it is best to use ComputeProgram, which is
 * a typedef for a shared_ptr to a ComputeProgramImpl.  This allows you to treat it as having
 * value semantics, and frees you from having to manage memory.  
 */

class OPENMM_EXPORT_COMMON ComputeProgramImpl {
public:
    virtual ~ComputeProgramImpl() {
    }
    /**
     * Create a ComputeKernel for one of the kernels in this program.
     * 
     * @param name    the name of the kernel to get
     */
    virtual ComputeKernel createKernel(const std::string& name) = 0;
};

typedef std::shared_ptr<ComputeProgramImpl> ComputeProgram;

} // namespace OpenMM

#endif /*OPENMM_COMPUTEPROGRAM_H_*/
