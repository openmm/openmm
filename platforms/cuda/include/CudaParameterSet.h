#ifndef OPENMM_CUDAPARAMETERSET_H_
#define OPENMM_CUDAPARAMETERSET_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2019 Stanford University and the Authors.      *
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
#include "CudaNonbondedUtilities.h"
#include "openmm/common/ComputeParameterSet.h"

namespace OpenMM {

class CudaNonbondedUtilities;

/**
 * This class exists for backward compatibility.  For most purposes you can use
 * ComputeParameterSet directly instead.
 */

class OPENMM_EXPORT_COMMON CudaParameterSet : public ComputeParameterSet {
public:
    /**
     * Create an CudaParameterSet.
     *
     * @param context          the context for which to create the parameter set
     * @param numParameters    the number of parameters for each object
     * @param numObjects       the number of objects to store parameter values for
     * @param name             the name of the parameter set
     * @param bufferPerParameter  if true, a separate buffer is created for each parameter.  If false,
     *                            multiple parameters may be combined into a single buffer.
     * @param useDoublePrecision  whether values should be stored as single or double precision
     */
    CudaParameterSet(CudaContext& context, int numParameters, int numObjects, const std::string& name, bool bufferPerParameter=false, bool useDoublePrecision=false);
    /**
     * Get a set of CudaNonbondedUtilities::ParameterInfo objects which describe the Buffers
     * containing the data.
     */
    std::vector<CudaNonbondedUtilities::ParameterInfo>& getBuffers() {
        return buffers;
    }
private:
    std::vector<CudaNonbondedUtilities::ParameterInfo> buffers;
};

} // namespace OpenMM

#endif /*OPENMM_CUDAPARAMETERSET_H_*/
