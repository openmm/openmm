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

/* ----------------------------------------------------------------------------- *
 *                                   AMD                                         *
 * ----------------------------------------------------------------------------- *
 * MIT License                                                                   *
 *                                                                               *
 * Copyright (c) 2020 Advanced Micro Devices, Inc.                               *
 *                                                                               *
 * Permission is hereby granted, free of charge, to any person obtaining a copy  *
 * of this software and associated documentation files (the "Software"), to deal *
 * in the Software without restriction, including without limitation the rights  *
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell     *
 * copies of the Software, and to permit persons to whom the Software is         *
 * furnished to do so, subject to the following conditions:                      *
 *                                                                               *
 * The above copyright notice and this permission notice shall be included in    *
 * all copies or substantial portions of the Software.                           *
 *                                                                               *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR    *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,      *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE   *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER        *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, *
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN     *
 * THE SOFTWARE.                                                                 *
 * ----------------------------------------------------------------------------- */

#include "HipContext.h"
#include "HipNonbondedUtilities.h"
#include "openmm/common/ComputeParameterSet.h"

namespace OpenMM {

class HipNonbondedUtilities;

/**
 * This class exists for backward compatibility.  For most purposes you can use
 * ComputeParameterSet directly instead.
 */

class OPENMM_EXPORT_COMMON HipParameterSet : public ComputeParameterSet {
public:
    /**
     * Create an HipParameterSet.
     *
     * @param context          the context for which to create the parameter set
     * @param numParameters    the number of parameters for each object
     * @param numObjects       the number of objects to store parameter values for
     * @param name             the name of the parameter set
     * @param bufferPerParameter  if true, a separate buffer is created for each parameter.  If false,
     *                            multiple parameters may be combined into a single buffer.
     * @param useDoublePrecision  whether values should be stored as single or double precision
     */
    HipParameterSet(HipContext& context, int numParameters, int numObjects, const std::string& name, bool bufferPerParameter=false, bool useDoublePrecision=false);
    /**
     * Get a set of HipNonbondedUtilities::ParameterInfo objects which describe the Buffers
     * containing the data.
     */
    std::vector<HipNonbondedUtilities::ParameterInfo>& getBuffers() {
        return buffers;
    }
private:
    std::vector<HipNonbondedUtilities::ParameterInfo> buffers;
};

} // namespace OpenMM

#endif /*OPENMM_CUDAPARAMETERSET_H_*/
