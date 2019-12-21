#ifndef OPENMM_CUDAEXPRESSIONUTILITIES_H_
#define OPENMM_CUDAEXPRESSIONUTILITIES_H_

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

#include "openmm/common/ExpressionUtilities.h"
#include "openmm/common/windowsExportCommon.h"

namespace OpenMM {

/**
 * This class exists only for backward compatibility.  It adds no features beyond
 * the base ExpressionUtilities class.
 */

class OPENMM_EXPORT_COMMON CudaExpressionUtilities : public ExpressionUtilities {
public:
    CudaExpressionUtilities(ComputeContext& context) : ExpressionUtilities(context) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_CUDAEXPRESSIONUTILITIES_H_*/
