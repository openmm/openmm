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
 * Portions copyright (c) 2009-2012 Stanford University and the Authors.      *
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

namespace OpenMM {

class CudaNonbondedUtilities;

/**
 * This class represents a set of floating point parameter values for a set of objects (particles, bonds, etc.).
 * It automatically creates an appropriate set of device buffers to hold the parameter values, based
 * on the number of parameters required.
 */

class OPENMM_EXPORT_CUDA CudaParameterSet {
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
    ~CudaParameterSet();
    /**
     * Get the number of parameters.
     */
    int getNumParameters() const {
        return numParameters;
    }
    /**
     * Get the number of objects.
     */
    int getNumObjects() const {
        return numObjects;
    }
    /**
     * Get the values of all parameters.
     *
     * @param values on exit, values[i][j] contains the value of parameter j for object i
     */
    template <class T>
    void getParameterValues(std::vector<std::vector<T> >& values);
    /**
     * Set the values of all parameters.
     *
     * @param values values[i][j] contains the value of parameter j for object i
     */
    template <class T>
    void setParameterValues(const std::vector<std::vector<T> >& values);
    /**
     * Get a set of CudaNonbondedUtilities::ParameterInfo objects which describe the Buffers
     * containing the data.
     */
    std::vector<CudaNonbondedUtilities::ParameterInfo>& getBuffers() {
        return buffers;
    }
    /**
     * Get a suffix to add to variable names when accessing a certain parameter.
     *
     * @param index         the index of the parameter
     * @param extraSuffix   an extra suffix to add to the variable name
     * @return the suffix to append
     */
    std::string getParameterSuffix(int index, const std::string& extraSuffix = "") const;
private:
    CudaContext& context;
    int numParameters, numObjects, elementSize;
    std::string name;
    std::vector<CudaNonbondedUtilities::ParameterInfo> buffers;
};

} // namespace OpenMM

#endif /*OPENMM_CUDAPARAMETERSET_H_*/
