#ifndef OPENMM_COMPUTEPARAMETERSET_H_
#define OPENMM_COMPUTEPARAMETERSET_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2024 Stanford University and the Authors.      *
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

#include "openmm/common/ArrayInterface.h"
#include "openmm/common/ComputeContext.h"
#include "openmm/common/ComputeParameterInfo.h"
#include <string>
#include <vector>

namespace OpenMM {

/**
 * This class represents a set of floating point parameter values for a set of objects (particles, bonds, etc.).
 * It automatically creates an appropriate set of arrays to hold the parameter values, based
 * on the number of parameters required.
 */

class OPENMM_EXPORT_COMMON ComputeParameterSet {
public:
    /**
     * Create an ComputeParameterSet.
     *
     * @param context          the context for which to create the parameter set
     * @param numParameters    the number of parameters for each object
     * @param numObjects       the number of objects to store parameter values for
     * @param name             the name of the parameter set
     * @param arrayPerParameter   if true, a separate array is created for each parameter.  If false,
     *                            multiple parameters may be combined into a single array for efficiency.
     * @param useDoublePrecision  whether values should be stored as single or double precision
     */
    ComputeParameterSet(ComputeContext& context, int numParameters, int numObjects, const std::string& name, bool arrayPerParameter=false, bool useDoublePrecision=false);
    ~ComputeParameterSet();
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
     * @param values    values[i][j] contains the value of parameter j for object i
     * @param convert   if true, automatic conversions between single and double
     *                  precision will be performed as necessary
     */
    template <class T>
    void setParameterValues(const std::vector<std::vector<T> >& values, bool convert=false);
    /**
     * Set the values of all parameters for a subset of objects.  They must be
     * contained in a continuous range.
     *
     * @param first     the index of the first object for which to set parameter values
     * @param values    values[i][j] contains the value of parameter j for object i+first
     * @param convert   if true, automatic conversions between single and double
     *                  precision will be performed as necessary
     */
    template <class T>
    void setParameterValuesSubset(int first, const std::vector<std::vector<T> >& values, bool convert=false);
    /**
     * Get a vector of ComputeParameterInfo objects which describe the arrays
     * containing the data.
     */
    std::vector<ComputeParameterInfo>& getParameterInfos() {
        return parameters;
    }
    /**
     * Get a suffix to add to variable names when accessing a certain parameter.
     *
     * @param index         the index of the parameter
     * @param extraSuffix   an extra suffix to add to the variable name
     * @return the suffix to append
     */
    std::string getParameterSuffix(int index, const std::string& extraSuffix="") const;
private:
    ComputeContext& context;
    int numParameters, numObjects, elementSize;
    std::string name;
    std::vector<ArrayInterface*> arrays;
    std::vector<ComputeParameterInfo> parameters;
};

} // namespace OpenMM

#endif /*OPENMM_COMPUTEPARAMETERSET_H_*/
