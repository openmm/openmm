#ifndef OPENMM_COMPUTEPARAMETERINFO_H_
#define OPENMM_COMPUTEPARAMETERINFO_H_

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

#include "openmm/common/ArrayInterface.h"
#include <sstream>
#include <string>

namespace OpenMM {

/**
 * This class stores information about a parameter that can be passed to a kernel.
 * It combines an ArrayInterface holding parameter values with additional information
 * describing how to represent it in kernels: the variable name, the data type, etc.
 * 
 * The array is assumed to contain a parameter value for each of many objects (atoms,
 * bonds, etc.).  Each value may in turn be a multi-component vector.  When creating
 * a ComputeParameterInfo, specify the number of components in the vector and the
 * type of each component.  For example, suppose you have an array of type float3
 * containing a dipole moment for each atom.  The ComputeParameterInfo would be
 * created like this:
 * 
 * ComputeParameterInfo parameter(dipoleArray, "dipole", "float", 3);
 */

class ComputeParameterInfo {
public:
    /**
     * Create a ComputeParameterInfo.
     *
     * @param array          the array containing the parameter values
     * @param name           the name of the variable to use for this parameter
     * @param type           the data type of the parameter's components
     * @param numComponents  the number of components in the parameter
     * @param constant       whether the array memory should be marked as constant
     */
    ComputeParameterInfo(ArrayInterface& array, const std::string& name, const std::string& componentType, int numComponents, bool constant=true) :
            array(array), name(name), componentType(componentType), numComponents(numComponents), constant(constant) {
        if (numComponents == 1)
            type = componentType;
        else {
            std::stringstream s;
            s << componentType << numComponents;
            type = s.str();
        }
    }
    virtual ~ComputeParameterInfo() {
    }
    /**
     * Get the array containing the parameter values.
     */
    ArrayInterface& getArray() {
        return array;
    }
    /**
     * Get the array containing the parameter values.
     */
    const ArrayInterface& getArray() const {
        return array;
    }
    /**
     * Get the name of the variable to use for this parameter.
     */
    const std::string& getName() const {
        return name;
    }
    /**
     * Get the data type of each component of the value.  For example, if getType() returns "float3",
     * this will return "float".
     */
    const std::string& getComponentType() const {
        return componentType;
    }
    /**
     * Get the data type of each value.
     */
    const std::string& getType() const {
        return type;
    }
    /**
     * Get the number of components in each value.  If the values are not a vector
     * type, this returns 1.
     */
    int getNumComponents() const {
        return numComponents;
    }
    /**
     * Get the size of each parameter value in bytes.
     */
    int getSize() const {
        return array.getElementSize();
    }
    /**
     * Get whether the array memory should be marked as constant.
     */
    bool isConstant() const {
        return constant;
    }
private:
    ArrayInterface& array;
    std::string name;
    std::string componentType;
    std::string type;
    int numComponents;
    bool constant;
};

} // namespace OpenMM

#endif /*OPENMM_COMPUTEPARAMETERINFO_H_*/
