#ifndef OPENMM_CUSTOMVOLUMEFORCE_H_
#define OPENMM_CUSTOMVOLUMEFORCE_H_

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
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "Force.h"
#include "internal/windowsExport.h"
#include <vector>

namespace OpenMM {

/**
 * This class computes an energy that depends only on the volume of the periodic box, or more generally
 * on the box shape as specified by the elements of the box vectors.  Because the energy does not
 * depend on particle positions, it does not apply any forces to particles.  It is primarily useful
 * for constant pressure simulations, where the volume-dependent energy can influence the behavior
 * of the barostat.  Energy terms of this sort are often used for pressure matching in coarse grained
 * force fields.
 *
 * To use this class, create a CustomVolumeForce object, passing an algebraic expression to the constructor
 * that defines the energy.  The expression may depend on the following variables.
 *
 * <ul>
 * <li>v: The volume of the periodic box in nm^3.</li>
 * <li>ax: The x component of the first box vector in nm.  (The y and z components are always zero.)</li>
 * <li>bx, by: The x and y components of the second box vector in nm.  (The z component is always zero.)</li>
 * <li>cx, cy, cz: The x, y and z components of the third box vector in nm.</li>
 * <li>Global parameters that you define by calling addGlobalParameter().</li>
 * </ul>
 * 
 * The initial value of a global parameter is specified in the call to addGlobalParameter().  Theire values
 * can be modified during a simulation by calling Context::setParameter().
 *
 * Expressions may involve the operators + (add), - (subtract), * (multiply), / (divide), and ^ (power), and the following
 * functions: sqrt, exp, log, sin, cos, sec, csc, tan, cot, asin, acos, atan, atan2, sinh, cosh, tanh, erf, erfc, min, max, abs, floor, ceil, step, delta, select.  All trigonometric functions
 * are defined in radians, and log is the natural logarithm.  step(x) = 0 if x is less than 0, 1 otherwise.  delta(x) = 1 if x is 0, 0 otherwise.
 * select(x,y,z) = z if x = 0, y otherwise.
 */

class OPENMM_EXPORT CustomVolumeForce : public Force {
public:
    /**
     * Create a CustomVolumeForce.
     *
     * @param energy    an algebraic expression giving the energy as a function of the box shape
     */
    explicit CustomVolumeForce(const std::string& energy);
    /**
     * Get the number of global parameters that the energy depends on.
     */
    int getNumGlobalParameters() const {
        return globalParameters.size();
    }
    /**
     * Get the algebraic expression that defines the energy.
     */
    const std::string& getEnergyFunction() const;
    /**
     * Set the algebraic expression that defines the energy.
     */
    void setEnergyFunction(const std::string& energy);
    /**
     * Add a new global parameter that the interaction may depend on.  The default value provided to
     * this method is the initial value of the parameter in newly created Contexts.  You can change
     * the value at any time by calling setParameter() on the Context.
     *
     * @param name             the name of the parameter
     * @param defaultValue     the default value of the parameter
     * @return the index of the parameter that was added
     */
    int addGlobalParameter(const std::string& name, double defaultValue);
    /**
     * Get the name of a global parameter.
     *
     * @param index     the index of the parameter for which to get the name
     * @return the parameter name
     */
    const std::string& getGlobalParameterName(int index) const;
    /**
     * Set the name of a global parameter.
     *
     * @param index          the index of the parameter for which to set the name
     * @param name           the name of the parameter
     */
    void setGlobalParameterName(int index, const std::string& name);
    /**
     * Get the default value of a global parameter.
     *
     * @param index     the index of the parameter for which to get the default value
     * @return the parameter default value
     */
    double getGlobalParameterDefaultValue(int index) const;
    /**
     * Set the default value of a global parameter.
     *
     * @param index          the index of the parameter for which to set the default value
     * @param defaultValue   the default value of the parameter
     */
    void setGlobalParameterDefaultValue(int index, double defaultValue);
    /**
     * Returns whether or not this force makes use of periodic boundary conditions.  Because this
     * class is only applicable to periodic systems, this always returns true.
     */
    bool usesPeriodicBoundaryConditions() const {
        return true;
    }
protected:
    ForceImpl* createImpl() const;
private:
    class GlobalParameterInfo;
    std::string energyExpression;
    std::vector<GlobalParameterInfo> globalParameters;
};

/**
 * This is an internal class used to record information about a global parameter.
 * @private
 */
class CustomVolumeForce::GlobalParameterInfo {
public:
    std::string name;
    double defaultValue;
    GlobalParameterInfo() {
    }
    GlobalParameterInfo(const std::string& name, double defaultValue) : name(name), defaultValue(defaultValue) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_CUSTOMVOLUMEFORCE_H_*/
