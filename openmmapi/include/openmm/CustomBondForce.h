#ifndef OPENMM_CUSTOMBONDEDFORCE_H_
#define OPENMM_CUSTOMBONDEDFORCE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
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
#include "Vec3.h"
#include <vector>
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This class implements bonded interactions between pairs of particles.  Unlike HarmonicBondForce, the functional form
 * of the interaction is completely customizable, and may involve arbitrary algebraic expressions.
 * It may depend on the distance between particles, as well as on arbitrary global and
 * per-bond parameters.
 *
 * To use this class, create a CustomBondForce object, passing an algebraic expression to the constructor
 * that defines the interaction energy between each pair of bonded particles.  The expression may depend on r, the distance
 * between the particles, as well as on any parameters you choose.  Then call addPerBondParameter() to define per-bond
 * parameters, and addGlobalParameter() to define global parameters.  The values of per-bond parameters are specified as
 * part of the system definition, while values of global parameters may be modified during a simulation by calling Context::setParameter().
 * Finally, call addBond() once for each bond.  After a bond has been added, you can modify its parameters by calling setBondParameters().
 * This will have no effect on Contexts that already exist unless you call updateParametersInContext().
 *
 * As an example, the following code creates a CustomBondForce that implements a harmonic potential:
 *
 * \verbatim embed:rst:leading-asterisk
 * .. code-block:: cpp
 *
 *    CustomBondForce* force = new CustomBondForce("0.5*k*(r-r0)^2");
 *
 * \endverbatim
 *
 * This force depends on two parameters: the spring constant k and equilibrium distance r0.  The following code defines these parameters:
 *
 * \verbatim embed:rst:leading-asterisk
 * .. code-block:: cpp
 *
 *    force->addPerBondParameter("k");
 *    force->addPerBondParameter("r0");
 *
 * \endverbatim
 * 
 * This class also has the ability to compute derivatives of the potential energy with respect to global parameters.
 * Call addEnergyParameterDerivative() to request that the derivative with respect to a particular parameter be
 * computed.  You can then query its value in a Context by calling getState() on it.
 *
 * Expressions may involve the operators + (add), - (subtract), * (multiply), / (divide), and ^ (power), and the following
 * functions: sqrt, exp, log, sin, cos, sec, csc, tan, cot, asin, acos, atan, atan2, sinh, cosh, tanh, erf, erfc, min, max, abs, floor, ceil, step, delta, select.  All trigonometric functions
 * are defined in radians, and log is the natural logarithm.  step(x) = 0 if x is less than 0, 1 otherwise.  delta(x) = 1 if x is 0, 0 otherwise.
 * select(x,y,z) = z if x = 0, y otherwise.
 */

class OPENMM_EXPORT CustomBondForce : public Force {
public:
    /**
     * Create a CustomBondForce.
     *
     * @param energy    an algebraic expression giving the interaction energy between two bonded particles as a function
     *                  of r, the distance between them
     */
    explicit CustomBondForce(const std::string& energy);
    /**
     * Get the number of bonds for which force field parameters have been defined.
     */
    int getNumBonds() const {
        return bonds.size();
    }
    /**
     * Get the number of per-bond parameters that the interaction depends on.
     */
    int getNumPerBondParameters() const {
        return parameters.size();
    }
    /**
     * Get the number of global parameters that the interaction depends on.
     */
    int getNumGlobalParameters() const {
        return globalParameters.size();
    }
    /**
     * Get the number of global parameters with respect to which the derivative of the energy
     * should be computed.
     */
    int getNumEnergyParameterDerivatives() const {
        return energyParameterDerivatives.size();
    }
    /**
     * Get the algebraic expression that gives the interaction energy for each bond
     */
    const std::string& getEnergyFunction() const;
    /**
     * Set the algebraic expression that gives the interaction energy for each bond
     */
    void setEnergyFunction(const std::string& energy);
    /**
     * Add a new per-bond parameter that the interaction may depend on.
     *
     * @param name             the name of the parameter
     * @return the index of the parameter that was added
     */
    int addPerBondParameter(const std::string& name);
    /**
     * Get the name of a per-bond parameter.
     *
     * @param index     the index of the parameter for which to get the name
     * @return the parameter name
     */
    const std::string& getPerBondParameterName(int index) const;
    /**
     * Set the name of a per-bond parameter.
     *
     * @param index          the index of the parameter for which to set the name
     * @param name           the name of the parameter
     */
    void setPerBondParameterName(int index, const std::string& name);
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
     * Request that this Force compute the derivative of its energy with respect to a global parameter.
     * The parameter must have already been added with addGlobalParameter().
     *
     * @param name             the name of the parameter
     */
    void addEnergyParameterDerivative(const std::string& name);
    /**
     * Get the name of a global parameter with respect to which this Force should compute the
     * derivative of the energy.
     *
     * @param index     the index of the parameter derivative, between 0 and getNumEnergyParameterDerivatives()
     * @return the parameter name
     */
    const std::string& getEnergyParameterDerivativeName(int index) const;
    /**
     * Add a bond term to the force field.
     *
     * @param particle1     the index of the first particle connected by the bond
     * @param particle2     the index of the second particle connected by the bond
     * @param parameters    the list of parameters for the new bond
     * @return the index of the bond that was added
     */
    int addBond(int particle1, int particle2, const std::vector<double>& parameters=std::vector<double>());
    /**
     * Get the force field parameters for a bond term.
     *
     * @param      index         the index of the bond for which to get parameters
     * @param[out] particle1     the index of the first particle connected by the bond
     * @param[out] particle2     the index of the second particle connected by the bond
     * @param[out] parameters    the list of parameters for the bond
     */
    void getBondParameters(int index, int& particle1, int& particle2, std::vector<double>& parameters) const;
    /**
     * Set the force field parameters for a bond term.
     *
     * @param index         the index of the bond for which to set parameters
     * @param particle1     the index of the first particle connected by the bond
     * @param particle2     the index of the second particle connected by the bond
     * @param parameters    the list of parameters for the bond
     */
    void setBondParameters(int index, int particle1, int particle2, const std::vector<double>& parameters=std::vector<double>());
    /**
     * Update the per-bond parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setBondParameters() to modify this object's parameters, then call updateParametersInContext()
     * to copy them over to the Context.
     *
     * This method has several limitations.  The only information it updates is the values of per-bond parameters.
     * All other aspects of the Force (such as the energy function) are unaffected and can only be changed by reinitializing
     * the Context.  The set of particles involved in a bond cannot be changed, nor can new bonds be added.
     */
    void updateParametersInContext(Context& context);
    /**
     * Set whether this force should apply periodic boundary conditions when calculating displacements.
     * Usually this is not appropriate for bonded forces, but there are situations when it can be useful.
     */
    void setUsesPeriodicBoundaryConditions(bool periodic);
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * @returns true if force uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const;
protected:
    ForceImpl* createImpl() const;
private:
    class BondInfo;
    class BondParameterInfo;
    class GlobalParameterInfo;
    std::string energyExpression;
    std::vector<BondParameterInfo> parameters;
    std::vector<GlobalParameterInfo> globalParameters;
    std::vector<BondInfo> bonds;
    std::vector<int> energyParameterDerivatives;
    bool usePeriodic;
    mutable int numContexts, firstChangedBond, lastChangedBond;
};

/**
 * This is an internal class used to record information about a bond.
 * @private
 */
class CustomBondForce::BondInfo {
public:
    int particle1, particle2;
    std::vector<double> parameters;
    BondInfo() : particle1(-1), particle2(-1) {
    }
    BondInfo(int particle1, int particle2, const std::vector<double>& parameters) :
            particle1(particle1), particle2(particle2), parameters(parameters) {
    }
};

/**
 * This is an internal class used to record information about a per-bond parameter.
 * @private
 */
class CustomBondForce::BondParameterInfo {
public:
    std::string name;
    BondParameterInfo() {
    }
    BondParameterInfo(const std::string& name) : name(name) {
    }
};

/**
 * This is an internal class used to record information about a global parameter.
 * @private
 */
class CustomBondForce::GlobalParameterInfo {
public:
    std::string name;
    double defaultValue;
    GlobalParameterInfo() {
    }
    GlobalParameterInfo(const std::string& name, double defaultValue) : name(name), defaultValue(defaultValue) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_CUSTOMBONDEDFORCE_H_*/
