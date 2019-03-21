#ifndef OPENMM_CUSTOMCVFORCE_H_
#define OPENMM_CUSTOMCVFORCE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2018 Stanford University and the Authors.      *
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
#include "TabulatedFunction.h"
#include "internal/windowsExport.h"
#include <string>
#include <vector>

namespace OpenMM {

/**
 * This class supports energy functions that depend on collective variables.  To use it,
 * you define a set of collective variables (scalar valued functions that depend on the
 * particle positions), and an algebraic expression for the energy as a function of the
 * collective variables.  The expression also may involve tabulated functions, and may
 * depend on arbitrary global parameters.
 * 
 * Each collective variable is defined by a Force object.  The Force's potential energy
 * is computed, and that becomes the value of the variable.  This provides enormous
 * flexibility in defining collective variables, especially by using custom forces.
 * Anything that can be computed as a potential function can also be used as a collective
 * variable.
 *
 * To use this class, create a CustomCVForce object, passing an algebraic expression to the
 * constructor that defines the potential energy.  Then call addCollectiveVariable() to define
 * collective variables and addGlobalParameter() to define global parameters.  The values
 * of global parameters may be modified during a simulation by calling Context::setParameter().
 *
 * This class also has the ability to compute derivatives of the potential energy with respect to global parameters.
 * Call addEnergyParameterDerivative() to request that the derivative with respect to a particular parameter be
 * computed.  You can then query its value in a Context by calling getState() on it.
 *
 * Expressions may involve the operators + (add), - (subtract), * (multiply), / (divide), and ^ (power), and the following
 * functions: sqrt, exp, log, sin, cos, sec, csc, tan, cot, asin, acos, atan, atan2, sinh, cosh, tanh, erf, erfc, min, max, abs, floor, ceil, step, delta, select.  All trigonometric functions
 * are defined in radians, and log is the natural logarithm.  step(x) = 0 if x is less than 0, 1 otherwise.  delta(x) = 1 if x is 0, 0 otherwise.
 * select(x,y,z) = z if x = 0, y otherwise.
 *
 * In addition, you can call addTabulatedFunction() to define a new function based on tabulated values.  You specify the function by
 * creating a TabulatedFunction object.  That function can then appear in the expression.
 */

class OPENMM_EXPORT CustomCVForce : public Force {
public:
    /**
     * Create a CustomCVForce.
     *
     * @param energy   an algebraic expression giving the energy of the system as a function
     *                 of the collective variables and global parameters
     */
    explicit CustomCVForce(const std::string& energy);
    ~CustomCVForce();
    /**
     * Get the number of collective variables that the interaction depends on.
     */
    int getNumCollectiveVariables() const {
        return variables.size();
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
     * Get the number of tabulated functions that have been defined.
     */
    int getNumTabulatedFunctions() const {
        return functions.size();
    }
    /**
     * Get the algebraic expression that gives the energy of the system
     */
    const std::string& getEnergyFunction() const;
    /**
     * Set the algebraic expression that gives the energy of the system
     */
    void setEnergyFunction(const std::string& energy);
    /**
     * Add a collective variable that the force may depend on.  The collective variable
     * is represented by a Force object, which should have been created on the heap with the
     * "new" operator.  The CustomCVForce takes over ownership of it, and deletes the Force when the
     * CustomCVForce itself is deleted.
     * 
     * @param name      the name of the collective variable, as it will appear in the energy expression
     * @param variable  the collective variable, represented by a Force object.  The value of the
     *                  variable is the energy computed by the Force.
     * @return the index within the Force of the variable that was added
     */
    int addCollectiveVariable(const std::string& name, Force* variable);
    /**
     * Get the name of a collective variable.
     *
     * @param index     the index of the collective variable for which to get the name
     * @return the variable name
     */
    const std::string& getCollectiveVariableName(int index) const;
    /**
     * Get a writable reference to the Force object that computes a collective variable.
     *
     * @param index     the index of the collective variable to get
     * @return the Force object
     */
    Force& getCollectiveVariable(int index);
    /**
     * Get a const reference to the Force object that computes a collective variable.
     *
     * @param index     the index of the collective variable to get
     * @return the Force object
     */
    const Force& getCollectiveVariable(int index) const;
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
     * Add a tabulated function that may appear in the energy expression.
     *
     * @param name           the name of the function as it appears in expressions
     * @param function       a TabulatedFunction object defining the function.  The TabulatedFunction
     *                       should have been created on the heap with the "new" operator.  The
     *                       Force takes over ownership of it, and deletes it when the Force itself is deleted.
     * @return the index of the function that was added
     */
    int addTabulatedFunction(const std::string& name, TabulatedFunction* function);
    /**
     * Get a const reference to a tabulated function that may appear in the energy expression.
     *
     * @param index     the index of the function to get
     * @return the TabulatedFunction object defining the function
     */
    const TabulatedFunction& getTabulatedFunction(int index) const;
    /**
     * Get a reference to a tabulated function that may appear in the energy expression.
     *
     * @param index     the index of the function to get
     * @return the TabulatedFunction object defining the function
     */
    TabulatedFunction& getTabulatedFunction(int index);
    /**
     * Get the name of a tabulated function that may appear in the energy expression.
     *
     * @param index     the index of the function to get
     * @return the name of the function as it appears in expressions
     */
    const std::string& getTabulatedFunctionName(int index) const;
    /**
     * Get the current values of the collective variables in a Context.
     *
     * @param context        the Context for which to get the values
     * @param[out] values    the values of the collective variables are computed and
     *                       stored into this
     */
    void getCollectiveVariableValues(Context& context, std::vector<double>& values);
    /**
     * Get the inner Context used for evaluating collective variables.
     * 
     * When you create a Context for a System that contains a CustomCVForce, internally
     * it creates a new System, adds the Forces that define the CVs to it, creates a new
     * Context for that System, and uses it to evaluate the variables.  In most cases you
     * can ignore all of this.  It is just an implementation detail.  However, there are
     * a few cases where you need to directly access that internal Context.  For example,
     * if you want to modify one of the Forces that defines a collective variable and
     * call updateParametersInContext() on it, you need to pass that inner Context to it.
     * This method returns a reference to it.
     * 
     * @param context    the Context containing the CustomCVForce
     * @return the inner Context used to evaluate the collective variables
     */
    Context& getInnerContext(Context& context);
    /**
     * Update the tabulated function parameters in a Context to match those stored in this Force object.  This method
     * provides an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call getTabulatedFunction(index).setFunctionParameters() to modify this object's parameters, then call
     * updateParametersInContext() to copy them over to the Context.
     *
     * This method is very limited.  The only information it updates is the parameters of tabulated functions.
     * All other aspects of the Force (the energy expression, the set of collective variables, etc.) are unaffected and can
     * only be changed by reinitializing the Context.
     */
    void updateParametersInContext(Context& context);
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
    class GlobalParameterInfo;
    class VariableInfo;
    class FunctionInfo;
    std::string energyExpression;
    std::vector<GlobalParameterInfo> globalParameters;
    std::vector<VariableInfo> variables;
    std::vector<FunctionInfo> functions;
    std::vector<int> energyParameterDerivatives;
};

/**
 * This is an internal class used to record information about a global parameter.
 * @private
 */
class CustomCVForce::GlobalParameterInfo {
public:
    std::string name;
    double defaultValue;
    GlobalParameterInfo() {
    }
    GlobalParameterInfo(const std::string& name, double defaultValue) : name(name), defaultValue(defaultValue) {
    }
};

/**
 * This is an internal class used to record information about a tabulated function.
 * @private
 */
class CustomCVForce::VariableInfo {
public:
    std::string name;
    Force* variable;
    VariableInfo() {
    }
    VariableInfo(const std::string& name, Force* variable) : name(name), variable(variable) {
    }
};

/**
 * This is an internal class used to record information about a tabulated function.
 * @private
 */
class CustomCVForce::FunctionInfo {
public:
    std::string name;
    TabulatedFunction* function;
    FunctionInfo() {
    }
    FunctionInfo(const std::string& name, TabulatedFunction* function) : name(name), function(function) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_CUSTOMCVFORCE_H_*/
