#ifndef OPENMM_CUSTOMCENTROIDBONDFORCE_H_
#define OPENMM_CUSTOMCENTROIDBONDFORCE_H_

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

#include "TabulatedFunction.h"
#include "Force.h"
#include "Vec3.h"
#include <vector>
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This class is similar to CustomCompoundBondForce, but instead of applying forces between individual particles,
 * it applies them between the centers of groups of particles.  This is useful for a variety of purposes, such as
 * restraints to keep two molecules from moving too far apart.
 *
 * When using this class, you define groups of particles, and the center of each group is calculated as a weighted
 * average of the particle positions.  By default, the particle masses are used as weights, so the center position
 * is the center of mass.  You can optionally specify different weights to use.  You then add bonds just as with
 * CustomCompoundBondForce, but instead of specifying the particles that make up a bond, you specify the groups.
 *
 * When creating a CustomCentroidBondForce, you specify the number of groups involved in a bond, and an expression
 * for the energy of each bond.  It may depend on the center positions of individual groups, the distances between
 * the centers of pairs of groups, the angles formed by sets of three groups, and the dihedral angles formed by
 * sets of four groups.
 *
 * We refer to the groups in a bond as g1, g2, g3, etc.  For each bond, CustomCentroidBondForce evaluates a
 * user supplied algebraic expression to determine the interaction energy.  The expression may depend on the
 * following variables and functions:
 *
 * <ul>
 * <li>x1, y1, z1, x2, y2, z2, etc.: The x, y, and z coordinates of the centers of the groups.  For example, x1
 * is the x coordinate of the center of group g1, and y3 is the y coordinate of the center of group g3.</li>
 * <li>distance(g1, g2): the distance between the centers of groups g1 and g2 (where "g1" and "g2" may be replaced
 * by the names of whichever groups you want to calculate the distance between).</li>
 * <li>angle(g1, g2, g3): the angle formed by the centers of the three specified groups.</li>
 * <li>dihedral(g1, g2, g3, g4): the dihedral angle formed by the centers of the four specified groups.</li>
 * </ul>
 *
 * The expression also may involve tabulated functions, and may depend on arbitrary global and per-bond parameters.
 *
 * To use this class, create a CustomCentroidBondForce object, passing an algebraic expression to the constructor
 * that defines the interaction energy of each bond.  Then call addPerBondParameter() to define per-bond
 * parameters and addGlobalParameter() to define global parameters.  The values of per-bond parameters are specified
 * as part of the system definition, while values of global parameters may be modified during a simulation by calling
 * Context::setParameter().
 *
 * Next call addGroup() to define the particle groups.  Each group is specified by the particles it contains, and
 * the weights to use when computing the center position.
 *
 * Then call addBond() to define bonds and specify their parameter values.  After a bond has been added, you can
 * modify its parameters by calling setBondParameters().  This will have no effect on Contexts that already exist unless
 * you call updateParametersInContext().
 *
 * As an example, the following code creates a CustomCentroidBondForce that implements a harmonic force between the
 * centers of mass of two groups of particles.
 *
 * <tt><pre>
 * CustomCentroidBondForce* force = new CustomCentroidBondForce(2, "0.5*k*distance(g1,g2)^2");
 * force->addPerBondParameter("k");
 * force->addGroup(particles1);
 * force->addGroup(particles2);
 * vector<int> bondGroups;
 * bondGroups.push_back(0);
 * bondGroups.push_back(1);
 * vector<double> bondParameters;
 * bondParameters.push_back(k);
 * force->addBond(bondGroups, bondParameters);
 * </pre></tt>
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

class OPENMM_EXPORT CustomCentroidBondForce : public Force {
public:
    /**
     * Create a CustomCentroidBondForce.
     *
     * @param numGroups     the number of groups used to define each bond
     * @param energy        an algebraic expression giving the interaction energy of each bond as a function
     *                      of particle positions, inter-particle distances, angles, and dihedrals, and any global
     *                      and per-bond parameters
     */
    explicit CustomCentroidBondForce(int numGroups, const std::string& energy);
    ~CustomCentroidBondForce();
    /**
     * Get the number of groups used to define each bond.
     */
    int getNumGroupsPerBond() const {
        return groupsPerBond;
    }
    /**
     * Get the number of particle groups that have been defined.
     */
    int getNumGroups() const {
        return groups.size();
    }
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
        return bondParameters.size();
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
     * Get the number of tabulated functions that have been defined.
     *
     * @deprecated This method exists only for backward compatibility.  Use getNumTabulatedFunctions() instead.
     */
    int getNumFunctions() const {
        return functions.size();
    }
    /**
     * Get the algebraic expression that gives the interaction energy of each bond
     */
    const std::string& getEnergyFunction() const;
    /**
     * Set the algebraic expression that gives the interaction energy of each bond
     */
    void setEnergyFunction(const std::string& energy);
    /**
     * Add a new per-bond parameter that the interaction may depend on.
     *
     * @param name     the name of the parameter
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
     * Add a particle group.
     *
     * @param particles   the indices of the particles to include in the group
     * @param weights     the weight to use for each particle when computing the center position.
     *                    If this is omitted, then particle masses will be used as weights.
     * @return the index of the group that was added
     */
    int addGroup(const std::vector<int>& particles, const std::vector<double>& weights=std::vector<double>());
    /**
     * Get the properties of a group.
     *
     * @param index            the index of the group to get
     * @param[out] particles   the indices of the particles in the group
     * @param[out] weights     the weight used for each particle when computing the center position.
     *                         If no weights were specified, this vector will be empty indicating that particle
     *                         masses should be used as weights.
     */
    void getGroupParameters(int index, std::vector<int>& particles, std::vector<double>& weights) const;
    /**
     * Set the properties of a group.
     *
     * @param index       the index of the group to set
     * @param particles   the indices of the particles in the group
     * @param weights     the weight to use for each particle when computing the center position.
     *                    If this is omitted, then particle masses will be used as weights.
     */
    void setGroupParameters(int index, const std::vector<int>& particles, const std::vector<double>& weights=std::vector<double>());
    /**
     * Add a bond to the force
     *
     * @param groups      the indices of the groups the bond depends on
     * @param parameters  the list of per-bond parameter values for the new bond
     * @return the index of the bond that was added
     */
    int addBond(const std::vector<int>& groups, const std::vector<double>& parameters=std::vector<double>());
    /**
     * Get the properties of a bond.
     *
     * @param      index       the index of the bond to get
     * @param[out] groups      the indices of the groups in the bond
     * @param[out] parameters  the list of per-bond parameter values for the bond
     */
    void getBondParameters(int index, std::vector<int>& groups, std::vector<double>& parameters) const;
    /**
     * Set the properties of a bond.
     *
     * @param index       the index of the bond to set
     * @param groups      the indices of the groups in the bond
     * @param parameters  the list of per-bond parameter values for the bond
     */
    void setBondParameters(int index, const std::vector<int>& groups, const std::vector<double>& parameters=std::vector<double>());
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
     * Update the per-bond parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setBondParameters() to modify this object's parameters, then call updateParametersInContext()
     * to copy them over to the Context.
     *
     * This method has several limitations.  The only information it updates is the values of per-bond parameters.
     * All other aspects of the Force (such as the energy function) are unaffected and can only be changed by reinitializing
     * the Context.  Neither the definitions of groups nor the set of groups involved in a bond can be changed, nor can new
     * bonds be added.
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
    class GroupInfo;
    class BondInfo;
    class BondParameterInfo;
    class GlobalParameterInfo;
    class FunctionInfo;
    int groupsPerBond;
    std::string energyExpression;
    std::vector<BondParameterInfo> bondParameters;
    std::vector<GlobalParameterInfo> globalParameters;
    std::vector<GroupInfo> groups;
    std::vector<BondInfo> bonds;
    std::vector<FunctionInfo> functions;
    std::vector<int> energyParameterDerivatives;
    bool usePeriodic;
};

/**
 * This is an internal class used to record information about a group.
 * @private
 */
class CustomCentroidBondForce::GroupInfo {
public:
    std::vector<int> particles;
    std::vector<double> weights;
    GroupInfo() {
    }
    GroupInfo(const std::vector<int>& particles, const std::vector<double>& weights) :
        particles(particles), weights(weights) {
    }
};

/**
 * This is an internal class used to record information about a bond.
 * @private
 */
class CustomCentroidBondForce::BondInfo {
public:
    std::vector<int> groups;
    std::vector<double> parameters;
    BondInfo() {
    }
    BondInfo(const std::vector<int>& groups, const std::vector<double>& parameters) :
        groups(groups), parameters(parameters) {
    }
};

/**
 * This is an internal class used to record information about a per-bond parameter.
 * @private
 */
class CustomCentroidBondForce::BondParameterInfo {
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
class CustomCentroidBondForce::GlobalParameterInfo {
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
class CustomCentroidBondForce::FunctionInfo {
public:
    std::string name;
    TabulatedFunction* function;
    FunctionInfo() {
    }
    FunctionInfo(const std::string& name, TabulatedFunction* function) : name(name), function(function) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_CUSTOMCENTROIDBONDFORCE_H_*/
