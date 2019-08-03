#ifndef OPENMM_CUSTOMHBONDFORCE_H_
#define OPENMM_CUSTOMHBONDFORCE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2014 Stanford University and the Authors.      *
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
#include <map>
#include <set>
#include <utility>
#include <vector>
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This class supports a wide variety of energy functions used to represent hydrogen bonding.  It computes
 * interactions between "donor" particle groups and "acceptor" particle groups, where each group may include
 * up to three particles.  Typically a donor group consists of a hydrogen atom and the atoms it is bonded to,
 * and an acceptor group consists of a negatively charged atom and the atoms it is bonded to.
 *
 * We refer to the particles in a donor group as d1, d2 and d3, and the particles in an acceptor group as
 * a1, a2, and a3.  For each donor and each acceptor, CustomHbondForce evaluates a user supplied algebraic
 * expression to determine the interaction energy.  The expression may depend on arbitrary distances, angles,
 * and dihedral angles defined by any of the six particles involved.  The function distance(p1, p2) is the distance
 * between the particles p1 and p2 (where "p1" and "p2" should be replaced by the names of the actual particles
 * to calculate the distance between), angle(p1, p2, p3) is the angle formed by the three specified particles,
 * and dihedral(p1, p2, p3, p4) is the dihedral angle formed by the four specified particles.
 *
 * The expression also may involve tabulated functions, and may depend on arbitrary
 * global, per-donor, and per-acceptor parameters.  It also optionally supports periodic boundary conditions
 * and cutoffs for long range interactions.
 *
 * To use this class, create a CustomHbondForce object, passing an algebraic expression to the constructor
 * that defines the interaction energy between each donor and acceptor.  Then call addPerDonorParameter() to define per-donor
 * parameters, addPerAcceptorParameter() to define per-acceptor parameters, and addGlobalParameter() to define
 * global parameters.  The values of per-donor and per-acceptor parameters are specified as part of the system
 * definition, while values of global parameters may be modified during a simulation by calling Context::setParameter().
 *
 * Next, call addDonor() and addAcceptor() to define donors and acceptors and specify their parameter values.
 * After a donor or acceptor has been added, you can modify its parameters by calling setDonorParameters() or
 * setAcceptorParameters().  This will have no effect on Contexts that already exist unless you call updateParametersInContext().
 *
 * CustomHbondForce also lets you specify "exclusions", particular combinations of donors and acceptors whose
 * interactions should be omitted from force and energy calculations.  This is most often used for particles
 * that are bonded to each other.
 *
 * As an example, the following code creates a CustomHbondForce that implements a simple harmonic potential
 * to keep the distance between a1 and d1, and the angle formed by a1-d1-d2, near ideal values:
 *
 * <tt>CustomHbondForce* force = new CustomHbondForce("k*(distance(a1,d1)-r0)^2*(angle(a1,d1,d2)-theta0)^2");</tt>
 *
 * This force depends on three parameters: k, r0, and theta0.  The following code defines these as per-donor parameters:
 *
 * <tt><pre>
 * force->addPerDonorParameter("k");
 * force->addPerDonorParameter("r0");
 * force->addPerDonorParameter("theta0");
 * </pre></tt>
 *
 * Expressions may involve the operators + (add), - (subtract), * (multiply), / (divide), and ^ (power), and the following
 * functions: sqrt, exp, log, sin, cos, sec, csc, tan, cot, asin, acos, atan, atan2, sinh, cosh, tanh, erf, erfc, min, max, abs, floor, ceil, step, delta, select.  All trigonometric functions
 * are defined in radians, and log is the natural logarithm.  step(x) = 0 if x is less than 0, 1 otherwise.  delta(x) = 1 if x is 0, 0 otherwise.
 * select(x,y,z) = z if x = 0, y otherwise.
 *
 * In addition, you can call addTabulatedFunction() to define a new function based on tabulated values.  You specify the function by
 * creating a TabulatedFunction object.  That function can then appear in the expression.
 */

class OPENMM_EXPORT CustomHbondForce : public Force {
public:
    /**
     * This is an enumeration of the different methods that may be used for handling long range nonbonded forces.
     */
    enum NonbondedMethod {
        /**
         * No cutoff is applied to nonbonded interactions.  The full set of N^2 interactions is computed exactly.
         * This necessarily means that periodic boundary conditions cannot be used.  This is the default.
         */
        NoCutoff = 0,
        /**
         * Interactions beyond the cutoff distance are ignored.
         */
        CutoffNonPeriodic = 1,
        /**
         * Periodic boundary conditions are used, so that each particle interacts only with the nearest periodic copy of
         * each other particle.  Interactions beyond the cutoff distance are ignored.
         */
        CutoffPeriodic = 2,
    };
    /**
     * Create a CustomHbondForce.
     *
     * @param energy    an algebraic expression giving the interaction energy between a donor and an acceptor as a function
     *                  of inter-particle distances, angles, and dihedrals, as well as any global, per-donor, and
     *                  per-acceptor parameters
     */
    explicit CustomHbondForce(const std::string& energy);
    ~CustomHbondForce();
    /**
     * Get the number of donors for which force field parameters have been defined.
     */
    int getNumDonors() const {
        return donors.size();
    }
    /**
     * Get the number of acceptors for which force field parameters have been defined.
     */
    int getNumAcceptors() const {
        return acceptors.size();
    }
    /**
     * Get the number of donor-acceptor pairs whose interactions should be excluded.
     */
    int getNumExclusions() const {
        return exclusions.size();
    }
    /**
     * Get the number of per-donor parameters that the interaction depends on.
     */
    int getNumPerDonorParameters() const {
        return donorParameters.size();
    }
    /**
     * Get the number of per-acceptor parameters that the interaction depends on.
     */
    int getNumPerAcceptorParameters() const {
        return acceptorParameters.size();
    }
    /**
     * Get the number of global parameters that the interaction depends on.
     */
    int getNumGlobalParameters() const {
        return globalParameters.size();
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
     * Get the algebraic expression that gives the interaction energy between a donor and an acceptor
     */
    const std::string& getEnergyFunction() const;
    /**
     * Set the algebraic expression that gives the interaction energy between a donor and an acceptor
     */
    void setEnergyFunction(const std::string& energy);
    /**
     * Get the method used for handling long range nonbonded interactions.
     */
    NonbondedMethod getNonbondedMethod() const;
    /**
     * Set the method used for handling long range nonbonded interactions.
     */
    void setNonbondedMethod(NonbondedMethod method);
    /**
     * Get the cutoff distance (in nm) being used.  All interactions for which the distance between d1 and a1
     * is greater than the cutoff will be ignored.  If the NonbondedMethod in use is NoCutoff, this value will have no effect.
     *
     * @return the cutoff distance, measured in nm
     */
    double getCutoffDistance() const;
    /**
     * Set the cutoff distance (in nm) being used.  All interactions for which the distance between d1 and a1
     * is greater than the cutoff will be ignored.  If the NonbondedMethod in use is NoCutoff, this value will have no effect.
     *
     * @param distance    the cutoff distance, measured in nm
     */
    void setCutoffDistance(double distance);
    /**
     * Add a new per-donor parameter that the interaction may depend on.
     *
     * @param name     the name of the parameter
     * @return the index of the parameter that was added
     */
    int addPerDonorParameter(const std::string& name);
    /**
     * Get the name of a per-donor parameter.
     *
     * @param index     the index of the parameter for which to get the name
     * @return the parameter name
     */
    const std::string& getPerDonorParameterName(int index) const;
    /**
     * Set the name of a per-donor parameter.
     *
     * @param index          the index of the parameter for which to set the name
     * @param name           the name of the parameter
     */
    void setPerDonorParameterName(int index, const std::string& name);
    /**
     * Add a new per-acceptor parameter that the interaction may depend on.
     *
     * @param name     the name of the parameter
     * @return the index of the parameter that was added
     */
    int addPerAcceptorParameter(const std::string& name);
    /**
     * Get the name of a per-acceptor parameter.
     *
     * @param index     the index of the parameter for which to get the name
     * @return the parameter name
     */
    const std::string& getPerAcceptorParameterName(int index) const;
    /**
     * Set the name of a per-acceptor parameter.
     *
     * @param index          the index of the parameter for which to set the name
     * @param name           the name of the parameter
     */
    void setPerAcceptorParameterName(int index, const std::string& name);
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
     * Add a donor group to the force
     *
     * @param d1          the index of the first particle for this donor group
     * @param d2          the index of the second particle for this donor group.  If the group only
     *                    includes one particle, this must be -1.
     * @param d3          the index of the third particle for this donor group.  If the group includes
     *                    less than three particles, this must be -1.
     * @param parameters  the list of per-donor parameter values for the new donor
     * @return the index of the donor that was added
     */
    int addDonor(int d1, int d2, int d3, const std::vector<double>& parameters=std::vector<double>());
    /**
     * Get the properties of a donor group.
     *
     * @param index            the index of the donor group to get
     * @param[out] d1          the index of the first particle for this donor group
     * @param[out] d2          the index of the second particle for this donor group.  If the group only
     *                         includes one particle, this will be -1.
     * @param[out] d3          the index of the third particle for this donor group.  If the group includes
     *                         less than three particles, this will be -1.
     * @param[out] parameters  the list of per-donor parameter values for the donor
     */
    void getDonorParameters(int index, int& d1, int& d2, int& d3, std::vector<double>& parameters) const;
    /**
     * Set the properties of a donor group.
     *
     * @param index       the index of the donor group to set
     * @param d1          the index of the first particle for this donor group
     * @param d2          the index of the second particle for this donor group.  If the group only
     *                    includes one particle, this must be -1.
     * @param d3          the index of the third particle for this donor group.  If the group includes
     *                    less than three particles, this must be -1.
     * @param parameters  the list of per-donor parameter values for the donor
     */
    void setDonorParameters(int index, int d1, int d2, int d3, const std::vector<double>& parameters=std::vector<double>());
    /**
     * Add an acceptor group to the force
     *
     * @param a1          the index of the first particle for this acceptor group
     * @param a2          the index of the second particle for this acceptor group.  If the group only
     *                    includes one particle, this must be -1.
     * @param a3          the index of the third particle for this acceptor group.  If the group includes
     *                    less than three particles, this must be -1.
     * @param parameters  the list of per-acceptor parameter values for the new acceptor
     * @return the index of the acceptor that was added
     */
    int addAcceptor(int a1, int a2, int a3, const std::vector<double>& parameters=std::vector<double>());
    /**
     * Get the properties of an acceptor group.
     *
     * @param index            the index of the acceptor group to get
     * @param[out] a1          the index of the first particle for this acceptor group
     * @param[out] a2          the index of the second particle for this acceptor group.  If the group only
     *                         includes one particle, this will be -1.
     * @param[out] a3          the index of the third particle for this acceptor group.  If the group includes
     *                         less than three particles, this will be -1.
     * @param[out] parameters  the list of per-acceptor parameter values for the acceptor
     */
    void getAcceptorParameters(int index, int& a1, int& a2, int& a3, std::vector<double>& parameters) const;
    /**
     * Set the properties of an acceptor group.
     *
     * @param index       the index of the acceptor group to set
     * @param a1          the index of the first particle for this acceptor group
     * @param a2          the index of the second particle for this acceptor group.  If the group only
     *                    includes one particle, this must be -1.
     * @param a3          the index of the third particle for this acceptor group.  If the group includes
     *                    less than three particles, this must be -1.
     * @param parameters  the list of per-acceptor parameter values for the acceptor
     */
    void setAcceptorParameters(int index, int a1, int a2, int a3, const std::vector<double>& parameters=std::vector<double>());
    /**
     * Add a donor-acceptor pair to the list of interactions that should be excluded.
     *
     * @param donor     the index of the donor to exclude
     * @param acceptor  the index of the acceptor to exclude
     * @return the index of the exclusion that was added
     */
    int addExclusion(int donor, int acceptor);
    /**
     * Get the donor and acceptor in a pair whose interaction should be excluded.
     *
     * @param index           the index of the exclusion for which to get donor and acceptor indices
     * @param[out] donor      the index of the donor
     * @param[out] acceptor   the index of the acceptor
     */
    void getExclusionParticles(int index, int& donor, int& acceptor) const;
    /**
     * Get the donor and acceptor in a pair whose interaction should be excluded.
     *
     * @param index      the index of the exclusion for which to get donor and acceptor indices
     * @param donor      the index of the donor
     * @param acceptor   the index of the acceptor
     */
    void setExclusionParticles(int index, int donor, int acceptor);
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
     * Add a tabulated function that may appear in the energy expression.
     *
     * @deprecated This method exists only for backward compatibility.  Use addTabulatedFunction() instead.
     */
    int addFunction(const std::string& name, const std::vector<double>& values, double min, double max);
    /**
     * Get the parameters for a tabulated function that may appear in the energy expression.
     *
     * @deprecated This method exists only for backward compatibility.  Use getTabulatedFunctionParameters() instead.
     * If the specified function is not a Continuous1DFunction, this throws an exception.
     */
    void getFunctionParameters(int index, std::string& name, std::vector<double>& values, double& min, double& max) const;
    /**
     * Set the parameters for a tabulated function that may appear in the energy expression.
     *
     * @deprecated This method exists only for backward compatibility.  Use setTabulatedFunctionParameters() instead.
     * If the specified function is not a Continuous1DFunction, this throws an exception.
     */
    void setFunctionParameters(int index, const std::string& name, const std::vector<double>& values, double min, double max);
    /**
     * Update the per-donor and per-acceptor parameters in a Context to match those stored in this Force object.  This method
     * provides an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setDonorParameters() and setAcceptorParameters() to modify this object's parameters, then call
     * updateParametersInContext() to copy them over to the Context.
     *
     * This method has several limitations.  The only information it updates is the values of per-donor and per-acceptor parameters.
     * All other aspects of the Force (the energy function, nonbonded method, cutoff distance, etc.) are unaffected and can only
     * be changed by reinitializing the Context.  The set of particles involved in a donor or acceptor cannot be changed, nor can
     * new donors or acceptors be added.
     */
    void updateParametersInContext(Context& context);
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * @returns true if force uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const {
        return nonbondedMethod == CustomHbondForce::CutoffPeriodic;
    }
protected:
    ForceImpl* createImpl() const;
private:
    class GroupInfo;
    class PerPairParameterInfo;
    class GlobalParameterInfo;
    class ExclusionInfo;
    class FunctionInfo;
    NonbondedMethod nonbondedMethod;
    double cutoffDistance;
    std::string energyExpression;
    std::vector<PerPairParameterInfo> donorParameters;
    std::vector<PerPairParameterInfo> acceptorParameters;
    std::vector<GlobalParameterInfo> globalParameters;
    std::vector<GroupInfo> donors;
    std::vector<GroupInfo> acceptors;
    std::vector<ExclusionInfo> exclusions;
    std::vector<FunctionInfo> functions;
};

/**
 * This is an internal class used to record information about a donor or acceptor.
 * @private
 */
class CustomHbondForce::GroupInfo {
public:
    std::vector<double> parameters;
    int p1, p2, p3;
    GroupInfo() : p1(-1), p2(-1), p3(-1) {
    }
    GroupInfo(int p1, int p2, int p3, const std::vector<double>& parameters) :
        parameters(parameters), p1(p1), p2(p2), p3(p3) {
    }
};

/**
 * This is an internal class used to record information about a per-donor or per-acceptor parameter.
 * @private
 */
class CustomHbondForce::PerPairParameterInfo {
public:
    std::string name;
    PerPairParameterInfo() {
    }
    PerPairParameterInfo(const std::string& name) : name(name) {
    }
};

/**
 * This is an internal class used to record information about a global parameter.
 * @private
 */
class CustomHbondForce::GlobalParameterInfo {
public:
    std::string name;
    double defaultValue;
    GlobalParameterInfo() {
    }
    GlobalParameterInfo(const std::string& name, double defaultValue) : name(name), defaultValue(defaultValue) {
    }
};

/**
 * This is an internal class used to record information about an exclusion.
 * @private
 */
class CustomHbondForce::ExclusionInfo {
public:
    int donor, acceptor;
    ExclusionInfo() {
        donor = acceptor = -1;
    }
    ExclusionInfo(int donor, int acceptor) :
        donor(donor), acceptor(acceptor) {
    }
};

/**
 * This is an internal class used to record information about a tabulated function.
 * @private
 */
class CustomHbondForce::FunctionInfo {
public:
    std::string name;
    TabulatedFunction* function;
    FunctionInfo() {
    }
    FunctionInfo(const std::string& name, TabulatedFunction* function) : name(name), function(function) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_CUSTOMHBONDFORCE_H_*/
