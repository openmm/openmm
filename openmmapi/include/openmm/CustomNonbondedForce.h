#ifndef OPENMM_CUSTOMNONBONDEDFORCE_H_
#define OPENMM_CUSTOMNONBONDEDFORCE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
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
#include <map>
#include <set>
#include <utility>
#include <vector>
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This class implements nonbonded interactions between particles.  Unlike NonbondedForce, the functional form
 * of the interaction is completely customizable, and may involve arbitrary algebraic expressions and tabulated
 * functions.  It may depend on the distance between particles, as well as on arbitrary global and
 * per-particle parameters.  It also optionally supports periodic boundary conditions and cutoffs for long range interactions.
 *
 * To use this class, create a CustomNonbondedForce object, passing an algebraic expression to the constructor
 * that defines the interaction energy between each pair of particles.  The expression may depend on r, the distance
 * between the particles, as well as on any parameters you choose.  Then call addParameter() to define per-particle
 * parameters, and addGlobalParameter() to define global parameters.  When defining a per-particle parameter, you
 * specify an arbitrary algebraic expression which serves as the combining rule for calculating the parameter value
 * based on the values for the two particles involved.  The values of global parameters may be modified during a
 * simulation by calling Context::setParameter().
 * 
 * Next, call addParticle() once for each particle in the System to set the values of its per-particle parameters.
 * The number of particles for which you set parameters must be exactly equal to the number of particles in the
 * System, or else an exception will be thrown when you try to create a Context.  After a particle has been added,
 * you can modify its parameters by calling setParticleParameters().
 *
 * CustomNonbondedForce also lets you specify "exceptions", particular pairs of particles whose interactions should be
 * computed based on different parameters than those defined for the individual particles.  This can be used to
 * completely exclude certain interactions from the force calculation, or to alter how they interact with each other.
 *
 * As an example, the following code creates a CustomNonbondedForce that implements a 12-6 Lennard-Jones potential:
 *
 * <tt>CustomNonbondedForce* force = new CustomNonbondedForce("4*epsilon*((sigma/r)^12-(sigma/r)^6)");</tt>
 *
 * This force depends on two parameters: sigma and epsilon.  The following code defines these parameters, and
 * specifies combining rules for them which correspond to the standard Lorentz-Bertelot combining rules:
 *
 * <tt><pre>
 * force->addParameter("sigma", "0.5*(sigma1*sigma2)");
 * force->addParameter("epsilon", "sqrt(epsilon1*epsilon2)");
 * </pre></tt>
 *
 * Expressions may involve the operators + (add), - (subtract), * (multiply), / (divide), and ^ (power), and the following
 * functions: sqrt, exp, log, sin, cos, sec, csc, tan, cot, asin, acos, atan, sinh, cosh, tanh.  All trigonometric functions
 * are defined in radians, and log is the natural logarithm.
 *
 * In addition, you can call addFunction() to define a new function based on tabulated values.  You specify a vector of
 * values, and an interpolating or approximating spline is created from them.  That function can then appear in expressions
 * that define energy or combining rules.
 */

class OPENMM_EXPORT CustomNonbondedForce : public Force {
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
     * Create a CustomNonbondedForce.
     *
     * @param energy    an algebraic expression giving the interaction energy between two particles as a function
     *                  of r, the distance between them
     */
    CustomNonbondedForce(const std::string& energy);
    /**
     * Get the number of particles for which force field parameters have been defined.
     */
    int getNumParticles() const {
        return particles.size();
    }
    /**
     * Get the number of special interactions that should be calculated differently from other interactions.
     */
    int getNumExceptions() const {
        return exceptions.size();
    }
    /**
     * Get the number of per-particle parameters that the interaction depends on.
     */
    int getNumParameters() const {
        return parameters.size();
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
    int getNumFunctions() const {
        return functions.size();
    }
    /**
     * Get the algebraic expression that gives the interaction energy between two particles
     */
    const std::string& getEnergyFunction() const;
    /**
     * Set the algebraic expression that gives the interaction energy between two particles
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
     * Get the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
     * is NoCutoff, this value will have no effect.
     */
    double getCutoffDistance() const;
    /**
     * Set the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
     * is NoCutoff, this value will have no effect.
     */
    void setCutoffDistance(double distance);
    /**
     * Add a new per-particle parmeter that the interaction may depend on.
     *
     * @param name             the name of the parameter
     * @param combiningRule    an algebraic expression giving the combining rule for this parameter
     * @return the index of the parameter that was added
     */
    int addParameter(const std::string& name, const std::string& combiningRule);
    /**
     * Get the name of a per-particle parameter.
     *
     * @param index     the index of the parameter for which to get the name
     * @return the parameter name
     */
    const std::string& getParameterName(int index) const;
    /**
     * Set the name of a per-particle parameter.
     *
     * @param index          the index of the parameter for which to set the name
     * @param name           the name of the parameter
     */
    void setParameterName(int index, const std::string& name);
    /**
     * Get the combining rule for a per-particle parameter.
     *
     * @param index     the index of the parameter for which to get the combining rule
     * @return an algebraic expression giving the combining rule for the parameter
     */
    const std::string& getParameterCombiningRule(int index) const;
    /**
     * Set the combining rule for a per-particle parameter.
     *
     * @param index          the index of the parameter for which to set the combining rule
     * @param combiningRule  an algebraic expression giving the combining rule for the parameter
     */
    void setParameterCombiningRule(int index, const std::string& combiningRule);
    /**
     * Add a new global parmeter that the interaction may depend on.
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
     * @param name           the default value of the parameter
     */
    void setGlobalParameterDefaultValue(int index, double defaultValue);
    /**
     * Add the nonbonded force parameters for a particle.  This should be called once for each particle
     * in the System.  When it is called for the i'th time, it specifies the parameters for the i'th particle.
     *
     * @param parameters    the list of parameters for the new particle
     * @return the index of the particle that was added
     */
    int addParticle(const std::vector<double>& parameters);
    /**
     * Get the nonbonded force parameters for a particle.
     *
     * @param index       the index of the particle for which to get parameters
     * @param parameters  the list of parameters for the specified particle
     */
    void getParticleParameters(int index, std::vector<double>& parameters) const;
    /**
     * Set the nonbonded force parameters for a particle.
     *
     * @param index       the index of the particle for which to set parameters
     * @param parameters  the list of parameters for the specified particle
     */
    void setParticleParameters(int index, const std::vector<double>& parameters);
    /**
     * Add an interaction to the list of exceptions that should be calculated differently from other interactions.
     *
     * @param particle1  the index of the first particle involved in the interaction
     * @param particle2  the index of the second particle involved in the interaction
     * @param parameters the list of parameters for the new interaction.  If this is an empty (zero length) vector, it
     *                   will cause the interaction to be completely omitted from force and energy calculations.
     * @param replace    determines the behavior if there is already an exception for the same two particles.  If true, the existing one is replaced.  If false,
     *                   an exception is thrown.
     * @return the index of the exception that was added
     */
    int addException(int particle1, int particle2, const std::vector<double>& parameters, bool replace = false);
    /**
     * Get the force field parameters for an interaction that should be calculated differently from others.
     *
     * @param index      the index of the interaction for which to get parameters
     * @param particle1  the index of the first particle involved in the interaction
     * @param particle2  the index of the second particle involved in the interaction
     * @param parameters the list of parameters for the interaction.  If this is an empty (zero length) vector, it means
     *                   the interaction will be completely omitted from force and energy calculations.
     */
    void getExceptionParameters(int index, int& particle1, int& particle2, std::vector<double>& parameters) const;
    /**
     * Set the force field parameters for an interaction that should be calculated differently from others.
     *
     * @param index      the index of the interaction for which to get parameters
     * @param particle1  the index of the first particle involved in the interaction
     * @param particle2  the index of the second particle involved in the interaction
     * @param parameters the list of parameters for the interaction.  If this is an empty (zero length) vector, it
     *                   will cause the interaction to be completely omitted from force and energy calculations.
     */
    void setExceptionParameters(int index, int particle1, int particle2, const std::vector<double>& parameters);
    /**
     * Add a tabulated function that may appear in algebraic expressions.
     *
     * @param name           the name of the function as it appears in expressions
     * @param values         the tabulated values of the function f(x) at uniformly spaced values of x between min and max.
     *                       The function is assumed to be zero for x &lt; min or x &gt; max.
     * @param min            the value of the independent variable corresponding to the first element of values
     * @param max            the value of the independent variable corresponding to the last element of values
     * @param interpolating  if true, an interpolating (Catmull-Rom) spline will be used to represent the function.
     *                       If false, an approximating spline (B-spline) will be used.
     * @return the index of the function that was added
     */
    int addFunction(const std::string& name, const std::vector<double>& values, double min, double max, bool interpolating);
    /**
     * Get the parameters for a tabulated function that may appear in algebraic expressions.
     *
     * @param index          the index of the function for which to get parameters
     * @param name           the name of the function as it appears in expressions
     * @param values         the tabulated values of the function f(x) at uniformly spaced values of x between min and max.
     *                       The function is assumed to be zero for x &lt; min or x &gt; max.
     * @param min            the value of the independent variable corresponding to the first element of values
     * @param max            the value of the independent variable corresponding to the last element of values
     * @param interpolating  if true, an interpolating (Catmull-Rom) spline will be used to represent the function.
     *                       If false, an approximating spline (B-spline) will be used.
     */
    void getFunctionParameters(int index, std::string& name, std::vector<double>& values, double& min, double& max, bool& interpolating) const;
    /**
     * Set the parameters for a tabulated function that may appear in algebraic expressions.
     *
     * @param index          the index of the function for which to set parameters
     * @param name           the name of the function as it appears in expressions
     * @param values         the tabulated values of the function f(x) at uniformly spaced values of x between min and max.
     *                       The function is assumed to be zero for x &lt; min or x &gt; max.
     * @param min            the value of the independent variable corresponding to the first element of values
     * @param max            the value of the independent variable corresponding to the last element of values
     * @param interpolating  if true, an interpolating (Catmull-Rom) spline will be used to represent the function.
     *                       If false, an approximating spline (B-spline) will be used.
     */
    void setFunctionParameters(int index, const std::string& name, const std::vector<double>& values, double min, double max, bool interpolating);
protected:
    ForceImpl* createImpl();
private:
    class ParticleInfo;
    class ParameterInfo;
    class GlobalParameterInfo;
    class ExceptionInfo;
    class FunctionInfo;
    NonbondedMethod nonbondedMethod;
    double cutoffDistance;
    std::string energyExpression;
    std::vector<ParameterInfo> parameters;
    std::vector<GlobalParameterInfo> globalParameters;
    std::vector<ParticleInfo> particles;
    std::vector<ExceptionInfo> exceptions;
    std::vector<FunctionInfo> functions;
    std::map<std::pair<int, int>, int> exceptionMap;
};

class CustomNonbondedForce::ParticleInfo {
public:
    std::vector<double> parameters;
    ParticleInfo() {
    }
    ParticleInfo(const std::vector<double>& parameters) : parameters(parameters) {
    }
};

class CustomNonbondedForce::ParameterInfo {
public:
    std::string name, combiningRule;
    ParameterInfo() {
    }
    ParameterInfo(const std::string& name, const std::string& combiningRule) : name(name), combiningRule(combiningRule) {
    }
};

class CustomNonbondedForce::GlobalParameterInfo {
public:
    std::string name;
    double defaultValue;
    GlobalParameterInfo() {
    }
    GlobalParameterInfo(const std::string& name, double defaultValue) : name(name), defaultValue(defaultValue) {
    }
};

class CustomNonbondedForce::ExceptionInfo {
public:
    int particle1, particle2;
    std::vector<double> parameters;
    ExceptionInfo() {
        particle1 = particle2 = -1;
    }
    ExceptionInfo(int particle1, int particle2, const std::vector<double>& parameters) :
        particle1(particle1), particle2(particle2), parameters(parameters) {
    }
};

class CustomNonbondedForce::FunctionInfo {
public:
    std::string name;
    std::vector<double> values;
    double min, max;
    bool interpolating;
    FunctionInfo() {
    }
    FunctionInfo(const std::string& name, const std::vector<double>& values, double min, double max, bool interpolating) :
        name(name), values(values), min(min), max(max), interpolating(interpolating) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_CUSTOMNONBONDEDFORCE_H_*/
