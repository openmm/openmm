#ifndef OPENMM_CUSTOMEXTERNALFORCE_H_
#define OPENMM_CUSTOMEXTERNALFORCE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2024 Stanford University and the Authors.      *
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
 * This class implements an "external" force on particles.  The force may be applied to any subset of the particles
 * in the System.  The force on each particle is specified by an arbitrary algebraic expression, which may depend
 * on the current position of the particle as well as on arbitrary global and per-particle parameters.
 *
 * To use this class, create a CustomExternalForce object, passing an algebraic expression to the constructor
 * that defines the potential energy of each affected particle.  The expression may depend on the particle's x, y, and
 * z coordinates, as well as on any parameters you choose.  Then call addPerParticleParameter() to define per-particle
 * parameters, and addGlobalParameter() to define global parameters.  The values of per-particle parameters are specified as
 * part of the system definition, while values of global parameters may be modified during a simulation by calling Context::setParameter().
 * Finally, call addParticle() once for each particle that should be affected by the force.  After a particle has been added,
 * you can modify its parameters by calling setParticleParameters().  This will have no effect on Contexts that already exist unless
 * you call updateParametersInContext().
 *
 * As an example, the following code creates a CustomExternalForce that attracts each particle to a target position (x0, y0, z0)
 * via a harmonic potential:
 *
 * \verbatim embed:rst:leading-asterisk
 * .. code-block:: cpp
 *
 *    CustomExternalForce* force = new CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)");
 *
 * \endverbatim
 *
 * This force depends on four parameters: the spring constant k and equilibrium coordinates x0, y0, and z0.  The following code defines these parameters:
 *
 * \verbatim embed:rst:leading-asterisk
 * .. code-block:: cpp
 *
 *    force->addGlobalParameter("k", 100.0);
 *    force->addPerParticleParameter("x0");
 *    force->addPerParticleParameter("y0");
 *    force->addPerParticleParameter("z0");
 *
 * \endverbatim
 *
 * Special care is needed in systems that use periodic boundary conditions.  In that case, each particle really represents
 * an infinite set of particles repeating through space.  The variables x, y, and z contain the coordinates of one of those
 * periodic copies, but there is no guarantee about which.  It might even change from one time step to the next.  You can handle
 * this situation by using the function periodicdistance(x1, y1, z1, x2, y2, z2), which returns the minimum distance between
 * periodic copies of the points (x1, y1, z1) and (x2, y2, z2).  For example, the force given above would be rewritten as
 *
 * \verbatim embed:rst:leading-asterisk
 * .. code-block:: cpp
 *
 *    CustomExternalForce* force = new CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2");
 *
 * \endverbatim
 *
 * Expressions may involve the operators + (add), - (subtract), * (multiply), / (divide), and ^ (power), and the following
 * functions: sqrt, exp, log, sin, cos, sec, csc, tan, cot, asin, acos, atan, atan2, sinh, cosh, tanh, erf, erfc, min, max, abs, floor, ceil, step, delta, select.  All trigonometric functions
 * are defined in radians, and log is the natural logarithm.  step(x) = 0 if x is less than 0, 1 otherwise.  delta(x) = 1 if x is 0, 0 otherwise.
 * select(x,y,z) = z if x = 0, y otherwise.
 */

class OPENMM_EXPORT CustomExternalForce : public Force {
public:
    /**
     * Create a CustomExternalForce.
     *
     * @param energy    an algebraic expression giving the potential energy of each particle as a function
     *                  of its x, y, and z coordinates
     */
    explicit CustomExternalForce(const std::string& energy);
    /**
     * Get the number of particles for which force field parameters have been defined.
     */
    int getNumParticles() const {
        return particles.size();
    }
    /**
     * Get the number of per-particle parameters that the force depends on
     */
    int getNumPerParticleParameters() const {
        return parameters.size();
    }
    /**
     * Get the number of global parameters that the force depends on.
     */
    int getNumGlobalParameters() const {
        return globalParameters.size();
    }
    /**
     * Get the algebraic expression that gives the potential energy of each particle
     */
    const std::string& getEnergyFunction() const;
    /**
     * Set the algebraic expression that gives the potential energy of each particle
     */
    void setEnergyFunction(const std::string& energy);
    /**
     * Add a new per-particle parameter that the force may depend on.
     *
     * @param name             the name of the parameter
     * @return the index of the parameter that was added
     */
    int addPerParticleParameter(const std::string& name);
    /**
     * Get the name of a per-particle parameter.
     *
     * @param index     the index of the parameter for which to get the name
     * @return the parameter name
     */
    const std::string& getPerParticleParameterName(int index) const;
    /**
     * Set the name of a per-particle parameter.
     *
     * @param index          the index of the parameter for which to set the name
     * @param name           the name of the parameter
     */
    void setPerParticleParameterName(int index, const std::string& name);
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
     * Add a particle term to the force field.
     *
     * @param particle     the index of the particle this term is applied to
     * @param parameters   the list of parameters for the new force term
     * @return the index of the particle term that was added
     */
    int addParticle(int particle, const std::vector<double>& parameters=std::vector<double>());
    /**
     * Get the force field parameters for a force field term.
     *
     * @param index              the index of the particle term for which to get parameters
     * @param[out] particle      the index of the particle this term is applied to
     * @param[out] parameters    the list of parameters for the force field term
     */
    void getParticleParameters(int index, int& particle, std::vector<double>& parameters) const;
    /**
     * Set the force field parameters for a force field term.
     *
     * @param index         the index of the particle term for which to set parameters
     * @param particle      the index of the particle this term is applied to
     * @param parameters    the list of parameters for the force field term
     */
    void setParticleParameters(int index, int particle, const std::vector<double>& parameters=std::vector<double>());
    /**
     * Update the per-particle parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setParticleParameters() to modify this object's parameters, then call updateParametersInContext()
     * to copy them over to the Context.
     * 
     * This method has several limitations.  The only information it updates is the values of per-particle parameters.
     * All other aspects of the Force (such as the energy function) are unaffected and can only be changed by reinitializing
     * the Context.  Also, this method cannot be used to add new particles, only to change the parameters of existing ones.
     */
    void updateParametersInContext(Context& context);
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * @returns false
     */
    bool usesPeriodicBoundaryConditions() const;
protected:
    ForceImpl* createImpl() const;
private:
    class ParticleInfo;
    class ParticleParameterInfo;
    class GlobalParameterInfo;
    std::string energyExpression;
    std::vector<ParticleParameterInfo> parameters;
    std::vector<GlobalParameterInfo> globalParameters;
    std::vector<ParticleInfo> particles;
    mutable int numContexts, firstChangedParticle, lastChangedParticle;
};

/**
 * This is an internal class used to record information about a particle.
 * @private
 */
class CustomExternalForce::ParticleInfo {
public:
    int particle;
    std::vector<double> parameters;
    ParticleInfo() : particle(-1) {
    }
    ParticleInfo(int particle, const std::vector<double>& parameters) : particle(particle), parameters(parameters) {
    }
};

/**
 * This is an internal class used to record information about a per-particle parameter.
 * @private
 */
class CustomExternalForce::ParticleParameterInfo {
public:
    std::string name;
    ParticleParameterInfo() {
    }
    ParticleParameterInfo(const std::string& name) : name(name) {
    }
};

/**
 * This is an internal class used to record information about a global parameter.
 * @private
 */
class CustomExternalForce::GlobalParameterInfo {
public:
    std::string name;
    double defaultValue;
    GlobalParameterInfo() {
    }
    GlobalParameterInfo(const std::string& name, double defaultValue) : name(name), defaultValue(defaultValue) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_CUSTOMEXTERNALFORCE_H_*/
