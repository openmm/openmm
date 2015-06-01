#ifndef OPENMM_CUSTOMINTEGRATOR_H_
#define OPENMM_CUSTOMINTEGRATOR_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011-2012 Stanford University and the Authors.      *
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

#include "Integrator.h"
#include "Vec3.h"
#include "openmm/Kernel.h"
#include "internal/windowsExport.h"
#include <string>
#include <vector>

namespace OpenMM {

/**
 * This is an Integrator that can be used to implemented arbitrary, user defined
 * integration algorithms.  It is flexible enough to support a wide range of
 * methods including both deterministic and stochastic integrators, Metropolized
 * integrators, and integrators that must integrate additional quantities along
 * with the particle positions and momenta.
 * 
 * To create an integration algorithm, you first define a set of variables the
 * integrator will compute.  Variables come in two types: <i>global</i> variables
 * have a single value, while <i>per-DOF</i> variables have a value for every
 * degree of freedom (x, y, or z coordinate of a particle).  You can define as
 * many variables as you want of each type.  The value of any variable can be
 * computed by the integration algorithm, or set directly by calling a method on
 * the CustomIntegrator.  All variables are persistent between integration
 * steps; once a value is set, it keeps that value until it is changed by the
 * user or recomputed in a later integration step.
 * 
 * Next, you define the algorithm as a series of computations.  To execute a
 * time step, the integrator performs the list of computations in order.  Each
 * computation updates the value of one global or per-DOF value.  There are
 * several types of computations that can be done:
 * 
 * <ul>
 * <li>Global: You provide a mathematical expression involving only global
 * variables.  It is evaluated and stored into a global variable.</li>
 * <li>Per-DOF: You provide a mathematical expression involving both global and
 * per-DOF variables.  It is evaluated once for every degree of freedom, and
 * the values are stored into a per-DOF variable.</li>
 * <li>Sum: You provide a mathematical expression involving both global and
 * per-DOF variables.  It is evaluated once for every degree of freedom.  All
 * of those values are then added together, and the sum is stored into a global
 * variable.</li>
 * <li>Constrain Positions: The particle positions are updated so that all
 * distance constraints are satisfied.</li>
 * <li>Constrain Velocities: The particle velocities are updated so the net
 * velocity along any constrained distance is 0.</li>
 * </ul>
 * 
 * Like all integrators, CustomIntegrator ignores any particle whose mass is 0.
 * It is skipped when doing per-DOF computations, and is not included when
 * computing sums over degrees of freedom.
 * 
 * In addition to the variables you define by calling addGlobalVariable() and
 * addPerDofVariable(), the integrator provides the following pre-defined
 * variables:
 * 
 * <ul>
 * <li>dt: (global) This is the step size being used by the integrator.</li>
 * <li>energy: (global, read-only) This is the current potential energy of the
 * system.</li>
 * <li>energy0, energy1, energy2, ...: (global, read-only) This is similar to
 * energy, but includes only the contribution from forces in one force group.
 * A single computation step may only depend on a single energy variable
 * (energy, energy0, energy1, etc.).</li>
 * <li>x: (per-DOF) This is the current value of the degree of freedom (the x,
 * y, or z coordinate of a particle).</li>
 * <li>v: (per-DOF) This is the current velocity associated with the degree of
 * freedom (the x, y, or z component of a particle's velocity).</li>
 * <li>f: (per-DOF, read-only) This is the current force acting on the degree of
 * freedom (the x, y, or z component of the force on a particle).</li>
 * <li>f0, f1, f2, ...: (per-DOF, read-only) This is similar to f, but includes
 * only the contribution from forces in one force group.  A single computation
 * step may only depend on a single force variable (f, f0, f1, etc.).</li>
 * <li>m: (per-DOF, read-only) This is the mass of the particle the degree of
 * freedom is associated with.</li>
 * <li>uniform: (either global or per-DOF, read-only) This is a uniformly
 * distributed random number between 0 and 1.  Every time an expression is
 * evaluated, a different value will be used.  When used in a per-DOF
 * expression, a different value will be used for every degree of freedom.
 * Note, however, that if this variable appears multiple times in a single
 * expression, the <i>same</i> value is used everywhere it appears in that
 * expression.</li>
 * <li>gaussian: (either global or per-DOF, read-only) This is a Gaussian
 * distributed random number with mean 0 and variance 1.  Every time an expression
 * is evaluated, a different value will be used.  When used in a per-DOF
 * expression, a different value will be used for every degree of freedom.
 * Note, however, that if this variable appears multiple times in a single
 * expression, the <i>same</i> value is used everywhere it appears in that
 * expression.</li>
 * <li>A global variable is created for every adjustable parameter defined
 * in the integrator's Context.</li>
 * </ul>
 * 
 * The following example uses a CustomIntegrator to implement a velocity Verlet
 * integrator:
 * 
 * <tt><pre>
 * CustomIntegrator integrator(0.001);
 * integrator.addComputePerDof("v", "v+0.5*dt*f/m");
 * integrator.addComputePerDof("x", "x+dt*v");
 * integrator.addComputePerDof("v", "v+0.5*dt*f/m");
 * </pre></tt>
 * 
 * The first step updates the velocities based on the current forces.
 * The second step updates the positions based on the new velocities, and the
 * third step updates the velocities again.  Although the first and third steps
 * look identical, the forces used in them are different.  You do not need to
 * tell the integrator that; it will recognize that the positions have changed
 * and know to recompute the forces automatically.
 * 
 * The above example has two problems.  First, it does not respect distance
 * constraints.  To make the integrator work with constraints, you need to add
 * extra steps to tell it when and how to apply them.  Second, it never gives
 * Forces an opportunity to update the context state.  This should be done every
 * time step so that, for example, an AndersenThermostat can randomize velocities
 * or a MonteCarloBarostat can scale particle positions.  You need to add a
 * step to tell the integrator when to do this.  The following example corrects
 * both these problems, using the RATTLE algorithm to apply constraints:
 * 
 * <tt><pre>
 * CustomIntegrator integrator(0.001);
 * integrator.addPerDofVariable("x1", 0);
 * integrator.addUpdateContextState();
 * integrator.addComputePerDof("v", "v+0.5*dt*f/m");
 * integrator.addComputePerDof("x", "x+dt*v");
 * integrator.addComputePerDof("x1", "x");
 * integrator.addConstrainPositions();
 * integrator.addComputePerDof("v", "v+0.5*dt*f/m+(x-x1)/dt");
 * integrator.addConstrainVelocities();
 * </pre></tt>
 * 
 * CustomIntegrator can be used to implement multiple time step integrators.  The
 * following example shows an r-RESPA integrator.  It assumes the quickly changing
 * forces are in force group 0 and the slowly changing ones are in force group 1.
 * It evaluates the "fast" forces four times as often as the "slow" forces.
 * 
 * <tt><pre>
 * CustomIntegrator integrator(0.004);
 * integrator.addComputePerDof("v", "v+0.5*dt*f1/m");
 * for (int i = 0; i &lt; 4; i++) {
 *     integrator.addComputePerDof("v", "v+0.5*(dt/4)*f0/m");
 *     integrator.addComputePerDof("x", "x+(dt/4)*v");
 *     integrator.addComputePerDof("v", "v+0.5*(dt/4)*f0/m");
 * }
 * integrator.addComputePerDof("v", "v+0.5*dt*f1/m");
 * </pre></tt>
 * 
 * An Integrator has one other job in addition to evolving the equations of motion:
 * it defines how to compute the kinetic energy of the system.  Depending on the
 * integration method used, simply summing mv<sup>2</sup>/2 over all degrees of
 * freedom may not give the correct answer.  For example, in a leapfrog integrator
 * the velocities are "delayed" by half a time step, so the above formula would
 * give the kinetic energy half a time step ago, not at the current time.
 * 
 * Call setKineticEnergyExpression() to set an expression for the kinetic energy.
 * It is computed for every degree of freedom (excluding ones whose mass is 0) and
 * the result is summed.  The default expression is "m*v*v/2", which is correct
 * for many integrators.
 * 
 * As example, the following line defines the correct way to compute kinetic energy
 * when using a leapfrog algorithm:
 * 
 * <tt><pre>
 * integrator.setKineticEnergyExpression("m*v1*v1/2; v1=v+0.5*dt*f/m");
 * </pre></tt>
 * 
 * The kinetic energy expression may depend on the following pre-defined variables:
 * x, v, f, m, dt.  It also may depend on user-defined global and per-DOF variables,
 * and on the values of adjustable parameters defined  in the integrator's Context.
 * It may <i>not</i> depend on any other variable, such as the potential energy,
 * the force from a single force group, or a random number.
 * 
 * Expressions may involve the operators + (add), - (subtract), * (multiply), / (divide), and ^ (power), and the following
 * functions: sqrt, exp, log, sin, cos, sec, csc, tan, cot, asin, acos, atan, sinh, cosh, tanh, erf, erfc, min, max, abs, floor, ceil, step, delta, select.  All trigonometric functions
 * are defined in radians, and log is the natural logarithm.  step(x) = 0 if x is less than 0, 1 otherwise.  delta(x) = 1 if x is 0, 0 otherwise.
 * select(x,y,z) = z if x = 0, y otherwise.  An expression may also involve intermediate quantities that are defined following the main expression, using ";" as a separator.
 */

class OPENMM_EXPORT CustomIntegrator : public Integrator {
public:
    /**
     * This is an enumeration of the different types of computations that may appear in an integration algorithm.
     */
    enum ComputationType {
        /**
         * Compute an expression and store it in a global variable.
         */
        ComputeGlobal = 0,
        /**
         * Compute an expression for every degree of freedom and store it in a per-DOF variable.
         */
        ComputePerDof = 1,
        /**
         * Compute an expression for every degree of freedom, sum the values, and store the result in a global variable.
         */
        ComputeSum = 2,
        /**
         * Update particle positions so all constraints are satisfied.
         */
        ConstrainPositions = 3,
        /**
         * Update particle velocities so the net velocity along all constraints is 0.
         */
        ConstrainVelocities = 4,
        /**
         * Allow Forces to update the context state.
         */
        UpdateContextState = 5,
        /**
         * Conditionally terminate the integration step based on an expression.
         */
        ConditionalTermination = 6
    };
    /**
     * Create a CustomIntegrator.
     * 
     * @param stepSize       the step size with which to integrate the system (in picoseconds)
     */
    CustomIntegrator(double stepSize);
    /**
     * Get the number of global variables that have been defined.
     */
    int getNumGlobalVariables() const {
        return globalNames.size();
    }
    /**
     * Get the number of per-DOF variables that have been defined.
     */
    int getNumPerDofVariables() const {
        return perDofNames.size();
    }
    /**
     * Get the number of computation steps that have been added.
     */
    int getNumComputations() const {
        return computations.size();
    }
    /**
     * Define a new global variable.
     * 
     * @param name          the name of the variable
     * @param initialValue  the variable will initially be set to this value
     * @return the index of the variable that was added
     */
    int addGlobalVariable(const std::string& name, double initialValue);
    /**
     * Get the name of a global variable.
     * 
     * @param index    the index of the variable to get
     * @return the name of the variable
     */
    const std::string& getGlobalVariableName(int index) const;
    /**
     * Define a new per-DOF variable.
     * 
     * @param name          the name of the variable
     * @param initialValue  the variable will initially be set to this value for
     *                      all degrees of freedom
     * @return the index of the variable that was added
     */
    int addPerDofVariable(const std::string& name, double initialValue);
    /**
     * Get the name of a per-DOF variable.
     * 
     * @param index    the index of the variable to get
     * @return the name of the variable
     */
    const std::string& getPerDofVariableName(int index) const;
    /**
     * Get the current value of a global variable.
     * 
     * @param index   the index of the variable to get
     * @return the current value of the variable
     */
    double getGlobalVariable(int index) const;
    /**
     * Get the current value of a global variable, specified by name.
     *
     * @param name    the name of the variable to get
     * @return the current value of the parameter
     */
    double getGlobalVariableByName(const std::string& name) const;
    /**
     * Set the value of a global variable.
     * 
     * @param index   the index of the variable to set
     * @param value   the new value of the variable
     */
    void setGlobalVariable(int index, double value);
    /**
     * Set the value of a global variable, specified by name.
     * 
     * @param name    the name of the variable to set
     * @param value   the new value of the variable
     */
    void setGlobalVariableByName(const std::string& name, double value);
    /**
     * Get the value of a per-DOF variable.
     * 
     * @param index   the index of the variable to get
     * @param values  the values of the variable for all degrees of freedom
     *                are stored into this
     */
    void getPerDofVariable(int index, std::vector<Vec3>& values) const;
    /**
     * Get the value of a per-DOF variable, specified by name.
     *
     * @param name    the name of the variable to get
     * @param values  the values of the variable for all degrees of freedom
     *                are stored into this
     */
    void getPerDofVariableByName(const std::string& name, std::vector<Vec3>& values) const;
    /**
     * Set the value of a per-DOF variable.
     * 
     * @param index   the index of the variable to set
     * @param values  the new values of the variable for all degrees of freedom
     */
    void setPerDofVariable(int index, const std::vector<Vec3>& values);
    /**
     * Set the value of a per-DOF variable, specified by name.
     * 
     * @param name    the name of the variable to set
     * @param values  the new values of the variable for all degrees of freedom
     */
    void setPerDofVariableByName(const std::string& name, const std::vector<Vec3>& values);
    /**
     * Add a step to the integration algorithm that computes a global value.
     * 
     * @param variable    the global variable to store the computed value into
     * @param expression  a mathematical expression involving only global variables.
     *                    In each integration step, its value is computed and
     *                    stored into the specified variable.
     * @return the index of the step that was added
     */
    int addComputeGlobal(const std::string& variable, const std::string& expression);
    /**
     * Add a step to the integration algorithm that computes a per-DOF value.
     * 
     * @param variable    the per-DOF variable to store the computed value into
     * @param expression  a mathematical expression involving both global and
     *                    per-DOF variables.  In each integration step, its value
     *                    is computed for every degree of freedom and stored into
     *                    the specified variable.
     * @return the index of the step that was added
     */
     int addConditionalTermination(const std::string& expression);
    /**
     * Add a step to the integration algorithm that conditionally terminates the step.
     *
     * @param expression  a mathematical expression involving both global and
     *                    per-DOF variables.  In each integration step, its value
     *                    is computed for every degree of freedom.
     * @return the index of the step that was added
     */
    int addComputePerDof(const std::string& variable, const std::string& expression);
    /**
     * Add a step to the integration algorithm that computes a sum over degrees of freedom.
     * 
     * @param variable    the global variable to store the computed value into
     * @param expression  a mathematical expression involving both global and
     *                    per-DOF variables.  In each integration step, its value
     *                    is computed for every degree of freedom.  Those values
     *                    are then added together, and the sum is stored in the
     *                    specified variable.
     * @return the index of the step that was added
     */
    int addComputeSum(const std::string& variable, const std::string& expression);
    /**
     * Add a step to the integration algorithm that updates particle positions so
     * all constraints are satisfied.
     * 
     * @return the index of the step that was added
     */
    int addConstrainPositions();
    /**
     * Add a step to the integration algorithm that updates particle velocities
     * so the net velocity along all constraints is 0.
     * 
     * @return the index of the step that was added
     */
    int addConstrainVelocities();
    /**
     * Add a step to the integration algorithm that allows Forces to update the
     * context state.
     * 
     * @return the index of the step that was added
     */
    int addUpdateContextState();
    /**
     * Get the details of a computation step that has been added to the integration algorithm.
     * 
     * @param index      the index of the computation step to get
     * @param type       on exit, the type of computation this step performs
     * @param variable   on exit, the variable into which this step stores its result.  If this
     *                   step does not store a result in a variable, this will be an
     *                   empty string.
     * @param expression on exit, the expression this step evaluates.  If this step does not
     *                   evaluate an expression, this will be an empty string.
     */
    void getComputationStep(int index, ComputationType& type, std::string& variable, std::string& expression) const;
    /**
     * Get the expression to use for computing the kinetic energy.  The expression is evaluated
     * for every degree of freedom.  Those values are then added together, and the sum
     * is reported as the current kinetic energy.
     */
    const std::string& getKineticEnergyExpression() const;
    /**
     * Set the expression to use for computing the kinetic energy.  The expression is evaluated
     * for every degree of freedom.  Those values are then added together, and the sum
     * is reported as the current kinetic energy.
     */
    void setKineticEnergyExpression(const std::string& expression);
    /**
     * Get the random number seed.  See setRandomNumberSeed() for details.
     */
    int getRandomNumberSeed() const {
        return randomNumberSeed;
    }
    /**
     * Set the random number seed.  The precise meaning of this parameter is undefined, and is left up
     * to each Platform to interpret in an appropriate way.  It is guaranteed that if two simulations
     * are run with different random number seeds, the sequence of random numbers will be different.  On
     * the other hand, no guarantees are made about the behavior of simulations that use the same seed.
     * In particular, Platforms are permitted to use non-deterministic algorithms which produce different
     * results on successive runs, even if those runs were initialized identically.
     *
     * If seed is set to 0 (which is the default value assigned), a unique seed is chosen when a Context
     * is created from this Force. This is done to ensure that each Context receives unique random seeds
     * without you needing to set them explicitly.
     */
    void setRandomNumberSeed(int seed) {
        randomNumberSeed = seed;
    }
    /**
     * Advance a simulation through time by taking a series of time steps.
     * 
     * @param steps   the number of time steps to take
     */
    void step(int steps);
protected:
    /**
     * This will be called by the Context when it is created.  It informs the Integrator
     * of what context it will be integrating, and gives it a chance to do any necessary initialization.
     * It will also get called again if the application calls reinitialize() on the Context.
     */
    void initialize(ContextImpl& context);
    /**
     * This will be called by the Context when it is destroyed to let the Integrator do any necessary
     * cleanup.  It will also get called again if the application calls reinitialize() on the Context.
     */
    void cleanup();
    /**
     * When the user modifies the state, we need to mark that the forces need to be recalculated.
     */
    void stateChanged(State::DataType changed);
    /**
     * Get the names of all Kernels used by this Integrator.
     */
    std::vector<std::string> getKernelNames();
    /**
     * Compute the kinetic energy of the system at the current time.
     */
    double computeKineticEnergy();
private:
    class ComputationInfo;
    std::vector<std::string> globalNames;
    std::vector<std::string> perDofNames;
    mutable std::vector<double> globalValues;
    std::vector<std::vector<Vec3> > perDofValues;
    std::vector<ComputationInfo> computations;
    std::string kineticEnergy;
    mutable bool globalsAreCurrent;
    int randomNumberSeed;
    bool forcesAreValid;
    Kernel kernel;
};

/**
 * This is an internal class used to record information about a computation step.
 * @private
 */
class CustomIntegrator::ComputationInfo {
public:
    ComputationType type;
    std::string variable, expression;
    ComputationInfo() {
    }
    ComputationInfo(ComputationType type, const std::string& variable, const std::string& expression) :
        type(type), variable(variable), expression(expression) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_CUSTOMINTEGRATOR_H_*/
