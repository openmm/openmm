#ifndef OPENMM_INTEGRATOR_H_
#define OPENMM_INTEGRATOR_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2020 Stanford University and the Authors.      *
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

#include "State.h"
#include "Vec3.h"
#include "openmm/serialization/SerializationNode.h"
#include <iosfwd>
#include <map>
#include <vector>
#include "internal/windowsExport.h"

namespace OpenMM {

class Context;
class ContextImpl;
class System;

/**
 * An Integrator defines a method for simulating a System by integrating the equations of motion.
 * This is an abstract class.  Subclasses define particular integration methods.
 * 
 * Each Integrator object is bound to a particular Context which it integrates.  This connection
 * is specified by passing the Integrator as an argument to the constructor of the Context.
 */

class OPENMM_EXPORT Integrator {
public:
    Integrator();
    virtual ~Integrator();
    /**
     * Get the size of each time step, in picoseconds.  If this integrator uses variable time steps,
     * the size of the most recent step is returned.
     *
     * @return the step size, measured in ps
     */
    virtual double getStepSize() const;
    /**
     * Set the size of each time step, in picoseconds.  If this integrator uses variable time steps,
     * the effect of calling this method is undefined, and it may simply be ignored.
     *
     * @param size    the step size, measured in ps
     */
    virtual void setStepSize(double size);
    /**
     * Get the distance tolerance within which constraints are maintained, as a fraction of the constrained distance.
     */
    virtual double getConstraintTolerance() const;
    /**
     * Set the distance tolerance within which constraints are maintained, as a fraction of the constrained distance.
     */
    virtual void setConstraintTolerance(double tol);
    /**
     * Advance a simulation through time by taking a series of time steps.
     * 
     * @param steps   the number of time steps to take
     */
    virtual void step(int steps) = 0;
    /**
     * Get which force groups to use for integration.  By default, all force groups
     * are included.  This is interpreted as a set of bit flags: the forces from group i
     * will be included if (groups&(1<<i)) != 0.
     */
    virtual int getIntegrationForceGroups() const;
    /**
     * Set which force groups to use for integration.  By default, all force groups
     * are included.  This is interpreted as a set of bit flags: the forces from group i
     * will be included if (groups&(1<<i)) != 0.
     */
    virtual void setIntegrationForceGroups(int groups);
protected:
    friend class Context;
    friend class ContextImpl;
    friend class CompoundIntegrator;
    ContextImpl* context;
    Context* owner;
    /**
     * This will be called by the Context when it is created.  It informs the Integrator
     * of what context it will be integrating, and gives it a chance to do any necessary initialization.
     * It will also get called again if the application calls reinitialize() on the Context.
     */
    virtual void initialize(ContextImpl& context) = 0;
    /**
     * This will be called by the Context when it is destroyed to let the Integrator do any necessary
     * cleanup.  It will also get called again if the application calls reinitialize() on the Context.
     */
    virtual void cleanup() {
    };
    /**
     * Get the names of all Kernels used by this Integrator.
     */
    virtual std::vector<std::string> getKernelNames() = 0;
    /**
     * This will be called by the Context when the user modifies aspects of the context state, such
     * as positions, velocities, or parameters.  This gives the Integrator a chance to discard cached
     * information.  This is <i>only</i> called when the user modifies information using methods of the Context
     * object.  It is <i>not</i> called when a ForceImpl object modifies state information in its updateContextState()
     * method (unless the ForceImpl calls a Context method to perform the modification).
     * 
     * @param changed     this specifies what aspect of the Context was changed
     */
    virtual void stateChanged(State::DataType changed) {
    }
    /**
     * Compute the kinetic energy of the system at the current time.  This may be different from simply
     * mv<sup>2</sup>/2.  For example, a leapfrog integrator will store velocities offset by half a step,
     * but the kinetic energy should be computed at the current time, not delayed by half a step.
     * 
     * If kineticEnergyRequiresForce() returns true, this method can assume that valid forces
     * have already been computed.
     */
    virtual double computeKineticEnergy() = 0;
    /**
     * Get whether computeKineticEnergy() expects forces to have been computed.  The default
     * implementation returns true to be safe.  Non-leapfrog integrators can override this to
     * return false, which makes calling getState() to query the energy less expensive.
     */
    virtual bool kineticEnergyRequiresForce() const {
        return true;
    }
    /**
     * Return a list of velocities normally distributed around a target temperature.  This may be
     * overridden by Drude integrators to ensure that Drude pairs have their center of mass velocity
     * assigned as a single entity, rather than treating both particles as being independent.
     *
     * @param system the system whose velocities are to be initialized.
     * @param temperature the target temperature in Kelvin.
     * @param randomSeed the random number seed to use when selecting velocities 
     */
    virtual std::vector<Vec3> getVelocitiesForTemperature(const System &system, double temperature, int randomSeed) const;
    /**
     * Get the time interval by which velocities are offset from positions.  This is used to
     * adjust velocities when setVelocitiesToTemperature() is called on a Context.
     */
    virtual double getVelocityTimeOffset() const {
        return 0.0;
    }
    /**
     * This is called while writing checkpoints.  It gives the integrator a chance to write
     * its own data.  The default implementation does nothing.
     */
    virtual void createCheckpoint(std::ostream& stream) const {
    }
    /**
     * This is called while loading a checkpoint.  The integrator should read in whatever
     * data it wrote in createCheckpoint() and update its internal state accordingly.
     */
    virtual void loadCheckpoint(std::istream& stream) {
    }
    /**
     * This is called while creating a State.  The Integrator should store the values
     * of all time-varying parameters into the SerializationNode so they can be saved
     * as part of the state.
     */
    virtual void serializeParameters(SerializationNode& node) const {
    }
    /**
     * This is called when loading a previously saved State.  The Integrator should
     * load the values of all time-varying parameters from the SerializationNode.  If
     * the node contains parameters that are not defined for this Integrator, it should
     * throw an exception.
     */
    virtual void deserializeParameters(const SerializationNode& node) {
    }
private:
    double stepSize, constraintTol;
    int forceGroups;
};

} // namespace OpenMM

#endif /*OPENMM_INTEGRATOR_H_*/
