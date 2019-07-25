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
 * Portions copyright (c) 2008-2015 Stanford University and the Authors.      *
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
#include <map>
#include <vector>
#include "internal/windowsExport.h"

namespace OpenMM {

class Context;
class ContextImpl;

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
     */
    virtual double computeKineticEnergy() = 0;
private:
    double stepSize, constraintTol;
};

} // namespace OpenMM

#endif /*OPENMM_INTEGRATOR_H_*/
