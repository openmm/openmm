#ifndef OPENMM_VARIABLEVERLETINTEGRATOR_H_
#define OPENMM_VARIABLEVERLETINTEGRATOR_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
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
#include "openmm/Kernel.h"
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This is an error contolled, variable time step Integrator that simulates a System using the
 * leap-frog Verlet algorithm.  It compares the result of the Verlet integrator to that of an
 * explicit Euler integrator, takes the difference between the two as a measure of the integration
 * error in each time step, and continuously adjusts the step size to keep the error below a
 * specified tolerance.  This both improves the stability of the integrator and allows it to take
 * larger steps on average, while still maintaining comparable accuracy to a fixed step size integrator.
 *
 * It is best not to think of the error tolerance as having any absolute meaning.  It is just an
 * adjustable parameter that affects the step size and integration accuracy.  You
 * should try different values to find the largest one that produces a trajectory sufficiently
 * accurate for your purposes.  0.001 is often a good starting point.
 *
 * Unlike a fixed step size Verlet integrator, variable step size Verlet is not symplectic.  This
 * means that at a given accuracy level, energy is not as precisely conserved over long time periods.
 * This makes it most appropriate for constant temperate simulations.  In constant energy simulations
 * where precise energy conservation over long time periods is important, a fixed step size Verlet
 * integrator may be more appropriate.
 */

class OPENMM_EXPORT VariableVerletIntegrator : public Integrator {
public:
    /**
     * Create a VariableVerletIntegrator.
     *
     * @param errorTol    the error tolerance
     */
    explicit VariableVerletIntegrator(double errorTol);
    /**
     * Get the error tolerance.
     */
    double getErrorTolerance() const {
        return errorTol;
    }
    /**
     * Set the error tolerance.
     */
    void setErrorTolerance(double tol) {
        errorTol = tol;
    }
    /**
     * Advance a simulation through time by taking a series of time steps.
     *
     * @param steps   the number of time steps to take
     */
    void step(int steps);
    /**
     * Advance a simulation through time by taking a series of steps until a specified time is
     * reached.  When this method returns, the simulation time will exactly equal the time which
     * was specified.  If you call this method and specify a time that is earlier than the
     * current time, it will return without doing anything.
     *
     * @param time   the time to which the simulation should be advanced
     */
    void stepTo(double time);
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
     * Get the names of all Kernels used by this Integrator.
     */
    std::vector<std::string> getKernelNames();
    /**
     * Compute the kinetic energy of the system at the current time.
     */
    double computeKineticEnergy();
private:
    double errorTol;
    Kernel kernel;
};

} // namespace OpenMM

#endif /*OPENMM_VARIABLEVERLETINTEGRATOR_H_*/
