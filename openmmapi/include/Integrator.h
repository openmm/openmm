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
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "Vec3.h"
#include <map>
#include <vector>

namespace OpenMM {

class OpenMMContext;

/**
 * An Integrator defines a method for simulating a System by integrating the equations of motion.
 * This is an abstract class.  Subclasses define particular integration methods.
 * 
 * Each Integrator object is bound to a particular OpenMMContext which it integrates.  This connection
 * is specified by passing the Integrator as an argument to the constructor of the OpenMMContext.
 */

class Integrator {
public:
    virtual ~Integrator() {
    }
    /**
     * Get the size of each time step, in picoseconds.
     */
    double getStepSize() const {
        return stepSize;
    }
    /**
     * Set the size of each time step, in picoseconds.
     */
    void setStepSize(double size) {
        stepSize = size;
    }
    /**
     * Advance a simulation through time by taking a series of time steps.
     * 
     * @param steps   the number of time steps to take
     */
    virtual void step(int steps) = 0;
protected:
    friend class OpenMMContextImpl;
    /**
     * This will be called by the OpenMMContext when it is created.  It informs the Integrator
     * of what context it will be integrating, and gives it a chance to do any necessary initialization.
     * It will also get called again if the application calls reinitialize() on the OpenMMContext.
     */
    virtual void initialize(OpenMMContextImpl& context) = 0;
    /**
     * Get the names of all Kernels used by this Integrator.
     */
    virtual std::vector<std::string> getKernelNames() = 0;
private:
    double stepSize;
};

} // namespace OpenMM

#endif /*OPENMM_INTEGRATOR_H_*/
