#ifndef OPENMM_RPMDINTEGRATOR_H_
#define OPENMM_RPMDINTEGRATOR_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2011 Stanford University and the Authors.      *
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

#include "openmm/Integrator.h"
#include "openmm/Kernel.h"
#include "openmm/State.h"
#include "openmm/Vec3.h"
#include "openmm/internal/windowsExport.h"

namespace OpenMM {

/**
 * This is an Integrator which simulates a System using ring polymer molecular dynamics (RPMD).
 * It simulates many copies of the System, with successive copies connected by harmonic
 * springs to form a ring.  This allows certain quantum mechanical effects to be efficiently
 * simulated.
 * 
 * Because this Integrator simulates many copies of the System at once, it must be used
 * differently from other Integrators.  Instead of setting positions and velocities by
 * calling methods of the Context, you should use the corresponding methods of the Integrator
 * to set them for specific copies of the System.  Similarly, you should retrieve state information
 * for particular copies by calling getState() on the Integrator.  Do not query the Context for
 * state information.
 */

class OPENMM_EXPORT RPMDIntegrator : public Integrator {
public:
    /**
     * Create a RPMDIntegrator.
     *
     * @param numCopies      the number of copies of the system that should be simulated
     * @param temperature    the temperature of the heat bath (in Kelvin)
     * @param frictionCoeff  the friction coefficient which couples the system to the heat bath (in inverse picoseconds)
     * @param stepSize       the step size with which to integrator the system (in picoseconds)
     */
    RPMDIntegrator(int numCopies, double temperature, double frictionCoeff, double stepSize);
    /**
     * Get the number of copies of the system being simulated.
     */
    int getNumCopies() const {
        return numCopies;
    }
    /**
     * Get the temperature of the heat bath (in Kelvin).
     *
     * @return the temperature of the heat bath, measured in Kelvin
     */
    double getTemperature() const {
        return temperature;
    }
    /**
     * Set the temperature of the heat bath (in Kelvin).
     *
     * @param temp    the temperature of the heat bath, measured in Kelvin
     */
    void setTemperature(double temp) {
        temperature = temp;
    }
    /**
     * Get the friction coefficient which determines how strongly the system is coupled to
     * the heat bath (in inverse ps).
     *
     * @return the friction coefficient, measured in 1/ps
     */
    double getFriction() const {
        return friction;
    }
    /**
     * Set the friction coefficient which determines how strongly the system is coupled to
     * the heat bath (in inverse ps).
     *
     * @param coeff    the friction coefficient, measured in 1/ps
     */
    void setFriction(double coeff) {
        friction = coeff;
    }
    /**
     * Get the random number seed.  See setRandomNumberSeed() for details.
     */
    int getRandomNumberSeed() const {
        return randomNumberSeed;
    }
    /**
     * Set the random number seed.  The precise meaning of this parameter is undefined, and is left up
     * to each Platform to interpret in an appropriate way.  It is guaranteed that if two simulations
     * are run with different random number seeds, the sequence of random forces will be different.  On
     * the other hand, no guarantees are made about the behavior of simulations that use the same seed.
     * In particular, Platforms are permitted to use non-deterministic algorithms which produce different
     * results on successive runs, even if those runs were initialized identically.
     */
    void setRandomNumberSeed(int seed) {
        randomNumberSeed = seed;
    }
    /**
     * Set the positions of all particles in one copy of the system.
     * 
     * @param copy      the index of the copy for which to set positions
     * @param positions the positions of all particles in the system
     */
    void setPositions(int copy, const std::vector<Vec3>& positions);
    /**
     * Get the velocities of all particles in one copy of the system.
     * 
     * @param copy       the index of the copy for which to set velocities
     * @param velocities the velocities of all particles in the system
     */
    void setVelocities(int copy, const std::vector<Vec3>& velocities);
    /**
     * Get a State object recording the current state information about one copy of the system.
     * 
     * @param copy  the index of the copy for which to retrieve state information
     * @param types the set of data types which should be stored in the State object.  This
     * should be a union of DataType values, e.g. (State::Positions | State::Velocities).
     */
    State getState(int copy, int types);
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
     * Get the names of all Kernels used by this Integrator.
     */
    std::vector<std::string> getKernelNames();
private:
    double temperature, friction;
    int numCopies, randomNumberSeed;
    ContextImpl* context;
    Context* owner;
    Kernel kernel;
};

} // namespace OpenMM

#endif /*OPENMM_RPMDINTEGRATOR_H_*/
