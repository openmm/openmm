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

#include "openmm/Integrator.h"
#include "openmm/Kernel.h"
#include "openmm/State.h"
#include "openmm/Vec3.h"
#include "openmm/internal/windowsExportRpmd.h"
#include <map>

namespace OpenMM {

/**
 * This is an Integrator which simulates a System using ring polymer molecular dynamics (RPMD).
 * It simulates many copies of the System, with successive copies connected by harmonic
 * springs to form a ring.  This allows certain quantum mechanical effects to be efficiently
 * simulated.
 * 
 * By default this Integrator applies a PILE thermostat to the system to simulate constant
 * temperature dynamics.  You can disable the thermostat by calling setApplyThermostat(false).
 * 
 * Because this Integrator simulates many copies of the System at once, it must be used
 * differently from other Integrators.  Instead of setting positions and velocities by
 * calling methods of the Context, you should use the corresponding methods of the Integrator
 * to set them for specific copies of the System.  Similarly, you should retrieve state information
 * for particular copies by calling getState() on the Integrator.  Do not query the Context for
 * state information.
 * 
 * You can optionally specify a set of "ring polymer contractions", by which different force
 * groups are evaluated on different numbers of copies, instead of computing every force on
 * every copy.  This can be much more efficient, since different forces may vary widely in
 * how many times they must be evaluated to produce sufficient accuracy.  For example, you
 * might simulate a 32 copy ring polymer and evaluate bonded forces on every copy, but contract
 * it down to only 6 copies for computing nonbonded interactions, and down to only a single
 * copy (the centroid) for computing the reciprocal space part of PME.
 */

class OPENMM_EXPORT_RPMD RPMDIntegrator : public Integrator {
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
     * Create a RPMDIntegrator.
     *
     * @param numCopies      the number of copies of the system that should be simulated
     * @param temperature    the temperature of the heat bath (in Kelvin)
     * @param frictionCoeff  the friction coefficient which couples the system to the heat bath (in inverse picoseconds)
     * @param stepSize       the step size with which to integrator the system (in picoseconds)
     * @param contractions   the ring polymer contractions to use for evaluating different force groups.  Each key in the
     *                       map is the index of a force group, and the corresponding value is the number of copies to evaluate
     *                       that force group on.  If no entry is provided for a force group (the default), it is evaluated
     *                       independently on every copy.
     */
    RPMDIntegrator(int numCopies, double temperature, double frictionCoeff, double stepSize, const std::map<int, int>& contractions);
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
     * Get whether a thermostat is applied to the system.
     */
    bool getApplyThermostat() const {
        return applyThermostat;
    }
    /**
     * Set whether a thermostat is applied to the system.
     */
    void setApplyThermostat(bool apply) {
        applyThermostat = apply;
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
     *
     * If seed is set to 0 (which is the default value assigned), a unique seed is chosen when a Context
     * is created from this Force. This is done to ensure that each Context receives unique random seeds
     * without you needing to set them explicitly.
     */
    void setRandomNumberSeed(int seed) {
        randomNumberSeed = seed;
    }
    /**
     * Get the ring polymer contractions to use for evaluating different force groups.  Each key in the
     * map is the index of a force group, and the corresponding value is the number of copies to evaluate
     * that force group on.  If no entry is provided for a force group, it is evaluated independently on
     * every copy.
     */
    const std::map<int, int>& getContractions() const {
        return contractions;
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
     * @param enforcePeriodicBox if false, the position of each particle will be whatever position
     * is stored by the integrator, regardless of periodic boundary conditions.  If true, particle
     * positions will be translated so the center of every molecule lies in the same periodic box.
     * @param groups a set of bit flags for which force groups to include when computing forces
     * and energies.  Group i will be included if (groups&(1<<i)) != 0.  The default value includes all groups.
     */
    State getState(int copy, int types, bool enforcePeriodicBox=false, int groups=0xFFFFFFFF);
    /**
     * Get the total energy of the ring polymer.  This includes the potential and kinetic energies of all copies,
     * plus the potential energy of the harmonic springs that link copies together.
     */
    double getTotalEnergy();
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
    double temperature, friction;
    int numCopies, randomNumberSeed;
    bool applyThermostat;
    std::map<int, int> contractions;
    bool forcesAreValid, hasSetPosition, hasSetVelocity, isFirstStep;
    Kernel kernel;
};

} // namespace OpenMM

#endif /*OPENMM_RPMDINTEGRATOR_H_*/
