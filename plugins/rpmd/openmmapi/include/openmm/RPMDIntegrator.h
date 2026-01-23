#ifndef OPENMM_RPMDINTEGRATOR_H_
#define OPENMM_RPMDINTEGRATOR_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
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
#include <set>

namespace OpenMM {

/**
 * This is an Integrator which simulates a System using ring polymer molecular dynamics (RPMD).
 * It simulates many copies of the System, with successive copies connected by harmonic
 * springs to form a ring.  This allows certain quantum mechanical effects to be efficiently
 * simulated.
 * 
 * By default this Integrator applies a PILE thermostat to the system to simulate constant
 * temperature dynamics.  You can disable the thermostat by calling setApplyThermostat(false)
 * or setThermostatType(RPMDIntegrator::NoneThermo).
 * 
 * Alternatively, you can use PILE_G mode which applies a Bussi stochastic velocity rescaling
 * (SVR) thermostat to the centroid mode while keeping PILE Langevin thermostat for internal
 * modes.  This follows the i-PI convention and can provide better ergodicity for the centroid
 * kinetic energy distribution.
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
     * This is an enumeration of the different thermostat types that can be used with RPMD.
     */
    enum ThermostatType {
        /**
         * PILE: Path Integral Langevin Equation thermostat applied to all normal modes.
         * This is the default and applies Langevin dynamics with mode-dependent
         * friction to all normal modes including the centroid.
         */
        Pile = 0,
        /**
         * PILE_G: PILE Global thermostat. Applies Bussi stochastic velocity rescaling
         * (SVR) to the centroid mode and PILE Langevin to internal modes.
         * This follows the i-PI convention and can improve ergodicity for
         * the centroid kinetic energy distribution.
         */
        PileG = 1,
        /**
         * No thermostat is applied. Useful for NVE dynamics or when an external
         * thermostat is used.
         */
        NoneThermo = 2
    };
    /**
     * This is an enumeration of the different thermostat types that can be used
     * for classical (non-RPMD) particles in hybrid quantum/classical simulations.
     */
    enum ClassicalThermostatType {
        /**
         * Bussi stochastic velocity rescaling for classical particles.
         * This maintains canonical ensemble efficiently.
         */
        BussiClassical = 0,
        /**
         * Langevin dynamics with friction matching RPMD friction coefficient.
         * Provides consistent friction with quantum sector.
         */
        LangevinClassical = 1,
        /**
         * No thermostat for classical particles (NVE dynamics).
         * Useful for microcanonical sampling.
         */
        NoneClassical = 2
    };
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
     * Get the type of thermostat being used.
     * 
     * @return the thermostat type (Pile, PileG, or NoneThermo)
     */
    ThermostatType getThermostatType() const {
        return thermostatType;
    }
    /**
     * Set the type of thermostat to use.
     * 
     * - Pile: Path Integral Langevin Equation thermostat on all normal modes (default)
     * - PileG: Bussi/SVR on centroid + PILE Langevin on internal modes (i-PI convention)
     * - NoneThermo: No thermostat (NVE dynamics)
     * 
     * Note: Setting this to NoneThermo has the same effect as setApplyThermostat(false).
     * Setting to Pile or PileG will also set applyThermostat to true.
     * 
     * @param type the thermostat type to use
     */
    void setThermostatType(ThermostatType type) {
        thermostatType = type;
        applyThermostat = (type != NoneThermo);
    }
    /**
     * Get the friction coefficient for the Bussi thermostat on the centroid mode
     * (only used when thermostatType is PILE_G), in inverse picoseconds.
     * If not explicitly set, defaults to the same value as getFriction().
     * 
     * @return the centroid thermostat friction coefficient, measured in 1/ps
     */
    double getCentroidFriction() const {
        return centroidFriction;
    }
    /**
     * Set the friction coefficient for the Bussi thermostat on the centroid mode
     * (only used when thermostatType is PILE_G), in inverse picoseconds.
     * This controls the coupling strength of the Bussi thermostat to the centroid.
     * Larger values give stronger coupling; smaller values give weaker coupling.
     * 
     * @param coeff the centroid thermostat friction coefficient, measured in 1/ps
     */
    void setCentroidFriction(double coeff) {
        centroidFriction = coeff;
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
     * Get a map whose keys are particle indices and whose values are particle types.
     * This contains only the particles that have been explicitly set with setParticleType().
     * Particles without explicit type assignment are treated as type 0.
     */
    const std::map<int, int>& getParticleTypes() const;
    /**
     * Set the type of a particle. This is an arbitrary integer that can be used to
     * group particles. Use setQuantumParticleTypes() to specify which types receive
     * quantum (RPMD) treatment.
     * 
     * @param index  the index of the particle
     * @param type   the particle type (arbitrary integer)
     */
    void setParticleType(int index, int type);
    /**
     * Get the set of particle types that receive quantum (RPMD) treatment.
     * Particles with types in this set will be propagated with multiple beads,
     * while other particles will be treated classically (single copy).
     */
    const std::set<int>& getQuantumParticleTypes() const;
    /**
     * Set which particle types receive quantum (RPMD) treatment.
     * Particles with types in this set will be propagated with multiple beads,
     * while other particles will be treated classically (single copy).
     * 
     * @param types  the set of particle types to treat as quantum
     */
    void setQuantumParticleTypes(const std::set<int>& types);
    /**
     * Get whether particles of type 0 (the default type) are treated as quantum.
     * 
     * @return true if type 0 particles get RPMD treatment, false if classical
     */
    bool getDefaultQuantum() const;
    /**
     * Set whether particles of type 0 (the default type) are treated as quantum.
     * This is convenient for setting the default behavior for particles without
     * explicit type assignment.
     * 
     * @param quantum  true to treat type 0 as quantum, false for classical
     */
    void setDefaultQuantum(bool quantum);
    /**
     * Get the type of thermostat used for classical (non-RPMD) particles.
     * 
     * @return the classical thermostat type
     */
    ClassicalThermostatType getClassicalThermostat() const;
    /**
     * Set the type of thermostat to use for classical (non-RPMD) particles.
     * Options are:
     * - BussiClassical: Bussi stochastic velocity rescaling (default)
     * - LangevinClassical: Langevin dynamics with RPMD friction coefficient
     * - NoneClassical: No thermostat (NVE dynamics)
     * 
     * @param type  the thermostat type for classical particles
     */
    void setClassicalThermostat(ClassicalThermostatType type);
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
    double temperature, friction, centroidFriction;
    int numCopies, randomNumberSeed;
    bool applyThermostat;
    ThermostatType thermostatType;
    ClassicalThermostatType classicalThermostat;
    std::map<int, int> contractions;
    std::map<int, int> particleType;
    std::set<int> quantumParticleTypes;
    bool defaultQuantum;
    bool forcesAreValid, hasSetPosition, hasSetVelocity, isFirstStep;
    Kernel kernel;
};

} // namespace OpenMM

#endif /*OPENMM_RPMDINTEGRATOR_H_*/
