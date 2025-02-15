#ifndef OPENMM_DPDINTEGRATOR_H_
#define OPENMM_DPDINTEGRATOR_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
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
#include <map>
#include <utility>
#include <vector>

namespace OpenMM {

/**
 * This integrator implements dissipative particle dynamics (DPD).  It is similar to a
 * LangevinIntegrator, but instead off applying the friction and noise forces to the
 * Cartesian coordinates of particles, they are applied to inter-particle distances.
 * 
 * Only particles within a cutoff distance apply friction to each other.  Most often
 * a single cutoff distance and friction coefficient are used for all particle pairs.
 * In some cases you may want to use different values for the interactions between
 * different types of particles.  In that case you can call setParticleType() to set
 * the types of particles and addTypePair() to set the parameters to use for their
 * interactions.
 * 
 * Particle types, type pairs, and default parameters must be set before you create
 * a Context.  Changing them will have no effect on an existing Context unless you
 * reinitialize it.
 * 
 * This integrator can be applied either to periodic or non-periodic systems.  It
 * applies periodic boundary conditions if system.usesPeriodicBoundaryConditions()
 * returns true, which happens if any force in the system uses periodic boundary
 * conditions.
 */

class OPENMM_EXPORT DPDIntegrator : public Integrator {
public:
    /**
     * Create a DPDIntegrator.  All particles default to having type 0.
     * 
     * @param temperature      the temperature of the heat bath (in Kelvin)
     * @param defaultFriction  the default friction coefficient (in inverse picoseconds).  This value is
     *                         used for interactions whose parameters have not been set with addTypePair().
     * @param defaultCutoff    the default cutoff distance (in nanometers).  This value is
     *                         used for interactions whose parameters have not been set with addTypePair().
     * @param stepSize         the step size with which to integrate the system (in picoseconds)
     */
    DPDIntegrator(double temperature, double defaultFriction, double defaultCutoff, double stepSize);
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
    void setTemperature(double temp);
    /**
     * Get the default friction coefficient for interactions that have not been specified with
     * addTypePair (in inverse ps).
     */
    double getDefaultFriction() const {
        return defaultFriction;
    }
    /**
     * Set the default friction coefficient for interactions that have not been specified with
     * addTypePair (in inverse ps).
     */
    void setDefaultFriction(double friction);
    /**
     * Get the default cutoff distance for interactions that have not been specified with
     * addTypePair (in nm).
     */
    double getDefaultCutoff() const {
        return defaultCutoff;
    }
    /**
     * Set the default cutoff distance for interactions that have not been specified with
     * addTypePair (in nm).
     */
    void setDefaultCutoff(double cutoff);
    /**
     * Get the type of a particle.  This is an arbitrary integer.  All particles initially
     * default to type 0.
     */
    int getParticleType(int index) const;
    /**
     * Set the type of a particle.  This is an arbitrary integer.  All particles initially
     * default to type 0.
     */
    void setParticleType(int index, int type);
    /**
     * Get a map whose keys are particle indices and whose values are particle types.  This
     * contains only the particles that have been specifically set with setParticleType().
     * All others have a default type of 0.
     */
    const std::map<int, int>& getParticleTypes() const;
    /**
     * Get the number of type pairs that have been added to the integrator.
     */
    int getNumTypePairs() const {
        return pairs.size();
    }
    /**
     * Add a type pair.  This overrides the default friction and cutoff distance for interactions
     * between particles of two particular types.
     * 
     * @param type1     the first particle type
     * @param type2     the second particle type
     * @param friction  the friction for interactions between particles of these two types
     * @param cutoff    the cutoff distance for interactions between particles of these two types
     * @return the index of the type pair that was just added.
     */
    int addTypePair(int type1, int type2, double friction, double cutoff);
    /**
     * Get the parameters of a type pair.  This overrides the default friction and cutoff distance
     * for interactions between particles of two particular types.
     * 
     * @param pairIndex      the index of the type pair
     * @param[out] type1     the index of the first particle type
     * @param[out] type2     the index of the second particle type
     * @param[out] friction  the friction for interactions between particles of these two types
     * @param[out] cutoff    the cutoff distance for interactions between particles of these two types
     */
    void getTypePairParameters(int pairIndex, int& type1, int& type2, double& friction, double& cutoff) const;
    /**
     * Set the parameters of a type pair.  This overrides the default friction and cutoff distance
     * for interactions between particles of two particular types.
     * 
     * @param pairIndex the index of the type pair
     * @param type1     the index of the first particle type
     * @param type2     the index of the second particle type
     * @param friction  the friction for interactions between particles of these two types
     * @param cutoff    the cutoff distance for interactions between particles of these two types
     */
    void setTypePairParameters(int pairIndex, int type1, int type2, double friction, double cutoff);
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
     * is created from this Integrator. This is done to ensure that each Context receives unique random seeds
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
     * Get the names of all Kernels used by this Integrator.
     */
    std::vector<std::string> getKernelNames();
    /**
     * Compute the kinetic energy of the system at the current time.
     */
    double computeKineticEnergy();
    /**
     * Computing kinetic energy for this integrator does not require forces.
     */
    bool kineticEnergyRequiresForce() const;
    /**
     * Get the time interval by which velocities are offset from positions.  This is used to
     * adjust velocities when setVelocitiesToTemperature() is called on a Context.
     */
    double getVelocityTimeOffset() const {
        return getStepSize()/2;
    }
private:
    class TypePairInfo;
    double temperature, defaultFriction, defaultCutoff;
    int randomNumberSeed;
    std::map<int, int> particleType;
    std::vector<TypePairInfo> pairs;
    Kernel kernel;
};

/**
 * This is an internal class used to record information about a type pair.
 * @private
 */
class DPDIntegrator::TypePairInfo {
public:
    int type1, type2;
    double friction, cutoff;
    TypePairInfo() : type1(-1), type2(-1), friction(1.0), cutoff(0.0) {
    }
    TypePairInfo(int type1, int type2, double friction, double cutoff) :
        type1(type1), type2(type2), friction(friction), cutoff(cutoff) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_DPDINTEGRATOR_H_*/
