#ifndef OPENMM_QTBINTEGRATOR_H_
#define OPENMM_QTBINTEGRATOR_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
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

namespace OpenMM {

/**
 * This integrator implements the Adaptive Quantum Thermal Bath (adQTB) algorithm
 * as described in https://doi.org/10.1021/acs.jctc.8b01164.  This is a fast method
 * for approximating nuclear quantum effects by applying a Langevin thermostat whose
 * random force varies with frequency to match the expected energy of a quantum
 * harmonic oscillator.
 * 
 * To compensate for zero point energy leakage, the spectrum of the random force
 * is adjusted automatically.  The trajectory is divided into short segments.  At
 * the end of each segment, the distribution of energy across frequencies is evaluated,
 * and the friction coefficients used in generating the random force are adjusted
 * to better match the target distribution.  You can monitor this process by calling
 * getAdaptedFriction().  It is important to equilibrate the simulation long enough
 * for the friction coefficients to converge before beginning the production part of
 * the simulation.
 * 
 * To make this process more robust, it is recommended to average the data over all
 * particles that are expected to behave identically.  To do this, you can optionally
 * call setParticleType() to define a particle as having a particular type.  The
 * data for all particles of the same type is averaged and they are adjusted together.
 * You also can set a different adaptation rate for each particle type.  The more
 * particles that are being averaged over, the higher the adaptation rate can reasonably
 * be set, leading to faster adaptation.
 * 
 * Most properties of this integrator are fixed at Context creation time and can
 * only be changed after that by reinitializing the Context.  That includes the
 * step size, friction coefficient, segment length, cutoff frequency, particle types,
 * and adaptation rates.  The only property of the integrator that can be freely
 * changed in the middle of a simulation is the temperature, and even that requires
 * some care.  The new temperature will not take effect until the beginning of the
 * next segment.  Furthermore, changing the temperature is a potentially expensive
 * operation, since it requires performing a calculation whose cost scales as the
 * cube of the number of time steps in a segment.
 * 
 * One must be very careful when trying to compute velocity dependent thermodynamic
 * quantities, such as the instantaneous temperature or instantaneous pressure.
 * The standard calculations for these quantities assume the velocities follow a
 * classical distribution.  They do not produce correct results for an adQTB
 * simulation, in which the velocities follow a quantum distribution.
 */
class OPENMM_EXPORT QTBIntegrator : public Integrator {
public:
    /**
     * Create a QTBIntegrator.
     * 
     * @param temperature    the temperature of the heat bath (in Kelvin)
     * @param frictionCoeff  the friction coefficient which couples the system to the heat bath (in inverse picoseconds)
     * @param stepSize       the step size with which to integrate the system (in picoseconds)
     */
    QTBIntegrator(double temperature, double frictionCoeff, double stepSize);
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
    void setFriction(double coeff);
    /**
     * Get the length of each segment used for adapting the noise spectrum.
     * 
     * @return the segment length, measured in ps
     */
    double getSegmentLength() const;
    /**
     * Set the length of each segment used for adapting the noise spectrum.
     * 
     * @param length    the segment length, measured in ps
     */
    void setSegmentLength(double length);
    /**
     * Get the cutoff frequency applied to the colored noise.
     * 
     * @return the cutoff frequency, measured in 1/ps
     */
    double getCutoffFrequency() const;
    /**
     * Set the cutoff frequency applied to the colored noise.
     * 
     * @param cutoff    the cutoff frequency, measured in 1/ps
     */
    void setCutoffFrequency(double cutoff);
    /**
     * Get a map whose keys are particle indices and whose values are particle types.  This
     * contains only the particles that have been specifically set with setParticleType().
     */
    const std::map<int, int>& getParticleTypes() const;
    /**
     * Set the type of a particle.  This is an arbitrary integer.
     */
    void setParticleType(int index, int type);
    /**
     * Get the default adaptation rate.  This is the rate used for particles whose
     * rate has not been otherwise specified with setTypeAdaptationRate().  It
     * determines how much the noise spectrum changes after each segment.
     */
    double getDefaultAdaptationRate() const;
    /**
     * Set the default adaptation rate.  This is the rate used for particles whose
     * rate has not been otherwise specified with setTypeAdaptationRate().  It
     * determines how much the noise spectrum changes after each segment.
     */
    void setDefaultAdaptationRate(double rate);
    /**
     * Get a map whose keys are particle types and whose values are adaptation rates.
     * These are the rates used for particles whose rates have been specified with
     * setTypeAdaptationRate().  The rate determines how much the noise spectrum
     * changes after each segment.
     */
    const std::map<int, double>& getTypeAdaptationRates() const;
    /**
     * Set the adaptation rate to use for particles of a given type.  The rate
     * determines how much the noise spectrum changes after each segment.
     */
    void setTypeAdaptationRate(int type, double rate);
    /**
     * Get the adapted friction coefficients for a particle.  The return value is
     * a vector of length numFreq = (3*n+1)/2, where n is the number of time steps
     * in a segment.  Element i is the friction coefficient used in generating
     * the random force with frequency i*pi/(numFreq*stepSize).
     * 
     * This method is guaranteed to return identical results for particles that
     * have the same type.
     * 
     * @param particle    the index of the particle for which to get the friction
     * @param friction    the adapted friction coefficients used in generating the
     *                    random force.
     */
    void getAdaptedFriction(int particle, std::vector<double>& friction) const;
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
    /**
     * This is called while writing checkpoints.  It gives the integrator a chance to write
     * its own data.
     */
    void createCheckpoint(std::ostream& stream) const;
    /**
     * This is called while loading a checkpoint.  The integrator should read in whatever
     * data it wrote in createCheckpoint() and update its internal state accordingly.
     */
    void loadCheckpoint(std::istream& stream);
    /**
     * This is called while creating a State.  The Integrator should store the values
     * of all time-varying parameters into the SerializationNode so they can be saved
     * as part of the state.
     */
    void serializeParameters(SerializationNode& node) const;
    /**
     * This is called when loading a previously saved State.  The Integrator should
     * load the values of all time-varying parameters from the SerializationNode.  If
     * the node contains parameters that are not defined for this Integrator, it should
     * throw an exception.
     */
    void deserializeParameters(const SerializationNode& node);
private:
    double temperature, friction, segmentLength, defaultAdaptationRate, cutoffFrequency;
    std::map<int, int> particleType;
    std::map<int, double> typeAdaptationRates;
    int randomNumberSeed;
    Kernel kernel;
};

} // namespace OpenMM

#endif /*OPENMM_QTBINTEGRATOR_H_*/
