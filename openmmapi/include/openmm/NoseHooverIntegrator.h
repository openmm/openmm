#ifndef OPENMM_NOSEHOOVERINTEGRATOR_H_
#define OPENMM_NOSEHOOVERINTEGRATOR_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2019 Stanford University and the Authors.           *
 * Authors: Andreas Krämer and Andrew C. Simmonett                            *
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
#include "openmm/State.h"
#include "openmm/Kernel.h"
#include "openmm/NoseHooverChain.h"
#include "internal/windowsExport.h"

namespace OpenMM {

class System;
/**
 * This is an Integrator which simulates a System using one or more Nose Hoover chain
 * thermostats, using the velocity Verlet propagation algorithm.
 */

class OPENMM_EXPORT NoseHooverIntegrator : public Integrator {
public:
    /**
     * Create a NoseHooverIntegrator.  This version creates a bare velocity Verlet integrator
     * with no thermostats; any thermostats should be added by calling addThermostat.
     * 
     * @param stepSize the step size with which to integrate the system (in picoseconds)
     */
    explicit NoseHooverIntegrator(double stepSize);
    /**
     * Create a NoseHooverIntegrator.
     *
     * @param stepSize the step size with which to integrate the system (in picoseconds)
     * @param temperature the target temperature for the system.
     * @param collisionFrequency the frequency of the interaction with the heat bath (in 1/ps).
     * @param chainLength the number of beads in the Nose-Hoover chain.
     * @param numMTS the number of step in the  multiple time step chain propagation algorithm.
     * @param numYoshidaSuzuki the number of terms in the Yoshida-Suzuki multi time step decomposition
     *        used in the chain propagation algorithm (must be 1, 3, or 5).
     */
    explicit NoseHooverIntegrator(double stepSize, double temperature, int collisionFrequnency,
                                  int chainLength = 3, int numMTS = 3, int numYoshidaSuzuki = 3);

    virtual ~NoseHooverIntegrator();
   /**
     * Advance a simulation through time by taking a series of time steps.
     * 
     * @param steps   the number of time steps to take
     */
    void step(int steps);
    /**
     * Add a simple Nose-Hoover Chain thermostat to control the temperature of the full system
     *
     * @param temperature the target temperature for the system.
     * @param collisionFrequency the frequency of the interaction with the heat bath (in 1/ps).
     * @param chainLength the number of beads in the Nose-Hoover chain.
     * @param numMTS the number of step in the  multiple time step chain propagation algorithm.
     * @param numYoshidaSuzuki the number of terms in the Yoshida-Suzuki multi time step decomposition
     *        used in the chain propagation algorithm (must be 1, 3, or 5).
     */
     int addThermostat(double temperature, double collisionFrequency,
                       int chainLength, int numMTS, int numYoshidaSuzuki);
    /**
     * Add a Nose-Hoover Chain thermostat to control the temperature of a collection of atoms and/or pairs of
     * connected atoms within the full system.  A list of atoms defining the atoms to be thermostated is
     * provided and the thermostat will only control members of that list.  Additionally a list of pairs of
     * connected atoms may be provided; in this case both the center of mass absolute motion of each pair is
     * controlled as well as their motion relative to each other, which is independently thermostated.
     *
     * @param thermostatedParticles list of particle ids to be thermostated.
     * @param thermostatedPairs a list of pairs of connected atoms whose absolute center of mass motion
     *        and motion relative to one another will be independently thermostated.
     * @param temperature the target temperature for each pair's absolute of center of mass motion.
     * @param relativeTemperature the target temperature for each pair's relative motion.
     * @param collisionFrequency the frequency of the interaction with the heat bath for the
     *        pairs' center of mass motion (in 1/ps).
     * @param relativeCollisionFrequency the frequency of the interaction with the heat bath for the
     *        pairs' relative motion (in 1/ps).
     * @param chainLength the number of beads in the Nose-Hoover chain.
     * @param numMTS the number of step in the  multiple time step chain propagation algorithm.
     * @param numYoshidaSuzuki the number of terms in the Yoshida-Suzuki multi time step decomposition
     *        used in the chain propagation algorithm (must be 1, 3, or 5).
     */
     int addSubsystemThermostat(const std::vector<int>& thermostatedParticles,
                                const std::vector< std::pair< int, int> >& thermostatedPairs,
                                double temperature, double relativeTemperature,
                                double collisionFrequency, double relativeCollisionFrequency,
                                int chainLength = 3, int numMTS = 3, int numYoshidaSuzuki = 3);
    /**
     * Get the temperature of the i-th chain for controling absolute particle motion (in Kelvin).
     * 
     * @param chainID the index of the Nose-Hoover chain (default=0).
     * 
     * @return the temperature.
     */
    double getTemperature(int chainID=0) const;
    /**
     * set the (absolute motion) temperature of the i-th chain.
     *
     * @param temperature the temperature for the Nose-Hoover chain thermostat (in Kelvin).
     * @param chainID The id of the Nose-Hoover chain for which the temperature is set (default=0).
     */
    void setTemperature(double temperature, int chainID=0);
    /**
     * Get the temperature of the i-th chain for controling pairs' relative particle motion (in Kelvin).
     *
     * @param chainID the index of the Nose-Hoover chain (default=0).
     *
     * @return the temperature.
     */
    double getRelativeTemperature(int chainID=0) const;
    /**
     * set the (relative pair motion) temperature of the i-th chain.
     *
     * @param temperature the temperature for the Nose-Hoover chain thermostat (in Kelvin).
     * @param chainID The id of the Nose-Hoover chain for which the temperature is set (default=0).
     */
    void setRelativeTemperature(double temperature, int chainID=0);
    /**
     * Get the collision frequency for absolute motion of the i-th chain (in 1/picosecond).
     * 
     * @param chainID the index of the Nose-Hoover chain (default=0).
     *
     * @return the collision frequency.
     */
    double getCollisionFrequency(int chainID=0) const;
    /**
     * Set the collision frequency for absolute motion of the i-th chain.
     *
     * @param frequency the collision frequency in picosecond.
     * @param chainID the index of the Nose-Hoover chain (default=0).
     */
    void setCollisionFrequency(double frequency, int chainID=0);
    /**
     * Get the collision frequency for pairs' relative motion of the i-th chain (in 1/picosecond).
     *
     * @param chainID the index of the Nose-Hoover chain (default=0).
     *
     * @return the collision frequency.
     */
    double getRelativeCollisionFrequency(int chainID=0) const;
    /**
     * Set the collision frequency for pairs' relative motion of the i-th chain.
     *
     * @param frequency the collision frequency in picosecond.
     * @param chainID the index of the Nose-Hoover chain (default=0).
     */
    void setRelativeCollisionFrequency(double frequency, int chainID=0);
    /**
     * Compute the total (potential + kinetic) heat bath energy for all heat baths
     * associated with this integrator, at the current time.
     */
    double computeHeatBathEnergy();
    /**
     * Get the number of Nose-Hoover chains registered with this integrator.
     */
    int getNumNoseHooverChains() const {
        return noseHooverChains.size();
    }
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
       if (State::Positions == changed) forcesAreValid = false;
    }

protected:
   /**
     * Advance any Nose-Hoover chains associated with this integrator and determine
     * scale factor for the velocities.
     * 
     * @param kineticEnergy  the {absolute, relative} kinetic energies of the system that the chain is thermostating
     * @param chainID        id of the Nose-Hoover-Chain
     * @return the scale factor to be applied to the velocities of the particles thermostated by the chain.
     */
    std::pair<double, double> propagateChain(std::pair<double, double> kineticEnergy, int chainID=0);
    /**
     * This will be called by the Context when it is created.  It informs the Integrator
     * of what context it will be integrating, and gives it a chance to do any necessary initialization.
     * It will also get called again if the application calls reinitialize() on the Context.
     */
    void initialize(ContextImpl& context);
    /**
     * Goes through the list of requested thermostat data and creates the thermostat chains.
     */
    void createThermostats(const System& system);
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
    virtual double computeKineticEnergy();

    struct ThermostatData {
        std::vector<int> thermostatedParticles;
        std::vector< std::pair< int, int> > thermostatedPairs;
        double temperature;
        double relativeTemperature;
        double collisionFrequency;
        double relativeCollisionFrequency;
        int chainLength;
        int numMTS;
        int numYoshidaSuzuki;
        ThermostatData(const std::vector<int>&particles, const std::vector<std::pair<int, int>> &pairs, double temp,
                       double relTemp, double freq, double relFreq, int length, int MTS, int YS) :
                               thermostatedParticles(particles), thermostatedPairs(pairs), temperature(temp),
                               relativeTemperature(relTemp), collisionFrequency(freq), relativeCollisionFrequency(relFreq),
                               chainLength(length), numMTS(MTS), numYoshidaSuzuki(YS) {}
    };

    std::vector<ThermostatData> thermostatData;

    std::vector<NoseHooverChain> noseHooverChains;
    bool forcesAreValid;
    Kernel vvKernel, nhcKernel;
};

} // namespace OpenMM

#endif /*OPENMM_NOSEHOOVERINTEGRATOR_H_*/
