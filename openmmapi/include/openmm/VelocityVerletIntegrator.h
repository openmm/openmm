#ifndef OPENMM_VELOCITYVERLETINTEGRATOR_H_
#define OPENMM_VELOCITYVERLETINTEGRATOR_H_

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

class OPENMM_EXPORT VelocityVerletIntegrator : public Integrator {
public:
    /**
     * Create a VelocityVerletIntegrator.
     * 
     * @param stepSize the step size with which to integrate the system (in picoseconds)
     */
    explicit VelocityVerletIntegrator(double stepSize);

    virtual ~VelocityVerletIntegrator(); 
   /**
     * Advance a simulation through time by taking a series of time steps.
     * 
     * @param steps   the number of time steps to take
     */
    void step(int steps);
   /**
     * Advance a simulation through time by taking a series of time steps.
     * 
     * @param kineticEnergy  the kinetic energy of the system that the chain is thermostating
     * @param chainID        id of the Nose-Hoover-Chain
     */
    double propagateChain(double kineticEnergy, int chainID=0);
    /**
     * Add a Nose-Hoover Chain thermostat to control the temperature of the system
     *
     * @param system the system to be thermostated.  Note: this must be setup, i.e. all
     *        particles should have been added, before calling this function.
     * @param temperature the target temperature for the system.
     * @param collisionFrequency the frequency of the interaction with the heat bath (in 1/ps).
     * @param chainLength the number of beads in the Nose-Hoover chain.
     * @param numMTS the number of step in the  multiple time step chain propagation algorithm.
     * @param numYoshidaSuzuki the number of terms in the Yoshida-Suzuki multi time step decomposition
     *        used in the chain propagation algorithm (must be 1, 3, or 5).
     */
     int addNoseHooverChainThermostat(System& system, double temperature, double collisionFrequency,
                                             int chainLength, int numMTS, int numYoshidaSuzuki);
    /**
     * Add a Nose-Hoover Chain thermostat to control the temperature of the system
     *
     * @param system the system to be thermostated.  Note: this must be setup, i.e. all
     *        particles should have been added, before calling this function.
     * @param mask list of particle ids to be thermostated.
     * @param parents either an empty list of a list describing the parent atoms that each thermostated
     *        atom is connected to.
     * @param temperature the target temperature for the system.
     * @param collisionFrequency the frequency of the interaction with the heat bath (in 1/ps).
     * @param chainLength the number of beads in the Nose-Hoover chain.
     * @param numMTS the number of step in the  multiple time step chain propagation algorithm.
     * @param numYoshidaSuzuki the number of terms in the Yoshida-Suzuki multi time step decomposition
     *        used in the chain propagation algorithm (must be 1, 3, or 5).
     */
     int addMaskedNoseHooverChainThermostat(System& system, const std::vector<int>& mask, const std::vector<int>& parents,
                                             double temperature, double collisionFrequency,
                                             int chainLength, int numMTS, int numYoshidaSuzuki);
    /**
     * Get the temperature of the i-th chain (in Kelvin).
     * 
     * @param chainID the index of the Nose-Hoover chain (default=0).
     * 
     * @return the temperature.
     */
    double getTemperature(int chainID=0) const;
    /**
     * set the temperature of the i-th chain.
     *
     * @param temperature the temperature for the Nose-Hoover chain thermostat (in Kelvin).
     * @param chainID The id of the Nose-Hoover chain for which the temperature is set (default=0).
     */
    void setTemperature(double temperature, int chainID=0);
    /**
     * Get the collision frequency of the i-th chain (in 1/picosecond).
     * 
     * @param chainID the index of the Nose-Hoover chain (default=0).
     *
     * @return the collision frequency.
     */
    double getCollisionFrequency(int chainID=0) const;
    /**
     * Set the collision frequency of the i-th chain.
     *
     * @param frequency the collision frequency in picosecond.
     * @param chainID the index of the Nose-Hoover chain (default=0).
     */
    void setCollisionFrequency(double frequency, int chainID=0);
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
    virtual double computeKineticEnergy();

    std::vector<NoseHooverChain> noseHooverChains;
    bool forcesAreValid;
    Kernel vvKernel, nhcKernel;
};

} // namespace OpenMM

#endif /*OPENMM_VELOCITYVERLETINTEGRATOR_H_*/
