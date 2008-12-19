#ifndef OPENMM_OPENMMCONTEXTIMPL_H_
#define OPENMM_OPENMMCONTEXTIMPL_H_

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

#include "Kernel.h"
#include "Platform.h"
#include "Stream.h"
#include <map>
#include <vector>

namespace OpenMM {

class ForceImpl;
class Integrator;
class OpenMMContext;
class System;

/**
 * This is the internal implementation of an OpenMMContext.
 */

class OPENMM_EXPORT OpenMMContextImpl {
public:
    /**
     * Create an OpenMMContextImpl for an OpenMMContext;
     */
    OpenMMContextImpl(OpenMMContext& owner, System& system, Integrator& integrator, Platform* platform);
    ~OpenMMContextImpl();
    /**
     * Get the OpenMMContext for which this is the implementation.
     */
    OpenMMContext& getOwner() {
        return owner;
    }
    /**
     * Get System being simulated in this context.
     */
    System& getSystem() {
        return system;
    }
    /**
     * Get Integrator being used to by this context.
     */
    Integrator& getIntegrator() {
        return integrator;
    }
    /**
     * Get the Platform implementation being used for computations.
     */
    Platform& getPlatform() {
        return *platform;
    }
    /**
     * Get the Stream containing the current position of each particle.
     */
    Stream& getPositions() {
        return positions;
    }
    /**
     * Get the Stream containing the current velocity of each particle.
     */
    Stream& getVelocities() {
        return velocities;
    }
    /**
     * Get the Stream containing the force on each particle that was calculated by
     * the most recent call to calcForces().
     */
    Stream& getForces() {
        return forces;
    }
    /**
     * Get the current time (in picoseconds).
     */
    double getTime() const {
        return time;
    }
    /**
     * Set the current time (in picoseconds).
     */
    void setTime(double t) {
        time = t;
    }
    /**
     * Get the value of an adjustable parameter.  If there is no parameter with the specified name, this
     * throws an exception.
     * 
     * @param name the name of the parameter to get
     */
    double getParameter(std::string name);
    /**
     * Set the value of an adjustable parameter.  If there is no parameter with the specified name, this
     * throws an exception.
     * 
     * @param name  the name of the parameter to set
     * @param value the value of the parameter
     */
    void setParameter(std::string name, double value);
    /**
     * Recalculate all of the forces in the system.  After calling this, use getForces() to retrieve
     * the forces that were calculated.
     */
    void calcForces();
    /**
     * Calculate the kinetic energy of the system (in kJ/mol).
     */
    double calcKineticEnergy();
    /**
     * Calculate the potential energy of the system (in kJ/mol).
     */
    double calcPotentialEnergy();
    /**
     * This should be called at the start of each time step.  It calls updateContextState() on each
     * ForceImpl in the system, allowing them to modify the values of state variables.
     */
    void updateContextState();
    /**
     * Get the platform-specific data stored in this context.
     */
    void* getPlatformData();
    /**
     * Set the platform-specific data stored in this context.
     */
    void setPlatformData(void* data);
private:
    friend class OpenMMContext;
    OpenMMContext& owner;
    System& system;
    Integrator& integrator;
    std::vector<ForceImpl*> forceImpls;
    double time;
    std::map<std::string, double> parameters;
    Platform* platform;
    Stream positions, velocities, forces;
    Kernel initializeForcesKernel, kineticEnergyKernel;
    void* platformData;
};

} // namespace OpenMM

#endif /*OPENMM_OPENMMCONTEXTIMPL_H_*/
