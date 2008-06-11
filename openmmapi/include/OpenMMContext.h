#ifndef OPENMM_OPENMMCONTEXT_H_
#define OPENMM_OPENMMCONTEXT_H_

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

#include "Integrator.h"
#include "State.h"
#include "System.h"
#include <string>
#include <vector>

namespace OpenMM {

class OpenMMContextImpl;
class Vec3;
class Platform;

/**
 * An OpenMMContext stores the complete state of a simulation.  More specifically, it includes:
 * 
 * <ul>
 * <li>The current time</li>
 * <li>The position of each atom</li>
 * <li>The velocity of each atom</li>
 * <li>The values of configurable parameters defined by Force objects in the System</li>
 * </ul>
 * 
 * You can retrieve a snapshot of the current state at any time by calling getState().  This
 * allows you to record the state of the simulation at various points, either for analysis
 * or for checkpointing.  getState() can also be used to retrieve the current forces on each
 * atom and the current energy of the System.
 */

class OpenMMContext {
public:
    /**
     * Construct a new OpenMMContext in which to run a simulation.
     * 
     * @param system      the System which will be simulated
     * @param integrator  the Integrator which will be used to simulate the System
     */
    OpenMMContext(System& system, Integrator& integrator);
    /**
     * Construct a new OpenMMContext in which to run a simulation, explicitly specifying what Platform should be used
     * to perform calculations.
     * 
     * @param system      the System which will be simulated
     * @param integrator  the Integrator which will be used to simulate the System
     * @param platform    the Platform to use for calculations
     */
    OpenMMContext(System& system, Integrator& integrator, Platform& platform);
    ~OpenMMContext();
    /**
     * Get System being simulated in this context.
     */
    const System& getSystem() const;
    /**
     * Get System being simulated in this context.
     */
    System& getSystem();
    /**
     * Get Integrator being used to by this context.
     */
    const Integrator& getIntegrator() const;
    /**
     * Get Integrator being used to by this context.
     */
    Integrator& getIntegrator();
    /**
     * Get a State object recording the current state information stored in this context.
     * 
     * @param types the set of data types which should be stored in the State object.  This
     * should be a union of DataType values, e.g. (State::Positions | State::Velocities).
     */
    State getState(int types) const;
    /**
     * Get the current time of the simulation (in picoseconds).
     */
    double getTime();
    /**
     * Set the current time of the simulation (in picoseconds).
     */
    void setTime(double time);
    /**
     * Set the positions of all atoms in the System (measured in nm).
     * 
     * @param a vector whose length equals the number of atoms in the System.  The i'th element
     * contains the position of the i'th atom.
     */
    void setPositions(const std::vector<Vec3>& positions);
    /**
     * Set the velocities of all atoms in the System (measured in nm/picosecond).
     * 
     * @param a vector whose length equals the number of atoms in the System.  The i'th element
     * contains the velocity of the i'th atom.
     */
    void setVelocities(const std::vector<Vec3>& velocities);
    /**
     * Get the value of an adjustable parameter defined by a Force object in the System.
     * 
     * @param name the name of the parameter to get
     */
    double getParameter(std::string name);
    /**
     * Set the value of an adjustable parameter defined by a Force object in the System.
     * 
     * @param name  the name of the parameter to set
     * @param value the value of the parameter
     */
    void setParameter(std::string name, double value);
    /**
     * When an OpenMMContext is created, it may cache information about the System being simulated
     * and the Force objects contained in it.  This means that, if the System or Forces are then
     * modified, the OpenMMContext might not see all of the changes.  Call reinitialize() to force
     * the OpenMMContext to rebuild its internal representation of the System and pick up any changes
     * that have been made.
     * 
     * This is an expensive operation, so you should try to avoid calling it too frequently.
     */
    void reinitialize();
private:
    OpenMMContextImpl* impl;
};

} // namespace OpenMM

#endif /*OPENMM_OPENMMCONTEXT_H_*/
