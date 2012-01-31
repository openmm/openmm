#ifndef OPENMM_CONTEXT_H_
#define OPENMM_CONTEXT_H_

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
#include <map>
#include <string>
#include <vector>
#include "internal/windowsExport.h"

namespace OpenMM {

class ContextImpl;
class Vec3;
class Platform;

/**
 * A Context stores the complete state of a simulation.  More specifically, it includes:
 * 
 * <ul>
 * <li>The current time</li>
 * <li>The position of each particle</li>
 * <li>The velocity of each particle</li>
 * <li>The values of configurable parameters defined by Force objects in the System</li>
 * </ul>
 * 
 * You can retrieve a snapshot of the current state at any time by calling getState().  This
 * allows you to record the state of the simulation at various points, either for analysis
 * or for checkpointing.  getState() can also be used to retrieve the current forces on each
 * particle and the current energy of the System.
 */

class OPENMM_EXPORT Context {
public:
    /**
     * Construct a new Context in which to run a simulation.
     * 
     * @param system      the System which will be simulated
     * @param integrator  the Integrator which will be used to simulate the System
     */
    Context(System& system, Integrator& integrator);
    /**
     * Construct a new Context in which to run a simulation, explicitly specifying what Platform should be used
     * to perform calculations.
     * 
     * @param system      the System which will be simulated
     * @param integrator  the Integrator which will be used to simulate the System
     * @param platform    the Platform to use for calculations
     */
    Context(System& system, Integrator& integrator, Platform& platform);
    /**
     * Construct a new Context in which to run a simulation, explicitly specifying what Platform should be used
     * to perform calculations and the values of platform-specific properties.
     *
     * @param system      the System which will be simulated
     * @param integrator  the Integrator which will be used to simulate the System
     * @param platform    the Platform to use for calculations
     * @param properties  a set of values for platform-specific properties.  Keys are the property names.
     */
    Context(System& system, Integrator& integrator, Platform& platform, const std::map<std::string, std::string>& properties);
    ~Context();
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
     * Get the Platform being used for calculations.
     */
    const Platform& getPlatform() const;
    /**
     * Get the Platform being used for calculations.
     */
    Platform& getPlatform();
    /**
     * Get a State object recording the current state information stored in this context.
     * 
     * @param types the set of data types which should be stored in the State object.  This
     * should be a union of DataType values, e.g. (State::Positions | State::Velocities).
     * @param enforcePeriodicBox if false, the position of each particle will be whatever position
     * is stored in the Context, regardless of periodic boundary conditions.  If true, particle
     * positions will be translated so the center of every molecule lies in the same periodic box.
     */
    State getState(int types, bool enforcePeriodicBox=false) const;
    /**
     * Set the current time of the simulation (in picoseconds).
     */
    void setTime(double time);
    /**
     * Set the positions of all particles in the System (measured in nm).  This method simply sets the positions
     * without checking to see whether they satisfy distance constraints.  If you want constraints to be
     * enforced, call applyConstraints() after setting the positions.
     * 
     * @param positions   a vector whose length equals the number of particles in the System.  The i'th element
     * contains the position of the i'th particle.
     */
    void setPositions(const std::vector<Vec3>& positions);
    /**
     * Set the velocities of all particles in the System (measured in nm/picosecond).
     * 
     * @param velocities  a vector whose length equals the number of particles in the System.  The i'th element
     * contains the velocity of the i'th particle.
     */
    void setVelocities(const std::vector<Vec3>& velocities);
    /**
     * Get the value of an adjustable parameter defined by a Force object in the System.
     * 
     * @param name the name of the parameter to get
     */
    double getParameter(const std::string& name);
    /**
     * Set the value of an adjustable parameter defined by a Force object in the System.
     * 
     * @param name  the name of the parameter to set
     * @param value the value of the parameter
     */
    void setParameter(const std::string& name, double value);
    /**
     * Set the vectors defining the axes of the periodic box (measured in nm).  They will affect
     * any Force that uses periodic boundary conditions.
     *
     * Currently, only rectangular boxes are supported.  This means that a, b, and c must be aligned with the
     * x, y, and z axes respectively.  Future releases may support arbitrary triclinic boxes.
     *
     * @param a      the vector defining the first edge of the periodic box
     * @param b      the vector defining the second edge of the periodic box
     * @param c      the vector defining the third edge of the periodic box
     */
    void setPeriodicBoxVectors(const Vec3& a, const Vec3& b, const Vec3& c);
    /**
     * Update the positions of particles so that all distance constraints are satisfied.  This also recomputes
     * the locations of all virtual sites.
     *
     * @param tol    the distance tolerance within which constraints must be satisfied.
     */
    void applyConstraints(double tol);
    /**
     * Recompute the locations of all virtual sites.  There is rarely a reason to call
     * this, since virtual sites are also updated by applyConstraints().  This is only
     * for the rare situations when you want to enforce virtual sites but <i>not</i>
     * constraints.
     */
    void computeVirtualSites();
    /**
     * When a Context is created, it may cache information about the System being simulated
     * and the Force objects contained in it.  This means that, if the System or Forces are then
     * modified, the Context might not see all of the changes.  Call reinitialize() to force
     * the Context to rebuild its internal representation of the System and pick up any changes
     * that have been made.
     * 
     * This is an expensive operation, so you should try to avoid calling it too frequently.
     */
    void reinitialize();
private:
    friend class Platform;
    ContextImpl* impl;
    std::map<std::string, std::string> properties;
};

} // namespace OpenMM

#endif /*OPENMM_CONTEXT_H_*/
