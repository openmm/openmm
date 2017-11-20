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
 * Portions copyright (c) 2008-2015 Stanford University and the Authors.      *
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
#include <iosfwd>
#include <map>
#include <string>
#include <vector>
#include "internal/windowsExport.h"
#include "internal/OSRngSeed.h"

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
    Context(const System& system, Integrator& integrator);
    /**
     * Construct a new Context in which to run a simulation, explicitly specifying what Platform should be used
     * to perform calculations.
     * 
     * @param system      the System which will be simulated
     * @param integrator  the Integrator which will be used to simulate the System
     * @param platform    the Platform to use for calculations
     */
    Context(const System& system, Integrator& integrator, Platform& platform);
    /**
     * Construct a new Context in which to run a simulation, explicitly specifying what Platform should be used
     * to perform calculations and the values of platform-specific properties.
     *
     * @param system      the System which will be simulated
     * @param integrator  the Integrator which will be used to simulate the System
     * @param platform    the Platform to use for calculations
     * @param properties  a set of values for platform-specific properties.  Keys are the property names.
     */
    Context(const System& system, Integrator& integrator, Platform& platform, const std::map<std::string, std::string>& properties);
    ~Context();
    /**
     * Get System being simulated in this context.
     */
    const System& getSystem() const;
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
     * @param groups a set of bit flags for which force groups to include when computing forces
     * and energies.  Group i will be included if (groups&(1<<i)) != 0.  The default value includes all groups.
     */
    State getState(int types, bool enforcePeriodicBox=false, int groups=0xFFFFFFFF) const;
    /**
     * Copy information from a State object into this Context.  This restores the Context to
     * approximately the same state it was in when the State was created.  If the State does not include
     * a piece of information (e.g. positions or velocities), that aspect of the Context is
     * left unchanged.
     * 
     * Even when all possible information is included in the State, the effect of calling this method
     * is still less complete than loadCheckpoint().  For example, it does not restore the internal
     * states of random number generators.  On the other hand, it has the advantage of not being hardware
     * specific.
     */
    void setState(const State& state);
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
     * Set the velocities of all particles in the System to random values chosen from a Boltzmann
     * distribution at a given temperature.
     * 
     * @param temperature    the temperature for which to select the velocities (measured in Kelvin)
     * @param randomSeed     the random number seed to use when selecting velocities
     */
    void setVelocitiesToTemperature(double temperature, int randomSeed=osrngseed());
    /**
     * Get all adjustable parameters that have been defined by Force objects in the System, along
     * with their current values.
     */
    const std::map<std::string, double>& getParameters() const;
    /**
     * Get the value of an adjustable parameter defined by a Force object in the System.
     * 
     * @param name the name of the parameter to get
     */
    double getParameter(const std::string& name) const;
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
     * Triclinic boxes are supported, but the vectors must satisfy certain requirements.  In particular,
     * a must point in the x direction, b must point "mostly" in the y direction, and c must point "mostly"
     * in the z direction.  See the documentation for details.
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
     * Update the velocities of particles so the net velocity of each constrained distance is zero.
     *
     * @param tol    the velocity tolerance within which constraints must be satisfied.
     */
    void applyVelocityConstraints(double tol);
    /**
     * Recompute the locations of all virtual sites.  There is rarely a reason to call
     * this, since virtual sites are also updated by applyConstraints().  This is only
     * for the rare situations when you want to enforce virtual sites but <i>not</i>
     * constraints.
     */
    void computeVirtualSites();
    /**
     * When a Context is created, it caches information about the System being simulated
     * and the Force objects contained in it.  This means that, if the System or Forces are then
     * modified, the Context does not see the changes.  Call reinitialize() to force
     * the Context to rebuild its internal representation of the System and pick up any changes
     * that have been made.
     * 
     * This is an expensive operation, so you should try to avoid calling it too frequently.
     * Most Force classes have an updateParametersInContext() method that provides a less expensive
     * way of updating certain types of information.  However, this method is the only way to
     * make some types of changes, so it is sometimes necessary to call it.
     * 
     * By default, reinitializing a Context causes all state information (positions, velocities,
     * etc.) to be discarded.  You can optionally tell it to try to preserve state information.
     * It does this by internally creating a checkpoint, then reinitializing the Context, then
     * loading the checkpoint.  Be aware that if the System has changed in a way that prevents
     * the checkpoint from being loaded (such as changing the number of particles), this will
     * throw an exception and the state information will be lost.
     */
    void reinitialize(bool preserveState=false);
    /**
     * Create a checkpoint recording the current state of the Context.  This should be treated
     * as an opaque block of binary data.  See loadCheckpoint() for more details.
     * 
     * @param stream    an output stream the checkpoint data should be written to
     */
    void createCheckpoint(std::ostream& stream);
    /**
     * Load a checkpoint that was written by createCheckpoint().
     * 
     * A checkpoint contains not only publicly visible data such as the particle positions and
     * velocities, but also internal data such as the states of random number generators.  Ideally,
     * loading a checkpoint should restore the Context to an identical state to when it was written,
     * such that continuing the simulation will produce an identical trajectory.  This is not strictly
     * guaranteed to be true, however, and should not be relied on.  For most purposes, however, the
     * internal state should be close enough to be reasonably considered equivalent.
     * 
     * A checkpoint contains data that is highly specific to the Context from which it was created.
     * It depends on the details of the System, the Platform being used, and the hardware and software
     * of the computer it was created on.  If you try to load it on a computer with different hardware,
     * or for a System that is different in any way, loading is likely to fail.  Checkpoints created
     * with different versions of OpenMM are also often incompatible.  If a checkpoint cannot be loaded,
     * that is signaled by throwing an exception.
     * 
     * @param stream    an input stream the checkpoint data should be read from
     */
    void loadCheckpoint(std::istream& stream);
    /**
     * Get a description of how the particles in the system are grouped into molecules.  Two particles are in the
     * same molecule if they are connected by constraints or bonds, where every Force object can define bonds
     * in whatever way are appropriate to that force.
     *
     * Each element lists the indices of all particles in a single molecule.  Every particle is guaranteed to
     * belong to exactly one molecule.
     */
    const std::vector<std::vector<int> >& getMolecules() const;
private:
    friend class ContextImpl;
    friend class Force;
    friend class ForceImpl;
    friend class Platform;
    Context(const System& system, Integrator& integrator, ContextImpl& linked);
    ContextImpl& getImpl();
    const ContextImpl& getImpl() const;
    ContextImpl* impl;
    std::map<std::string, std::string> properties;
};

} // namespace OpenMM

#endif /*OPENMM_CONTEXT_H_*/
