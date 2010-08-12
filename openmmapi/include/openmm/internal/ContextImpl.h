#ifndef OPENMM_CONTEXTIMPL_H_
#define OPENMM_CONTEXTIMPL_H_

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

#include "openmm/Kernel.h"
#include "openmm/Platform.h"
#include "openmm/Vec3.h"
#include <map>
#include <vector>

namespace OpenMM {

class ForceImpl;
class Integrator;
class Context;
class System;

/**
 * This is the internal implementation of a Context.
 */

class OPENMM_EXPORT ContextImpl {
public:
    /**
     * Create an ContextImpl for a Context;
     */
    ContextImpl(Context& owner, System& system, Integrator& integrator, Platform* platform, const std::map<std::string, std::string>& properties);
    ~ContextImpl();
    /**
     * Get the Context for which this is the implementation.
     */
    Context& getOwner() {
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
     * Get the current time (in picoseconds).
     */
    double getTime() const;
    /**
     * Set the current time (in picoseconds).
     */
    void setTime(double t);
    /**
     * Get the positions of all particles.
     *
     * @param positions  on exit, this contains the particle positions
     */
    void getPositions(std::vector<Vec3>& positions);
    /**
     * Set the positions of all particles.
     *
     * @param positions  a vector containg the particle positions
     */
    void setPositions(const std::vector<Vec3>& positions);
    /**
     * Get the velocities of all particles.
     *
     * @param velocities  on exit, this contains the particle velocities
     */
    void getVelocities(std::vector<Vec3>& velocities);
    /**
     * Set the velocities of all particles.
     *
     * @param velocities  a vector containg the particle velocities
     */
    void setVelocities(const std::vector<Vec3>& velocities);
    /**
     * Get the current forces on all particles.
     *
     * @param forces  on exit, this contains the forces
     */
    void getForces(std::vector<Vec3>& forces);
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
     * Get the vectors defining the axes of the periodic box (measured in nm).  They will affect
     * any Force that uses periodic boundary conditions.
     *
     * Currently, only rectangular boxes are supported.  This means that a, b, and c must be aligned with the
     * x, y, and z axes respectively.  Future releases may support arbitrary triclinic boxes.
     *
     * @param a      the vector defining the first edge of the periodic box
     * @param b      the vector defining the second edge of the periodic box
     * @param c      the vector defining the third edge of the periodic box
     */
    void getPeriodicBoxVectors(Vec3& a, Vec3& b, Vec3& c);
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
     * Update the positions of particles so that all distance constraints are satisfied.
     *
     * @param tol    the distance tolerance within which constraints must be satisfied.
     */
    void applyConstraints(double tol);
    /**
     * Recalculate all of the forces in the system and/or the potential energy of the system (in kJ/mol).
     * After calling this, use getForces() to retrieve the forces that were calculated.
     *
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy of the system, or 0 if includeEnergy is false
     */
    double calcForcesAndEnergy(bool includeForces, bool includeEnergy);
    /**
     * Calculate the kinetic energy of the system (in kJ/mol).
     */
    double calcKineticEnergy();
    /**
     * This should be called at the start of each time step.  It calls updateContextState() on each
     * ForceImpl in the system, allowing them to modify the values of state variables.
     */
    void updateContextState();
    /**
     * Get the list of ForceImpls belonging to this ContextImpl.
     */
    const std::vector<ForceImpl*>& getForceImpls() const;
    /**
     * Get the platform-specific data stored in this context.
     */
    void* getPlatformData();
    /**
     * Get the platform-specific data stored in this context.
     */
    const void* getPlatformData() const;
    /**
     * Set the platform-specific data stored in this context.
     */
    void setPlatformData(void* data);
    /**
     * Get a list of the particles in each molecules in the system.  Two particles are in the
     * same molecule if they are connected by constraints or bonds.
     */
    const std::vector<std::vector<int> >& getMolecules() const;
private:
    friend class Context;
    static void tagParticlesInMolecule(int particle, int molecule, std::vector<int>& particleMolecule, std::vector<std::vector<int> >& particleBonds);
    Context& owner;
    System& system;
    Integrator& integrator;
    std::vector<ForceImpl*> forceImpls;
    std::map<std::string, double> parameters;
    mutable std::vector<std::vector<int> > molecules;
    bool hasInitializedForces;
    Platform* platform;
    Kernel initializeForcesKernel, kineticEnergyKernel, updateStateDataKernel, applyConstraintsKernel;
    void* platformData;
};

} // namespace OpenMM

#endif /*OPENMM_CONTEXTIMPL_H_*/
