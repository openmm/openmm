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
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
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
#include <iosfwd>
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
    ContextImpl(Context& owner, const System& system, Integrator& integrator, Platform* platform, const std::map<std::string, std::string>& properties,
            ContextImpl* originalContext=NULL);
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
    const System& getSystem() const {
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
     * Get the set of all adjustable parameters and their values
     */
    const std::map<std::string, double>& getParameters() const;
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
     * Get the derivatives of the energy with respect to parameters.
     */
    void getEnergyParameterDerivatives(std::map<std::string, double>& derivs);
    /**
     * Get the vectors defining the axes of the periodic box (measured in nm).  They will affect
     * any Force that uses periodic boundary conditions.
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
     * Recalculate all of the forces in the system and/or the potential energy of the system (in kJ/mol).
     * After calling this, use getForces() to retrieve the forces that were calculated.
     *
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @param groups         a set of bit flags for which force groups to include.  Group i will be included
     *                       if (groups&(1<<i)) != 0.  The default value includes all groups.
     * @return the potential energy of the system, or 0 if includeEnergy is false
     */
    double calcForcesAndEnergy(bool includeForces, bool includeEnergy, int groups=0xFFFFFFFF);
    /**
     * Get the set of force group flags that were passed to the most recent call to calcForcesAndEnergy().
     * 
     * Note that this returns a reference, so it's possible to modify it.  Be very very cautious about
     * doing that!  Only do it if you're also modifying forces stored inside the context.
     */
    int& getLastForceGroups();
    /**
     * Calculate the kinetic energy of the system (in kJ/mol).
     */
    double calcKineticEnergy();
    /**
     * This should be called at the start of each time step.  It calls updateContextState() on each
     * ForceImpl in the system, allowing them to modify the values of state variables.
     * 
     * @return true if the state was modified in any way that would cause the forces on particles
     * to change, false otherwise
     */
    bool updateContextState();
    /**
     * Get the list of ForceImpls belonging to this ContextImpl.
     */
    const std::vector<ForceImpl*>& getForceImpls() const;
    /**
     * Get the list of ForceImpls belonging to this ContextImpl.
     */
    std::vector<ForceImpl*>& getForceImpls();
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
    /**
     * Create a checkpoint recording the current state of the Context.
     * 
     * @param stream    an output stream the checkpoint data should be written to
     */
    void createCheckpoint(std::ostream& stream);
    /**
     * Load a checkpoint that was written by createCheckpoint().
     * 
     * @param stream    an input stream the checkpoint data should be read from
     */
    void loadCheckpoint(std::istream& stream);
    /**
     * This is invoked by the Integrator when it is deleted.  This is needed to ensure the cleanup process
     * is done correctly, since we don't know whether the Integrator or Context will be deleted first.
     */
    void integratorDeleted() {
        integratorIsDeleted = true;
    }
    /**
     * Notify the integrator that some aspect of the system has changed, and cached information should be discarded.
     */
    void systemChanged();
    /**
     * This is the routine that actually computes the list of molecules returned by getMolecules().  Normally
     * you should never call it.  It is exposed here because the same logic is useful to other classes too.
     */
    static std::vector<std::vector<int> > findMolecules(int numParticles, std::vector<std::vector<int> >& particleBonds);
    /**
     * Create a new Context based on this one.  The new context will use the same Platform, device, and property
     * values as this one.  With the CUDA and OpenCL platforms, it also shares the same GPU context, allowing data
     * to be transferred between them without leaving the GPU.
     * 
     * This method exists for very specialized purposes.  If you aren't certain whether you should use it, that probably
     * means you shouldn't.
     */
    Context* createLinkedContext(const System& system, Integrator& integrator);
private:
    friend class Context;
    void initialize();
    Context& owner;
    const System& system;
    Integrator& integrator;
    std::vector<ForceImpl*> forceImpls;
    std::map<std::string, double> parameters;
    mutable std::vector<std::vector<int> > molecules;
    bool hasInitializedForces, hasSetPositions, integratorIsDeleted;
    int lastForceGroups;
    Platform* platform;
    Kernel initializeForcesKernel, updateStateDataKernel, applyConstraintsKernel, virtualSitesKernel;
    void* platformData;
};

} // namespace OpenMM

#endif /*OPENMM_CONTEXTIMPL_H_*/
