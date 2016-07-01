#ifndef OPENMM_SYSTEM_H_
#define OPENMM_SYSTEM_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
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

#include "Vec3.h"
#include <vector>
#include "internal/windowsExport.h"

namespace OpenMM {

class OPENMM_EXPORT Force;
class OPENMM_EXPORT VirtualSite;

/**
 * This class represents a molecular system.  The definition of a System involves
 * four elements:
 *
 * <ol>
 * <li>The set of particles in the system</li>
 * <li>The forces acting on them</li>
 * <li>Pairs of particles whose separation should be constrained to a fixed value</li>
 * <li>For periodic systems, the dimensions of the periodic box</li>
 * </ol>
 *
 * The particles and constraints are defined directly by the System object, while
 * forces are defined by objects that extend the Force class.  After creating a
 * System, call addParticle() once for each particle, addConstraint() for each constraint,
 * and addForce() for each Force.
 *
 * In addition, particles may be designated as "virtual sites".  These are particles
 * whose positions are computed automatically based on the positions of other particles.
 * To define a virtual site, call setVirtualSite(), passing in a VirtualSite object
 * that defines the rules for computing its position.
 */

class OPENMM_EXPORT System {
public:
    /**
     * Create a new System.
     */
    System();
    ~System();
    /**
     * Get the number of particles in this System.
     */
    int getNumParticles() const {
        return masses.size();
    }
    /**
     * Add a particle to the System.  If the mass is 0, Integrators will ignore
     * the particle and not modify its position or velocity.  This is most often
     * used for virtual sites, but can also be used as a way to prevent a particle
     * from moving.
     *
     * @param mass   the mass of the particle (in atomic mass units)
     * @return the index of the particle that was added
     */
    int addParticle(double mass) {
        masses.push_back(mass);
        return masses.size()-1;
    }
    /**
     * Get the mass (in atomic mass units) of a particle.  If the mass is 0, Integrators will ignore
     * the particle and not modify its position or velocity.  This is most often
     * used for virtual sites, but can also be used as a way to prevent a particle
     * from moving.
     *
     * @param index the index of the particle for which to get the mass
     */
    double getParticleMass(int index) const;
    /**
     * Set the mass (in atomic mass units) of a particle.  If the mass is 0, Integrators will ignore
     * the particle and not modify its position or velocity.  This is most often
     * used for virtual sites, but can also be used as a way to prevent a particle
     * from moving.
     *
     * @param index  the index of the particle for which to set the mass
     * @param mass   the mass of the particle
     */
    void setParticleMass(int index, double mass);
    /**
     * Set a particle to be a virtual site.  The VirtualSite object should have
     * been created on the heap with the "new" operator.  The System takes over
     * ownership of it, and deletes it when the System itself is deleted.
     *
     * @param index        the index of the particle that should be treated as a
     *                     virtual site
     * @param virtualSite  a pointer to the VirtualSite object describing it
     */
    void setVirtualSite(int index, VirtualSite* virtualSite);
    /**
     * Get whether a particle is a VirtualSite.
     *
     * @param index  the index of the particle to check
     */
    bool isVirtualSite(int index) const {
        return (index < (int) virtualSites.size() && virtualSites[index] != NULL);
    }
    /**
     * Get VirtualSite object for a particle.  If the particle is not a virtual
     * site, this throws an exception.
     *
     * @param index  the index of the particle to get
     */
    const VirtualSite& getVirtualSite(int index) const;
    /**
     * Get the number of distance constraints in this System.
     */
    int getNumConstraints() const {
        return constraints.size();
    }
    /**
     * Add a constraint to the System.  Particles whose mass is 0 cannot participate
     * in constraints.
     *
     * @param particle1 the index of the first particle involved in the constraint
     * @param particle2 the index of the second particle involved in the constraint
     * @param distance  the required distance between the two particles, measured in nm
     * @return the index of the constraint that was added
     */
    int addConstraint(int particle1, int particle2, double distance);
    /**
     * Get the parameters defining a distance constraint.
     *
     * @param index     the index of the constraint for which to get parameters
     * @param[out] particle1 the index of the first particle involved in the constraint
     * @param[out] particle2 the index of the second particle involved in the constraint
     * @param[out] distance  the required distance between the two particles, measured in nm
     */
    void getConstraintParameters(int index, int& particle1, int& particle2, double& distance) const;
    /**
     * Set the parameters defining a distance constraint.  Particles whose mass is 0 cannot participate
     * in constraints.
     *
     * @param index     the index of the constraint for which to set parameters
     * @param particle1 the index of the first particle involved in the constraint
     * @param particle2 the index of the second particle involved in the constraint
     * @param distance  the required distance between the two particles, measured in nm
     */
    void setConstraintParameters(int index, int particle1, int particle2, double distance);
    /**
     * Remove a constraint from the System.
     *
     * @param index    the index of the constraint to remove
     */
    void removeConstraint(int index);
    /**
     * Add a Force to the System.  The Force should have been created on the heap with the
     * "new" operator.  The System takes over ownership of it, and deletes the Force when the
     * System itself is deleted.
     *
     * @param force   a pointer to the Force object to be added
     * @return        the index within the System of the Force that was added
     */
    int addForce(Force* force) {
        forces.push_back(force);
        return forces.size()-1;
    }
    /**
     * Get the number of Force objects that have been added to the System.
     */
    int getNumForces() const {
        return forces.size();
    }
    /**
     * Get a const reference to one of the Forces in this System.
     *
     * @param index  the index of the Force to get
     */
    const Force& getForce(int index) const;
    /**
     * Get a writable reference to one of the Forces in this System.
     *
     * @param index  the index of the Force to get
     */
    Force& getForce(int index);
    /**
     * Remove a Force from the System.  The memory associated with the removed Force
     * object is deleted.
     *
     * @param index   the index of the Force to remove
     */
    void removeForce(int index);
    /**
     * Get the default values of the vectors defining the axes of the periodic box (measured in nm).  Any newly
     * created Context will have its box vectors set to these.  They will affect
     * any Force added to the System that uses periodic boundary conditions.
     *
     * @param[out] a     the vector defining the first edge of the periodic box
     * @param[out] b     the vector defining the second edge of the periodic box
     * @param[out] c     the vector defining the third edge of the periodic box
     */
    void getDefaultPeriodicBoxVectors(Vec3& a, Vec3& b, Vec3& c) const;
    /**
     * Set the default values of the vectors defining the axes of the periodic box (measured in nm).  Any newly
     * created Context will have its box vectors set to these.  They will affect
     * any Force added to the System that uses periodic boundary conditions.
     *
     * Triclinic boxes are supported, but the vectors must satisfy certain requirements.  In particular,
     * a must point in the x direction, b must point "mostly" in the y direction, and c must point "mostly"
     * in the z direction.  See the documentation for details.
     *
     * @param a      the vector defining the first edge of the periodic box
     * @param b      the vector defining the second edge of the periodic box
     * @param c      the vector defining the third edge of the periodic box
     */
    void setDefaultPeriodicBoxVectors(const Vec3& a, const Vec3& b, const Vec3& c);
    /**
     * Returns whether or not any forces in this System use periodic boundaries.
     *
     * If a force in this System does not implement usesPeriodicBoundaryConditions
     * a OpenMM::OpenMMException is thrown
     *
     * @return true if at least one force uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const;
private:
    class ConstraintInfo;
    Vec3 periodicBoxVectors[3];
    std::vector<double> masses;
    std::vector<ConstraintInfo> constraints;
    std::vector<Force*> forces;
    std::vector<VirtualSite*> virtualSites;
};

/**
 * This is an internal class used to record information about a constraint.
 * @private
 */
class System::ConstraintInfo {
public:
    int particle1, particle2;
    double distance;
    ConstraintInfo() {
        particle1 = particle2 = -1;
        distance = 0.0;
    }
    ConstraintInfo(int particle1, int particle2, double distance) :
        particle1(particle1), particle2(particle2), distance(distance) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_SYSTEM_H_*/
