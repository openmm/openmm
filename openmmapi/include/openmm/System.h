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

#include <vector>
#include "internal/windowsExport.h"

namespace OpenMM {

class OPENMM_EXPORT Force;

/**
 * This class represents a molecular system.  The definition of a System involves
 * three elements:
 * 
 * <ol>
 * <li>The set of particles in the system</li>
 * <li>The forces acting on them</li>
 * <li>Pairs of particles whose separation should be constrained to a fixed value</li>
 * </ol>
 * 
 * The particles and constraints are defined directly by the System object, while
 * forces are defined by objects that extend the Force class.  After creating a
 * System, call addParticle() once for each particle, addConstraint() for each constraint,
 * and addForce() for each Force.
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
     * Add a particle to the System.
     *
     * @param mass   the mass of the particle (in atomic mass units)
     * @return the index of the particle that was added
     */
    int addParticle(double mass) {
        masses.push_back(mass);
        return masses.size()-1;
    }
    /**
     * Get the mass (in atomic mass units) of a particle.
     * 
     * @param index the index of the particle for which to get the mass
     */
    double getParticleMass(int index) const {
        return masses[index];
    }
    /**
     * Set the mass (in atomic mass units) of a particle.
     * 
     * @param index  the index of the particle for which to set the mass
     * @param mass   the mass of the particle
     */
    void setParticleMass(int index, double mass) {
        masses[index] = mass;
    }
    /**
     * Get the number of distance constraints in this System.
     */
    int getNumConstraints() const {
        return constraints.size();
    }
    /**
     * Add a constraint to the System.
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
     * @param particle1 the index of the first particle involved in the constraint
     * @param particle2 the index of the second particle involved in the constraint
     * @param distance  the required distance between the two particles, measured in nm
     */
    void getConstraintParameters(int index, int& particle1, int& particle2, double& distance) const;
    /**
     * Set the parameters defining a distance constraint.
     * 
     * @param index     the index of the constraint for which to set parameters
     * @param particle1 the index of the first particle involved in the constraint
     * @param particle2 the index of the second particle involved in the constraint
     * @param distance  the required distance between the two particles, measured in nm
     */
    void setConstraintParameters(int index, int particle1, int particle2, double distance);
    /**
     * Add a Force to the System.  The Force should have been created on the heap with the
     * "new" operator.  The System takes over ownership of it, and deletes the Force when the
     * System itself is deleted.
     */
    void addForce(Force* force) {
        forces.push_back(force);
    }
    /**
     * Get the number of Force objects that have been added to the System.
     */
    int getNumForces() const {
        return forces.size();
    }
    /**
     * Get a reference to one of the Forces in this System.
     * 
     * @param index  the index of the Force to get
     */
    const Force& getForce(int index) const {
        return *forces[index];
    }
    /**
     * Get a reference to one of the Forces in this System.
     *
     * @param index  the index of the Force to get
     */
    Force& getForce(int index) {
        return *forces[index];
    }
private:
    class ConstraintInfo;

// Retarded visual studio compiler complains about being unable to 
// export private stl class members.
// This stanza explains that it should temporarily shut up.
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4251)
#endif

    std::vector<double> masses;
    std::vector<ConstraintInfo> constraints;
    std::vector<Force*> forces;

#if defined(_MSC_VER)
#pragma warning(pop)
#endif

};

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
