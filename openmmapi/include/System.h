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

#include <vector>

namespace OpenMM {

class Force;

/**
 * This class represents a molecular system.  The definition of a System involves
 * three elements:
 * 
 * <ol>
 * <li>The set of atoms in the system</li>
 * <li>The forces acting on them</li>
 * <li>Pairs of atoms whose separation should be connstrained to a fixed value</li>
 * </ol>
 * 
 * The atoms and constraints are defined directly by the System object.
 * The forces are defined by objects that extend the Force class.  The System
 * stores a list of Force objects that determine the motion of the atoms.
 */

class System {
public:
    /**
     * Create a new System.
     * 
     * @param numAtoms       the number of atoms in the System
     * @param numConstraints the number of distance constraints in the System.
     */
    System(int numAtoms, int numConstraints);
    ~System();
    /**
     * Get the number of atoms in this System.
     */
    int getNumAtoms() const {
        return masses.size();
    }
    /**
     * Get the mass (in atomic mass units) of an atom.
     * 
     * @param index the index of the atom for which to get the mass
     */
    double getAtomMass(int index) const {
        return masses[index];
    }
    /**
     * Set the mass (in atomic mass units) of an atom.
     * 
     * @param index  the index of the atom for which to set the mass
     * @param mass   the mass of the atom
     */
    void setAtomMass(int index, double mass) {
        masses[index] = mass;
    }
    /**
     * Get the number of distance constraints in this System.
     */
    int getNumConstraints() const {
        return constraints.size();
    }
    /**
     * Get the parameters defining a distance constraint.
     * 
     * @param index     the index of the constraint for which to get parameters
     * @param atom1     the index of the first atom involved in the constraint
     * @param atom2     the index of the second atom involved in the constraint
     * @param distance  the required distance between the two atoms, measured in nm
     */
    void getConstraintParameters(int index, int& atom1, int& atom2, double& distance) const;
    /**
     * Set the parameters defining a distance constraint.
     * 
     * @param index     the index of the constraint for which to set parameters
     * @param atom1     the index of the first atom involved in the constraint
     * @param atom2     the index of the second atom involved in the constraint
     * @param distance  the required distance between the two atoms, measured in nm
     */
    void setConstraintParameters(int index, int atom1, int atom2, double distance);
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
    Force& getForce(int index) {
        return *forces[index];
    }
private:
    class ConstraintInfo;
    std::vector<int> masses;
    std::vector<ConstraintInfo> constraints;
    std::vector<Force*> forces;
};

class System::ConstraintInfo {
public:
    int atom1, atom2;
    double distance;
    ConstraintInfo() {
        atom1 = atom2 = -1;
        distance = 0.0;
    }
};

} // namespace OpenMM

#endif /*OPENMM_SYSTEM_H_*/
