#ifndef OPENMM_NONBONDEDFORCE_H_
#define OPENMM_NONBONDEDFORCE_H_

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

#include "Force.h"
#include "Vec3.h"
#include <map>
#include <vector>
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This class implements nonbonded interactions between particles, including a Coulomb force to represent
 * electrostatics and a Lennard-Jones force to represent van der Waals interactions.  It optionally supports
 * periodic boundary conditions and cutoffs for long range interactions.
 * 
 * If the System also contains a HarmonicBondForce, nonbonded interactions are automatically excluded between
 * particles which are separated by three or fewer bonds.  Most molecular force fields omit Coulomb and
 * Lennard-Jones interactions between particles separated by one or two bonds, while using modified parameters
 * for those separated by three bonds (known as "1-4 interactions").  This class lets you provide a list of
 * 1-4 interactions to include in the potential, along with the parameters to use for each one.
 * 
 * When creating a NonbondedForce, you specify the number of atoms and the number of 1-4 interactions as
 * arguments to the constructor.  You then loop over them and call setAtomParameters() for each atom and
 * setNonbond14Parameters() for each 1-4 interaction.
 */

class OPENMM_EXPORT NonbondedForce : public Force {
public:
    /**
     * This is an enumeration of the different methods that may be used for handling long range nonbonded forces.
     */
    enum NonbondedMethod {
        /**
         * No cutoff is applied to nonbonded interactions.  The full set of N^2 interactions is computed exactly.
         * This necessarily means that periodic boundary conditions cannot be used.  This is the default.
         */
        NoCutoff = 0,
        /**
         * Interactions beyond the cutoff distance are ignored.  Coulomb interactions closer than the cutoff distance
         * are modified using the reaction field method.
         */
        CutoffNonPeriodic = 1,
        /**
         * Periodic boundary conditions are used, so that each atom interacts only with the nearest periodic copy of
         * each other atom.  Interactions beyond the cutoff distance are ignored.  Coulomb interactions closer than the
         * cutoff distance are modified using the reaction field method.
         */
        CutoffPeriodic = 2
    };
    /**
     * Create a NonbondedForce.
     * 
     * @param numAtoms            the number of atoms in the system
     * @param numNonbonded14      the number of nonbonded 1-4 terms in the potential function
     */
    NonbondedForce(int numAtoms, int numNonbonded14);
    /**
     * Get the number of atoms in the system.
     */
    int getNumAtoms() const {
        return atoms.size();
    }
    /**
     * Get the number of nonbonded 1-4 terms in the potential function
     */
    int getNumNonbonded14() const {
        return nb14s.size();
    }
    /**
     * Get the method used for handling long range nonbonded interactions.
     */
    NonbondedMethod getNonbondedMethod() const;
    /**
     * Set the method used for handling long range nonbonded interactions.
     */
    void setNonbondedMethod(NonbondedMethod method);
    /**
     * Get the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
     * does not use cutoffs, this value will have no effect.
     */
    double getCutoffDistance() const;
    /**
     * Set the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
     * does not use cutoffs, this value will have no effect.
     */
    void setCutoffDistance(double distance);
    /**
     * Get the vectors which define the axes of the periodic box (measured in nm).  If the NonbondedMethod
     * in use does not use periodic boundary conditions, these values will have no effect.
     *
     * Currently, only rectangular boxes are supported.  This means that a, b, and c must be aligned with the
     * x, y, and z axes respectively.  Future releases may support arbitrary triclinic boxes.
     *
     * @param a      on exit, this contains the vector defining the first edge of the periodic box
     * @param b      on exit, this contains the vector defining the second edge of the periodic box
     * @param c      on exit, this contains the vector defining the third edge of the periodic box
     */
    void getPeriodicBoxVectors(Vec3& a, Vec3& b, Vec3& c) const;
    /**
     * Set the vectors which define the axes of the periodic box (measured in nm).  If the NonbondedMethod
     * in use does not use periodic boundary conditions, these values will have no effect.
     *
     * Currently, only rectangular boxes are supported.  This means that a, b, and c must be aligned with the
     * x, y, and z axes respectively.  Future releases may support arbitrary triclinic boxes.
     *
     * @param a      the vector defining the first edge of the periodic box
     * @param b      the vector defining the second edge of the periodic box
     * @param c      the vector defining the third edge of the periodic box
     */
    void setPeriodicBoxVectors(Vec3 a, Vec3 b, Vec3 c);
    /**
     * Get the nonbonded force parameters for an atom.
     * 
     * @param index     the index of the atom for which to get parameters
     * @param charge    the charge of the atom, measured in units of the proton charge
     * @param radius    the van der Waals radius of the atom, measured in nm
     * @param depth     the well depth of the van der Waals interaction, measured in kJ/mol
     */
    void getAtomParameters(int index, double& charge, double& radius, double& depth) const;
    /**
     * Set the nonbonded force parameters for an atom.
     * 
     * @param index     the index of the atom for which to set parameters
     * @param charge    the charge of the atom, measured in units of the proton charge
     * @param radius    the van der Waals radius of the atom (sigma in the Lennard Jones potential), measured in nm
     * @param depth     the well depth of the van der Waals interaction (epsilon in the Lennard Jones potential), measured in kJ/mol
     */
    void setAtomParameters(int index, double charge, double radius, double depth);
    /**
     * Get the force field parameters for a nonbonded 1-4 term.
     * 
     * @param index     the index of the interaction for which to get parameters
     * @param atom1     the index of the first atom involved in the interaction
     * @param atom2     the index of the second atom involved in the interaction
     * @param charge    the scaled product of the atomic charges (i.e. the strength of the Coulomb interaction), measured in units of the proton charge
     * @param radius    the van der Waals radius of the atom (sigma in the Lennard Jones potential), measured in nm
     * @param depth     the well depth of the van der Waals interaction (epsilon in the Lennard Jones potential), measured in kJ/mol
     */
    void getNonbonded14Parameters(int index, int& atom1, int& atom2, double& charge, double& radius, double& depth) const;
    /**
     * Set the force field parameters for a nonbonded 1-4 term.
     * 
     * @param index     the index of the interaction for which to get parameters
     * @param atom1     the index of the first atom involved in the interaction
     * @param atom2     the index of the second atom involved in the interaction
     * @param charge    the scaled product of the atomic charges (i.e. the strength of the Coulomb interaction), measured in units of the proton charge
     * @param radius    the van der Waals radius of the atom (sigma in the Lennard Jones potential), measured in nm
     * @param depth     the well depth of the van der Waals interaction (epsilon in the Lennard Jones potential), measured in kJ/mol
     */
    void setNonbonded14Parameters(int index, int atom1, int atom2, double charge, double radius, double depth);
protected:
    ForceImpl* createImpl();
private:
    class AtomInfo;
    class NB14Info;
    NonbondedMethod nonbondedMethod;
    double cutoffDistance;
    Vec3 periodicBoxVectors[3];

// Retarded visual studio compiler complains about being unable to 
// export private stl class members.
// This stanza explains that it should temporarily shut up.
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4251)
#endif
    std::vector<AtomInfo> atoms;
    std::vector<NB14Info> nb14s;
#if defined(_MSC_VER)
#pragma warning(pop)
#endif

};

class NonbondedForce::AtomInfo {
public:
    double charge, radius, depth;
    AtomInfo() {
        charge = radius = depth = 0.0;
    }
};

class NonbondedForce::NB14Info {
public:
    int atom1, atom2;
    double charge, radius, depth;
    NB14Info() {
        atom1 = atom2 = -1;
        charge = radius = depth = 0.0;
    }
};

} // namespace OpenMM

#endif /*OPENMM_NONBONDEDFORCE_H_*/
