#ifndef OPENMM_STANDARDMMFORCEFIELD_H_
#define OPENMM_STANDARDMMFORCEFIELD_H_

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
#include <map>
#include <vector>
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This class can be used for a variety of standard molecular mechanics force fields.  It includes
 * terms for Coulomb and Lennard-Jones nonbonded interactions, and harmonic terms for the following
 * types of bonded interactions: bond length, bond angle, and torsion (both periodic and Ryckaert-Bellemans).
 * 
 * To create a StandardMMForceField, you specify the number of atoms and the number of each type
 * of bonded term to the constructor.  Then loop over them and call the appropriate setXXXParameter()
 * method to set the force field parameters for each one.
 */

class OPENMM_EXPORT StandardMMForceField : public Force {
public:
    /**
     * Create a StandardMMForceField.
     * 
     * @param numAtoms            the number of atoms in the system
     * @param numBonds            the number of harmonic bond stretch terms in the potential function
     * @param numAngles           the number of harmonic bond angle terms in the potential function
     * @param numPeriodicTorsions the number of periodic torsion terms in the potential function
     * @param numRBTorsions       the number of Ryckaert-Bellemans torsion terms in the potential function
     */
    StandardMMForceField(int numAtoms, int numBonds, int numAngles, int numPeriodicTorsions, int numRBTorsions);
    /**
     * Get the number of atoms in the system.
     */
    int getNumAtoms() const {
        return atoms.size();
    }
    /**
     * Get the number of harmonic bond stretch terms in the potential function
     */
    int getNumBonds() const {
        return bonds.size();
    }
    /**
     * Get the number of harmonic bond angle terms in the potential function
     */
    int getNumAngles() const {
        return angles.size();
    }
    /**
     * Get the number of periodic torsion terms in the potential function
     */
    int getNumPeriodicTorsions() const {
        return periodicTorsions.size();
    }
    /**
     * Get the number of Ryckaert-Bellemans torsion terms in the potential function
     */
    int getNumRBTorsions() const {
        return rbTorsions.size();
    }
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
     * Get the force field parameters for a bond term.
     * 
     * @param index     the index of the bond for which to get parameters
     * @param atom1     the index of the first atom connected by the bond
     * @param atom2     the index of the second atom connected by the bond
     * @param length    the equilibrium length of the bond, measured in nm
     * @param k         the harmonic force constant for the bond
     */
    void getBondParameters(int index, int& atom1, int& atom2, double& length, double& k) const;
    /**
     * Set the force field parameters for a bond term.
     * 
     * @param index     the index of the bond for which to set parameters
     * @param atom1     the index of the first atom connected by the bond
     * @param atom2     the index of the second atom connected by the bond
     * @param length    the equilibrium length of the bond, measured in nm
     * @param k         the harmonic force constant for the bond
     */
    void setBondParameters(int index, int atom1, int atom2, double length, double k);
    /**
     * Get the force field parameters for an angle term.
     * 
     * @param index     the index of the angle for which to get parameters
     * @param atom1     the index of the first atom forming the angle
     * @param atom2     the index of the second atom forming the angle
     * @param atom3     the index of the third atom forming the angle
     * @param length    the equilibrium angle, measured in radians
     * @param k         the harmonic force constant for the angle
     */
    void getAngleParameters(int index, int& atom1, int& atom2, int& atom3, double& angle, double& k) const;
    /**
     * Set the force field parameters for an angle term.
     * 
     * @param index     the index of the angle for which to set parameters
     * @param atom1     the index of the first atom forming the angle
     * @param atom2     the index of the second atom forming the angle
     * @param atom3     the index of the third atom forming the angle
     * @param length    the equilibrium angle, measured in radians
     * @param k         the harmonic force constant for the angle
     */
    void setAngleParameters(int index, int atom1, int atom2, int atom3, double angle, double k);
    /**
     * Get the force field parameters for a periodic torsion term.
     * 
     * @param index        the index of the torsion for which to get parameters
     * @param atom1        the index of the first atom forming the torsion
     * @param atom2        the index of the second atom forming the torsion
     * @param atom3        the index of the third atom forming the torsion
     * @param atom3        the index of the fourth atom forming the torsion
     * @param periodicity  the periodicity of the torsion
     * @param phase        the phase offset of the torsion, measured in radians
     * @param k            the force constant for the torsion
     */
    void getPeriodicTorsionParameters(int index, int& atom1, int& atom2, int& atom3, int& atom4, int& periodicity, double& phase, double& k) const;
    /**
     * Set the force field parameters for a periodic torsion term.
     * 
     * @param index        the index of the torsion for which to set parameters
     * @param atom1        the index of the first atom forming the torsion
     * @param atom2        the index of the second atom forming the torsion
     * @param atom3        the index of the third atom forming the torsion
     * @param atom3        the index of the fourth atom forming the torsion
     * @param periodicity  the periodicity of the torsion
     * @param phase        the phase offset of the torsion, measured in radians
     * @param k            the force constant for the torsion
     */
    void setPeriodicTorsionParameters(int index, int atom1, int atom2, int atom3, int atom4, int periodicity, double phase, double k);
    /**
     * Get the force field parameters for a Ryckaert-Bellemans torsion term.
     * 
     * @param index        the index of the torsion for which to get parameters
     * @param atom1        the index of the first atom forming the torsion
     * @param atom2        the index of the second atom forming the torsion
     * @param atom3        the index of the third atom forming the torsion
     * @param atom3        the index of the fourth atom forming the torsion
     * @param c0           the coefficient of the constant term
     * @param c1           the coefficient of the 1st order term
     * @param c2           the coefficient of the 2nd order term
     * @param c3           the coefficient of the 3rd order term
     * @param c4           the coefficient of the 4th order term
     * @param c5           the coefficient of the 5th order term
     */
    void getRBTorsionParameters(int index, int& atom1, int& atom2, int& atom3, int& atom4, double& c0, double& c1, double& c2, double& c3, double& c4, double& c5) const;
    /**
     * Set the force field parameters for a Ryckaert-Bellemans torsion term.
     * 
     * @param index        the index of the torsion for which to set parameters
     * @param atom1        the index of the first atom forming the torsion
     * @param atom2        the index of the second atom forming the torsion
     * @param atom3        the index of the third atom forming the torsion
     * @param atom3        the index of the fourth atom forming the torsion
     * @param c0           the coefficient of the constant term
     * @param c1           the coefficient of the 1st order term
     * @param c2           the coefficient of the 2nd order term
     * @param c3           the coefficient of the 3rd order term
     * @param c4           the coefficient of the 4th order term
     * @param c5           the coefficient of the 5th order term
     */
    void setRBTorsionParameters(int index, int atom1, int atom2, int atom3, int atom4, double c0, double c1, double c2, double c3, double c4, double c5);
protected:
    ForceImpl* createImpl();
private:
    class AtomInfo;
    class BondInfo;
    class AngleInfo;
    class PeriodicTorsionInfo;
    class RBTorsionInfo;

// Retarded visual studio compiler complains about being unable to 
// export private stl class members.
// This stanza explains that it should temporarily shut up.
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4251)
#endif

    std::vector<AtomInfo> atoms;
    std::vector<BondInfo> bonds;
    std::vector<AngleInfo> angles;
    std::vector<PeriodicTorsionInfo> periodicTorsions;
    std::vector<RBTorsionInfo> rbTorsions;

#if defined(_MSC_VER)
#pragma warning(pop)
#endif

};

class StandardMMForceField::AtomInfo {
public:
    double charge, radius, depth;
    AtomInfo() {
        charge = radius = depth = 0.0;
    }
};

class StandardMMForceField::AngleInfo {
public:
    int atom1, atom2, atom3;
    double angle, k;
    AngleInfo() {
        atom1 = atom2 = atom3 = -1;
        angle = k = 0.0;
    }
};

class StandardMMForceField::PeriodicTorsionInfo {
public:
    int atom1, atom2, atom3, atom4, periodicity;
    double phase, k;
    PeriodicTorsionInfo() {
        atom1 = atom2 = atom3 = atom4 = -1;
        periodicity = 1;
        phase = k = 0.0;
    }
};

class StandardMMForceField::RBTorsionInfo {
public:
    int atom1, atom2, atom3, atom4;
    double c[6];
    RBTorsionInfo() {
        atom1 = atom2 = atom3 = atom4 = -1;
        c[0] = c[1] = c[2] = c[3] = c[4] = c[5] = 0.0;
    }
};

class StandardMMForceField::BondInfo {
public:
    int atom1, atom2;
    double length, k;
    BondInfo() {
        atom1 = atom2 = -1;
        length = k = 0.0;
    }
};

} // namespace OpenMM

#endif /*OPENMM_STANDARDMMFORCEFIELD_H_*/
