#ifndef OPENMM_HARMONICBONDFORCE_H_
#define OPENMM_HARMONICBONDFORCE_H_

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
 * This class implements an interaction between pairs of atoms that varies harmonically with the distance
 * between them.  When creating a HarmonicBondForce, you specify the number of bonds as an argument to the
 * constructor, then loop over them and call setBondParameters() to set the force field parameters for each one.
 */

class OPENMM_EXPORT HarmonicBondForce : public Force {
public:
    /**
     * Create a HarmonicBondForce.
     * 
     * @param numBonds            the number of harmonic bond stretch terms in the potential function
     */
    HarmonicBondForce(int numBonds);
    /**
     * Get the number of harmonic bond stretch terms in the potential function
     */
    int getNumBonds() const {
        return bonds.size();
    }
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
protected:
    ForceImpl* createImpl();
private:
    class BondInfo;

// Retarded visual studio compiler complains about being unable to 
// export private stl class members.
// This stanza explains that it should temporarily shut up.
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4251)
#endif
    std::vector<BondInfo> bonds;
#if defined(_MSC_VER)
#pragma warning(pop)
#endif

};

class HarmonicBondForce::BondInfo {
public:
    int atom1, atom2;
    double length, k;
    BondInfo() {
        atom1 = atom2 = -1;
        length = k = 0.0;
    }
};

} // namespace OpenMM

#endif /*OPENMM_HARMONICBONDFORCE_H_*/
