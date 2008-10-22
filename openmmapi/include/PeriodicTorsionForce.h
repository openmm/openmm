#ifndef OPENMM_PERIODICTORSIONFORCE_H_
#define OPENMM_PERIODICTORSIONFORCE_H_

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
 * This class implements an interaction between groups of four atoms that varies periodically with the torsion angle
 * between them.  When creating a PeriodicTorsionForce, you specify the number of torsions as an argument to the
 * constructor, then loop over them and call setTorsionParameters() to set the force field parameters for each one.
 */

class OPENMM_EXPORT PeriodicTorsionForce : public Force {
public:
    /**
     * Create a PeriodicTorsionForce.
     * 
     * @param numTorsions the number of periodic torsion terms in the potential function
     */
    PeriodicTorsionForce(int numTorsions);
    /**
     * Get the number of periodic torsion terms in the potential function
     */
    int getNumTorsions() const {
        return periodicTorsions.size();
    }
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
    void getTorsionParameters(int index, int& atom1, int& atom2, int& atom3, int& atom4, int& periodicity, double& phase, double& k) const;
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
    void setTorsionParameters(int index, int atom1, int atom2, int atom3, int atom4, int periodicity, double phase, double k);
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
protected:
    ForceImpl* createImpl();
private:
    class PeriodicTorsionInfo;

// Retarded visual studio compiler complains about being unable to 
// export private stl class members.
// This stanza explains that it should temporarily shut up.
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4251)
#endif
    std::vector<PeriodicTorsionInfo> periodicTorsions;
#if defined(_MSC_VER)
#pragma warning(pop)
#endif

};

class PeriodicTorsionForce::PeriodicTorsionInfo {
public:
    int atom1, atom2, atom3, atom4, periodicity;
    double phase, k;
    PeriodicTorsionInfo() {
        atom1 = atom2 = atom3 = atom4 = -1;
        periodicity = 1;
        phase = k = 0.0;
    }
};

} // namespace OpenMM

#endif /*OPENMM_PERIODICTORSIONFORCE_H_*/
