#ifndef OPENMM_AMOEBA_PI_TORSION_FORCE_H_
#define OPENMM_AMOEBA_PI_TORSION_FORCE_H_

/* -------------------------------------------------------------------------- *
 *                              AmoebaOpenMM                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
 * Authors:                                                                   *
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

#include "openmm/Force.h"
#include <vector>
#include "openmm/internal/windowsExport.h"

namespace OpenMM {

/**
 * This class implements the Amoeba pi-torsion interaction
 * To use it, create a PiTorsionForce object then call addPiTorsion() once for each torsion.  After
 * a torsion has been added, you can modify its force field parameters by calling setPiTorsionParameters().
 */

 class OPENMM_EXPORT AmoebaPiTorsionForce : public Force {
public:
    /**
     * Create a Amoeba PiTorsionForce.
     */
    AmoebaPiTorsionForce();

    /**
     * Get the number of pi torsion terms in the potential function
     */
    int getNumPiTorsions() const {
        return piTorsions.size();
    }

    /**
     * Add a torsion term to the force field.
     *
     * @param particle1     the index of the first particle connected by the torsion
     * @param particle2     the index of the second particle connected by the torsion
     * @param particle3     the index of the third particle connected by the torsion
     * @param particle4     the index of the fourth particle connected by the torsion
     * @param particle5     the index of the fifth particle connected by the torsion
     * @param particle6     the index of the sixth particle connected by the torsion
     * @param k             the force constant for the torsion
     * @return the index of the torsion that was added
     */
    int addPiTorsion(int particle1, int particle2, int particle3, int particle4, int particle5, int particle6, double k );

    /**
     * Get the force field parameters for a torsion term.
     * 
     * @param index         the index of the torsion for which to get parameters
     * @param particle1     the index of the first particle connected by the torsion
     * @param particle2     the index of the second particle connected by the torsion
     * @param particle3     the index of the third particle connected by the torsion
     * @param particle4     the index of the fourth particle connected by the torsion
     * @param particle5     the index of the fifth particle connected by the torsion
     * @param particle6     the index of the sixth particle connected by the torsion
     * @param k             the force constant for the torsion
     */
    void getPiTorsionParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4, int& particle5, int& particle6, double& k ) const;

    /**
     * Set the force field parameters for a pi torsion term.
     * 
     * @param index         the index of the torsion for which to set parameters
     * @param particle1     the index of the first particle connected by the torsion
     * @param particle2     the index of the second particle connected by the torsion
     * @param particle3     the index of the third particle connected by the torsion
     * @param particle4     the index of the fourth particle connected by the torsion
     * @param particle5     the index of the fifth particle connected by the torsion
     * @param particle6     the index of the sixth particle connected by the torsion
     * @param k             the force constant for the torsion
     */
    void setPiTorsionParameters(int index, int particle1, int particle2, int particle3, int particle4, int particle5, int particle6, double k );

protected:
    ForceImpl* createImpl();
private:

    class PiTorsionInfo;

// Retarded visual studio compiler complains about being unable to 
// export private stl class members.
// This stanza explains that it should temporarily shut up.
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4251)
#endif
    std::vector<PiTorsionInfo> piTorsions;
#if defined(_MSC_VER)
#pragma warning(pop)
#endif

};

class AmoebaPiTorsionForce::PiTorsionInfo {
public:
    int particle1, particle2, particle3, particle4, particle5, particle6;
    double k;
    PiTorsionInfo() {
        particle1 = particle2  = particle3 = particle4 = particle5 = particle6 = -1;
        k = 0.0;
    }
    PiTorsionInfo(int particle1, int particle2, int particle3, int particle4, int particle5, int particle6, double k ) :
        particle1(particle1), particle2(particle2), particle3(particle3), particle4(particle4), particle5(particle5), particle6(particle6), k(k) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_AMOEBA_PI_TORSION_FORCE_H_*/
