#ifndef OPENMM_AMOEBA_STRETCH_BEND_FORCE_H_
#define OPENMM_AMOEBA_STRETCH_BEND_FORCE_H_

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
 * This class implements the Amoeba stretch bend interaction
 * To use it, create a StretchBendForce object then call addStretchBend() once for each stretchBend.  After
 * a stretchBend has been added, you can modify its force field parameters by calling setStretchBendParameters().
 */

 class OPENMM_EXPORT AmoebaStretchBendForce : public Force {
public:
    /**
     * Create a Amoeba StretchBendForce.
     */
    AmoebaStretchBendForce();

    /**
     * Get the number of stretchBend terms in the potential function
     */
    int getNumStretchBends() const {
        return stretchBends.size();
    }

    /**
     * Add a stretchBend term to the force field.
     *
     * @param particle1     the index of the first particle connected by the stretchBend
     * @param particle2     the index of the second particle connected by the stretchBend
     * @param particle3     the index of the third particle connected by the stretchBend
     * @param lengthAB      the equilibrium length of the stretchBend in bond ab [particle1, particle2], measured in nm
     * @param lengthCB      the equilibrium length of the stretchBend in bond cb [particle3, particle2], measured in nm
     * @param angle         the equilibrium angle in radians
     * @param k             the force constant for the stretchBend
     * @return the index of the stretchBend that was added
     */
    int addStretchBend(int particle1, int particle2, int particle3, double lengthAB,  double lengthCB, double angle,
                       double k );

    /**
     * Get the force field parameters for a stretchBend term.
     * 
     * @param index         the index of the stretchBend for which to get parameters
     * @param particle1     the index of the first particle connected by the stretchBend
     * @param particle2     the index of the second particle connected by the stretchBend
     * @param particle3     the index of the third particle connected by the stretchBend
     * @param lengthAB      the equilibrium length of the stretchBend in bond ab [particle1, particle2], measured in nm
     * @param lengthCB      the equilibrium length of the stretchBend in bond cb [particle3, particle2], measured in nm
     * @param angle         the equilibrium angle in radians
     * @param k             the force constant for the stretchBend
     */
    void getStretchBendParameters(int index, int& particle1, int& particle2, int& particle3,
                                  double& lengthAB, double& lengthCB, double& angle, double& k ) const;

    /**
     * Set the force field parameters for a stretchBend term.
     * 
     * @param index         the index of the stretchBend for which to set parameters
     * @param particle1     the index of the first particle connected by the stretchBend
     * @param particle2     the index of the second particle connected by the stretchBend
     * @param particle3     the index of the third particle connected by the stretchBend
     * @param lengthAB      the equilibrium length of the stretchBend in bond ab [particle1, particle2], measured in nm
     * @param lengthCB      the equilibrium length of the stretchBend in bond cb [particle3, particle2], measured in nm
     * @param angle         the equilibrium angle in radians
     * @param k             the force constant for the stretchBend
     */
    void setStretchBendParameters(int index, int particle1, int particle2, int particle3, 
                                  double lengthAB,  double lengthCB, double angle, double k );

protected:
    ForceImpl* createImpl();
private:

    class StretchBendInfo;

// Retarded visual studio compiler complains about being unable to 
// export private stl class members.
// This stanza explains that it should temporarily shut up.
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4251)
#endif
    std::vector<StretchBendInfo> stretchBends;
#if defined(_MSC_VER)
#pragma warning(pop)
#endif

};

class AmoebaStretchBendForce::StretchBendInfo {
public:
    int particle1, particle2, particle3;
    double lengthAB, lengthCB, angle, k;
    StretchBendInfo() {
        particle1 = particle2  = particle3 = -1;
        lengthAB  = lengthCB   = angle     = k   = 0.0;
    }
    StretchBendInfo(int particle1, int particle2, int particle3, 
                    double lengthAB,  double lengthCB, double angle, double k ) :
                    particle1(particle1), particle2(particle2), particle3(particle3), 
                    lengthAB(lengthAB), lengthCB(lengthCB), angle(angle), k(k) {
     
    }
};

} // namespace OpenMM

#endif /*OPENMM_AMOEBA_STRETCH_BEND_FORCE_H_*/
