#ifndef OPENMM_AMOEBA_OUT_OF_PLANE_BEND_FORCE_H_
#define OPENMM_AMOEBA_OUT_OF_PLANE_BEND_FORCE_H_

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
 * This class implements the Amoeba Out-of-plane bend interaction
 * To use it, create a OutOfPlaneBendForce object then call addOutOfPlaneBend() once for each outOfPlaneBend.  After
 * a outOfPlaneBend has been added, you can modify its force field parameters by calling setOutOfPlaneBendParameters().
 */

class OPENMM_EXPORT AmoebaOutOfPlaneBendForce : public Force {
public:
    /**
     * Create a Amoeba OutOfPlaneBendForce.
     */
    AmoebaOutOfPlaneBendForce();

    /**
     * Get the number of outOfPlaneBend terms in the potential function
     */
    int getNumOutOfPlaneBends() const {
        return outOfPlaneBends.size();
    }

    /** 
     * Set the global cubic term
     * 
     * @param cubicK        the cubic harmonic force constant for the angle
     */
    void setAmoebaGlobalOutOfPlaneBendCubic( double cubicK );

    /** 
     * Get the global cubic term
     * 
     * @return global cubicK term
     */
    double getAmoebaGlobalOutOfPlaneBendCubic( void ) const;

    /** 
     * Set the global cubic term
     * 
     * @param quarticK       the quartic harmonic force constant for the angle
     */
    void setAmoebaGlobalOutOfPlaneBendQuartic( double quarticK );

    /** 
     * Get the global quartic term
     * 
     * @return global  quartic term
     */
    double getAmoebaGlobalOutOfPlaneBendQuartic( void ) const;

    /** 
     * Set the global pentic term
     * 
     * @param penticK the pentic harmonic force constant for the angle
     */
    void setAmoebaGlobalOutOfPlaneBendPentic( double penticK );

    /** 
     * Get the global pentic term
     * 
     * @return global penticK term
     */
    double getAmoebaGlobalOutOfPlaneBendPentic( void ) const;

    /** 
     * Set the global sextic term
     * 
     * @param sexticK       the sextic harmonic force constant for the angle
     */
    void setAmoebaGlobalOutOfPlaneBendSextic( double sexticK );

    /** 
     * Get the global sextic term
     * 
     * @return global sexticK term
     */
    double getAmoebaGlobalOutOfPlaneBendSextic( void ) const;

    /**
     * Add a outOfPlaneBend term to the force field.
     *
     * @param particle1     the index of the first particle connected by the outOfPlaneBend
     * @param particle2     the index of the second particle connected by the outOfPlaneBend
     * @param particle3     the index of the third particle connected by the outOfPlaneBend
     * @param particle4     the index of the fourth particle connected by the outOfPlaneBend
     * @param k             the force constant for the outOfPlaneBend
     * @return the index of the outOfPlaneBend that was added
     */
    int addOutOfPlaneBend(int particle1, int particle2, int particle3, int particle4, double k );

    /**
     * Get the force field parameters for a outOfPlaneBend term.
     * 
     * @param index         the index of the outOfPlaneBend for which to get parameters
     * @param particle1     the index of the first particle connected by the outOfPlaneBend
     * @param particle2     the index of the second particle connected by the outOfPlaneBend
     * @param particle3     the index of the third particle connected by the outOfPlaneBend
     * @param particle4     the index of the fourth particle connected by the outOfPlaneBend
     * @param k             the force constant for the outOfPlaneBend
     */
    void getOutOfPlaneBendParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4, double& k ) const;

    /**
     * Set the force field parameters for a outOfPlaneBend term.
     * 
     * @param index         the index of the outOfPlaneBend for which to set parameters
     * @param particle1     the index of the first particle connected by the outOfPlaneBend
     * @param particle2     the index of the second particle connected by the outOfPlaneBend
     * @param particle3     the index of the third particle connected by the outOfPlaneBend
     * @param particle4     the index of the fourth particle connected by the outOfPlaneBend
     * @param k             the force constant for the outOfPlaneBend
     */
    void setOutOfPlaneBendParameters(int index, int particle1, int particle2, int particle3, int particle4, double k );

protected:
    ForceImpl* createImpl();
    double _globalCubicK, _globalQuarticK, _globalPenticK, _globalSexticK;
private:

    class OutOfPlaneBendInfo;

// Retarded visual studio compiler complains about being unable to 
// export private stl class members.
// This stanza explains that it should temporarily shut up.
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4251)
#endif
    std::vector<OutOfPlaneBendInfo> outOfPlaneBends;
#if defined(_MSC_VER)
#pragma warning(pop)
#endif

};

class AmoebaOutOfPlaneBendForce::OutOfPlaneBendInfo {
public:
    int particle1, particle2, particle3, particle4;
    double k;
    OutOfPlaneBendInfo() {
        particle1 = particle2  = particle3 = particle4 = -1;
        k   = 0.0;
    }
    OutOfPlaneBendInfo(int particle1, int particle2, int particle3, int particle4,
                           double k ) :
                    particle1(particle1), particle2(particle2), particle3(particle3), particle4(particle4), k(k) {
     
    }
};

} // namespace OpenMM

#endif /*OPENMM_AMOEBA_OUT_OF_PLANE_BEND_FORCE_H_*/
