#ifndef OPENMM_AMOEBA_HARMONIC_ANGLE_FORCE_H_
#define OPENMM_AMOEBA_HARMONIC_ANGLE_FORCE_H_

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
 * This class implements an interaction between pairs of particles that varies harmonically with the distance
 * between them.  To use it, create a HarmonicAngleForce object then call addAngle() once for each angle.  After
 * a angle has been added, you can modify its force field parameters by calling setAngleParameters().
 */

 class OPENMM_EXPORT AmoebaHarmonicAngleForce : public Force {
public:
    /**
     * Create a Amoeba HarmonicAngleForce.
     */
    AmoebaHarmonicAngleForce();

    /**
     * Get the number of harmonic angle stretch terms in the potential function
     */
    int getNumAngles() const {
        return angles.size();
    }

    /**
     * Set the global cubic term
     * 
     * @param cubicK        the cubic harmonic force constant for the angle
     */
    void setAmoebaGlobalHarmonicAngleCubic( double cubicK );

    /**
     * Get the global cubic term
     * 
     * @return global cubicK term
     */
    double getAmoebaGlobalHarmonicAngleCubic( void ) const;

    /**
     * Set the global cubic term
     * 
     * @param quarticK       the quartic harmonic force constant for the angle
     */
    void setAmoebaGlobalHarmonicAngleQuartic( double quarticK );

    /**
     * Get the global quartic term
     * 
     * @return global  quartic term
     */
    double getAmoebaGlobalHarmonicAngleQuartic( void ) const;

    /**
     * Set the global pentic term
     * 
     * @param penticK the pentic harmonic force constant for the angle
     */
    void setAmoebaGlobalHarmonicAnglePentic( double penticK );

    /**
     * Get the global pentic term
     * 
     * @return global penticK term
     */
    double getAmoebaGlobalHarmonicAnglePentic( void ) const;

    /**
     * Set the global sextic term
     * 
     * @param sexticK       the sextic harmonic force constant for the angle
     */
    void setAmoebaGlobalHarmonicAngleSextic( double sexticK );

    /**
     * Get the global sextic term
     * 
     * @return global  sextic term
     */
    double getAmoebaGlobalHarmonicAngleSextic( void ) const;

    /**
     * Add a angle term to the force field.
     *
     * @param particle1     the index of the first particle connected by the angle
     * @param particle2     the index of the second particle connected by the angle
     * @param particle3     the index of the third particle connected by the angle
     * @param length        the equilibrium length of the angle, measured in nm
     * @param quadratic k   the quadratic harmonic force constant for the angle
     * @return the index of the angle that was added
     */
    int addAngle(int particle1, int particle2, int particle3, double length, 
                 double quadraticK );

    /**
     * Get the force field parameters for a angle term.
     * 
     * @param index         the index of the angle for which to get parameters
     * @param particle1     the index of the first particle connected by the angle
     * @param particle2     the index of the second particle connected by the angle
     * @param particle3     the index of the third particle connected by the angle
     * @param length        the equilibrium length of the angle, measured in nm
     * @param quadratic k   the quadratic harmonic force constant for the angle
     */
    void getAngleParameters(int index, int& particle1, int& particle2, int& particle3, double& length,
                            double& quadraticK ) const;

    /**
     * Set the force field parameters for a angle term.
     * 
     * @param index         the index of the angle for which to set parameters
     * @param particle1     the index of the first particle connected by the angle
     * @param particle2     the index of the second particle connected by the angle
     * @param particle3     the index of the third particle connected by the angle
     * @param length        the equilibrium length of the angle, measured in nm
     * @param quadratic k   the quadratic harmonic force constant for the angle
     */
    void setAngleParameters(int index, int particle1, int particle2, int particle3, double length, double quadraticK );

protected:
    ForceImpl* createImpl();
    double _globalCubicK, _globalQuarticK, _globalPenticK, _globalSexticK;
private:

    class AngleInfo;

// Retarded visual studio compiler complains about being unable to 
// export private stl class members.
// This stanza explains that it should temporarily shut up.
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4251)
#endif
    std::vector<AngleInfo> angles;
#if defined(_MSC_VER)
#pragma warning(pop)
#endif

};

class AmoebaHarmonicAngleForce::AngleInfo {
public:
    int particle1, particle2, particle3;
    double length, quadraticK;
    AngleInfo() {
        particle1 = particle2  = particle3 = -1;
        length    = quadraticK = 0.0;
    }
    AngleInfo(int particle1, int particle2, int particle3, double length, double  quadraticK ) :
        particle1(particle1), particle2(particle2), particle3(particle3), length(length), quadraticK(quadraticK)  {
    }
};

} // namespace OpenMM

#endif /*OPENMM_AMOEBA_HARMONIC_ANGLE_FORCE_H_*/
