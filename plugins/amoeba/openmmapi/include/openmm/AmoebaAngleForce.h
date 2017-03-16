#ifndef OPENMM_AMOEBA_ANGLE_FORCE_H_
#define OPENMM_AMOEBA_ANGLE_FORCE_H_

/* -------------------------------------------------------------------------- *
 *                              OpenMMAmoeba                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
 * Authors: Mark Friedrichs, Peter Eastman                                    *
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
#include "internal/windowsExportAmoeba.h"
#include <vector>

namespace OpenMM {

/**
 * This class implements an interaction between triplets of particles that varies with the angle
 * between them.  The interaction is defined by a 6th order polynomial.  Only the quadratic term
 * is set per-angle.  The coefficients of the higher order terms each have a single value that
 * is set globally.
 * 
 * To use it, create an AmoebaAngleForce object then call addAngle() once for each angle.  After
 * an angle has been added, you can modify its force field parameters by calling setAngleParameters().
 * This will have no effect on Contexts that already exist unless you call updateParametersInContext().
 */

class OPENMM_EXPORT_AMOEBA AmoebaAngleForce : public Force {

public:

    /**
     * Create an AmoebaAngleForce.
     */
    AmoebaAngleForce();

    /**
     * Get the number of angle stretch terms in the potential function
     */
    int getNumAngles() const {
        return angles.size();
    }

    /**
     * Set the global cubic term
     * 
     * @param cubicK        the cubic force constant for the angle
     */
    void setAmoebaGlobalAngleCubic(double cubicK);

    /**
     * Get the global cubic term
     * 
     * @return global cubicK term
     */
    double getAmoebaGlobalAngleCubic() const;

    /**
     * Set the global quartic term
     * 
     * @param quarticK       the quartic force constant for the angle
     */
    void setAmoebaGlobalAngleQuartic(double quarticK);

    /**
     * Get the global quartic term
     * 
     * @return global  quartic term
     */
    double getAmoebaGlobalAngleQuartic() const;

    /**
     * Set the global pentic term
     * 
     * @param penticK the pentic force constant for the angle
     */
    void setAmoebaGlobalAnglePentic(double penticK);

    /**
     * Get the global pentic term
     * 
     * @return global penticK term
     */
    double getAmoebaGlobalAnglePentic() const;

    /**
     * Set the global sextic term
     * 
     * @param sexticK       the sextic force constant for the angle
     */
    void setAmoebaGlobalAngleSextic(double sexticK);

    /**
     * Get the global sextic term
     * 
     * @return global sextic term
     */
    double getAmoebaGlobalAngleSextic() const;

    /**
     * Add an angle term to the force field.
     *
     * @param particle1     the index of the first particle connected by the angle
     * @param particle2     the index of the second particle connected by the angle
     * @param particle3     the index of the third particle connected by the angle
     * @param length        the equilibrium angle, measured in degrees
     * @param quadraticK    the quadratic force constant for the angle, measured in kJ/mol/radian^2
     * @return              the index of the angle that was added
     */
    int addAngle(int particle1, int particle2, int particle3, double length, double quadraticK);

    /**
     * Get the force field parameters for an angle term.
     * 
     * @param index              the index of the angle for which to get parameters
     * @param[out] particle1     the index of the first particle connected by the angle
     * @param[out] particle2     the index of the second particle connected by the angle
     * @param[out] particle3     the index of the third particle connected by the angle
     * @param[out] length        the equilibrium angle, measured in degrees
     * @param[out] quadraticK    the quadratic force constant for the angle, measured in kJ/mol/radian^2
     */
    void getAngleParameters(int index, int& particle1, int& particle2, int& particle3, double& length, double& quadraticK) const;

    /**
     * Set the force field parameters for an angle term.
     * 
     * @param index         the index of the angle for which to set parameters
     * @param particle1     the index of the first particle connected by the angle
     * @param particle2     the index of the second particle connected by the angle
     * @param particle3     the index of the third particle connected by the angle
     * @param length        the equilibrium angle, measured in degrees
     * @param quadraticK    the quadratic force constant for the angle, measured in kJ/mol/radian^2
     */
    void setAngleParameters(int index, int particle1, int particle2, int particle3, double length, double quadraticK);
    /**
     * Update the per-angle parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setAngleParameters() to modify this object's parameters, then call updateParametersInContext()
     * to copy them over to the Context.
     * 
     * The only information this method updates is the values of per-angle parameters.  The set of particles involved
     * in an angle cannot be changed, nor can new angles be added.
     */
    void updateParametersInContext(Context& context);
    /**
     * Set whether this force should apply periodic boundary conditions when calculating displacements.
     * Usually this is not appropriate for bonded forces, but there are situations when it can be useful.
     */
    void setUsesPeriodicBoundaryConditions(bool periodic);
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * @returns true if force uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const;
protected:
    ForceImpl* createImpl() const;
    double _globalCubicK, _globalQuarticK, _globalPenticK, _globalSexticK;
private:
    class AngleInfo;
    std::vector<AngleInfo> angles;
    bool usePeriodic;
};

/**
 * This is an internal class used to record information about an angle.
 * @private
 */
class AmoebaAngleForce::AngleInfo {
public:
    int particle1, particle2, particle3;
    double length, quadraticK;
    AngleInfo() {
        particle1 = particle2  = particle3 = -1;
        length    = quadraticK = 0.0;
    }
    AngleInfo(int particle1, int particle2, int particle3, double length, double  quadraticK) :
        particle1(particle1), particle2(particle2), particle3(particle3), length(length), quadraticK(quadraticK)  {
    }
};

} // namespace OpenMM

#endif /*OPENMM_AMOEBA_ANGLE_FORCE_H_*/
