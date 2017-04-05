#ifndef OPENMM_AMOEBA_BOND_FORCE_H_
#define OPENMM_AMOEBA_BOND_FORCE_H_

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
#include "openmm/Vec3.h"
#include "internal/windowsExportAmoeba.h"
#include <map>
#include <vector>

namespace OpenMM {

/**
 * This class implements an interaction between pairs of particles that varies with the distance
 * between them.  The interaction is defined by a 4th order polynomial.  Only the quadratic term
 * is set per-bond.  The coefficients of the higher order terms each have a single value that
 * is set globally.
 * 
 * To use it, create an AmoebaBondForce object then call addBond() once for each bond.  After
 * a bond has been added, you can modify its force field parameters by calling setBondParameters().
 * This will have no effect on Contexts that already exist unless you call updateParametersInContext().
 */

class OPENMM_EXPORT_AMOEBA AmoebaBondForce : public Force {

public:

    /**
     * Create an AmoebaBondForce.
     */
    AmoebaBondForce();
    /**
     * Get the number of bond stretch terms in the potential function
     */
    int getNumBonds() const {
        return bonds.size();
    }

    /**
     * Set the global cubic term
     * 
     * @param cubicK        the cubic force constant for the bond
     */
    void setAmoebaGlobalBondCubic(double cubicK);

    /**
     * Get the global cubic term
     * 
     * @return global cubicK term
     */
    double getAmoebaGlobalBondCubic() const;

    /**
     * Set the global quartic term
     * 
     * @param quarticK       the quartic force constant for the bond
     */
    void setAmoebaGlobalBondQuartic(double quarticK);

    /**
     * Get the global quartic term
     * 
     * @return global  quartic term
     */
    double getAmoebaGlobalBondQuartic() const;

    /**
     * Add a bond term to the force field.
     *
     * @param particle1     the index of the first particle connected by the bond
     * @param particle2     the index of the second particle connected by the bond
     * @param length        the equilibrium length of the bond, measured in nm
     * @param quadraticK    the quadratic force constant for the bond
     * @return the index of the bond that was added
     */

    int addBond(int particle1, int particle2, double length, double quadraticK);

    /**
     * Get the force field parameters for a bond term.
     * 
     * @param index              the index of the bond for which to get parameters
     * @param[out] particle1     the index of the first particle connected by the bond
     * @param[out] particle2     the index of the second particle connected by the bond
     * @param[out] length        the equilibrium length of the bond, measured in nm
     * @param[out] quadraticK    the quadratic force constant for the bond
     */

    void getBondParameters(int index, int& particle1, int& particle2, double& length, double& quadraticK) const;

    /**
     * Set the force field parameters for a bond term.
     * 
     * @param index       the index of the bond for which to set parameters
     * @param particle1   the index of the first particle connected by the bond
     * @param particle2   the index of the second particle connected by the bond
     * @param length      the equilibrium length of the bond, measured in nm
     * @param quadraticK  the quadratic force constant for the bond
     */
    void setBondParameters(int index, int particle1, int particle2, double length, double quadraticK);
    /**
     * Update the per-bond parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setBondParameters() to modify this object's parameters, then call updateParametersInContext()
     * to copy them over to the Context.
     * 
     * The only information this method updates is the values of per-bond parameters.  The set of particles involved
     * in a bond cannot be changed, nor can new bonds be added.
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
    double _globalQuarticK, _globalCubicK;
    ForceImpl* createImpl() const;
private:
    class BondInfo;
    std::vector<BondInfo> bonds;
    bool usePeriodic;
};

/**
 * This is an internal class used to record information about a bond.
 * @private
 */
class AmoebaBondForce::BondInfo {
public:
    int particle1, particle2;
    double length, quadraticK;
    BondInfo() {
        particle1 = particle2    = -1;
        length    = quadraticK   = 0.0;
    }
    BondInfo(int particle1, int particle2, double length, double  quadraticK) :
        particle1(particle1), particle2(particle2), length(length), quadraticK(quadraticK) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_AMOEBA_BOND_FORCE_H_*/
