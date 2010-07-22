#ifndef OPENMM_AMOEBA_HARMONIC_BOND_FORCE_H_
#define OPENMM_AMOEBA_HARMONIC_BOND_FORCE_H_

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
#include "openmm/Vec3.h"
#include <map>
#include <vector>
#include "openmm/internal/windowsExport.h"

//namespace AmoebaOpenMM {
namespace OpenMM {

/**
 * This class implements an interaction between pairs of particles that varies harmonically with the distance
 * between them.  To use it, create a HarmonicBondForce object then call addBond() once for each bond.  After
 * a bond has been added, you can modify its force field parameters by calling setBondParameters().
 */

 class OPENMM_EXPORT AmoebaHarmonicBondForce : public Force {
public:
    /**
     * Create a Amoeba HarmonicBondForce.
     */
    AmoebaHarmonicBondForce();
    /**
     * Get the number of harmonic bond stretch terms in the potential function
     */
    int getNumBonds() const {
        return bonds.size();
    }

    /**
     * Set the global cubic term
     * 
     * @param cubicK        the cubic harmonic force constant for the bond
     */
    void setAmoebaGlobalHarmonicBondCubic( double cubicK );

    /**
     * Get the global cubic term
     * 
     * @return global cubicK term
     */
    double getAmoebaGlobalHarmonicBondCubic( void ) const;

    /**
     * Set the global cubic term
     * 
     * @param quarticK       the quartic harmonic force constant for the bond
     */
    void setAmoebaGlobalHarmonicBondQuartic( double quarticK );

    /**
     * Get the global quartic term
     * 
     * @return global  quartic term
     */
    double getAmoebaGlobalHarmonicBondQuartic( void ) const;

    /**
     * Add a bond term to the force field.
     *
     * @param particle1     the index of the first particle connected by the bond
     * @param particle2     the index of the second particle connected by the bond
     * @param length        the equilibrium length of the bond, measured in nm
     * @param cubic k       the cubic harmonic force constant for the bond
     * @param quartic k     the quartic harmonic force constant for the bond
     * @param quadratic k   the quadratic harmonic force constant for the bond
     * @return the index of the bond that was added
     */
    int addBond(int particle1, int particle2, double length, double quadraticK, double cubicK = DEFAULT_GLOBAL_K, double quarticK = DEFAULT_GLOBAL_K );

    /**
     * Get the force field parameters for a bond term.
     * 
     * @param index     the index of the bond for which to get parameters
     * @param particle1 the index of the first particle connected by the bond
     * @param particle2 the index of the second particle connected by the bond
     * @param length    the equilibrium length of the bond, measured in nm
     * @param quadratic k   the quadratic harmonic force constant for the bond
     * @param cubic k       the cubic harmonic force constant for the bond
     * @param quartic k     the quartic harmonic force constant for the bond
     */
    void getBondParameters(int index, int& particle1, int& particle2, double& length, double& quadraticK, double& cubicK, double& quarticK ) const;

    /**
     * Set the force field parameters for a bond term.
     * 
     * @param index     the index of the bond for which to set parameters
     * @param particle1 the index of the first particle connected by the bond
     * @param particle2 the index of the second particle connected by the bond
     * @param length    the equilibrium length of the bond, measured in nm
     * @param cubic k       the cubic harmonic force constant for the bond
     * @param quartic k     the quartic harmonic force constant for the bond
     * @param quadratic k   the quadratic harmonic force constant for the bond
     */
    void setBondParameters(int index, int particle1, int particle2, double length, double quadraticK, double cubicK = DEFAULT_GLOBAL_K, double quarticK = DEFAULT_GLOBAL_K);

protected:
    double _globalQuarticK, _globalCubicK;
    ForceImpl* createImpl();
private:

    static const double DEFAULT_GLOBAL_K;
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

class AmoebaHarmonicBondForce::BondInfo {
public:
    int particle1, particle2;
    double length, quadraticK, cubicK, quarticK;
    BondInfo() {
        particle1 = particle2    = -1;
        length    = quadraticK   = 0.0;
        cubicK    = quarticK     = DEFAULT_GLOBAL_K;
    }
    BondInfo(int particle1, int particle2, double length, double  quadraticK, double cubicK, double quarticK) :
        particle1(particle1), particle2(particle2), length(length), quadraticK(quadraticK), cubicK(cubicK), quarticK(quarticK) {
    }
};

} // namespace AmoebaOpenMM

#endif /*OPENMM_AMOEBA_HARMONIC_BOND_FORCE_H_*/
