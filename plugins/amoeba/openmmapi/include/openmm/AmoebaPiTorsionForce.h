#ifndef OPENMM_AMOEBA_PI_TORSION_FORCE_H_
#define OPENMM_AMOEBA_PI_TORSION_FORCE_H_

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
 * This class implements the Amoeba pi-torsion interaction.
 * 
 * To use it, create an AmoebaPiTorsionForce object then call addPiTorsion() once for each torsion.  After
 * a torsion has been added, you can modify its force field parameters by calling setPiTorsionParameters().
 * This will have no effect on Contexts that already exist unless you call updateParametersInContext().
 */

class OPENMM_EXPORT_AMOEBA AmoebaPiTorsionForce : public Force {

public:

    /**
     * Create an AmoebaPiTorsionForce.
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
    int addPiTorsion(int particle1, int particle2, int particle3, int particle4, int particle5, int particle6, double k);

    /**
     * Get the force field parameters for a torsion term.
     * 
     * @param index              the index of the torsion for which to get parameters
     * @param[out] particle1     the index of the first particle connected by the torsion
     * @param[out] particle2     the index of the second particle connected by the torsion
     * @param[out] particle3     the index of the third particle connected by the torsion
     * @param[out] particle4     the index of the fourth particle connected by the torsion
     * @param[out] particle5     the index of the fifth particle connected by the torsion
     * @param[out] particle6     the index of the sixth particle connected by the torsion
     * @param[out] k             the force constant for the torsion
     */
    void getPiTorsionParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4, int& particle5, int& particle6, double& k) const;

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
    void setPiTorsionParameters(int index, int particle1, int particle2, int particle3, int particle4, int particle5, int particle6, double k);
    /**
     * Update the per-torsion parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setPiTorsionParameters() to modify this object's parameters, then call updateParametersInContext()
     * to copy them over to the Context.
     * 
     * The only information this method updates is the values of per-torsion parameters.  The set of particles involved
     * in a torsion cannot be changed, nor can new torsions be added.
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
private:
    class PiTorsionInfo;
    std::vector<PiTorsionInfo> piTorsions;
    bool usePeriodic;
};

/**
 * This is an internal class used to record information about a torsion.
 * @private
 */
class AmoebaPiTorsionForce::PiTorsionInfo {
public:
    int particle1, particle2, particle3, particle4, particle5, particle6;
    double k;
    PiTorsionInfo() {
        particle1 = particle2  = particle3 = particle4 = particle5 = particle6 = -1;
        k = 0.0;
    }
    PiTorsionInfo(int particle1, int particle2, int particle3, int particle4, int particle5, int particle6, double k) :
        particle1(particle1), particle2(particle2), particle3(particle3), particle4(particle4), particle5(particle5), particle6(particle6), k(k) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_AMOEBA_PI_TORSION_FORCE_H_*/
