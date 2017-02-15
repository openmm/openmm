#ifndef OPENMM_AMOEBA_VDW_FORCE_H_
#define OPENMM_AMOEBA_VDW_FORCE_H_

/* -------------------------------------------------------------------------- *
 *                              OpenMMAmoeba                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
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
 * This class implements a buffered 14-7 potential used to model van der Waals forces.
 *
 * To use it, create an AmoebaVdwForce object then call addParticle() once for each particle.  Then call
 * addVdwprByOldTypes() to use special vdwpr parameters.  After particles have been added,
 * you can modify its force field parameters by calling setParticleParameters() and setVdwprParameters().  This
 * will have no effect on Contexts that already exist unless you call updateParametersInContext().
 *
 * A unique feature of this class is that the interaction site for a particle does not need to be
 * exactly at the particle's location.  Instead, it can be placed a fraction of the distance from that
 * particle to another one.  This is typically done for hydrogens to place the interaction site slightly
 * closer to the parent atom.  The fraction is known as the "reduction factor", since it reduces the distance
 * from the parent atom to the interaction site.
 */

class OPENMM_EXPORT_AMOEBA AmoebaVdwForce : public Force {
public:
    /**
     * This is an enumeration of the different methods that may be used for handling long range nonbonded forces.
     */
    enum NonbondedMethod {
        /**
         * No cutoff is applied to nonbonded interactions.  The full set of N^2 interactions is computed exactly.
         * This necessarily means that periodic boundary conditions cannot be used.  This is the default.
         */
        NoCutoff = 0,
        /**
         * Periodic boundary conditions are used, so that each particle interacts only with the nearest periodic copy of
         * each other particle.  Interactions beyond the cutoff distance are ignored.
         */
        CutoffPeriodic = 1,
    };

    /**
     * Create an Amoeba VdwForce.
     */
    AmoebaVdwForce();

    /**
     * Get the number of particles.
     */
    int getNumParticles() const {
        return parameters.size();
    }

    /**
     * Set the force field parameters for a vdw particle.
     *
     * @param particleIndex   the particle index
     * @param parentIndex     the index of the parent particle
     * @param vdwprType       the vdw class/type number
     * @param sigma           vdw sigma
     * @param epsilon         vdw epsilon
     * @param reductionFactor the fraction of the distance along the line from the parent particle to this particle
     *                        at which the interaction site should be placed
     * @param lambda          the alchemical order
     */
    void setParticleParameters(int particleIndex, int parentIndex, int vdwprType, double sigma, double epsilon, double reductionFactor, double lambda);

    /**
     * Get the force field parameters for a vdw particle.
     *
     * @param particleIndex        the particle index
     * @param[out] parentIndex     the index of the parent particle
     * @param[out] vdwprType       the vdw class/type number
     * @param[out] sigma           vdw sigma
     * @param[out] epsilon         vdw epsilon
     * @param[out] reductionFactor the fraction of the distance along the line from the parent particle to this particle
     *                             at which the interaction site should be placed
     * @param[out] lambda          the alchemical order
     */
    void getParticleParameters(int particleIndex, int& parentIndex, int& vdwprType, double& sigma, double& epsilon, double& reductionFactor, double& lambda) const;


    /**
     * Add the force field parameters for a vdw particle.
     *
     * @param parentIndex     the index of the parent particle
     * @param vdwprType       the vdw class/type number
     * @param sigma           vdw sigma
     * @param epsilon         vdw epsilon
     * @param reductionFactor the fraction of the distance along the line from the parent particle to this particle
     *                        at which the interaction site should be placed
     * @param lambda          the alchemical order
     * @return index of added particle
     */
    int addParticle(int parentIndex, int vdwprType, double sigma, double epsilon, double reductionFactor, double lambda);

    /**
     * Compute the combined sigmas and epsilons using combining rules.
     */
    void computeCombinedSigmaEpsilon();

    /**
     * Get the number of vdw classes/types.
     */
    int getNumVdwprTypes() const {
        return numVdwprTypes;
    }

    /**
     * Set the number of vdw classes/types.
     *
     * @param newNum the new number of vdw classes/types
     */
    void setNumVdwprTypes(int newNum);

    /**
     * Get the new vdw class/type.
     *
     * @param oldType the old vdw class/type
     * @return the new vdw class/type number
     */
    int getNewVdwprType(int oldType) const;

    /**
     * Get the old vdw class/type number from a new vdw class/type number.
     *
     * @param newType the new vdw class/type number
     * @return the old vdw class/type number
     */
    int getOldVdwprType(int newType) const;

    /**
     * Set an old vdw class/type number with a new vdw class/type.
     *
     * @param newType the new vdw class/type number
     * @param oldType the old vdw class/type number
     */
    void setOldVdwprType(int newType, int oldType);

    /**
     * Resize some internal storing variables.
     *
     * @param newSize the new number of vdw classes/types.
     */
    void resize(int newSize);

    /**
     * Set vdwpr parameters by old vdw class/type numbers.
     *
     * @param oldtype1        the old vdw class/type number 1
     * @param oldtype2        the old vdw class/type number 2
     * @param combinedSigma   combined sigma
     * @param combinedEpsilon combined epsilon
     */
    void setVdwprParametersByOldTypes(int oldtype1, int oldtype2, double combinedSigma, double combinedEpsilon);

    /**
     * Add vdwpr parameters by old vdw class/type numbers.
     *
     * @param oldtype1        the old vdw class/type number 1
     * @param oldtype2        the old vdw class/type number 2
     * @param combinedSigma   combined sigma
     * @param combinedEpsilon combined epsilon
     * @return the indexed of the added pair
     */
    int addVdwprByOldTypes(int oldtype1, int oldtype2, double combinedSigma, double combinedEpsilon);

    /**
     * Set vdwpr parameters by new vdw class/type numbers.
     *
     * @param ntype1          the new vdw class/type number 1
     * @param ntype2          the new vdw class/type number 2
     * @param combinedSigma   combined sigma
     * @param combinedEpsilon combined epsilon
     */
    void setVdwprParameters(int ntype1, int ntype2, double combinedSigma, double combinedEpsilon);

    /**
     * Get vdwpr parameters by new vdw class/type numbers.
     *
     * @param ntype1          the new vdw class/vdw number 1
     * @param ntype2          the new vdw class/vdw number 2
     * @param combinedSigma   combined sigma
     * @param combinedEpsilon combined epsilon
     */
    void getVdwprParameters(int ntype1, int ntype2, double& combinedSigma, double& combinedEpsilon) const;

    /**
     * Add vdwpr parameters by new vdw class/type numbers.
     *
     * @param ntype1          the new vdw class/type number 1
     * @param ntype2          the new vdw class/type number 2
     * @param combinedSigma   combined sigma
     * @param combinedEpsilon combined epsilon
     */
    int addVdwpr(int ntype1, int ntype2, double combinedSigma, double combinedEpsilon);

    /**
     * Set sigma combining rule
     *
     * @param sigmaCombiningRule   sigma combining rule:  'ARITHMETIC', 'GEOMETRIC'. 'CUBIC-MEAN'
     */
    void setSigmaCombiningRule(const std::string& sigmaCombiningRule);

    /**
     * Get sigma combining rule
     *
     * @return sigmaCombiningRule   sigma combining rule:  'ARITHMETIC', 'GEOMETRIC'. 'CUBIC-MEAN'
     */
    const std::string& getSigmaCombiningRule(void) const;

    /**
     * Set epsilon combining rule
     *
     * @param epsilonCombiningRule   epsilon combining rule:   'ARITHMETIC', 'GEOMETRIC'. 'HARMONIC', 'HHG'
     */
    void setEpsilonCombiningRule(const std::string& epsilonCombiningRule);

    /**
     * Get epsilon combining rule
     *
     * @return epsilonCombiningRule   epsilon combining rule:  'ARITHMETIC', 'GEOMETRIC'. 'HARMONIC', 'HHG'
     */
    const std::string& getEpsilonCombiningRule(void) const;

    /**
     * Set vdw functional form
     *
     * @param functionalForm   functional form:  'BUFFERED-14-7', 'LENNARD-JONES'
     */
    void setFunctionalForm(const std::string& functionalForm);

    /**
     * Get vdw functional form
     *
     * @return functionalForm   functional form:  'BUFFERED-14-7', 'LENNARD-JONES'
     */
    const std::string& getFunctionalForm(void) const;

    /**
     * Get whether to add a contribution to the energy that approximately represents the effect of VdW
     * interactions beyond the cutoff distance.  The energy depends on the volume of the periodic box, and is only
     * applicable when periodic boundary conditions are used.  When running simulations at constant pressure, adding
     * this contribution can improve the quality of results.
     */
    bool getUseDispersionCorrection() const {
        return useDispersionCorrection;
    }

    /**
     * Set whether to add a contribution to the energy that approximately represents the effect of VdW
     * interactions beyond the cutoff distance.  The energy depends on the volume of the periodic box, and is only
     * applicable when periodic boundary conditions are used.  When running simulations at constant pressure, adding
     * this contribution can improve the quality of results.
     */
    void setUseDispersionCorrection(bool useCorrection) {
        useDispersionCorrection = useCorrection;
    }

    /**
     * Set exclusions for specified particle
     *
     * @param particleIndex particle index
     * @param exclusions vector of exclusions
     */
    void setParticleExclusions(int particleIndex, const std::vector<int>& exclusions);

    /**
     * Get exclusions for specified particle
     *
     * @param particleIndex   particle index
     * @param[out] exclusions vector of exclusions
     */
    void getParticleExclusions(int particleIndex, std::vector<int>& exclusions) const;

    /**
     * Set the cutoff distance.
     */
    void setCutoff(double cutoff);

    /**
     * Get the cutoff distance.
     */
    double getCutoff() const;

    /**
     * Get the method used for handling long range nonbonded interactions.
     */
    NonbondedMethod getNonbondedMethod() const;

    /**
     * Set the method used for handling long range nonbonded interactions.
     */
    void setNonbondedMethod(NonbondedMethod method);
    /**
     * Update the per-particle parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setParticleParameters() to modify this object's parameters, then call updateParametersInContext()
     * to copy them over to the Context.
     *
     * The only information this method updates is the values of per-particle parameters.  All other aspects of the Force
     * (the nonbonded method, the cutoff distance, etc.) are unaffected and can only be changed by reinitializing the Context.
     */
    void updateParametersInContext(Context& context);
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * @returns true if nonbondedMethod uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const {
        return nonbondedMethod == AmoebaVdwForce::CutoffPeriodic;
    }

protected:
    ForceImpl* createImpl() const;

private:
    class VdwInfo;
    class VdwprInfo;
    NonbondedMethod nonbondedMethod;
    double cutoff;
    bool useDispersionCorrection;
    int numVdwprTypes;

    std::string sigmaCombiningRule;
    std::string epsilonCombiningRule;
    std::string functionalForm;

    std::vector< std::vector<int> > exclusions; // size = number of atoms
    std::vector<VdwInfo> parameters; // size = number of atoms
    std::vector<VdwprInfo> arguments; // size = (number of vdw classes/types)^2
    std::vector<int> typeMap; // size = number of vdw classes/types
};

/**
 * This is an internal class used to record information about a particle.
 * @private
 */
class AmoebaVdwForce::VdwInfo {
public:
    int parentIndex, vdwprType;
    double reductionFactor, sigma, epsilon, lambda;

    VdwInfo()
        : parentIndex(-1)
        , vdwprType(-1)
        , reductionFactor(0.0)
        , sigma(1.0)
        , epsilon(0.0)
        , lambda(1.0) {}

    VdwInfo(int parentIndex, int vdwprType, double sigma, double epsilon, double reductionFactor, double lambda)
        : parentIndex(parentIndex)
        , vdwprType(vdwprType)
        , reductionFactor(reductionFactor)
        , sigma(sigma)
        , epsilon(epsilon)
        , lambda(lambda) {}
};

/**
 * This is an internal class used to record information about vdw pair.
 * @private
 */
class AmoebaVdwForce::VdwprInfo {
public:
    double combinedSigma;
    double combinedEpsilon;

    VdwprInfo()
        : combinedSigma(1.0)
        , combinedEpsilon(0.0) {}

    VdwprInfo(double combinedSigma, double combinedEpsilon)
        : combinedSigma(combinedSigma)
        , combinedEpsilon(combinedEpsilon) {}
};
} // namespace OpenMM

#endif /*OPENMM_AMOEBA_VDW_FORCE_H_*/

