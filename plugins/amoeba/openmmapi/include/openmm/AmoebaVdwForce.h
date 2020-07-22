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
 * Portions copyright (c) 2008-2020 Stanford University and the Authors.      *
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
#include <string>
#include <vector>

namespace OpenMM {

/**
 * This class models van der Waals forces in the AMOEBA force field.  It can use
 * either buffered 14-7 potential or a Lennard-Jones 12-6 potential.
 *
 * This class can operate in two different modes.  In one mode, force field parameters
 * are defined for each particle.  When two particles interact, a combining rule is
 * used to calculate the interaction parameters based on the parameters for the two
 * particles.  To use the class in this mode, call the version of addParticle() that
 * takes sigma and epsilon values.  It should be called once for each particle in the
 * System.
 * 
 * In the other mode, each particle has a type index, and parameters are specified for
 * each type rather than each individual particle.  By default this mode also uses a
 * combining rule, but you can override it by defining alternate parameters to use for
 * specific pairs of particle types.  To use the class in this mode, call the version of
 * addParticle() that takes a type index.  It should be called once for each particle
 * in the System.  You also must call addParticleType() once for each type.  If you
 * wish to override the combining for particular pairs of types, do so by calling
 * addTypePair().
 * 
 * A unique feature of this class is that the interaction site for a particle does not need to be
 * exactly at the particle's location.  Instead, it can be placed a fraction of the distance from that
 * particle to another one.  This is typically done for hydrogens to place the interaction site slightly
 * closer to the parent atom.  The fraction is known as the "reduction factor", since it reduces the distance
 * from the parent atom to the interaction site.
 *
 * Support is also available for softcore interactions based on setting a per particle alchemical flag and
 * setting the AmoebaVdwForce to use an "AlchemicalMethod" -- either Decouple or Annihilate.
 * For Decouple, two alchemical atoms interact normally. For Annihilate, all interactions involving an 
 * alchemical atom are influenced. The softcore state is specified by setting a single 
 * Context parameter "AmoebaVdwLambda" between 0.0 and 1.0.
 *
 * The softcore functional form can be modified by setting the softcore power (default of 5) and the softcore
 * alpha (default of 0,7). For more information on the softcore functional form see Eq. 2 from:
 * Jiao, D.;  Golubkov, P. A.;  Darden, T. A.; Ren, P., 
 * Calculation of protein-ligand binding free energy by using a polarizable potential.
 * Proc. Natl. Acad. Sci. U.S.A. 2008, 105 (17), 6290-6295.
 * https://www.pnas.org/content/105/17/6290.
 */

class OPENMM_EXPORT_AMOEBA AmoebaVdwForce : public Force {
public:
    /**
     * This is the name of the parameter which stores the current Amoeba vdW lambda value.
     */
    static const std::string& Lambda() {
        static const std::string key = "AmoebaVdwLambda";
        return key;
    }

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
     * This is an enumeration of the different potential functions that can be used.
     */
    enum PotentialFunction {
        /**
         * Use a buffered 14-7 potential.  This is the default.
         */
        Buffered147 = 0,
        /**
         * Use a Lennard-Jones 12-6 potential.
         */
        LennardJones = 1,
    };

    /**
     * This is an enumeration of the different alchemical methods used when applying softcore interactions.
     */
    enum AlchemicalMethod {
        /**
         * All vdW interactions are treated normally. This is the default.
         */
        None = 0,
        /**
         * Maintain full strength vdW interactions between two alchemical particles.
         */
        Decouple = 1,
        /**
         * Interactions between two alchemical particles are turned off at lambda=0.
         */
        Annihilate = 2,
    };

    /**
     * Create an Amoeba VdwForce.
     */
    AmoebaVdwForce();

    /**
     * Get the number of particles
     */
    int getNumParticles() const {
        return parameters.size();
    }

    /**
     * Get the number of particle types.
     */
    int getNumParticleTypes() const {
        return types.size();
    }

    /**
     * Get the number of type pairs.
     */
    int getNumTypePairs() const {
        return pairs.size();
    }

    /**
     * Set the force field parameters for a vdw particle.
     *
     * @param particleIndex   the particle index
     * @param parentIndex     the index of the parent particle
     * @param sigma           vdw sigma
     * @param epsilon         vdw epsilon
     * @param reductionFactor the fraction of the distance along the line from the parent particle to this particle
     *                        at which the interaction site should be placed
     * @param isAlchemical    if true, this vdW particle is undergoing an alchemical change.
     * @param typeIndex       the index of the particle type for this particle
     */
    void setParticleParameters(int particleIndex, int parentIndex, double sigma, double epsilon, 
                               double reductionFactor, bool isAlchemical=false, int typeIndex=-1);

    /**
     * Get the force field parameters for a vdw particle.
     *
     * @param particleIndex        the particle index
     * @param[out] parentIndex     the index of the parent particle
     * @param[out] sigma           vdw sigma
     * @param[out] epsilon         vdw epsilon
     * @param[out] reductionFactor the fraction of the distance along the line from the parent particle to this particle
     *                             at which the interaction site should be placed
     * @param[out] isAlchemical    if true, this vdW particle is undergoing an alchemical change.
     * @param[out] typeIndex       the index of the particle type for this particle
     */
    void getParticleParameters(int particleIndex, int& parentIndex, double& sigma, double& epsilon, 
                               double& reductionFactor, bool& isAlchemical, int& typeIndex) const;

    /**
     * Add the force field parameters for a vdw particle.  This version is used when parameters
     * are defined for each particle.
     *
     * @param parentIndex     the index of the parent particle
     * @param sigma           vdw sigma
     * @param epsilon         vdw epsilon
     * @param reductionFactor the fraction of the distance along the line from the parent particle to this particle
     *                        at which the interaction site should be placed
     * @param isAlchemical    if true, this vdW particle is undergoing an alchemical change.
     * @return index of added particle
     */
    int addParticle(int parentIndex, double sigma, double epsilon, double reductionFactor, bool isAlchemical = false);

    /**
     * Add the force field parameters for a vdw particle. This version is used when parameters
     * are defined by particle type.
     *
     * @param parentIndex     the index of the parent particle
     * @param typeIndex       the index of the particle type for this particle
     * @param reductionFactor the fraction of the distance along the line from the parent particle to this particle
     *                        at which the interaction site should be placed
     * @param isAlchemical    if true, this vdW particle is undergoing an alchemical change.
     * @return index of added particle
     */
    int addParticle(int parentIndex, int typeIndex, double reductionFactor, bool isAlchemical = false);

    /**
     * Add a particle type.
     * 
     * @param sigma     the sigma value for particles of this type
     * @param epsilon   the epsilon value for particles of this type
     * @return the index of the particle type that was just added.
     */
    int addParticleType(double sigma, double epsilon);

    /**
     * Get the force field parameters for a particle type.
     * 
     * @param typeIndex      the index of the particle type
     * @param[out] sigma     the sigma value for particles of this type
     * @param[out] epsilon   the epsilon value for particles of this type
     */
    void getParticleTypeParameters(int typeIndex, double& sigma, double& epsilon) const;

    /**
     * Set the force field parameters for a particle type.
     * 
     * @param typeIndex the index of the particle type
     * @param sigma     the sigma value for particles of this type
     * @param epsilon   the epsilon value for particles of this type
     */
    void setParticleTypeParameters(int typeIndex, double sigma, double epsilon);

    /**
     * Add a type pair.  This overrides the standard combining rule for interactions
     * between particles of two particular types.
     * 
     * @param type1     the index of the first particle type
     * @param type2     the index of the second particle type
     * @param sigma     the sigma value for interactions between particles of these two types
     * @param epsilon   the epsilon  value for interactions between particles of these two types
     * @return the index of the type pair that was just added.
     */
    int addTypePair(int type1, int type2, double sigma, double epsilon);

    /**
     * Get the force field parameters for a type pair.  This overrides the standard
     * combining rule for interactions between particles of two particular types.
     * 
     * @param pairIndex      the index of the type pair
     * @param[out] type1     the index of the first particle type
     * @param[out] type2     the index of the second particle type
     * @param[out] sigma     the sigma value for interactions between particles of these two types
     * @param[out] epsilon   the epsilon  value for interactions between particles of these two types
     */
    void getTypePairParameters(int pairIndex, int& type1, int& type2, double& sigma, double& epsilon) const;

    /**
     * Set the force field parameters for a type pair.  This overrides the standard
     * combining rule for interactions between particles of two particular types.
     * 
     * @param pairIndex the index of the type pair
     * @param type1     the index of the first particle type
     * @param type2     the index of the second particle type
     * @param sigma     the sigma value for interactions between particles of these two types
     * @param epsilon   the epsilon  value for interactions between particles of these two types
     */
    void setTypePairParameters(int pairIndex, int type1, int type2, double sigma, double epsilon);
    
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
     * @param epsilonCombiningRule   epsilon combining rule:   'ARITHMETIC', 'GEOMETRIC'. 'HARMONIC', 'W-H', 'HHG'
     */
    void setEpsilonCombiningRule(const std::string& epsilonCombiningRule);

    /**
     * Get epsilon combining rule
     *
     * @return epsilonCombiningRule   epsilon combining rule:  'ARITHMETIC', 'GEOMETRIC'. 'HARMONIC', 'W-H', 'HHG'
     */
    const std::string& getEpsilonCombiningRule(void) const;

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
     * Get whether parameters were specified by particle or by particle type.
     */
    bool getUseParticleTypes() const {
        return useTypes;
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
     * Get the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
     * is NoCutoff, this value will have no effect.
     *
     * @return the cutoff distance, measured in nm
     */

    double getCutoffDistance() const;
    
    /**
     * Set the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
     * is NoCutoff, this value will have no effect.
     *
     * @param distance    the cutoff distance, measured in nm
     */
    void setCutoffDistance(double distance);

    /**
     * Set the cutoff distance.
     * 
     * @deprecated This method exists only for backward compatibility.  Use setCutoffDistance() instead.
     */
    void setCutoff(double cutoff);

    /**
     * Get the cutoff distance.
     * 
     * @deprecated This method exists only for backward compatibility.  Use getCutoffDistance() instead.
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
     * Get the potential function to use.
     */
    PotentialFunction getPotentialFunction() const {
        return potentialFunction;
    }

    /**
     * Set the potential function to use.
     */
    void setPotentialFunction(PotentialFunction potential) {
        potentialFunction = potential;
    }

    /**
     * Set the softcore power on lambda (default = 5).
     */
    void setSoftcorePower(int n);

    /**
     * Get the softcore power on lambda.
     */
    int getSoftcorePower() const;

    /**
     * Set the softcore alpha value (default = 0.7).
     */
    void setSoftcoreAlpha(double alpha);

    /**
     * Get the softcore alpha value.
     */
    double getSoftcoreAlpha() const;

    /**
     * Get the method used for alchemical interactions.
     */
    AlchemicalMethod getAlchemicalMethod() const;

    /**
     * Set the method used for handling long range nonbonded interactions.
     */
    void setAlchemicalMethod(AlchemicalMethod method);

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
    class ParticleTypeInfo;
    class TypePairInfo;
    NonbondedMethod nonbondedMethod;
    PotentialFunction potentialFunction;
    double cutoff;
    bool useDispersionCorrection, useTypes;
    AlchemicalMethod alchemicalMethod;
    int n;
    double alpha;

    std::string sigmaCombiningRule;
    std::string epsilonCombiningRule;

    std::vector< std::vector<int> > exclusions;
    std::vector<VdwInfo> parameters;
    std::vector<ParticleTypeInfo> types;
    std::vector<TypePairInfo> pairs;
};

/**
 * This is an internal class used to record information about a particle.
 * @private
 */
class AmoebaVdwForce::VdwInfo {
public:
    int parentIndex, typeIndex;
    double reductionFactor, sigma, epsilon, cutoff;
    bool isAlchemical;
    VdwInfo() {
        parentIndex = -1;
        reductionFactor      = 0.0;
        sigma                = 1.0;
        epsilon              = 0.0;
        isAlchemical         = false;
    }
    VdwInfo(int parentIndex, double sigma, double epsilon, int typeIndex, double reductionFactor, bool isAlchemical) :
        parentIndex(parentIndex), reductionFactor(reductionFactor), sigma(sigma), epsilon(epsilon), typeIndex(typeIndex), isAlchemical(isAlchemical) {
    }
};

/**
 * This is an internal class used to record information about a particle type.
 * @private
 */
class AmoebaVdwForce::ParticleTypeInfo {
public:
    double sigma, epsilon;
    ParticleTypeInfo() : sigma(1.0), epsilon(0.0) {
    }
    ParticleTypeInfo(double sigma, double epsilon) : sigma(sigma), epsilon(epsilon) {
    }
};

/**
 * This is an internal class used to record information about a type pair.
 * @private
 */
class AmoebaVdwForce::TypePairInfo {
public:
    int type1, type2;
    double sigma, epsilon;
    TypePairInfo() : type1(-1), type2(-1), sigma(1.0), epsilon(0.0) {
    }
    TypePairInfo(int type1, int type2, double sigma, double epsilon) :
        type1(type1), type2(type2), sigma(sigma), epsilon(epsilon) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_AMOEBA_VDW_FORCE_H_*/

