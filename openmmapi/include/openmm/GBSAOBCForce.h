#ifndef OPENMM_GBSAOBCFORCEFIELD_H_
#define OPENMM_GBSAOBCFORCEFIELD_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
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

#include "Force.h"
#include <vector>
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This class implements an implicit solvation force using the GBSA-OBC model.
 *
 * To use this class, create a GBSAOBCForce object, then call addParticle() once for each particle in the
 * System to define its parameters.  The number of particles for which you define GBSA parameters must
 * be exactly equal to the number of particles in the System, or else an exception will be thrown when you
 * try to create a Context.  After a particle has been added, you can modify its force field parameters
 * by calling setParticleParameters().  This will have no effect on Contexts that already exist unless you
 * call updateParametersInContext().
 *
 * When using this Force, the System should also include a NonbondedForce, and both objects must specify
 * identical charges for all particles.  Otherwise, the results will not be correct.  Furthermore, if the
 * nonbonded method is set to CutoffNonPeriodic or CutoffPeriodic, you should call setReactionFieldDielectric(1.0)
 * on the NonbondedForce to turn off the reaction field approximation, which does not produce correct results
 * when combined with GBSA.
 */

class OPENMM_EXPORT GBSAOBCForce : public Force {
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
         * Interactions beyond the cutoff distance are ignored.
         */
        CutoffNonPeriodic = 1,
        /**
         * Periodic boundary conditions are used, so that each particle interacts only with the nearest periodic copy of
         * each other particle.  Interactions beyond the cutoff distance are ignored.
         */
        CutoffPeriodic = 2,
    };
    /**
     * Create a GBSAOBCForce.
     */
    GBSAOBCForce();
    /**
     * Get the number of particles in the system.
     */
    int getNumParticles() const {
        return particles.size();
    }
    /**
     * Add the GBSA parameters for a particle.  This should be called once for each particle
     * in the System.  When it is called for the i'th time, it specifies the parameters for the i'th particle.
     *
     * @param charge         the charge of the particle, measured in units of the proton charge
     * @param radius         the GBSA radius of the particle, measured in nm
     * @param scalingFactor  the OBC scaling factor for the particle
     * @return the index of the particle that was added
     */
    int addParticle(double charge, double radius, double scalingFactor);
    /**
     * Get the force field parameters for a particle.
     *
     * @param index               the index of the particle for which to get parameters
     * @param[out] charge         the charge of the particle, measured in units of the proton charge
     * @param[out] radius         the GBSA radius of the particle, measured in nm
     * @param[out] scalingFactor  the OBC scaling factor for the particle
     */
    void getParticleParameters(int index, double& charge, double& radius, double& scalingFactor) const;
    /**
     * Set the force field parameters for a particle.
     *
     * @param index          the index of the particle for which to set parameters
     * @param charge         the charge of the particle, measured in units of the proton charge
     * @param radius         the GBSA radius of the particle, measured in nm
     * @param scalingFactor  the OBC scaling factor for the particle
     */
    void setParticleParameters(int index, double charge, double radius, double scalingFactor);
    /**
     * Get the dielectric constant for the solvent.
     */
    double getSolventDielectric() const {
        return solventDielectric;
    }
    /**
     * Set the dielectric constant for the solvent.
     */
    void setSolventDielectric(double dielectric) {
        solventDielectric = dielectric;
    }
    /**
     * Get the dielectric constant for the solute.
     */
    double getSoluteDielectric() const {
        return soluteDielectric;
    }
    /**
     * Set the dielectric constant for the solute.
     */
    void setSoluteDielectric(double dielectric) {
        soluteDielectric = dielectric;
    }
    /**
     * Get the energy scale for the surface energy term, measured in kJ/mol/nm^2.
     */
    double getSurfaceAreaEnergy() const {
        return surfaceAreaEnergy;
    }
    /**
     * Set the energy scale for the surface energy term, measured in kJ/mol/nm^2.
     */
    void setSurfaceAreaEnergy(double energy) {
        surfaceAreaEnergy = energy;
    }
    /**
     * Get the method used for handling long range nonbonded interactions.
     */
    NonbondedMethod getNonbondedMethod() const;
    /**
     * Set the method used for handling long range nonbonded interactions.
     */
    void setNonbondedMethod(NonbondedMethod method);
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
     * Update the particle parameters in a Context to match those stored in this Force object.  This method
     * provides an efficient method to update certain parameters in an existing Context without needing to
     * reinitialize it.  Simply call setParticleParameters() to modify this object's parameters, then call
     * updateParametersInContext() to copy them over to the Context.
     *
     * The only information this method updates is the values of per-particle parameters.  All other aspects
     * of the Force (the nonbonded method, the cutoff distance, etc.) are unaffected and can only be changed
     * by reinitializing the Context.  Furthermore, this method cannot be used to add new particles, only to
     * change the parameters of existing ones.
     */
    void updateParametersInContext(Context& context);
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * @returns true if force uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const {
        return nonbondedMethod == GBSAOBCForce::CutoffPeriodic;
    }
protected:
    ForceImpl* createImpl() const;
private:
    class ParticleInfo;
    NonbondedMethod nonbondedMethod;
    double cutoffDistance, solventDielectric, soluteDielectric, surfaceAreaEnergy;
    std::vector<ParticleInfo> particles;
};

/**
 * This is an internal class used to record information about a particle.
 * @private
 */
class GBSAOBCForce::ParticleInfo {
public:
    double charge, radius, scalingFactor;
    ParticleInfo() {
        charge = radius = scalingFactor = 0.0;
    }
    ParticleInfo(double charge, double radius, double scalingFactor) :
        charge(charge), radius(radius), scalingFactor(scalingFactor) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_GBSAOBCFORCEFIELD_H_*/
