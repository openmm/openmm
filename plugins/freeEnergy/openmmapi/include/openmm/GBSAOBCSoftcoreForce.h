#ifndef OPENMM_FREE_ENERGY_GBSA_OBC_FORCE_FIELD_H_
#define OPENMM_FREE_ENERGY_GBSA_OBC_FORCE_FIELD_H_

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

#include "openmm/Force.h"
#include <vector>
#include "openmm/internal/windowsExport.h"

namespace OpenMM {

/**
 * This class implements an implicit solvation force using the GBSA-OBC model.
 * <p>
 * To use this class, create a GBSAOBCSoftcoreForce object, then call addParticle() once for each particle in the
 * System to define its parameters.  The number of particles for which you define GBSA parameters must
 * be exactly equal to the number of particles in the System, or else an exception will be thrown when you
 * try to create a Context.  After a particle has been added, you can modify its force field parameters
 * by calling setParticleParameters().
 * <p>
 * If the System also contains a NonbondedForce, this force will use the cutoffs
 * and periodic boundary conditions specified in it.
 */

class OPENMM_EXPORT GBSAOBCSoftcoreForce : public Force {

public:

    /** 
     * This is an enumeration of the different methods that may be used for handling long range nonbonded forces.
     */
    enum NonbondedSoftcoreMethod {
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

    /*
     * Create a GBSAOBCSoftcoreForce.
     */
    GBSAOBCSoftcoreForce();
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
     * @param charge                        the charge of the particle, measured in units of the proton charge
     * @param radius                        the GBSA radius of the particle, measured in nm
     * @param scalingFactor                 the OBC scaling factor for the particle
     * @param nonPolarScalingFactor         the nonpolar scaling factor for the particle (default: 1.0)
     * @return the index of the particle that was added
     * @return the index of the particle that was added
     */
    int addParticle(double charge, double radius, double scalingFactor, double nonPolarScalingFactor = 1.0);
    /**
     * Get the force field parameters for a particle.
     * 
     * @param index                  the index of the particle for which to get parameters
     * @param charge                 the charge of the particle, measured in units of the proton charge
     * @param radius                 the GBSA radius of the particle, measured in nm
     * @param scalingFactor          the OBC scaling factor for the particle
     */
    void getParticleParameters(int index, double& charge, double& radius, double& scalingFactor) const;
    /**
     * Get the force field parameters for a particle.
     * 
     * @param index                  the index of the particle for which to get parameters
     * @param charge                 the charge of the particle, measured in units of the proton charge
     * @param radius                 the GBSA radius of the particle, measured in nm
     * @param scalingFactor          the OBC scaling factor for the particle
     * @param nonPolarScalingFactor  the nonpolar scaling factor for the particle
     */
    void getParticleParameters(int index, double& charge, double& radius, double& scalingFactor, double& nonPolarScalingFactor) const;
    /**
     * Set the force field parameters for a particle.
     * 
     * @param index          the index of the particle for which to set parameters
     * @param charge         the charge of the particle, measured in units of the proton charge
     * @param radius         the GBSA radius of the particle, measured in nm
     * @param scalingFactor  the OBC scaling factor for the particle
     * @param nonPolarScalingFactor  the nonpolar scaling factor for the particle (default: 1.0)
     */
    void setParticleParameters(int index, double charge, double radius, double scalingFactor, double nonPolarScalingFactor = 1.0);
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
     * Get the nonPolarPrefactor.
     */
    double getNonPolarPrefactor() const {
        return nonPolarPrefactor;
    }

    /**
     * Set the nonPolarPrefactor; units are kJ/mol/nm^2
     */
    void setNonPolarPrefactor(double inputNonPolarPrefactor) {
        nonPolarPrefactor = inputNonPolarPrefactor;
    }

    /** 
     * Get the method used for handling long range nonbonded interactions.
     */
    NonbondedSoftcoreMethod getNonbondedMethod() const;

    /** 
     * Set the method used for handling long range nonbonded interactions.
     */
    void setNonbondedMethod(NonbondedSoftcoreMethod method);

    /** 
     * Get the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
     * is NoCutoff, this value will have no effect.
     */
    double getCutoffDistance() const;

    /** 
     * Set the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
     * is NoCutoff, this value will have no effect.
     */
    void setCutoffDistance(double distance);
protected:
    ForceImpl* createImpl();
private:
    class ParticleInfo;
    NonbondedSoftcoreMethod nonbondedMethod;
    double cutoffDistance, solventDielectric, soluteDielectric, nonPolarPrefactor;

// Retarded visual studio compiler complains about being unable to 
// export private stl class members.
// This stanza explains that it should temporarily shut up.
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4251)
#endif

    std::vector<ParticleInfo> particles;

#if defined(_MSC_VER)
#pragma warning(pop)
#endif

};

class GBSAOBCSoftcoreForce::ParticleInfo {
public:
    double charge, radius, scalingFactor, nonPolarScalingFactor;
    ParticleInfo() {
        charge = radius = scalingFactor = nonPolarScalingFactor = 0.0;
    }
    ParticleInfo(double charge, double radius, double scalingFactor, double nonPolarScalingFactor) :
        charge(charge), radius(radius), scalingFactor(scalingFactor), nonPolarScalingFactor(nonPolarScalingFactor) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_FREE_ENERGY_GBSA_OBC_FORCE_FIELD_H_*/
