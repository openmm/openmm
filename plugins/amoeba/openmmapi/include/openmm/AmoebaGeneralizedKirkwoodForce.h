#ifndef AMOEBA_OPENMM_GK_FORCE_FIELD_H_
#define AMOEBA_OPENMM_GK_FORCE_FIELD_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
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
 * This class implements an implicit solvation force using the generalized Kirkwood/Grycuk model.
 * <p>
 * To use this class, create an AmoebaGeneralizedKirkwoodForce object, then call addParticle() once for each particle in the
 * System to define its parameters.  The number of particles for which you define parameters must
 * be equal to the number of particles in the System, or else an exception will be thrown when you
 * try to create a Context.  After a particle has been added, you can modify its force field parameters
 * by calling setParticleParameters().  This will have no effect on Contexts that already exist unless you
 * call updateParametersInContext().
 */

class OPENMM_EXPORT_AMOEBA AmoebaGeneralizedKirkwoodForce : public Force {

public:

    /*
     * Create an AmoebaGeneralizedKirkwoodForce.
     */
    AmoebaGeneralizedKirkwoodForce();

    /**
     * Get the number of particles in the system.
     */
    int getNumParticles() const {
        return particles.size();
    }

    /**
     * Add the  parameters for a particle.  This should be called once for each particle
     * in the System.  When it is called for the i'th time, it specifies the parameters for the i'th particle.
     *
     * @param charge         the charge of the particle, measured in units of the proton charge
     * @param radius         the atomic radius of the particle, measured in nm
     * @param scalingFactor  the scaling factor for the particle
     * @return the index of the particle that was added
     */
    int addParticle(double charge, double radius, double scalingFactor);

    /**
     * Get the force field parameters for a particle.
     * 
     * @param      index          the index of the particle for which to get parameters
     * @param[out] charge         the charge of the particle, measured in units of the proton charge
     * @param[out] radius         the atomic radius of the particle, measured in nm
     * @param[out] scalingFactor  the scaling factor for the particle
     */
    void getParticleParameters(int index, double& charge, double& radius, double& scalingFactor) const;

    /**
     * Set the force field parameters for a particle.
     * 
     * @param index          the index of the particle for which to set parameters
     * @param charge         the charge of the particle, measured in units of the proton charge
     * @param radius         the atomic radius of the particle, measured in nm
     * @param scalingFactor  the scaling factor for the particle
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
     * Get the flag signaling whether the cavity term should be included
     */
    int getIncludeCavityTerm() const;

    /**
     * Set the flag signaling whether the cavity term should be included
     */
    void setIncludeCavityTerm(int includeCavityTerm);

    /**
     * Get the probe radius (nm) used in SASA contribution
     */
    double getProbeRadius() const;

    /**
     * Set the probe radius (nm) used in SASA contribution
     */
    void setProbeRadius(double probeRadius);

    /**
     * Get the surface area factor kJ/(nm*nm) used in SASA contribution
     */
    double getSurfaceAreaFactor() const;

    /**
     * Set the surface area factor kJ/(nm*nm) used in SASA contribution
     */
    void setSurfaceAreaFactor(double surfaceAreaFactor);
    /**
     * Update the per-particle parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setParticleParameters() to modify this object's parameters, then call updateParametersInContext()
     * to copy them over to the Context.
     * 
     * The only information this method updates is the values of per-particle parameters.  All other aspects of the Force
     * (the probe radius, the surface area factor, etc.) are unaffected and can only be changed by reinitializing the Context.
     */
    void updateParametersInContext(Context& context);
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * @returns true if nonbondedMethod uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const {
        return false;
    }
protected:
    ForceImpl* createImpl() const;
private:
    class ParticleInfo;
    int includeCavityTerm;
    double solventDielectric, soluteDielectric, dielectricOffset,
           probeRadius, surfaceAreaFactor;
    std::vector<ParticleInfo> particles;
};

/**
 * This is an internal class used to record information about a particle.
 * @private
 */
class AmoebaGeneralizedKirkwoodForce::ParticleInfo {
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

#endif /*AMOEBA_OPENMM_GK_FORCE_FIELD_H_*/
