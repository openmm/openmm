#ifndef OPENMM_DRUDEFORCE_H_
#define OPENMM_DRUDEFORCE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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

#include "openmm/Context.h"
#include "openmm/Force.h"
#include <vector>
#include "internal/windowsExportDrude.h"

namespace OpenMM {

/**
 * This class implements forces that are specific to Drude oscillators.  There are two distinct forces
 * it applies: an anisotropic harmonic force connecting each Drude particle to its parent particle; and
 * a screened Coulomb interaction between specific pairs of dipoles.  The latter is typically used between
 * closely bonded particles whose Coulomb interaction would otherwise be fully excluded.
 *
 * To use this class, create a DrudeForce object, then call addParticle() once for each Drude particle in the
 * System to define its parameters.  After a particle has been added, you can modify its force field parameters
 * by calling setParticleParameters().  This will have no effect on Contexts that already exist unless you
 * call updateParametersInContext().  Likewise, call addScreenedPair() for each pair of dipoles (each dipole
 * consisting of a Drude particle and its parent) that should be computed.
 */

class OPENMM_EXPORT_DRUDE DrudeForce : public Force {
public:
    /**
     * Create a DrudeForce.
     */
    DrudeForce();
    /**
     * Get the number of particles for which force field parameters have been defined.
     */
    int getNumParticles() const {
        return particles.size();
    }
    /**
     * Get the number of special interactions that should be calculated differently from other interactions.
     */
    int getNumScreenedPairs() const {
        return screenedPairs.size();
    }
    /**
     * Add a Drude particle to which forces should be applied.
     *
     * @param particle        the index within the System of the Drude particle
     * @param particle1       the index within the System of the particle to which the Drude particle is attached
     * @param particle2       the index within the System of the second particle used for defining anisotropic polarizability.
     *                        This may be set to -1, in which case aniso12 will be ignored.
     * @param particle3       the index within the System of the third particle used for defining anisotropic polarizability.
     *                        This may be set to -1, in which case aniso34 will be ignored.
     * @param particle4       the index within the System of the fourth particle used for defining anisotropic polarizability.
     *                        This may be set to -1, in which case aniso34 will be ignored.
     * @param charge          The charge on the Drude particle
     * @param polarizability  The isotropic polarizability
     * @param aniso12         The scale factor for the polarizability along the direction defined by particle1 and particle2
     * @param aniso34         The scale factor for the polarizability along the direction defined by particle3 and particle4
     * @return the index of the particle that was added
     */
    int addParticle(int particle, int particle1, int particle2, int particle3, int particle4, double charge, double polarizability, double aniso12, double aniso34);
    /**
     * Get the parameters for a Drude particle.
     *
     * @param index                the index of the Drude particle for which to get parameters
     * @param[out] particle        the index within the System of the Drude particle
     * @param[out] particle1       the index within the System of the particle to which the Drude particle is attached
     * @param[out] particle2       the index within the System of the second particle used for defining anisotropic polarizability.
     *                             This may be set to -1, in which case aniso12 will be ignored.
     * @param[out] particle3       the index within the System of the third particle used for defining anisotropic polarizability.
     *                             This may be set to -1, in which case aniso34 will be ignored.
     * @param[out] particle4       the index within the System of the fourth particle used for defining anisotropic polarizability.
     *                             This may be set to -1, in which case aniso34 will be ignored.
     * @param[out] charge          The charge on the Drude particle
     * @param[out] polarizability  The isotropic polarizability
     * @param[out] aniso12         The scale factor for the polarizability along the direction defined by particle1 and particle2
     * @param[out] aniso34         The scale factor for the polarizability along the direction defined by particle3 and particle4
     */
    void getParticleParameters(int index, int& particle, int& particle1, int& particle2, int& particle3, int& particle4, double& charge, double& polarizability, double& aniso12, double& aniso34) const;
    /**
     * Set the parameters for a Drude particle.
     *
     * @param index           the index of the Drude particle for which to set parameters
     * @param particle        the index within the System of the Drude particle
     * @param particle1       the index within the System of the particle to which the Drude particle is attached
     * @param particle2       the index within the System of the second particle used for defining anisotropic polarizability.
     *                        This may be set to -1, in which case aniso12 will be ignored.
     * @param particle3       the index within the System of the third particle used for defining anisotropic polarizability.
     *                        This may be set to -1, in which case aniso34 will be ignored.
     * @param particle4       the index within the System of the fourth particle used for defining anisotropic polarizability.
     *                        This may be set to -1, in which case aniso34 will be ignored.
     * @param charge          The charge on the Drude particle
     * @param polarizability  The isotropic polarizability
     * @param aniso12         The scale factor for the polarizability along the direction defined by particle1 and particle2
     * @param aniso34         The scale factor for the polarizability along the direction defined by particle3 and particle4
     */
    void setParticleParameters(int index, int particle, int particle1, int particle2, int particle3, int particle4, double charge, double polarizability, double aniso12, double aniso34);
    /**
     * Add an interaction to the list of screened pairs.
     *
     * @param particle1  the index within this Force of the first particle involved in the interaction
     * @param particle2  the index within this Force of the second particle involved in the interaction
     * @param thole      the Thole screening factor
     * @return the index of the screenedPair that was added
     */
    int addScreenedPair(int particle1, int particle2, double thole);
    /**
     * Get the force field parameters for screened pair.
     *
     * @param index           the index of the pair for which to get parameters
     * @param[out] particle1  the index within this Force of the first particle involved in the interaction
     * @param[out] particle2  the index within this Force of the second particle involved in the interaction
     * @param[out] thole      the Thole screening factor
     */
    void getScreenedPairParameters(int index, int& particle1, int& particle2, double& thole) const;
    /**
     * Set the force field parameters for screened pair.
     *
     * @param index      the index of the pair for which to get parameters
     * @param particle1  the index within this Force of the first particle involved in the interaction
     * @param particle2  the index within this Force of the second particle involved in the interaction
     * @param thole      the Thole screening factor
     */
    void setScreenedPairParameters(int index, int particle1, int particle2, double thole);
    /**
     * Update the particle and screened pair parameters in a Context to match those stored in this Force object.  This method
     * provides an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setParticleParameters() and setScreenedPairParameters() to modify this object's parameters, then call
     * updateParametersInContext() to copy them over to the Context.
     *
     * This method has several limitations.  It can be used to modify the numeric parameters associated with a particle or
     * screened pair (polarizability, thole, etc.), but not the identities of the particles they involve.  It also cannot
     * be used to add new particles or screenedPairs, only to change the parameters of existing ones.
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
    class ScreenedPairInfo;
    std::vector<ParticleInfo> particles;
    std::vector<ScreenedPairInfo> screenedPairs;
};

/**
 * This is an internal class used to record information about a particle.
 * @private
 */
class DrudeForce::ParticleInfo {
public:
    int particle, particle1, particle2, particle3, particle4;
    double charge, polarizability, aniso12, aniso34;
    ParticleInfo() {
        particle = particle1 = particle2 = particle3 = particle4 = -1;
        charge = polarizability = aniso12 = aniso34 = 0.0;
    }
    ParticleInfo(int particle, int particle1, int particle2, int particle3, int particle4, double charge, double polarizability, double aniso12, double aniso34) :
        particle(particle), particle1(particle1), particle2(particle2), particle3(particle3), particle4(particle4), charge(charge), polarizability(polarizability), aniso12(aniso12), aniso34(aniso34) {
    }
};

/**
 * This is an internal class used to record information about a screenedPair.
 * @private
 */
class DrudeForce::ScreenedPairInfo {
public:
    int particle1, particle2;
    double thole;
    ScreenedPairInfo() {
        particle1 = particle2 = -1;
        thole = 0.0;
    }
    ScreenedPairInfo(int particle1, int particle2, double thole) :
        particle1(particle1), particle2(particle2), thole(thole) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_DRUDEFORCE_H_*/
