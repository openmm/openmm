#ifndef AMOEBA_OPENMM_SASA_FORCE_FIELD_H_
#define AMOEBA_OPENMM_SASA_FORCE_FIELD_H_

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
 * This class implements an SASA force
 * by performing an analytical computation of the weighted
 * solvent accessible surface area of each atom and the first
 * derivatives of the area with respect to Cartesian coordinates
 * 
 * Literature references:
 * 
 * T. J. Richmond, "Solvent Accessible Surface Area and
 * Excluded Volume in Proteins", Journal of Molecular Biology,
 * 178, 63-89 (1984)
 * 
 * L. Wesson and D. Eisenberg, "Atomic Solvation Parameters
 * Applied to Molecular Dynamics of Proteins in Solution",
 * Protein Science, 1, 227-235 (1992)
 *
 * <p>
 * To use this class, create a AmoebaSASAForce object, then call addParticle() once for each particle in the
 * System to define its parameters.  The number of particles for which you define SASA parameters must
 * be exactly equal to the number of particles in the System, or else an exception will be thrown when you
 * try to create a Context.  After a particle has been added, you can modify its force field parameters
 * by calling setParticleParameters().
 */

class OPENMM_EXPORT AmoebaSASAForce : public Force {
public:
    /*
     * Create a AmoebaSASAForce.
     */
    AmoebaSASAForce();

    /**
     * Get the number of particles in the system.
     */
    int getNumParticles() const {
        return particles.size();
    }

    /**
     * Add the SASA parameters for a particle.  This should be called once for each particle
     * in the System.  When it is called for the i'th time, it specifies the parameters for the i'th particle.
     *
     * @param radius         the radius of the particle, measured in nm
     * @param weight         the sasa weight (Tinker asolv())
     * @return the index of the particle that was added
     */
    int addParticle(double radius, double weight);

    /**
     * Get the force field parameters for a particle.
     * 
     * @param index          the index of the particle for which to get parameters
     * @param radius         the radius of the particle, measured in nm
     * @param weight         the sasa weight(Tinker asolv())
     */
    void getParticleParameters(int index, double& radius, double& weight) const;

    /**
     * Set the force field parameters for a particle.
     * 
     * @param index          the index of the particle for which to set parameters
     * @param radius         the radius of the particle, measured in nm
     * @param weight         the sasa weight(Tinker asolv())
     */
    void setParticleParameters(int index, double radius, double weight);

    /**
     * Get the probe radius
     */
    double getProbeRadius() const;

    /**
     * Set the probe radius
     */
    void setProbeRadius(double distance);

protected:
    ForceImpl* createImpl();
private:
    class ParticleInfo;
    double probeRadius;
    std::vector<ParticleInfo> particles;
};

class AmoebaSASAForce::ParticleInfo {
public:
    double radius, weight;
    ParticleInfo() {
        radius = weight = 0.0;
    }
    ParticleInfo(double radius, double weight) :
        radius(radius), weight(weight) {
    }
};

} // namespace OpenMM

#endif /*AMOEBA_OPENMM_GBSA_OBC_FORCE_FIELD_H_*/
