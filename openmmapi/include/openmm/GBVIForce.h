#ifndef OPENMM_GBVIFORCEFIELD_H_
#define OPENMM_GBVIFORCEFIELD_H_

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
 * This class implements an implicit solvation force using the GB/VI model.
 * <p>
 * To use this class, create a GBVIForce object, then call addParticle() once for each particle in the
 * System to define its parameters.  The number of particles for which you define GB/VI parameters must
 * be exactly equal to the number of particles in the System, or else an exception will be thrown when you
 * try to create an OpenMMContext.  After a particle has been added, you can modify its force field parameters
 * by calling setParticleParameters().
 * <p>
 * If the System also contains a NonbondedForce, this force will use the cutoffs
 * and periodic boundary conditions specified in it.
 */

class OPENMM_EXPORT GBVIForce : public Force {
public:
    /*
     * Create a GBVIForce.
     */
    GBVIForce();
    /**
     * Get the number of particles in the system.
     */
    int getNumParticles() const {
        return particles.size();
    }
    /**
     * Add the GB/VI parameters for a particle.  This should be called once for each particle
     * in the System.  When it is called for the i'th time, it specifies the parameters for the i'th particle.
     *
     * @param charge         the charge of the particle, measured in units of the proton charge
     * @param radius         the GB/VI radius of the particle, measured in nm
     * @param gamma          the gamma parameter
     * @return the index of the particle that was added
     */
    int addParticle(double charge, double radius, double gamma);
    /**
     * Get the force field parameters for a particle.
     * 
     * @param index          the index of the particle for which to get parameters
     * @param charge         the charge of the particle, measured in units of the proton charge
     * @param radius         the GBSA radius of the particle, measured in nm
     * @param gamma          the gamma parameter
     */
    void getParticleParameters(int index, double& charge, double& radius, double& gamma) const;
    /**
     * Set the force field parameters for a particle.
     * 
     * @param index          the index of the particle for which to set parameters
     * @param charge         the charge of the particle, measured in units of the proton charge
     * @param radius         the GB/VI radius of the particle, measured in nm
     * @param gamma          the gamma parameter
     */
    void setParticleParameters(int index, double charge, double radius, double gamma);
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
protected:
    ForceImpl* createImpl();
private:
    class ParticleInfo;
    double solventDielectric, soluteDielectric;

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

class GBVIForce::ParticleInfo {
public:
    double charge, radius, gamma;
    ParticleInfo() {
        charge = radius = gamma = 0.0;
    }
    ParticleInfo(double charge, double radius, double gamma) :
        charge(charge), radius(radius), gamma(gamma) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_GBVIFORCEFIELD_H_*/
