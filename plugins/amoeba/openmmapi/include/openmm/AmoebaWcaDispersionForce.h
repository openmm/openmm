#ifndef OPENMM_AMOEBA_WCA_DISPERSION_FORCE_H_
#define OPENMM_AMOEBA_WCA_DISPERSION_FORCE_H_

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
 * This class implements a nonbonded interaction between pairs of particles typically used along with
 * AmoebaGeneralizedKirkwoodForce as part of an implicit solvent model.
 * 
 * To use it, create an AmoebaWcaDispersionForce object then call addParticle() once for each particle.  After
 * a particle has been added, you can modify its force field parameters by calling setParticleParameters().
 * This will have no effect on Contexts that already exist unless you call updateParametersInContext().
 */

class OPENMM_EXPORT_AMOEBA AmoebaWcaDispersionForce : public Force {

public:

    /**
     * Create an AmoebaWcaDispersionForce.
     */
    AmoebaWcaDispersionForce();

    /**
     * Get the number of particles
     */
    int getNumParticles() const {
        return parameters.size();
    }

    /**
     * Set the force field parameters for a WCA dispersion particle.
     * 
     * @param particleIndex   the particle index
     * @param radius          radius
     * @param epsilon         epsilon 
     */
    void setParticleParameters(int particleIndex, double radius, double epsilon);

    /**
     * Get the force field parameters for a WCA dispersion particle.
     * 
     * @param particleIndex        the particle index
     * @param[out] radius          radius
     * @param[out] epsilon         epsilon 
     */
    void getParticleParameters(int particleIndex, double& radius, double& epsilon) const;

    /**
     * Set the force field parameters for a WCA dispersion particle.
     * 
     * @param radius          radius
     * @param epsilon         epsilon 
     * @return index of added particle
     */
    int addParticle(double radius, double epsilon);
    /**
     * Update the per-particle parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setParticleParameters() to modify this object's parameters, then call updateParametersInContext()
     * to copy them over to the Context.
     * 
     * The only information this method updates is the values of per-particle parameters.  All other aspects of the Force
     * are unaffected and can only be changed by reinitializing the Context.
     */
    void updateParametersInContext(Context& context);

    /* 
     * Constants
     */

    double getEpso() const;
    double getEpsh() const;
    double getRmino() const;
    double getRminh() const;
    double getAwater() const;
    double getShctd() const;
    double getDispoff() const;
    double getSlevy() const;

    void setEpso(double inputValue);
    void setEpsh(double inputValue);
    void setRmino(double inputValue);
    void setRminh(double inputValue);
    void setAwater(double inputValue);
    void setShctd(double inputValue);
    void setDispoff(double inputValue);
    void setSlevy(double inputValue);
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
    class WcaDispersionInfo;
    double epso;
    double epsh;
    double rmino;
    double rminh;
    double awater;
    double slevy;
    double shctd;
    double dispoff;
    std::vector<WcaDispersionInfo> parameters;
};

/**
 * This is an internal class used to record information about a particle.
 * @private
 */
class AmoebaWcaDispersionForce::WcaDispersionInfo {
public:
    double radius, epsilon;
    WcaDispersionInfo() {
        radius              = 1.0;
        epsilon             = 0.0;
    }
    WcaDispersionInfo(double radius, double epsilon) : radius(radius), epsilon(epsilon) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_AMOEBA_WCA_DISPERSION_FORCE_H_*/
