#ifndef OPENMM_AMOEBA_WCA_DISPERSION_FORCE_H_
#define OPENMM_AMOEBA_WCA_DISPERSION_FORCE_H_

/* -------------------------------------------------------------------------- *
 *                              AmoebaOpenMM                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
 * Authors:                                                                   *
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
 * This class implements an interaction between pairs of particles that varies harmonically with the distance
 * between them.  To use it, create a WcaDispersionForce object then call addAngle() once for each angle.  After
 * a angle has been added, you can modify its force field parameters by calling setAngleParameters().
 */

 class OPENMM_EXPORT AmoebaWcaDispersionForce : public Force {
public:
    /**
     * Create a Amoeba WcaDispersionForce.
     */
    AmoebaWcaDispersionForce();

    /**
     * Get the number of particles
     */
    int getNumParticles() const {
        return parameters.size();
    }

    /**
     * Set the force field parameters for a wca dispersion particle.
     * 
     * @param particleIndex   the particle index
     * @param radius          radius
     * @param epsilon         epsilon 
     */
    void setParticleParameters(int particleIndex, double radius, double epsilon );

    /**
     * Get the force field parameters for a wca dispersion particle.
     * 
     * @param particleIndex   the particle index
     * @param radius          radius
     * @param epsilon         epsilon 
     */
    void getParticleParameters(int particleIndex, double& radius, double& epsilon) const;

    /**
     * Set the force field parameters for a wca dispersion particle.
     * 
     * @param radius          radius
     * @param epsilon         epsilon 
     * @return index of added particle
     */
    int addParticle( double radius, double epsilon );

    /* 
     * Constants
     */

    double getEpso(     void ) const;
    double getEpsh(     void ) const;
    double getRmino(    void ) const;
    double getRminh(    void ) const;
    double getAwater(   void ) const;
    double getShctd(    void ) const;
    double getDispoff(  void ) const;
    double getSlevy(    void ) const;

    void setEpso(     double inputValue );
    void setEpsh(     double inputValue );
    void setRmino(    double inputValue );
    void setRminh(    double inputValue );
    void setAwater(   double inputValue );
    void setShctd(    double inputValue );
    void setDispoff(  double inputValue );
    void setSlevy(    double inputValue );

protected:
    ForceImpl* createImpl();
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

// Retarded visual studio compiler complains about being unable to 
// export private stl class members.
// This stanza explains that it should temporarily shut up.
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4251)
#endif
    std::vector<WcaDispersionInfo> parameters;
#if defined(_MSC_VER)
#pragma warning(pop)
#endif

};

class AmoebaWcaDispersionForce::WcaDispersionInfo {
public:
    double radius, epsilon;
    WcaDispersionInfo() {
        radius              = 1.0;
        epsilon             = 0.0;
    }
    WcaDispersionInfo( double radius, double epsilon ) : radius(radius), epsilon(epsilon) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_AMOEBA_WCA_DISPERSION_FORCE_H_*/
