#ifndef OPENMM_LCPOFORCE_H_
#define OPENMM_LCPOFORCE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2008-2025 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Evan Pretti                                        *
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
#include "internal/windowsExport.h"
#include <vector>

namespace OpenMM {

/**
 * This class implements the functional form of the potential used in the LCPO
 * (linear combinations of pairwise overlaps) method for estimating
 * solvent-accessible surface areas.  Specifically, it implements the three-body
 * potential described in Weiser, Shenkin, and Still, J. Comput. Chem. 20, 217
 * (1999).
 *
 * To use this class, create an LCPOForce object, then call addParticle() once
 * for each particle in the System to define its parameters.  The number of
 * particles for which you define nonbonded parameters must be exactly equal to
 * the number of particles in the System, or else an exception will be thrown
 * when you try to create a Context.  After a particle has been added, you can
 * modify its parameters by calling setParticleParameters().  This will have no
 * effect on Contexts that already exist unless you call
 * updateParametersInContext().
 *
 * To exclude a particle from the force entirely, you can either set its radius
 * to zero or set all of its LCPO parameters to zero (or both).
 */

class OPENMM_EXPORT LCPOForce : public Force {
public:
    /**
     * Create an LCPOForce.
     */
    LCPOForce();
    /**
     * Set whether this force should apply periodic boundary conditions when calculating distances between particles.
     */
    void setUsesPeriodicBoundaryConditions(bool periodic);
    /**
     * Get the number of particles for which force field parameters have been defined.
     */
    int getNumParticles() const;
    /**
     * Add parameters for a particle.  This should be called once for each particle in the System.
     * When it is called for the i'th time, it specifies the parameters for the i'th particle.
     *
     * @param radius  the radius of the particle for the LCPO method, including the solvent probe radius if applicable
     * @param p1      the product of the surface energy with the value of the LCPO parameter P1 for the particle
     * @param p2      the product of the surface energy with the value of the LCPO parameter P2 for the particle
     * @param p3      the product of the surface energy with the value of the LCPO parameter P3 for the particle
     * @param p4      the product of the surface energy with the value of the LCPO parameter P4 for the particle
     * @return        the index of the particle that was added
     */
    int addParticle(double radius, double p1, double p2, double p3, double p4);
    /**
     * Get the nonbonded force parameters for a particle.
     *
     * @param index        the index of the particle for which to get parameters
     * @param[out] radius  the radius of the particle for the LCPO method, including the solvent probe radius if applicable
     * @param[out] p1      the product of the surface energy with the value of the LCPO parameter P1 for the particle
     * @param[out] p2      the product of the surface energy with the value of the LCPO parameter P2 for the particle
     * @param[out] p3      the product of the surface energy with the value of the LCPO parameter P3 for the particle
     * @param[out] p4      the product of the surface energy with the value of the LCPO parameter P4 for the particle
     */
    void getParticleParameters(int index, double& radius, double& p1, double& p2, double& p3, double& p4) const;
    /**
     * Set the nonbonded force parameters for a particle.
     *
     * @param index   the index of the particle for which to set parameters
     * @param radius  the radius of the particle for the LCPO method, including the solvent probe radius if applicable
     * @param p1      the product of the surface energy with the value of the LCPO parameter P1 for the particle
     * @param p2      the product of the surface energy with the value of the LCPO parameter P2 for the particle
     * @param p3      the product of the surface energy with the value of the LCPO parameter P3 for the particle
     * @param p4      the product of the surface energy with the value of the LCPO parameter P4 for the particle
     */
    void setParticleParameters(int index, double radius, double p1, double p2, double p3, double p4);
    /**
     * Update the parameters in a Context to match those stored in this Force
     * object.  This method provides an efficient method to update certain
     * parameters in an existing Context without needing to reinitialize it.
     * Simply call setParticleParameters() to modify this object's parameters,
     * then call updateParametersInContext() to copy them over to the Context.
     * Only the values of particle parameters can be updated by this method;
     * other changes made will not be reflected in the Context after calling it.
     */
    void updateParametersInContext(Context& context);
    /**
     * Returns whether or not this force makes use of periodic boundary conditions.
     *
     * @returns true if force uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const;
protected:
    ForceImpl* createImpl() const;
private:
    class ParticleInfo;
    bool usePeriodic;
    std::vector<ParticleInfo> particles;
    mutable int numContexts, firstChangedParticle, lastChangedParticle;
};

/**
 * This is an internal class used to record information about a particle.
 * @private
 */
class LCPOForce::ParticleInfo {
public:
    double radius, p1, p2, p3, p4;
    ParticleInfo(double radius, double p1, double p2, double p3, double p4) : radius(radius), p1(p1), p2(p2), p3(p3), p4(p4) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_LCPOFORCE_H_*/
