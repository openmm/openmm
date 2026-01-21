#ifndef OPENMM_CAVITYPARTICLEDISPLACER_H_
#define OPENMM_CAVITYPARTICLEDISPLACER_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
 * Authors: Muhammad Hasyim                                                   *
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
#include "Vec3.h"
#include <string>
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This class implements the finite-Q displacement for cavity molecular dynamics.
 * When the cavity coupling is switched on, the cavity particle needs to be displaced
 * to its equilibrium position to avoid sudden energy injection into the system.
 * 
 * The equilibrium position is calculated as:
 *   q_eq = -(lambda / omega_c) * d_xy
 * 
 * where:
 *   - lambda is the dimensionless coupling parameter
 *   - omega_c is the cavity frequency
 *   - d_xy is the molecular dipole moment (only x,y components)
 * 
 * This displacement preserves the z-coordinate of the cavity particle and does not
 * modify the velocity, ensuring smooth energy transitions when coupling is activated.
 * 
 * IMPORTANT: This force should be used in conjunction with CavityForce when using
 * time-varying coupling schedules.
 */
class OPENMM_EXPORT CavityParticleDisplacer : public Force {
public:
    /**
     * Create a CavityParticleDisplacer.
     * 
     * @param cavityParticleIndex  the index of the cavity (photon) particle
     * @param omegac               the cavity frequency in atomic units
     * @param photonMass           the effective photon mass in atomic units (default: 1.0)
     */
    CavityParticleDisplacer(int cavityParticleIndex, double omegac, double photonMass = 1.0);
    /**
     * Get the index of the cavity (photon) particle.
     */
    int getCavityParticleIndex() const {
        return cavityParticleIndex;
    }
    /**
     * Set the index of the cavity (photon) particle.
     */
    void setCavityParticleIndex(int index) {
        cavityParticleIndex = index;
    }
    /**
     * Get the cavity frequency omega_c (in atomic units).
     */
    double getOmegac() const {
        return omegac;
    }
    /**
     * Set the cavity frequency omega_c.
     */
    void setOmegac(double omegac) {
        this->omegac = omegac;
    }
    /**
     * Get the photon mass (in atomic units).
     */
    double getPhotonMass() const {
        return photonMass;
    }
    /**
     * Set the photon mass.
     */
    void setPhotonMass(double mass) {
        photonMass = mass;
    }
    /**
     * Get the timestep at which the coupling will be switched on.
     * The displacer will trigger at this timestep.
     */
    int getSwitchOnStep() const {
        return switchOnStep;
    }
    /**
     * Set the timestep at which the coupling will be switched on.
     * 
     * @param step  the timestep number
     */
    void setSwitchOnStep(int step) {
        switchOnStep = step;
    }
    /**
     * Get the lambda coupling value that will be applied when switched on.
     */
    double getSwitchOnLambda() const {
        return switchOnLambda;
    }
    /**
     * Set the lambda coupling value for the displacement calculation.
     * 
     * @param lambda  the dimensionless coupling parameter
     */
    void setSwitchOnLambda(double lambda) {
        switchOnLambda = lambda;
    }
    /**
     * Get whether the displacement has already been triggered.
     */
    bool hasTriggered() const {
        return triggered;
    }
    /**
     * Manually displace the cavity particle to its equilibrium position.
     * This is typically called internally when the switch-on step is reached.
     * 
     * @param context        the Context to modify
     * @param lambdaCoupling the coupling value to use for displacement calculation
     */
    void displaceToEquilibrium(Context& context, double lambdaCoupling);
    /**
     * Returns whether this force uses periodic boundary conditions.
     */
    bool usesPeriodicBoundaryConditions() const {
        return true;
    }
protected:
    ForceImpl* createImpl() const;
private:
    int cavityParticleIndex;
    double omegac;
    double photonMass;
    int switchOnStep;
    double switchOnLambda;
    mutable bool triggered;
};

} // namespace OpenMM

#endif /*OPENMM_CAVITYPARTICLEDISPLACER_H_*/
