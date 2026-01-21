#ifndef OPENMM_CAVITYFORCE_H_
#define OPENMM_CAVITYFORCE_H_

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
#include <vector>
#include <utility>
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This class implements the cavity-molecule interaction force for cavity molecular 
 * dynamics simulations. The cavity is modeled as a single fictitious particle (photon)
 * that couples to the molecular dipole moment.
 * 
 * The cavity Hamiltonian is:
 *   H = (1/2) * K * q^2 + epsilon * q . d + (epsilon^2 / 2K) * d^2
 * 
 * where:
 *   - q is the cavity mode position (unwrapped coordinates of the cavity particle)
 *   - d is the molecular dipole moment (sum of charge * position for all molecular particles)
 *   - K = photonMass * omegac^2 is the cavity spring constant (converted to kJ/(mol·nm²))
 *   - epsilon = lambdaCoupling * omegac is the effective coupling strength (converted to kJ/(mol·nm²·e))
 *   - lambdaCoupling is the dimensionless coupling parameter
 *   - omegac is the cavity frequency in atomic units (Hartree)
 *   - photonMass is in amu, converted internally to atomic units (1 amu = 1822.888 mₑ)
 *   
 * Note: All atomic unit parameters are converted internally to OpenMM units (kJ/mol, nm, ps).
 * 
 * The energy is decomposed into three components for analysis:
 *   - Harmonic energy: (1/2) * K * q^2
 *   - Coupling energy: epsilon * q . d (only x,y components)
 *   - Dipole self-energy: (epsilon^2 / 2K) * d^2 (only x,y components)
 * 
 * IMPORTANT: All position calculations use UNWRAPPED coordinates. This means
 * the actual position = wrapped_position + image_flags * box_vectors. This is
 * critical for correct dipole moment calculations across periodic boundaries.
 * 
 * The coupling strength can be time-varying using setLambdaCouplingSchedule().
 * This allows for sudden or gradual switching of the coupling.
 */
class OPENMM_EXPORT CavityForce : public Force {
public:
    /**
     * Parameter name for the cavity frequency (in atomic units).
     */
    static const std::string& OmegaC() {
        static const std::string key = "CavityOmegaC";
        return key;
    }
    /**
     * Parameter name for the dimensionless coupling parameter lambda.
     */
    static const std::string& LambdaCoupling() {
        static const std::string key = "CavityLambdaCoupling";
        return key;
    }
    /**
     * Parameter name for the photon mass (in atomic units).
     */
    static const std::string& PhotonMass() {
        static const std::string key = "CavityPhotonMass";
        return key;
    }
    /**
     * Create a CavityForce.
     * 
     * @param cavityParticleIndex  the index of the cavity (photon) particle
     * @param omegac               the cavity frequency in atomic units (Hartree)
     * @param lambdaCoupling       the dimensionless coupling parameter
     * @param photonMass           the effective photon mass in amu (default: 1.0)
     */
    CavityForce(int cavityParticleIndex, double omegac, double lambdaCoupling, double photonMass = 1.0);
    /**
     * Get the index of the cavity (photon) particle.
     */
    int getCavityParticleIndex() const {
        return cavityParticleIndex;
    }
    /**
     * Set the index of the cavity (photon) particle.
     * 
     * @param index  the particle index
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
     * 
     * @param omegac  the cavity frequency in atomic units (Hartree)
     */
    void setOmegac(double omegac) {
        this->omegac = omegac;
    }
    /**
     * Get the dimensionless coupling parameter lambda.
     */
    double getLambdaCoupling() const {
        return lambdaCoupling;
    }
    /**
     * Set the dimensionless coupling parameter lambda.
     * 
     * @param lambda  the dimensionless coupling parameter
     */
    void setLambdaCoupling(double lambda) {
        lambdaCoupling = lambda;
    }
    /**
     * Get the photon mass (in amu).
     */
    double getPhotonMass() const {
        return photonMass;
    }
    /**
     * Set the photon mass.
     * 
     * @param mass  the photon mass in amu
     */
    void setPhotonMass(double mass) {
        photonMass = mass;
    }
    /**
     * Get the spring constant K = photonMass * omegac^2.
     */
    double getSpringConstant() const {
        return photonMass * omegac * omegac;
    }
    /**
     * Get the effective coupling strength epsilon = lambdaCoupling * omegac.
     */
    double getEffectiveCoupling() const {
        return lambdaCoupling * omegac;
    }
    /**
     * Set the step at which coupling should be switched on and the coupling value.
     * This is a simpler alternative to setLambdaCouplingSchedule.
     * 
     * @param step   the timestep to start coupling
     * @param value  the lambda coupling value to use after that step
     */
    void setCouplingOnStep(int step, double value);
    /**
     * Get the coupling start step.
     */
    int getCouplingOnStep() const {
        return couplingOnStep;
    }
    /**
     * Get the coupling value that will be used when coupling is turned on.
     */
    double getCouplingOnValue() const {
        return couplingOnValue;
    }
    /// @cond INTERNAL
#ifndef SWIG
    /**
     * Set a schedule for time-varying coupling. The coupling value will change
     * at the specified timesteps. Between schedule points, the coupling remains
     * constant (step function).
     * 
     * @param schedule  vector of (timestep, lambda_value) pairs, sorted by timestep
     */
    void setLambdaCouplingSchedule(const std::vector<std::pair<int, double>>& schedule);
    /**
     * Get the coupling schedule.
     */
    const std::vector<std::pair<int, double>>& getLambdaCouplingSchedule() const {
        return couplingSchedule;
    }
#endif
    /// @endcond
    /**
     * Get the lambda coupling value at a specific timestep.
     * If no schedule is set, returns the default lambdaCoupling value.
     * 
     * @param step  the timestep
     * @return the lambda coupling value at that timestep
     */
    double getLambdaCouplingAtStep(int step) const;
    /**
     * Get the harmonic energy component: (1/2) * K * q^2
     * This must be called after forces have been computed.
     * 
     * @param context  the Context to query
     * @return the harmonic energy in kJ/mol
     */
    double getHarmonicEnergy(const Context& context) const;
    /**
     * Get the coupling energy component: epsilon * q . d
     * This must be called after forces have been computed.
     * 
     * @param context  the Context to query
     * @return the coupling energy in kJ/mol
     */
    double getCouplingEnergy(const Context& context) const;
    /**
     * Get the dipole self-energy component: (epsilon^2 / 2K) * d^2
     * This must be called after forces have been computed.
     * 
     * @param context  the Context to query
     * @return the dipole self-energy in kJ/mol
     */
    double getDipoleSelfEnergy(const Context& context) const;
    /**
     * Get the total cavity energy (sum of all three components).
     * 
     * @param context  the Context to query
     * @return the total cavity energy in kJ/mol
     */
    double getTotalCavityEnergy(const Context& context) const;
    /**
     * Update the parameters in a Context to match those stored in this Force object.
     * 
     * @param context  the Context to update
     */
    void updateParametersInContext(Context& context);
    /**
     * Returns whether this force uses periodic boundary conditions.
     * 
     * @return true (cavity force uses PBC for unwrapped position calculations)
     */
    bool usesPeriodicBoundaryConditions() const {
        return true;
    }
protected:
    ForceImpl* createImpl() const;
private:
    int cavityParticleIndex;
    double omegac;
    double lambdaCoupling;
    double photonMass;
    int couplingOnStep;
    double couplingOnValue;
    std::vector<std::pair<int, double>> couplingSchedule;
};

} // namespace OpenMM

#endif /*OPENMM_CAVITYFORCE_H_*/
