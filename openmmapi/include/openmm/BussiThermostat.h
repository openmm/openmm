#ifndef OPENMM_BUSSITHERMOSTAT_H_
#define OPENMM_BUSSITHERMOSTAT_H_

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
#include <string>
#include <set>
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This class implements the Bussi-Donadio-Parrinello stochastic velocity rescaling 
 * thermostat (also known as the "Bussi thermostat" or "canonical sampling through 
 * velocity rescaling" - CSVR).
 * 
 * Unlike simple velocity rescaling, this method correctly samples the canonical 
 * ensemble by introducing stochastic noise that satisfies detailed balance. It can
 * be applied to a subset of particles, allowing different thermostats to be used 
 * for different parts of the system.
 * 
 * The algorithm rescales velocities according to:
 *   v_new = alpha * v_old
 * 
 * where alpha is computed stochastically to maintain the canonical distribution:
 *   alpha^2 = c + (1-c) * (sum_i R_i^2 + 2*R_1*sqrt(c*K/((1-c)*K_target)))
 * 
 * where c = exp(-dt/tau), K is the current kinetic energy, K_target is the target
 * kinetic energy, and R_i are independent normal random variates.
 * 
 * The rescaling factor alpha may be positive or negative (signed alpha per
 * Bussi et al. 2009 Eq. A8) for correct canonical sampling. When used with
 * VerletIntegrator, the thermostat is applied after the first half-kick
 * (HOOMD-style order). Zero kinetic energy in the thermostatted group throws.
 *
 * The thermostat tracks the cumulative energy transferred to/from the reservoir,
 * which is useful for computing work and heat in non-equilibrium simulations.
 *
 * References:
 *   Bussi, Donadio, and Parrinello, J. Chem. Phys. 126, 014101 (2007)
 *   Bussi, Zykova-Timan, and Parrinello, J. Chem. Phys. 130, 074101 (2009)
 */
class OPENMM_EXPORT BussiThermostat : public Force {
public:
    /**
     * This is the name of the parameter which stores the current temperature of the
     * heat bath (in Kelvin).
     */
    static const std::string& Temperature() {
        static const std::string key = "BussiTemperature";
        return key;
    }
    /**
     * This is the name of the parameter which stores the current time constant tau (in ps).
     */
    static const std::string& Tau() {
        static const std::string key = "BussiTau";
        return key;
    }
    /**
     * This is the name of the parameter which stores the cumulative reservoir energy
     * from translational degrees of freedom (in kJ/mol).
     */
    static const std::string& ReservoirEnergyTranslational() {
        static const std::string key = "BussiReservoirEnergyTranslational";
        return key;
    }
    /**
     * This is the name of the parameter which stores the cumulative reservoir energy
     * from rotational degrees of freedom (in kJ/mol).
     */
    static const std::string& ReservoirEnergyRotational() {
        static const std::string key = "BussiReservoirEnergyRotational";
        return key;
    }
    /**
     * Create a BussiThermostat.
     * 
     * @param defaultTemperature  the default temperature of the heat bath (in Kelvin)
     * @param defaultTau          the default time constant (in picoseconds)
     */
    BussiThermostat(double defaultTemperature, double defaultTau);
    /**
     * Get the default temperature of the heat bath (in Kelvin).
     *
     * @return the default temperature of the heat bath, measured in Kelvin.
     */
    double getDefaultTemperature() const {
        return defaultTemp;
    }
    /**
     * Set the default temperature of the heat bath. This will affect any new Contexts you create,
     * but not ones that already exist.
     *
     * @param temperature  the default temperature of the heat bath (in Kelvin)
     */
    void setDefaultTemperature(double temperature) {
        defaultTemp = temperature;
    }
    /**
     * Get the default time constant tau (in picoseconds).
     *
     * @return the default time constant, measured in picoseconds.
     */
    double getDefaultTau() const {
        return defaultTau;
    }
    /**
     * Set the default time constant. This will affect any new Contexts you create,
     * but not ones that already exist.
     *
     * @param tau  the default time constant (in picoseconds)
     */
    void setDefaultTau(double tau) {
        defaultTau = tau;
    }
    /**
     * Get the number of particles this thermostat is applied to.
     */
    int getNumParticles() const {
        return particles.size();
    }
    /**
     * Add a particle to the thermostat.
     * 
     * @param particle  the index of the particle to add
     * @return the index of the particle that was added
     */
    int addParticle(int particle);
    /**
     * Get the index of a particle this thermostat is applied to.
     * 
     * @param index  the index of the particle within the thermostat's particle list
     * @return the index of the particle within the System
     */
    int getParticle(int index) const;
    /// @cond INTERNAL
#ifndef SWIG
    /**
     * Set the particle indices this thermostat applies to.
     * 
     * @param particles  the set of particle indices
     */
    void setParticles(const std::set<int>& particles);
    /**
     * Get the set of particle indices this thermostat applies to.
     * 
     * @return the set of particle indices
     */
    const std::set<int>& getParticles() const {
        return particles;
    }
#endif
    /// @endcond
    /**
     * Get whether this thermostat applies to all particles. If true, the particle
     * list is ignored and the thermostat applies to all particles in the system
     * (except those with zero mass).
     */
    bool getApplyToAllParticles() const {
        return applyToAll;
    }
    /**
     * Set whether this thermostat applies to all particles. If true, the particle
     * list is ignored and the thermostat applies to all particles in the system
     * (except those with zero mass).
     * 
     * @param apply  whether to apply to all particles
     */
    void setApplyToAllParticles(bool apply) {
        applyToAll = apply;
    }
    /**
     * Get the random number seed. See setRandomNumberSeed() for details.
     */
    int getRandomNumberSeed() const {
        return randomNumberSeed;
    }
    /**
     * Set the random number seed. The precise meaning of this parameter is undefined, 
     * and is left up to each Platform to interpret in an appropriate way. It is 
     * guaranteed that if two simulations are run with different random number seeds, 
     * the sequence of velocity rescalings will be different.
     *
     * If seed is set to 0 (which is the default value assigned), a unique seed is 
     * chosen when a Context is created from this Force.
     */
    void setRandomNumberSeed(int seed) {
        randomNumberSeed = seed;
    }
    /**
     * Get the frequency (in time steps) at which the thermostat is applied. Default is 1.
     */
    int getFrequency() const {
        return frequency;
    }
    /**
     * Set the frequency (in time steps) at which the thermostat is applied.
     * 
     * @param freq  the frequency (1 = every step, 2 = every other step, etc.)
     */
    void setFrequency(int freq) {
        frequency = freq;
    }
    /**
     * Returns whether or not this force makes use of periodic boundary conditions.
     *
     * @returns false - this force does not use PBC
     */
    bool usesPeriodicBoundaryConditions() const {
        return false;
    }
protected:
    ForceImpl* createImpl() const;
private:
    double defaultTemp;
    double defaultTau;
    std::set<int> particles;
    bool applyToAll;
    int randomNumberSeed;
    int frequency;
};

} // namespace OpenMM

#endif /*OPENMM_BUSSITHERMOSTAT_H_*/
