#ifndef OPENMM_NOSEHOOVERCHAIN_H_
#define OPENMM_NOSEHOOVERCHAIN_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2019 Stanford University and the Authors.           *
 * Authors: Andreas Krämer and Andrew C. Simmonett                            *
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
#include <vector>
#include "openmm/OpenMMException.h"
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This class defines a chain of Nose-Hoover particles to be used as a heat bath to
 * scale the velocities of a collection of particles subject to thermostating.  The
 * heat bath is propagated using the multi time step approach detailed in
 *
 * G. J. Martyna, M. E. Tuckerman, D. J. Tobias and M. L. Klein, Mol. Phys. 87, 1117 (1996).
 *
 * where the total number of timesteps used to propagate the chain in each step is
 * the number of MTS steps multiplied by the number of terms in the Yoshida-Suzuki decomposition.
 *
 * Two types of NHC may be created.  The first is a simple thermostat that couples with a given subset
 * of the atoms within a system, controling their absolute motion.  The second is more elaborate and
 * can thermostat tethered pairs of atoms and in this case two thermostats are created: one that controls
 * the absolute center of mass velocity of each pair and another that controls their motion relative to
 * one another.
 */

class OPENMM_EXPORT NoseHooverChain {
public:
    /**
     * Create a NoseHooverChain.
     * 
     * @param defaultTemperature                the default temperature of the heat bath for absolute motion (in Kelvin)
     * @param defaultRelativeTemperature        the default temperature of the heat bath for relative motion(in Kelvin).
     *                                          This is only used if the list of thermostated pairs is not empty.
     * @param defaultCollisionFrequency         the default collision frequency for absolute motion (in 1/ps)
     * @param defaultRelativeCollisionFrequency the default collision frequency for relative motion(in 1/ps).
     *                                          This is only used if the list of thermostated pairs is not empty.
     * @param defaultNumDOFs                    the default number of degrees of freedom in the particles that
     *                                          interact with this chain
     * @param defaultChainLength                the default length of (number of particles in) this heat bath
     * @param defaultNumMTS                     the default number of multi time steps used to propagate this chain
     * @param defaultNumYoshidaSuzuki           the default number of Yoshida Suzuki steps used to propagate this chain (1, 3, or 5).
     * @param defaultChainID                    the default chain id used to distinguish this Nose-Hoover chain from others that may
     *                                          be used to control a different set of particles, e.g. for Drude oscillators
     * @param thermostatedAtoms                 the list of atoms to be handled by this thermostat
     * @param thermostatedPairs                 the list of connected pairs to be thermostated; their absolute center of mass motion will
     *                                          be thermostated independently from their motion relative to one another.
     */
    NoseHooverChain(double defaultTemperature, double defaultRelativeTemperature, double defaultCollisionFrequency,
                    double defaultRelativeCollisionFrequency, int defaultNumDOFs, int defaultChainLength,
                    int defaultNumMTS, int defaultNumYoshidaSuzuki, int defaultChainID,
                    const std::vector<int>& thermostatedAtoms, const std::vector< std::pair< int, int > > &thermostatedPairs);
    /**
     * Get the default temperature of the heat bath for treating absolute particle motion (in Kelvin).
     *
     * @return the default temperature of the heat bath, measured in Kelvin.
     */
    double getDefaultTemperature() const {
        return defaultTemp;
    }
    /**
     * Set the default temperature of the heat bath for treating absolute particle motion.
     * This will affect any new Contexts you create, but not ones that already exist.
     *
     * @param temperature the default temperature of the heat bath (in Kelvin)
     */
    void setDefaultTemperature(double temperature) {
        defaultTemp = temperature;
    }
    /**
     * Get the default temperature of the heat bath for treating relative particle motion (in Kelvin).
     *
     * @return the default temperature of the heat bath, measured in Kelvin.
     */
    double getDefaultRelativeTemperature() const {
        return defaultRelativeTemp;
    }
    /**
     * Set the default temperature of the heat bath for treating relative motion if this thermostat has
     * been set up to treat connected pairs of atoms.  This will affect any new Contexts you create,
     * but not ones that already exist.
     *
     * @param temperature the default temperature of the heat bath for relative motion (in Kelvin)
     */
    void setDefaultRelativeTemperature(double temperature) {
        defaultRelativeTemp = temperature;
    }
    /**
     * Get the default collision frequency for treating absolute particle motion (in 1/ps).
     *
     * @return the default collision frequency, measured in 1/ps.
     */
    double getDefaultCollisionFrequency() const {
        return defaultFreq;
    }
    /**
     * Set the default collision frequency for treating absolute particle motion.
     * This will affect any new Contexts you create, but not those that already exist.
     *
     * @param frequency the default collision frequency (in 1/ps)
     */
    void setDefaultCollisionFrequency(double frequency) {
        defaultFreq = frequency;
    }
    /**
     * Get the default collision frequency for treating relative particle motion (in 1/ps).
     *
     * @return the default collision frequency, measured in 1/ps.
     */
    double getDefaultRelativeCollisionFrequency() const {
        return defaultRelativeFreq;
    }
    /**
     * Set the default collision frequency for treating relative particle motion if this thermostat has
     * been set up to handle connected pairs of atoms.  This will affect any new Contexts you create,
     * but not those that already exist.
     *
     * @param frequency the default collision frequency (in 1/ps)
     */
    void setDefaultRelativeCollisionFrequency(double frequency) {
        defaultRelativeFreq = frequency;
    }
    /**
     * Get the default number of degrees of freedom in the particles controled by this heat bath.
     *
     * @return the default number of degrees of freedom.
     */
    int getDefaultNumDegreesOfFreedom() const {
        return defaultNumDOFs;
    }
    /**
     * Set the default number of degrees of freedom in the particles controled by this heat bath.
     * This will affect any new Contexts you create, but not those that already exist.
     *
     * @param numDOFs 
     */
    void setDefaultNumDegreesOfFreedom(int numDOFs) {
        defaultNumDOFs = numDOFs;
    }
    /**
     * Get the default chain length of this heat bath.
     *
     * @return the default chain length.
     */
    int getDefaultChainLength() const {
        return defaultChainLength;
    }
    /**
     * Set the default chain length of this heat bath.
     * This will affect any new Contexts you create, but not those that already exist.
     *
     * @param numDOFs 
     */
    void setDefaultChainLength(int chainLength) {
        defaultChainLength = chainLength;
    }
    /**
     * Get the default number of steps used in the multi time step propagation.
     *
     * @returns the default number of multi time steps.
     */
    int getDefaultNumMultiTimeSteps() const {
        return defaultNumMTS;
    }
    /**
     * Set the default number of steps used in the multi time step propagation.
     * This will affect any new Contexts you create, but not those that already exist.
     *
     * @param numSteps 
     */
    void setDefaultNumMultiTimeSteps(int numSteps) {
        defaultNumMTS = numSteps;
    }
    /**
     * Get the default number of steps used in the Yoshida-Suzuki decomposition for
     * multi time step propagation.
     *
     * @returns the default number of multi time steps in the Yoshida-Suzuki decomposition.
     */
    int getDefaultNumYoshidaSuzukiTimeSteps() const {
        return defaultNumYS;
    }
    /**
     * Set the default number of steps used in the Yoshida-Suzuki decomposition for
     * multi time step propagation. This will affect any new Contexts you create,
     * but not those that already exist.
     * 
     * @param the default number of multi time steps in the Yoshida-Suzuki decomposition.
     */
    void setDefaultNumYoshidaSuzukiTimeSteps(int numSteps) {
        defaultNumYS = numSteps;
    }
    /**
     * Get the default chain id used to identify this chain
     *
     * @returns the default chain id 
     */
    int getDefaultChainID() const {
        return defaultChainID;
    }
    /**
     * Set the default chain id used to identify this chain
     * This will affect any new Contexts you create, but not those that already exist.
     *
     * @param the default chain id
     */
    void setDefaultChainID(int chainID) {
        defaultChainID = chainID;
    }
    /**
     * Get the atom ids of all atoms that are thermostated 
     *
     * @returns ids of all atoms that are being handled by this thermostat
     */
    const std::vector<int>& getThermostatedAtoms() const {
        return thermostatedAtoms;
    }
    /**
     * Set list of atoms that are handled by this thermostat
     *
     * @param atomIDs 
     */
    void setThermostatedAtoms(const std::vector<int>& atomIDs){
        thermostatedAtoms = atomIDs;
    }
    /**
     * Get the list of any connected pairs to be handled by this thermostat.
     * If this is a regular thermostat, returns an empty vector.
     *
     * @returns list of connected pairs.
     */
    const std::vector< std::pair< int, int > >& getThermostatedPairs() const {
        return thermostatedPairs;
    }
    /**
     * In case this thermostat handles the kinetic energy of Drude particles 
     * set the atom IDs of all parent atoms. 
     *
     * @param pairIDs the list of connected pairs to thermostat.
     */
    void setThermostatedPairs(const std::vector< std::pair< int, int > >& pairIDs){
        thermostatedPairs = pairIDs;
    }
    /**
     * Get the default weights used in the Yoshida Suzuki multi time step decomposition (dimensionless) 
     *
     * @returns the weights for the Yoshida-Suzuki integration
     */
    std::vector<double> getDefaultYoshidaSuzukiWeights() const;
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * @returns true if force uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const {
        return false;
    }
private:
    double defaultTemp, defaultFreq, defaultRelativeTemp, defaultRelativeFreq;
    int defaultNumDOFs, defaultChainLength, defaultNumMTS, defaultNumYS;
    // The suffix used to distinguish NH chains, e.g. for Drude particles vs. regular particles.
    int defaultChainID;
    std::vector<int> thermostatedAtoms;
    std::vector<std::pair<int, int>> thermostatedPairs;
};

} // namespace OpenMM

#endif /*OPENMM_NOSEHOOVERCHAIN_H_*/