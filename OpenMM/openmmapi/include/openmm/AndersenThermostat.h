#ifndef OPENMM_ANDERSENTHERMOSTAT_H_
#define OPENMM_ANDERSENTHERMOSTAT_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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
#include <string>
#include "internal/windowsExport.h"
#include "openmm/internal/OSRngSeed.h"

namespace OpenMM {

/**
 * This class uses the Andersen method to maintain constant temperature.
 */

class OPENMM_EXPORT AndersenThermostat : public Force {
public:
    /**
     * This is the name of the parameter which stores the current temperature of the
     * heat bath (in Kelvin).
     */
    static const std::string& Temperature() {
        static const std::string key = "AndersenTemperature";
        return key;
    }
    /**
     * This is the name of the parameter which store the current collision frequency (in 1/ps).
     */
    static const std::string& CollisionFrequency() {
        static const std::string key = "AndersenCollisionFrequency";
        return key;
    }
    /**
     * Create an AndersenThermostat.
     * 
     * @param defaultTemperature         the default temperature of the heat bath (in Kelvin)
     * @param defaultCollisionFrequency  the default collision frequency (in 1/ps)
     */
    AndersenThermostat(double defaultTemperature, double defaultCollisionFrequency);
    /**
     * Get the default temperature of the heat bath (in Kelvin).
     *
     * @return the default temperature of the heat bath, measured in Kelvin.
     */
    double getDefaultTemperature() const {
        return defaultTemp;
    }
    /**
     * Set the default temperature of the heat bath.  This will affect any new Contexts you create,
     * but not ones that already exist.
     *
     * @param temperature         the default temperature of the heat bath (in Kelvin)
     */
    void setDefaultTemperature(double temperature) {
        defaultTemp = temperature;
    }
    /**
     * Get the default collision frequency (in 1/ps).
     *
     * @return the default collision frequency, measured in 1/ps.
     */
    double getDefaultCollisionFrequency() const {
        return defaultFreq;
    }
    /**
     * Set the default collision frequency.  This will affect any new Contexts you create,
     * but not ones that already exist.
     *
     * @param frequency         the default collision frequency (in 1/ps)
     */
    void setDefaultCollisionFrequency(double frequency) {
        defaultFreq = frequency;
    }
    /**
     * Get the random number seed.  See setRandomNumberSeed() for details.
     */
    int getRandomNumberSeed() const {
        return randomNumberSeed;
    }
    /**
     * Set the random number seed.  The precise meaning of this parameter is undefined, and is left up
     * to each Platform to interpret in an appropriate way.  It is guaranteed that if two simulations
     * are run with different random number seeds, the sequence of collisions will be different.  On
     * the other hand, no guarantees are made about the behavior of simulations that use the same seed.
     * In particular, Platforms are permitted to use non-deterministic algorithms which produce different
     * results on successive runs, even if those runs were initialized identically.
     *
     * If seed is set to 0 (which is the default value assigned), a unique seed is chosen when a Context
     * is created from this Force. This is done to ensure that each Context receives unique random seeds
     * without you needing to set them explicitly.
     */
    void setRandomNumberSeed(int seed) {
        randomNumberSeed = seed;
    }
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * @returns true if force uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const {
        return false;
    }
protected:
    ForceImpl* createImpl() const;
private:
    double defaultTemp, defaultFreq;
    int randomNumberSeed;
};

} // namespace OpenMM

#endif /*OPENMM_ANDERSENTHERMOSTAT_H_*/
