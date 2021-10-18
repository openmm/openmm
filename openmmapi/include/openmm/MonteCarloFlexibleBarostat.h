#ifndef OPENMM_MONTECARLOFLEXIBLEBAROSTAT_H_
#define OPENMM_MONTECARLOFLEXIBLEBAROSTAT_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2021 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Sander Vandenhaute                                 *
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

namespace OpenMM {

/**
 * This class uses a Monte Carlo algorithm to adjust the size and shape of the periodic box, simulating the
 * effect of constant pressure.
 *
 * This class is similar to MonteCarloBarostat, but it simulates a fully flexible periodic box
 * in which all three lengths and all three angles are free to change independently.  It is primarily
 * useful for simulations of bulk materials where the shape of a crystal's unit cell may not
 * be known in advance, or could even change with time as it transitions between phases.
 * 
 * Like MonteCarloBarostat, the default behavior of this class is to scale the centroid position
 * of each molecule while holding it rigid.  In simulations of materials where all atoms are
 * covalently bonded to each other, this behavior will not work well since the entire system
 * then consists of a single molecule.  You can use setScaleMoleculesAsRigid() to disable this
 * behavior and instead have it scale the position of every atom independently.
 * 
 * This class assumes the simulation is also being run at constant temperature, and requires you
 * to specify the system temperature (since it affects the acceptance probability for Monte Carlo
 * moves).  It does not actually perform temperature regulation, however.  You must use another
 * mechanism along with it to maintain the temperature, such as LangevinIntegrator or AndersenThermostat.
 */

class OPENMM_EXPORT MonteCarloFlexibleBarostat : public Force {
public:
    /**
     * This is the name of the parameter which stores the current pressure acting on
     * the system (in bar).
     */
    static const std::string& Pressure() {
        static const std::string key = "MonteCarloPressure";
        return key;
    }
    /**
     * This is the name of the parameter which stores the current temperature at which the
     * system is being maintained (in Kelvin)
     */
    static const std::string& Temperature() {
        static const std::string key = "MonteCarloTemperature";
        return key;
    }
    /**
     * Create a MonteCarloFlexibleBarostat.
     *
     * @param defaultPressure         the default pressure acting on the system (in bar)
     * @param defaultTemperature      the default temperature at which the system is being maintained (in Kelvin)
     * @param frequency               the frequency at which Monte Carlo pressure changes should be attempted (in time steps)
     * @param scaleMoleculesAsRigid   if true, coordinate scaling keeps molecules rigid, scaling only the center of mass
     *                                of each one.  If false, every atom is scaled independently.
     */
    MonteCarloFlexibleBarostat(double defaultPressure, double defaultTemperature, int frequency = 25, bool scaleMoleculesAsRigid = true);
    /**
     * Get the default pressure acting on the system (in bar).
     *
     * @return the default pressure acting on the system, measured in bar.
     */
    double getDefaultPressure() const {
        return defaultPressure;
    }
    /**
     * Set the default pressure acting on the system.  This will affect any new Contexts you create,
     * but not ones that already exist.
     *
     * @param pressure   the default pressure acting on the system, measured in bar.
     */
    void setDefaultPressure(double pressure);
    /**
     * Get the frequency (in time steps) at which Monte Carlo pressure changes should be attempted.  If this is set to
     * 0, the barostat is disabled.
     */
    int getFrequency() const {
        return frequency;
    }
    /**
     * Set the frequency (in time steps) at which Monte Carlo pressure changes should be attempted.  If this is set to
     * 0, the barostat is disabled.
     */
    void setFrequency(int freq);
    /**
     * Get the default temperature at which the system is being maintained, measured in Kelvin.
     */
    double getDefaultTemperature() const {
        return defaultTemperature;
    }
    /**
     * Set the default temperature at which the system is being maintained.  This will affect any new Contexts you create,
     * but not ones that already exist.
     *
     * @param temp     the system temperature, measured in Kelvin.
     */
    void setDefaultTemperature(double temp);
    /**
     * Get the random number seed.  See setRandomNumberSeed() for details.
     */
    int getRandomNumberSeed() const {
        return randomNumberSeed;
    }
    /**
     * Set the random number seed.  It is guaranteed that if two simulations are run
     * with different random number seeds, the sequence of Monte Carlo steps will be different.  On
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
    /**
     * Get whether scaling is applied to the centroid of each molecule while keeping
     * the molecules rigid, or to each atom independently.
     *
     * @returns true if scaling is applied to molecule centroids, false if it is applied to each atom independently.
     */
    bool getScaleMoleculesAsRigid() const {
        return scaleMoleculesAsRigid;
    }
    /**
     * Set whether scaling is applied to the centroid of each molecule while keeping
     * the molecules rigid, or to each atom independently.
     */
    void setScaleMoleculesAsRigid(bool rigid) {
        scaleMoleculesAsRigid = rigid;
    }
protected:
    ForceImpl* createImpl() const;
private:
    bool scaleMoleculesAsRigid;
    double defaultPressure, defaultTemperature;
    int frequency, randomNumberSeed;
};

} // namespace OpenMM

#endif /*OPENMM_MONTECARLOFLEXIBLEBAROSTAT_H_*/
