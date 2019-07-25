#ifndef OPENMM_MONTECARLOMEMBRANEBAROSTAT_H_
#define OPENMM_MONTECARLOMEMBRANEBAROSTAT_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2016 Stanford University and the Authors.      *
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

namespace OpenMM {

/**
 * This is a Monte Carlo barostat designed specifically for membrane simulations.  It assumes the
 * membrane lies in the XY plane.  The Monte Carlo acceptance criterion includes a term to model
 * isotropic pressure, which depends on the volume of the periodic box, and a second term to model
 * surface tension, which depends on the cross sectional area of the box in the XY plane.  Note
 * that pressure and surface tension are defined with opposite senses: a larger pressure tends to
 * make the box smaller, but a larger surface tension tends to make the box larger.
 * 
 * There are options for configuring exactly how the various box dimensions are allowed to change:
 * 
 * <ul>
 * <li>The X and Y axes may be treated isotropically, in which case they always scale by the same
 * amount and remain in proportion to each other; or they may be treated anisotropically, in which
 * case they can vary independently of each other.</li>
 * <li>The Z axis can be allowed to vary independently of the other axes; or held fixed; or constrained
 * to vary in inverse proportion to the other two axes, so that the total box volume remains fixed.</li>
 * </ul>
 *
 * This class assumes the simulation is also being run at constant temperature, and requires you
 * to specify the system temperature (since it affects the acceptance probability for Monte Carlo
 * moves).  It does not actually perform temperature regulation, however.  You must use another
 * mechanism along with it to maintain the temperature, such as LangevinIntegrator or AndersenThermostat.
 */

class OPENMM_EXPORT MonteCarloMembraneBarostat : public Force {
public:
    /**
     * This is an enumeration of the different behaviors for the X and Y axes.
     */
    enum XYMode {
        /**
         * The X and Y axes are always scaled by the same amount, so the ratio of their lengths remains constant.
         */
        XYIsotropic = 0,
        /**
         * The X and Y axes are allowed to vary independently of each other.
         */
        XYAnisotropic = 1
    };
    /**
     * This is an enumeration of the different behaviors for Z axis.
     */
    enum ZMode {
        /**
         * The Z axis is allowed to vary freely, independent of the other two axes.
         */
        ZFree = 0,
        /**
         * The Z axis is held fixed and does not change.
         */
        ZFixed = 1,
        /**
         * The Z axis is always scaled in inverse proportion to the other two axes so the box volume remains
         * fixed.  Note that in this mode pressure has no effect on the system, only surface tension.
         */
        ConstantVolume = 2
    };
    /**
     * This is the name of the parameter which stores the current pressure acting on
     * the system (in bar).
     */
    static const std::string& Pressure() {
        static const std::string key = "MembraneMonteCarloPressure";
        return key;
    }
    /**
     * This is the name of the parameter which stores the current surface tension acting on
     * the system (in bar*nm).
     */
    static const std::string& SurfaceTension() {
        static const std::string key = "MembraneMonteCarloSurfaceTension";
        return key;
    }
    /**
     * This is the name of the parameter which stores the current temperature at which the
     * system is being maintained (in Kelvin)
     */
    static const std::string& Temperature() {
        static const std::string key = "MembraneMonteCarloTemperature";
        return key;
    }
    /**
     * Create a MonteCarloMembraneBarostat.
     *
     * @param defaultPressure        the default pressure acting on the system (in bar)
     * @param defaultSurfaceTension  the default surface tension acting on the system (in bar*nm)
     * @param defaultTemperature     the default temperature at which the system is being maintained (in Kelvin)
     * @param xymode                 the mode specifying the behavior of the X and Y axes
     * @param zmode                  the mode specifying the behavior of the Z axis
     * @param frequency              the frequency at which Monte Carlo volume changes should be attempted (in time steps)
     */
    MonteCarloMembraneBarostat(double defaultPressure, double defaultSurfaceTension, double defaultTemperature, XYMode xymode, ZMode zmode, int frequency = 25);
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
    void setDefaultPressure(double pressure) {
        defaultPressure = pressure;
    }
    /**
     * Get the default surface tension acting on the system (in bar*nm).
     *
     * @return the default surface tension acting on the system, measured in bar*nm.
     */
    double getDefaultSurfaceTension() const {
        return defaultSurfaceTension;
    }
    /**
     * Set the default surface tension acting on the system.  This will affect any new Contexts you create,
     * but not ones that already exist.
     *
     * @param surfaceTension   the default surface tension acting on the system, measured in bar.
     */
    void setDefaultSurfaceTension(double surfaceTension) {
        defaultSurfaceTension = surfaceTension;
    }
    /**
     * Get the frequency (in time steps) at which Monte Carlo volume changes should be attempted.  If this is set to
     * 0, the barostat is disabled.
     */
    int getFrequency() const {
        return frequency;
    }
    /**
     * Set the frequency (in time steps) at which Monte Carlo volume changes should be attempted.  If this is set to
     * 0, the barostat is disabled.
     */
    void setFrequency(int freq) {
        frequency = freq;
    }
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
    void setDefaultTemperature(double temp) {
        defaultTemperature = temp;
    }
    /**
     * Get the mode specifying the behavior of the X and Y axes.
     */
    XYMode getXYMode() const {
        return xymode;
    }
    /**
     * Set the mode specifying the behavior of the X and Y axes.
     */
    void setXYMode(XYMode mode) {
        xymode = mode;
    }
    /**
     * Get the mode specifying the behavior of the Z axis.
     */
    ZMode getZMode() const {
        return zmode;
    }
    /**
     * Set the mode specifying the behavior of the Z axis.
     */
    void setZMode(ZMode mode) {
        zmode = mode;
    }
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
protected:
    ForceImpl* createImpl() const;
private:
    double defaultPressure, defaultSurfaceTension, defaultTemperature;
    XYMode xymode;
    ZMode zmode;
    int frequency, randomNumberSeed;
};

} // namespace OpenMM

#endif /*OPENMM_MONTECARLOMEMBRANEBAROSTAT_H_*/
