#ifndef OPENMM_MONTECARLOANISOTROPICBAROSTAT_H_
#define OPENMM_MONTECARLOANISOTROPICBAROSTAT_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2013 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Lee-Ping Wang                                      *
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
#include "OpenMMException.h"
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This class uses a Monte Carlo algorithm to adjust the size of the periodic box, simulating the
 * effect of constant pressure.
 *
 * This class is similar to MonteCarloBarostat, but each Monte Carlo move is applied to only one axis
 * of the periodic box (unlike MonteCarloBarostat, which scales the entire box isotropically).  This
 * means that the box may change shape as well as size over the course of the simulation.  It also
 * allows you to specify a different pressure for each axis of the box, or to keep the box size fixed
 * along certain axes while still allowing it to change along others.
 * 
 * This class assumes the simulation is also being run at constant temperature, and requires you
 * to specify the system temperature (since it affects the acceptance probability for Monte Carlo
 * moves).  It does not actually perform temperature regulation, however.  You must use another
 * mechanism along with it to maintain the temperature, such as LangevinIntegrator or AndersenThermostat.
 */

class OPENMM_EXPORT MonteCarloAnisotropicBarostat : public Force {
public:
    /**
     * This is the name of the parameter which stores the current pressure acting on
     * the X-axis (in bar).
     */
    static const std::string& PressureX() {
        static const std::string key = "MonteCarloPressureX";
        return key;
    }
    /**
     * This is the name of the parameter which stores the current pressure acting on
     * the Y-axis (in bar).
     */
    static const std::string& PressureY() {
        static const std::string key = "MonteCarloPressureY";
        return key;
    }
    /**
     * This is the name of the parameter which stores the current pressure acting on
     * the Z-axis (in bar).
     */
    static const std::string& PressureZ() {
        static const std::string key = "MonteCarloPressureZ";
        return key;
    }
    /**
     * Create a MonteCarloAnisotropicBarostat. You can only specify one of coupleXY, coupleXZ, or coupleYZ. And
     * each of the dimensions coupled must be able to be scaled.
     *
     * @param defaultPressure   The default pressure acting on each axis (in bar)
     * @param temperature       the temperature at which the system is being maintained (in Kelvin)
     * @param scaleX            whether to allow the X dimension of the periodic box to change size
     * @param scaleY            whether to allow the Y dimension of the periodic box to change size
     * @param scaleZ            whether to allow the Z dimension of the periodic box to change size
     * @param frequency         the frequency at which Monte Carlo pressure changes should be attempted (in time steps)
     * @param coupleXY          whether to change the X and Y dimensions of the periodic box together (semi-isotropic)
     * @param coupleXZ          whether to change the X and Z dimensions of the periodic box together (semi-isotropic)
     * @param coupleYZ          whether to change the Y and Z dimensions of the periodic box together (semi-isotropic)
     */
    MonteCarloAnisotropicBarostat(const Vec3& defaultPressure, double temperature, bool scaleX = true, bool scaleY = true, bool scaleZ = true,
                                  int frequency = 25, bool coupleXY = false, bool coupleXZ = false, bool coupleYZ = false);
    /**
     * Get the default pressure (in bar).
     *
     * @return the default pressure acting along each axis, measured in bar.
     */
    const Vec3& getDefaultPressure() const {
        return defaultPressure;
    }
    /**
     * Set the default pressure acting on the system.  This will affect any new Contexts you create,
     * but not ones that already exist.
     *
     * @param pressure   the default pressure acting on the system, measured in bar.
     */
    void setDefaultPressure(const Vec3& pressure) {
        defaultPressure = pressure;
    }
    /**
     * Get whether to allow the X dimension of the periodic box to change size.
     */
    bool getScaleX() const {
      return scaleX;
    }
    /**
     * Get whether to allow the Y dimension of the periodic box to change size.
     */
    bool getScaleY() const {
      return scaleY;
    }
    /**
     * Get whether to allow the Z dimension of the periodic box to change size.
     */
    bool getScaleZ() const {
      return scaleZ;
    }
    /**
     * Get whether to change the X- and Y-dimensions of the periodic box the same amount
     */
    bool getCoupleXY() const {
        return coupleXY;
    }
    /**
     * Sets whether to change the X- and Y-dimensions of the periodic box the same amount
     * throws OpenMMException if two other axes are coupled already
     */
    void setCoupleXY(bool set) {
        if (set && (coupleXZ || coupleYZ))
            throw OpenMMException("MonteCarloAnisotropicBarostat: Can only couple one pair of axes.");
        coupleXY = set;
    }
    /**
     * Get whether to change the X- and Z-dimensions of the periodic box the same amount
     */
    bool getCoupleXZ() const {
        return coupleXZ;
    }
    /**
     * Sets whether to change the X- and Z-dimensions of the periodic box the same amount
     * throws OpenMMException if two other axes are coupled already
     */
    void setCoupleXZ(bool set) {
        if (set && (coupleXY || coupleYZ))
            throw OpenMMException("MonteCarloAnisotropicBarostat: Can only couple one pair of axes.");
        coupleXZ = set;
    }
    /**
     * Get whether to change the Y- and Z-dimensions of the periodic box the same amount
     */
    bool getCoupleYZ() const {
        return coupleYZ;
    }
    /**
     * Sets whether to change the Y- and Z-dimensions of the periodic box the same amount
     * throws OpenMMException if two other axes are coupled already
     */
    void setCoupleYZ(bool set) {
        if (set && (coupleXY || coupleXZ))
            throw OpenMMException("MonteCarloAnisotropicBarostat: Can only couple one pair of axes.");
        coupleYZ = set;
    }
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
    void setFrequency(int freq) {
        frequency = freq;
    }
    /**
     * Get the temperature at which the system is being maintained, measured in Kelvin.
     */
    double getTemperature() const {
        return temperature;
    }
    /**
     * Set the temperature at which the system is being maintained.
     *
     * @param temp     the system temperature, measured in Kelvin.
     */
    void setTemperature(double temp) {
        temperature = temp;
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
     */
    void setRandomNumberSeed(int seed) {
        randomNumberSeed = seed;
    }
protected:
    ForceImpl* createImpl() const;
private:
    Vec3 defaultPressure;
    double temperature;
    bool scaleX, scaleY, scaleZ;
    bool coupleXY, coupleXZ, coupleYZ;
    int frequency, randomNumberSeed;
};

} // namespace OpenMM

#endif /*OPENMM_MONTECARLOANISOTROPICBAROSTAT_H_*/
