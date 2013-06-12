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
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
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
#include "Vec3.h"
#include <string>
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This class uses a Monte Carlo algorithm to adjust the size of the periodic box, simulating the
 * effect of constant pressure.
 *
 * Compared to MonteCarloBarostat, this class scales the three axes of the simulation cell independently.
 * The user supplies three doubles to specify the pressure along each axis.
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
     * Create a MonteCarloAnisotropicBarostat.
     *
     * @param defaultPressure   The default pressure acting on each axis (in bar)
     * @param temperature       the temperature at which the system is being maintained (in Kelvin)
     * @param frequency         the frequency at which Monte Carlo pressure changes should be attempted (in time steps)
     * @param scaleX            on/off switch for whether to scale the X axis
     * @param scaleY            on/off switch for whether to scale the Y axis
     * @param scaleZ            on/off switch for whether to scale the Z axis
     */
    MonteCarloAnisotropicBarostat(const Vec3& defaultPressure, double temperature, int frequency = 25, bool scaleX = 1, bool scaleY = 1, bool scaleZ = 1);
    /**
     * Get the default pressure (in bar).
     *
     * @return the default pressure acting on the system, measured in bar.
     */
    Vec3 getDefaultPressure() const {
        return defaultPressure;
    }
    /**
     * Get the true/false flag for scaling the X-axis.
     *
     * @return the true/false flag for scaling the X-axis.
     */
    bool getScaleX() const {
      return scaleX;
    }
    /**
     * Get the true/false flag for scaling the Y-axis.
     *
     * @return the true/false flag for scaling the Y-axis.
     */
    bool getScaleY() const {
      return scaleY;
    }
    /**
     * Get the true/false flag for scaling the Z-axis.
     *
     * @return the true/false flag for scaling the Z-axis.
     */
    bool getScaleZ() const {
      return scaleZ;
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
    int frequency, randomNumberSeed;

        double GetTemperature() const {
            return temperature;
        }

        void SetTemperature(double temperature) {
            this->temperature = temperature;
        }
    };

} // namespace OpenMM

#endif /*OPENMM_MONTECARLOANISOTROPICBAROSTAT_H_*/
