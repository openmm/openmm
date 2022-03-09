#ifndef OPENMM_DRUDEINTEGRATOR_H_
#define OPENMM_DRUDEINTEGRATOR_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2022 Stanford University and the Authors.      *
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

#include "openmm/Integrator.h"
#include "openmm/Kernel.h"
#include "openmm/internal/windowsExportDrude.h"

namespace OpenMM {

/**
 * A base class to encapsulate features common to Drude integrators.
 */

class OPENMM_EXPORT_DRUDE DrudeIntegrator : public Integrator {
public:
    /**
     * Create a DrudeSCFIntegrator.
     *
     * @param stepSize       the step size with which to integrator the system (in picoseconds)
     */
    DrudeIntegrator(double stepSize) {};
    /**
     * Advance a simulation through time by taking a series of time steps.
     *
     * @param steps   the number of time steps to take
     */
    virtual void step(int steps) override {};
    /**
     * Get the temperature of the heat bath applied to internal coordinates of Drude particles (in Kelvin).
     *
     * @return the temperature of the heat bath, measured in Kelvin
     */
    double getDrudeTemperature() const {
        return drudeTemperature;
    }
    /**
     * Set the temperature of the heat bath applied to internal coordinates of Drude particles (in Kelvin).
     *
     * @param temp    the temperature of the heat bath, measured in Kelvin
     */
    void setDrudeTemperature(double temp);
    /**
     * Get the maximum distance a Drude particle can ever move from its parent particle, measured in nm.  This is implemented
     * with a hard wall constraint.  The default value is 0.02.  If this distance is set to 0, the hard wall constraint is omitted.
     */
    double getMaxDrudeDistance() const;
    /**
     * Set the maximum distance a Drude particle can ever move from its parent particle, measured in nm.  This is implemented
     * with a hard wall constraint.  The default value is 0.02.  If this distance is set to 0, the hard wall constraint is omitted.
     */
    void setMaxDrudeDistance(double distance);
    /**
     * Set the random number seed.  The precise meaning of this parameter is undefined, and is left up
     * to each Platform to interpret in an appropriate way.  It is guaranteed that if two simulations
     * are run with different random number seeds, the sequence of random forces will be different.  On
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
     * Get the random number seed.  See setRandomNumberSeed() for details.
     */
    int getRandomNumberSeed() const {
        return randomNumberSeed;
    }
protected:
    /**
     * This will be called by the Context when it is created.  It informs the Integrator
     * of what context it will be integrating, and gives it a chance to do any necessary initialization.
     * It will also get called again if the application calls reinitialize() on the Context.
     */
    virtual void initialize(ContextImpl& context) override {};
    /**
     * This will be called by the Context when it is destroyed to let the Integrator do any necessary
     * cleanup.  It will also get called again if the application calls reinitialize() on the Context.
     */
    virtual void cleanup() override {};
    /**
     * Get the names of all Kernels used by this Integrator.
     */
    virtual std::vector<std::string> getKernelNames() override { return std::vector<std::string>(); }
    /**
     * Compute the kinetic energy of the system at the current time.
     */
    virtual double computeKineticEnergy() override { return 0; }
    /**
     * Return a list of velocities normally distributed around a target temperature, with the Drude
     * temperatures assigned according to the Drude temperature assigned to the integrator.
     *
     * @param system the system whose velocities are to be initialized.
     * @param temperature the target temperature in Kelvin.
     * @param randomSeed the random number seed to use when selecting velocities 
     */
    virtual std::vector<Vec3> getVelocitiesForTemperature(const System &system, double temperature,
                                                          int randomSeed) const override;

    int randomNumberSeed;
    double drudeTemperature, maxDrudeDistance;
};

} // namespace OpenMM

#endif /*OPENMM_DRUDEINTEGRATOR_H_*/
