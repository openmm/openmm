#ifndef OPENMM_DRUDENOSEHOOVERINTEGRATOR_H_
#define OPENMM_DRUDENOSEHOOVERINTEGRATOR_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2019-2022 Stanford University and the Authors.      *
 * Authors: Andreas Krämer and Andrew C. Simmonett                            *
 * Contributors: Peter Eastman                                                *
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

#include "openmm/NoseHooverIntegrator.h"
#include "openmm/Kernel.h"
#include "openmm/internal/windowsExportDrude.h"

namespace OpenMM {

/**
 * This Integrator simulates systems that include Drude particles.  It applies two different Nose-Hoover
 * chain thermostats to the different parts of the system.  The first is applied to ordinary particles (ones
 * that are not part of a Drude particle pair), as well as to the center of mass of each Drude particle pair.
 * A second thermostat, typically with a much lower temperature, is applied to the relative internal
 * displacement of each pair.
 *
 * This integrator can optionally set an upper limit on how far any Drude particle is ever allowed to
 * get from its parent particle.  This can sometimes help to improve stability.  The limit is enforced
 * with a hard wall constraint.  By default the limit is set to 0.02 nm.
 *
 * This Integrator requires the System to include a DrudeForce, which it uses to identify the Drude
 * particles.
 */

class OPENMM_EXPORT_DRUDE DrudeNoseHooverIntegrator : public NoseHooverIntegrator {
public:
    /**
     * Create a DrudeNoseHooverIntegrator.
     *
     * @param temperature the target temperature for the system (in Kelvin).
     * @param collisionFrequency the frequency of the system's interaction with the heat bath (in inverse picoseconds).
     * @param drudeTemperature the target temperature for the Drude particles, relative to their parent atom (in Kelvin).
     * @param drudeCollisionFrequency the frequency of the drude particles' interaction with the heat bath (in inverse picoseconds).
     * @param stepSize       the step size with which to integrator the system (in picoseconds)
     * @param chainLength the number of beads in the Nose-Hoover chain.
     * @param numMTS the number of step in the  multiple time step chain propagation algorithm.
     * @param numYoshidaSuzuki the number of terms in the Yoshida-Suzuki multi time step decomposition
     *        used in the chain propagation algorithm (must be 1, 3, or 5).
     */
    DrudeNoseHooverIntegrator(double temperature, double collisionFrequency, 
                              double drudeTemperature, double drudeCollisionFrequency, double stepSize, 
                              int chainLength = 3, int numMTS = 3, int numYoshidaSuzuki = 7);

    virtual ~DrudeNoseHooverIntegrator();
    /**
     * This will be called by the Context when it is created.  It informs the Integrator
     * of what context it will be integrating, and gives it a chance to do any necessary initialization.
     * It will also get called again if the application calls reinitialize() on the Context.
     */
    void initialize(ContextImpl& context) override;
    /**
     * Get the maximum distance a Drude particle can ever move from its parent particle, measured in nm.  This is implemented
     * with a hard wall constraint.  If this distance is set to 0 (the default), the hard wall constraint is omitted.
     */
    double getMaxDrudeDistance() const;
    /**
     * Set the maximum distance a Drude particle can ever move from its parent particle, measured in nm.  This is implemented
     * with a hard wall constraint.  If this distance is set to 0 (the default), the hard wall constraint is omitted.
     */
    void setMaxDrudeDistance(double distance);
    /**
     * Compute the kinetic energy of the drude particles at the current time.
     */
    double computeDrudeKineticEnergy();
    /**
     * Compute the kinetic energy of all (real and drude) particles at the current time.
     */
    double computeTotalKineticEnergy();
    /**
     * Compute the instantaneous temperature of the System, measured in Kelvin.
     * This is calculated based on the kinetic energy of the ordinary particles (ones
     * not attached to a Drude particle), as well as the center of mass motion of the
     * Drude particle pairs.  It does not include the internal motion of the pairs.
     * On average, this should be approximately equal to the value returned by
     * getTemperature().
     */
    double computeSystemTemperature();
    /**
     * Compute the instantaneous temperature of the Drude system, measured in Kelvin.
     * This is calculated based on the kinetic energy of the internal motion of Drude pairs
     * and should remain close to the prescribed Drude temperature.
     */ 
     double computeDrudeTemperature();
protected:
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
    double drudeTemperature;
};

} // namespace OpenMM

#endif /*OPENMM_DRUDENOSEHOOVERINTEGRATOR_H_*/
