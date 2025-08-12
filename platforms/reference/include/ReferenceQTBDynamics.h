/* Portions copyright (c) 2025 Stanford University and the Authors.
 * Authors: Peter Eastman
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef __ReferenceQTBDynamics_H__
#define __ReferenceQTBDynamics_H__

#include "ReferenceDynamics.h"
#include "openmm/QTBIntegrator.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/ThreadPool.h"
#include "openmm/internal/windowsExport.h"
#include <complex>

namespace OpenMM {

class OPENMM_EXPORT ReferenceQTBDynamics : public ReferenceDynamics {
protected:
    double friction, lastTemperature;
    int segmentLength, stepIndex, numFreq;
    std::vector<OpenMM::Vec3> xPrime, oldx, randomForce, segmentVelocity;
    std::vector<double> inverseMasses, typeAdaptationRate, typeMass;
    std::vector<double> noise, theta, thetad, cutoffFunction;
    std::vector<int> particleType;
    std::vector<std::vector<int> > typeParticles;
    std::vector<std::vector<double> > adaptedFriction;

public:
    /**
     * Constructor
     * 
     * @param system       the system to simulate
     * @param integrator   the integrator being used
     */
    ReferenceQTBDynamics(const System& system, const QTBIntegrator& integrator);

    /**
     * Destructor
     */
    ~ReferenceQTBDynamics();

    /**
     * Perform a time step, updating the positions and velocities.
     * 
     * @param context             the context this integrator is updating
     * @param system              the System to be integrated
     * @param atomCoordinates     atom coordinates
     * @param velocities          velocities
     * @param masses              atom masses
     * @param tolerance           the constraint tolerance
     * @param boxVectors          the current periodic box vectors
     * @param threads      a ThreadPool to use for parallelization
     */
    void update(OpenMM::ContextImpl& context, std::vector<OpenMM::Vec3>& atomCoordinates,
                std::vector<OpenMM::Vec3>& velocities, std::vector<double>& masses, double tolerance, const Vec3* boxVectors, ThreadPool& threads);

    /**
     * The first stage of the update algorithm.
     */
    virtual void updatePart1(int numParticles, std::vector<OpenMM::Vec3>& velocities, std::vector<OpenMM::Vec3>& forces);

    /**
     * The second stage of the update algorithm.
     */
    virtual void updatePart2(int numParticles, std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<OpenMM::Vec3>& velocities,
                             std::vector<OpenMM::Vec3>& xPrime);

    /**
     * The third stage of the update algorithm.
     */
    virtual void updatePart3(OpenMM::ContextImpl& context, int numParticles, std::vector<OpenMM::Vec3>& atomCoordinates,
                             std::vector<OpenMM::Vec3>& velocities, std::vector<OpenMM::Vec3>& xPrime);
    /**
     * Get the adapted friction coefficients for a particle.
     * 
     * @param particle   the index of the particle for which to get the friction
     * @param friction   the adapted friction coefficients used in generating the
     *                   random force.
     */
    void getAdaptedFriction(int particle, std::vector<double>& friction) const;
    /**
     * Set the adapted friction coefficients for a particle.  This affects the
     * specified particle, and all others that have the same type.
     * 
     * @param particle   the index of the particle for which to get the friction
     * @param friction   the adapted friction coefficients used in generating the
     *                   random force.
     */
    void setAdaptedFriction(int particle, const std::vector<double>& friction);
    /**
     * Write the adapted friction to a checkpoint.
     */
    void createCheckpoint(std::ostream& stream) const;
    /**
     * Load the adapted friction from a checkpoint.
     */
    void loadCheckpoint(std::istream& stream);

private:
    /**
     * Generate noise for the next segment.
     */
    void generateNoise(int numParticles, std::vector<double>& masses, ThreadPool& threads);
    /**
     * Update the friction rates used for generating noise.
     */
    void adaptFriction(ThreadPool& threads);
};

} // namespace OpenMM

#endif // __ReferenceQTBDynamics_H__
