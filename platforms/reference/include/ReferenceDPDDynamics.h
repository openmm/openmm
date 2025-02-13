
/* Portions copyright (c) 2006-2025 Stanford University and Simbios.
 * Contributors: Pande Group
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

#ifndef __ReferenceDPDDynamics_H__
#define __ReferenceDPDDynamics_H__

#include "ReferenceDynamics.h"
#include "ReferenceNeighborList.h"
#include "openmm/DPDIntegrator.h"
#include "openmm/Vec3.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/windowsExport.h"
#include <vector>

namespace OpenMM {

class OPENMM_EXPORT ReferenceDPDDynamics : public ReferenceDynamics {
protected:
      std::vector<Vec3> xPrime, oldx;
      std::vector<double> masses, inverseMasses;
      bool periodic;
      Vec3 periodicBoxVectors[3];
      NeighborList neighborList;
      std::vector<int> particleType;
      std::vector<std::vector<double> > frictionTable, cutoffTable;
      double maxCutoff;
public:
    /**
     * Constructor
     * 
     * @param system       the system to simulate
     * @param integrator   the integrator being used
     */
    ReferenceDPDDynamics(const System& system, const DPDIntegrator& integrator);

    /**
     * Destructor
     */
    ~ReferenceDPDDynamics();

    /**
     * Set the periodic box vectors.
     */
    void setPeriodicBoxVectors(OpenMM::Vec3* vectors);
    /**
     * Get the maximum cutoff distance for any pair of types.
     */
    double getMaxCutoff() const;
    /**
     * Perform a time step, updating the positions and velocities.
     * 
     * @param context             the context this integrator is updating
     * @param system              the System to be integrated
     * @param atomCoordinates     atom coordinates
     * @param velocities          velocities
     * @param masses              atom masses
     * @param tolerance           the constraint tolerance
     */
    void update(OpenMM::ContextImpl& context, std::vector<OpenMM::Vec3>& atomCoordinates,
                std::vector<OpenMM::Vec3>& velocities, std::vector<double>& masses, double tolerance);

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
};

} // namespace OpenMM

#endif // __ReferenceDPDDynamics_H__
