
/* Portions copyright (c) 2008-2010 Stanford University and Simbios.
 * Contributors: Peter Eastman
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

#include <cmath>
#include <string.h>
#include <sstream>
#include <exception>
#include "SimTKOpenMMUtilities.h"
#include "ReferenceNoseHooverChain.h"

using std::vector;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   Constructor

   --------------------------------------------------------------------------------------- */

 ReferenceNoseHooverChain::ReferenceNoseHooverChain() {
 }

/**---------------------------------------------------------------------------------------

   Destructor

   --------------------------------------------------------------------------------------- */

 ReferenceNoseHooverChain::~ReferenceNoseHooverChain() {
 }

double ReferenceNoseHooverChain::propagate(double kineticEnergy, vector<double>& chainVelocities,
                                           vector<double>& chainPositions, int numDOFs,
                                           double temperature, double collisionFrequency, double timeStep,
                                           int numMTS, const vector<double>& YSWeights) const {
     double scale = 1;
     const double kT = BOLTZ * temperature;
     const size_t chainLength = chainPositions.size();
     std::vector<double> chainForces(chainLength, 0);
     std::vector<double> chainMasses(chainLength, kT/(collisionFrequency*collisionFrequency));
     chainMasses[0] *= numDOFs;
     double KE2 = 2 * kineticEnergy;
     chainForces[0] = (KE2 - numDOFs * kT) / chainMasses[0];
     for (int bead = 0; bead < chainLength - 1; ++bead) {
         chainForces[bead + 1] = (chainMasses[bead] * chainVelocities[bead] * chainVelocities[bead] - kT) / chainMasses[bead + 1];
     }
     for (int mts = 0; mts < numMTS; ++mts) {
         for (const auto &ys : YSWeights) {
             double wdt = ys * timeStep / numMTS;
             chainVelocities.back() += 0.25 * wdt * chainForces.back();
             for (int bead = chainLength - 2; bead >= 0; --bead) {
                 double aa = exp(-0.125 * wdt * chainVelocities[bead + 1]);
                 chainVelocities[bead] = aa * (chainVelocities[bead] * aa + 0.25 * wdt * chainForces[bead]);
             }
             // update particle velocities
             double aa = exp(-0.5 * wdt * chainVelocities[0]);
             scale *= aa;
             // update the thermostat positions
             for (int bead = 0; bead < chainLength; ++bead) {
                 chainPositions[bead] += 0.5 * chainVelocities[bead] * wdt;
             }
             // update the forces
             chainForces[0] = (scale * scale * KE2 - numDOFs * kT) / chainMasses[0];
             // update thermostat velocities
             for (int bead = 0; bead < chainLength - 1; ++bead) {
                 double aa = exp(-0.125 * wdt * chainVelocities[bead + 1]);
                 chainVelocities[bead] = aa * (aa * chainVelocities[bead] + 0.25 * wdt * chainForces[bead]);
                 chainForces[bead + 1] = (chainMasses[bead] * chainVelocities[bead] * chainVelocities[bead] - kT) / chainMasses[bead + 1];
             }
             chainVelocities[chainLength-1] += 0.25 * wdt * chainForces.back();
         }  // YS loop
     } // MTS loop
     return scale;
}
