/* Portions copyright (c) 2025 Stanford University and Simbios.
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

#include "ReferenceRGForce.h"
#include <cmath>

using namespace OpenMM;
using namespace std;

ReferenceRGForce::ReferenceRGForce(vector<int>& particles) : particles(particles) {
}

double ReferenceRGForce::calculateIxn(vector<Vec3>& atomCoordinates, vector<Vec3>& forces) const {
    // Compute the center position.

    int numParticles = particles.size();
    Vec3 center;
    for (int i : particles)
        center += atomCoordinates[i];
    center /= numParticles;

    // Compute the radius of gyration.

    double sum = 0.0;
    for (int i = 0; i < numParticles; i++) {
        Vec3 delta = atomCoordinates[particles[i]]-center;
        sum += delta.dot(delta);
    }
    double rg = sqrt(sum/numParticles);

    // Compute the forces.

    double scale = 1.0/(rg*numParticles);
    for (int i = 0; i < numParticles; i++) {
        Vec3 delta = atomCoordinates[particles[i]]-center;
        forces[particles[i]] -= scale*delta;
    }
    return rg;
}
