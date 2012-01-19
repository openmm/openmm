/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
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

#include "ReferenceVirtualSites.h"
#include "openmm/VirtualSite.h"

using namespace OpenMM;
using namespace std;

void ReferenceVirtualSites::computePositions(const OpenMM::System& system, vector<OpenMM::RealVec>& atomCoordinates) {
    for (int i = 0; i < system.getNumParticles(); i++)
        if (system.isVirtualSite(i)) {
            if (dynamic_cast<const VirtualSite::TwoParticleAverage*>(&system.getVirtualSite(i)) != NULL) {
                // A two particle average.
                
                const VirtualSite::TwoParticleAverage& site = dynamic_cast<const VirtualSite::TwoParticleAverage&>(system.getVirtualSite(i));
                int p1 = site.getParticle(0), p2 = site.getParticle(1);
                RealOpenMM w1 = site.getWeight(0), w2 = site.getWeight(1);
                atomCoordinates[i] = atomCoordinates[p1]*w1 + atomCoordinates[p2]*w2;
            }
            else if (dynamic_cast<const VirtualSite::ThreeParticleAverage*>(&system.getVirtualSite(i)) != NULL) {
                // A three particle average.
                
                const VirtualSite::ThreeParticleAverage& site = dynamic_cast<const VirtualSite::ThreeParticleAverage&>(system.getVirtualSite(i));
                int p1 = site.getParticle(0), p2 = site.getParticle(1), p3 = site.getParticle(2);
                RealOpenMM w1 = site.getWeight(0), w2 = site.getWeight(1), w3 = site.getWeight(2);
                atomCoordinates[i] = atomCoordinates[p1]*w1 + atomCoordinates[p2]*w2 + atomCoordinates[p3]*w3;
            }
            else if (dynamic_cast<const VirtualSite::OutOfPlane*>(&system.getVirtualSite(i)) != NULL) {
                // An out of plane site.
                
                const VirtualSite::OutOfPlane& site = dynamic_cast<const VirtualSite::OutOfPlane&>(system.getVirtualSite(i));
                int p1 = site.getParticle(0), p2 = site.getParticle(1), p3 = site.getParticle(2);
                RealOpenMM w12 = site.getWeight12(), w13 = site.getWeight13(), wcross = site.getWeightCross();
                RealVec v12 = atomCoordinates[p2]-atomCoordinates[p1];
                RealVec v13 = atomCoordinates[p3]-atomCoordinates[p1];
                RealVec cross = v12.cross(v13);
                atomCoordinates[i] = atomCoordinates[p1] + v12*w12 + v13*w13 + cross*wcross;
            }
        }

}

void ReferenceVirtualSites::distributeForces(const OpenMM::System& system, const vector<OpenMM::RealVec>& atomCoordinates, vector<OpenMM::RealVec>& forces) {
    for (int i = 0; i < system.getNumParticles(); i++)
        if (system.isVirtualSite(i)) {
            RealVec f = forces[i];
            if (dynamic_cast<const VirtualSite::TwoParticleAverage*>(&system.getVirtualSite(i)) != NULL) {
                // A two particle average.
                
                const VirtualSite::TwoParticleAverage& site = dynamic_cast<const VirtualSite::TwoParticleAverage&>(system.getVirtualSite(i));
                int p1 = site.getParticle(0), p2 = site.getParticle(1);
                RealOpenMM w1 = site.getWeight(0), w2 = site.getWeight(1);
                forces[p1] += f*w1;
                forces[p2] += f*w2;
            }
            else if (dynamic_cast<const VirtualSite::ThreeParticleAverage*>(&system.getVirtualSite(i)) != NULL) {
                // A three particle average.
                
                const VirtualSite::ThreeParticleAverage& site = dynamic_cast<const VirtualSite::ThreeParticleAverage&>(system.getVirtualSite(i));
                int p1 = site.getParticle(0), p2 = site.getParticle(1), p3 = site.getParticle(2);
                RealOpenMM w1 = site.getWeight(0), w2 = site.getWeight(1), w3 = site.getWeight(2);
                forces[p1] += f*w1;
                forces[p2] += f*w2;
                forces[p3] += f*w3;
            }
            else if (dynamic_cast<const VirtualSite::OutOfPlane*>(&system.getVirtualSite(i)) != NULL) {
                // An out of plane site.
                
                const VirtualSite::OutOfPlane& site = dynamic_cast<const VirtualSite::OutOfPlane&>(system.getVirtualSite(i));
                int p1 = site.getParticle(0), p2 = site.getParticle(1), p3 = site.getParticle(2);
                RealOpenMM w12 = site.getWeight12(), w13 = site.getWeight13(), wcross = site.getWeightCross();
                RealVec v12 = atomCoordinates[p2]-atomCoordinates[p1];
                RealVec v13 = atomCoordinates[p3]-atomCoordinates[p1];
                RealVec f2(w12*f[0] - wcross*v13[2]*f[1] + wcross*v13[1]*f[2],
                           wcross*v13[2]*f[0] + w12*f[1] - wcross*v13[0]*f[2],
                          -wcross*v13[1]*f[0] + wcross*v13[0]*f[1] + w12*f[2]);
                RealVec f3(w13*f[0] + wcross*v12[2]*f[1] - wcross*v12[1]*f[2],
                          -wcross*v12[2]*f[0] + w13*f[1] + wcross*v12[0]*f[2],
                           wcross*v12[1]*f[0] - wcross*v12[0]*f[1] + w13*f[2]);
                forces[p1] += f-f2-f3;
                forces[p2] += f2;
                forces[p3] += f3;
            }
        }
}
