/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2012-2017 Stanford University and the Authors.      *
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
#include <cmath>

using namespace OpenMM;
using namespace std;

void ReferenceVirtualSites::computePositions(const OpenMM::System& system, vector<OpenMM::Vec3>& atomCoordinates) {
    for (int i = 0; i < system.getNumParticles(); i++)
        if (system.isVirtualSite(i)) {
            if (dynamic_cast<const TwoParticleAverageSite*>(&system.getVirtualSite(i)) != NULL) {
                // A two particle average.
                
                const TwoParticleAverageSite& site = dynamic_cast<const TwoParticleAverageSite&>(system.getVirtualSite(i));
                int p1 = site.getParticle(0), p2 = site.getParticle(1);
                double w1 = site.getWeight(0), w2 = site.getWeight(1);
                atomCoordinates[i] = atomCoordinates[p1]*w1 + atomCoordinates[p2]*w2;
            }
            else if (dynamic_cast<const ThreeParticleAverageSite*>(&system.getVirtualSite(i)) != NULL) {
                // A three particle average.
                
                const ThreeParticleAverageSite& site = dynamic_cast<const ThreeParticleAverageSite&>(system.getVirtualSite(i));
                int p1 = site.getParticle(0), p2 = site.getParticle(1), p3 = site.getParticle(2);
                double w1 = site.getWeight(0), w2 = site.getWeight(1), w3 = site.getWeight(2);
                atomCoordinates[i] = atomCoordinates[p1]*w1 + atomCoordinates[p2]*w2 + atomCoordinates[p3]*w3;
            }
            else if (dynamic_cast<const OutOfPlaneSite*>(&system.getVirtualSite(i)) != NULL) {
                // An out of plane site.
                
                const OutOfPlaneSite& site = dynamic_cast<const OutOfPlaneSite&>(system.getVirtualSite(i));
                int p1 = site.getParticle(0), p2 = site.getParticle(1), p3 = site.getParticle(2);
                double w12 = site.getWeight12(), w13 = site.getWeight13(), wcross = site.getWeightCross();
                Vec3 v12 = atomCoordinates[p2]-atomCoordinates[p1];
                Vec3 v13 = atomCoordinates[p3]-atomCoordinates[p1];
                Vec3 cross = v12.cross(v13);
                atomCoordinates[i] = atomCoordinates[p1] + v12*w12 + v13*w13 + cross*wcross;
            }
            else if (dynamic_cast<const LocalCoordinatesSite*>(&system.getVirtualSite(i)) != NULL) {
                // A local coordinates site.
                
                const LocalCoordinatesSite& site = dynamic_cast<const LocalCoordinatesSite&>(system.getVirtualSite(i));
                int numParticles = site.getNumParticles();
                vector<double> originWeights, xWeights, yWeights;
                site.getOriginWeights(originWeights);
                site.getXWeights(xWeights);
                site.getYWeights(yWeights);
                Vec3 origin, xdir, ydir;
                for (int j = 0; j < numParticles; j++) {
                    Vec3 pos = atomCoordinates[site.getParticle(j)];
                    origin += pos*originWeights[j];
                    xdir += pos*xWeights[j];
                    ydir += pos*yWeights[j];
                }
                Vec3 localPosition = site.getLocalPosition();
                Vec3 zdir = xdir.cross(ydir);
                double normXdir = sqrt(xdir.dot(xdir));
                double normZdir = sqrt(zdir.dot(zdir));
                if (normXdir > 0.0)
                    xdir /= normXdir;
                if (normZdir > 0.0)
                    zdir /= normZdir;
                ydir = zdir.cross(xdir);
                atomCoordinates[i] = origin + xdir*localPosition[0] + ydir*localPosition[1] + zdir*localPosition[2];
            }
        }
}

void ReferenceVirtualSites::distributeForces(const OpenMM::System& system, const vector<OpenMM::Vec3>& atomCoordinates, vector<OpenMM::Vec3>& forces) {
    for (int i = 0; i < system.getNumParticles(); i++)
        if (system.isVirtualSite(i)) {
            Vec3 f = forces[i];
            if (dynamic_cast<const TwoParticleAverageSite*>(&system.getVirtualSite(i)) != NULL) {
                // A two particle average.
                
                const TwoParticleAverageSite& site = dynamic_cast<const TwoParticleAverageSite&>(system.getVirtualSite(i));
                int p1 = site.getParticle(0), p2 = site.getParticle(1);
                double w1 = site.getWeight(0), w2 = site.getWeight(1);
                forces[p1] += f*w1;
                forces[p2] += f*w2;
            }
            else if (dynamic_cast<const ThreeParticleAverageSite*>(&system.getVirtualSite(i)) != NULL) {
                // A three particle average.
                
                const ThreeParticleAverageSite& site = dynamic_cast<const ThreeParticleAverageSite&>(system.getVirtualSite(i));
                int p1 = site.getParticle(0), p2 = site.getParticle(1), p3 = site.getParticle(2);
                double w1 = site.getWeight(0), w2 = site.getWeight(1), w3 = site.getWeight(2);
                forces[p1] += f*w1;
                forces[p2] += f*w2;
                forces[p3] += f*w3;
            }
            else if (dynamic_cast<const OutOfPlaneSite*>(&system.getVirtualSite(i)) != NULL) {
                // An out of plane site.
                
                const OutOfPlaneSite& site = dynamic_cast<const OutOfPlaneSite&>(system.getVirtualSite(i));
                int p1 = site.getParticle(0), p2 = site.getParticle(1), p3 = site.getParticle(2);
                double w12 = site.getWeight12(), w13 = site.getWeight13(), wcross = site.getWeightCross();
                Vec3 v12 = atomCoordinates[p2]-atomCoordinates[p1];
                Vec3 v13 = atomCoordinates[p3]-atomCoordinates[p1];
                Vec3 f2(w12*f[0] - wcross*v13[2]*f[1] + wcross*v13[1]*f[2],
                        wcross*v13[2]*f[0] + w12*f[1] - wcross*v13[0]*f[2],
                       -wcross*v13[1]*f[0] + wcross*v13[0]*f[1] + w12*f[2]);
                Vec3 f3(w13*f[0] + wcross*v12[2]*f[1] - wcross*v12[1]*f[2],
                       -wcross*v12[2]*f[0] + w13*f[1] + wcross*v12[0]*f[2],
                        wcross*v12[1]*f[0] - wcross*v12[0]*f[1] + w13*f[2]);
                forces[p1] += f-f2-f3;
                forces[p2] += f2;
                forces[p3] += f3;
            }
            else if (dynamic_cast<const LocalCoordinatesSite*>(&system.getVirtualSite(i)) != NULL) {
                // A local coordinates site.
                
                const LocalCoordinatesSite& site = dynamic_cast<const LocalCoordinatesSite&>(system.getVirtualSite(i));
                int numParticles = site.getNumParticles();
                vector<double> originWeights, wx, wy;
                site.getOriginWeights(originWeights);
                site.getXWeights(wx);
                site.getYWeights(wy);
                Vec3 xdir, ydir;
                for (int j = 0; j < numParticles; j++) {
                    Vec3 pos = atomCoordinates[site.getParticle(j)];
                    xdir += pos*wx[j];
                    ydir += pos*wy[j];
                }
                Vec3 localPosition = site.getLocalPosition();
                Vec3 zdir = xdir.cross(ydir);
                double normXdir = sqrt(xdir.dot(xdir));
                double normZdir = sqrt(zdir.dot(zdir));
                double invNormXdir = (normXdir > 0.0 ? 1.0/normXdir : 0.0);
                double invNormZdir = (normZdir > 0.0 ? 1.0/normZdir : 0.0);
                Vec3 dx = xdir*invNormXdir;
                Vec3 dz = zdir*invNormZdir;
                Vec3 dy = dz.cross(dx);
                
                // The derivatives for this case are very complicated.  They were computed with SymPy then simplified by hand.
                
                vector<double> wxScaled(numParticles);
                for (int j = 0; j < numParticles; j++)
                    wxScaled[j] = wx[j]*invNormXdir;
                Vec3 fp1 = localPosition*f[0];
                Vec3 fp2 = localPosition*f[1];
                Vec3 fp3 = localPosition*f[2];
                for (int j = 0; j < numParticles; j++) {
                    double t1 = (wx[j]*ydir[0]-wy[j]*xdir[0])*invNormZdir;
                    double t2 = (wx[j]*ydir[1]-wy[j]*xdir[1])*invNormZdir;
                    double t3 = (wx[j]*ydir[2]-wy[j]*xdir[2])*invNormZdir;
                    double sx = t3*dz[1]-t2*dz[2];
                    double sy = t1*dz[2]-t3*dz[0];
                    double sz = t2*dz[0]-t1*dz[1];
                    int p = site.getParticle(j);
                    forces[p][0] += fp1[0]*wxScaled[j]*(1-dx[0]*dx[0]) + fp1[2]*(dz[0]*sx   ) + fp1[1]*((-dx[0]*dy[0]      )*wxScaled[j] + dy[0]*sx - dx[1]*t2 - dx[2]*t3) + f[0]*originWeights[j];
                    forces[p][1] += fp1[0]*wxScaled[j]*( -dx[0]*dx[1]) + fp1[2]*(dz[0]*sy+t3) + fp1[1]*((-dx[1]*dy[0]-dz[2])*wxScaled[j] + dy[0]*sy + dx[1]*t1);
                    forces[p][2] += fp1[0]*wxScaled[j]*( -dx[0]*dx[2]) + fp1[2]*(dz[0]*sz-t2) + fp1[1]*((-dx[2]*dy[0]+dz[1])*wxScaled[j] + dy[0]*sz + dx[2]*t1);
                    forces[p][0] += fp2[0]*wxScaled[j]*( -dx[1]*dx[0]) + fp2[2]*(dz[1]*sx-t3) - fp2[1]*(( dx[0]*dy[1]-dz[2])*wxScaled[j] - dy[1]*sx - dx[0]*t2);
                    forces[p][1] += fp2[0]*wxScaled[j]*(1-dx[1]*dx[1]) + fp2[2]*(dz[1]*sy   ) - fp2[1]*(( dx[1]*dy[1]      )*wxScaled[j] - dy[1]*sy + dx[0]*t1 + dx[2]*t3) + f[1]*originWeights[j];
                    forces[p][2] += fp2[0]*wxScaled[j]*( -dx[1]*dx[2]) + fp2[2]*(dz[1]*sz+t1) - fp2[1]*(( dx[2]*dy[1]+dz[0])*wxScaled[j] - dy[1]*sz - dx[2]*t2);
                    forces[p][0] += fp3[0]*wxScaled[j]*( -dx[2]*dx[0]) + fp3[2]*(dz[2]*sx+t2) + fp3[1]*((-dx[0]*dy[2]-dz[1])*wxScaled[j] + dy[2]*sx + dx[0]*t3);
                    forces[p][1] += fp3[0]*wxScaled[j]*( -dx[2]*dx[1]) + fp3[2]*(dz[2]*sy-t1) + fp3[1]*((-dx[1]*dy[2]+dz[0])*wxScaled[j] + dy[2]*sy + dx[1]*t3);
                    forces[p][2] += fp3[0]*wxScaled[j]*(1-dx[2]*dx[2]) + fp3[2]*(dz[2]*sz   ) + fp3[1]*((-dx[2]*dy[2]      )*wxScaled[j] + dy[2]*sz - dx[0]*t1 - dx[1]*t2) + f[2]*originWeights[j];
                }
           }
        }
}
