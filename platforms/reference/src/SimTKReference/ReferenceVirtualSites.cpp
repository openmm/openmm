/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2012-2014 Stanford University and the Authors.      *
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
                int p1 = site.getParticle(0), p2 = site.getParticle(1), p3 = site.getParticle(2);
                Vec3 originWeights = site.getOriginWeights();
                Vec3 xWeights = site.getXWeights();
                Vec3 yWeights = site.getYWeights();
                Vec3 localPosition = site.getLocalPosition();
                Vec3 origin = atomCoordinates[p1]*originWeights[0] + atomCoordinates[p2]*originWeights[1] + atomCoordinates[p3]*originWeights[2];
                Vec3 xdir = atomCoordinates[p1]*xWeights[0] + atomCoordinates[p2]*xWeights[1] + atomCoordinates[p3]*xWeights[2];
                Vec3 ydir = atomCoordinates[p1]*yWeights[0] + atomCoordinates[p2]*yWeights[1] + atomCoordinates[p3]*yWeights[2];
                Vec3 zdir = xdir.cross(ydir);
                xdir /= sqrt(xdir.dot(xdir));
                zdir /= sqrt(zdir.dot(zdir));
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
                int p1 = site.getParticle(0), p2 = site.getParticle(1), p3 = site.getParticle(2);
                Vec3 originWeights = site.getOriginWeights();
                Vec3 wx = site.getXWeights();
                Vec3 wy = site.getYWeights();
                Vec3 localPosition = site.getLocalPosition();
                Vec3 xdir = atomCoordinates[p1]*wx[0] + atomCoordinates[p2]*wx[1] + atomCoordinates[p3]*wx[2];
                Vec3 ydir = atomCoordinates[p1]*wy[0] + atomCoordinates[p2]*wy[1] + atomCoordinates[p3]*wy[2];
                Vec3 zdir = xdir.cross(ydir);
                double invNormXdir = 1.0/sqrt(xdir.dot(xdir));
                double invNormZdir = 1.0/sqrt(zdir.dot(zdir));
                Vec3 dx = xdir*invNormXdir;
                Vec3 dz = zdir*invNormZdir;
                Vec3 dy = dz.cross(dx);
                
                // The derivatives for this case are very complicated.  They were computed with SymPy then simplified by hand.
                
                double t11 = (wx[0]*ydir[0]-wy[0]*xdir[0])*invNormZdir;
                double t12 = (wx[0]*ydir[1]-wy[0]*xdir[1])*invNormZdir;
                double t13 = (wx[0]*ydir[2]-wy[0]*xdir[2])*invNormZdir;
                double t21 = (wx[1]*ydir[0]-wy[1]*xdir[0])*invNormZdir;
                double t22 = (wx[1]*ydir[1]-wy[1]*xdir[1])*invNormZdir;
                double t23 = (wx[1]*ydir[2]-wy[1]*xdir[2])*invNormZdir;
                double t31 = (wx[2]*ydir[0]-wy[2]*xdir[0])*invNormZdir;
                double t32 = (wx[2]*ydir[1]-wy[2]*xdir[1])*invNormZdir;
                double t33 = (wx[2]*ydir[2]-wy[2]*xdir[2])*invNormZdir;
                double sx1 = t13*dz[1]-t12*dz[2];
                double sy1 = t11*dz[2]-t13*dz[0];
                double sz1 = t12*dz[0]-t11*dz[1];
                double sx2 = t23*dz[1]-t22*dz[2];
                double sy2 = t21*dz[2]-t23*dz[0];
                double sz2 = t22*dz[0]-t21*dz[1];
                double sx3 = t33*dz[1]-t32*dz[2];
                double sy3 = t31*dz[2]-t33*dz[0];
                double sz3 = t32*dz[0]-t31*dz[1];
                Vec3 wxScaled = wx*invNormXdir;
                Vec3 fp1 = localPosition*f[0];
                Vec3 fp2 = localPosition*f[1];
                Vec3 fp3 = localPosition*f[2];
                forces[p1][0] += fp1[0]*wxScaled[0]*(1-dx[0]*dx[0]) + fp1[2]*(dz[0]*sx1    ) + fp1[1]*((-dx[0]*dy[0]      )*wxScaled[0] + dy[0]*sx1 - dx[1]*t12 - dx[2]*t13) + f[0]*originWeights[0];
                forces[p1][1] += fp1[0]*wxScaled[0]*( -dx[0]*dx[1]) + fp1[2]*(dz[0]*sy1+t13) + fp1[1]*((-dx[1]*dy[0]-dz[2])*wxScaled[0] + dy[0]*sy1 + dx[1]*t11);
                forces[p1][2] += fp1[0]*wxScaled[0]*( -dx[0]*dx[2]) + fp1[2]*(dz[0]*sz1-t12) + fp1[1]*((-dx[2]*dy[0]+dz[1])*wxScaled[0] + dy[0]*sz1 + dx[2]*t11);
                forces[p2][0] += fp1[0]*wxScaled[1]*(1-dx[0]*dx[0]) + fp1[2]*(dz[0]*sx2    ) + fp1[1]*((-dx[0]*dy[0]      )*wxScaled[1] + dy[0]*sx2 - dx[1]*t22 - dx[2]*t23) + f[0]*originWeights[1];
                forces[p2][1] += fp1[0]*wxScaled[1]*( -dx[0]*dx[1]) + fp1[2]*(dz[0]*sy2+t23) + fp1[1]*((-dx[1]*dy[0]-dz[2])*wxScaled[1] + dy[0]*sy2 + dx[1]*t21);
                forces[p2][2] += fp1[0]*wxScaled[1]*( -dx[0]*dx[2]) + fp1[2]*(dz[0]*sz2-t22) + fp1[1]*((-dx[2]*dy[0]+dz[1])*wxScaled[1] + dy[0]*sz2 + dx[2]*t21);
                forces[p3][0] += fp1[0]*wxScaled[2]*(1-dx[0]*dx[0]) + fp1[2]*(dz[0]*sx3    ) + fp1[1]*((-dx[0]*dy[0]      )*wxScaled[2] + dy[0]*sx3 - dx[1]*t32 - dx[2]*t33) + f[0]*originWeights[2];
                forces[p3][1] += fp1[0]*wxScaled[2]*( -dx[0]*dx[1]) + fp1[2]*(dz[0]*sy3+t33) + fp1[1]*((-dx[1]*dy[0]-dz[2])*wxScaled[2] + dy[0]*sy3 + dx[1]*t31);
                forces[p3][2] += fp1[0]*wxScaled[2]*( -dx[0]*dx[2]) + fp1[2]*(dz[0]*sz3-t32) + fp1[1]*((-dx[2]*dy[0]+dz[1])*wxScaled[2] + dy[0]*sz3 + dx[2]*t31);
                forces[p1][0] += fp2[0]*wxScaled[0]*( -dx[1]*dx[0]) + fp2[2]*(dz[1]*sx1-t13) - fp2[1]*(( dx[0]*dy[1]-dz[2])*wxScaled[0] - dy[1]*sx1 - dx[0]*t12);
                forces[p1][1] += fp2[0]*wxScaled[0]*(1-dx[1]*dx[1]) + fp2[2]*(dz[1]*sy1    ) - fp2[1]*(( dx[1]*dy[1]      )*wxScaled[0] - dy[1]*sy1 + dx[0]*t11 + dx[2]*t13) + f[1]*originWeights[0];
                forces[p1][2] += fp2[0]*wxScaled[0]*( -dx[1]*dx[2]) + fp2[2]*(dz[1]*sz1+t11) - fp2[1]*(( dx[2]*dy[1]+dz[0])*wxScaled[0] - dy[1]*sz1 - dx[2]*t12);
                forces[p2][0] += fp2[0]*wxScaled[1]*( -dx[1]*dx[0]) + fp2[2]*(dz[1]*sx2-t23) - fp2[1]*(( dx[0]*dy[1]-dz[2])*wxScaled[1] - dy[1]*sx2 - dx[0]*t22);
                forces[p2][1] += fp2[0]*wxScaled[1]*(1-dx[1]*dx[1]) + fp2[2]*(dz[1]*sy2    ) - fp2[1]*(( dx[1]*dy[1]      )*wxScaled[1] - dy[1]*sy2 + dx[0]*t21 + dx[2]*t23) + f[1]*originWeights[1];
                forces[p2][2] += fp2[0]*wxScaled[1]*( -dx[1]*dx[2]) + fp2[2]*(dz[1]*sz2+t21) - fp2[1]*(( dx[2]*dy[1]+dz[0])*wxScaled[1] - dy[1]*sz2 - dx[2]*t22);
                forces[p3][0] += fp2[0]*wxScaled[2]*( -dx[1]*dx[0]) + fp2[2]*(dz[1]*sx3-t33) - fp2[1]*(( dx[0]*dy[1]-dz[2])*wxScaled[2] - dy[1]*sx3 - dx[0]*t32);
                forces[p3][1] += fp2[0]*wxScaled[2]*(1-dx[1]*dx[1]) + fp2[2]*(dz[1]*sy3    ) - fp2[1]*(( dx[1]*dy[1]      )*wxScaled[2] - dy[1]*sy3 + dx[0]*t31 + dx[2]*t33) + f[1]*originWeights[2];
                forces[p3][2] += fp2[0]*wxScaled[2]*( -dx[1]*dx[2]) + fp2[2]*(dz[1]*sz3+t31) - fp2[1]*(( dx[2]*dy[1]+dz[0])*wxScaled[2] - dy[1]*sz3 - dx[2]*t32);
                forces[p1][0] += fp3[0]*wxScaled[0]*( -dx[2]*dx[0]) + fp3[2]*(dz[2]*sx1+t12) + fp3[1]*((-dx[0]*dy[2]-dz[1])*wxScaled[0] + dy[2]*sx1 + dx[0]*t13);
                forces[p1][1] += fp3[0]*wxScaled[0]*( -dx[2]*dx[1]) + fp3[2]*(dz[2]*sy1-t11) + fp3[1]*((-dx[1]*dy[2]+dz[0])*wxScaled[0] + dy[2]*sy1 + dx[1]*t13);
                forces[p1][2] += fp3[0]*wxScaled[0]*(1-dx[2]*dx[2]) + fp3[2]*(dz[2]*sz1    ) + fp3[1]*((-dx[2]*dy[2]      )*wxScaled[0] + dy[2]*sz1 - dx[0]*t11 - dx[1]*t12) + f[2]*originWeights[0];
                forces[p2][0] += fp3[0]*wxScaled[1]*( -dx[2]*dx[0]) + fp3[2]*(dz[2]*sx2+t22) + fp3[1]*((-dx[0]*dy[2]-dz[1])*wxScaled[1] + dy[2]*sx2 + dx[0]*t23);
                forces[p2][1] += fp3[0]*wxScaled[1]*( -dx[2]*dx[1]) + fp3[2]*(dz[2]*sy2-t21) + fp3[1]*((-dx[1]*dy[2]+dz[0])*wxScaled[1] + dy[2]*sy2 + dx[1]*t23);
                forces[p2][2] += fp3[0]*wxScaled[1]*(1-dx[2]*dx[2]) + fp3[2]*(dz[2]*sz2    ) + fp3[1]*((-dx[2]*dy[2]      )*wxScaled[1] + dy[2]*sz2 - dx[0]*t21 - dx[1]*t22) + f[2]*originWeights[1];
                forces[p3][0] += fp3[0]*wxScaled[2]*( -dx[2]*dx[0]) + fp3[2]*(dz[2]*sx3+t32) + fp3[1]*((-dx[0]*dy[2]-dz[1])*wxScaled[2] + dy[2]*sx3 + dx[0]*t33);
                forces[p3][1] += fp3[0]*wxScaled[2]*( -dx[2]*dx[1]) + fp3[2]*(dz[2]*sy3-t31) + fp3[1]*((-dx[1]*dy[2]+dz[0])*wxScaled[2] + dy[2]*sy3 + dx[1]*t33);
                forces[p3][2] += fp3[0]*wxScaled[2]*(1-dx[2]*dx[2]) + fp3[2]*(dz[2]*sz3    ) + fp3[1]*((-dx[2]*dy[2]      )*wxScaled[2] + dy[2]*sz3 - dx[0]*t31 - dx[1]*t32) + f[2]*originWeights[2];
           }
        }
}
