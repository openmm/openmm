/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2016 Stanford University and the Authors.           *
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

#ifdef _MSC_VER
    // Prevent Windows from defining macros that interfere with other code.
    #define NOMINMAX
#endif
#include "ReferenceGayBerneForce.h"
#include "ReferenceForce.h"
#include "openmm/OpenMMException.h"
#include <algorithm>
#include <cmath>

using namespace OpenMM;
using namespace std;

ReferenceGayBerneForce::ReferenceGayBerneForce(const GayBerneForce& force) {
    // Record the force parameters.

    int numParticles = force.getNumParticles();
    particles.resize(numParticles);
    for (int i = 0; i < numParticles; i++) {
        ParticleInfo& p = particles[i];
        double sx, sy, sz;
        force.getParticleParameters(i, p.sigma, p.epsilon, p.xparticle, p.yparticle, sx, sy, sz, p.ex, p.ey, p.ez);
        p.rx = 0.5*sx;
        p.ry = 0.5*sy;
        p.rz = 0.5*sz;
    }
    int numExceptions = force.getNumExceptions();
    exceptions.resize(numExceptions);
    for (int i = 0; i < numExceptions; i++) {
        ExceptionInfo& e = exceptions[i];
        force.getExceptionParameters(i, e.particle1, e.particle2, e.sigma, e.epsilon);
        exclusions.insert(make_pair(min(e.particle1, e.particle2), max(e.particle1, e.particle2)));
    }
    nonbondedMethod = force.getNonbondedMethod();
    cutoffDistance = force.getCutoffDistance();
    switchingDistance = force.getSwitchingDistance();
    useSwitchingFunction = force.getUseSwitchingFunction();

    // Allocate workspace for calculations.

    s.resize(numParticles);
    A.resize(numParticles);
    B.resize(numParticles);
    G.resize(numParticles);

    // We can precompute the shape factors.

    for (int i = 0; i < numParticles; i++) {
        ParticleInfo& p = particles[i];
        s[i] = (p.rx*p.ry + p.rz*p.rz)*sqrt(p.rx*p.ry);
    }
}

double ReferenceGayBerneForce::calculateForce(const vector<Vec3>& positions, vector<Vec3>& forces, const Vec3* boxVectors) {
    if (nonbondedMethod == GayBerneForce::CutoffPeriodic) {
        double minAllowedSize = 1.999999*cutoffDistance;
        if (boxVectors[0][0] < minAllowedSize || boxVectors[1][1] < minAllowedSize || boxVectors[2][2] < minAllowedSize)
            throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
    }

    // First find the orientations of the particles and compute the matrices we'll be needing.

    computeEllipsoidFrames(positions);

    // Compute standard interactions.

    double energy = 0;
    int numParticles = particles.size();
    vector<Vec3> torques(numParticles, Vec3());
    for (int i = 1; i < numParticles; i++) {
        if (particles[i].epsilon == 0.0)
            continue;
        for (int j = 0; j < i; j++) {
            if (particles[j].epsilon == 0.0)
                continue;
            if (exclusions.find(make_pair(j, i)) != exclusions.end())
                continue; // This interaction will be handled by an exception.
            double sigma = 0.5*(particles[i].sigma+particles[j].sigma);
            double epsilon = sqrt(particles[i].epsilon*particles[j].epsilon);
            energy += computeOneInteraction(i, j, sigma, epsilon, positions, forces, torques, boxVectors);
        }
    }

    // Compute exceptions.

    int numExceptions = exceptions.size();
    for (int i = 0; i < numExceptions; i++) {
        ExceptionInfo& e = exceptions[i];
        energy += computeOneInteraction(e.particle1, e.particle2, e.sigma, e.epsilon, positions, forces, torques, boxVectors);
    }
    
    // Apply torques.
    
    applyTorques(positions, forces, torques);
    return energy;
}

void ReferenceGayBerneForce::computeEllipsoidFrames(const vector<Vec3>& positions) {
    int numParticles = particles.size();
    for (int particle = 0; particle < numParticles; particle++) {
        ParticleInfo& p = particles[particle];

        // Compute the local coordinate system of the ellipsoid;

        Vec3 xdir, ydir, zdir;
        if (p.xparticle == -1) {
            xdir = Vec3(1, 0, 0);
            ydir = Vec3(0, 1, 0);
        }
        else {
            xdir = positions[particle]-positions[p.xparticle];
            xdir /= sqrt(xdir.dot(xdir));
            if (p.yparticle == -1) {
                if (xdir[1] > -0.5 && xdir[1] < 0.5)
                    ydir = Vec3(0, 1, 0);
                else
                    ydir = Vec3(1, 0, 0);
            }
            else
                ydir = positions[particle]-positions[p.yparticle];
            ydir -= xdir*(xdir.dot(ydir));
            ydir /= sqrt(ydir.dot(ydir));
        }
        zdir = xdir.cross(ydir);

        // Compute matrices we will need later.

        double (&a)[3][3] = A[particle].v;
        double (&b)[3][3] = B[particle].v;
        double (&g)[3][3] = G[particle].v;
        a[0][0] = xdir[0];
        a[0][1] = xdir[1];
        a[0][2] = xdir[2];
        a[1][0] = ydir[0];
        a[1][1] = ydir[1];
        a[1][2] = ydir[2];
        a[2][0] = zdir[0];
        a[2][1] = zdir[1];
        a[2][2] = zdir[2];
        Vec3 r2(p.rx*p.rx, p.ry*p.ry, p.rz*p.rz);
        Vec3 e2(1/sqrt(p.ex), 1/sqrt(p.ey), 1/sqrt(p.ez));
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
                b[i][j] = 0;
                g[i][j] = 0;
                for (int k = 0; k < 3; k++) {
                    b[i][j] += a[k][i]*e2[k]*a[k][j];
                    g[i][j] += a[k][i]*r2[k]*a[k][j];
                }
            }
    }
}

void ReferenceGayBerneForce::applyTorques(const vector<Vec3>& positions, vector<Vec3>& forces, const vector<Vec3>& torques) {
    int numParticles = particles.size();
    for (int particle = 0; particle < numParticles; particle++) {
        ParticleInfo& p = particles[particle];
        Vec3 pos = positions[particle];
        if (p.xparticle != -1) {
            // Apply a force to the x particle.
            
            Vec3 dx = positions[p.xparticle]-pos;
            double dx2 = dx.dot(dx);
            Vec3 f = torques[particle].cross(dx)/dx2;
            forces[p.xparticle] += f;
            forces[particle] -= f;
            if (p.yparticle != -1) {
                // Apply a force to the y particle.  This is based on the component of the torque
                // that was not already applied to the x particle.
                
                Vec3 dy = positions[p.yparticle]-pos;
                double dy2 = dy.dot(dy);
                Vec3 torque = dx*(torques[particle].dot(dx)/dx2);
                f = torque.cross(dy)/dy2;
                forces[p.yparticle] += f;
                forces[particle] -= f;
            }
        }
    }
}

double ReferenceGayBerneForce::computeOneInteraction(int particle1, int particle2, double sigma, double epsilon, const vector<Vec3>& positions,
        vector<Vec3>& forces, vector<Vec3>& torques, const Vec3* boxVectors) {
    // Compute the displacement and check against the cutoff.

    double deltaR[ReferenceForce::LastDeltaRIndex];
    if (nonbondedMethod == GayBerneForce::CutoffPeriodic)
        ReferenceForce::getDeltaRPeriodic(positions[particle2], positions[particle1], boxVectors, deltaR);
    else
        ReferenceForce::getDeltaR(positions[particle2], positions[particle1], deltaR);
    double r = deltaR[ReferenceForce::RIndex];
    if (nonbondedMethod != GayBerneForce::NoCutoff && r >= cutoffDistance)
        return 0;

    // Compute vectors and matrices we'll be needing.

    double rInv = 1/r;
    Vec3 dr(deltaR[ReferenceForce::XIndex], deltaR[ReferenceForce::YIndex], deltaR[ReferenceForce::ZIndex]);
    Vec3 drUnit = dr*rInv;
    Matrix B12 = B[particle1]+B[particle2];
    Matrix G12 = G[particle1]+G[particle2];
    Matrix B12inv = B12.inverse();
    Matrix G12inv = G12.inverse();
    double detG12 = G12.determinant();
    
    // Compute the switching function.

    double switchValue = 1, switchDeriv = 0;
    if (useSwitchingFunction && r > switchingDistance) {
        double t = (r-switchingDistance)/(cutoffDistance-switchingDistance);
        switchValue = 1+t*t*t*(-10+t*(15-t*6));
        switchDeriv = t*t*(-30+t*(60-t*30))/(cutoffDistance-switchingDistance);
    }

    // Estimate the distance between the ellipsoids and compute the first terms needed for the energy.

    double sigma12 = 1/sqrt(0.5*drUnit.dot(G12inv*drUnit));
    double h12 = r - sigma12;
    double rho = sigma/(h12+sigma);
    double rho2 = rho*rho;
    double rho6 = rho2*rho2*rho2;
    double u = 4*epsilon*(rho6*rho6-rho6);
    double eta = sqrt(2*s[particle1]*s[particle2]/detG12);
    double chi = 2*drUnit.dot(B12inv*drUnit);
    chi *= chi;
    double energy = u*eta*chi;
    
    // Compute the terms needed for the force.

    Vec3 kappa = G12inv*dr;
    Vec3 iota = B12inv*dr;
    double rInv2 = rInv*rInv;
    double dUSLJdr = 24*epsilon*(2*rho6-1)*rho6*rho/sigma;
    double temp = 0.5*sigma12*sigma12*sigma12*rInv2;
    Vec3 dudr = (drUnit + (kappa-drUnit*kappa.dot(drUnit))*temp)*dUSLJdr;
    Vec3 dchidr = (iota-drUnit*iota.dot(drUnit))*(-8*rInv2*sqrt(chi));
    Vec3 force = (dchidr*u + dudr*chi)*(eta*switchValue) - drUnit*(energy*switchDeriv);
    forces[particle1] += force;
    forces[particle2] -= force;

    // Compute the terms needed for the torque.

    for (int j = 0; j < 2; j++) {
        int particle = (j == 0 ? particle1 : particle2);
        Vec3 dudq = (kappa*G[particle]).cross(kappa*(temp*dUSLJdr));
        Vec3 dchidq = (iota*B[particle]).cross(iota)*(-4*rInv2);
        double (&g12)[3][3] = G12.v;
        double (&a)[3][3] = A[particle].v;
        ParticleInfo& p = particles[particle];
        Vec3 scale = Vec3(p.rx*p.rx, p.ry*p.ry, p.rz*p.rz)*(-0.5*eta/detG12);
        Matrix D;
        double (&d)[3][3] = D.v;
        d[0][0] = scale[0]*(2*a[0][0]*(g12[1][1]*g12[2][2] - g12[1][2]*g12[2][1]) +
                              a[0][2]*(g12[1][2]*g12[0][1] + g12[1][0]*g12[2][1] - g12[1][1]*(g12[0][2] + g12[2][0])) +
                              a[0][1]*(g12[0][2]*g12[2][1] + g12[2][0]*g12[1][2] - g12[2][2]*(g12[0][1] + g12[1][0])));
        d[0][1] = scale[0]*(  a[0][0]*(g12[0][2]*g12[2][1] + g12[2][0]*g12[1][2] - g12[2][2]*(g12[0][1] + g12[1][0])) +
                            2*a[0][1]*(g12[0][0]*g12[2][2] - g12[2][0]*g12[0][2]) +
                              a[0][2]*(g12[1][0]*g12[0][2] + g12[2][0]*g12[0][1] - g12[0][0]*(g12[1][2] + g12[2][1])));
        d[0][2] = scale[0]*(  a[0][0]*(g12[0][1]*g12[1][2] + g12[1][0]*g12[2][1] - g12[1][1]*(g12[0][2] + g12[2][0])) +
                              a[0][1]*(g12[1][0]*g12[0][2] + g12[2][0]*g12[0][1] - g12[0][0]*(g12[1][2] + g12[2][1])) +
                            2*a[0][2]*(g12[1][1]*g12[0][0] - g12[1][0]*g12[0][1]));
        d[1][0] = scale[1]*(2*a[1][0]*(g12[1][1]*g12[2][2] - g12[1][2]*g12[2][1]) +
                              a[1][1]*(g12[0][2]*g12[2][1] + g12[2][0]*g12[1][2] - g12[2][2]*(g12[0][1] + g12[1][0])) +
                              a[1][2]*(g12[1][2]*g12[0][1] + g12[1][0]*g12[2][1] - g12[1][1]*(g12[0][2] + g12[2][0])));
        d[1][1] = scale[1]*(  a[1][0]*(g12[0][2]*g12[2][1] + g12[2][0]*g12[1][2] - g12[2][2]*(g12[0][1] + g12[1][0])) +
                            2*a[1][1]*(g12[2][2]*g12[0][0] - g12[2][0]*g12[0][2]) +
                              a[1][2]*(g12[1][0]*g12[0][2] + g12[0][1]*g12[2][0] - g12[0][0]*(g12[1][2] + g12[2][1])));
        d[1][2] = scale[1]*(  a[1][0]*(g12[0][1]*g12[1][2] + g12[1][0]*g12[2][1] - g12[1][1]*(g12[0][2] + g12[2][0])) +
                              a[1][1]*(g12[1][0]*g12[0][2] + g12[0][1]*g12[2][0] - g12[0][0]*(g12[1][2] + g12[2][1])) +
                            2*a[1][2]*(g12[1][1]*g12[0][0] - g12[1][0]*g12[0][1]));
        d[2][0] = scale[2]*(2*a[2][0]*(g12[1][1]*g12[2][2] - g12[2][1]*g12[1][2]) +
                              a[2][1]*(g12[0][2]*g12[2][1] + g12[1][2]*g12[2][0] - g12[2][2]*(g12[0][1] + g12[1][0])) +
                              a[2][2]*(g12[0][1]*g12[1][2] + g12[2][1]*g12[1][0] - g12[1][1]*(g12[0][2] + g12[2][0])));
        d[2][1] = scale[2]*(  a[2][0]*(g12[0][2]*g12[2][1] + g12[1][2]*g12[2][0] - g12[2][2]*(g12[0][1] + g12[1][0])) +
                            2*a[2][1]*(g12[0][0]*g12[2][2] - g12[0][2]*g12[2][0]) +
                              a[2][2]*(g12[1][0]*g12[0][2] + g12[0][1]*g12[2][0] - g12[0][0]*(g12[1][2] + g12[2][1])));
        d[2][2] = scale[2]*(  a[2][0]*(g12[0][1]*g12[1][2] + g12[2][1]*g12[1][0] - g12[1][1]*(g12[0][2] + g12[2][0])) +
                              a[2][1]*(g12[1][0]*g12[0][2] + g12[2][0]*g12[0][1] - g12[0][0]*(g12[1][2] + g12[2][1])) +
                            2*a[2][2]*(g12[1][1]*g12[0][0] - g12[1][0]*g12[0][1]));
        Vec3 detadq;
        for (int i = 0; i < 3; i++)
            detadq += Vec3(a[i][0], a[i][1], a[i][2]).cross(Vec3(d[i][0], d[i][1], d[i][2]));
        Vec3 torque = (dchidq*(u*eta) + detadq*(u*chi) + dudq*(eta*chi))*switchValue;
        torques[particle] -= torque;
    }
    return switchValue*energy;
}
