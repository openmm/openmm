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

#include "ReferenceOrientationRestraintForce.h"
#include "jama_eig.h"

using namespace OpenMM;
using namespace std;

ReferenceOrientationRestraintForce::ReferenceOrientationRestraintForce(double k, vector<Vec3>& referencePos, vector<int>& particles) :
        k(k), referencePos(referencePos), particles(particles) {
}

ReferenceOrientationRestraintForce::~ReferenceOrientationRestraintForce() {
}

double ReferenceOrientationRestraintForce::calculateIxn(vector<Vec3>& atomCoordinates, vector<Vec3>& forces) const {
    // Find the optimal transformation using the algorithm described in Coutsias et al,
    // "Using quaternions to calculate RMSD" (doi: 10.1002/jcc.20110).  First subtract
    // the centroid from the atom positions.  The reference positions have already been centered.

    int numParticles = particles.size();
    Vec3 center;
    for (int i : particles)
        center += atomCoordinates[i];
    center /= numParticles;
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numParticles; i++)
        positions[i] = atomCoordinates[particles[i]]-center;

    // Compute the correlation matrix.

    double R[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < numParticles; k++) {
                int index = particles[k];
                R[i][j] += positions[k][j]*referencePos[index][i];
            }

    // Compute the F matrix.

    Array2D<double> F(4, 4);
    F[0][0] =  R[0][0] + R[1][1] + R[2][2];
    F[1][0] =  R[1][2] - R[2][1];
    F[2][0] =  R[2][0] - R[0][2];
    F[3][0] =  R[0][1] - R[1][0];

    F[0][1] =  R[1][2] - R[2][1];
    F[1][1] =  R[0][0] - R[1][1] - R[2][2];
    F[2][1] =  R[0][1] + R[1][0];
    F[3][1] =  R[0][2] + R[2][0];

    F[0][2] =  R[2][0] - R[0][2];
    F[1][2] =  R[0][1] + R[1][0];
    F[2][2] = -R[0][0] + R[1][1] - R[2][2];
    F[3][2] =  R[1][2] + R[2][1];

    F[0][3] =  R[0][1] - R[1][0];
    F[1][3] =  R[0][2] + R[2][0];
    F[2][3] =  R[1][2] + R[2][1];
    F[3][3] = -R[0][0] - R[1][1] + R[2][2];

    // Find the maximum eigenvalue and eigenvector.

    JAMA::Eigenvalue<double> eigen(F);
    Array1D<double> values;
    eigen.getRealEigenvalues(values);
    Array2D<double> vectors;
    eigen.getV(vectors);

    // Construct the quaternion and use it to compute the energy.

    double q[] = {vectors[0][3], vectors[1][3], vectors[2][3], vectors[3][3]};
    double energy = 2*k*(1.0-q[0]*q[0]);

    // Compute the forces.  The algorithm is modelled after the calculation of the
    // orientationAngle CV in Colvars.  I can't find any published source for it.

    if (q[0]*q[0] < 1.0) {
        double theta = 2*asin(sqrt(1.0-q[0]*q[0]));
        double dxdq = 4.0*k*sin(theta/2)*cos(theta/2)/sqrt(1.0-q[0]*q[0]);
        if (vectors[0][3] > 0)
            dxdq = -dxdq;
        for (int index = 0; index < numParticles; index++) {
            const Vec3& p = referencePos[particles[index]];
            Vec3 ds[4][4] = {
                {Vec3(p[0], p[1], p[2]), Vec3(0.0, -p[2], p[1]), Vec3(p[2], 0.0, -p[0]), Vec3(-p[1], p[0], 0.0)},
                {Vec3(0.0, -p[2], p[1]), Vec3(p[0], -p[1], -p[2]), Vec3(p[1], p[0], 0.0), Vec3(p[2], 0.0, p[0])},
                {Vec3(p[2], 0.0, -p[0]), Vec3(p[1], p[0], 0.0), Vec3(-p[0], p[1], -p[2]), Vec3(0.0, p[2], p[1])},
                {Vec3(-p[1], p[0], 0.0), Vec3(p[2], 0.0, p[0]), Vec3(0.0, p[2], p[1]), Vec3(-p[0], -p[1], p[2])}
            };
            Vec3 dq;
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    dq += ((vectors[i][2]*vectors[0][2]) / (values[3]-values[2]) +
                           (vectors[i][1]*vectors[0][1]) / (values[3]-values[1]) +
                           (vectors[i][0]*vectors[0][0]) / (values[3]-values[0])) * vectors[j][3] * ds[i][j];
                }
            }
            forces[particles[index]] -= dxdq*dq;
        }
    }
    return energy;
}
