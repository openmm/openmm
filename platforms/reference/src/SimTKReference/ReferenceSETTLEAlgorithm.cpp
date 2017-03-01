/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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

#include "ReferenceSETTLEAlgorithm.h"

using namespace OpenMM;
using namespace std;

ReferenceSETTLEAlgorithm::ReferenceSETTLEAlgorithm(const vector<int>& atom1, const vector<int>& atom2, const vector<int>& atom3,
        const vector<double>& distance1, const vector<double>& distance2, vector<double>& masses) :
        atom1(atom1), atom2(atom2), atom3(atom3), distance1(distance1), distance2(distance2), masses(masses) {
}

int ReferenceSETTLEAlgorithm::getNumClusters() const {
    return atom1.size();
}

void ReferenceSETTLEAlgorithm::getClusterParameters(int index, int& atom1, int& atom2, int& atom3, double& distance1, double& distance2) const {
    atom1 = this->atom1[index];
    atom2 = this->atom2[index];
    atom3 = this->atom3[index];
    distance1 = this->distance1[index];
    distance2 = this->distance2[index];
}

void ReferenceSETTLEAlgorithm::apply(vector<OpenMM::Vec3>& atomCoordinates, vector<OpenMM::Vec3>& atomCoordinatesP, vector<double>& inverseMasses, double tolerance) {
    for (int index = 0; index < (int) atom1.size(); ++index) {
        Vec3 apos0 = atomCoordinates[atom1[index]];
        Vec3 xp0 = atomCoordinatesP[atom1[index]]-apos0;
        Vec3 apos1 = atomCoordinates[atom2[index]];
        Vec3 xp1 = atomCoordinatesP[atom2[index]]-apos1;
        Vec3 apos2 = atomCoordinates[atom3[index]];
        Vec3 xp2 = atomCoordinatesP[atom3[index]]-apos2;
        double m0 = masses[atom1[index]];
        double m1 = masses[atom2[index]];
        double m2 = masses[atom3[index]];

        // Apply the SETTLE algorithm.

        double xb0 = apos1[0]-apos0[0];
        double yb0 = apos1[1]-apos0[1];
        double zb0 = apos1[2]-apos0[2];
        double xc0 = apos2[0]-apos0[0];
        double yc0 = apos2[1]-apos0[1];
        double zc0 = apos2[2]-apos0[2];

        double invTotalMass = 1/(m0+m1+m2);
        double xcom = (xp0[0]*m0 + (xb0+xp1[0])*m1 + (xc0+xp2[0])*m2) * invTotalMass;
        double ycom = (xp0[1]*m0 + (yb0+xp1[1])*m1 + (yc0+xp2[1])*m2) * invTotalMass;
        double zcom = (xp0[2]*m0 + (zb0+xp1[2])*m1 + (zc0+xp2[2])*m2) * invTotalMass;

        double xa1 = xp0[0] - xcom;
        double ya1 = xp0[1] - ycom;
        double za1 = xp0[2] - zcom;
        double xb1 = xb0 + xp1[0] - xcom;
        double yb1 = yb0 + xp1[1] - ycom;
        double zb1 = zb0 + xp1[2] - zcom;
        double xc1 = xc0 + xp2[0] - xcom;
        double yc1 = yc0 + xp2[1] - ycom;
        double zc1 = zc0 + xp2[2] - zcom;

        double xaksZd = yb0*zc0 - zb0*yc0;
        double yaksZd = zb0*xc0 - xb0*zc0;
        double zaksZd = xb0*yc0 - yb0*xc0;
        double xaksXd = ya1*zaksZd - za1*yaksZd;
        double yaksXd = za1*xaksZd - xa1*zaksZd;
        double zaksXd = xa1*yaksZd - ya1*xaksZd;
        double xaksYd = yaksZd*zaksXd - zaksZd*yaksXd;
        double yaksYd = zaksZd*xaksXd - xaksZd*zaksXd;
        double zaksYd = xaksZd*yaksXd - yaksZd*xaksXd;

        double axlng = sqrt(xaksXd*xaksXd + yaksXd*yaksXd + zaksXd*zaksXd);
        double aylng = sqrt(xaksYd*xaksYd + yaksYd*yaksYd + zaksYd*zaksYd);
        double azlng = sqrt(xaksZd*xaksZd + yaksZd*yaksZd + zaksZd*zaksZd);
        double trns11 = xaksXd / axlng;
        double trns21 = yaksXd / axlng;
        double trns31 = zaksXd / axlng;
        double trns12 = xaksYd / aylng;
        double trns22 = yaksYd / aylng;
        double trns32 = zaksYd / aylng;
        double trns13 = xaksZd / azlng;
        double trns23 = yaksZd / azlng;
        double trns33 = zaksZd / azlng;

        double xb0d = trns11*xb0 + trns21*yb0 + trns31*zb0;
        double yb0d = trns12*xb0 + trns22*yb0 + trns32*zb0;
        double xc0d = trns11*xc0 + trns21*yc0 + trns31*zc0;
        double yc0d = trns12*xc0 + trns22*yc0 + trns32*zc0;
        double za1d = trns13*xa1 + trns23*ya1 + trns33*za1;
        double xb1d = trns11*xb1 + trns21*yb1 + trns31*zb1;
        double yb1d = trns12*xb1 + trns22*yb1 + trns32*zb1;
        double zb1d = trns13*xb1 + trns23*yb1 + trns33*zb1;
        double xc1d = trns11*xc1 + trns21*yc1 + trns31*zc1;
        double yc1d = trns12*xc1 + trns22*yc1 + trns32*zc1;
        double zc1d = trns13*xc1 + trns23*yc1 + trns33*zc1;

        //                                        --- Step2  A2' ---

        double rc = 0.5*distance2[index];
        double rb = sqrt(distance1[index]*distance1[index]-rc*rc);
        double ra = rb*(m1+m2)*invTotalMass;
        rb -= ra;
        double sinphi = za1d / ra;
        double cosphi = sqrt(1 - sinphi*sinphi);
        double sinpsi = (zb1d - zc1d) / (2*rc*cosphi);
        double cospsi = sqrt(1 - sinpsi*sinpsi);

        double ya2d =   ra*cosphi;
        double xb2d = - rc*cospsi;
        double yb2d = - rb*cosphi - rc*sinpsi*sinphi;
        double yc2d = - rb*cosphi + rc*sinpsi*sinphi;
        double xb2d2 = xb2d*xb2d;
        double hh2 = 4.0f*xb2d2 + (yb2d-yc2d)*(yb2d-yc2d) + (zb1d-zc1d)*(zb1d-zc1d);
        double deltx = 2.0f*xb2d + sqrt(4.0f*xb2d2 - hh2 + distance2[index]*distance2[index]);
        xb2d -= deltx*0.5;

        //                                        --- Step3  al,be,ga ---

        double alpha = (xb2d*(xb0d-xc0d) + yb0d*yb2d + yc0d*yc2d);
        double beta = (xb2d*(yc0d-yb0d) + xb0d*yb2d + xc0d*yc2d);
        double gamma = xb0d*yb1d - xb1d*yb0d + xc0d*yc1d - xc1d*yc0d;

        double al2be2 = alpha*alpha + beta*beta;
        double sintheta = (alpha*gamma - beta*sqrt(al2be2 - gamma*gamma)) / al2be2;

        //                                        --- Step4  A3' ---

        double costheta = sqrt(1 - sintheta*sintheta);
        double xa3d = - ya2d*sintheta;
        double ya3d =   ya2d*costheta;
        double za3d = za1d;
        double xb3d =   xb2d*costheta - yb2d*sintheta;
        double yb3d =   xb2d*sintheta + yb2d*costheta;
        double zb3d = zb1d;
        double xc3d = - xb2d*costheta - yc2d*sintheta;
        double yc3d = - xb2d*sintheta + yc2d*costheta;
        double zc3d = zc1d;

        //                                        --- Step5  A3 ---

        double xa3 = trns11*xa3d + trns12*ya3d + trns13*za3d;
        double ya3 = trns21*xa3d + trns22*ya3d + trns23*za3d;
        double za3 = trns31*xa3d + trns32*ya3d + trns33*za3d;
        double xb3 = trns11*xb3d + trns12*yb3d + trns13*zb3d;
        double yb3 = trns21*xb3d + trns22*yb3d + trns23*zb3d;
        double zb3 = trns31*xb3d + trns32*yb3d + trns33*zb3d;
        double xc3 = trns11*xc3d + trns12*yc3d + trns13*zc3d;
        double yc3 = trns21*xc3d + trns22*yc3d + trns23*zc3d;
        double zc3 = trns31*xc3d + trns32*yc3d + trns33*zc3d;

        xp0[0] = xcom + xa3;
        xp0[1] = ycom + ya3;
        xp0[2] = zcom + za3;
        xp1[0] = xcom + xb3 - xb0;
        xp1[1] = ycom + yb3 - yb0;
        xp1[2] = zcom + zb3 - zb0;
        xp2[0] = xcom + xc3 - xc0;
        xp2[1] = ycom + yc3 - yc0;
        xp2[2] = zcom + zc3 - zc0;

        // Record the new positions.

        atomCoordinatesP[atom1[index]] = xp0+apos0;
        atomCoordinatesP[atom2[index]] = xp1+apos1;
        atomCoordinatesP[atom3[index]] = xp2+apos2;
    }
}

void ReferenceSETTLEAlgorithm::applyToVelocities(vector<OpenMM::Vec3>& atomCoordinates, vector<OpenMM::Vec3>& velocities, vector<double>& inverseMasses, double tolerance) {
    for (int index = 0; index < (int) atom1.size(); ++index) {
        Vec3 apos0 = atomCoordinates[atom1[index]];
        Vec3 apos1 = atomCoordinates[atom2[index]];
        Vec3 apos2 = atomCoordinates[atom3[index]];
        Vec3 v0 = velocities[atom1[index]];
        Vec3 v1 = velocities[atom2[index]];
        Vec3 v2 = velocities[atom3[index]];
        
        // Compute intermediate quantities: the atom masses, the bond directions, the relative velocities,
        // and the angle cosines and sines.
        
        double mA = masses[atom1[index]];
        double mB = masses[atom2[index]];
        double mC = masses[atom3[index]];
        Vec3 eAB = apos1-apos0;
        Vec3 eBC = apos2-apos1;
        Vec3 eCA = apos0-apos2;
        eAB /= sqrt(eAB[0]*eAB[0] + eAB[1]*eAB[1] + eAB[2]*eAB[2]);
        eBC /= sqrt(eBC[0]*eBC[0] + eBC[1]*eBC[1] + eBC[2]*eBC[2]);
        eCA /= sqrt(eCA[0]*eCA[0] + eCA[1]*eCA[1] + eCA[2]*eCA[2]);
        double vAB = (v1[0]-v0[0])*eAB[0] + (v1[1]-v0[1])*eAB[1] + (v1[2]-v0[2])*eAB[2];
        double vBC = (v2[0]-v1[0])*eBC[0] + (v2[1]-v1[1])*eBC[1] + (v2[2]-v1[2])*eBC[2];
        double vCA = (v0[0]-v2[0])*eCA[0] + (v0[1]-v2[1])*eCA[1] + (v0[2]-v2[2])*eCA[2];
        double cA = -(eAB[0]*eCA[0] + eAB[1]*eCA[1] + eAB[2]*eCA[2]);
        double cB = -(eAB[0]*eBC[0] + eAB[1]*eBC[1] + eAB[2]*eBC[2]);
        double cC = -(eBC[0]*eCA[0] + eBC[1]*eCA[1] + eBC[2]*eCA[2]);
        double s2A = 1-cA*cA;
        double s2B = 1-cB*cB;
        double s2C = 1-cC*cC;
        
        // Solve the equations.  These are different from those in the SETTLE paper (JCC 13(8), pp. 952-962, 1992), because
        // in going from equations B1 to B2, they make the assumption that mB=mC (but don't bother to mention they're
        // making that assumption).  We allow all three atoms to have different masses.
        
        double mABCinv = 1/(mA*mB*mC);
        double denom = (((s2A*mB+s2B*mA)*mC+(s2A*mB*mB+2*(cA*cB*cC+1)*mA*mB+s2B*mA*mA))*mC+s2C*mA*mB*(mA+mB))*mABCinv;
        double tab = ((cB*cC*mA-cA*mB-cA*mC)*vCA + (cA*cC*mB-cB*mC-cB*mA)*vBC + (s2C*mA*mA*mB*mB*mABCinv+(mA+mB+mC))*vAB)/denom;
        double tbc = ((cA*cB*mC-cC*mB-cC*mA)*vCA + (s2A*mB*mB*mC*mC*mABCinv+(mA+mB+mC))*vBC + (cA*cC*mB-cB*mA-cB*mC)*vAB)/denom;
        double tca = ((s2B*mA*mA*mC*mC*mABCinv+(mA+mB+mC))*vCA + (cA*cB*mC-cC*mB-cC*mA)*vBC + (cB*cC*mA-cA*mB-cA*mC)*vAB)/denom;
        v0 += (eAB*tab - eCA*tca)*inverseMasses[atom1[index]];
        v1 += (eBC*tbc - eAB*tab)*inverseMasses[atom2[index]];
        v2 += (eCA*tca - eBC*tbc)*inverseMasses[atom3[index]];
        velocities[atom1[index]] = v0;
        velocities[atom2[index]] = v1;
        velocities[atom3[index]] = v2;
    }
}
