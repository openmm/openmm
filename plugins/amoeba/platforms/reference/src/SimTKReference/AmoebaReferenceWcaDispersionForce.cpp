
/* Portions copyright (c) 2006 Stanford University and Simbios.
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

#ifdef WIN32
#define _USE_MATH_DEFINES // Needed to get M_PI
#endif

#include "AmoebaReferenceForce.h"
#include "AmoebaReferenceWcaDispersionForce.h"
#include <cmath>

using std::vector;
using namespace OpenMM;

AmoebaReferenceWcaDispersionForce::AmoebaReferenceWcaDispersionForce(double epso, double epsh, double rmino,
                                                                     double rminh, double awater, double shctd,
                                                                     double dispoff, double slevy) :
        _epso(epso), _epsh(epsh), _rmino(rmino), _rminh(rminh), _awater(awater), _shctd(shctd), _dispoff(dispoff),
        _slevy(slevy) {
}

static double integralBeforeRMin(double eps, double r, double r2, double sk2,
                                 double lik2, double lik3, double lik4,
                                 double uik2, double uik3, double uik4) {
    return -eps * (4.0e0 * M_PI / (48.0e0 * r) *
                   (3.0e0 * (lik4 - uik4) - 8.0e0 * r * (lik3 - uik3) + 6.0e0 * (r2 - sk2) * (lik2 - uik2)));
}

static double
integralBeforeRminDerivative(double ri, double eps, double rmin, double r, double r2, double r3, double sk, double sk2,
                             double lik, double lik2, double lik3, double uik, double uik2, double uik3) {
    double dl;
    if (ri > r - sk) {
        dl = (-lik2 + 2.0 * r2 + 2.0 * sk2) * lik2;
    } else {
        dl = (-lik3 + 4.0 * lik2 * r - 6.0 * lik * r2 + 2.0 * lik * sk2 + 4.0 * r3 - 4.0 * r * sk2) * lik;
    }
    double du;
    if (r + sk > rmin) {
        du = -(-uik2 + 2.0 * r2 + 2.0 * sk2) * uik2;
    } else {
        du = -(-uik3 + 4.0 * uik2 * r - 6.0 * uik * r2 + 2.0 * uik * sk2 + 4.0 * r3 - 4.0 * r * sk2) * uik;
    }
    return -eps * M_PI * (dl + du) / (4.0 * r2);
}

static double integratlAfterRmin(double eps, double rmin7, double r, double r2, double sk2,
                                 double lik, double lik2, double lik3, double lik4, double lik5, double lik10,
                                 double lik11, double lik12, double uik, double uik2, double uik3, double uik4,
                                 double uik5, double uik10, double uik11, double uik12) {
    double er7 = eps * rmin7;
    double term = 4.0 * M_PI / (120.0 * r * lik5 * uik5)
                  * (15.0 * uik * lik * r * (uik4 - lik4)
                     - 10.0 * uik2 * lik2 * (uik3 - lik3)
                     + 6.0 * (sk2 - r2) * (uik5 - lik5));
    double term2 = 4.0 * M_PI / (2640.0 * r * lik12 * uik12)
                   * (120.0 * uik * lik * r * (uik11 - lik11)
                      - 66.0 * uik2 * lik2 * (uik10 - lik10)
                      + 55.0 * (sk2 - r2) * (uik12 - lik12));
    double idisp = -2.0 * er7 * term;
    double irep = er7 * rmin7 * term2;
    return irep + idisp;
}

static double integratlAfterRminDerivative(double ri, double eps, double rmin, double rmin7, double rmax,
                                           double r, double r2, double r3, double sk, double sk2, double lik,
                                           double lik2, double lik3, double lik5, double lik6, double lik12,
                                           double lik13, double uik, double uik2, double uik3, double uik6,
                                           double uik13) {
    double er7 = eps * rmin7;
    double lowerTerm = lik2 * r + r3 - r * sk2;
    double upperTerm = uik2 * r + r3 - r * sk2;

    double dl;
    if (ri > r - sk || rmax < rmin) {
        dl = -(-5.0 * lik2 + 3.0 * r2 + 3.0 * sk2) / lik5;
    } else {
        dl = (5.0 * lik3 - 33.0 * lik * r2 - 3.0 * lik * sk2 + 15.0 * lowerTerm) / lik6;
    }
    double du = -(5.0 * uik3 - 33.0 * uik * r2 - 3.0 * uik * sk2 + 15.0 * upperTerm) / uik6;
    double de = -2.0 * M_PI * er7 * (dl + du) / (15.0 * r2);

    if (ri > r - sk || rmax < rmin) {
        dl = -(-6.0 * lik2 + 5.0 * r2 + 5.0 * sk2) / lik12;
    } else {
        dl = (6.0 * lik3 - 125.0 * lik * r2 - 5.0 * lik * sk2 + 60.0 * lowerTerm) / lik13;
    }
    du = -(6.0 * uik3 - 125.0 * uik * r2 - 5.0 * uik * sk2 + 60.0 * upperTerm) / uik13;
    de += M_PI * er7 * rmin7 * (dl + du) / (60.0 * r2);

    return de;
}

static double interact(double factor, double ri, double sk, double rmix, double emix,
                       double r, double r2, double r3, Vec3 &force) {
    double sum = 0.0;
    // Nothing to do if the integral begins beyond r + sk (i.e. atom k does not exclude solvent)
    if (ri < r + sk) {
        // Zero out the derivative contribution of atom k.
        double de = 0.0;
        double sk2 = sk * sk;
        // Compute the maximum of 1) the beginning of the integral and 2) closest edge of atom K.
        double iStart = ri > r - sk ? ri : r - sk;
        // Use this as the lower limit for integrating the constant eps value below Rmin.
        double lik = iStart;
        // Interaction with water from lik to Rmin; nothing to do if the lower limit is greater than Rmin.
        if (lik < rmix) {
            double lik2 = lik * lik;
            double lik3 = lik2 * lik;
            double lik4 = lik3 * lik;
            // Upper limit is the minimum of Rmin and the farthest edge of atom K.
            double uik = r + sk < rmix ? r + sk : rmix;
            double uik2 = uik * uik;
            double uik3 = uik2 * uik;
            double uik4 = uik3 * uik;
            sum = integralBeforeRMin(emix, r, r2, sk2, lik2, lik3, lik4, uik2, uik3, uik4);
            de = integralBeforeRminDerivative(ri, emix, rmix, r, r2, r3, sk, sk2, lik, lik2, lik3, uik, uik2, uik3);
        }
        // Upper limit the variable part of Uwca always the farthest edge of atom K.
        double uik = r + sk;
        // Interaction with water beyond Rmin, from lik to uik = r + sk.
        if (uik > rmix) {
            // Start the integral at the max of 1) iStart and 2) Rmin.
            lik = iStart > rmix ? iStart : rmix;
            double lik2 = lik * lik;
            double lik3 = lik2 * lik;
            double lik4 = lik3 * lik;
            double lik5 = lik4 * lik;
            double lik6 = lik5 * lik;
            double lik10 = lik5 * lik5;
            double lik11 = lik10 * lik;
            double lik12 = lik11 * lik;
            double uik2 = uik * uik;
            double uik3 = uik2 * uik;
            double uik4 = uik3 * uik;
            double uik5 = uik4 * uik;
            double uik10 = uik5 * uik5;
            double uik11 = uik10 * uik;
            double uik12 = uik11 * uik;
            double rmix3 = rmix * rmix * rmix;
            double rmix7 = rmix3 * rmix3 * rmix;
            sum += integratlAfterRmin(emix, rmix7, r, r2, sk2, lik, lik2, lik3, lik4, lik5, lik10, lik11, lik12, uik,
                                      uik2, uik3, uik4, uik5, uik10, uik11, uik12);
            double lik13 = lik12 * lik;
            double uik6 = uik5 * uik;
            double uik13 = uik12 * uik;
            de += integratlAfterRminDerivative(ri, emix, rmix, rmix7, iStart, r, r2, r3, sk, sk2, lik, lik2, lik3,
                                               lik5, lik6, lik12, lik13, uik, uik2, uik3, uik6, uik13);
        }
        // Increment the individual dispersion gradient components.
        de *= factor / r;
        force[0] += de;
        force[1] += de;
        force[2] += de;
    }
    return factor * sum;
}

double AmoebaReferenceWcaDispersionForce::calculatePairIxn(double radiusJ,
                                                           const Vec3 &particleIPosition,
                                                           const Vec3 &particleJPosition,
                                                           const double *const intermediateValues,
                                                           Vec3 &force) const {
    // Atomic separation
    double xr = particleIPosition[0] - particleJPosition[0];
    double yr = particleIPosition[1] - particleJPosition[1];
    double zr = particleIPosition[2] - particleJPosition[2];
    double r2 = xr * xr + yr * yr + zr * zr;
    double r = sqrt(r2);
    double r3 = r2 * r;

    // Parameters for atom i with water oxygen.
    double rmixo = intermediateValues[RMIXO];
    double emixo = intermediateValues[EMIXO];

    // Parameters for atom i with water hydrogen.
    double rmixh = intermediateValues[RMIXH];
    double emixh = intermediateValues[EMIXH];

    // Start of integration of dispersion for atom i with water oxygen.
    double riO = rmixo / 2.0e0 + _dispoff;
    double nO = 1.0e0;

    // Start of integration of dispersion for atom i with water hydrogen.
    double riH = rmixh / 2.0e0 + _dispoff;
    double nH = 2.0e0;

    // Atom k blocks the interaction of atom i with solvent.
    double sJ = radiusJ * _shctd;

    double sum = interact(nO, riO, sJ, rmixo, emixo, r, r2, r3, force) +
                 interact(nH, riH, sJ, rmixh, emixh, r, r2, r3, force);

    force[0] *= xr;
    force[1] *= yr;
    force[2] *= zr;

    return sum;
}

double AmoebaReferenceWcaDispersionForce::calculateForceAndEnergy(int numParticles,
                                                                  const vector<Vec3> &particlePositions,
                                                                  const std::vector<double> &radii,
                                                                  const std::vector<double> &epsilons,
                                                                  double totalMaximumDispersionEnergy,
                                                                  vector<Vec3> &forces) const {
    double energy = 0.0;
    double rmino2 = _rmino * _rmino;
    double rmino3 = rmino2 * _rmino;
    double rminh2 = _rminh * _rminh;
    double rminh3 = rminh2 * _rminh;
    double intermediateValues[LastIntermediateValueIndex];

    // loop over all ixns
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(numParticles); ii++) {

        double epsi = epsilons[ii];
        double rmini = radii[ii];
        double rminI2 = rmini * rmini;
        double rminI3 = rminI2 * rmini;

        double denominator = sqrt(_epso) + sqrt(epsi);
        double emixo = 4.0 * _epso * epsi / (denominator * denominator);
        intermediateValues[EMIXO] = emixo;

        double rmixo = 2.0 * (rmino3 + rminI3) / (rmino2 + rminI2);
        intermediateValues[RMIXO] = rmixo;

        denominator = sqrt(_epsh) + sqrt(epsi);
        double emixh = 4.0 * _epsh * epsi / (denominator * denominator);
        intermediateValues[EMIXH] = emixh;

        double rmixh = 2.0 * (rminh3 + rminI3) / (rminh2 + rminI2);
        intermediateValues[RMIXH] = rmixh;

        // Remove dispersion for atom i by atom j.
        for (unsigned int jj = 0; jj < static_cast<unsigned int>(numParticles); jj++) {

            if (ii == jj) continue;

            Vec3 force = Vec3();
            energy += calculatePairIxn(radii[jj], particlePositions[ii], particlePositions[jj],
                                       intermediateValues, force);

            force[0] *= _slevy * _awater;
            force[1] *= _slevy * _awater;
            force[2] *= _slevy * _awater;

            forces[ii][0] += force[0];
            forces[ii][1] += force[1];
            forces[ii][2] += force[2];

            forces[jj][0] -= force[0];
            forces[jj][1] -= force[1];
            forces[jj][2] -= force[2];

        }

    }

    energy = totalMaximumDispersionEnergy - _slevy * _awater * energy;

    return energy;
}
