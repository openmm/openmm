
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

AmoebaReferenceWcaDispersionForce::AmoebaReferenceWcaDispersionForce(double epso, double epsh, double rmino, double rminh, 
                                                                     double awater, double shctd, double dispoff, double slevy) :
                               _epso(epso), _epsh(epsh), _rmino(rmino), _rminh(rminh), _awater(awater), _shctd(shctd), _dispoff(dispoff), _slevy(slevy) {
}   


double AmoebaReferenceWcaDispersionForce::calculatePairIxn(double radiusI, double radiusK, 
                                                           const Vec3& particleIPosition,
                                                           const Vec3& particleJPosition,
                                                           const double* const intermediateValues,
                                                           Vec3& force) const {

    static const double PI = M_PI;

    double xr           = particleIPosition[0] - particleJPosition[0];
    double yr           = particleIPosition[1] - particleJPosition[1];
    double zr           = particleIPosition[2] - particleJPosition[2];

    double r2           = xr*xr + yr*yr + zr*zr;
    double r            = sqrt(r2);
    double r3           = r2*r;

    double sK           = radiusK*_shctd;
    double sK2          = sK*sK;

    double rmixo        = intermediateValues[RMIXO];
    double rmixo7       = intermediateValues[RMIXO7];

    double emixo        = intermediateValues[EMIXO];

    double rmixh        = intermediateValues[RMIXH];
    double rmixh7       = intermediateValues[RMIXH7];

    double emixh        = intermediateValues[EMIXH];

    double ao           = intermediateValues[AO];
    double ah           = intermediateValues[AH];

    double sum          = 0.0;
    double de           = 0.0;

    if (radiusI < (r + sK)) {

        double rmax     = (radiusI > (r - sK)) ? radiusI : (r - sK);

        double lik      = rmax;
        double lik2     = lik*lik;
        double lik3     = lik2*lik;
        double lik4     = lik2*lik2;

        if (lik < rmixo) { 

            double uik  = (r + sK) < rmixo ? (r + sK) : rmixo;
            double uik2 = uik*uik;
            double uik3 = uik2*uik;
            double uik4 = uik2*uik2;

            double term = 4.0*PI/(48.0*r)* (3.0*(lik4-uik4) - 8.0*r*(lik3-uik3) + 6.0*(r2-sK2)*(lik2-uik2));

            double dl;
            if (radiusI >  (r - sK)) {
                dl  = -lik2 + 2.0*(r2 + sK2);
                dl *= lik2;
            } else {
                dl  = -lik3 + 4.0*lik2*r - 6.0*lik*r2 + 2.0*lik*sK2 + 4.0*r*(r2 - sK2);
                dl *= lik;
            }

            double du;
            if ((r+sK) > rmixo) {
                du  = -uik2 + 2.0*(r2 + sK2);
                du *= -uik2;
            } else {
                du  = -uik3 + 4.0*uik2*r - 6.0*uik*r2 + 2.0*uik*sK2 + 4.0*r*(r2 - sK2);
                du *= -uik;
            }
            de     = -emixo*PI*(dl+du)/(4.0*r2);
            sum   += -emixo*term;
        }

        if (lik < rmixh) {

            double uik  = (r + sK) < rmixh ? (r + sK) : rmixh;
            double uik2 = uik*uik;
            double uik3 = uik2*uik;
            double uik4 = uik2*uik2;
            double term = 4.0*PI / (48.0*r)*(3.0*(lik4-uik4) - 8.0*r*(lik3-uik3) + 6.0*(r2-sK2)*(lik2-uik2));
            double dl;
            if (radiusI > (r-sK)) {
                dl  = -lik2 + 2.0*(r2 + sK2);
                dl *= lik2;
            } else {
                dl  = -lik3 + 4.0*lik2*r - 6.0*lik*r2 + 2.0*lik*sK2 + 4.0*r*(r2 -sK2);
                dl *= lik;
            }
      
            double du;
            if (r+sK > rmixh) {
                du  = -uik2 + 2.0*(r2 + sK2);
                du *= -uik2;
            } else {
                du  = -uik3 +4.0*uik2*r - 6.0*uik*r2 + 2.0*uik*sK2 +4.0*r*(r2 - sK2);
                du *= -uik;
            }
            de   -= 2.0*emixh*PI*(dl+du)/(4.0*r2);
            sum  -= 2.0*emixh*term;
        }

        double uik   = r + sK;
        double uik2  = uik   * uik;
        double uik3  = uik2  * uik;
        double uik4  = uik2  * uik2;
        double uik5  = uik3  * uik2;
        double uik6  = uik3  * uik3;
        double uik10 = uik5  * uik5;
        double uik11 = uik5  * uik6;
        double uik12 = uik6  * uik6;
        double uik13 = uik10 * uik3;

        if (uik > rmixo) {

            double lik   = rmax > rmixo ? rmax : rmixo;
            double lik2  = lik   * lik;
            double lik3  = lik2  * lik;
            double lik4  = lik2  * lik2;
            double lik5  = lik2  * lik3;
            double lik6  = lik3  * lik3;
            double lik10 = lik5  * lik5;
            double lik11 = lik5  * lik6;
            double lik12 = lik6  * lik6;
            double lik13 = lik10 * lik3;

            double term  = 4.0*PI/(120.0*r*lik5*uik5)*(15.0*uik*lik*r*(uik4-lik4) - 10.0*uik2*lik2*(uik3-lik3) + 6.0*(sK2-r2)*(uik5-lik5));
            double dl;
            if (radiusI > (r-sK) || rmax < rmixo) {
                dl  = -5.0*lik2 + 3.0*(r2 + sK2);
                dl /= -lik5;
            } else {
                dl  = 5.0*lik3 - 33.0*lik*r2 - 3.0*lik*sK2 + 15.0*(lik2*r+r3-r*sK2);
                dl /= lik6;
            }
            double du     = 5.0*uik3 - 33.0*uik*r2 - 3.0*uik*sK2 + 15.0*(uik2*r+r3-r*sK2);
                   du    /= -uik6;
            double idisp  = -2.0*ao*term;
                   de    -= 2.0*ao*PI*(dl + du)/(15.0*r2);
                   term   = 4.0*PI/(2640.0*r*lik12*uik12) * (120.0*uik*lik*r*(uik11-lik11) - 66.0*uik2*lik2*(uik10-lik10) + 55.0*(sK2-r2)*(uik12-lik12));
            if (radiusI > (r-sK) || rmax < rmixo) {
                dl  = -6.0*lik2 + 5.0*r2 + 5.0*sK2;
                dl /= -lik12;
            } else {
                dl  = 6.0*lik3 - 125.0*lik*r2 - 5.0*lik*sK2 + 60.0*(lik2*r+r3-r*sK2);
                dl /= lik13;
            }
            du               = 6.0*uik3 - 125.0*uik*r2 - 5.0*uik*sK2 + 60.0*(uik2*r+r3-r*sK2);
            du              /= -uik13;
            de              += ao*rmixo7*PI*(dl + du)/(60.0*r2);
            sum             += ao*rmixo7*term + idisp;
        }
        if (uik > rmixh) {
                   lik   = rmax > rmixh ? rmax : rmixh;
                   lik2  = lik  * lik;
                   lik3  = lik2 * lik;
                   lik4  = lik2 * lik2;
            double lik5  = lik2 * lik3;
            double lik6  = lik3 * lik3;
            double lik10 = lik5 * lik5;
            double lik11 = lik5 * lik6;
            double lik12 = lik6 * lik6;
            double lik13 = lik3 * lik10;

            double term  = 4.0*PI / (120.0*r*lik5*uik5)*(15.0*uik*lik*r*(uik4-lik4) - 10.0*uik2*lik2*(uik3-lik3) + 6.0*(sK2-r2)*(uik5-lik5));
            double dl;
            if (radiusI > (r-sK) || rmax < rmixh) {
                dl  = -5.0*lik2 + 3.0*(r2 + sK2);
                dl /= -lik5;
            } else {
                dl  = 5.0*lik3 - 33.0*lik*r2 - 3.0*lik*sK2 + 15.0*(lik2*r+r3-r*sK2);
                dl /= lik6;
            }
            double du     = 5.0*uik3 - 33.0*uik*r2 - 3.0*uik*sK2 + 15.0*(uik2*r+r3-r*sK2);
                   du    /= -uik6;
            double idisp  = -4.0*ah*term;
            de                = de - 4.0*ah*PI*(dl + du)/(15.0*r2);
            term              = 4.0*PI / (2640.0*r*lik12*uik12)* (120.0*uik*lik*r*(uik11-lik11)- 66.0*uik2*lik2*(uik10-lik10)+ 55.0*(sK2-r2)*(uik12-lik12));
            if (radiusI > (r-sK) || rmax < rmixh) {
                dl = -6.0*lik2 + 5.0*r2 + 5.0*sK2;
                dl = -dl / lik12;
            } else {;
                dl = 6.0*lik3 - 125.0*lik*r2 - 5.0*lik*sK2 + 60.0*(lik2*r+r3-r*sK2);
                dl = dl / lik13;
            }
            du            = 6.0*uik3 - 125.0*uik*r2 - 5.0*uik*sK2 + 60.0*(uik2*r+r3-r*sK2);
            du           /= -uik13;
            double irep   = 2.0*ah*rmixh7*term;
            de           += ah*rmixh7*PI*(dl+du)/(30.0*r2);
            sum          += irep + idisp;
        }
    }

    // increment the individual dispersion gradient components

    de          *= (_slevy*_awater)/r;

    force[0]     = de*xr;
    force[1]     = de*yr;
    force[2]     = de*zr;

    return sum;

}

double AmoebaReferenceWcaDispersionForce::calculateForceAndEnergy(int numParticles,
                                                                  const vector<Vec3>& particlePositions,
                                                                  const std::vector<double>& radii,
                                                                  const std::vector<double>& epsilons,
                                                                  double totalMaximumDispersionEnergy,
                                                                  vector<Vec3>& forces) const {

    // loop over all ixns

    double energy     = 0.0;

    double rmino2     = _rmino*_rmino;
    double rmino3     = rmino2*_rmino;

    double rminh2     = _rminh*_rminh;
    double rminh3     = rminh2*_rminh;

    double intermediateValues[LastIntermediateValueIndex];

    for (unsigned int ii = 0; ii < static_cast<unsigned int>(numParticles); ii++) {
 
        double epsi              = epsilons[ii];
        double rmini             = radii[ii];

        double denominator       = sqrt(_epso) + sqrt(epsi);
        double emixo             = 4.0*_epso*epsi/(denominator*denominator);
        intermediateValues[EMIXO]    = emixo;

        double rminI2            = rmini*rmini;
        double rminI3            = rminI2*rmini;
 
        double rmixo             = 2.0*(rmino3 + rminI3) / (rmino2 + rminI2);
        intermediateValues[RMIXO] = rmixo;

        double rmixo7            = rmixo*rmixo*rmixo;
               rmixo7            = rmixo7*rmixo7*rmixo;
        intermediateValues[RMIXO7] = rmixo7;

        intermediateValues[AO]     = emixo*rmixo7;

                   denominator     = sqrt(_epsh) + sqrt(epsi);

        double emixh              = 4.0*_epsh*epsi/ (denominator*denominator);
        intermediateValues[EMIXH] = emixh;

        double rmixh              = 2.0 * (rminh3 + rminI3) / (rminh2 + rminI2);
        intermediateValues[RMIXH] = rmixh;

        double rmixh7              = rmixh*rmixh*rmixh;
                   rmixh7          = rmixh7*rmixh7*rmixh;
        intermediateValues[RMIXH7] = rmixh7;

        intermediateValues[AH]     = emixh*rmixh7;

        for (unsigned int jj = 0; jj < static_cast<unsigned int>(numParticles); jj++) {

            if (ii == jj)continue;

            Vec3 force;
            energy += calculatePairIxn(rmini, radii[jj],
                                       particlePositions[ii], particlePositions[jj],
                                       intermediateValues, force);
            
            forces[ii][0] += force[0];
            forces[ii][1] += force[1];
            forces[ii][2] += force[2];

            forces[jj][0] -= force[0];
            forces[jj][1] -= force[1];
            forces[jj][2] -= force[2];
        }

    }

    energy = totalMaximumDispersionEnergy - _slevy*_awater*energy;

    return energy;
}
