
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

#include "AmoebaReferenceForce.h"
#include "AmoebaReferenceWcaDispersionForce.h"

using std::vector;
using OpenMM::RealVec;

AmoebaReferenceWcaDispersionForce::AmoebaReferenceWcaDispersionForce( RealOpenMM epso, RealOpenMM epsh, RealOpenMM rmino, RealOpenMM rminh, 
                                                                      RealOpenMM awater, RealOpenMM shctd, RealOpenMM dispoff, RealOpenMM slevy ) :
                               _epso(epso), _epsh(epsh), _rmino(rmino), _rminh(rminh), _awater(awater), _shctd(shctd), _dispoff(dispoff), _slevy(slevy) {
}   


RealOpenMM AmoebaReferenceWcaDispersionForce::calculatePairIxn( RealOpenMM radiusI, RealOpenMM radiusK, 
                                                                const RealVec& particleIPosition,
                                                                const RealVec& particleJPosition,
                                                                const RealOpenMM* const intermediateValues,
                                                                Vec3& force ) const {

   // ---------------------------------------------------------------------------------------

    static const RealOpenMM zero          =  0.0;
    static const RealOpenMM one           =  1.0;
    static const RealOpenMM two           =  2.0;
    static const RealOpenMM three         =  3.0;
    static const RealOpenMM four          =  4.0;
    static const RealOpenMM five          =  5.0;
    static const RealOpenMM six           =  6.0;
    static const RealOpenMM seven         =  7.0;
    static const RealOpenMM eight         =  8.0;
    static const RealOpenMM ten           = 10.0;
    static const RealOpenMM fortyEight    = 48.0;
    static const RealOpenMM PI            = 3.1415926535897932384;

   // ---------------------------------------------------------------------------------------

    RealOpenMM xr           = particleIPosition[0] - particleJPosition[0];
    RealOpenMM yr           = particleIPosition[1] - particleJPosition[1];
    RealOpenMM zr           = particleIPosition[2] - particleJPosition[2];

    RealOpenMM r2           = xr*xr + yr*yr + zr*zr;
    RealOpenMM r            = SQRT( r2 );
    RealOpenMM r3           = r2*r;

    RealOpenMM sK           = radiusK*_shctd;
    RealOpenMM sK2          = sK*sK;

    RealOpenMM rmixo        = intermediateValues[RMIXO];
    RealOpenMM rmixo7       = intermediateValues[RMIXO7];

    RealOpenMM emixo        = intermediateValues[EMIXO];

    RealOpenMM rmixh        = intermediateValues[RMIXH];
    RealOpenMM rmixh7       = intermediateValues[RMIXH7];

    RealOpenMM emixh        = intermediateValues[EMIXH];

    RealOpenMM ao           = intermediateValues[AO];
    RealOpenMM ah           = intermediateValues[AH];

    RealOpenMM sum          = zero;
    RealOpenMM de           = zero;

    if( radiusI < (r + sK) ){

        RealOpenMM rmax     = (radiusI > (r - sK)) ? radiusI : (r - sK);

        RealOpenMM lik      = rmax;
        RealOpenMM lik2     = lik*lik;
        RealOpenMM lik3     = lik2*lik;
        RealOpenMM lik4     = lik2*lik2;

        if( lik < rmixo ){ 

            RealOpenMM uik  = (r + sK) < rmixo ? (r + sK) : rmixo;
            RealOpenMM uik2 = uik*uik;
            RealOpenMM uik3 = uik2*uik;
            RealOpenMM uik4 = uik2*uik2;

            RealOpenMM term = four*PI/(fortyEight*r)* (three*(lik4-uik4) - eight*r*(lik3-uik3) + six*(r2-sK2)*(lik2-uik2));

            RealOpenMM dl;
            if( radiusI >  (r - sK) ){
                dl  = -lik2 + two*(r2 + sK2);
                dl *= lik2;
            } else {
                dl  = -lik3 + four*lik2*r - six*lik*r2 + two*lik*sK2 + four*r*(r2 - sK2);
                dl *= lik;
            }

            RealOpenMM du;
            if( (r+sK) > rmixo ){
                du  = -uik2 + two*(r2 + sK2);
                du *= -uik2;
            } else {
                du  = -uik3 + four*uik2*r - six*uik*r2 + two*uik*sK2 + four*r*(r2 - sK2);
                du *= -uik;
            }
            de     = -emixo*PI*(dl+du)/(four*r2);
            sum   += -emixo*term;
        }

        if( lik < rmixh ){

            RealOpenMM uik  = (r + sK) < rmixh ? (r + sK) : rmixh;
            RealOpenMM uik2 = uik*uik;
            RealOpenMM uik3 = uik2*uik;
            RealOpenMM uik4 = uik2*uik2;
            RealOpenMM term = four*PI / (fortyEight*r)*(three*(lik4-uik4) - eight*r*(lik3-uik3) + six*(r2-sK2)*(lik2-uik2));
            RealOpenMM dl;
            if( radiusI > (r-sK) ){
                dl  = -lik2 + two*(r2 + sK2);
                dl *= lik2;
            } else {
                dl  = -lik3 + four*lik2*r - six*lik*r2 + two*lik*sK2 + four*r*(r2 -sK2);
                dl *= lik;
            }
      
            RealOpenMM du;
            if (r+sK > rmixh){
                du  = -uik2 + two*(r2 + sK2);
                du *= -uik2;
            } else {
                du  = -uik3 +four*uik2*r - six*uik*r2 + two*uik*sK2 +four*r*(r2 - sK2);
                du *= -uik;
            }
            de   -= two*emixh*PI*(dl+du)/(four*r2);
            sum  -= two*emixh*term;
        }

        RealOpenMM uik   = r + sK;
        RealOpenMM uik2  = uik   * uik;
        RealOpenMM uik3  = uik2  * uik;
        RealOpenMM uik4  = uik2  * uik2;
        RealOpenMM uik5  = uik3  * uik2;
        RealOpenMM uik6  = uik3  * uik3;
        RealOpenMM uik10 = uik5  * uik5;
        RealOpenMM uik11 = uik5  * uik6;
        RealOpenMM uik12 = uik6  * uik6;
        RealOpenMM uik13 = uik10 * uik3;

        if( uik > rmixo ){

            RealOpenMM lik   = rmax > rmixo ? rmax : rmixo;
            RealOpenMM lik2  = lik   * lik;
            RealOpenMM lik3  = lik2  * lik;
            RealOpenMM lik4  = lik2  * lik2;
            RealOpenMM lik5  = lik2  * lik3;
            RealOpenMM lik6  = lik3  * lik3;
            RealOpenMM lik10 = lik5  * lik5;
            RealOpenMM lik11 = lik5  * lik6;
            RealOpenMM lik12 = lik6  * lik6;
            RealOpenMM lik13 = lik10 * lik3;

            RealOpenMM term  = four*PI/(120.0*r*lik5*uik5)*(15.0*uik*lik*r*(uik4-lik4) - ten*uik2*lik2*(uik3-lik3) + six*(sK2-r2)*(uik5-lik5));
            RealOpenMM dl;
            if( radiusI > (r-sK) || rmax < rmixo ){
                dl  = -five*lik2 + three*(r2 + sK2);
                dl /= -lik5;
            } else {
                dl  = five*lik3 - 33.0*lik*r2 - three*lik*sK2 + 15.0*(lik2*r+r3-r*sK2);
                dl /= lik6;
            }
            RealOpenMM du     = five*uik3 - 33.0*uik*r2 - three*uik*sK2 + 15.0*(uik2*r+r3-r*sK2);
                       du    /= -uik6;
            RealOpenMM idisp  = -two*ao*term;
                       de    -= two*ao*PI*(dl + du)/(15.0*r2);
                       term   = four*PI/(2640.0*r*lik12*uik12) * (120.0*uik*lik*r*(uik11-lik11) - 66.0*uik2*lik2*(uik10-lik10) + 55.0*(sK2-r2)*(uik12-lik12));
            if( radiusI > (r-sK) || rmax < rmixo ){
                dl  = -six*lik2 + five*r2 + five*sK2;
                dl /= -lik12;
            } else {
                dl  = six*lik3 - 125.0*lik*r2 - five*lik*sK2 + 60.0*(lik2*r+r3-r*sK2);
                dl /= lik13;
            }
            du               = six*uik3 - 125.0*uik*r2 - five*uik*sK2 + 60.0*(uik2*r+r3-r*sK2);
            du              /= -uik13;
            de              += ao*rmixo7*PI*(dl + du)/(60.0*r2);
            sum             += ao*rmixo7*term + idisp;
        }
        if (uik > rmixh){
                       lik   = rmax > rmixh ? rmax : rmixh;
                       lik2  = lik  * lik;
                       lik3  = lik2 * lik;
                       lik4  = lik2 * lik2;
            RealOpenMM lik5  = lik2 * lik3;
            RealOpenMM lik6  = lik3 * lik3;
            RealOpenMM lik10 = lik5 * lik5;
            RealOpenMM lik11 = lik5 * lik6;
            RealOpenMM lik12 = lik6 * lik6;
            RealOpenMM lik13 = lik3 * lik10;

            RealOpenMM term  = four*PI / (120.0*r*lik5*uik5)*(15.0*uik*lik*r*(uik4-lik4) - ten*uik2*lik2*(uik3-lik3) + six*(sK2-r2)*(uik5-lik5));
            RealOpenMM dl;
            if( radiusI > (r-sK) || rmax < rmixh ){
                dl  = -five*lik2 + three*(r2 + sK2);
                dl /= -lik5;
            } else {
                dl  = five*lik3 - 33.0*lik*r2 - three*lik*sK2 + 15.0*(lik2*r+r3-r*sK2);
                dl /= lik6;
            }
            RealOpenMM du     = five*uik3 - 33.0*uik*r2 - 3.0*uik*sK2 + 15.0*(uik2*r+r3-r*sK2);
                       du    /= -uik6;
            RealOpenMM idisp  = -four*ah*term;
            de                = de - four*ah*PI*(dl + du)/(15.0*r2);
            term              = four*PI / (2640.0*r*lik12*uik12)* (120.0*uik*lik*r*(uik11-lik11)- 66.0*uik2*lik2*(uik10-lik10)+ 55.0*(sK2-r2)*(uik12-lik12));
            if( radiusI > (r-sK) || rmax < rmixh){
                dl = -six*lik2 + five*r2 + five*sK2;
                dl = -dl / lik12;
            } else {;
                dl = six*lik3 - 125.0*lik*r2 - five*lik*sK2 + 60.0*(lik2*r+r3-r*sK2);
                dl = dl / lik13;
            }
            du                = six*uik3 - 125.0*uik*r2 - five*uik*sK2 + 60.0*(uik2*r+r3-r*sK2);
            du               /= -uik13;
            RealOpenMM irep   = two*ah*rmixh7*term;
            de               += ah*rmixh7*PI*(dl+du)/(30.0*r2);
            sum              += irep + idisp;
        }
    }

    // increment the individual dispersion gradient components

    de          *= (_slevy*_awater)/r;

    force[0]     = de*xr;
    force[1]     = de*yr;
    force[2]     = de*zr;

    return sum;

}

RealOpenMM AmoebaReferenceWcaDispersionForce::calculateForceAndEnergy( int numParticles,
                                                                       const vector<RealVec>& particlePositions,
                                                                       const std::vector<RealOpenMM>& radii,
                                                                       const std::vector<RealOpenMM>& epsilons,
                                                                       RealOpenMM totalMaximumDispersionEnergy,
                                                                       vector<RealVec>& forces ) const {

    // ---------------------------------------------------------------------------------------

    //static const std::string methodName = "AmoebaReferenceWcaDispersionForce::calculateForceAndEnergy";

    static const RealOpenMM zero          = 0.0;
    static const RealOpenMM one           = 1.0;
    static const RealOpenMM two           = 2.0;
    static const RealOpenMM four          = 4.0;

    // ---------------------------------------------------------------------------------------

    // loop over all ixns

    RealOpenMM energy     = zero;

    RealOpenMM rmino2     = _rmino*_rmino;
    RealOpenMM rmino3     = rmino2*_rmino;

    RealOpenMM rminh2     = _rminh*_rminh;
    RealOpenMM rminh3     = rminh2*_rminh;

    RealOpenMM intermediateValues[LastIntermediateValueIndex];

    for( unsigned int ii = 0; ii < static_cast<unsigned int>(numParticles); ii++ ){
 
        RealOpenMM epsi              = epsilons[ii];
        RealOpenMM rmini             = radii[ii];

        RealOpenMM denominator       = SQRT(_epso) + SQRT(epsi);
        RealOpenMM emixo             = four*_epso*epsi/(denominator*denominator);
        intermediateValues[EMIXO]    = emixo;

        RealOpenMM rminI2            = rmini*rmini;
        RealOpenMM rminI3            = rminI2*rmini;
 
        RealOpenMM rmixo             = two*(rmino3 + rminI3 ) / (rmino2 + rminI2);
        intermediateValues[RMIXO]    = rmixo;

        RealOpenMM rmixo7            = rmixo*rmixo*rmixo;
                   rmixo7            = rmixo7*rmixo7*rmixo;
        intermediateValues[RMIXO7]   = rmixo7;

        intermediateValues[AO]       = emixo*rmixo7;

                   denominator       = SQRT(_epsh) + SQRT(epsi);

        RealOpenMM emixh             = four*_epsh*epsi/ (denominator*denominator);
        intermediateValues[EMIXH]    = emixh;

        RealOpenMM rmixh             = two * (rminh3 + rminI3) / (rminh2 + rminI2);
        intermediateValues[RMIXH]    = rmixh;

        RealOpenMM rmixh7            = rmixh*rmixh*rmixh;
                   rmixh7            = rmixh7*rmixh7*rmixh;
        intermediateValues[RMIXH7]   = rmixh7;

        intermediateValues[AH]       = emixh*rmixh7;

        for( unsigned int jj = 0; jj < static_cast<unsigned int>(numParticles); jj++ ){

            if( ii == jj )continue;

            Vec3 force;
            energy += calculatePairIxn( rmini, radii[jj],
                                        particlePositions[ii], particlePositions[jj],
                                        intermediateValues, force );
            
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
