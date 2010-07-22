/* -------------------------------------------------------------------------- *
 *                                AmoebaOpenMM                                *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
 * Authors:                                                                   *
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

#include "openmm/Force.h"
#include "openmm/OpenMMException.h"
#include "AmoebaWcaDispersionForce.h"
#include "internal/AmoebaWcaDispersionForceImpl.h"
#include <cmath>

using namespace OpenMM;

AmoebaWcaDispersionForce::AmoebaWcaDispersionForce() {
    epso      = 0.1100;
    epsh      = 0.0135;
    rmino     = 1.7025;
    rminh     = 1.3275;
    awater    = 0.033428;
    slevy     = 1.0;
    shctd     = 0.81;
    dispoff   = 0.26;
}

int AmoebaWcaDispersionForce::addParticle( double radius, double epsilon ) {
    parameters.push_back(WcaDispersionInfo( radius, epsilon));
    return parameters.size()-1;
}

void AmoebaWcaDispersionForce::getParticleParameters(int particleIndex, double& radius, double& epsilon ) const {
    radius          = parameters[particleIndex].radius;
    epsilon         = parameters[particleIndex].epsilon;
}

void AmoebaWcaDispersionForce::setParticleParameters(int particleIndex, double radius, double epsilon ) {
    parameters[particleIndex].radius          = radius;
    parameters[particleIndex].epsilon         = epsilon;
}

void AmoebaWcaDispersionForce::getMaximumDispersionEnergy(int particleIndex, double& maxDispersionEnergy ) const {

    const double pi  = 3.1415926535897;

    // from last loop in subroutine knp in ksolv.f

    double rdisp, epsi;
    getParticleParameters( particleIndex, rdisp, epsi );
    if( epsi <= 0.0 || rdisp <= 0.0 ){
        maxDispersionEnergy = 0.0;
        return;
    }
    double rmini     = rdisp;
    rdisp           += getDispoff();
    double epso      = getEpso();
    double emixo     = std::sqrt(epso) + std::sqrt(epsi);
           emixo     = 4.0*epso*epsi/(emixo*emixo);

    double rmino     = getRmino();
    double rmino2    = rmino*rmino;
    double rmini2    = rmini*rmini;
    double rmixo     = 2.0*(rmino2*rmino + rmini2*rmini) / (rmino2 + rmini2);

    double rmixo3    = rmixo*rmixo*rmixo;
    double rmixo7    = rmixo*rmixo3*rmixo3;
    double ao        = emixo*rmixo7;

    double epsh      = getEpsh();
    double emixh     = std::sqrt(epsh) + std::sqrt(epsi);
           emixh     = 4.0*epsh*epsi/(emixh*emixh);
    double rminh     = getRminh();
    double rminh2    = rminh*rminh;
    double rmixh     = rminh*rminh + rmini2;
           rmixh     = 2.0 * (rminh2*rminh + rmini2*rmini) / (rminh2 + rmini2);
    double rmixh3    = rmixh*rmixh*rmixh;
    double rmixh7    = rmixh3*rmixh3*rmixh;
    double ah        = emixh*rmixh7;

    double rdisp3    = rdisp*rdisp*rdisp;
    double rdisp7    = rdisp*rdisp3*rdisp3;
    double rdisp11   = rdisp7*rdisp3*rdisp;

    double cdisp;
    if( rdisp < rmixh) {
        cdisp = -4.0*pi*emixh*(rmixh3-rdisp3)/3.0 - emixh*18.0/11.0*rmixh3*pi;
    } else {
        cdisp = 2.0*pi*(2.0*rmixh7-11.0*rdisp7)*ah/ (11.0*rdisp11);
    }
    cdisp *= 2.0;
    if (rdisp < rmixo ) {
        cdisp -= 4.0*pi*emixo*(rmixo3-rdisp3)/3.0;
        cdisp -= emixo*18.0/11.0*rmixo3*pi;
    } else { 
        cdisp += 2.0*pi*(2.0*rmixo7-11.0*rdisp7) * ao/(11.0*rdisp11);
    }
    maxDispersionEnergy = getSlevy()*getAwater()*cdisp;
//    (void) fprintf( stderr,"Wca %5d %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e\n",
//                    particleIndex, rdisp,rmini,epsi, emixh,rmixh,emixo,rmixo,cdisp );

    return;
}

double AmoebaWcaDispersionForce::getTotalMaximumDispersionEnergy( void ) const {

    double totalMaximumDispersionEnergy = 0.0;
    for( int ii = 0; ii < getNumParticles(); ii++ ){
        double maximumDispersionEnergy;
        getMaximumDispersionEnergy( ii, maximumDispersionEnergy );
        totalMaximumDispersionEnergy += maximumDispersionEnergy;
    }
    return totalMaximumDispersionEnergy;
}

double AmoebaWcaDispersionForce::getEpso( void ) const {
    return epso;
}

double AmoebaWcaDispersionForce::getEpsh( void ) const {
    return epsh;
}

double AmoebaWcaDispersionForce::getRmino( void ) const {
    return rmino;
}

double AmoebaWcaDispersionForce::getRminh( void ) const {
    return rminh;
}

double AmoebaWcaDispersionForce::getAwater( void ) const {
    return awater;
}

double AmoebaWcaDispersionForce::getShctd( void ) const {
    return shctd;
}

double AmoebaWcaDispersionForce::getDispoff( void ) const {
    return dispoff;
}

double AmoebaWcaDispersionForce::getSlevy( void ) const {
    return slevy;
}

void AmoebaWcaDispersionForce::setEpso( double inputEpso ){
    epso = inputEpso;
}

void AmoebaWcaDispersionForce::setEpsh( double inputEpsh ){
    epsh = inputEpsh;
}

void AmoebaWcaDispersionForce::setRmino( double inputRmino ){
    rmino = inputRmino;
}

void AmoebaWcaDispersionForce::setRminh( double inputRminh ){
    rminh = inputRminh;
}

void AmoebaWcaDispersionForce::setAwater( double inputAwater ){
    awater = inputAwater;
}

void AmoebaWcaDispersionForce::setShctd( double inputShctd ){
    shctd = inputShctd;
}

void AmoebaWcaDispersionForce::setDispoff( double inputDispoff ){
    dispoff = inputDispoff;
}

void AmoebaWcaDispersionForce::setSlevy( double inputSlevy ){
    slevy = inputSlevy;
}

ForceImpl* AmoebaWcaDispersionForce::createImpl() {
    return new AmoebaWcaDispersionForceImpl(*this);
}
