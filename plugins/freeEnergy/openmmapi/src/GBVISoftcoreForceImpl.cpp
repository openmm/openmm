/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "openmm/internal/GBVISoftcoreForceImpl.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/freeEnergyKernels.h"
#include "openmm/kernels.h"
#include <vector>
#include <cmath>
#include <cstdio>
#include <sstream>

using namespace OpenMM;
using std::vector;

GBVISoftcoreForceImpl::GBVISoftcoreForceImpl(GBVISoftcoreForce& owner) : owner(owner) {
}

void GBVISoftcoreForceImpl::initialize(ContextImpl& context) {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "GBVISoftcoreForceImpl::initialize";

// ---------------------------------------------------------------------------------------

    kernel                    = context.getPlatform().createKernel(CalcGBVISoftcoreForceKernel::Name(), context);
    System& system            = context.getSystem();
    int     numberOfParticles = owner.getNumParticles();

    // validate number of particles consistent

    if( numberOfParticles != system.getNumParticles() ){
        throw OpenMMException("GBVISoftcoreForce must have exactly as many particles as the System it belongs to.");
    }

    // load 1-2 atom pairs along w/ bond distance using HarmonicBondForce & constraints
    // if numberOfBonds < 1, indicating they were not set by the user

    if( owner.getNumBonds() < 1 ){
        getBondsFromForces( context );
    }   
    int numberOfBonds = owner.getNumBonds();

    std::vector< std::vector<int> > bondIndices;
    bondIndices.resize( numberOfBonds );

    std::vector<double> bondLengths;
    bondLengths.resize( numberOfBonds );
//double maxBond = -1.0;
    for (int i = 0; i < numberOfBonds; i++) {
        int particle1, particle2;
        double bondLength;
        owner.getBondParameters(i, particle1, particle2, bondLength);
        if (particle1 < 0 || particle1 >= owner.getNumParticles()) {
            std::stringstream msg;
            msg << "GBVISoftcoreForce: Illegal particle index: ";
            msg << particle1;
            throw OpenMMException(msg.str());
        }
        if (particle2 < 0 || particle2 >= owner.getNumParticles()) {
            std::stringstream msg;
            msg << "GBVISoftcoreForce: Illegal particle index: ";
            msg << particle2;
            throw OpenMMException(msg.str());
        }
//if( maxBond < bondLength ) maxBond = bondLength;
        if (bondLength < 0.0 ) {
            std::stringstream msg;
            msg << "GBVISoftcoreForce: negative bondlength: ";
            msg << bondLength;
            throw OpenMMException(msg.str());
        }
        bondIndices[i].push_back( particle1 ); 
        bondIndices[i].push_back( particle2 );
        bondLengths[i] = bondLength;
    }

//fprintf( stderr, "%s MaxBond=%14.6e\n", methodName.c_str(), maxBond );

    vector<double> scaledRadii;
    scaledRadii.resize(numberOfParticles);
    findScaledRadii( numberOfParticles, bondIndices, bondLengths, scaledRadii);

    dynamic_cast<CalcGBVISoftcoreForceKernel&>(kernel.getImpl()).initialize(context.getSystem(), owner, scaledRadii);
}

int GBVISoftcoreForceImpl::getBondsFromForces(ContextImpl& context) {

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "GBVISoftcoreForceImpl::getBondsFromForces";

// ---------------------------------------------------------------------------------------

    // load 1-2 atom pairs along w/ bond distance using HarmonicBondForce & constraints
    
    System& system            = context.getSystem();
    for (int i = 0; i < system.getNumForces(); i++) {
        if (dynamic_cast<HarmonicBondForce*>(&system.getForce(i)) != NULL) {
            const HarmonicBondForce& force = dynamic_cast<const HarmonicBondForce&>(system.getForce(i));
            for (int j = 0; j < force.getNumBonds(); ++j) {
                int particle1, particle2;
                double length, k;
                force.getBondParameters(j, particle1, particle2, length, k); 
                owner.addBond( particle1, particle2, length );
            }
            break;
        }
    }

    // Also treat constrained distances as bonds if mass of one particle is < (2.0 + epsilon) (~2=deuterium) 

    for (int j = 0; j < system.getNumConstraints(); j++) {
        int particle1, particle2;
        double distance;
        system.getConstraintParameters(j, particle1, particle2, distance);
        if( system.getParticleMass( particle1 ) < 2.1 || system.getParticleMass( particle2 ) < 2.1 ){
            owner.addBond( particle1, particle2, distance );
        }
    }

    return 0;

}

#define GBVIDebug 0

void GBVISoftcoreForceImpl::findScaledRadii( int numberOfParticles, const std::vector<std::vector<int> >& bondIndices,
                                             const std::vector<double> & bondLengths, std::vector<double> & scaledRadii) const {

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "GBVISoftcoreForceImpl::findScaledRadii";

// ---------------------------------------------------------------------------------------

    // load 1-2 indicies for each atom 

    std::vector<std::vector<int> > bonded12(numberOfParticles);

    for (int i = 0; i < (int) bondIndices.size(); ++i) {
        bonded12[bondIndices[i][0]].push_back(i);
        bonded12[bondIndices[i][1]].push_back(i);
    }

    int errors = 0;

    // compute scaled radii (Eq. 5 of Labute paper [JCC 29 p. 1693-1698 2008])

    for (int j = 0; j < (int) bonded12.size(); ++j){

        double charge;
        double gamma;
        double radiusJ;
        double scaledRadiusJ;
     
        owner.getParticleParameters(j, charge, radiusJ, gamma); 

        if(  bonded12[j].size() == 0 ){
//            (void) fprintf( stderr, "Warning GBVISoftcoreForceImpl::findScaledRadii atom %d has no covalent bonds; using atomic radius=%.3f.\n", j, radiusJ );
            scaledRadiusJ = radiusJ;
//             errors++;
        } else {

            double rJ2    = radiusJ*radiusJ;
    
            // loop over bonded neighbors of atom j, applying Eq. 5 in Labute

            scaledRadiusJ = 0.0;
            for (int i = 0; i < (int) bonded12[j].size(); ++i){
    
               int index            = bonded12[j][i];
               int bondedAtomIndex  = (j == bondIndices[index][0]) ? bondIndices[index][1] : bondIndices[index][0];
              
               double radiusI;
               owner.getParticleParameters(bondedAtomIndex, charge, radiusI, gamma); 
               double rI2           = radiusI*radiusI;
    
               double a_ij          = (radiusI - bondLengths[index]);
                      a_ij         *= a_ij;
                      a_ij          = (rJ2 - a_ij)/(2.0*bondLengths[index]);
    
               double a_ji          = radiusJ - bondLengths[index];
                      a_ji         *= a_ji;
                      a_ji          = (rI2 - a_ji)/(2.0*bondLengths[index]);
    
               scaledRadiusJ       += a_ij*a_ij*(3.0*radiusI - a_ij) + a_ji*a_ji*( 3.0*radiusJ - a_ji );
            }
    
#if GBVIDebug
//(void) fprintf( stderr, "j=%d rJ=%12.4f scaledRadiusJRaw==%12.4f", j, radiusJ, scaledRadiusJ );
#endif
            scaledRadiusJ  = (radiusJ*radiusJ*radiusJ) - 0.125*scaledRadiusJ; 

#if GBVIDebug
//(void) fprintf( stderr, " scRJInt=%12.4f\n", scaledRadiusJ );
#endif
            if( scaledRadiusJ > 0.0 ){
                scaledRadiusJ  = 0.95*pow( scaledRadiusJ, (1.0/3.0) );
            } else {
                scaledRadiusJ  = 0.0;
            }
        }
        scaledRadii[j] = scaledRadiusJ;

    }

    // abort if errors

    if( errors ){
        throw OpenMMException("GBVISoftcoreForceImpl::findScaledRadii errors -- aborting");
    }

#if GBVIDebug
    (void) fprintf( stderr, "                  R              q          gamma   scaled radii no. bnds\n" );
    double totalQ = 0.0;
    for( int i = 0; i < (int) scaledRadii.size(); i++ ){

        double charge;
        double gamma;
        double radiusI;
     
        owner.getParticleParameters(i, charge, radiusI, gamma); 
        totalQ += charge;
        (void) fprintf( stderr, "%4d %14.5e %14.5e %14.5e %14.5e %d\n", i, radiusI, charge, gamma, scaledRadii[i], (int) bonded12[i].size() );
    }
    (void) fprintf( stderr, "Total charge=%e\n", totalQ );
    (void) fflush( stderr );
#endif
#undef GBVIDebug

}

void GBVISoftcoreForceImpl::calcForces(ContextImpl& context) {
    dynamic_cast<CalcGBVISoftcoreForceKernel&>(kernel.getImpl()).executeForces(context);
}

double GBVISoftcoreForceImpl::calcEnergy(ContextImpl& context) {
    return dynamic_cast<CalcGBVISoftcoreForceKernel&>(kernel.getImpl()).executeEnergy(context);
}

std::vector<std::string> GBVISoftcoreForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcGBVISoftcoreForceKernel::Name());
    return names;
}
