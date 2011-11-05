/* -------------------------------------------------------------------------- *
 *                                OpenMMAmoeba                                *
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
#include "openmm/AmoebaVdwForce.h"
#include "openmm/internal/AmoebaVdwForceImpl.h"

using namespace OpenMM;
using std::string;
using std::vector;

AmoebaVdwForce::AmoebaVdwForce() : sigmaCombiningRule("CUBIC-MEAN"), epsilonCombiningRule("HHG"), usePBC(0), cutoff(1.0e+10), useNeighborList(0) {
}

int AmoebaVdwForce::addParticle(int ivIndex, int classIndex, double sigma, double epsilon, double reductionFactor ) {
    parameters.push_back(VdwInfo(ivIndex, classIndex, sigma, epsilon, reductionFactor));
    return parameters.size()-1;
}

void AmoebaVdwForce::getParticleParameters(int particleIndex, int& ivIndex, int& classIndex,
                                           double& sigma, double& epsilon, double& reductionFactor ) const {
    ivIndex         = parameters[particleIndex].ivIndex;
    classIndex      = parameters[particleIndex].classIndex;
    sigma           = parameters[particleIndex].sigma;
    epsilon         = parameters[particleIndex].epsilon;
    reductionFactor = parameters[particleIndex].reductionFactor;
}

void AmoebaVdwForce::setParticleParameters(int particleIndex, int ivIndex, int classIndex,
                                           double sigma, double epsilon, double reductionFactor ) {
    parameters[particleIndex].ivIndex         = ivIndex;
    parameters[particleIndex].classIndex      = classIndex;
    parameters[particleIndex].sigma           = sigma;
    parameters[particleIndex].epsilon         = epsilon;
    parameters[particleIndex].reductionFactor = reductionFactor;
}

void AmoebaVdwForce::setSigmaCombiningRule( const std::string& inputSigmaCombiningRule ) {
    sigmaCombiningRule = inputSigmaCombiningRule;
}

const std::string& AmoebaVdwForce::getSigmaCombiningRule( void ) const {
    return sigmaCombiningRule;
}

void AmoebaVdwForce::setEpsilonCombiningRule( const std::string& inputEpsilonCombiningRule ) {
    epsilonCombiningRule = inputEpsilonCombiningRule;
}

const std::string& AmoebaVdwForce::getEpsilonCombiningRule( void ) const {
    return epsilonCombiningRule;
}

void AmoebaVdwForce::setParticleExclusions( int particleIndex, const std::vector< int >& inputExclusions ) {

   if( exclusions.size() < parameters.size() ){
       exclusions.resize( parameters.size() );
   }
   if(  exclusions.size() < particleIndex ){
       exclusions.resize( particleIndex + 10 );
   }
   for( unsigned int ii = 0; ii < inputExclusions.size(); ii++ ){
       exclusions[particleIndex].push_back( inputExclusions[ii] );
   }
}

void AmoebaVdwForce::getParticleExclusions( int particleIndex, std::vector< int >& outputExclusions ) const {

   if( particleIndex < exclusions.size() ){
       outputExclusions.resize( exclusions[particleIndex].size() );
       for( unsigned int ii = 0; ii < exclusions[particleIndex].size(); ii++ ){
           outputExclusions[ii] = exclusions[particleIndex][ii];
       }
   }

}

void AmoebaVdwForce::setCutoff( double inputCutoff ){
    cutoff = inputCutoff;
}

double AmoebaVdwForce::getCutoff( void ) const {
    return cutoff;
}

void AmoebaVdwForce::setUseNeighborList( int useNeighborListFlag ){
    useNeighborList = useNeighborListFlag;
}

int AmoebaVdwForce::getUseNeighborList( void ) const {
    return useNeighborList;
}

void AmoebaVdwForce::setPBC( int pbcFlag ){
    usePBC = pbcFlag;
}

int AmoebaVdwForce::getPBC( void ) const {
    return usePBC;
}

ForceImpl* AmoebaVdwForce::createImpl() {
    return new AmoebaVdwForceImpl(*this);
}
