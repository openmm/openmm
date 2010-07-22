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
#include "AmoebaVdwForce.h"
#include "internal/AmoebaVdwForceImpl.h"

using namespace OpenMM;

AmoebaVdwForce::AmoebaVdwForce() {
}

int AmoebaVdwForce::addParticle(int ivIndex, int classIndex,  double sigma, double sigma4, double epsilon, double epsilon4, double reductionFactor ) {
    parameters.push_back(VdwInfo(ivIndex, classIndex,  sigma, sigma4, epsilon, epsilon4, reductionFactor));
    return parameters.size()-1;
}

void AmoebaVdwForce::getParticleParameters(int particleIndex, int& ivIndex, int& classIndex,
                                           double& sigma, double& sigma4, double& epsilon, double& epsilon4,  double& reductionFactor ) const {
    ivIndex         = parameters[particleIndex].ivIndex;
    classIndex      = parameters[particleIndex].classIndex;
    sigma           = parameters[particleIndex].sigma;
    sigma4          = parameters[particleIndex].sigma4;
    epsilon         = parameters[particleIndex].epsilon;
    epsilon4        = parameters[particleIndex].epsilon4;
    reductionFactor = parameters[particleIndex].reductionFactor;
}

void AmoebaVdwForce::setParticleParameters(int particleIndex, int ivIndex, int classIndex,
                                           double sigma, double sigma4, double epsilon, double epsilon4, double reductionFactor ) {
    parameters[particleIndex].ivIndex         = ivIndex;
    parameters[particleIndex].classIndex      = classIndex;
    parameters[particleIndex].sigma           = sigma;
    parameters[particleIndex].sigma4          = sigma4;
    parameters[particleIndex].epsilon         = epsilon;
    parameters[particleIndex].epsilon4        = epsilon4;
    parameters[particleIndex].reductionFactor = reductionFactor;
}

void AmoebaVdwForce::setSigEpsTableSize(int tableSize ) {
    sigEpsTable.resize( tableSize );
    for( unsigned int ii = 0; ii < tableSize; ii++ ){
        sigEpsTable[ii].resize( tableSize );
    }
}

int AmoebaVdwForce::getSigEpsTableSize(void ) const {
    return static_cast<int>(sigEpsTable.size( ));
}

void AmoebaVdwForce::setSigEpsTableEntry(int indexI, int indexJ, double combinedSigma, double combinedEpsilon, double combinedSigma4, double combinedEpsilon4 ) {
    sigEpsTable[indexI][indexJ].resize( 4 );
    sigEpsTable[indexI][indexJ][0] = combinedSigma;
    sigEpsTable[indexI][indexJ][1] = combinedEpsilon;
    sigEpsTable[indexI][indexJ][2] = combinedSigma4;
    sigEpsTable[indexI][indexJ][3] = combinedEpsilon4;
}

void AmoebaVdwForce::getSigEpsTableEntry(int indexI, int indexJ, double& combinedSigma, double& combinedEpsilon, double& combinedSigma4, double& combinedEpsilon4 ) const {
    combinedSigma        = sigEpsTable[indexI][indexJ][0];
    combinedEpsilon      = sigEpsTable[indexI][indexJ][1];
    combinedSigma4       = sigEpsTable[indexI][indexJ][2];
    combinedEpsilon4     = sigEpsTable[indexI][indexJ][3];
}

void AmoebaVdwForce::setSigmaCombiningRule(std::string& inputSigmaCombiningRule ) {
    sigmaCombiningRule = inputSigmaCombiningRule;
}

std::string AmoebaVdwForce::getSigmaCombiningRule( void ) const {
    return sigmaCombiningRule;
}

void AmoebaVdwForce::setEpsilonCombiningRule(std::string& inputEpsilonCombiningRule ) {
    epsilonCombiningRule = inputEpsilonCombiningRule;
}

std::string AmoebaVdwForce::getEpsilonCombiningRule( void ) const {
    return epsilonCombiningRule;
}

void AmoebaVdwForce::setParticleExclusions( int particleIndex, std::vector< int >& inputExclusions ) {

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

   outputExclusions.resize( exclusions[particleIndex].size() );
   for( unsigned int ii = 0; ii < exclusions[particleIndex].size(); ii++ ){
       outputExclusions[ii] = exclusions[particleIndex][ii];
   }

}


ForceImpl* AmoebaVdwForce::createImpl() {
    return new AmoebaVdwForceImpl(*this);
}
