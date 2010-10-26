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
#include "AmoebaTorsionForce.h"
#include "internal/AmoebaTorsionForceImpl.h"

using namespace OpenMM;

AmoebaTorsionForce::AmoebaTorsionForce() {
}

int AmoebaTorsionForce::addTorsion(int particle1, int particle2, int particle3, int particle4,
                                   std::vector<double>& torsion1, std::vector<double>& torsion2, std::vector<double>& torsion3 ) {
    torsions.push_back(TorsionInfo(particle1, particle2, particle3, particle4, torsion1, torsion2, torsion3));
    return torsions.size()-1;
}

void AmoebaTorsionForce::getTorsionParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4,
                                               std::vector<double>& torsion1, std::vector<double>& torsion2, std::vector<double>& torsion3 ) const {
    particle1       = torsions[index].particle1;
    particle2       = torsions[index].particle2;
    particle3       = torsions[index].particle3;
    particle4       = torsions[index].particle4;

    torsion1.resize( AmoebaTorsionForce::ParametersPerTorsion );
    torsion2.resize( AmoebaTorsionForce::ParametersPerTorsion );
    torsion3.resize( AmoebaTorsionForce::ParametersPerTorsion );

    for( unsigned int ii = 0; ii < AmoebaTorsionForce::ParametersPerTorsion; ii++ ){
       torsion1[ii]  =  torsions[index].torsionParameters[0][ii];
       torsion2[ii]  =  torsions[index].torsionParameters[1][ii];
       torsion3[ii]  =  torsions[index].torsionParameters[2][ii];
    }   
}

void AmoebaTorsionForce::setTorsionParameters(int index, int particle1, int particle2, int particle3, int particle4,
                                              std::vector<double>& torsion1, std::vector<double>& torsion2, std::vector<double>& torsion3 ) {

    torsions[index].particle1  = particle1;
    torsions[index].particle2  = particle2;
    torsions[index].particle3  = particle3;
    torsions[index].particle4  = particle4;

    torsions[index].copyTorsionParameter( 0, torsion1 );
    torsions[index].copyTorsionParameter( 1, torsion2 );
    torsions[index].copyTorsionParameter( 2, torsion3 );
}

ForceImpl* AmoebaTorsionForce::createImpl() {
    return new AmoebaTorsionForceImpl(*this);
}
