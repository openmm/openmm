/* -------------------------------------------------------------------------- *
 *                                OpenMMAmoeba                                *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
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

#include "../../../tests/AssertionUtilities.h"
#include "openmm/AmoebaVdwForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

void testSerialization() {
    // Create a Force.

    AmoebaVdwForce force1;
    force1.setSigmaCombiningRule(   "GEOMETRIC" );
    force1.setEpsilonCombiningRule( "GEOMETRIC" );
    force1.setCutoff( 0.9 );
    force1.setUseNeighborList( 1 );
    force1.setPBC( 1 );

    force1.addParticle(0, 1, 1.0, 2.0, 0.9);
    force1.addParticle(1, 0, 1.1, 2.1, 0.9);
    force1.addParticle(2, 3, 1.3, 4.1, 0.9);
    for( unsigned int ii = 0; ii < 3; ii++ ){
        std::vector< int > exclusions;
        exclusions.push_back( ii );
        exclusions.push_back( ii + 1 );
        exclusions.push_back( ii + 10 );
        force1.setParticleExclusions( ii, exclusions );
    }

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<AmoebaVdwForce>(&force1, "Force", buffer);
#ifdef AMOEBA_DEBUG
    if( 0 ){
        FILE* filePtr = fopen("Vdw.xml", "w" );
        (void) fprintf( filePtr, "%s", buffer.str().c_str() );
        (void) fclose( filePtr );
    }
#endif

    AmoebaVdwForce* copy = XmlSerializer::deserialize<AmoebaVdwForce>(buffer);

    // Compare the two forces to see if they are identical.  
    AmoebaVdwForce& force2 = *copy;

    ASSERT_EQUAL(force1.getSigmaCombiningRule(),    force2.getSigmaCombiningRule());
    ASSERT_EQUAL(force1.getEpsilonCombiningRule(),  force2.getEpsilonCombiningRule());
    ASSERT_EQUAL(force1.getCutoff(),                force2.getCutoff());
    ASSERT_EQUAL(force1.getUseNeighborList(),       force2.getUseNeighborList());
    ASSERT_EQUAL(force1.getPBC(),                   force2.getPBC());

    ASSERT_EQUAL(force1.getNumParticles(),          force2.getNumParticles());

    for (unsigned int ii = 0; ii < force1.getNumParticles(); ii++) {

        int ivIndex1, classIndex1;
        int ivIndex2, classIndex2;

        double sigma1, epsilon1, reductionFactor1;
        double sigma2, epsilon2, reductionFactor2;

        force1.getParticleParameters( ii, ivIndex1, classIndex1, sigma1, epsilon1, reductionFactor1 );
        force2.getParticleParameters( ii, ivIndex2, classIndex2, sigma2, epsilon2, reductionFactor2 );

        ASSERT_EQUAL(ivIndex1,          ivIndex2 );
        ASSERT_EQUAL(classIndex1,       classIndex2 );
        ASSERT_EQUAL(sigma1,            sigma2);
        ASSERT_EQUAL(epsilon1,          epsilon2);
        ASSERT_EQUAL(reductionFactor1,  reductionFactor2);
    }
    for (unsigned int ii = 0; ii < force1.getNumParticles(); ii++) {

        std::vector< int > exclusions1;
        std::vector< int > exclusions2;

        force1.getParticleExclusions( ii, exclusions1 );
        force2.getParticleExclusions( ii, exclusions2 );

        ASSERT_EQUAL(exclusions1.size(), exclusions2.size());
        for (unsigned int jj = 0; jj < exclusions1.size(); jj++) {
            int hit = 0;
            for (unsigned int kk = 0; kk < exclusions2.size(); kk++) {
                if( exclusions2[jj] == exclusions1[kk] )hit++;
            }
            ASSERT_EQUAL(hit, 1);
        }
    }
}

int main() {
    try {
        testSerialization();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

