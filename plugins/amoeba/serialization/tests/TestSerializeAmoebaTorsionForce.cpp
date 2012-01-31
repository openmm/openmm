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

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/AmoebaTorsionForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

void loadTorsion( std::vector<double>& torsion, double offset ){
    torsion.push_back( offset + 1.0 );
    torsion.push_back( offset + 2.0 );
//    torsion.push_back( offset + 3.0 );
//    torsion.push_back( offset + 4.0 );
}

void compareTorsion( std::vector<double> torsionA, std::vector<double> torsionB ){
    ASSERT_EQUAL(torsionA.size(), torsionB.size());
    for (unsigned int ii = 0; ii < torsionA.size(); ii++) {
        ASSERT_EQUAL(torsionA[ii], torsionB[ii]);
    }
}

void testSerialization() {

    // Create a Force.

    AmoebaTorsionForce force1;
    for( unsigned int ii = 0; ii < 5; ii++ ){
        std::vector<double> torsion1;
        std::vector<double> torsion2;
        std::vector<double> torsion3;
        loadTorsion( torsion1, static_cast<double>(5*ii) + 11.1);
        loadTorsion( torsion2, static_cast<double>(5*ii) + 21.2);
        loadTorsion( torsion3, static_cast<double>(5*ii) + 31.3);
        force1.addTorsion( ii, ii+1,ii+3, ii+4, torsion1, torsion2, torsion3 );
    }

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<AmoebaTorsionForce>(&force1, "Force", buffer);

#ifdef AMOEBA_DEBUG
    if( 0 ){
        FILE* filePtr = fopen("Torsion.xml", "w" );
        (void) fprintf( filePtr, "%s", buffer.str().c_str() );
        (void) fclose( filePtr );
    }
#endif

    AmoebaTorsionForce* copy = XmlSerializer::deserialize<AmoebaTorsionForce>(buffer);

    // Compare the two forces to see if they are identical.  
    AmoebaTorsionForce& force2 = *copy;
    ASSERT_EQUAL(force1.getNumTorsions(), force2.getNumTorsions());
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(force1.getNumTorsions()); ii++) {

        int a1, a2, a3, a4, b1, b2, b3, b4;

        std::vector<double> torsion1a;
        std::vector<double> torsion2a;
        std::vector<double> torsion3a;
        std::vector<double> torsion1b;
        std::vector<double> torsion2b;
        std::vector<double> torsion3b;

        force1.getTorsionParameters( ii, a1, a2, a3, a4, torsion1a, torsion2a, torsion3a);
        force2.getTorsionParameters( ii, b1, b2, b3, b4, torsion1b, torsion2b, torsion3b);

        ASSERT_EQUAL(a1, b1);
        ASSERT_EQUAL(a2, b2);
        ASSERT_EQUAL(a3, b3);
        ASSERT_EQUAL(a4, b4);

        compareTorsion( torsion1a, torsion1b );
        compareTorsion( torsion2a, torsion2b );
        compareTorsion( torsion3a, torsion3b );
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

