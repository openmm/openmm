/* -------------------------------------------------------------------------- *
 *                                OpenMMAmoeba                                *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2016 Stanford University and the Authors.      *
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

#include "openmm/Platform.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/AmoebaTorsionTorsionForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>
#include <stdlib.h>

using namespace OpenMM;
using namespace std;

extern "C" void registerAmoebaSerializationProxies();

static void loadTorsionTorsionGrid(std::vector< std::vector< std::vector<double> > >& gridVector) {
 
    static const int gridSize = 25;
    gridVector.resize(gridSize);
    for (unsigned int ii = 0; ii < gridSize; ii++) {
        gridVector[ii].resize(gridSize);
        for (unsigned int jj = 0; jj < gridSize; jj++) {
            gridVector[ii][jj].resize(6);
            for (unsigned int kk = 0; kk < 6; kk++) {
                gridVector[ii][jj][0] = -180.0 + 15.0*static_cast<double>(ii);
                gridVector[ii][jj][1] = -180.0 + 15.0*static_cast<double>(jj);
                gridVector[ii][jj][2] = static_cast<double>(rand());
                gridVector[ii][jj][3] = static_cast<double>(rand());
                gridVector[ii][jj][4] = static_cast<double>(rand());
                gridVector[ii][jj][5] = static_cast<double>(rand());
            }
        }
     }
}

static void compareGrids(const std::vector< std::vector< std::vector<double> > >& grid1, const std::vector< std::vector< std::vector<double> > >& grid2) {

    ASSERT_EQUAL(grid1.size(), grid2.size());
    for (unsigned int ii = 0; ii < grid1.size(); ii++) {
        ASSERT_EQUAL(grid1[ii].size(), grid2[ii].size());
        for (unsigned int jj = 0; jj < grid1[ii].size(); jj++) {
            ASSERT_EQUAL(grid1[ii][jj].size(), grid2[ii][jj].size());
            for (unsigned int kk = 0; kk < grid1[ii][jj].size(); kk++) {
                ASSERT_EQUAL(grid1[ii][jj][kk], grid2[ii][jj][kk]);
            }
        }
    }
}

void testSerialization() {
    // Create a Force.

    AmoebaTorsionTorsionForce force1;

    force1.setForceGroup(3);
    for (unsigned int ii = 0; ii < 5; ii++) {
        std::vector< std::vector< std::vector<double> > > gridVector;
        loadTorsionTorsionGrid(gridVector);
        force1.setTorsionTorsionGrid(ii, gridVector);
    }
    for (unsigned int ii = 0; ii < 5; ii++) {
        force1.addTorsionTorsion(ii, ii+1,ii+3, ii+4, ii+5, ((ii % 2) ? 1 : 0), (ii % 4));
    }
    force1.setUsesPeriodicBoundaryConditions(true);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<AmoebaTorsionTorsionForce>(&force1, "Force", buffer);
    AmoebaTorsionTorsionForce* copy = XmlSerializer::deserialize<AmoebaTorsionTorsionForce>(buffer);

    // Compare the two force1s to see if they are identical.  

    AmoebaTorsionTorsionForce & force2 = *copy;
    ASSERT_EQUAL(force1.getForceGroup(), force2.getForceGroup());
    ASSERT_EQUAL(force1.usesPeriodicBoundaryConditions(), force2.usesPeriodicBoundaryConditions());
    ASSERT_EQUAL(force1.getNumTorsionTorsions(), force2.getNumTorsionTorsions());
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(force1.getNumTorsionTorsions()); ii++) {

        int a1, a2, a3, a4, a5, aChiral, aGridIndex, b1, b2, b3, b4, b5, bChiral, bGridIndex;

        force1.getTorsionTorsionParameters(ii, a1, a2, a3, a4, a5, aChiral, aGridIndex);
        force2.getTorsionTorsionParameters(ii, b1, b2, b3, b4, b5, bChiral, bGridIndex);

        ASSERT_EQUAL(a1, b1);
        ASSERT_EQUAL(a2, b2);
        ASSERT_EQUAL(a3, b3);
        ASSERT_EQUAL(a4, b4);
        ASSERT_EQUAL(a5, b5);
        ASSERT_EQUAL(aChiral, bChiral);
        ASSERT_EQUAL(aGridIndex, bGridIndex);
    }

    ASSERT_EQUAL(force1.getNumTorsionTorsionGrids(), force2.getNumTorsionTorsionGrids());
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(force1.getNumTorsionTorsionGrids()); ii++) {
        const std::vector< std::vector< std::vector<double> > >& grid1 = force1.getTorsionTorsionGrid(ii);
        const std::vector< std::vector< std::vector<double> > >& grid2 = force2.getTorsionTorsionGrid(ii);
        compareGrids(grid1, grid2);
    }

}

int main() {
    try {
        registerAmoebaSerializationProxies();
        testSerialization();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

