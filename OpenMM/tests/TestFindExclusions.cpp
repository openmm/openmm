/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
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
#include "openmm/CustomNonbondedForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/OpenMMException.h"
#include <iostream>
#include <set>
#include <vector>

using namespace OpenMM;
using namespace std;

static const int NUM_ATOMS = 20;

/**
 * Add a pair of atoms to the list of exclusions.
 */

void addAtomsToExclusions(int atom1, int atom2, vector<set<int> >& exclusions, int& totalExclusions) {
    if (atom2 < NUM_ATOMS) {
        exclusions[atom1].insert(atom2);
        exclusions[atom2].insert(atom1);
        totalExclusions++;
    }
}

/**
 * This tests the createExceptionsFromBonds() method of NonbondedForce, which identifies pairs of atoms
 * whose nonbonded atoms are either excluded or decreased.  The test system is a chain with branches:
 *
 * 1  3  5  7  9  11  13  15  17  19
 * |  |  |  |  |  |   |   |   |   |
 * 0--2--4--6--8--10--12--14--16--18
 */

void testFindExceptions() {
    NonbondedForce nonbonded;
    vector<pair<int, int> > bonds;
    for (int i = 0; i < NUM_ATOMS; i++)
        nonbonded.addParticle(1.0, 1.0, 2.0);
    // loop over all main-chain atoms (even numbered atoms)
    for (int i = 0; i < NUM_ATOMS-1; i += 2)
    {
        // side-chain bonds
        bonds.push_back(pair<int, int>(i, i+1));
        // main-chain bonds
        if (i < NUM_ATOMS-2) // penultimate atom (NUM_ATOMS-2) has no subsequent main-chain atom
            bonds.push_back(pair<int, int>(i, i+2));
    }
    nonbonded.createExceptionsFromBonds(bonds, 0.2, 0.4);

    // Build lists of the expected exclusions and 1-4s.

    vector<set<int> > expectedExclusions(NUM_ATOMS);
    int totalExclusions = 0;
    for (int i = 0; i < NUM_ATOMS; i += 2) {
        addAtomsToExclusions(i, i+1, expectedExclusions, totalExclusions);
        addAtomsToExclusions(i, i+2, expectedExclusions, totalExclusions);
        addAtomsToExclusions(i, i+3, expectedExclusions, totalExclusions);
        addAtomsToExclusions(i, i+4, expectedExclusions, totalExclusions);
        addAtomsToExclusions(i+1, i+2, expectedExclusions, totalExclusions);
    }
    vector<set<int> > expected14(NUM_ATOMS);
    int total14 = 0;
    for (int i = 0; i < NUM_ATOMS; i += 2) {
        addAtomsToExclusions(i, i+5, expected14, total14);
        addAtomsToExclusions(i, i+6, expected14, total14);
        addAtomsToExclusions(i+1, i+3, expected14, total14);
        addAtomsToExclusions(i+1, i+4, expected14, total14);
    }

    // Compare them to the exceptions that were generated.

    ASSERT_EQUAL(totalExclusions+total14, nonbonded.getNumExceptions());
    for (int i = 0; i < nonbonded.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        nonbonded.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
        if (chargeProd == 0) {
            // This is an exclusion.

            ASSERT_EQUAL(0.0, epsilon);
            ASSERT(expectedExclusions[particle1].find(particle2) != expectedExclusions[particle1].end());
        }
        else {
            // This is a 1-4.

            ASSERT_EQUAL_TOL(0.2, chargeProd, 1e-10);
            ASSERT_EQUAL_TOL(1.0, sigma, 1e-10);
            ASSERT_EQUAL_TOL(0.8, epsilon, 1e-10);
            ASSERT(expected14[particle1].find(particle2) != expected14[particle1].end());
        }
    }
}

/**
 * Test replacing existing exclusions.
 */

void testReplaceExceptions() {
    NonbondedForce nonbonded;
    for (int i = 0; i < 10; i++)
        nonbonded.addParticle(1.0, 1.0, 0.0);
    nonbonded.addException(0, 1, 1.0, 0.0, 0.0);
    nonbonded.addException(0, 2, 2.0, 0.0, 0.0);
    nonbonded.addException(0, 3, 3.0, 0.0, 0.0);
    try {
        nonbonded.addException(0, 3, 4.0, 0.0, 0.0);
        throw std::exception();
    }
    catch (const OpenMMException ex) {
        // This should have thrown an exception.
    }
    int p1, p2;
    double charge, sigma, epsilon;
    nonbonded.getExceptionParameters(2, p1, p2, charge, sigma, epsilon);
    ASSERT(p1 == 0);
    ASSERT(p2 == 3);
    ASSERT(charge == 3.0);
    ASSERT(nonbonded.addException(0, 3, 4.0, 0.0, 0.0, true) == 2);
    nonbonded.getExceptionParameters(2, p1, p2, charge, sigma, epsilon);
    ASSERT(p1 == 0);
    ASSERT(p2 == 3);
    ASSERT(charge == 4.0);
    ASSERT(nonbonded.addException(1, 3, 4.0, 0.0, 0.0) == 3);
    try {
        nonbonded.addException(3, 0, 5.0, 0.0, 0.0);
        throw std::exception();
    }
    catch (const OpenMMException ex) {
        // This should have thrown an exception.
    }
    ASSERT(nonbonded.addException(3, 0, 5.0, 0.0, 0.0, true) == 2);
    nonbonded.getExceptionParameters(2, p1, p2, charge, sigma, epsilon);
    ASSERT(p1 == 3);
    ASSERT(p2 == 0);
    ASSERT(charge == 5.0);
}

/**
 * This is the same as testFindExceptions(), except it tests adding exclusions to a CustomNonbondedForce.
 */

void testFindCustomExclusions() {
    CustomNonbondedForce nonbonded("r");
    vector<pair<int, int> > bonds;
    vector<double> params;
    for (int i = 0; i < NUM_ATOMS; i++)
        nonbonded.addParticle(params);
    // loop over all main-chain atoms (even numbered atoms)
    for (int i = 0; i < NUM_ATOMS-1; i += 2)
    {
        // side-chain bonds
        bonds.push_back(pair<int, int>(i, i+1));
        // main-chain bonds
        if (i < NUM_ATOMS-2) // penultimate atom (NUM_ATOMS-2) has no subsequent main-chain atom
            bonds.push_back(pair<int, int>(i, i+2));
    }
    nonbonded.createExclusionsFromBonds(bonds, 3);

    // Build lists of the expected exclusions.

    vector<set<int> > expectedExclusions(NUM_ATOMS);
    int totalExclusions = 0;
    for (int i = 0; i < NUM_ATOMS; i += 2) {
        addAtomsToExclusions(i, i+1, expectedExclusions, totalExclusions);
        addAtomsToExclusions(i, i+2, expectedExclusions, totalExclusions);
        addAtomsToExclusions(i, i+3, expectedExclusions, totalExclusions);
        addAtomsToExclusions(i, i+4, expectedExclusions, totalExclusions);
        addAtomsToExclusions(i+1, i+2, expectedExclusions, totalExclusions);
        addAtomsToExclusions(i, i+5, expectedExclusions, totalExclusions);
        addAtomsToExclusions(i, i+6, expectedExclusions, totalExclusions);
        addAtomsToExclusions(i+1, i+3, expectedExclusions, totalExclusions);
        addAtomsToExclusions(i+1, i+4, expectedExclusions, totalExclusions);
    }
    for (int i = 0; i < nonbonded.getNumExclusions(); i++) {
        int particle1, particle2;
        nonbonded.getExclusionParticles(i, particle1, particle2);
    }

    // Compare them to the exceptions that were generated.

    ASSERT_EQUAL(totalExclusions, nonbonded.getNumExclusions());
    for (int i = 0; i < nonbonded.getNumExclusions(); i++) {
        int particle1, particle2;
        nonbonded.getExclusionParticles(i, particle1, particle2);
        ASSERT(expectedExclusions[particle1].find(particle2) != expectedExclusions[particle1].end());
    }
}

int main() {
    try {
        testFindExceptions();
        testReplaceExceptions();
        testFindCustomExclusions();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
