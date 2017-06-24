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

#include "openmm/Platform.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/AmoebaVdwForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

extern "C" void registerAmoebaSerializationProxies();

void testSerialization() {
    // Create a Force.

    AmoebaVdwForce force1;
    force1.setSigmaCombiningRule("GEOMETRIC");
    force1.setEpsilonCombiningRule("GEOMETRIC");
    force1.setCutoff(0.9);
    force1.setNonbondedMethod(AmoebaVdwForce::CutoffPeriodic);

    force1.addParticle(0, 34, 1.0, 2.0, 0.9, 1.0);
    force1.addParticle(1, 35, 1.1, 2.1, 0.9, 1.0);
    force1.addParticle(2, 35, 1.3, 4.1, 0.9, 1.0);
    for (unsigned int ii = 0; ii < 3; ii++) {
        std::vector<int> exclusions;
        exclusions.push_back(ii);
        exclusions.push_back(ii + 1);
        exclusions.push_back(ii + 10);
        force1.setParticleExclusions(ii, exclusions);
    }

    force1.computeCombinedSigmaEpsilon();
    force1.addVdwprByOldTypes(34, 35, 1.7, 2.9);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<AmoebaVdwForce>(&force1, "Force", buffer);
    AmoebaVdwForce* copy = XmlSerializer::deserialize<AmoebaVdwForce>(buffer);
    // Compare the two forces to see if they are identical.
    AmoebaVdwForce& force2 = *copy;

    ASSERT_EQUAL(force1.getSigmaCombiningRule(), force2.getSigmaCombiningRule());
    ASSERT_EQUAL(force1.getEpsilonCombiningRule(), force2.getEpsilonCombiningRule());
    ASSERT_EQUAL(force1.getCutoff(), force2.getCutoff());
    ASSERT_EQUAL(force1.getNonbondedMethod(), force2.getNonbondedMethod());
    ASSERT_EQUAL(force1.getNumParticles(), force2.getNumParticles());
    ASSERT_EQUAL(force1.getNumVdwprTypes(), force2.getNumVdwprTypes());

    for (unsigned int ii = 0; ii < static_cast<unsigned int>(force1.getNumParticles()); ii++) {
        int ivIndex1;
        int ivIndex2;
        int vType1;
        int vType2;

        double sigma1, epsilon1, reductionFactor1, lambda1;
        double sigma2, epsilon2, reductionFactor2, lambda2;

        force1.getParticleParameters(ii, ivIndex1, vType1, sigma1, epsilon1, reductionFactor1, lambda1);
        force2.getParticleParameters(ii, ivIndex2, vType2, sigma2, epsilon2, reductionFactor2, lambda2);

        ASSERT_EQUAL(ivIndex1, ivIndex2);
        ASSERT_EQUAL(vType1, vType2);
        ASSERT_EQUAL(sigma1, sigma2);
        ASSERT_EQUAL(epsilon1, epsilon2);
        ASSERT_EQUAL(reductionFactor1, reductionFactor2);
        ASSERT_EQUAL(lambda1, lambda2);
    }
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(force1.getNumParticles()); ii++) {

        std::vector<int> exclusions1;
        std::vector<int> exclusions2;

        force1.getParticleExclusions(ii, exclusions1);
        force2.getParticleExclusions(ii, exclusions2);

        ASSERT_EQUAL(exclusions1.size(), exclusions2.size());
        for (unsigned int jj = 0; jj < exclusions1.size(); jj++) {
            int hit = 0;
            for (unsigned int kk = 0; kk < exclusions2.size(); kk++) {
                if (exclusions2[jj] == exclusions1[kk])
                    hit++;
            }
            ASSERT_EQUAL(hit, 1);
        }
    }
    int numVdwprTypes = force1.getNumVdwprTypes();
    for (int ii = 0; ii < numVdwprTypes; ++ii) {
        for (int jj = 0; jj < numVdwprTypes; ++jj) {
            double combinedSig1, combinedSig2;
            double combinedEps1, combinedEps2;

            force1.getVdwprParameters(ii, jj, combinedSig1, combinedEps1);
            force2.getVdwprParameters(ii, jj, combinedSig2, combinedEps2);

            ASSERT_EQUAL(combinedSig1, combinedSig2);
            ASSERT_EQUAL(combinedEps1, combinedEps2);
        }
    }
}

int main() {
    try {
        registerAmoebaSerializationProxies();
        testSerialization();
    }
    catch (const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
