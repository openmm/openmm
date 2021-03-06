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
    force1.setForceGroup(3);
    force1.setName("custom name");
    force1.setSigmaCombiningRule("GEOMETRIC");
    force1.setEpsilonCombiningRule("GEOMETRIC");
    force1.setCutoff(0.9);
    force1.setNonbondedMethod(AmoebaVdwForce::CutoffPeriodic);
    force1.setAlchemicalMethod(AmoebaVdwForce::None);
    force1.setPotentialFunction(AmoebaVdwForce::Buffered147);

    force1.addParticle(0, 1.0, 2.0, 0.9, false);
    force1.addParticle(1, 1.1, 2.1, 0.9, true);
    force1.addParticle(2, 1.3, 4.1, 0.9, false);
    for (int i = 0; i < 3; i++) {
        std::vector< int > exclusions;
        exclusions.push_back(i);
        exclusions.push_back(i + 1);
        exclusions.push_back(i + 10);
        force1.setParticleExclusions(i, exclusions);
    }

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<AmoebaVdwForce>(&force1, "Force", buffer);
    AmoebaVdwForce* copy = XmlSerializer::deserialize<AmoebaVdwForce>(buffer);

    // Compare the two forces to see if they are identical.  
    AmoebaVdwForce& force2 = *copy;

    ASSERT_EQUAL(force1.getForceGroup(), force2.getForceGroup());
    ASSERT_EQUAL(force1.getName(), force2.getName());
    ASSERT_EQUAL(force1.getSigmaCombiningRule(),    force2.getSigmaCombiningRule());
    ASSERT_EQUAL(force1.getEpsilonCombiningRule(),  force2.getEpsilonCombiningRule());
    ASSERT_EQUAL(force1.getCutoff(),                force2.getCutoff());
    ASSERT_EQUAL(force1.getNonbondedMethod(),       force2.getNonbondedMethod());
    ASSERT_EQUAL(force1.getAlchemicalMethod(),      force2.getAlchemicalMethod());
    ASSERT_EQUAL(force1.getPotentialFunction(),     force2.getPotentialFunction());

    ASSERT_EQUAL(force1.getNumParticles(),          force2.getNumParticles());

    for (int i = 0; i < force1.getNumParticles(); i++) {

        int ivIndex1, type1;
        int ivIndex2, type2;

        double sigma1, epsilon1, reductionFactor1;
        double sigma2, epsilon2, reductionFactor2;

        bool isAlchemical1;
        bool isAlchemical2;

        force1.getParticleParameters(i, ivIndex1, sigma1, epsilon1, reductionFactor1, isAlchemical1, type1);
        force2.getParticleParameters(i, ivIndex2, sigma2, epsilon2, reductionFactor2, isAlchemical2, type2);

        ASSERT_EQUAL(ivIndex1,          ivIndex2);
        ASSERT_EQUAL(sigma1,            sigma2);
        ASSERT_EQUAL(epsilon1,          epsilon2);
        ASSERT_EQUAL(reductionFactor1,  reductionFactor2);
        ASSERT_EQUAL(isAlchemical1,     isAlchemical2);
        ASSERT_EQUAL(type1,             type2);
    }
    for (int i = 0; i < force1.getNumParticles(); i++) {

        std::vector< int > exclusions1;
        std::vector< int > exclusions2;

        force1.getParticleExclusions(i, exclusions1);
        force2.getParticleExclusions(i, exclusions2);

        ASSERT_EQUAL(exclusions1.size(), exclusions2.size());
        for (int j = 0; j < exclusions1.size(); j++) {
            int hit = 0;
            for (int kk = 0; kk < exclusions2.size(); kk++) {
                if (exclusions2[j] == exclusions1[kk])hit++;
            }
            ASSERT_EQUAL(hit, 1);
        }
    }
}

void testSerializeTypes() {
    // Create a Force that specifies parameters by type.

    AmoebaVdwForce force1;
    force1.setPotentialFunction(AmoebaVdwForce::LennardJones);

    force1.addParticle(0, 2, 1.0, false);
    force1.addParticle(1, 2, 0.9, true);
    force1.addParticle(2, 0, 1.0, false);
    force1.addParticle(3, 1, 0.9, false);
    force1.addParticleType(1.1, 2.0);
    force1.addParticleType(1.2, 2.1);
    force1.addParticleType(1.3, 2.2);
    force1.addTypePair(0, 2, 1.5, 2.5);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<AmoebaVdwForce>(&force1, "Force", buffer);
    AmoebaVdwForce* copy = XmlSerializer::deserialize<AmoebaVdwForce>(buffer);

    // Compare the two forces to see if they are identical.  
    AmoebaVdwForce& force2 = *copy;

    ASSERT_EQUAL(force1.getPotentialFunction(), force2.getPotentialFunction());
    ASSERT_EQUAL(force1.getNumParticles(),      force2.getNumParticles());
    ASSERT_EQUAL(force1.getNumParticleTypes(),  force2.getNumParticleTypes());
    ASSERT_EQUAL(force1.getNumTypePairs(),      force2.getNumTypePairs());

    for (int i = 0; i < force1.getNumParticles(); i++) {
        int ivIndex1, type1;
        int ivIndex2, type2;

        double sigma1, epsilon1, reductionFactor1;
        double sigma2, epsilon2, reductionFactor2;

        bool isAlchemical1;
        bool isAlchemical2;

        force1.getParticleParameters(i, ivIndex1, sigma1, epsilon1, reductionFactor1, isAlchemical1, type1);
        force2.getParticleParameters(i, ivIndex2, sigma2, epsilon2, reductionFactor2, isAlchemical2, type2);

        ASSERT_EQUAL(ivIndex1,          ivIndex2);
        ASSERT_EQUAL(sigma1,            sigma2);
        ASSERT_EQUAL(epsilon1,          epsilon2);
        ASSERT_EQUAL(reductionFactor1,  reductionFactor2);
        ASSERT_EQUAL(isAlchemical1,     isAlchemical2);
        ASSERT_EQUAL(type1,             type2);
    }
    for (int i = 0; i < force1.getNumParticleTypes(); i++) {
        double sigma1, epsilon1;
        double sigma2, epsilon2;
        force1.getParticleTypeParameters(i, sigma1, epsilon1);
        force2.getParticleTypeParameters(i, sigma2, epsilon2);
        ASSERT_EQUAL(sigma1, sigma2);
        ASSERT_EQUAL(epsilon1, epsilon2);
    }
    for (int i = 0; i < force1.getNumTypePairs(); i++) {
        int type11, type21;
        int type12, type22;
        double sigma1, epsilon1;
        double sigma2, epsilon2;
        force1.getTypePairParameters(i, type11, type21, sigma1, epsilon1);
        force2.getTypePairParameters(i, type12, type22, sigma2, epsilon2);
        ASSERT_EQUAL(type11, type12);
        ASSERT_EQUAL(type21, type22);
        ASSERT_EQUAL(sigma1, sigma2);
        ASSERT_EQUAL(epsilon1, epsilon2);
    }
}

int main() {
    try {
        registerAmoebaSerializationProxies();
        testSerialization();
        testSerializeTypes();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

