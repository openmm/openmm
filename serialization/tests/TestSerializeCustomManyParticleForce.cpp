/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2014 Stanford University and the Authors.      *
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
#include "openmm/CustomManyParticleForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

void testSerialization() {
    // Create a Force.

    CustomManyParticleForce force(3, "C*(a1+a2+a3)*(distance(p1,p2)+distance(p1,p3))");
    force.setForceGroup(3);
    force.setNonbondedMethod(CustomManyParticleForce::CutoffPeriodic);
    force.setPermutationMode(CustomManyParticleForce::UniqueCentralParticle);
    force.setCutoffDistance(2.1);
    force.addGlobalParameter("C", 1.3);
    force.addPerParticleParameter("a");
    vector<double> params(1);
    params[0] = 1.0;
    force.addParticle(params, 0);
    params[0] = -3.3;
    force.addParticle(params, 2);
    params[0] = 2.1;
    force.addParticle(params, 3);
    set<int> types1, types2;
    types1.insert(0);
    types2.insert(2);
    types2.insert(3);
    force.setTypeFilter(0, types1);
    force.setTypeFilter(2, types2);
    force.addExclusion(0, 1);
    force.addExclusion(1, 2);
    vector<double> values(10);
    for (int i = 0; i < 10; i++)
        values[i] = sin((double) i);
    force.addTabulatedFunction("f", new Continuous1DFunction(values, 0.5, 1.5));

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<CustomManyParticleForce>(&force, "Force", buffer);
    CustomManyParticleForce* copy = XmlSerializer::deserialize<CustomManyParticleForce>(buffer);

    // Compare the two forces to see if they are identical.

    CustomManyParticleForce& force2 = *copy;
    ASSERT_EQUAL(force.getForceGroup(), force2.getForceGroup());
    ASSERT_EQUAL(force.getNumParticlesPerSet(), force2.getNumParticlesPerSet());
    ASSERT_EQUAL(force.getEnergyFunction(), force2.getEnergyFunction());
    ASSERT_EQUAL(force.getNonbondedMethod(), force2.getNonbondedMethod());
    ASSERT_EQUAL(force.getCutoffDistance(), force2.getCutoffDistance());
    ASSERT_EQUAL(force.getNumPerParticleParameters(), force2.getNumPerParticleParameters());
    for (int i = 0; i < force.getNumPerParticleParameters(); i++)
        ASSERT_EQUAL(force.getPerParticleParameterName(i), force2.getPerParticleParameterName(i));
    ASSERT_EQUAL(force.getNumGlobalParameters(), force2.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        ASSERT_EQUAL(force.getGlobalParameterName(i), force2.getGlobalParameterName(i));
        ASSERT_EQUAL(force.getGlobalParameterDefaultValue(i), force2.getGlobalParameterDefaultValue(i));
    }
    ASSERT_EQUAL(force.getNumParticles(), force2.getNumParticles());
    for (int i = 0; i < force.getNumParticles(); i++) {
        vector<double> params1, params2;
        int type1, type2;
        force.getParticleParameters(i, params1, type1);
        force2.getParticleParameters(i, params2, type2);
        ASSERT_EQUAL(type1, type2);
        ASSERT_EQUAL(params1.size(), params2.size());
        for (int j = 0; j < (int) params1.size(); j++)
            ASSERT_EQUAL(params1[j], params2[j]);
    }
    ASSERT_EQUAL(force.getNumExclusions(), force2.getNumExclusions());
    for (int i = 0; i < force.getNumExclusions(); i++) {
        int a1, a2, b1, b2;
        force.getExclusionParticles(i, a1, b1);
        force2.getExclusionParticles(i, a2, b2);
        ASSERT_EQUAL(a1, a2);
        ASSERT_EQUAL(b1, b2);
    }
    for (int i = 0; i < force.getNumParticlesPerSet(); i++) {
        set<int> set1, set2;
        force.getTypeFilter(i, set1);
        force2.getTypeFilter(i, set2);
        ASSERT_EQUAL_CONTAINERS(set1, set2);
    }
    ASSERT_EQUAL(force.getNumTabulatedFunctions(), force2.getNumTabulatedFunctions());
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        double min1, min2, max1, max2;
        vector<double> val1, val2;
        dynamic_cast<Continuous1DFunction&>(force.getTabulatedFunction(i)).getFunctionParameters(val1, min1, max1);
        dynamic_cast<Continuous1DFunction&>(force2.getTabulatedFunction(i)).getFunctionParameters(val2, min2, max2);
        ASSERT_EQUAL(force.getTabulatedFunctionName(i), force2.getTabulatedFunctionName(i));
        ASSERT_EQUAL(min1, min2);
        ASSERT_EQUAL(max1, max2);
        ASSERT_EQUAL(val1.size(), val2.size());
        for (int j = 0; j < (int) val1.size(); j++)
            ASSERT_EQUAL(val1[j], val2[j]);
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
