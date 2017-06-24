/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
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

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/CustomGBForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

void testSerialization() {
    // Create a Force.

    CustomGBForce force;
    force.setForceGroup(3);
    force.setNonbondedMethod(CustomGBForce::CutoffPeriodic);
    force.setCutoffDistance(2.1);
    force.addGlobalParameter("x", 1.3);
    force.addGlobalParameter("y", 2.221);
    force.addPerParticleParameter("z");
    force.addEnergyParameterDerivative("y");
    force.addComputedValue("a", "x+1", CustomGBForce::ParticlePairNoExclusions);
    force.addComputedValue("b", "y-1", CustomGBForce::SingleParticle);
    force.addEnergyTerm("a*b", CustomGBForce::SingleParticle);
    force.addEnergyTerm("a*b+1", CustomGBForce::ParticlePair);
    vector<double> params(1);
    params[0] = 1.0;
    force.addParticle(params);
    params[0] = -3.3;
    force.addParticle(params);
    params[0] = 2.1;
    force.addParticle(params);
    force.addExclusion(0, 1);
    force.addExclusion(1, 2);
    vector<double> values(10);
    for (int i = 0; i < 10; i++)
        values[i] = sin((double) i);
    force.addTabulatedFunction("f", new Discrete1DFunction(values));

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<CustomGBForce>(&force, "Force", buffer);
    CustomGBForce* copy = XmlSerializer::deserialize<CustomGBForce>(buffer);

    // Compare the two forces to see if they are identical.

    CustomGBForce& force2 = *copy;
    ASSERT_EQUAL(force.getForceGroup(), force2.getForceGroup());
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
    ASSERT_EQUAL(force.getNumEnergyParameterDerivatives(), force2.getNumEnergyParameterDerivatives());
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++)
        ASSERT_EQUAL(force.getEnergyParameterDerivativeName(i), force2.getEnergyParameterDerivativeName(i));
    ASSERT_EQUAL(force.getNumComputedValues(), force2.getNumComputedValues());
    for (int i = 0; i < force.getNumComputedValues(); i++) {
        string name1, name2, expression1, expression2;
        CustomGBForce::ComputationType type1, type2;
        force.getComputedValueParameters(i, name1, expression1, type1);
        force2.getComputedValueParameters(i, name2, expression2, type2);
        ASSERT_EQUAL(name1, name2);
        ASSERT_EQUAL(expression1, expression2);
        ASSERT_EQUAL(type1, type2);
    }
    ASSERT_EQUAL(force.getNumEnergyTerms(), force2.getNumEnergyTerms());
    for (int i = 0; i < force.getNumEnergyTerms(); i++) {
        string expression1, expression2;
        CustomGBForce::ComputationType type1, type2;
        force.getEnergyTermParameters(i, expression1, type1);
        force2.getEnergyTermParameters(i, expression2, type2);
        ASSERT_EQUAL(expression1, expression2);
        ASSERT_EQUAL(type1, type2);
    }
    ASSERT_EQUAL(force.getNumParticles(), force2.getNumParticles());
    for (int i = 0; i < force.getNumParticles(); i++) {
        vector<double> params1, params2;
        force.getParticleParameters(i, params1);
        force2.getParticleParameters(i, params2);
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
    ASSERT_EQUAL(force.getNumTabulatedFunctions(), force2.getNumTabulatedFunctions());
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        vector<double> val1, val2;
        dynamic_cast<Discrete1DFunction&>(force.getTabulatedFunction(i)).getFunctionParameters(val1);
        dynamic_cast<Discrete1DFunction&>(force2.getTabulatedFunction(i)).getFunctionParameters(val2);
        ASSERT_EQUAL(force.getTabulatedFunctionName(i), force2.getTabulatedFunctionName(i));
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
