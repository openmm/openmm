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
#include "openmm/CustomCentroidBondForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

void testSerialization() {
    // Create a Force.

    CustomCentroidBondForce force(3, "5*sin(distance(g1,g2))^2+y*z");
    force.setForceGroup(3);
    force.addGlobalParameter("x", 1.3);
    force.addGlobalParameter("y", 2.221);
    force.addPerBondParameter("z");
    force.addEnergyParameterDerivative("y");
    for (int i = 0; i < 3; i++) {
        vector<int> particles;
        vector<double> weights;
        for (int j = 0; j < i+1; j++) {
            particles.push_back(i+j);
            if (i < 2)
                weights.push_back(1.0/(i+1));
        }
        force.addGroup(particles, weights);
    }
    vector<int> groups(3);
    vector<double> params(1);
    groups[0] = 0;
    groups[1] = 1;
    groups[2] = 2;
    params[0] = 1.0;
    force.addBond(groups, params);
    groups[0] = 2;
    groups[1] = 3;
    groups[2] = 4;
    params[0] = -3.3;
    force.addBond(groups, params);
    groups[0] = 3;
    groups[1] = 5;
    groups[2] = 1;
    params[0] = 2.1;
    force.addBond(groups, params);
    vector<double> values(10);
    for (int i = 0; i < 10; i++)
        values[i] = sin((double) i);
    force.addTabulatedFunction("f", new Continuous1DFunction(values, 0.5, 1.5));
    force.setUsesPeriodicBoundaryConditions(true);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<CustomCentroidBondForce>(&force, "Force", buffer);
    CustomCentroidBondForce* copy = XmlSerializer::deserialize<CustomCentroidBondForce>(buffer);

    // Compare the two forces to see if they are identical.

    CustomCentroidBondForce& force2 = *copy;
    ASSERT_EQUAL(force.getForceGroup(), force2.getForceGroup());
    ASSERT_EQUAL(force.getNumGroupsPerBond(), force2.getNumGroupsPerBond());
    ASSERT_EQUAL(force.getEnergyFunction(), force2.getEnergyFunction());
    ASSERT_EQUAL(force.getNumPerBondParameters(), force2.getNumPerBondParameters());
    for (int i = 0; i < force.getNumPerBondParameters(); i++)
        ASSERT_EQUAL(force.getPerBondParameterName(i), force2.getPerBondParameterName(i));
    ASSERT_EQUAL(force.getNumGlobalParameters(), force2.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        ASSERT_EQUAL(force.getGlobalParameterName(i), force2.getGlobalParameterName(i));
        ASSERT_EQUAL(force.getGlobalParameterDefaultValue(i), force2.getGlobalParameterDefaultValue(i));
    }
    ASSERT_EQUAL(force.getNumEnergyParameterDerivatives(), force2.getNumEnergyParameterDerivatives());
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++)
        ASSERT_EQUAL(force.getEnergyParameterDerivativeName(i), force2.getEnergyParameterDerivativeName(i));
    ASSERT_EQUAL(force.usesPeriodicBoundaryConditions(), force2.usesPeriodicBoundaryConditions());
    ASSERT_EQUAL(force.getNumGroups(), force2.getNumGroups());
    for (int i = 0; i < force.getNumGroups(); i++) {
        vector<int> particles1, particles2;
        vector<double> weights1, weights2;
        force.getGroupParameters(i, particles1, weights1);
        force2.getGroupParameters(i, particles2, weights2);
        ASSERT_EQUAL(weights1.size(), weights2.size());
        for (int j = 0; j < (int) weights1.size(); j++)
            ASSERT_EQUAL(weights1[j], weights2[j]);
        ASSERT_EQUAL(particles1.size(), particles2.size());
        for (int j = 0; j < (int) particles1.size(); j++)
            ASSERT_EQUAL(particles1[j], particles2[j]);
    }
    ASSERT_EQUAL(force.getNumBonds(), force2.getNumBonds());
    for (int i = 0; i < force.getNumBonds(); i++) {
        vector<int> groups1, groups2;
        vector<double> params1, params2;
        force.getBondParameters(i, groups1, params1);
        force2.getBondParameters(i, groups2, params2);
        ASSERT_EQUAL(params1.size(), params2.size());
        for (int j = 0; j < (int) params1.size(); j++)
            ASSERT_EQUAL(params1[j], params2[j]);
        ASSERT_EQUAL(groups1.size(), groups2.size());
        for (int j = 0; j < (int) groups1.size(); j++)
            ASSERT_EQUAL(groups1[j], groups2[j]);
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
