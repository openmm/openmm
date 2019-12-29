/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2014 Stanford University and the Authors.      *
 * Authors: Andrew C. Simmonett and Andreas Kraemer
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
#include "openmm/NoseHooverIntegrator.h"
#include "openmm/serialization/XmlSerializer.h"
#include "openmm/System.h"
#include "openmm/Context.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

void assertIntegratorsEqual(const NoseHooverIntegrator& integrator1, const NoseHooverIntegrator& integrator2){
    ASSERT_EQUAL(integrator1.getStepSize(), integrator2.getStepSize());
    ASSERT_EQUAL(integrator1.getConstraintTolerance(), integrator2.getConstraintTolerance());
    ASSERT_EQUAL(integrator1.getMaximumPairDistance(), integrator2.getMaximumPairDistance());
    ASSERT_EQUAL(integrator1.getNumThermostats(), integrator2.getNumThermostats());
    for (int i = 0; i < integrator1.getNumThermostats(); i++) {
        const auto &thermostat1 = integrator1.getThermostat(i);
        const auto &thermostat2 = integrator2.getThermostat(i);
        ASSERT_EQUAL(thermostat1.getTemperature(), thermostat2.getTemperature());
        ASSERT_EQUAL(thermostat1.getCollisionFrequency(), thermostat2.getCollisionFrequency());
        ASSERT_EQUAL(thermostat1.getRelativeTemperature(), thermostat2.getRelativeTemperature());
        ASSERT_EQUAL(thermostat1.getRelativeCollisionFrequency(), thermostat2.getRelativeCollisionFrequency());
        ASSERT_EQUAL(thermostat1.getChainLength(), thermostat2.getChainLength());
        ASSERT_EQUAL(thermostat1.getNumMultiTimeSteps(), thermostat2.getNumMultiTimeSteps());
        ASSERT_EQUAL(thermostat1.getNumYoshidaSuzukiTimeSteps(), thermostat2.getNumYoshidaSuzukiTimeSteps());
        ASSERT_EQUAL(thermostat1.getChainID(), thermostat2.getChainID());
        const auto &thermostat1Atoms = thermostat1.getThermostatedAtoms();
        const auto &thermostat2Atoms = thermostat2.getThermostatedAtoms();
        ASSERT_EQUAL(thermostat1Atoms.size(), thermostat2Atoms.size());
        for (int j = 0; j < thermostat1Atoms.size(); ++j) {
            ASSERT_EQUAL(thermostat1Atoms[j], thermostat2Atoms[j]);
        }
        const auto &thermostat1Pairs = thermostat1.getThermostatedPairs();
        const auto &thermostat2Pairs = thermostat2.getThermostatedPairs();
        ASSERT_EQUAL(thermostat1Pairs.size(), thermostat2Pairs.size());
        for (int j = 0; j < thermostat1Pairs.size(); ++j) {
            ASSERT_EQUAL(thermostat1Pairs[j].first, thermostat2Pairs[j].first);
            ASSERT_EQUAL(thermostat1Pairs[j].second, thermostat2Pairs[j].second);
        }
    }
}

void testSerialization() {

    // Check with custom subsystem thermostats

    NoseHooverIntegrator integrator_sub (0.0006);
    integrator_sub.setConstraintTolerance(0.0404);
    integrator_sub.setMaximumPairDistance(0.0051);
    integrator_sub.addSubsystemThermostat(
        {0,1,2,3,4,7}, {{0,7}}, 301.1, 1.1, 1.2, 1.3, 9, 2, 5
    );

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<NoseHooverIntegrator>(&integrator_sub, "Integrator", buffer);
    NoseHooverIntegrator* copy = XmlSerializer::deserialize<NoseHooverIntegrator>(buffer);
    assertIntegratorsEqual(integrator_sub, *copy);

    // Check with default constructor

    System system;
    for (int i=0; i<10; i++) system.addParticle(1.0);
    NoseHooverIntegrator integrator(331, 1.1, 0.004, 5, 5, 5);
    Context context(system, integrator);

    // Serialize and then deserialize it.
    stringstream buffer2;
    XmlSerializer::serialize<NoseHooverIntegrator>(&integrator, "Integrator", buffer2);
    copy = XmlSerializer::deserialize<NoseHooverIntegrator>(buffer2);
    // for thermostats that apply to the whole system, the particles are not serialized ...
    ASSERT_EQUAL(copy->getThermostat(0).getThermostatedAtoms().size(), 0);

    // ... but assigned when creating a context.
    Context context2(system, *copy);
    assertIntegratorsEqual(integrator, *copy);

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

