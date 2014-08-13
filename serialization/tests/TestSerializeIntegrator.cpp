/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Yutong Zhao                                        *
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

#include "openmm/BrownianIntegrator.h"
#include "openmm/CustomIntegrator.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/VariableLangevinIntegrator.h"
#include "openmm/VariableVerletIntegrator.h"
#include "openmm/VerletIntegrator.h"

#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>
#include <stdlib.h>

#include <fstream>

using namespace OpenMM;
using namespace std;

void testSerializeVerletIntegrator() {
    VerletIntegrator *intg = new VerletIntegrator(0.00342);
    stringstream ss;
    XmlSerializer::serialize<Integrator>(intg, "VerletIntegrator", ss);
    VerletIntegrator *intg2 = dynamic_cast<VerletIntegrator*>(XmlSerializer::deserialize<Integrator>(ss));
    ASSERT_EQUAL(intg->getConstraintTolerance(), intg2->getConstraintTolerance());
    ASSERT_EQUAL(intg->getStepSize(), intg2->getStepSize());
    delete intg;
    delete intg2;
}

void testSerializeLangevinIntegrator() {
    LangevinIntegrator *intg = new LangevinIntegrator(372.4, 1.234, 0.0018);
    stringstream ss;
    XmlSerializer::serialize<Integrator>(intg, "LangevinIntegrator", ss);
    LangevinIntegrator *intg2 = dynamic_cast<LangevinIntegrator*>(XmlSerializer::deserialize<Integrator>(ss));
    ASSERT_EQUAL(intg->getConstraintTolerance(), intg2->getConstraintTolerance());
    ASSERT_EQUAL(intg->getStepSize(), intg2->getStepSize());
    ASSERT_EQUAL(intg->getTemperature(), intg2->getTemperature());
    ASSERT_EQUAL(intg->getFriction(), intg2->getFriction());
    ASSERT_EQUAL(intg->getRandomNumberSeed(), intg2->getRandomNumberSeed());
    delete intg;
    delete intg2;
}

void testSerializeBrownianIntegrator() {
    BrownianIntegrator *intg = new BrownianIntegrator(243.1, 3.234, 0.0021);
    stringstream ss;
    XmlSerializer::serialize<Integrator>(intg, "BrownianIntegrator", ss);
    BrownianIntegrator *intg2 = dynamic_cast<BrownianIntegrator*>(XmlSerializer::deserialize<Integrator>(ss));
    ASSERT_EQUAL(intg->getConstraintTolerance(), intg2->getConstraintTolerance());
    ASSERT_EQUAL(intg->getStepSize(), intg2->getStepSize());
    ASSERT_EQUAL(intg->getTemperature(), intg2->getTemperature());
    ASSERT_EQUAL(intg->getFriction(), intg2->getFriction());
    ASSERT_EQUAL(intg->getRandomNumberSeed(), intg2->getRandomNumberSeed());
    delete intg;
    delete intg2;
}

void testSerializeVariableVerletIntegrator() {
    VariableVerletIntegrator *intg = new VariableVerletIntegrator(0.04234);
    stringstream ss;
    XmlSerializer::serialize<Integrator>(intg, "VariableVerletIntegrator", ss);
    VariableVerletIntegrator *intg2 = dynamic_cast<VariableVerletIntegrator*>(XmlSerializer::deserialize<Integrator>(ss));
    ASSERT_EQUAL(intg->getConstraintTolerance(), intg2->getConstraintTolerance());
    ASSERT_EQUAL(intg->getStepSize(), intg2->getStepSize());
    ASSERT_EQUAL(intg->getErrorTolerance(), intg2->getErrorTolerance());
    delete intg;
    delete intg2;
}

void testSerializeVariableLangevinIntegrator() {
    VariableLangevinIntegrator *intg = new VariableLangevinIntegrator(243.1, 3.234, 0.0021);
    stringstream ss;
    XmlSerializer::serialize<Integrator>(intg, "VariableLangevinIntegrator", ss);
    VariableLangevinIntegrator *intg2 = dynamic_cast<VariableLangevinIntegrator*>(XmlSerializer::deserialize<Integrator>(ss));
    ASSERT_EQUAL(intg->getConstraintTolerance(), intg2->getConstraintTolerance());
    ASSERT_EQUAL(intg->getStepSize(), intg2->getStepSize());
    ASSERT_EQUAL(intg->getErrorTolerance(), intg2->getErrorTolerance());
    ASSERT_EQUAL(intg->getFriction(), intg2->getFriction());
    ASSERT_EQUAL(intg->getTemperature(), intg2->getTemperature());
    ASSERT_EQUAL(intg->getRandomNumberSeed(), intg2->getRandomNumberSeed());
    delete intg;
    delete intg2;
}

void testSerializeCustomIntegrator() {
    CustomIntegrator *intg = new CustomIntegrator(0.002234);
    intg->addPerDofVariable("temp",0);
    vector<Vec3> initialValues(123);
    for(int i = 0; i < 123; i++)
        initialValues[i] = Vec3(i+0.1, i+0.2, i+0.3);
    intg->setPerDofVariable(0, initialValues);
    intg->addPerDofVariable("oldx", 0);
    intg->addComputePerDof("v", "v+dt*f/m");
    intg->addComputePerDof("oldx", "x");
    intg->addComputePerDof("x", "x+dt*v");
    intg->addConstrainPositions();
    intg->addComputePerDof("v", "(x-oldx)/dt");
    intg->addUpdateContextState();
    intg->addConstrainVelocities();
    intg->addComputeSum("summand", "x*x+v*v");
    intg->addPerDofVariable("outf", 0);
    intg->addPerDofVariable("outf1", 0);
    intg->addPerDofVariable("outf2", 0);
    intg->addGlobalVariable("oute", 0);
    intg->addGlobalVariable("oute1", 0);
    intg->addGlobalVariable("oute2", 0);
    intg->addComputePerDof("outf", "f");
    intg->addComputePerDof("outf1", "f1");
    intg->addComputePerDof("outf2", "f2");
    intg->addComputeGlobal("oute", "energy");
    intg->addComputeGlobal("oute1", "energy1");
    intg->addComputeGlobal("oute2", "energy2");
    intg->addUpdateContextState();
    intg->addConstrainVelocities();
    intg->addComputeSum("summand2", "v*v+f*f");
    intg->setConstraintTolerance(1e-5);
    intg->setKineticEnergyExpression("m*v1*v1/2; v1=v+0.5*dt*f/m");
    stringstream ss;
    XmlSerializer::serialize<Integrator>(intg, "CustomIntegrator", ss);
    CustomIntegrator *intg2 = dynamic_cast<CustomIntegrator*>(XmlSerializer::deserialize<Integrator>(ss));
    ASSERT_EQUAL(intg->getNumGlobalVariables(), intg->getNumGlobalVariables());
    for (int i = 0; i < intg->getNumGlobalVariables(); i++) {
        ASSERT_EQUAL(intg->getGlobalVariable(i), intg2->getGlobalVariable(i));
        ASSERT_EQUAL(intg->getGlobalVariableName(i), intg2->getGlobalVariableName(i));
    }
    ASSERT_EQUAL(intg->getNumPerDofVariables(), intg2->getNumPerDofVariables());
    for(int i = 0; i < intg->getNumPerDofVariables(); i++) {
        vector<Vec3> vars1; intg->getPerDofVariable(i, vars1);
        vector<Vec3> vars2; intg2->getPerDofVariable(i, vars2);
        ASSERT_EQUAL(vars1.size(),vars2.size());
        for (int j = 0; j < (int) vars1.size(); j++) {
            ASSERT_EQUAL(vars1[j][0], vars2[j][0]);
            ASSERT_EQUAL(vars1[j][1], vars2[j][1]);
            ASSERT_EQUAL(vars1[j][2], vars2[j][2]);
        }
    }
    ASSERT_EQUAL(intg->getNumComputations(), intg2->getNumComputations());
    for(int i=0; i<intg->getNumComputations(); i++) {
        CustomIntegrator::ComputationType type1, type2;
        string variable1, variable2;
        string expression1, expression2;
        intg->getComputationStep(i, type1, variable1, expression1);
        intg2->getComputationStep(i, type2, variable2, expression2);
        ASSERT_EQUAL(type1, type2);
        ASSERT_EQUAL(variable1, variable2);
        ASSERT_EQUAL(expression1, expression2);
    }
    ASSERT_EQUAL(intg->getKineticEnergyExpression(), intg2->getKineticEnergyExpression());
    ASSERT_EQUAL(intg->getRandomNumberSeed(), intg2->getRandomNumberSeed());
    ASSERT_EQUAL(intg->getStepSize(), intg2->getStepSize());
    ASSERT_EQUAL(intg->getConstraintTolerance(), intg2->getConstraintTolerance());
    delete intg;
    delete intg2;
}

int main() {
    try {
        testSerializeBrownianIntegrator();
        testSerializeCustomIntegrator();
        testSerializeVerletIntegrator();
        testSerializeVariableLangevinIntegrator();
        testSerializeVariableVerletIntegrator();
        testSerializeLangevinIntegrator(); 
    }
    catch(const exception& e) {
		return 1;
    }
    return 0;
}


