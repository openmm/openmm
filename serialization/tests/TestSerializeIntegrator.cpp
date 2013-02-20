/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
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
#include "openmm/LangevinIntegrator.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>
#include <stdlib.h>

#include <fstream>

using namespace OpenMM;
using namespace std;

void testSerialization() {

	{
	VerletIntegrator *vInt = new VerletIntegrator(0.00342);
	stringstream ss;
	XmlSerializer::serialize<Integrator>(vInt, "VerletIntegrator", ss);
	VerletIntegrator *vInt2 = dynamic_cast<VerletIntegrator*>(XmlSerializer::deserialize<Integrator>(ss));
	ASSERT_EQUAL(vInt->getConstraintTolerance(), vInt2->getConstraintTolerance());
	ASSERT_EQUAL(vInt->getStepSize(), vInt2->getStepSize());
	delete vInt;
	delete vInt2;
	}
	{
	LangevinIntegrator *lInt = new LangevinIntegrator(372.4, 1.234, 0.0018);
	stringstream ss;
	XmlSerializer::serialize<Integrator>(lInt, "LangevinIntegrator", ss);
	LangevinIntegrator *lInt2 = dynamic_cast<LangevinIntegrator*>(XmlSerializer::deserialize<Integrator>(ss));
	ASSERT_EQUAL(lInt->getConstraintTolerance(), lInt2->getConstraintTolerance());
	ASSERT_EQUAL(lInt->getStepSize(), lInt2->getStepSize());
	ASSERT_EQUAL(lInt->getTemperature(), lInt2->getTemperature());
	ASSERT_EQUAL(lInt->getFriction(), lInt2->getFriction());
	delete lInt;
	delete lInt2;
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


