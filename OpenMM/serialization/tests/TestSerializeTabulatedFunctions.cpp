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

#include "openmm/TabulatedFunction.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

void testContinuous1DFunction() {
    // Create a function.

    double min = 0.5, max = 1.5;
    vector<double> values(60);
    for (int i = 0; i < (int) values.size(); i++)
        values[i] = sin((double) i);
    Continuous1DFunction function(values, min, max);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<Continuous1DFunction>(&function, "Function", buffer);
    Continuous1DFunction* copy = XmlSerializer::deserialize<Continuous1DFunction>(buffer);

    // Compare the two forces to see if they are identical.

    double min2, max2;
    vector<double> values2;
    copy->getFunctionParameters(values2, min2, max2);
    ASSERT_EQUAL(min, min2);
    ASSERT_EQUAL(max, max2);
    ASSERT_EQUAL(values.size(), values2.size());
    for (int j = 0; j < (int) values.size(); j++)
        ASSERT_EQUAL(values[j], values2[j]);
}

void testContinuous2DFunction() {
    // Create a function.

    int xsize = 5, ysize = 12;
    double xmin = 0.5, xmax = 1.5, ymin = 0.1, ymax = 5.0;
    vector<double> values(xsize*ysize);
    for (int i = 0; i < (int) values.size(); i++)
        values[i] = sin((double) i);
    Continuous2DFunction function(xsize, ysize, values, xmin, xmax, ymin, ymax);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<Continuous2DFunction>(&function, "Function", buffer);
    Continuous2DFunction* copy = XmlSerializer::deserialize<Continuous2DFunction>(buffer);

    // Compare the two forces to see if they are identical.

    int xsize2, ysize2;
    double xmin2, xmax2, ymin2, ymax2;
    vector<double> values2;
    copy->getFunctionParameters(xsize2, ysize2, values2, xmin2, xmax2, ymin2, ymax2);
    ASSERT_EQUAL(xsize, xsize2);
    ASSERT_EQUAL(ysize, ysize2);
    ASSERT_EQUAL(xmin, xmin2);
    ASSERT_EQUAL(xmax, xmax2);
    ASSERT_EQUAL(ymin, ymin2);
    ASSERT_EQUAL(ymax, ymax2);
    ASSERT_EQUAL(values.size(), values2.size());
    for (int j = 0; j < (int) values.size(); j++)
        ASSERT_EQUAL(values[j], values2[j]);
}

void testContinuous3DFunction() {
    // Create a function.

    int xsize = 5, ysize = 4, zsize = 3;
    double xmin = 0.5, xmax = 1.5, ymin = 0.1, ymax = 5.0, zmin = 0.3, zmax = 0.9;
    vector<double> values(xsize*ysize*zsize);
    for (int i = 0; i < (int) values.size(); i++)
        values[i] = sin((double) i);
    Continuous3DFunction function(xsize, ysize, zsize, values, xmin, xmax, ymin, ymax, zmin, zmax);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<Continuous3DFunction>(&function, "Function", buffer);
    Continuous3DFunction* copy = XmlSerializer::deserialize<Continuous3DFunction>(buffer);

    // Compare the two forces to see if they are identical.

    int xsize2, ysize2, zsize2;
    double xmin2, xmax2, ymin2, ymax2, zmin2, zmax2;
    vector<double> values2;
    copy->getFunctionParameters(xsize2, ysize2, zsize2, values2, xmin2, xmax2, ymin2, ymax2, zmin2, zmax2);
    ASSERT_EQUAL(xsize, xsize2);
    ASSERT_EQUAL(ysize, ysize2);
    ASSERT_EQUAL(zsize, zsize2);
    ASSERT_EQUAL(xmin, xmin2);
    ASSERT_EQUAL(xmax, xmax2);
    ASSERT_EQUAL(ymin, ymin2);
    ASSERT_EQUAL(ymax, ymax2);
    ASSERT_EQUAL(zmin, zmin2);
    ASSERT_EQUAL(zmax, zmax2);
    ASSERT_EQUAL(values.size(), values2.size());
    for (int j = 0; j < (int) values.size(); j++)
        ASSERT_EQUAL(values[j], values2[j]);
}

void testDiscrete1DFunction() {
    // Create a function.

    vector<double> values(60);
    for (int i = 0; i < (int) values.size(); i++)
        values[i] = sin((double) i);
    Discrete1DFunction function(values);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<Discrete1DFunction>(&function, "Function", buffer);
    Discrete1DFunction* copy = XmlSerializer::deserialize<Discrete1DFunction>(buffer);

    // Compare the two forces to see if they are identical.

    vector<double> values2;
    copy->getFunctionParameters(values2);
    ASSERT_EQUAL(values.size(), values2.size());
    for (int j = 0; j < (int) values.size(); j++)
        ASSERT_EQUAL(values[j], values2[j]);
}

void testDiscrete2DFunction() {
    // Create a function.

    int xsize = 5, ysize = 12;
    vector<double> values(xsize*ysize);
    for (int i = 0; i < (int) values.size(); i++)
        values[i] = sin((double) i);
    Discrete2DFunction function(xsize, ysize, values);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<Discrete2DFunction>(&function, "Function", buffer);
    Discrete2DFunction* copy = XmlSerializer::deserialize<Discrete2DFunction>(buffer);

    // Compare the two forces to see if they are identical.

    int xsize2, ysize2;
    vector<double> values2;
    copy->getFunctionParameters(xsize2, ysize2, values2);
    ASSERT_EQUAL(xsize, xsize2);
    ASSERT_EQUAL(ysize, ysize2);
    ASSERT_EQUAL(values.size(), values2.size());
    for (int j = 0; j < (int) values.size(); j++)
        ASSERT_EQUAL(values[j], values2[j]);
}

void testDiscrete3DFunction() {
    // Create a function.

    int xsize = 5, ysize = 4, zsize = 3;
    vector<double> values(xsize*ysize*zsize);
    for (int i = 0; i < (int) values.size(); i++)
        values[i] = sin((double) i);
    Discrete3DFunction function(xsize, ysize, zsize, values);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<Discrete3DFunction>(&function, "Function", buffer);
    Discrete3DFunction* copy = XmlSerializer::deserialize<Discrete3DFunction>(buffer);

    // Compare the two forces to see if they are identical.

    int xsize2, ysize2, zsize2;
    vector<double> values2;
    copy->getFunctionParameters(xsize2, ysize2, zsize2, values2);
    ASSERT_EQUAL(xsize, xsize2);
    ASSERT_EQUAL(ysize, ysize2);
    ASSERT_EQUAL(zsize, zsize2);
    ASSERT_EQUAL(values.size(), values2.size());
    for (int j = 0; j < (int) values.size(); j++)
        ASSERT_EQUAL(values[j], values2[j]);
}

int main() {
    try {
        testContinuous1DFunction();
        testContinuous2DFunction();
        testContinuous3DFunction();
        testDiscrete1DFunction();
        testDiscrete2DFunction();
        testDiscrete3DFunction();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
