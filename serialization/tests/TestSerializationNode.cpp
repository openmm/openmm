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
#include "openmm/serialization/SerializationNode.h"
#include <iostream>

using namespace OpenMM;
using namespace std;

void testProperties() {
    SerializationNode node;
    ASSERT_EQUAL(false, node.hasProperty("prop1"));
    ASSERT_EQUAL(false, node.hasProperty("prop2"));
    bool exists = false;
    try {
        node.getStringProperty("prop1");
        exists = true;
    }
    catch (const exception& ex) {
    }
    ASSERT_EQUAL(false, exists);
    try {
        node.getIntProperty("prop1");
        exists = true;
    }
    catch (const exception& ex) {
    }
    ASSERT_EQUAL(false, exists);
    try {
        node.getDoubleProperty("prop1");
        exists = true;
    }
    catch (const exception& ex) {
    }
    ASSERT_EQUAL(false, exists);
    ASSERT_EQUAL(3, node.getIntProperty("prop1", 3));
    ASSERT_EQUAL(3.5, node.getDoubleProperty("prop1", 3.5));
    ASSERT_EQUAL("abc", node.getStringProperty("prop1", "abc"));
    node.setIntProperty("prop1", 1);
    ASSERT_EQUAL(1, node.getIntProperty("prop1"));
    node.setDoubleProperty("prop1", 1.5);
    ASSERT_EQUAL(1.5, node.getDoubleProperty("prop1"));
    node.setStringProperty("prop1", "hello");
    ASSERT_EQUAL("hello", node.getStringProperty("prop1"));
    ASSERT_EQUAL(true, node.hasProperty("prop1"));
    ASSERT_EQUAL(false, node.hasProperty("prop2"));
}

int main() {
    try {
        testProperties();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

