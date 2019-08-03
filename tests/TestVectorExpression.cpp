/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2018 Stanford University and the Authors.           *
 * Authors: Robert McGibbon                                                   *
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
#include "openmm/internal/VectorExpression.h"
#include <iostream>

using namespace OpenMM;
using namespace std;

void verifyEvaluation(const string& expression, Vec3 expectedValue, Vec3 x=Vec3(), Vec3 y=Vec3()) {
    map<string, Vec3> variables;
    variables["x"] = x;
    variables["y"] = y;
    map<string, Lepton::CustomFunction*> customFunctions;
    VectorExpression expr(expression, customFunctions);
    Vec3 value = expr.evaluate(variables);
    ASSERT_EQUAL_VEC(expectedValue, value, 1e-10);
}

void testExpressions() {
    verifyEvaluation("5", Vec3(5, 5, 5));
    verifyEvaluation("10*x", Vec3(50, 100, 150), Vec3(5, 10, 15));
    verifyEvaluation("2*dot(x, 3)", Vec3(180, 180, 180), Vec3(5, 10, 15));
    Vec3 a(1, 1.5, 2);
    Vec3 b(3, -1, 4);
    verifyEvaluation("cross(y, x)", b.cross(a), a, b);
    verifyEvaluation("_x(x)", Vec3(a[0], a[0], a[0]), a);
    verifyEvaluation("_y(x)", Vec3(a[1], a[1], a[1]), a);
    verifyEvaluation("_z(x)", Vec3(a[2], a[2], a[2]), a);
    verifyEvaluation("vector(x, 5, y)", Vec3(a[0], 5, b[2]), a, b);
}

int main(int argc, char* argv[]) {
    try {
        testExpressions();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
