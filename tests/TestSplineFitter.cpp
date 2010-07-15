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

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "AssertionUtilities.h"
#include "openmm/internal/SplineFitter.h"
#include <cmath>
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

void testNaturalSpline() {
    vector<double> x(20);
    vector<double> y(20);
    for (int i = 0; i < x.size(); i++) {
        x[i] = 0.5*i+0.1*sin(i);
        y[i] = sin(x[i]);
    }
    vector<double> deriv;
    SplineFitter::createNaturalSpline(x, y, deriv);
    ASSERT_EQUAL_TOL(deriv[0], 0.0, 1e-6);
    ASSERT_EQUAL_TOL(deriv[deriv.size()-1], 0.0, 1e-6);
    for (int i = 0; i < x.size(); i++) {
        ASSERT_EQUAL_TOL(y[i], SplineFitter::evaluateSpline(x, y, deriv, x[i]), 1e-6);
        ASSERT_EQUAL_TOL(cos(x[i]), SplineFitter::evaluateSplineDerivative(x, y, deriv, x[i]), 0.05);
    }
    for (int i = 1; i < 9; i++) {
        ASSERT_EQUAL_TOL(sin(i), SplineFitter::evaluateSpline(x, y, deriv, i), 0.05);
        ASSERT_EQUAL_TOL(cos(i), SplineFitter::evaluateSplineDerivative(x, y, deriv, i), 0.05);
    }
}

void testPeriodicSpline() {
    vector<double> x(26);
    vector<double> y(26);
    for (int i = 0; i < x.size()-1; i++) {
        x[i] = 0.5*i+0.1*sin(i);
        y[i] = sin(x[i]);
    }
    x[x.size()-1] = 4*M_PI;
    y[y.size()-1] = y[0];
    vector<double> deriv;
    SplineFitter::createPeriodicSpline(x, y, deriv);
    ASSERT_EQUAL_TOL(deriv[0], deriv[deriv.size()-1], 1e-6);
    for (int i = 0; i < x.size(); i++) {
        ASSERT_EQUAL_TOL(y[i], SplineFitter::evaluateSpline(x, y, deriv, x[i]), 1e-6);
        ASSERT_EQUAL_TOL(cos(x[i]), SplineFitter::evaluateSplineDerivative(x, y, deriv, x[i]), 0.05);
    }
    for (int i = 1; i < 9; i++) {
        ASSERT_EQUAL_TOL(sin(i), SplineFitter::evaluateSpline(x, y, deriv, i), 0.05);
        ASSERT_EQUAL_TOL(cos(i), SplineFitter::evaluateSplineDerivative(x, y, deriv, i), 0.05);
    }
}

int main() {
    try {
        testNaturalSpline();
        testPeriodicSpline();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

