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
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/SplineFitter.h"
#include <cmath>
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

void testNaturalSpline() {
    vector<double> x(20);
    vector<double> y(20);
    for (unsigned int i = 0; i < x.size(); i++) {
        x[i] = 0.5*i+0.1*sin(double(i));
        y[i] = sin(double(x[i]));
    }
    vector<double> deriv;
    SplineFitter::createNaturalSpline(x, y, deriv);
    ASSERT_EQUAL_TOL(deriv[0], 0.0, 1e-6);
    ASSERT_EQUAL_TOL(deriv[deriv.size()-1], 0.0, 1e-6);
    for (unsigned int i = 0; i < x.size(); i++) {
        ASSERT_EQUAL_TOL(y[i], SplineFitter::evaluateSpline(x, y, deriv, x[i]), 1e-6);
        ASSERT_EQUAL_TOL(cos(x[i]), SplineFitter::evaluateSplineDerivative(x, y, deriv, x[i]), 0.05);
    }
    for (int i = 1; i < 9; i++) {
        ASSERT_EQUAL_TOL(sin((double)i), SplineFitter::evaluateSpline(x, y, deriv, i), 0.05);
        ASSERT_EQUAL_TOL(cos((double)i), SplineFitter::evaluateSplineDerivative(x, y, deriv, i), 0.05);
    }
}

void testPeriodicSpline() {
    vector<double> x(26);
    vector<double> y(26);
    for (unsigned int i = 0; i < x.size()-1; i++) {
        x[i] = 0.5*i+0.1*sin((double)i);
        y[i] = sin((double)x[i]);
    }
    x[x.size()-1] = 4*M_PI;
    y[y.size()-1] = y[0];
    vector<double> deriv;
    SplineFitter::createPeriodicSpline(x, y, deriv);
    ASSERT_EQUAL_TOL(deriv[0], deriv[deriv.size()-1], 1e-6);
    for (unsigned int i = 0; i < x.size(); i++) {
        ASSERT_EQUAL_TOL(y[i], SplineFitter::evaluateSpline(x, y, deriv, x[i]), 1e-6);
        ASSERT_EQUAL_TOL(cos(x[i]), SplineFitter::evaluateSplineDerivative(x, y, deriv, x[i]), 0.05);
    }
    for (int i = 1; i < 9; i++) {
        ASSERT_EQUAL_TOL(sin((double)i), SplineFitter::evaluateSpline(x, y, deriv, i), 0.05);
        ASSERT_EQUAL_TOL(cos((double)i), SplineFitter::evaluateSplineDerivative(x, y, deriv, i), 0.05);
    }
}

void test2DSpline() {
    const int xsize = 15;
    const int ysize = 17;
    vector<double> x(xsize);
    vector<double> y(ysize);
    vector<double> f(xsize*ysize);
    for (int i = 0; i < xsize; i++)
        x[i] = 0.5*i+0.1*sin(double(i));
    for (int i = 0; i < ysize; i++)
        y[i] = 0.6*i+0.1*sin(double(i));
    for (int i = 0; i < xsize; i++)
        for (int j = 0; j < ysize; j++)
            f[i+j*xsize] = sin(x[i])*cos(0.4*y[j]);
    vector<vector<double> > c;
    SplineFitter::create2DNaturalSpline(x, y, f, c);
    for (int i = 0; i < xsize; i++)
        for (int j = 0; j < ysize; j++) {
            double value = SplineFitter::evaluate2DSpline(x, y, f, c, x[i], y[j]);
            ASSERT_EQUAL_TOL(f[i+j*xsize], value, 1e-6);
        }
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            double s = x[0]+(i+1)*(x[xsize-1]-x[0])/12.0;
            double t = y[0]+(j+1)*(y[ysize-1]-y[0])/12.0;
            double value = SplineFitter::evaluate2DSpline(x, y, f, c, s, t);
            ASSERT_EQUAL_TOL(sin(s)*cos(0.4*t), value, 0.02);
            double dx, dy;
            SplineFitter::evaluate2DSplineDerivatives(x, y, f, c, s, t, dx, dy);
            ASSERT_EQUAL_TOL(cos(s)*cos(0.4*t), dx, 0.05);
            ASSERT_EQUAL_TOL(-0.4*sin(s)*sin(0.4*t), dy, 0.05);
        }
    }
}

void test3DSpline() {
    const int xsize = 8;
    const int ysize = 9;
    const int zsize = 10;
    vector<double> x(xsize);
    vector<double> y(ysize);
    vector<double> z(zsize);
    vector<double> f(xsize*ysize*zsize);
    for (int i = 0; i < xsize; i++)
        x[i] = 0.2*i+0.02*sin(0.4*double(i));
    for (int i = 0; i < ysize; i++)
        y[i] = 0.2*i+0.02*sin(0.45*double(i));
    for (int i = 0; i < zsize; i++)
        z[i] = 0.2*i+0.02*sin(0.5*double(i));
    for (int i = 0; i < xsize; i++)
        for (int j = 0; j < ysize; j++)
            for (int k = 0; k < zsize; k++)
                f[i+j*xsize+k*xsize*ysize] = sin(x[i])*cos(0.4*y[j])*(1+z[k]);
    vector<vector<double> > c;
    SplineFitter::create3DNaturalSpline(x, y, z, f, c);
    for (int i = 0; i < xsize; i++)
        for (int j = 0; j < ysize; j++) {
            for (int k = 0; k < zsize; k++) {
                double value = SplineFitter::evaluate3DSpline(x, y, z, f, c, x[i], y[j], z[k]);
                ASSERT_EQUAL_TOL(f[i+j*xsize+k*xsize*ysize], value, 1e-6);
            }
        }
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            for (int k = 0; k < 10; k++) {
                double s = x[0]+(i+1)*(x[xsize-1]-x[0])/12.0;
                double t = y[0]+(j+1)*(y[ysize-1]-y[0])/12.0;
                double u = z[0]+(k+1)*(z[zsize-1]-z[0])/12.0;
                double value = SplineFitter::evaluate3DSpline(x, y, z, f, c, s, t, u);
                ASSERT_EQUAL_TOL(sin(s)*cos(0.4*t)*(1+u), value, 0.02);
                double dx, dy, dz;
                SplineFitter::evaluate3DSplineDerivatives(x, y, z, f, c, s, t, u, dx, dy, dz);
                ASSERT_EQUAL_TOL(cos(s)*cos(0.4*t)*(1+u), dx, 0.1);
                ASSERT_EQUAL_TOL(-0.4*sin(s)*sin(0.4*t)*(1+u), dy, 0.1);
                ASSERT_EQUAL_TOL(sin(s)*cos(0.4*t), dz, 0.1);
            }
        }
    }
}

int main() {
    try {
        testNaturalSpline();
        testPeriodicSpline();
        test2DSpline();
        test3DSpline();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

