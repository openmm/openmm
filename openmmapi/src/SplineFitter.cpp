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

#include <vector>

#include "openmm/internal/SplineFitter.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;
using namespace std;

void SplineFitter::createNaturalSpline(const vector<double>& x, const vector<double>& y, vector<double>& deriv) {
    int n = x.size();
    if (y.size() != n)
        throw OpenMMException("createNaturalSpline: x and y vectors must have same length");
    if (n < 2)
        throw OpenMMException("createNaturalSpline: the length of the input array must be at least 2");
    deriv.resize(n);
    if (n == 2) {
        // This is just a straight line.

        deriv[0] = 0;
        deriv[1] = 0;
    }

    // Create the system of equations to solve.

    vector<double> a(n), b(n), c(n), rhs(n);
    a[0] = 0.0;
    b[0] = 1.0;
    c[0] = 0.0;
    rhs[0] = 0.0;
    for (int i = 1; i < n-1; i++) {
        a[i] = x[i]-x[i-1];
        b[i] = 2.0*(x[i+1]-x[i-1]);
        c[i] = x[i+1]-x[i];
        rhs[i] = 6.0*((y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]));
    }
    a[n-1] = 0.0;
    b[n-1] = 1.0;
    c[n-1] = 0.0;
    rhs[n-1] = 0.0;

    // Solve them.

    solveTridiagonalMatrix(a, b, c, rhs, deriv);
}

void SplineFitter::createPeriodicSpline(const vector<double>& x, const vector<double>& y, vector<double>& deriv) {
    int n = x.size();
    if (y.size() != n)
        throw OpenMMException("createPeriodicSpline: x and y vectors must have same length");
    if (n < 3)
        throw OpenMMException("createPeriodicSpline: the length of the input array must be at least 3");
    if (y[0] != y[n-1])
        throw OpenMMException("createPeriodicSpline: the first and last points must have the same value");
    deriv.resize(n);

    // Create the system of equations to solve.

    vector<double> a(n), b(n), c(n), rhs(n);
    a[0] = 0.0;
    b[0] = 2.0*(x[1]-x[0]+x[n-1]-x[n-2]);
    c[0] = x[1]-x[0];
    rhs[0] = 6.0*((y[1]-y[0])/(x[1]-x[0]) - (y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
    for (int i = 1; i < n-1; i++) {
        a[i] = x[i]-x[i-1];
        b[i] = 2.0*(x[i+1]-x[i-1]);
        c[i] = x[i+1]-x[i];
        rhs[i] = 6.0*((y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]));
    }
    a[n-1] = 0.0;
    b[n-1] = 1.0;
    c[n-1] = 0.0;
    rhs[n-1] = 0.0;
    double beta = x[n-1]-x[n-2];
    double alpha = -1.0;
    double gamma = -b[0];

    // This is a cyclic tridiagonal matrix.  We solve it using the Sherman-Morrison method,
    // which involves solving two tridiagonal systems.

    b[0] -= gamma;
    b[n-1] -= alpha*beta/gamma;
    solveTridiagonalMatrix(a, b, c, rhs, deriv);
    vector<double> u(n, 0.0), z(n);
    u[0] = gamma;
    u[n-1] = alpha;
    solveTridiagonalMatrix(a, b, c, u, z);
    double scale = (deriv[0]+beta*deriv[n-1]/gamma)/(1.0+z[0]+beta*z[n-1]/gamma);
    for (int i = 0; i < n; i++)
        deriv[i] -= scale*z[i];
}

double SplineFitter::evaluateSpline(const vector<double>& x, const vector<double>& y, const vector<double>& deriv, double t) {
    int n = x.size();
    if (t < x[0] || t > x[n-1])
        throw OpenMMException("evaluateSpline: specified point is outside the range defined by the spline");

    // Perform a binary search to identify the interval containing the point to evaluate.

    int lower = 0;
    int upper = n-1;
    while (upper-lower > 1) {
        int middle = (upper+lower)/2;
        if (x[middle] > t)
            upper = middle;
        else
            lower = middle;
    }

    // Evaluate the spline.

    double dx = x[upper]-x[lower];
    double a = (x[upper]-t)/dx;
    double b = 1.0-a;
    return a*y[lower]+b*y[upper]+((a*a*a-a)*deriv[lower] + (b*b*b-b)*deriv[upper])*dx*dx/6.0;
}

double SplineFitter::evaluateSplineDerivative(const vector<double>& x, const vector<double>& y, const vector<double>& deriv, double t) {
    int n = x.size();
    if (t < x[0] || t > x[n-1])
        throw OpenMMException("evaluateSplineDerivative: specified point is outside the range defined by the spline");

    // Perform a binary search to identify the interval containing the point to evaluate.

    int lower = 0;
    int upper = n-1;
    while (upper-lower > 1) {
        int middle = (upper+lower)/2;
        if (x[middle] > t)
            upper = middle;
        else
            lower = middle;
    }

    // Evaluate the spline.

    double dx = x[upper]-x[lower];
    double a = (x[upper]-t)/dx;
    double b = 1.0-a;
    double dadx = -1.0/dx;
    return dadx*y[lower]-dadx*y[upper]+((1.0-3.0*a*a)*deriv[lower] + (3.0*b*b-1.0)*deriv[upper])*dx/6.0;
}

void SplineFitter::solveTridiagonalMatrix(const vector<double>& a, const vector<double>& b, const vector<double>& c, const vector<double>& rhs, vector<double>& sol) {
    int n = a.size();
    vector<double> gamma(n);

    // Decompose the matrix.

    sol[0] = rhs[0]/b[0];
    double beta = b[0];
    for (int i = 1; i < n; i++) {
        gamma[i] = c[i-1]/beta;
        beta = b[i]-a[i]*gamma[i];
        sol[i] = (rhs[i]-a[i]*sol[i-1])/beta;
    }

    // Perform backsubstitation.

    for (int i = n-2; i >= 0; i--)
        sol[i] -= gamma[i+1]*sol[i+1];
}

void SplineFitter::create2DNaturalSpline(const vector<double>& x, const vector<double>& y, const vector<double>& values, vector<vector<double> >& c) {
    int xsize = x.size(), ysize = y.size();
    if (xsize < 2 || ysize < 2)
        throw OpenMMException("create2DNaturalSpline: must have at least two points along each axis");
    if (values.size() != xsize*ysize)
        throw OpenMMException("create2DNaturalSpline: incorrect number of values");
    vector<double> d1(xsize*ysize), d2(xsize*ysize), d12(xsize*ysize);
    vector<double> t(xsize), deriv(xsize);

    // Compute derivatives with respect to x.

    for (int i = 0; i < ysize; i++) {
        for (int j = 0; j < xsize; j++)
            t[j] = values[j+xsize*i];
        SplineFitter::createNaturalSpline(x, t, deriv);
        for (int j = 0; j < xsize; j++)
            d1[j+xsize*i] = SplineFitter::evaluateSplineDerivative(x, t, deriv, x[j]);
    }

    // Compute derivatives with respect to y.

    t.resize(ysize);
    deriv.resize(ysize);
    for (int i = 0; i < xsize; i++) {
        for (int j = 0; j < ysize; j++)
            t[j] = values[i+xsize*j];
        SplineFitter::createNaturalSpline(y, t, deriv);
        for (int j = 0; j < ysize; j++)
            d2[i+xsize*j] = SplineFitter::evaluateSplineDerivative(y, t, deriv, y[j]);
    }

    // Compute cross derivatives.

    t.resize(xsize);
    deriv.resize(xsize);
    for (int i = 0; i < ysize; i++) {
        for (int j = 0; j < xsize; j++)
            t[j] = d2[j+xsize*i];
        SplineFitter::createNaturalSpline(x, t, deriv);
        for (int j = 0; j < xsize; j++)
            d12[j+xsize*i] = SplineFitter::evaluateSplineDerivative(x, t, deriv, x[j]);
    }

    // Now compute the coefficients.

    const int wt[] = {
        1, 0, -3, 2, 0, 0, 0, 0, -3, 0, 9, -6, 2, 0, -6, 4,
        0, 0, 0, 0, 0, 0, 0, 0, 3, 0, -9, 6, -2, 0, 6, -4,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, -6, 0, 0, -6, 4,
        0, 0, 3, -2, 0, 0, 0, 0, 0, 0, -9, 6, 0, 0, 6, -4,
        0, 0, 0, 0, 1, 0, -3, 2, -2, 0, 6, -4, 1, 0, -3, 2,
        0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 3, -2, 1, 0, -3, 2,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 2, 0, 0, 3, -2,
        0, 0, 0, 0, 0, 0, 3, -2, 0, 0, -6, 4, 0, 0, 3, -2,
        0, 1, -2, 1, 0, 0, 0, 0, 0, -3, 6, -3, 0, 2, -4, 2,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 3, -6, 3, 0, -2, 4, -2,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, 2, -2,
        0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 3, -3, 0, 0, -2, 2,
        0, 0, 0, 0, 0, 1, -2, 1, 0, -2, 4, -2, 0, 1, -2, 1,
        0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 2, -1, 0, 1, -2, 1,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, -1, 1,
        0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 2, -2, 0, 0, -1, 1
    };
    vector<double> rhs(16);
    c.resize((xsize-1)*(ysize-1));
    for (int i = 0; i < xsize-1; i++) {
        for (int j = 0; j < ysize-1; j++) {
            // Compute the 16 coefficients for patch (i, j).

            int nexti = i+1;
            int nextj = j+1;
            double deltax = x[nexti]-x[i];
            double deltay = y[nextj]-y[j];
            double e[] = {values[i+j*xsize], values[nexti+j*xsize], values[nexti+nextj*xsize], values[i+nextj*xsize]};
            double e1[] = {d1[i+j*xsize], d1[nexti+j*xsize], d1[nexti+nextj*xsize], d1[i+nextj*xsize]};
            double e2[] = {d2[i+j*xsize], d2[nexti+j*xsize], d2[nexti+nextj*xsize], d2[i+nextj*xsize]};
            double e12[] = {d12[i+j*xsize], d12[nexti+j*xsize], d12[nexti+nextj*xsize], d12[i+nextj*xsize]};

            for (int k = 0; k < 4; k++) {
                rhs[k] = e[k];
                rhs[k+4] = e1[k]*deltax;
                rhs[k+8] = e2[k]*deltay;
                rhs[k+12] = e12[k]*deltax*deltay;
            }
            vector<double>& coeff = c[i+j*(xsize-1)];
            coeff.resize(16);
            for (int k = 0; k < 16; k++) {
                double sum = 0.0;
                for (int m = 0; m < 16; m++)
                    sum += wt[k+16*m]*rhs[m];
                coeff[k] = sum;
            }
        }
    }
}

double SplineFitter::evaluate2DSpline(const vector<double>& x, const vector<double>& y, const vector<double>& values, const vector<vector<double> >& c, double u, double v) {
    int xsize = x.size();
    int ysize = y.size();
    if (u < x[0] || u > x[xsize-1] || v < y[0] || v > y[ysize-1])
        throw OpenMMException("evaluate2DSpline: specified point is outside the range defined by the spline");

    // Perform a binary search to identify the interval containing the point to evaluate.

    int lowerx = 0;
    int upperx = xsize-1;
    while (upperx-lowerx > 1) {
        int middle = (upperx+lowerx)/2;
        if (x[middle] > u)
            upperx = middle;
        else
            lowerx = middle;
    }
    int lowery = 0;
    int uppery = ysize-1;
    while (uppery-lowery > 1) {
        int middle = (uppery+lowery)/2;
        if (y[middle] > v)
            uppery = middle;
        else
            lowery = middle;
    }
    double deltax = x[upperx]-x[lowerx];
    double deltay = y[uppery]-y[lowery];
    double da = (u-x[lowerx])/deltax;
    double db = (v-y[lowery])/deltay;
    const vector<double>& coeff = c[lowerx+(xsize-1)*lowery];

    // Evaluate the spline to determine the value and gradients.

    double value = 0;
    for (int i = 3; i >= 0; i--)
        value = da*value + ((coeff[i*4+3]*db + coeff[i*4+2])*db + coeff[i*4+1])*db + coeff[i*4+0];
    return value;
}

void SplineFitter::evaluate2DSplineDerivatives(const vector<double>& x, const vector<double>& y, const vector<double>& values, const vector<vector<double> >& c, double u, double v, double& dx, double &dy) {
    int xsize = x.size();
    int ysize = y.size();
    if (u < x[0] || u > x[xsize-1] || v < y[0] || v > y[ysize-1])
        throw OpenMMException("evaluate2DSplineDerivatives: specified point is outside the range defined by the spline");

    // Perform a binary search to identify the interval containing the point to evaluate.

    int lowerx = 0;
    int upperx = xsize-1;
    while (upperx-lowerx > 1) {
        int middle = (upperx+lowerx)/2;
        if (x[middle] > u)
            upperx = middle;
        else
            lowerx = middle;
    }
    int lowery = 0;
    int uppery = ysize-1;
    while (uppery-lowery > 1) {
        int middle = (uppery+lowery)/2;
        if (y[middle] > v)
            uppery = middle;
        else
            lowery = middle;
    }
    double deltax = x[upperx]-x[lowerx];
    double deltay = y[uppery]-y[lowery];
    double da = (u-x[lowerx])/deltax;
    double db = (v-y[lowery])/deltay;
    const vector<double>& coeff = c[lowerx+(xsize-1)*lowery];

    // Evaluate the spline to determine the value and gradients.

    dx = 0;
    dy = 0;
    for (int i = 3; i >= 0; i--) {
        dx = db*dx + (3.0*coeff[i+3*4]*da + 2.0*coeff[i+2*4])*da + coeff[i+1*4];
        dy = da*dy + (3.0*coeff[i*4+3]*db + 2.0*coeff[i*4+2])*db + coeff[i*4+1];
    }
    dx /= deltax;
    dy /= deltay;
}
