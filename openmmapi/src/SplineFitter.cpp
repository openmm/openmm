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

    // Evaluate the spline to determine the value.

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

    // Evaluate the spline to determine the derivatives.

    dx = 0;
    dy = 0;
    for (int i = 3; i >= 0; i--) {
        dx = db*dx + (3.0*coeff[i+3*4]*da + 2.0*coeff[i+2*4])*da + coeff[i+1*4];
        dy = da*dy + (3.0*coeff[i*4+3]*db + 2.0*coeff[i*4+2])*db + coeff[i*4+1];
    }
    dx /= deltax;
    dy /= deltay;
}

void SplineFitter::create3DNaturalSpline(const vector<double>& x, const vector<double>& y, const vector<double>& z, const vector<double>& values, vector<vector<double> >& c) {
    int xsize = x.size(), ysize = y.size(), zsize = z.size();
    int xysize = xsize*ysize;
    if (xsize < 2 || ysize < 2 || zsize < 2)
        throw OpenMMException("create2DNaturalSpline: must have at least two points along each axis");
    if (values.size() != xsize*ysize*zsize)
        throw OpenMMException("create2DNaturalSpline: incorrect number of values");
    vector<double> d1(xsize*ysize*zsize), d2(xsize*ysize*zsize), d3(xsize*ysize*zsize);
    vector<double> d12(xsize*ysize*zsize), d13(xsize*ysize*zsize), d23(xsize*ysize*zsize), d123(xsize*ysize*zsize);
    vector<double> t(xsize), deriv(xsize);

    // Compute derivatives with respect to x.

    for (int i = 0; i < ysize; i++) {
        for (int j = 0; j < zsize; j++) {
            for (int k = 0; k < xsize; k++)
                t[k] = values[k+xsize*i+xysize*j];
            SplineFitter::createNaturalSpline(x, t, deriv);
            for (int k = 0; k < xsize; k++)
                d1[k+xsize*i+xysize*j] = SplineFitter::evaluateSplineDerivative(x, t, deriv, x[k]);
        }
    }

    // Compute derivatives with respect to y.

    t.resize(ysize);
    deriv.resize(ysize);
    for (int i = 0; i < xsize; i++) {
        for (int j = 0; j < zsize; j++) {
            for (int k = 0; k < ysize; k++)
                t[k] = values[i+xsize*k+xysize*j];
            SplineFitter::createNaturalSpline(y, t, deriv);
            for (int k = 0; k < ysize; k++)
                d2[i+xsize*k+xysize*j] = SplineFitter::evaluateSplineDerivative(y, t, deriv, y[k]);
        }
    }

    // Compute derivatives with respect to z.

    t.resize(zsize);
    deriv.resize(zsize);
    for (int i = 0; i < xsize; i++) {
        for (int j = 0; j < ysize; j++) {
            for (int k = 0; k < zsize; k++)
                t[k] = values[i+xsize*j+xysize*k];
            SplineFitter::createNaturalSpline(z, t, deriv);
            for (int k = 0; k < zsize; k++)
                d3[i+xsize*j+xysize*k] = SplineFitter::evaluateSplineDerivative(z, t, deriv, z[k]);
        }
    }

    // Compute second derivatives with respect to x and y.

    t.resize(xsize);
    deriv.resize(xsize);
    for (int i = 0; i < ysize; i++) {
        for (int j = 0; j < zsize; j++) {
            for (int k = 0; k < xsize; k++)
                t[k] = d2[k+xsize*i+xysize*j];
            SplineFitter::createNaturalSpline(x, t, deriv);
            for (int k = 0; k < xsize; k++)
                d12[k+xsize*i+xysize*j] = SplineFitter::evaluateSplineDerivative(x, t, deriv, x[k]);
        }
    }

    // Compute second derivatives with respect to y and z.

    t.resize(ysize);
    deriv.resize(ysize);
    for (int i = 0; i < zsize; i++) {
        for (int j = 0; j < xsize; j++) {
            for (int k = 0; k < ysize; k++)
                t[k] = d3[j+xsize*k+xysize*i];
            SplineFitter::createNaturalSpline(y, t, deriv);
            for (int k = 0; k < ysize; k++)
                d23[j+xsize*k+xysize*i] = SplineFitter::evaluateSplineDerivative(y, t, deriv, y[k]);
        }
    }

    // Compute second derivatives with respect to x and z.

    t.resize(zsize);
    deriv.resize(zsize);
    for (int i = 0; i < xsize; i++) {
        for (int j = 0; j < ysize; j++) {
            for (int k = 0; k < zsize; k++)
                t[k] = d1[i+xsize*j+xysize*k];
            SplineFitter::createNaturalSpline(z, t, deriv);
            for (int k = 0; k < zsize; k++)
                d13[i+xsize*j+xysize*k] = SplineFitter::evaluateSplineDerivative(z, t, deriv, z[k]);
        }
    }

    // Compute third derivatives with respect to x, y, and z.

    t.resize(xsize);
    deriv.resize(xsize);
    for (int i = 0; i < ysize; i++) {
        for (int j = 0; j < zsize; j++) {
            for (int k = 0; k < xsize; k++)
                t[k] = d23[k+xsize*i+xysize*j];
            SplineFitter::createNaturalSpline(x, t, deriv);
            for (int k = 0; k < xsize; k++)
                d123[k+xsize*i+xysize*j] = SplineFitter::evaluateSplineDerivative(x, t, deriv, x[k]);
        }
    }

    // Now compute the coefficients.  This involves multiplying by a sparse 64x64 matrix, given
    // here in packed form.

    const int wt[] = {
        1,0,1,
        1,8,1,
        4,0,-3,1,3,8,-2,9,-1,
        4,0,2,1,-2,8,1,9,1,
        1,16,1,
        1,32,1,
        4,16,-3,17,3,32,-2,33,-1,
        4,16,2,17,-2,32,1,33,1,
        4,0,-3,2,3,16,-2,18,-1,
        4,8,-3,10,3,32,-2,34,-1,
        16,0,9,1,-9,2,-9,3,9,8,6,9,3,10,-6,11,-3,16,6,17,-6,18,3,19,-3,32,4,33,2,34,2,35,1,
        16,0,-6,1,6,2,6,3,-6,8,-3,9,-3,10,3,11,3,16,-4,17,4,18,-2,19,2,32,-2,33,-2,34,-1,35,-1,
        4,0,2,2,-2,16,1,18,1,
        4,8,2,10,-2,32,1,34,1,
        16,0,-6,1,6,2,6,3,-6,8,-4,9,-2,10,4,11,2,16,-3,17,3,18,-3,19,3,32,-2,33,-1,34,-2,35,-1,
        16,0,4,1,-4,2,-4,3,4,8,2,9,2,10,-2,11,-2,16,2,17,-2,18,2,19,-2,32,1,33,1,34,1,35,1,
        1,24,1,
        1,40,1,
        4,24,-3,25,3,40,-2,41,-1,
        4,24,2,25,-2,40,1,41,1,
        1,48,1,
        1,56,1,
        4,48,-3,49,3,56,-2,57,-1,
        4,48,2,49,-2,56,1,57,1,
        4,24,-3,26,3,48,-2,50,-1,
        4,40,-3,42,3,56,-2,58,-1,
        16,24,9,25,-9,26,-9,27,9,40,6,41,3,42,-6,43,-3,48,6,49,-6,50,3,51,-3,56,4,57,2,58,2,59,1,
        16,24,-6,25,6,26,6,27,-6,40,-3,41,-3,42,3,43,3,48,-4,49,4,50,-2,51,2,56,-2,57,-2,58,-1,59,-1,
        4,24,2,26,-2,48,1,50,1,
        4,40,2,42,-2,56,1,58,1,
        16,24,-6,25,6,26,6,27,-6,40,-4,41,-2,42,4,43,2,48,-3,49,3,50,-3,51,3,56,-2,57,-1,58,-2,59,-1,
        16,24,4,25,-4,26,-4,27,4,40,2,41,2,42,-2,43,-2,48,2,49,-2,50,2,51,-2,56,1,57,1,58,1,59,1,
        4,0,-3,4,3,24,-2,28,-1,
        4,8,-3,12,3,40,-2,44,-1,
        16,0,9,1,-9,4,-9,5,9,8,6,9,3,12,-6,13,-3,24,6,25,-6,28,3,29,-3,40,4,41,2,44,2,45,1,
        16,0,-6,1,6,4,6,5,-6,8,-3,9,-3,12,3,13,3,24,-4,25,4,28,-2,29,2,40,-2,41,-2,44,-1,45,-1,
        4,16,-3,20,3,48,-2,52,-1,
        4,32,-3,36,3,56,-2,60,-1,
        16,16,9,17,-9,20,-9,21,9,32,6,33,3,36,-6,37,-3,48,6,49,-6,52,3,53,-3,56,4,57,2,60,2,61,1,
        16,16,-6,17,6,20,6,21,-6,32,-3,33,-3,36,3,37,3,48,-4,49,4,52,-2,53,2,56,-2,57,-2,60,-1,61,-1,
        16,0,9,2,-9,4,-9,6,9,16,6,18,3,20,-6,22,-3,24,6,26,-6,28,3,30,-3,48,4,50,2,52,2,54,1,
        16,8,9,10,-9,12,-9,14,9,32,6,34,3,36,-6,38,-3,40,6,42,-6,44,3,46,-3,56,4,58,2,60,2,62,1,
        64,0,-27,1,27,2,27,3,-27,4,27,5,-27,6,-27,7,27,8,-18,9,-9,10,18,11,9,12,18,13,9,14,-18,15,-9,16,-18,17,18,18,-9,19,9,20,18,21,-18,22,9,23,-9,24,-18,25,18,26,18,27,-18,28,-9,29,9,30,9,31,-9,32,-12,33,-6,34,-6,35,-3,36,12,37,6,38,6,39,3,40,-12,41,-6,42,12,43,6,44,-6,45,-3,46,6,47,3,48,-12,49,12,50,-6,51,6,52,-6,53,6,54,-3,55,3,56,-8,57,-4,58,-4,59,-2,60,-4,61,-2,62,-2,63,-1,
        64,0,18,1,-18,2,-18,3,18,4,-18,5,18,6,18,7,-18,8,9,9,9,10,-9,11,-9,12,-9,13,-9,14,9,15,9,16,12,17,-12,18,6,19,-6,20,-12,21,12,22,-6,23,6,24,12,25,-12,26,-12,27,12,28,6,29,-6,30,-6,31,6,32,6,33,6,34,3,35,3,36,-6,37,-6,38,-3,39,-3,40,6,41,6,42,-6,43,-6,44,3,45,3,46,-3,47,-3,48,8,49,-8,50,4,51,-4,52,4,53,-4,54,2,55,-2,56,4,57,4,58,2,59,2,60,2,61,2,62,1,63,1,
        16,0,-6,2,6,4,6,6,-6,16,-3,18,-3,20,3,22,3,24,-4,26,4,28,-2,30,2,48,-2,50,-2,52,-1,54,-1,
        16,8,-6,10,6,12,6,14,-6,32,-3,34,-3,36,3,38,3,40,-4,42,4,44,-2,46,2,56,-2,58,-2,60,-1,62,-1,
        64,0,18,1,-18,2,-18,3,18,4,-18,5,18,6,18,7,-18,8,12,9,6,10,-12,11,-6,12,-12,13,-6,14,12,15,6,16,9,17,-9,18,9,19,-9,20,-9,21,9,22,-9,23,9,24,12,25,-12,26,-12,27,12,28,6,29,-6,30,-6,31,6,32,6,33,3,34,6,35,3,36,-6,37,-3,38,-6,39,-3,40,8,41,4,42,-8,43,-4,44,4,45,2,46,-4,47,-2,48,6,49,-6,50,6,51,-6,52,3,53,-3,54,3,55,-3,56,4,57,2,58,4,59,2,60,2,61,1,62,2,63,1,
        64,0,-12,1,12,2,12,3,-12,4,12,5,-12,6,-12,7,12,8,-6,9,-6,10,6,11,6,12,6,13,6,14,-6,15,-6,16,-6,17,6,18,-6,19,6,20,6,21,-6,22,6,23,-6,24,-8,25,8,26,8,27,-8,28,-4,29,4,30,4,31,-4,32,-3,33,-3,34,-3,35,-3,36,3,37,3,38,3,39,3,40,-4,41,-4,42,4,43,4,44,-2,45,-2,46,2,47,2,48,-4,49,4,50,-4,51,4,52,-2,53,2,54,-2,55,2,56,-2,57,-2,58,-2,59,-2,60,-1,61,-1,62,-1,63,-1,
        4,0,2,4,-2,24,1,28,1,
        4,8,2,12,-2,40,1,44,1,
        16,0,-6,1,6,4,6,5,-6,8,-4,9,-2,12,4,13,2,24,-3,25,3,28,-3,29,3,40,-2,41,-1,44,-2,45,-1,
        16,0,4,1,-4,4,-4,5,4,8,2,9,2,12,-2,13,-2,24,2,25,-2,28,2,29,-2,40,1,41,1,44,1,45,1,
        4,16,2,20,-2,48,1,52,1,
        4,32,2,36,-2,56,1,60,1,
        16,16,-6,17,6,20,6,21,-6,32,-4,33,-2,36,4,37,2,48,-3,49,3,52,-3,53,3,56,-2,57,-1,60,-2,61,-1,
        16,16,4,17,-4,20,-4,21,4,32,2,33,2,36,-2,37,-2,48,2,49,-2,52,2,53,-2,56,1,57,1,60,1,61,1,
        16,0,-6,2,6,4,6,6,-6,16,-4,18,-2,20,4,22,2,24,-3,26,3,28,-3,30,3,48,-2,50,-1,52,-2,54,-1,
        16,8,-6,10,6,12,6,14,-6,32,-4,34,-2,36,4,38,2,40,-3,42,3,44,-3,46,3,56,-2,58,-1,60,-2,62,-1,
        64,0,18,1,-18,2,-18,3,18,4,-18,5,18,6,18,7,-18,8,12,9,6,10,-12,11,-6,12,-12,13,-6,14,12,15,6,16,12,17,-12,18,6,19,-6,20,-12,21,12,22,-6,23,6,24,9,25,-9,26,-9,27,9,28,9,29,-9,30,-9,31,9,32,8,33,4,34,4,35,2,36,-8,37,-4,38,-4,39,-2,40,6,41,3,42,-6,43,-3,44,6,45,3,46,-6,47,-3,48,6,49,-6,50,3,51,-3,52,6,53,-6,54,3,55,-3,56,4,57,2,58,2,59,1,60,4,61,2,62,2,63,1,
        64,0,-12,1,12,2,12,3,-12,4,12,5,-12,6,-12,7,12,8,-6,9,-6,10,6,11,6,12,6,13,6,14,-6,15,-6,16,-8,17,8,18,-4,19,4,20,8,21,-8,22,4,23,-4,24,-6,25,6,26,6,27,-6,28,-6,29,6,30,6,31,-6,32,-4,33,-4,34,-2,35,-2,36,4,37,4,38,2,39,2,40,-3,41,-3,42,3,43,3,44,-3,45,-3,46,3,47,3,48,-4,49,4,50,-2,51,2,52,-4,53,4,54,-2,55,2,56,-2,57,-2,58,-1,59,-1,60,-2,61,-2,62,-1,63,-1,
        16,0,4,2,-4,4,-4,6,4,16,2,18,2,20,-2,22,-2,24,2,26,-2,28,2,30,-2,48,1,50,1,52,1,54,1,
        16,8,4,10,-4,12,-4,14,4,32,2,34,2,36,-2,38,-2,40,2,42,-2,44,2,46,-2,56,1,58,1,60,1,62,1,
        64,0,-12,1,12,2,12,3,-12,4,12,5,-12,6,-12,7,12,8,-8,9,-4,10,8,11,4,12,8,13,4,14,-8,15,-4,16,-6,17,6,18,-6,19,6,20,6,21,-6,22,6,23,-6,24,-6,25,6,26,6,27,-6,28,-6,29,6,30,6,31,-6,32,-4,33,-2,34,-4,35,-2,36,4,37,2,38,4,39,2,40,-4,41,-2,42,4,43,2,44,-4,45,-2,46,4,47,2,48,-3,49,3,50,-3,51,3,52,-3,53,3,54,-3,55,3,56,-2,57,-1,58,-2,59,-1,60,-2,61,-1,62,-2,63,-1,
        64,0,8,1,-8,2,-8,3,8,4,-8,5,8,6,8,7,-8,8,4,9,4,10,-4,11,-4,12,-4,13,-4,14,4,15,4,16,4,17,-4,18,4,19,-4,20,-4,21,4,22,-4,23,4,24,4,25,-4,26,-4,27,4,28,4,29,-4,30,-4,31,4,32,2,33,2,34,2,35,2,36,-2,37,-2,38,-2,39,-2,40,2,41,2,42,-2,43,-2,44,2,45,2,46,-2,47,-2,48,2,49,-2,50,2,51,-2,52,2,53,-2,54,2,55,-2,56,1,57,1,58,1,59,1,60,1,61,1,62,1,63,1
    };
    vector<vector<int> > weight(64);
    int index = 0;
    for (int i = 0; i < 64; i++) {
        int numElements = wt[index++];
        for (int j = 0; j < numElements; j++) {
            weight[i].push_back(wt[index++]);
            weight[i].push_back(wt[index++]);
        }
    }
    vector<double> rhs(64);
    c.resize((xsize-1)*(ysize-1)*(zsize-1));
    for (int i = 0; i < xsize-1; i++) {
        for (int j = 0; j < ysize-1; j++) {
            for (int k = 0; k < zsize-1; k++) {
                // Compute the 64 coefficients for patch (i, j, k).

                int nexti = i+1;
                int nextj = j+1;
                int nextk = k+1;
                double deltax = x[nexti]-x[i];
                double deltay = y[nextj]-y[j];
                double deltaz = z[nextk]-z[k];
                double e[] = {values[i+j*xsize+k*xysize], values[nexti+j*xsize+k*xysize], values[i+nextj*xsize+k*xysize], values[nexti+nextj*xsize+k*xysize], values[i+j*xsize+nextk*xysize], values[nexti+j*xsize+nextk*xysize], values[i+nextj*xsize+nextk*xysize], values[nexti+nextj*xsize+nextk*xysize]};
                double e1[] = {d1[i+j*xsize+k*xysize], d1[nexti+j*xsize+k*xysize], d1[i+nextj*xsize+k*xysize], d1[nexti+nextj*xsize+k*xysize], d1[i+j*xsize+nextk*xysize], d1[nexti+j*xsize+nextk*xysize], d1[i+nextj*xsize+nextk*xysize], d1[nexti+nextj*xsize+nextk*xysize]};
                double e2[] = {d2[i+j*xsize+k*xysize], d2[nexti+j*xsize+k*xysize], d2[i+nextj*xsize+k*xysize], d2[nexti+nextj*xsize+k*xysize], d2[i+j*xsize+nextk*xysize], d2[nexti+j*xsize+nextk*xysize], d2[i+nextj*xsize+nextk*xysize], d2[nexti+nextj*xsize+nextk*xysize]};
                double e3[] = {d3[i+j*xsize+k*xysize], d3[nexti+j*xsize+k*xysize], d3[i+nextj*xsize+k*xysize], d3[nexti+nextj*xsize+k*xysize], d3[i+j*xsize+nextk*xysize], d3[nexti+j*xsize+nextk*xysize], d3[i+nextj*xsize+nextk*xysize], d3[nexti+nextj*xsize+nextk*xysize]};
                double e12[] = {d12[i+j*xsize+k*xysize], d12[nexti+j*xsize+k*xysize], d12[i+nextj*xsize+k*xysize], d12[nexti+nextj*xsize+k*xysize], d12[i+j*xsize+nextk*xysize], d12[nexti+j*xsize+nextk*xysize], d12[i+nextj*xsize+nextk*xysize], d12[nexti+nextj*xsize+nextk*xysize]};
                double e13[] = {d13[i+j*xsize+k*xysize], d13[nexti+j*xsize+k*xysize], d13[i+nextj*xsize+k*xysize], d13[nexti+nextj*xsize+k*xysize], d13[i+j*xsize+nextk*xysize], d13[nexti+j*xsize+nextk*xysize], d13[i+nextj*xsize+nextk*xysize], d13[nexti+nextj*xsize+nextk*xysize]};
                double e23[] = {d23[i+j*xsize+k*xysize], d23[nexti+j*xsize+k*xysize], d23[i+nextj*xsize+k*xysize], d23[nexti+nextj*xsize+k*xysize], d23[i+j*xsize+nextk*xysize], d23[nexti+j*xsize+nextk*xysize], d23[i+nextj*xsize+nextk*xysize], d23[nexti+nextj*xsize+nextk*xysize]};
                double e123[] = {d123[i+j*xsize+k*xysize], d123[nexti+j*xsize+k*xysize], d123[i+nextj*xsize+k*xysize], d123[nexti+nextj*xsize+k*xysize], d123[i+j*xsize+nextk*xysize], d123[nexti+j*xsize+nextk*xysize], d123[i+nextj*xsize+nextk*xysize], d123[nexti+nextj*xsize+nextk*xysize]};
                for (int m = 0; m < 8; m++) {
                    rhs[m] = e[m];
                    rhs[m+8] = e1[m]*deltax;
                    rhs[m+16] = e2[m]*deltay;
                    rhs[m+24] = e3[m]*deltaz;
                    rhs[m+32] = e12[m]*deltax*deltay;
                    rhs[m+40] = e13[m]*deltax*deltaz;
                    rhs[m+48] = e23[m]*deltay*deltaz;
                    rhs[m+56] = e123[m]*deltax*deltay*deltaz;
                }
                vector<double>& coeff = c[i+j*(xsize-1)+k*(xsize-1)*(ysize-1)];
                coeff.resize(64);
                for (int m = 0; m < 64; m++) {
                    double sum = 0.0;
                    int numElements = weight[m].size();
                    for (int n = 0; n < numElements; n += 2)
                        sum += weight[m][n+1]*rhs[weight[m][n]];
                    coeff[m] = sum;
                }
            }
        }
    }
}

double SplineFitter::evaluate3DSpline(const vector<double>& x, const vector<double>& y, const vector<double>& z, const vector<double>& values, const vector<vector<double> >& c, double u, double v, double w) {
    int xsize = x.size();
    int ysize = y.size();
    int zsize = z.size();
    if (u < x[0] || u > x[xsize-1] || v < y[0] || v > y[ysize-1] || w < z[0] || w > z[zsize-1])
        throw OpenMMException("evaluate3DSpline: specified point is outside the range defined by the spline");

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
    int lowerz = 0;
    int upperz = zsize-1;
    while (upperz-lowerz > 1) {
        int middle = (upperz+lowerz)/2;
        if (z[middle] > w)
            upperz = middle;
        else
            lowerz = middle;
    }
    double deltax = x[upperx]-x[lowerx];
    double deltay = y[uppery]-y[lowery];
    double deltaz = z[upperz]-z[lowerz];
    double da = (u-x[lowerx])/deltax;
    double db = (v-y[lowery])/deltay;
    double dc = (w-z[lowerz])/deltaz;
    const vector<double>& coeff = c[lowerx+(xsize-1)*lowery+(xsize-1)*(ysize-1)*lowerz];

    // Evaluate the spline to determine the value and gradients.

    double value[] = {0, 0, 0, 0};
    for (int i = 3; i >= 0; i--) {
        for (int j = 0; j < 4; j++) {
            int base = 4*i + 16*j;
            value[j] = db*value[j] + ((coeff[base+3]*da + coeff[base+2])*da + coeff[base+1])*da + coeff[base];
        }
    }
    return value[0] + dc*(value[1] + dc*(value[2] + dc*value[3]));
}

void SplineFitter::evaluate3DSplineDerivatives(const vector<double>& x, const vector<double>& y, const vector<double>& z, const vector<double>& values, const vector<vector<double> >& c, double u, double v, double w, double& dx, double& dy, double& dz) {
    int xsize = x.size();
    int ysize = y.size();
    int zsize = z.size();
    if (u < x[0] || u > x[xsize-1] || v < y[0] || v > y[ysize-1] || w < z[0] || w > z[zsize-1])
        throw OpenMMException("evaluate3DSpline: specified point is outside the range defined by the spline");

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
    int lowerz = 0;
    int upperz = zsize-1;
    while (upperz-lowerz > 1) {
        int middle = (upperz+lowerz)/2;
        if (z[middle] > w)
            upperz = middle;
        else
            lowerz = middle;
    }
    double deltax = x[upperx]-x[lowerx];
    double deltay = y[uppery]-y[lowery];
    double deltaz = z[upperz]-z[lowerz];
    double da = (u-x[lowerx])/deltax;
    double db = (v-y[lowery])/deltay;
    double dc = (w-z[lowerz])/deltaz;
    const vector<double>& coeff = c[lowerx+(xsize-1)*lowery+(xsize-1)*(ysize-1)*lowerz];

    // Evaluate the spline to determine the derivatives.

    double derivx[] = {0, 0, 0, 0};
    double derivy[] = {0, 0, 0, 0};
    double derivz[] = {0, 0, 0, 0};
    for (int i = 3; i >= 0; i--) {
        for (int j = 0; j < 4; j++) {
            int base = 4*i + 16*j;
            derivx[j] = db*derivx[j] + (3.0*coeff[base+3]*da + 2.0*coeff[base+2])*da + coeff[base+1];
            derivz[j] = db*derivz[j] + ((coeff[base+3]*da + coeff[base+2])*da + coeff[base+1])*da + coeff[base];
            base = i + 16*j;
            derivy[j] = da*derivy[j] + (3.0*coeff[base+12]*db + 2.0*coeff[base+8])*db + coeff[base+4];
        }
    }
    dx = derivx[0] + dc*(derivx[1] + dc*(derivx[2] + dc*derivx[3]));
    dy = derivy[0] + dc*(derivy[1] + dc*(derivy[2] + dc*derivy[3]));
    dz = derivz[1] + dc*(2.0*derivz[2] + 3.0*dc*derivz[3]);
    dx /= deltax;
    dy /= deltay;
    dz /= deltaz;
}
