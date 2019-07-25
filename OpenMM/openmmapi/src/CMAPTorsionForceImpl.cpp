/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2015 Stanford University and the Authors.      *
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
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/CMAPTorsionForceImpl.h"
#include "openmm/internal/SplineFitter.h"
#include "openmm/kernels.h"
#include <cmath>

using namespace OpenMM;
using namespace std;

CMAPTorsionForceImpl::CMAPTorsionForceImpl(const CMAPTorsionForce& owner) : owner(owner) {
}

CMAPTorsionForceImpl::~CMAPTorsionForceImpl() {
}

void CMAPTorsionForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcCMAPTorsionForceKernel::Name(), context);
    kernel.getAs<CalcCMAPTorsionForceKernel>().initialize(context.getSystem(), owner);
}

double CMAPTorsionForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0)
        return kernel.getAs<CalcCMAPTorsionForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

vector<string> CMAPTorsionForceImpl::getKernelNames() {
    vector<string> names;
    names.push_back(CalcCMAPTorsionForceKernel::Name());
    return names;
}

void CMAPTorsionForceImpl::calcMapDerivatives(int size, const vector<double>& energy, vector<vector<double> >& c) {
    vector<double> d1(size*size), d2(size*size), d12(size*size);
    vector<double> x(size+1), y(size+1), deriv(size+1);
    for (int i = 0; i < size+1; i++)
        x[i] = i*2*M_PI/size;

    // Compute derivatives with respect to the first angle.

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++)
            y[j] = energy[j+size*i];
        y[size] = energy[size*i];
        SplineFitter::createPeriodicSpline(x, y, deriv);
        for (int j = 0; j < size; j++) {
            d1[j+size*i] = SplineFitter::evaluateSplineDerivative(x, y, deriv, x[j]);
        }
    }

    // Compute derivatives with respect to the second angle.

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++)
            y[j] = energy[i+size*j];
        y[size] = energy[i];
        SplineFitter::createPeriodicSpline(x, y, deriv);
        for (int j = 0; j < size; j++)
            d2[i+size*j] = SplineFitter::evaluateSplineDerivative(x, y, deriv, x[j]);
    }

    // Compute cross derivatives.

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++)
            y[j] = d2[j+size*i];
        y[size] = d2[size*i];
        SplineFitter::createPeriodicSpline(x, y, deriv);
        for (int j = 0; j < size; j++)
            d12[j+size*i] = SplineFitter::evaluateSplineDerivative(x, y, deriv, x[j]);
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
    double delta = 2*M_PI/size;
    c.resize(size*size);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            // Compute the 16 coefficients for patch (i, j).

            int nexti = (i+1)%size;
            int nextj = (j+1)%size;
            double e[] = {energy[i+j*size], energy[nexti+j*size], energy[nexti+nextj*size], energy[i+nextj*size]};
            double e1[] = {d1[i+j*size], d1[nexti+j*size], d1[nexti+nextj*size], d1[i+nextj*size]};
            double e2[] = {d2[i+j*size], d2[nexti+j*size], d2[nexti+nextj*size], d2[i+nextj*size]};
            double e12[] = {d12[i+j*size], d12[nexti+j*size], d12[nexti+nextj*size], d12[i+nextj*size]};

            for (int k = 0; k < 4; k++) {
                rhs[k] = e[k];
                rhs[k+4] = e1[k]*delta;
                rhs[k+8] = e2[k]*delta;
                rhs[k+12] = e12[k]*delta*delta;
            }
            vector<double>& coeff = c[i+j*size];
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

void CMAPTorsionForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcCMAPTorsionForceKernel>().copyParametersToContext(context, owner);
    context.systemChanged();
}
