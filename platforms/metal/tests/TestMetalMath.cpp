/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2026 Stanford University and the Authors.           *
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

/**
 * This tests the accuracy of math operations on the Metal platform.
 */

#include "MetalTests.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/CustomIntegrator.h"
#include "openmm/System.h"
#include <cfloat>
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

void testFunctions() {
    const int numParticles = 50;
    System system;
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        positions[i] = Vec3(0.01*pow(-1.5, i), 0.0, 0.0);
    }
    CustomIntegrator integrator(0.1);
    integrator.addPerDofVariable("sqrtx", 0);
    integrator.addPerDofVariable("cosx", 0);
    integrator.addPerDofVariable("sinx", 0);
    integrator.addPerDofVariable("tanx", 0);
    integrator.addPerDofVariable("expx", 0);
    integrator.addPerDofVariable("logx", 0);
    integrator.addComputePerDof("sqrtx", "sqrt(abs(x))");
    integrator.addComputePerDof("cosx", "cos(x)");
    integrator.addComputePerDof("sinx", "sin(x)");
    integrator.addComputePerDof("tanx", "tan(x)");
    integrator.addComputePerDof("expx", "exp(x)");
    integrator.addComputePerDof("logx", "log(abs(x))");
    Context context(system, integrator, platform);
    context.setPositions(positions);
    integrator.step(1);

    // See if it computed the values correctly.

    vector<Vec3> sqrtx, cosx, sinx, tanx, expx, logx;
    integrator.getPerDofVariable(0, sqrtx);
    integrator.getPerDofVariable(1, cosx);
    integrator.getPerDofVariable(2, sinx);
    integrator.getPerDofVariable(3, tanx);
    integrator.getPerDofVariable(4, expx);
    integrator.getPerDofVariable(5, logx);
    bool isSingle = (platform.getPropertyValue(context, "Precision") == "single");
    double tol = (isSingle ? 1e-6 : 1e-12);
    for (int i = 0; i < numParticles; i++) {
        double x = positions[i][0];
        if (isSingle)
            x = (float) x;
        ASSERT_EQUAL_TOL(sqrt(abs(x)), sqrtx[i][0], tol);
        ASSERT_EQUAL_TOL(cos(x), cosx[i][0], tol);
        ASSERT_EQUAL_TOL(sin(x), sinx[i][0], tol);
        ASSERT_EQUAL_TOL(tan(x), tanx[i][0], tol);
        if (exp(x) < FLT_MAX)
            ASSERT_EQUAL_TOL(exp(x), expx[i][0], tol);
        ASSERT_EQUAL_TOL(log(abs(x)), logx[i][0], tol);
    }
}

int main(int argc, char* argv[]) {
    try {
        if (argc > 1)
            platform.setPropertyDefaultValue("Precision", string(argv[1]));
        testFunctions();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
