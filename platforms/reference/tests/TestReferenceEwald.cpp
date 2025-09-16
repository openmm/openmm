/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2015-2025 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Evan Prettt                                        *
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

#include "ReferenceTests.h"
#include "TestEwald.h"
#include "ReferencePME.h"

void testReferencePmeDerivatives() {
    // Ensures that derivatives reported by the reference PME implementation
    // match those estimated by finite differences.  Checks both forces and
    // charge derivatives for a system with a non-zero net charge.

    static const double DELTA = 1e-5;
    static const double EPSILON = 1e-5;

    Vec3 boxVectors[3] = {
        Vec3(10, 0, 0),
        Vec3(1, 9, 0),
        Vec3(2, 3, 8)
    };

    int numParticles = 10;
    vector<Vec3> positions;
    vector<double> charges;
    vector<int> indices;
    for (int i = 0; i < numParticles; i++) {
        double f = (double)i / numParticles;
        positions.push_back(f * boxVectors[0] + f * f * boxVectors[1] + f * f * f * boxVectors[2]);
        charges.push_back((i % 2 ? -1 : 1) + f * f);
        indices.push_back(i);
    }

    double ewaldAlpha = 1.752;
    int gridSize[3] = {54, 49, 43};

    pme_t pme;
    pme_init(&pme, ewaldAlpha, numParticles, gridSize, 5, 1);

    double dummyEnergy=0;
    vector<Vec3> dummyForces(numParticles);

    vector<Vec3> testForces(numParticles);
    pme_exec(pme, positions, testForces, charges, boxVectors, &dummyEnergy);

    for (int i = 0; i < numParticles; i++) {
        for (int j = 0; j < 3; j++) {
            double referencePosition = positions[i][j];

            double energyLess = 0.0;
            positions[i][j] = referencePosition - DELTA;
            pme_exec(pme, positions, dummyForces, charges, boxVectors, &energyLess);

            double energyMore = 0.0;
            positions[i][j] = referencePosition + DELTA;
            pme_exec(pme, positions, dummyForces, charges, boxVectors, &energyMore);

            positions[i][j] = referencePosition;

            ASSERT_EQUAL_TOL((energyMore - energyLess) / (2 * DELTA), -testForces[i][j], EPSILON);
        }
    }

    vector<double> testDerivatives(numParticles);
    pme_exec_charge_derivatives(pme, positions, testDerivatives, indices, charges, boxVectors);

    for (int i = 0; i < numParticles; i++) {
        double referenceCharge = charges[i];

        double energyLess = 0.0;
        charges[i] = referenceCharge - DELTA;
        pme_exec(pme, positions, dummyForces, charges, boxVectors, &energyLess);

        double energyMore = 0.0;
        charges[i] = referenceCharge + DELTA;
        pme_exec(pme, positions, dummyForces, charges, boxVectors, &energyMore);

        charges[i] = referenceCharge;

        ASSERT_EQUAL_TOL((energyMore - energyLess) / (2 * DELTA), testDerivatives[i], EPSILON);
    }

    pme_destroy(pme);
}

void runPlatformTests() {
    testReferencePmeDerivatives();
}
