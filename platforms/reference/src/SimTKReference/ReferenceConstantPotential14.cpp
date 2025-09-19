/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2006-2025 Stanford University and the Authors.      *
 * Authors: Pande Group, Evan Pretti                                          *
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

#include "SimTKOpenMMUtilities.h"
#include "ReferenceConstantPotential14.h"
#include "ReferenceForce.h"

using std::vector;
using namespace OpenMM;

ReferenceConstantPotential14::ReferenceConstantPotential14() : periodic(false) {
}

ReferenceConstantPotential14::~ReferenceConstantPotential14() {
}

void ReferenceConstantPotential14::setPeriodic(OpenMM::Vec3* vectors) {
    periodic = true;
    periodicBoxVectors[0] = vectors[0];
    periodicBoxVectors[1] = vectors[1];
    periodicBoxVectors[2] = vectors[2];
}

void ReferenceConstantPotential14::calculateBondIxn(
    vector<int>& atomIndices, vector<Vec3>& atomCoordinates,
    vector<double>& parameters, vector<Vec3>& forces,
    double* totalEnergy, double* energyParamDerivs
) {
    double deltaR[ReferenceForce::LastDeltaRIndex];

    int atomAIndex = atomIndices[0];
    int atomBIndex = atomIndices[1];
    if (periodic) {
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[atomBIndex], atomCoordinates[atomAIndex], periodicBoxVectors, deltaR);
    }
    else {
        ReferenceForce::getDeltaR(atomCoordinates[atomBIndex], atomCoordinates[atomAIndex], deltaR);
    }

    double inverseR = 1.0 / deltaR[ReferenceForce::RIndex];
    double energy = ONE_4PI_EPS0 * parameters[0] * inverseR;
    double dEdR = energy * inverseR * inverseR;

    for (int ii = 0; ii < 3; ii++) {
        double force = dEdR * deltaR[ii];
        forces[atomAIndex][ii] += force;
        forces[atomBIndex][ii] -= force;
    }

    if (totalEnergy != NULL) {
        *totalEnergy += energy;
    }
}
