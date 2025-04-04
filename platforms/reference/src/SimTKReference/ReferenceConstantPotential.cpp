/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
 * Authors: Evan Pretti                                                       *
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

#include "ReferenceConstantPotential.h"
#include "ReferenceForce.h"
#include "ReferencePME.h"
#include "SimTKOpenMMUtilities.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/MSVC_erfc.h"

using namespace OpenMM;

ReferenceConstantPotentialMatrix::ReferenceConstantPotentialMatrix() : valid(false) {
}

void ReferenceConstantPotentialMatrix::invalidate() {
    valid = false;
}

ReferenceConstantPotential::ReferenceConstantPotential(
    double nonbondedCutoff,
    const NeighborList* neighborList,
    const Vec3* boxVectors,
    bool exceptionsArePeriodic,
    double ewaldAlpha,
    const int* gridSize,
    double cgErrorTol,
    bool useChargeConstraint,
    double chargeTarget,
    Vec3 externalField
) : nonbondedCutoff(nonbondedCutoff),
    neighborList(neighborList),
    exceptionsArePeriodic(exceptionsArePeriodic),
    ewaldAlpha(ewaldAlpha),
    cgErrorTol(cgErrorTol),
    useChargeConstraint(useChargeConstraint),
    chargeTarget(chargeTarget),
    externalField(externalField)
{
    this->boxVectors[0] = boxVectors[0];
    this->boxVectors[1] = boxVectors[1];
    this->boxVectors[2] = boxVectors[2];
    this->gridSize[0] = gridSize[0];
    this->gridSize[1] = gridSize[1];
    this->gridSize[2] = gridSize[2];
}

void ReferenceConstantPotential::updateMatrix(
    int numParticles,
    int numElectrodeParticles,
    const std::vector<Vec3>& posData,
    const std::vector<double>& charges,
    const std::vector<std::set<int>>& exclusions,
    const std::vector<int>& electrodeIndexMap,
    const std::vector<int>& electrodeIndices,
    const std::vector<std::array<double, 3> >& electrodeParamArray,
    ReferenceConstantPotentialMatrix* matrix
) {
    if (matrix->valid) {
        if (matrix->boxVectors[0] != boxVectors[0] || matrix->boxVectors[1] != boxVectors[1] || matrix->boxVectors[2] != boxVectors[2]) {
            throw OpenMMException("Box vectors changed since electrode matrix precomputed (use conjugate gradient instead)");
        }
        return;
    }

    matrix->valid = true;
    matrix->boxVectors[0] = boxVectors[0];
    matrix->boxVectors[1] = boxVectors[1];
    matrix->boxVectors[2] = boxVectors[2];

    TNT::Array2D<double> A(numElectrodeParticles, numElectrodeParticles);
    
    pme_t pmeData;
    pme_init(&pmeData, ewaldAlpha, numParticles, gridSize, 5, 1);

    // Get derivatives when all electrode charges are zeroed.
    std::vector<double> q(charges);
    for (int iElectrode = 0; iElectrode < numElectrodeParticles; iElectrode++) {
        q[electrodeIndices[iElectrode]] = 0.0;
    }
    std::vector<double> dUdQ0(numElectrodeParticles, 0.0);
    getDerivatives(numParticles, numElectrodeParticles, posData, q, exclusions, electrodeIndexMap, electrodeIndices, electrodeParamArray, dUdQ0, pmeData);

    for (int iElectrode = 0; iElectrode < numElectrodeParticles; iElectrode++) {
        int i = electrodeIndices[iElectrode];

        // Get derivatives when one electrode charge is set.
        q[i] = 1.0;
        std::vector<double> dUdQ(numElectrodeParticles, 0.0);
        getDerivatives(numParticles, numElectrodeParticles, posData, q, exclusions, electrodeIndexMap, electrodeIndices, electrodeParamArray, dUdQ, pmeData);
        q[i] = 0.0;

        // Set matrix elements (subtract zero charge derivatives so that the
        // matrix elements will not depend on non-electrode particle positions).
        for (int jElectrode = 0; jElectrode < numElectrodeParticles; jElectrode++) {
            A[iElectrode][jElectrode] = dUdQ[jElectrode] - dUdQ0[jElectrode];
        }
    }

    // Compute Cholesky decomposition representation of the inverse.
    matrix->Ainv = JAMA::Cholesky<double>(A);
    if (!matrix->Ainv.is_spd()) {
        throw OpenMMException("Electrode matrix not symmetric positive definite");
    }

    pme_destroy(pmeData);
}

void ReferenceConstantPotential::execute(
    int numParticles,
    int numElectrodeParticles,
    const std::vector<Vec3>& posData,
    std::vector<Vec3>& forceData,
    std::vector<double>& charges,
    const std::vector<std::set<int>>& exclusions,
    const std::vector<int>& electrodeIndexMap,
    const std::vector<int>& electrodeIndices,
    const std::vector<std::array<double, 3> >& electrodeParamArray,
    double* energy,
    ReferenceConstantPotentialMatrix* matrix
) {
    pme_t pmeData;
    pme_init(&pmeData, ewaldAlpha, numParticles, gridSize, 5, 1);
    
    solve(numParticles, numElectrodeParticles, posData, charges, exclusions, electrodeIndexMap, electrodeIndices, electrodeParamArray, matrix, pmeData);
    getEnergiesForces(numParticles, numElectrodeParticles, posData, forceData, charges, exclusions, electrodeIndexMap, electrodeIndices, electrodeParamArray, energy, pmeData);

    pme_destroy(pmeData);
}

void ReferenceConstantPotential::solve(
    int numParticles,
    int numElectrodeParticles,
    const std::vector<Vec3>& posData,
    std::vector<double>& charges,
    const std::vector<std::set<int>>& exclusions,
    const std::vector<int>& electrodeIndexMap,
    const std::vector<int>& electrodeIndices,
    const std::vector<std::array<double, 3> >& electrodeParamArray,
    ReferenceConstantPotentialMatrix* matrix
) {
    pme_t pmeData;
    pme_init(&pmeData, ewaldAlpha, numParticles, gridSize, 5, 1);

    solve(numParticles, numElectrodeParticles, posData, charges, exclusions, electrodeIndexMap, electrodeIndices, electrodeParamArray, matrix, pmeData);

    pme_destroy(pmeData);
}

void ReferenceConstantPotential::solve(
    int numParticles,
    int numElectrodeParticles,
    const std::vector<Vec3>& posData,
    std::vector<double>& charges,
    const std::vector<std::set<int>>& exclusions,
    const std::vector<int>& electrodeIndexMap,
    const std::vector<int>& electrodeIndices,
    const std::vector<std::array<double, 3> >& electrodeParamArray,
    ReferenceConstantPotentialMatrix* matrix,
    pme_t pmeData
) {
    if(matrix != NULL) {
        // Precomputed matrix solver.

        std::vector<double> q(charges);
        for (int iElectrode = 0; iElectrode < numElectrodeParticles; iElectrode++) {
            q[electrodeIndices[iElectrode]] = 0.0;
        }
        std::vector<double> dUdQ0(numElectrodeParticles, 0.0);
        getDerivatives(numParticles, numElectrodeParticles, posData, q, exclusions, electrodeIndexMap, electrodeIndices, electrodeParamArray, dUdQ0, pmeData);

        std::vector<double> dUdQ(numElectrodeParticles, 0.0);
        getDerivatives(numParticles, numElectrodeParticles, posData, charges, exclusions, electrodeIndexMap, electrodeIndices, electrodeParamArray, dUdQ, pmeData);

        TNT::Array1D<double> b(numElectrodeParticles);
        for (int iElectrode = 0; iElectrode < numElectrodeParticles; iElectrode++) {
            b[iElectrode] = dUdQ[iElectrode] - dUdQ0[iElectrode];
        }

        TNT::Array1D<double> qOpt = matrix->Ainv.solve(b);
        for (int iElectrode = 0; iElectrode < numElectrodeParticles; iElectrode++) {
            charges[electrodeIndices[iElectrode]] = qOpt[iElectrode];
        }
    } else {
        // Conjugate gradient solver.
        
        // TODO
        throw OpenMMException("CG solver not yet supported!");
    }
}

void ReferenceConstantPotential::getEnergiesForces(
    int numParticles,
    int numElectrodeParticles,
    const std::vector<Vec3>& posData,
    std::vector<Vec3>& forceData,
    std::vector<double>& charges,
    const std::vector<std::set<int>>& exclusions,
    const std::vector<int>& electrodeIndexMap,
    const std::vector<int>& electrodeIndices,
    const std::vector<std::array<double, 3> >& electrodeParamArray,
    double* energy,
    pme_t pmeData
) {
    const double SQRT_PI = sqrt(PI_M);
    const double TWO_BY_SQRT_PI = 2.0 / SQRT_PI;
    const double SELF_ALPHA_SCALE = ONE_4PI_EPS0 / SQRT_PI;
    const double SELF_ETA_SCALE = SELF_ALPHA_SCALE / sqrt(2.0);

    // Compute total energy and forces with final charge vector.
    double energyAccum = 0.0;

    // 1. Direct space.
    for (auto& pair : *neighborList) {
        int i = pair.first;
        int j = pair.second;
        int iElectrode = electrodeIndexMap[i];
        int jElectrode = electrodeIndexMap[j];

        double iWidth = iElectrode == -1 ? 0.0 : electrodeParamArray[iElectrode][GaussianWidthIndex];
        double jWidth = jElectrode == -1 ? 0.0 : electrodeParamArray[jElectrode][GaussianWidthIndex];
        double width = sqrt(iWidth * iWidth + jWidth * jWidth);
        double eta = 1.0 / width;

        double deltaR[ReferenceForce::LastDeltaRIndex];
        ReferenceForce::getDeltaRPeriodic(posData[i], posData[j], boxVectors, deltaR);
        double r = deltaR[ReferenceForce::RIndex];
        double inverseR = 1.0 / r;

        double alphaR = ewaldAlpha * r;
        double etaR = eta * r;
        double erfcAlphaR = erfc(alphaR);
        double erfcEtaR = erfc(etaR);

        double qqFactor = ONE_4PI_EPS0 * charges[i] * charges[j];
        double forceFactor = qqFactor * inverseR * inverseR * inverseR * (erfcAlphaR - erfcEtaR + TWO_BY_SQRT_PI * (alphaR * exp(-alphaR * alphaR) + etaR * exp(-etaR * etaR)));
        for (int k = 0; k < 3; k++) {
            double force = forceFactor * deltaR[k];
            forceData[i][k] -= force;
            forceData[j][k] += force;
        }

        energyAccum += qqFactor * inverseR * (erfcAlphaR - erfcEtaR);
    }

    // 2. Exceptions.
    for (int i = 0; i < numParticles; i++) {
        for (int j : exclusions[i]) {
            if (j <= i) {
                continue;
            }

            int iElectrode = electrodeIndexMap[i];
            int jElectrode = electrodeIndexMap[j];

            double iWidth = iElectrode == -1 ? 0.0 : electrodeParamArray[iElectrode][GaussianWidthIndex];
            double jWidth = jElectrode == -1 ? 0.0 : electrodeParamArray[jElectrode][GaussianWidthIndex];
            double width = sqrt(iWidth * iWidth + jWidth * jWidth);
            double eta = 1.0 / width;

            double deltaR[ReferenceForce::LastDeltaRIndex];
            ReferenceForce::getDeltaRPeriodic(posData[i], posData[j], boxVectors, deltaR);
            double r = deltaR[ReferenceForce::RIndex];
            double inverseR = 1.0 / r;

            double alphaR = ewaldAlpha * r;
            double etaR = eta * r;
            double erfAlphaR = erf(alphaR);
            double erfcEtaR = erfc(etaR);

            double qqFactor = ONE_4PI_EPS0 * charges[i] * charges[j];
            double forceFactor = qqFactor * inverseR * inverseR * inverseR * (erfAlphaR + erfcEtaR + TWO_BY_SQRT_PI * (alphaR * exp(-alphaR * alphaR) + etaR * exp(-etaR * etaR)));
            for (int k = 0; k < 3; k++) {
                double force = forceFactor * deltaR[k];
                forceData[i][k] += force;
                forceData[j][k] -= force;
            }

            energyAccum -= qqFactor * inverseR * (erfAlphaR + erfcEtaR);
        }
    }
    
    // 3. Reciprocal space.
    double pmeEnergy = 0.0;
    pme_exec(pmeData, posData, forceData, charges, boxVectors, &pmeEnergy);
    energyAccum += pmeEnergy;

    // 4. Self-interaction.
    for (int i = 0; i < numParticles; i++) {
        energyAccum -= SELF_ALPHA_SCALE * charges[i] * charges[i] * ewaldAlpha;
    }
    for (int iElectrode = 0; iElectrode < numElectrodeParticles; iElectrode++) {
        int i = electrodeIndices[iElectrode];
        energyAccum += SELF_ETA_SCALE * charges[i] * charges[i] * electrodeParamArray[iElectrode][GaussianWidthIndex];
    }

    // 5. Potential and field.
    // TODO.

    if (energy) {
        *energy += energyAccum;
    }
}

void ReferenceConstantPotential::getDerivatives(
    int numParticles,
    int numElectrodeParticles,
    const std::vector<Vec3>& posData,
    std::vector<double>& charges,
    const std::vector<std::set<int>>& exclusions,
    const std::vector<int>& electrodeIndexMap,
    const std::vector<int>& electrodeIndices,
    const std::vector<std::array<double, 3> >& electrodeParamArray,
    std::vector<double>& chargeDerivatives,
    pme_t pmeData
) {
    const double SQRT_PI = sqrt(PI_M);
    const double SELF_ALPHA_SCALE = ONE_4PI_EPS0 / SQRT_PI;
    const double SELF_ETA_SCALE = SELF_ALPHA_SCALE / sqrt(2.0);

    for (int iElectrode = 0; iElectrode < numElectrodeParticles; iElectrode++) {
        chargeDerivatives[iElectrode] = 0.0;
    }

    // 1. Direct space.
    for (auto& pair : *neighborList) {
        int i = pair.first;
        int j = pair.second;
        int iElectrode = electrodeIndexMap[i];
        int jElectrode = electrodeIndexMap[j];

        if (iElectrode == -1 && jElectrode == -1) {
            return;
        }

        double iWidth = iElectrode == -1 ? 0.0 : electrodeParamArray[iElectrode][GaussianWidthIndex];
        double jWidth = jElectrode == -1 ? 0.0 : electrodeParamArray[jElectrode][GaussianWidthIndex];
        double width = sqrt(iWidth * iWidth + jWidth * jWidth);
        double eta = 1.0 / width;

        double deltaR[ReferenceForce::LastDeltaRIndex];
        ReferenceForce::getDeltaRPeriodic(posData[i], posData[j], boxVectors, deltaR);
        double r = deltaR[ReferenceForce::RIndex];
        double inverseR = 1.0 / r;

        double alphaR = ewaldAlpha * r;
        double etaR = eta * r;
        double erfcAlphaR = erfc(alphaR);
        double erfcEtaR = erfc(etaR);

        double factor = ONE_4PI_EPS0 * inverseR * (erfcAlphaR - erfcEtaR);
        if (iElectrode != -1) {
            chargeDerivatives[iElectrode] += factor;
        }
        if (jElectrode != -1) {
            chargeDerivatives[jElectrode] += factor;
        }
    }

    // 2. Exceptions.
    for (int i = 0; i < numParticles; i++) {
        for (int j : exclusions[i]) {
            if (j <= i) {
                continue;
            }

            int iElectrode = electrodeIndexMap[i];
            int jElectrode = electrodeIndexMap[j];

            if (iElectrode == -1 && jElectrode == -1) {
                return;
            }

            double iWidth = iElectrode == -1 ? 0.0 : electrodeParamArray[iElectrode][GaussianWidthIndex];
            double jWidth = jElectrode == -1 ? 0.0 : electrodeParamArray[jElectrode][GaussianWidthIndex];
            double width = sqrt(iWidth * iWidth + jWidth * jWidth);
            double eta = 1.0 / width;

            double deltaR[ReferenceForce::LastDeltaRIndex];
            ReferenceForce::getDeltaRPeriodic(posData[i], posData[j], boxVectors, deltaR);
            double r = deltaR[ReferenceForce::RIndex];
            double inverseR = 1.0 / r;

            double alphaR = ewaldAlpha * r;
            double etaR = eta * r;
            double erfAlphaR = erf(alphaR);
            double erfcEtaR = erfc(etaR);

            double factor = ONE_4PI_EPS0 * inverseR * (erfAlphaR + erfcEtaR);
            if (iElectrode != -1) {
                chargeDerivatives[iElectrode] -= factor;
            }
            if (jElectrode != -1) {
                chargeDerivatives[jElectrode] -= factor;
            }
        }
    }
    
    // 3. Reciprocal space.
    pme_exec_charge_derivatives(pmeData, posData, chargeDerivatives, electrodeIndices, charges, boxVectors);

    // 4. Self-interaction.
    for (int iElectrode = 0; iElectrode < numElectrodeParticles; iElectrode++) {
        chargeDerivatives[iElectrode] += 2.0 * charges[electrodeIndices[iElectrode]] * (SELF_ETA_SCALE * electrodeParamArray[iElectrode][GaussianWidthIndex] - SELF_ALPHA_SCALE * ewaldAlpha);
    }

    // 5. Potential and field.
    // TODO.
}
