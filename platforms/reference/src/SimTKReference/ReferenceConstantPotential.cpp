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

ReferenceConstantPotentialMatrix::ReferenceConstantPotentialMatrix(int numElectrodeParticles) : valid(false), constraintVector(numElectrodeParticles) {
}

void ReferenceConstantPotentialMatrix::invalidate() {
    valid = false;
}

ReferenceConstantPotentialCG::ReferenceConstantPotentialCG(int numElectrodeParticles) : valid(false), precondVector(numElectrodeParticles) {
}

void ReferenceConstantPotentialCG::invalidate() {
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
    // Initializes or updates the precomputed capacitance matrix if this is its
    // first use or electrode parameters have changed since its initialization.

    if (matrix->valid) {
        // If the matrix has already been initialized, ensure that the box size
        // has not changed.  No further action is required.
        if (matrix->boxVectors[0] != boxVectors[0] || matrix->boxVectors[1] != boxVectors[1] || matrix->boxVectors[2] != boxVectors[2]) {
            throw OpenMMException("Box vectors changed since electrode matrix precomputed");
        }
        return;
    }

    matrix->valid = true;
    matrix->boxVectors[0] = boxVectors[0];
    matrix->boxVectors[1] = boxVectors[1];
    matrix->boxVectors[2] = boxVectors[2];

    TNT::Array2D<double> A(numElectrodeParticles, numElectrodeParticles);
    std::vector<double> dUdQ0(numElectrodeParticles);
    std::vector<double> dUdQ(numElectrodeParticles);
    
    pme_t pmeData;
    pme_init(&pmeData, ewaldAlpha, numParticles, gridSize, 5, 1);

    // Get derivatives when all electrode charges are zeroed.
    std::vector<double> q(charges);
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        q[electrodeIndices[ii]] = 0.0;
    }
    getDerivatives(numParticles, numElectrodeParticles, posData, q, exclusions, electrodeIndexMap, electrodeIndices, electrodeParamArray, dUdQ0, pmeData);

    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        int i = electrodeIndices[ii];

        // Get derivatives when one electrode charge is set.
        q[i] = 1.0;
        getDerivatives(numParticles, numElectrodeParticles, posData, q, exclusions, electrodeIndexMap, electrodeIndices, electrodeParamArray, dUdQ, pmeData);
        q[i] = 0.0;

        // Set matrix elements, subtracting zero charge derivatives so that the
        // matrix will end up being the (charge-independent) Hessian.
        for (int jj = 0; jj < ii; jj++) {
            A[ii][jj] = A[jj][ii] = dUdQ[jj] - dUdQ0[jj];
        }
        A[ii][ii] = dUdQ[ii] - dUdQ0[ii];
    }

    pme_destroy(pmeData);

    // Compute Cholesky decomposition representation of the inverse.
    matrix->capacitance = JAMA::Cholesky<double>(A);
    if (!matrix->capacitance.is_spd()) {
        throw OpenMMException("Electrode matrix not positive definite");
    }

    // Precompute the appropriate scaling vector to enforce constant total
    // charge if requested.  The vector is parallel to one obtained by solving
    // Aq = b for all q_i = 1 (ensuring the constrained charges will actually be
    // the correct constrained minimum of the quadratic form for the energy),
    // and is scaled so that adding it to a vector of charges increases the
    // total charge by 1 (making it easy to calculate the necessary offset).
    if (useChargeConstraint) {
        TNT::Array1D<double> constraintVector = matrix->capacitance.solve(TNT::Array1D<double>(numElectrodeParticles, 1.0));
        matrix->constraintVector.assign(static_cast<double*>(constraintVector), static_cast<double*>(constraintVector) + numElectrodeParticles);

        double constraintScaleInv = 0.0;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            constraintScaleInv += matrix->constraintVector[ii];
        }
        double constraintScale = 1.0 / constraintScaleInv;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            matrix->constraintVector[ii] *= constraintScale;
        }
    }
}

void ReferenceConstantPotential::updateCG(
    int numParticles,
    int numElectrodeParticles,
    const std::vector<int>& electrodeIndexMap,
    const std::vector<int>& electrodeIndices,
    const std::vector<std::array<double, 3> >& electrodeParamArray,
    ReferenceConstantPotentialCG* cg
) {
    // Initializes or updates information for a preconditioner for the conjugate
    // gradient method if this is its first use or electrode parameters have
    // changed since its initialization.
    
    const double SQRT_PI = sqrt(PI_M);
    const double SELF_ALPHA_SCALE = 2.0 * ONE_4PI_EPS0 / SQRT_PI;
    const double SELF_ETA_SCALE = SELF_ALPHA_SCALE / sqrt(2.0);
    const double TF_SCALE = 1.0 / EPSILON0;

    // No action is required if the box vectors have not changed.
    if (cg->valid && cg->boxVectors[0] == boxVectors[0] && cg->boxVectors[1] == boxVectors[1] && cg->boxVectors[2] == boxVectors[2]) {
        return;
    }

    cg->valid = true;
    cg->boxVectors[0] = boxVectors[0];
    cg->boxVectors[1] = boxVectors[1];
    cg->boxVectors[2] = boxVectors[2];

    // Perform a reference PME calculation with a single charge at the origin to
    // find the constant offset on the preconditioner diagonal due to the PME
    // calculation.  This will actually vary slightly with position but only due
    // to finite accuracy of the PME splines, so it is fine to assume it will be
    // constant for the preconditioner.
    pme_t pmeData;
    pme_init(&pmeData, ewaldAlpha, 1, gridSize, 5, 1);
    std::vector<Vec3> pmePosData(1);
    std::vector<double> pmeChargeDerivatives(1);
    std::vector<int> pmeElectrodeIndices(1);
    std::vector<double> pmeCharges(1, 1.0);
    pme_exec_charge_derivatives(pmeData, pmePosData, pmeChargeDerivatives, pmeElectrodeIndices, pmeCharges, boxVectors);
    pme_destroy(pmeData);

    // The diagonal has a contribution from reciprocal space, Ewald
    // self-interaction, Gaussian self-interaction, and Thomas-Fermi
    // contributions.
    cg->precondScale = 0.0;
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        cg->precondVector[ii] = 1.0 / (pmeChargeDerivatives[0] + SELF_ETA_SCALE / electrodeParamArray[ii][GaussianWidthIndex]
            - SELF_ALPHA_SCALE * ewaldAlpha + TF_SCALE * electrodeParamArray[ii][ThomasFermiScaleIndex]);
        cg->precondScale += cg->precondVector[ii];
    }
    cg->precondScale = 1.0 / cg->precondScale;
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
    ReferenceConstantPotentialMatrix* matrix,
    ReferenceConstantPotentialCG* cg
) {
    // Solves for charges, then calculates energies and forces.

    pme_t pmeData;
    pme_init(&pmeData, ewaldAlpha, numParticles, gridSize, 5, 1);
    
    solve(numParticles, numElectrodeParticles, posData, charges, exclusions, electrodeIndexMap, electrodeIndices, electrodeParamArray, matrix, cg, pmeData);
    getEnergyForces(numParticles, numElectrodeParticles, posData, forceData, charges, exclusions, electrodeIndexMap, electrodeIndices, electrodeParamArray, energy, pmeData);

    pme_destroy(pmeData);
}

void ReferenceConstantPotential::getCharges(
    int numParticles,
    int numElectrodeParticles,
    const std::vector<Vec3>& posData,
    std::vector<double>& charges,
    const std::vector<std::set<int>>& exclusions,
    const std::vector<int>& electrodeIndexMap,
    const std::vector<int>& electrodeIndices,
    const std::vector<std::array<double, 3> >& electrodeParamArray,
    ReferenceConstantPotentialMatrix* matrix,
    ReferenceConstantPotentialCG* cg
) {
    // Like execute(), but only solves for the charges rather than also
    // computing the energy and forces with the optimized charges.

    pme_t pmeData;
    pme_init(&pmeData, ewaldAlpha, numParticles, gridSize, 5, 1);

    solve(numParticles, numElectrodeParticles, posData, charges, exclusions, electrodeIndexMap, electrodeIndices, electrodeParamArray, matrix, cg, pmeData);

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
    ReferenceConstantPotentialCG* cg,
    pme_t pmeData
) {
    // Solves for electrode charges that minimize the system energy.  It is
    // assumed that either matrix != NULL or cg != NULL.

    if (matrix != NULL) {
        // Precomputed matrix solver.
        
        // Zero electrode charges and get derivatives at zero charge.
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            charges[electrodeIndices[ii]] = 0.0;
        }
        std::vector<double> b(numElectrodeParticles);
        getDerivatives(numParticles, numElectrodeParticles, posData, charges, exclusions, electrodeIndexMap, electrodeIndices, electrodeParamArray, b, pmeData);
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            b[ii] = -b[ii];
        }

        // Solve for electrode charges directly using capacitance matrix and
        // calculated derivatives.
        TNT::Array1D<double> q = matrix->capacitance.solve(TNT::Array1D<double>(numElectrodeParticles, b.data()));
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            charges[electrodeIndices[ii]] = q[ii];
        }

        // Enforce total charge constraint if requested.
        if (useChargeConstraint) {
            double chargeOffset = chargeTarget;
            for (int ii = 0; ii < numElectrodeParticles; ii++) {
                chargeOffset -= q[ii];
            }
            for (int ii = 0; ii < numElectrodeParticles; ii++) {
                charges[electrodeIndices[ii]] += chargeOffset * matrix->constraintVector[ii];
            }
        }

    }
    else {
        // Conjugate gradient solver.

        std::vector<double> q(numElectrodeParticles);
        std::vector<double> grad(numElectrodeParticles);
        std::vector<double> projGrad(numElectrodeParticles);
        std::vector<double> qStep(numElectrodeParticles);
        std::vector<double> gradStep(numElectrodeParticles);
        std::vector<double> grad0(numElectrodeParticles);
        double gamma, gammaNext, offset, error;

        // Ensure that initial guess charges satisfy the constraint.
        if (useChargeConstraint) {
            double offset = chargeTarget;
            for (int ii = 0; ii < numElectrodeParticles; ii++) {
                offset -= charges[electrodeIndices[ii]];
            }
            offset /= numElectrodeParticles;
            for (int ii = 0; ii < numElectrodeParticles; ii++) {
                charges[electrodeIndices[ii]] += offset;
            }
        }

        // Evaluate the initial gradient Aq - b.
        getDerivatives(numParticles, numElectrodeParticles, posData, charges, exclusions, electrodeIndexMap, electrodeIndices, electrodeParamArray, grad, pmeData);

        // Project the initial gradient with preconditioning.
        offset = 0.0;
        if (useChargeConstraint) {
            for (int ii = 0; ii < numElectrodeParticles; ii++) {
                offset += cg->precondVector[ii] * grad[ii];
            }
            offset *= cg->precondScale;
        }
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            projGrad[ii] = cg->precondVector[ii] * (grad[ii] - offset);
        }

        // Save the current charges, then evaluate the gradient with zero
        // charges (-b) so that we can later compute Aq as (Aq - b) - (-b).
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            int i = electrodeIndices[ii];
            q[ii] = charges[i];
            charges[i] = 0.0;
        }
        getDerivatives(numParticles, numElectrodeParticles, posData, charges, exclusions, electrodeIndexMap, electrodeIndices, electrodeParamArray, grad0, pmeData);

        // Initialize parameters for conjugate gradient iterations.
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            qStep[ii] = -projGrad[ii];
        }

        gamma = 0.0;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            gamma += grad[ii] * projGrad[ii];
        }

        // Perform conjugate gradient iterations.
        bool converged = false;
        for (int iter = 0; iter < numElectrodeParticles; iter++) {
            // Evaluate the matrix-vector product A qStep.
            for (int ii = 0; ii < numElectrodeParticles; ii++) {
                charges[electrodeIndices[ii]] = qStep[ii];
            }
            getDerivatives(numParticles, numElectrodeParticles, posData, charges, exclusions, electrodeIndexMap, electrodeIndices, electrodeParamArray, gradStep, pmeData);
            for (int ii = 0; ii < numElectrodeParticles; ii++) {
                gradStep[ii] -= grad0[ii];
            }
            
            // Evaluate the conjugate gradient parameter alpha and update the
            // charge and gradient vectors.
            double alphaNumerator = 0.0;
            double alphaDenominator = 0.0;
            for (int ii = 0; ii < numElectrodeParticles; ii++) {
                alphaNumerator -= grad[ii] * qStep[ii];
                alphaDenominator += qStep[ii] * gradStep[ii];
            }
            double alpha = alphaNumerator / alphaDenominator;

            for (int ii = 0; ii < numElectrodeParticles; ii++) {
                q[ii] += alpha * qStep[ii];
                grad[ii] += alpha * gradStep[ii];
            }

            // Check for convergence.  We have to check the step size (see
            // Shariff, Computers Math. Applic. 30, 11, 25-37, 1995); the
            // projected gradient with preconditioning may not give meaningful
            // convergence information.
            error = 0.0;
            for (int ii = 0; ii < numElectrodeParticles; ii++) {
                error += qStep[ii] * qStep[ii];
            }
            error *= alpha * alpha;
            if (error <= cgErrorTol * cgErrorTol * numElectrodeParticles) {
                converged = true;
                break;
            }

            // Project the current gradient with preconditioning.
            offset = 0.0;
            if (useChargeConstraint) {
                for (int ii = 0; ii < numElectrodeParticles; ii++) {
                    offset += cg->precondVector[ii] * grad[ii];
                }
                offset *= cg->precondScale;
            }
            for (int ii = 0; ii < numElectrodeParticles; ii++) {
                projGrad[ii] = cg->precondVector[ii] * (grad[ii] - offset);
            }

            // Evaluate the conjugate gradient parameter beta and update
            // parameters for the following iteration.
            gammaNext = 0.0;
            for (int ii = 0; ii < numElectrodeParticles; ii++) {
                gammaNext += grad[ii] * projGrad[ii];
            }
            double beta = gammaNext / gamma;
            gamma = gammaNext;

            for (int ii = 0; ii < numElectrodeParticles; ii++) {
                qStep[ii] = beta * qStep[ii] - projGrad[ii];
            }

            if (useChargeConstraint) {
                // Analytically, when using the total charge constraint, qStep
                // should always sum to 0 so that the constraint will never be
                // violated as long as it is initially satisfied.  However,
                // numerically, the algorithm is unstable with respect to a
                // catastrophic buildup of constraint error!  We thus explicitly
                // subtract off any accumulated drift after each iteration.
                offset = 0.0;
                for (int ii = 0; ii < numElectrodeParticles; ii++) {
                    offset += qStep[ii];
                }
                offset /= numElectrodeParticles;
                for (int ii = 0; ii < numElectrodeParticles; ii++) {
                    qStep[ii] -= offset;
                }
            }
        }

        if (!converged) {
            throw OpenMMException("Constant potential conjugate gradient iterations not converged");
        }

        // Store the final charges.
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            charges[electrodeIndices[ii]] = q[ii];
        }
    }
}

void ReferenceConstantPotential::getEnergyForces(
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
    // Computes energies and forces given optimized charges.

    const double SQRT_PI = sqrt(PI_M);
    const double TWO_OVER_SQRT_PI = 2.0 / SQRT_PI;
    const double SELF_ALPHA_SCALE = ONE_4PI_EPS0 / SQRT_PI;
    const double SELF_ETA_SCALE = SELF_ALPHA_SCALE / sqrt(2.0);
    const double TF_SCALE = 1.0 / (2.0 * EPSILON0);

    double energyAccum = 0.0;

    // Direct space.
    for (auto& pair : *neighborList) {
        int i = pair.first;
        int j = pair.second;
        int ii = electrodeIndexMap[i];
        int jj = electrodeIndexMap[j];

        double iWidth = ii == -1 ? 0.0 : electrodeParamArray[ii][GaussianWidthIndex];
        double jWidth = jj == -1 ? 0.0 : electrodeParamArray[jj][GaussianWidthIndex];
        double width = sqrt(iWidth * iWidth + jWidth * jWidth);

        double deltaR[ReferenceForce::LastDeltaRIndex];
        ReferenceForce::getDeltaRPeriodic(posData[i], posData[j], boxVectors, deltaR);
        double r = deltaR[ReferenceForce::RIndex];
        double inverseR = 1.0 / r;

        double alphaR = ewaldAlpha * r;
        double erfcAlphaR = erfc(alphaR);

        double etaR = 0.0;
        double erfcEtaR = 0.0;
        double expEtaRTerm = 0.0;
        if (width != 0.0) {
            etaR = r / width;
            erfcEtaR = erfc(etaR);
            expEtaRTerm = etaR * exp(-etaR * etaR);
        }

        double qqFactor = ONE_4PI_EPS0 * charges[i] * charges[j];
        double forceFactor = qqFactor * inverseR * inverseR * inverseR * (erfcAlphaR - erfcEtaR + TWO_OVER_SQRT_PI * (alphaR * exp(-alphaR * alphaR) - expEtaRTerm));
        for (int k = 0; k < 3; k++) {
            double force = forceFactor * deltaR[k];
            forceData[i][k] -= force;
            forceData[j][k] += force;
        }

        energyAccum += qqFactor * inverseR * (erfcAlphaR - erfcEtaR);
    }

    // Exceptions.
    for (int i = 0; i < numParticles; i++) {
        for (int j : exclusions[i]) {
            if (j <= i) {
                continue;
            }

            double deltaR[ReferenceForce::LastDeltaRIndex];
            ReferenceForce::getDeltaRPeriodic(posData[i], posData[j], boxVectors, deltaR);
            double r = deltaR[ReferenceForce::RIndex];
            double inverseR = 1.0 / r;

            double qqFactor = ONE_4PI_EPS0 * charges[i] * charges[j];
            double alphaR = ewaldAlpha * r;

            if (alphaR > 1e-6) {
                double erfAlphaR = erf(alphaR);

                double forceFactor = qqFactor * inverseR * inverseR * inverseR * (erfAlphaR + TWO_OVER_SQRT_PI * alphaR * exp(-alphaR * alphaR));
                for (int k = 0; k < 3; k++) {
                    double force = forceFactor * deltaR[k];
                    forceData[i][k] += force;
                    forceData[j][k] -= force;
                }

                energyAccum -= qqFactor * inverseR * erfAlphaR;
            }
            else {
                energyAccum -= qqFactor * TWO_OVER_SQRT_PI * ewaldAlpha;
            }
        }
    }
    
    // Reciprocal space.
    double pmeEnergy = 0.0;
    pme_exec(pmeData, posData, forceData, charges, boxVectors, &pmeEnergy);
    energyAccum += pmeEnergy;

    // Ewald self-interaction and external field contributions (all particles).
    for (int i = 0; i < numParticles; i++) {
        double q = charges[i];
        energyAccum -= q * (SELF_ALPHA_SCALE * q * ewaldAlpha + posData[i].dot(externalField));
        forceData[i] += q * externalField;
    }

    // Gaussian self-interaction, potential, and Thomas-Fermi contributions
    // (electrode particles only).
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        int i = electrodeIndices[ii];
        double q = charges[i];
        energyAccum += q * (q * (SELF_ETA_SCALE / electrodeParamArray[ii][GaussianWidthIndex] + TF_SCALE * electrodeParamArray[ii][ThomasFermiScaleIndex]) - electrodeParamArray[ii][PotentialIndex]);
    }

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
    // Computes energy derivatives with respect to electrode charges.

    const double SQRT_PI = sqrt(PI_M);
    const double SELF_ALPHA_SCALE = 2.0 * ONE_4PI_EPS0 / SQRT_PI;
    const double SELF_ETA_SCALE = SELF_ALPHA_SCALE / sqrt(2.0);
    const double TF_SCALE = 1.0 / EPSILON0;

    // chargeDerivatives is to be overwritten by this function, not incremented,
    // so zero all derivatives initially.
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        chargeDerivatives[ii] = 0.0;
    }

    // Direct space (both particles in each exception will be non-electrode, so
    // we do not need to loop over exceptions when computing derivatives with
    // respect to electrode charges).
    for (auto& pair : *neighborList) {
        int i = pair.first;
        int j = pair.second;
        int ii = electrodeIndexMap[i];
        int jj = electrodeIndexMap[j];

        if (ii == -1 && jj == -1) {
            continue;
        }

        double iWidth = ii == -1 ? 0.0 : electrodeParamArray[ii][GaussianWidthIndex];
        double jWidth = jj == -1 ? 0.0 : electrodeParamArray[jj][GaussianWidthIndex];
        double width = sqrt(iWidth * iWidth + jWidth * jWidth);

        double deltaR[ReferenceForce::LastDeltaRIndex];
        ReferenceForce::getDeltaRPeriodic(posData[i], posData[j], boxVectors, deltaR);
        double r = deltaR[ReferenceForce::RIndex];

        double erfcEtaR = width == 0.0 ? 0.0 : erfc(r / width);
        double factor = ONE_4PI_EPS0 * (erfc(ewaldAlpha * r) - erfcEtaR) / r;

        if (ii != -1) {
            chargeDerivatives[ii] += charges[j] * factor;
        }
        if (jj != -1) {
            chargeDerivatives[jj] += charges[i] * factor;
        }
    }

    // Reciprocal space.
    pme_exec_charge_derivatives(pmeData, posData, chargeDerivatives, electrodeIndices, charges, boxVectors);

    // Ewald self-interaction, Gaussian self-interaction, potential, external
    // field, and Thomas-Fermi contributions.
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        int i = electrodeIndices[ii];
        double q = charges[i];
        chargeDerivatives[ii] += q * (SELF_ETA_SCALE / electrodeParamArray[ii][GaussianWidthIndex] - SELF_ALPHA_SCALE * ewaldAlpha + TF_SCALE * electrodeParamArray[ii][ThomasFermiScaleIndex])
            - electrodeParamArray[ii][PotentialIndex] - posData[i].dot(externalField);
    }
}
