/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
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

ReferenceConstantPotentialSolver::ReferenceConstantPotentialSolver() : valid(false) {
}

ReferenceConstantPotentialSolver::~ReferenceConstantPotentialSolver() {
}

void ReferenceConstantPotentialSolver::invalidate() {
    valid = false;
}

ReferenceConstantPotentialMatrixSolver::ReferenceConstantPotentialMatrixSolver(int numElectrodeParticles) : ReferenceConstantPotentialSolver(),
    electrodePosData(numElectrodeParticles), constraintVector(numElectrodeParticles) {
}

void ReferenceConstantPotentialMatrixSolver::update(
    int numParticles,
    int numElectrodeParticles,
    const std::vector<Vec3>& posData,
    const std::vector<double>& charges,
    const std::vector<std::set<int>>& exclusions,
    const std::vector<int>& sysToElec,
    const std::vector<int>& elecToSys,
    const std::vector<std::array<double, 3> >& electrodeParamArray,
    ReferenceConstantPotential& conp
) {
    // Initializes or updates the precomputed capacitance matrix if this is its
    // first use or electrode parameters have changed since its initialization.

    // Check for changes to box vectors or electrode positions that might
    // invalidate a matrix that is currently marked valid.
    if (valid) {
        if (boxVectors[0] != conp.boxVectors[0] || boxVectors[1] != conp.boxVectors[1] || boxVectors[2] != conp.boxVectors[2]) {
            valid = false;
        }
    }
    if (valid) {
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            if (electrodePosData[ii] != posData[elecToSys[ii]]) {
                valid = false;
                break;
            }
        }
    }
    if (valid) {
        return;
    }

    // Store the current box vectors and electrode positions before updating the
    // capacitance matrix.
    valid = true;
    boxVectors[0] = conp.boxVectors[0];
    boxVectors[1] = conp.boxVectors[1];
    boxVectors[2] = conp.boxVectors[2];
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        electrodePosData[ii] = posData[elecToSys[ii]];
    }

    TNT::Array2D<double> A(numElectrodeParticles, numElectrodeParticles);
    std::vector<double> dUdQ0(numElectrodeParticles);
    std::vector<double> dUdQ(numElectrodeParticles);
    
    pme_t pmeData;
    pme_init(&pmeData, conp.ewaldAlpha, numParticles, conp.gridSize, 5, 1);

    // Get derivatives when all electrode charges are zeroed.
    std::vector<double> q(charges);
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        q[elecToSys[ii]] = 0.0;
    }
    conp.getDerivatives(numParticles, numElectrodeParticles, posData, q, exclusions, sysToElec, elecToSys, electrodeParamArray, dUdQ0, pmeData);

    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        int i = elecToSys[ii];

        // Get derivatives when one electrode charge is set.
        q[i] = 1.0;
        conp.getDerivatives(numParticles, numElectrodeParticles, posData, q, exclusions, sysToElec, elecToSys, electrodeParamArray, dUdQ, pmeData);
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
    capacitance = JAMA::Cholesky<double>(A);
    if (!capacitance.is_spd()) {
        throw OpenMMException("Electrode matrix not positive definite");
    }

    // Precompute the appropriate scaling vector to enforce constant total
    // charge if requested.  The vector is parallel to one obtained by solving
    // Aq = b for all q_i = 1 (ensuring the constrained charges will actually be
    // the correct constrained minimum of the quadratic form for the energy),
    // and is scaled so that adding it to a vector of charges increases the
    // total charge by 1 (making it easy to calculate the necessary offset).
    if (conp.useChargeConstraint) {
        TNT::Array1D<double> solution = capacitance.solve(TNT::Array1D<double>(numElectrodeParticles, 1.0));
        constraintVector.assign(static_cast<double*>(solution), static_cast<double*>(solution) + numElectrodeParticles);

        double constraintScaleInv = 0.0;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            constraintScaleInv += constraintVector[ii];
        }
        double constraintScale = 1.0 / constraintScaleInv;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            constraintVector[ii] *= constraintScale;
        }
    }
}

void ReferenceConstantPotentialMatrixSolver::solve(
    int numParticles,
    int numElectrodeParticles,
    const std::vector<Vec3>& posData,
    std::vector<double>& charges,
    const std::vector<std::set<int>>& exclusions,
    const std::vector<int>& sysToElec,
    const std::vector<int>& elecToSys,
    const std::vector<std::array<double, 3> >& electrodeParamArray,
    ReferenceConstantPotential& conp,
    pme_t pmeData
) {
    // Solves for charges using the matrix method.

    // Zero electrode charges and get derivatives at zero charge.
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        charges[elecToSys[ii]] = 0.0;
    }
    std::vector<double> b(numElectrodeParticles);
    conp.getDerivatives(numParticles, numElectrodeParticles, posData, charges, exclusions, sysToElec, elecToSys, electrodeParamArray, b, pmeData);
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        b[ii] = -b[ii];
    }

    // Solve for electrode charges directly using capacitance matrix and
    // calculated derivatives.
    TNT::Array1D<double> q = capacitance.solve(TNT::Array1D<double>(numElectrodeParticles, b.data()));
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        charges[elecToSys[ii]] = q[ii];
    }

    // Enforce total charge constraint if requested.
    if (conp.useChargeConstraint) {
        double chargeOffset = conp.chargeTarget;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            chargeOffset -= q[ii];
        }
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            charges[elecToSys[ii]] += chargeOffset * constraintVector[ii];
        }
    }
}

ReferenceConstantPotentialCGSolver::ReferenceConstantPotentialCGSolver(int numElectrodeParticles, bool precond) : ReferenceConstantPotentialSolver(),
    precond(precond), precondVector(numElectrodeParticles) {
}

void ReferenceConstantPotentialCGSolver::update(
    int numParticles,
    int numElectrodeParticles,
    const std::vector<Vec3>& posData,
    const std::vector<double>& charges,
    const std::vector<std::set<int>>& exclusions,
    const std::vector<int>& sysToElec,
    const std::vector<int>& elecToSys,
    const std::vector<std::array<double, 3> >& electrodeParamArray,
    ReferenceConstantPotential& conp
) {
    // Initializes or updates information for a preconditioner for the conjugate
    // gradient method if this is its first use or electrode parameters have
    // changed since its initialization.
    
    const double SQRT_PI = sqrt(PI_M);
    const double SELF_ALPHA_SCALE = 2.0 * ONE_4PI_EPS0 / SQRT_PI;
    const double SELF_ETA_SCALE = SELF_ALPHA_SCALE / sqrt(2.0);
    const double TF_SCALE = 1.0 / EPSILON0;
    const double PLASMA_SCALE = TF_SCALE / 4.0;

    // No action is required if the box vectors have not changed.
    if (valid && boxVectors[0] == conp.boxVectors[0] && boxVectors[1] == conp.boxVectors[1] && boxVectors[2] == conp.boxVectors[2]) {
        return;
    }

    valid = true;
    boxVectors[0] = conp.boxVectors[0];
    boxVectors[1] = conp.boxVectors[1];
    boxVectors[2] = conp.boxVectors[2];

    if (precond) {
        // Perform a reference PME calculation with a single charge at the origin to
        // find the constant offset on the preconditioner diagonal due to the PME
        // calculation.  This will actually vary slightly with position but only due
        // to finite accuracy of the PME splines, so it is fine to assume it will be
        // constant for the preconditioner.
        pme_t pmeData;
        pme_init(&pmeData, conp.ewaldAlpha, 1, conp.gridSize, 5, 1);
        std::vector<Vec3> pmePosData(1);
        std::vector<double> pmeChargeDerivatives(1);
        std::vector<int> pmeElectrodeIndices(1);
        std::vector<double> pmeCharges(1, 1.0);
        pme_exec_charge_derivatives(pmeData, pmePosData, pmeChargeDerivatives, pmeElectrodeIndices, pmeCharges, boxVectors);
        pme_destroy(pmeData);

        // The diagonal has a contribution from reciprocal space, Ewald
        // self-interaction, Ewald neutralizing plasma, Gaussian self-interaction,
        // and Thomas-Fermi contributions.
        double volume = boxVectors[0][0] * boxVectors[1][1] * boxVectors[2][2];
        double ewaldTerm = pmeChargeDerivatives[0] - SELF_ALPHA_SCALE * conp.ewaldAlpha - PLASMA_SCALE / (volume * conp.ewaldAlpha * conp.ewaldAlpha);

        double precondScaleInv = 0.0;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            precondVector[ii] = 1.0 / (SELF_ETA_SCALE / electrodeParamArray[ii][conp.GaussianWidthIndex]
                + TF_SCALE * electrodeParamArray[ii][conp.ThomasFermiScaleIndex] + ewaldTerm);
            precondScaleInv += precondVector[ii];
        }
        double precondScale = 1.0 / precondScaleInv;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            precondVector[ii] *= precondScale;
        }
    }
}

void ReferenceConstantPotentialCGSolver::solve(
    int numParticles,
    int numElectrodeParticles,
    const std::vector<Vec3>& posData,
    std::vector<double>& charges,
    const std::vector<std::set<int>>& exclusions,
    const std::vector<int>& sysToElec,
    const std::vector<int>& elecToSys,
    const std::vector<std::array<double, 3> >& electrodeParamArray,
    ReferenceConstantPotential& conp,
    pme_t pmeData
) {
    // Solves for charges using the conjugate gradient method.
    
    std::vector<double> q(numElectrodeParticles);
    std::vector<double> grad(numElectrodeParticles);
    std::vector<double> projGrad(numElectrodeParticles);
    std::vector<double> precGrad(numElectrodeParticles);
    std::vector<double> qStep(numElectrodeParticles);
    std::vector<double> gradStep(numElectrodeParticles);
    std::vector<double> grad0(numElectrodeParticles);
    double offset, error, paramScale, alpha, beta;

    const double errorTarget = conp.cgErrorTol * conp.cgErrorTol * numElectrodeParticles;

    // Ensure that initial guess charges satisfy the constraint.
    if (conp.useChargeConstraint) {
        offset = conp.chargeTarget;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            offset -= charges[elecToSys[ii]];
        }
        offset /= numElectrodeParticles;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            charges[elecToSys[ii]] += offset;
        }
    }

    // Evaluate the initial gradient Aq - b.
    conp.getDerivatives(numParticles, numElectrodeParticles, posData, charges, exclusions, sysToElec, elecToSys, electrodeParamArray, grad, pmeData);

    // Project the initial gradient without preconditioning.
    offset = 0.0;
    if (conp.useChargeConstraint) {
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            offset += grad[ii];
        }
        offset /= numElectrodeParticles;
    }
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        projGrad[ii] = grad[ii] - offset;
    }

    // Check for convergence at the initial guess charges.
    error = 0.0;
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        error += projGrad[ii] * projGrad[ii];
    }
    if (error <= errorTarget) {
        return;
    }

    // Save the current charges, then evaluate the gradient with zero
    // charges (-b) so that we can later compute Ap as (Ap - b) - (-b).
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        int i = elecToSys[ii];
        q[ii] = charges[i];
        charges[i] = 0.0;
    }
    conp.getDerivatives(numParticles, numElectrodeParticles, posData, charges, exclusions, sysToElec, elecToSys, electrodeParamArray, grad0, pmeData);

    // Project the initial gradient with preconditioning.
    if (precond) {
        offset = 0.0;
        if (conp.useChargeConstraint) {
            for (int ii = 0; ii < numElectrodeParticles; ii++) {
                offset += precondVector[ii] * grad[ii];
            }
        }
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            precGrad[ii] = precondVector[ii] * (grad[ii] - offset);
        }
    }
    else {
        precGrad.assign(projGrad.begin(), projGrad.end());
    }

    // Initialize step vector for conjugate gradient iterations.
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        qStep[ii] = -precGrad[ii];
    }

    // Perform conjugate gradient iterations.
    bool converged = false;
    for (int iter = 0; iter <= numElectrodeParticles; iter++) {
        // Evaluate the matrix-vector product A qStep.
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            charges[elecToSys[ii]] = qStep[ii];
        }
        conp.getDerivatives(numParticles, numElectrodeParticles, posData, charges, exclusions, sysToElec, elecToSys, electrodeParamArray, gradStep, pmeData);
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            gradStep[ii] -= grad0[ii];
        }

        // If A qStep is small enough, stop to prevent, e.g., division by
        // zero in the calculation of alpha, or too large step sizes.
        error = 0.0;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            error += gradStep[ii] * gradStep[ii];
        }
        if (error <= errorTarget) {
            converged = true;
            break;
        }
        
        // Evaluate the scalar 1 / (qStep^T A qStep).
        paramScale = 0.0;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            paramScale += qStep[ii] * gradStep[ii];
        }
        paramScale = 1.0 / paramScale;

        // Evaluate the conjugate gradient parameter alpha.
        alpha = 0.0;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            alpha -= qStep[ii] * grad[ii];
        }
        alpha *= paramScale;

        // Update the charge vector.
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            q[ii] += alpha * qStep[ii];
        }

        if (conp.useChargeConstraint) {
            // Remove any accumulated drift from the charge vector.  This
            // would be zero in exact arithmetic, but error can accumulate
            // over time in finite precision.
            offset = conp.chargeTarget;
            for (int ii = 0; ii < numElectrodeParticles; ii++) {
                offset -= q[ii];
            }
            offset /= numElectrodeParticles;
            for (int ii = 0; ii < numElectrodeParticles; ii++) {
                q[ii] += offset;
            }
        }

        // Update the gradient vector (but periodically recompute it instead
        // of updating to reduce the accumulation of roundoff error).
        if (iter != 0 && iter % 32 == 0) {
            for (int ii = 0; ii < numElectrodeParticles; ii++) {
                charges[elecToSys[ii]] = q[ii];
            }
            conp.getDerivatives(numParticles, numElectrodeParticles, posData, charges, exclusions, sysToElec, elecToSys, electrodeParamArray, grad, pmeData);
        }
        else {
            for (int ii = 0; ii < numElectrodeParticles; ii++) {
                grad[ii] += alpha * gradStep[ii];
            }
        }

        // Project the current gradient without preconditioning.
        offset = 0.0;
        if (conp.useChargeConstraint) {
            for (int ii = 0; ii < numElectrodeParticles; ii++) {
                offset += grad[ii];
            }
            offset /= numElectrodeParticles;
        }
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            projGrad[ii] = grad[ii] - offset;
        }

        // Check for convergence.
        error = 0.0;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            error += projGrad[ii] * projGrad[ii];
        }
        if (error <= errorTarget) {
            converged = true;
            break;
        }

        // Project the current gradient with preconditioning.
        if (precond) {
            offset = 0.0;
            if (conp.useChargeConstraint) {
                for (int ii = 0; ii < numElectrodeParticles; ii++) {
                    offset += precondVector[ii] * grad[ii];
                }
            }
            for (int ii = 0; ii < numElectrodeParticles; ii++) {
                precGrad[ii] = precondVector[ii] * (grad[ii] - offset);
            }
        }
        else {
            precGrad.assign(projGrad.begin(), projGrad.end());
        }

        // Evaluate the conjugate gradient parameter beta.
        beta = 0.0;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            beta += precGrad[ii] * gradStep[ii];
        }
        beta *= paramScale;

        // Update the step vector.
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            qStep[ii] = beta * qStep[ii] - precGrad[ii];
        }

        if (conp.useChargeConstraint) {
            // Project out any deviation off of the constraint plane from
            // the step vector.  This would be zero in exact arithmetic, but
            // error can accumulate over time in finite precision.
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
        charges[elecToSys[ii]] = q[ii];
    }
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


void ReferenceConstantPotential::execute(
    int numParticles,
    int numElectrodeParticles,
    const std::vector<Vec3>& posData,
    std::vector<Vec3>& forceData,
    std::vector<double>& charges,
    const std::vector<std::set<int>>& exclusions,
    const std::vector<int>& sysToElec,
    const std::vector<int>& elecToSys,
    const std::vector<std::array<double, 3> >& electrodeParamArray,
    double* energy,
    ReferenceConstantPotentialSolver* solver
) {
    pme_t pmeData;
    pme_init(&pmeData, ewaldAlpha, numParticles, gridSize, 5, 1);

    solver->solve(numParticles, numElectrodeParticles, posData, charges, exclusions, sysToElec, elecToSys, electrodeParamArray, *this, pmeData);
    getEnergyForces(numParticles, numElectrodeParticles, posData, forceData, charges, exclusions, sysToElec, elecToSys, electrodeParamArray, energy, pmeData);

    pme_destroy(pmeData);
}

void ReferenceConstantPotential::getCharges(
    int numParticles,
    int numElectrodeParticles,
    const std::vector<Vec3>& posData,
    std::vector<double>& charges,
    const std::vector<std::set<int>>& exclusions,
    const std::vector<int>& sysToElec,
    const std::vector<int>& elecToSys,
    const std::vector<std::array<double, 3> >& electrodeParamArray,
    ReferenceConstantPotentialSolver* solver
) {
    pme_t pmeData;
    pme_init(&pmeData, ewaldAlpha, numParticles, gridSize, 5, 1);

    solver->solve(numParticles, numElectrodeParticles, posData, charges, exclusions, sysToElec, elecToSys, electrodeParamArray, *this, pmeData);

    pme_destroy(pmeData);
}

void ReferenceConstantPotential::getEnergyForces(
    int numParticles,
    int numElectrodeParticles,
    const std::vector<Vec3>& posData,
    std::vector<Vec3>& forceData,
    std::vector<double>& charges,
    const std::vector<std::set<int>>& exclusions,
    const std::vector<int>& sysToElec,
    const std::vector<int>& elecToSys,
    const std::vector<std::array<double, 3> >& electrodeParamArray,
    double* energy,
    pme_t pmeData
) {
    const double SQRT_PI = sqrt(PI_M);
    const double TWO_OVER_SQRT_PI = 2.0 / SQRT_PI;
    const double SELF_ALPHA_SCALE = ONE_4PI_EPS0 / SQRT_PI;
    const double SELF_ETA_SCALE = SELF_ALPHA_SCALE / sqrt(2.0);
    const double TF_SCALE = 1.0 / (2.0 * EPSILON0);
    const double PLASMA_SCALE = TF_SCALE / 4.0;

    double energyAccum = 0.0;

    // Direct space.
    for (auto& pair : *neighborList) {
        int i = pair.first;
        int j = pair.second;
        int ii = sysToElec[i];
        int jj = sysToElec[j];

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
            if (exceptionsArePeriodic) {
                ReferenceForce::getDeltaRPeriodic(posData[i], posData[j], boxVectors, deltaR);
            }
            else {
                ReferenceForce::getDeltaR(posData[i], posData[j], deltaR);
            }
            double r = deltaR[ReferenceForce::RIndex];
            double inverseR = 1.0 / r;

            double qqFactor = ONE_4PI_EPS0 * charges[i] * charges[j];
            double alphaR = ewaldAlpha * r;

            if (alphaR > 1e-6) {
                double erfAlphaR = erf(alphaR);

                double forceFactor = qqFactor * inverseR * inverseR * inverseR * (erfAlphaR - TWO_OVER_SQRT_PI * alphaR * exp(-alphaR * alphaR));
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
    double qTotal = 0.0;
    for (int i = 0; i < numParticles; i++) {
        double q = charges[i];
        qTotal += q;
        energyAccum -= q * (SELF_ALPHA_SCALE * q * ewaldAlpha + posData[i].dot(externalField));
        forceData[i] += q * externalField;
    }

    // Ewald neutralizing plasma.
    double volume = boxVectors[0][0] * boxVectors[1][1] * boxVectors[2][2];
    energyAccum -= PLASMA_SCALE * qTotal * qTotal / (volume * ewaldAlpha * ewaldAlpha);

    // Gaussian self-interaction, potential, and Thomas-Fermi contributions
    // (electrode particles only).
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        int i = elecToSys[ii];
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
    const std::vector<int>& sysToElec,
    const std::vector<int>& elecToSys,
    const std::vector<std::array<double, 3> >& electrodeParamArray,
    std::vector<double>& chargeDerivatives,
    pme_t pmeData
) {
    const double SQRT_PI = sqrt(PI_M);
    const double SELF_ALPHA_SCALE = 2.0 * ONE_4PI_EPS0 / SQRT_PI;
    const double SELF_ETA_SCALE = SELF_ALPHA_SCALE / sqrt(2.0);
    const double TF_SCALE = 1.0 / EPSILON0;
    const double PLASMA_SCALE = TF_SCALE / 4.0;

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
        int ii = sysToElec[i];
        int jj = sysToElec[j];

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
    pme_exec_charge_derivatives(pmeData, posData, chargeDerivatives, elecToSys, charges, boxVectors);

    // Ewald neutralizing plasma precalculation.
    double qTotal = 0.0;
    for (int i = 0; i < numParticles; i++) {
        qTotal += charges[i];
    }
    double volume = boxVectors[0][0] * boxVectors[1][1] * boxVectors[2][2];
    double plasmaTerm = PLASMA_SCALE * qTotal / (volume * ewaldAlpha * ewaldAlpha);

    // Ewald self-interaction, Ewald neutralizing plasma, Gaussian
    // self-interaction, potential, external field, and Thomas-Fermi
    // contributions.
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        int i = elecToSys[ii];
        double q = charges[i];
        chargeDerivatives[ii] += q * (SELF_ETA_SCALE / electrodeParamArray[ii][GaussianWidthIndex] - SELF_ALPHA_SCALE * ewaldAlpha + TF_SCALE * electrodeParamArray[ii][ThomasFermiScaleIndex])
            - plasmaTerm - electrodeParamArray[ii][PotentialIndex] - posData[i].dot(externalField);
    }
}
