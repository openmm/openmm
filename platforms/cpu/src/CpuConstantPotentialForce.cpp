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

#include "CpuConstantPotentialForce.h"
#include "SimTKOpenMMRealType.h"

#include <algorithm>

// In case we're using some primitive version of Visual Studio this will
// make sure that erf() and erfc() are defined.
#include "openmm/internal/MSVC_erfc.h"

using namespace std;
using namespace OpenMM;

CpuConstantPotentialSolver::CpuConstantPotentialSolver(int numParticles, int numElectrodeParticles) :
    numParticles(numParticles),
    numElectrodeParticles(numElectrodeParticles),
    valid(false),
    hasSavedSolution(false),
    savedPositions(numParticles),
    savedCharges(numElectrodeParticles)
{
}

CpuConstantPotentialSolver::~CpuConstantPotentialSolver() {
}

void CpuConstantPotentialSolver::invalidate() {
    valid = false;
    hasSavedSolution = false;
}

void CpuConstantPotentialSolver::discardSavedSolution() {
    hasSavedSolution = false;
}

void CpuConstantPotentialSolver::solve(CpuConstantPotentialForce& conp, ThreadPool& threads, Kernel& pmeKernel) {
    // There is nothing to do if all particles have fixed charges.
    if(!numElectrodeParticles) {
        return;
    }

    // If box vectors or positions have not changed, and there is a solution
    // already computed, we can simply reload it instead of solving again.
    if (hasSavedSolution) {
        if (savedBoxVectors[0] != conp.boxVectors[0] || savedBoxVectors[1] != conp.boxVectors[1] || savedBoxVectors[2] != conp.boxVectors[2]) {
            hasSavedSolution = false;
        }
    }
    if (hasSavedSolution) {
        for (int i = 0; i < numParticles; i++) {
            if (savedPositions[i] != conp.posData[i]) {
                hasSavedSolution = false;
                break;
            }
        }
    }
    if (hasSavedSolution) {
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            conp.posq[4 * conp.elecToSys[ii] + 3] = savedCharges[ii];
        }
        return;
    }

    solveImpl(conp, threads, pmeKernel);

    hasSavedSolution = true;
    savedBoxVectors[0] = conp.boxVectors[0];
    savedBoxVectors[1] = conp.boxVectors[1];
    savedBoxVectors[2] = conp.boxVectors[2];
    savedPositions.assign(static_cast<const Vec3*>(conp.posData), static_cast<const Vec3*>(conp.posData + numParticles));
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        savedCharges[ii] = conp.posq[4 * conp.elecToSys[ii] + 3];
    }
}

CpuConstantPotentialMatrixSolver::CpuConstantPotentialMatrixSolver(int numParticles, int numElectrodeParticles) : CpuConstantPotentialSolver(numParticles, numElectrodeParticles),
        electrodePosData(numElectrodeParticles), constraintVector(numElectrodeParticles), b(numElectrodeParticles) {
}

void CpuConstantPotentialMatrixSolver::solveImpl(CpuConstantPotentialForce& conp, ThreadPool& threads, Kernel& pmeKernel) {
    ensureValid(conp, threads, pmeKernel);

    // Zero electrode charges and get derivatives at zero charge.
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        conp.posq[4 * conp.elecToSys[ii] + 3] = 0.0f;
    }
    conp.getDerivatives(threads, pmeKernel);
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        b[ii] = -conp.chargeDerivatives[ii];
    }

    // Solve for electrode charges directly using capacitance matrix and
    // calculated derivatives.
    TNT::Array1D<float> q = capacitance.solve(b);

    // Enforce total charge constraint if requested.
    if (conp.useChargeConstraint) {
        float chargeOffset = conp.chargeTarget;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            chargeOffset -= q[ii];
        }
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            q[ii] += chargeOffset * constraintVector[ii];
        }
    }

    // Copy solution out of solver memory into posq.
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        conp.posq[4 * conp.elecToSys[ii] + 3] = q[ii];
    }
}

void CpuConstantPotentialMatrixSolver::ensureValid(CpuConstantPotentialForce& conp, ThreadPool& threads, Kernel& pmeKernel) {
    // Initializes or updates the precomputed capacitance matrix if this is its
    // first use or electrode parameters have changed since its initialization.

    // In this method, we modify the current electrode charges in posq.  This
    // is fine since it will only get called at the start of solve(), which
    // overwrites the electrode charges in posq after solving for them.

    // Check for changes to box vectors or electrode positions that might
    // invalidate a matrix that is currently marked valid.
    if (valid) {
        if (boxVectors[0] != conp.boxVectors[0] || boxVectors[1] != conp.boxVectors[1] || boxVectors[2] != conp.boxVectors[2]) {
            valid = false;
        }
    }
    if (valid) {
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            if (electrodePosData[ii] != conp.posData[conp.elecToSys[ii]]) {
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
        electrodePosData[ii] = conp.posData[conp.elecToSys[ii]];
    }

    TNT::Array2D<float> A(numElectrodeParticles, numElectrodeParticles);
    vector<float> dUdQ0(numElectrodeParticles);
    vector<float> dUdQ(numElectrodeParticles);

    // Get derivatives when all electrode charges are zeroed.
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        conp.posq[4 * conp.elecToSys[ii] + 3] = 0.0f;
    }
    conp.getDerivatives(threads, pmeKernel);
    dUdQ0 = conp.chargeDerivatives;

    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        int i = conp.elecToSys[ii];

        // Get derivatives when one electrode charge is set.
        conp.posq[4 * i + 3] = 1.0f;
        conp.getDerivatives(threads, pmeKernel);
        dUdQ = conp.chargeDerivatives;
        conp.posq[4 * i + 3] = 0.0f;

        // Set matrix elements, subtracting zero charge derivatives so that the
        // matrix will end up being the (charge-independent) Hessian.
        for (int jj = 0; jj < ii; jj++) {
            A[ii][jj] = A[jj][ii] = dUdQ[jj] - dUdQ0[jj];
        }
        A[ii][ii] = dUdQ[ii] - dUdQ0[ii];
    }

    // Compute Cholesky decomposition representation of the inverse.
    capacitance = JAMA::Cholesky<float>(A);
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
        TNT::Array1D<float> solution = capacitance.solve(TNT::Array1D<float>(numElectrodeParticles, 1.0));
        constraintVector.assign(static_cast<float*>(solution), static_cast<float*>(solution) + numElectrodeParticles);

        float constraintScaleInv = 0.0;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            constraintScaleInv += constraintVector[ii];
        }
        float constraintScale = 1.0 / constraintScaleInv;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            constraintVector[ii] *= constraintScale;
        }
    }
}

CpuConstantPotentialCGSolver::CpuConstantPotentialCGSolver(int numParticles,int numElectrodeParticles, bool precond) : CpuConstantPotentialSolver(numParticles, numElectrodeParticles),
    precondRequested(precond),
    precondVector(numElectrodeParticles),
    q(numElectrodeParticles),
    grad(numElectrodeParticles),
    projGrad(numElectrodeParticles),
    precGrad(numElectrodeParticles),
    qStep(numElectrodeParticles),
    gradStep(numElectrodeParticles),
    grad0(numElectrodeParticles),
    qLast(numElectrodeParticles)
{
}

void CpuConstantPotentialCGSolver::solveImpl(CpuConstantPotentialForce& conp, ThreadPool& threads, Kernel& pmeKernel) {
    ensureValid(conp, threads, pmeKernel);
    
    double offset;
    float error, paramScale, alpha, beta;
    const float errorTarget = conp.cgErrorTol * conp.cgErrorTol * numElectrodeParticles;

    // Set initial guess charges as linear extrapolations from the current and
    // previous charges fed through the solver, and save the current charges as
    // the previous charges.
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        int i = conp.elecToSys[ii];
        float qGuess = conp.posq[4 * i + 3];
        conp.posq[4 * i + 3] = 2.0f * qGuess - qLast[ii];
        qLast[ii] = qGuess;
    }

    // Ensure that initial guess charges satisfy the constraint.
    if (conp.useChargeConstraint) {
        offset = conp.chargeTarget;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            offset -= conp.posq[4 * conp.elecToSys[ii] + 3];
        }
        offset /= numElectrodeParticles;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            conp.posq[4 * conp.elecToSys[ii] + 3] += offset;
        }
    }

    // Evaluate the initial gradient Aq - b.
    conp.getDerivatives(threads, pmeKernel);
    grad.assign(&conp.chargeDerivatives[0], &conp.chargeDerivatives[numElectrodeParticles]);

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
    error = 0.0f;
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        error += projGrad[ii] * projGrad[ii];
    }
    if (error <= errorTarget) {
        return;
    }

    // Save the current charges, then evaluate the gradient with zero
    // charges (-b) so that we can later compute Ap as (Ap - b) - (-b).
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        int i = conp.elecToSys[ii];
        q[ii] = conp.posq[4 * i + 3];
        conp.posq[4 * i + 3] = 0.0f;
    }
    conp.getDerivatives(threads, pmeKernel);
    grad0.assign(&conp.chargeDerivatives[0], &conp.chargeDerivatives[numElectrodeParticles]);

    // Project the initial gradient with preconditioning.
    if (precondActivated) {
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
            conp.posq[4 * conp.elecToSys[ii] + 3] = qStep[ii];
        }
        conp.getDerivatives(threads, pmeKernel);
        gradStep.assign(&conp.chargeDerivatives[0], &conp.chargeDerivatives[numElectrodeParticles]);
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            gradStep[ii] -= grad0[ii];
        }

        // If A qStep is small enough, stop to prevent, e.g., division by
        // zero in the calculation of alpha, or too large step sizes.
        error = 0.0f;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            error += gradStep[ii] * gradStep[ii];
        }
        if (error <= errorTarget) {
            converged = true;
            break;
        }

        // Evaluate the scalar 1 / (qStep^T A qStep).
        paramScale = 0.0f;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            paramScale += qStep[ii] * gradStep[ii];
        }
        paramScale = 1.0f / paramScale;

        // Evaluate the conjugate gradient parameter alpha.
        alpha = 0.0f;
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
                conp.posq[4 * conp.elecToSys[ii] + 3] = q[ii];
            }
            conp.getDerivatives(threads, pmeKernel);
            grad.assign(&conp.chargeDerivatives[0], &conp.chargeDerivatives[numElectrodeParticles]);
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
        error = 0.0f;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            error += projGrad[ii] * projGrad[ii];
        }
        if (error <= errorTarget) {
            converged = true;
            break;
        }

        // Project the current gradient with preconditioning.
        if (precondActivated) {
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
        beta = 0.0f;
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

    // Store the final charges in posq.
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        conp.posq[4 * conp.elecToSys[ii] + 3] = q[ii];
    }
}

void CpuConstantPotentialCGSolver::ensureValid(CpuConstantPotentialForce& conp, ThreadPool& threads, Kernel& pmeKernel) {
    // Initializes or updates information for a preconditioner for the conjugate
    // gradient method if this is its first use or electrode parameters have
    // changed since its initialization.

    // No action is required if the box vectors have not changed.
    if (valid && boxVectors[0] == conp.boxVectors[0] && boxVectors[1] == conp.boxVectors[1] && boxVectors[2] == conp.boxVectors[2]) {
        return;
    }

    valid = true;
    boxVectors[0] = conp.boxVectors[0];
    boxVectors[1] = conp.boxVectors[1];
    boxVectors[2] = conp.boxVectors[2];

    precondActivated = false;
    if (precondRequested && numElectrodeParticles) {
        // If electrode self-energy contributions differ between electrodes,
        // a preconditioner may help convergence; otherwise, it provides no
        // benefit and may slow convergence due to roundoff error.
        for (int ie = 1; ie < conp.numElectrodes; ie++) {
            // Note electrodeSelfScales[1] has the scale for electrode 0, etc.
            if (conp.electrodeSelfScales[ie + 1] != conp.electrodeSelfScales[ie]) {
                precondActivated = true;
                break;
            }
        }
    }

    if (precondActivated) {
        // Save the position of the first electrode particle.
        int i0 = conp.elecToSys[0];
        float x0 = conp.posq[4 * i0];
        float y0 = conp.posq[4 * i0 + 1];
        float z0 = conp.posq[4 * i0 + 2];

        // Save the charges of all particles, and then zero them.
        vector<float> qSave(numParticles);
        for (int i = 0; i < numParticles; i++) {
            qSave[i] = conp.posq[4 * i + 3];
            conp.posq[4 * i + 3] = 0.0f;
        }

        // Place a unit charge at the origin.
        conp.posq[4 * i0] = 0.0f;
        conp.posq[4 * i0 + 1] = 0.0f;
        conp.posq[4 * i0 + 2] = 0.0f;
        conp.posq[4 * i0 + 3] = 1.0f;

        // Perform a reference PME calculation with a single charge at the
        // origin to find the constant offset on the preconditioner diagonal due
        // to the PME calculation.  This will actually vary slightly with
        // position but only due to finite accuracy of the PME splines, so it is
        // fine to assume it will be constant for the preconditioner.
        vector<float> derivatives(numElectrodeParticles);
        CpuConstantPotentialPmeIO io(conp.posq, NULL, &derivatives[0], numParticles, numElectrodeParticles);
        pmeKernel.getAs<CalcPmeReciprocalForceKernel>().beginComputation(io, boxVectors, false, false, true);
        pmeKernel.getAs<CalcPmeReciprocalForceKernel>().finishComputation(io);
        float pmeTerm = derivatives[0];

        // Restore particle positions and charges.
        conp.posq[4 * i0] = x0;
        conp.posq[4 * i0 + 1] = y0;
        conp.posq[4 * i0 + 2] = z0;
        for (int i = 0; i < numParticles; i++) {
            conp.posq[4 * i + 3] = qSave[i];
        }

        // The diagonal has a contribution from reciprocal space, Ewald
        // self-interaction, Ewald neutralizing plasma, Gaussian self-interaction,
        // and Thomas-Fermi contributions.
        double precondScaleInv = 0.0;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            precondVector[ii] = 1.0 / (2.0 * (conp.electrodeSelfScales[conp.elecElec[ii] + 1] - conp.plasmaScale) + pmeTerm);
            precondScaleInv += precondVector[ii];
        }
        double precondScale = 1.0 / precondScaleInv;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            precondVector[ii] *= precondScale;
        }
    }
}

const int CpuConstantPotentialForce::PotentialIndex = 0;
const int CpuConstantPotentialForce::GaussianWidthIndex = 1;
const int CpuConstantPotentialForce::ThomasFermiScaleIndex = 2;

const int CpuConstantPotentialForce::NUM_TABLE_POINTS = 2048;

const double CpuConstantPotentialForce::TWO_OVER_SQRT_PI = 2.0 / sqrt(PI_M);
const double CpuConstantPotentialForce::SELF_ALPHA_SCALE = ONE_4PI_EPS0 / sqrt(PI_M);
const double CpuConstantPotentialForce::SELF_ETA_SCALE = ONE_4PI_EPS0 / sqrt(2.0 * PI_M);
const double CpuConstantPotentialForce::SELF_TF_SCALE = 1.0 / (2.0 * EPSILON0);

CpuConstantPotentialForce::CpuConstantPotentialForce() {
}

CpuConstantPotentialForce::~CpuConstantPotentialForce() {
}

void CpuConstantPotentialForce::initialize(
    int numParticles,
    int numElectrodeParticles,
    int posqIndex,
    float nonbondedCutoff,
    float ewaldAlpha,
    float cgErrorTol,
    const int* gridSize,
    bool exceptionsArePeriodic,
    bool useChargeConstraint,
    const CpuNeighborList& neighborList,
    CpuConstantPotentialSolver* solver,
    const vector<set<int> >& exclusions,
    const vector<int>& sysToElec,
    const vector<int>& elecToSys,
    const vector<int>& sysElec,
    const vector<int>& elecElec,
    const vector<array<double, 3> >& electrodeParams,
    float chargeTarget,
    const float* externalField
) {
    this->numParticles = numParticles;
    this->numElectrodeParticles = numElectrodeParticles;
    this->posqIndex = posqIndex;
    this->nonbondedCutoff = nonbondedCutoff;
    this->ewaldAlpha = ewaldAlpha;
    this->cgErrorTol = cgErrorTol;
    this->gridSize[0] = gridSize[0];
    this->gridSize[1] = gridSize[1];
    this->gridSize[2] = gridSize[2];
    this->exceptionsArePeriodic = exceptionsArePeriodic;
    this->useChargeConstraint = useChargeConstraint;
    this->neighborList = &neighborList;
    this->solver = solver;

    // We store pointers to exclusions, sysToElec, elecToSys, sysElec, and
    // elecElec, which should not change, and to electrodeParams, whose values
    // might be updated.  If this happens, update() will be called with the
    // electrode indices that have new parameters.
    this->exclusions = &exclusions[0];
    this->sysToElec = &sysToElec[0];
    this->elecToSys = &elecToSys[0];
    this->sysElec = &sysElec[0];
    this->elecElec = &elecElec[0];
    this->electrodeParams = &electrodeParams[0];

    numElectrodes = electrodeParams.size();
    electrodePotentials.resize(numElectrodes + 1);
    electrodeSelfScales.resize(numElectrodes + 1);
    for (int ie = -1; ie < numElectrodes; ie++) {
        // Electrode indices of particles range from -1 to numElectrodes - 1,
        // with -1 indicating a non-electrode particle.  Precompute electrode
        // parameters, storing values for [-1, numElectrodes - 1] as [0,
        // numElectrodes] for convenient lookup.
        double electrodePotential = 0.0;
        double electrodeSelfScale = -SELF_ALPHA_SCALE * ewaldAlpha;
        if (ie != -1) {
            electrodePotential = electrodeParams[ie][PotentialIndex];
            electrodeSelfScale += SELF_ETA_SCALE / electrodeParams[ie][GaussianWidthIndex] + SELF_TF_SCALE * electrodeParams[ie][ThomasFermiScaleIndex];
        }
        electrodePotentials[ie + 1] = (float) electrodePotential;
        electrodeSelfScales[ie + 1] = (float) electrodeSelfScale;
    }

    this->chargeTarget = chargeTarget;
    this->externalField = fvec4(externalField[0], externalField[1], externalField[2], 0.0f);

    chargeDerivatives.resize(numElectrodeParticles);
    energyLookupTable.resize((numElectrodes + 1) * (numElectrodes + 1) * (NUM_TABLE_POINTS + 4));
    forceLookupTable.resize((numElectrodes + 1) * (numElectrodes + 1) * (NUM_TABLE_POINTS + 4));
    tableScale = NUM_TABLE_POINTS / nonbondedCutoff;

    for (int ie = -1; ie < numElectrodes; ie++) {
        for (int je = -1; je <= ie; je++) {
            initializeLookupTables(ie, je);
        }
    }
}

void CpuConstantPotentialForce::update(float chargeTarget, const float* externalField, int firstElectrode, int lastElectrode) {
    solver->discardSavedSolution();
    if (firstElectrode <= lastElectrode) {
        solver->invalidate();
    }

    for (int ie = firstElectrode; ie <= lastElectrode; ie++) {
        electrodePotentials[ie + 1] = (float) electrodeParams[ie][PotentialIndex];
        electrodeSelfScales[ie + 1] = (float) (-SELF_ALPHA_SCALE * ewaldAlpha + SELF_ETA_SCALE / electrodeParams[ie][GaussianWidthIndex] + SELF_TF_SCALE * electrodeParams[ie][ThomasFermiScaleIndex]);
    }

    this->chargeTarget = chargeTarget;
    this->externalField = fvec4(externalField[0], externalField[1], externalField[2], 0.0f);

    for (int ie = firstElectrode; ie <= lastElectrode; ie++) {
        for (int je = -1; je < numElectrodes; je++) {
            if (je >= firstElectrode && je <= lastElectrode && je > ie) {
                continue;
            }
            initializeLookupTables(ie, je);
        }
    }
}

void CpuConstantPotentialForce::execute(
    const Vec3* boxVectors,
    const vector<Vec3>& posData,
    vector<float>& charges,
    float* posq,
    vector<AlignedArray<float> >& threadForce,
    double* energy,
    ThreadPool& threads,
    Kernel& pmeKernel
) {
    setThreadData(boxVectors, posData, posq, threadForce);
    if (numElectrodeParticles) {
        solver->solve(*this, threads, pmeKernel);
    }
    getEnergyForces(threads, pmeKernel, energy);
    saveCharges(charges);
}

void CpuConstantPotentialForce::getCharges(
    const Vec3* boxVectors,
    const vector<Vec3>& posData,
    vector<float>& charges,
    float* posq,
    vector<AlignedArray<float> >& threadForce,
    ThreadPool& threads,
    Kernel& pmeKernel
) {
    setThreadData(boxVectors, posData, posq, threadForce);
    if (numElectrodeParticles) {
        solver->solve(*this, threads, pmeKernel);
    }
    saveCharges(charges);
}

void CpuConstantPotentialForce::setThreadData(const Vec3* boxVectors, const vector<Vec3>& posData, float* posq, vector<AlignedArray<float> >& threadForce) {
    this->boxVectors[0] = boxVectors[0];
    this->boxVectors[1] = boxVectors[1];
    this->boxVectors[2] = boxVectors[2];

    boxVectorsVec4.resize(3);
    boxVectorsVec4[0] = fvec4((float) boxVectors[0][0], (float) boxVectors[0][1], (float) boxVectors[0][2], 0.0f);
    boxVectorsVec4[1] = fvec4((float) boxVectors[1][0], (float) boxVectors[1][1], (float) boxVectors[1][2], 0.0f);
    boxVectorsVec4[2] = fvec4((float) boxVectors[2][0], (float) boxVectors[2][1], (float) boxVectors[2][2], 0.0f);

    boxSize = fvec4(boxVectorsVec4[0][0], boxVectorsVec4[1][1], boxVectorsVec4[2][2], 0.0f);
    recipBoxSize = fvec4((float) (1.0 / boxVectors[0][0]), (float) (1.0 / boxVectors[1][1]), (float) (1.0 / boxVectors[2][2]), 0.0f);

    triclinic = boxVectors[0][1] != 0.0 || boxVectors[0][2] != 0.0 ||
        boxVectors[1][0] != 0.0 || boxVectors[1][2] != 0.0 ||
        boxVectors[2][0] != 0.0 || boxVectors[2][1] != 0.0;

    plasmaScale = (float) (1.0 / (8.0 * EPSILON0 * boxVectors[0][0] * boxVectors[1][1] * boxVectors[2][2] * ewaldAlpha * ewaldAlpha));

    this->posData = &posData[0];
    this->posq = posq;
    this->threadForce = &threadForce;
}

void CpuConstantPotentialForce::saveCharges(vector<float>& charges) {
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        int i = elecToSys[ii];
        charges[i] = posq[4 * i + 3];
    }
}

void CpuConstantPotentialForce::getEnergyForces(ThreadPool& threads, Kernel& pmeKernel, double* energy) {
    // Everything except reciprocal space.
    int numThreads = threads.getNumThreads();
    threadEnergy.resize(numThreads);
    atomicBlockCounter = 0;
    atomicParticleCounter = 0;
    threads.execute([&] (ThreadPool& threads, int threadIndex) { getEnergyForcesThread(threads, threadIndex, energy != NULL); });
    threads.waitForThreads();

    // Accumulate thread energies.
    if (energy != NULL) {
        double energyAccum = 0.0;
        for (int i = 0; i < numThreads; i++) {
            energyAccum += threadEnergy[i];
        }
        *energy += energyAccum;
    }

    // Reciprocal space.
    CpuConstantPotentialPmeIO io(posq, &(*threadForce)[0][0], NULL, numParticles, numElectrodeParticles);
    pmeKernel.getAs<CalcPmeReciprocalForceKernel>().beginComputation(io, boxVectors, energy != NULL, true, false);
    double pmeEnergy = pmeKernel.getAs<CalcPmeReciprocalForceKernel>().finishComputation(io);
    if (energy != NULL) {
        // Ewald neutralizing plasma.
        float qTotal = 0.0;
        for (int i = 0; i < numParticles; i++) {
            qTotal += posq[4 * i + 3];
        }
        pmeEnergy -= plasmaScale * qTotal * qTotal;

        // Accumulate energy.
        *energy += pmeEnergy;
    }
}

void CpuConstantPotentialForce::getEnergyForcesThread(ThreadPool& threads, int threadIndex, bool includeEnergy) {
    int numThreads = threads.getNumThreads();
    int groupSize = max(1, numParticles / (10 * numThreads));
    float* forces = &(*threadForce)[threadIndex][0];

    if (includeEnergy) {
        threadEnergy[threadIndex] = 0.0;
    }

    // Direct space.
    while (true) {
        int nextBlock = atomicBlockCounter++;
        if (nextBlock >= neighborList->getNumBlocks()) {
            break;
        }

        getEnergyForcesBlock(nextBlock, forces, includeEnergy ? &threadEnergy[threadIndex] : NULL);
    }

    while (true) {
        int start = atomicParticleCounter.fetch_add(groupSize);
        if (start >= numParticles) {
            break;
        }
        int end = min(start + groupSize, numParticles);

        for (int i = start; i < end; i++) {
            fvec4 posI((float) posData[i][0], (float) posData[i][1], (float) posData[i][2], 0.0f);
            float chargeI = posq[4 * i + 3];
            float chargeIScaled = ONE_4PI_EPS0 * chargeI;

            // Exceptions.
            for (int j : exclusions[i]) {
                if (j <= i) {
                    continue;
                }

                fvec4 posJ((float) posData[j][0], (float) posData[j][1], (float) posData[j][2], 0.0f);
                float chargeJ = posq[4 * j + 3];

                fvec4 deltaR;
                float r2;
                getExceptionDeltaR(posI, posJ, deltaR, r2);
                float r = sqrt(r2);
                float inverseR = 1.0f / r;

                float qqFactor = chargeIScaled * chargeJ;
                float alphaR = ewaldAlpha * r;

                if (alphaR > 1e-6) {
                    float erfAlphaR = erf(alphaR);

                    float qqFactorInverseR = qqFactor * inverseR;
                    float forceFactor = qqFactorInverseR * inverseR * inverseR * (erfAlphaR - (float) TWO_OVER_SQRT_PI * alphaR * exp(-alphaR * alphaR));
                    fvec4 force = forceFactor * deltaR;
                    (fvec4(forces + 4 * i) + force).store(forces + 4 * i);
                    (fvec4(forces + 4 * j) - force).store(forces + 4 * j);

                    if (includeEnergy) {
                        threadEnergy[threadIndex] -= qqFactorInverseR * erfAlphaR;
                    }
                }
                else if (includeEnergy) {
                    threadEnergy[threadIndex] -= qqFactor * (float) TWO_OVER_SQRT_PI * ewaldAlpha;
                }
            }

            if (includeEnergy) {
                // Ewald self-interaction and external field contributions (all
                // particles), along with Gaussian self-interaction, potential,
                // and Thomas-Fermi contributions (electrode particles only).
                int iePlus1 = sysElec[i] + 1;
                threadEnergy[threadIndex] += chargeI * (chargeI * electrodeSelfScales[iePlus1] - electrodePotentials[iePlus1] - dot3(posI, externalField));
            }

            // External field force.
            (fvec4(forces + 4 * i) + chargeI * externalField).store(forces + 4 * i);
        }
    }
}

void CpuConstantPotentialForce::getDerivatives(ThreadPool& threads, Kernel& pmeKernel) {
    // Ewald neutralizing plasma.
    float qTotal = 0.0f;
    for (int i = 0; i < numParticles; i++) {
        qTotal += posq[4 * i + 3];
    }
    float plasmaTerm = -2.0f * plasmaScale * qTotal;

    // Everything except reciprocal space (but including neutralizing plasma).
    int numThreads = threads.getNumThreads();
    atomicBlockCounter = 0;
    threads.execute([&] (ThreadPool& threads, int threadIndex) { getDerivativesThread(threads, threadIndex, plasmaTerm); });
    threads.waitForThreads();

    // Accumulate thread derivatives.
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        chargeDerivatives[ii] = 0.0f;
        for (int j = 0; j < numThreads; j++) {
            chargeDerivatives[ii] += (*threadForce)[j][4 * ii + 3];
        }
    }

    // Reciprocal space.
    CpuConstantPotentialPmeIO io(posq, NULL, &chargeDerivatives[0], numParticles, numElectrodeParticles);
    pmeKernel.getAs<CalcPmeReciprocalForceKernel>().beginComputation(io, boxVectors, false, false, true);
    pmeKernel.getAs<CalcPmeReciprocalForceKernel>().finishComputation(io);
}

void CpuConstantPotentialForce::getDerivativesThread(
    ThreadPool& threads,
    int threadIndex,
    float plasmaTerm
) {
    int numThreads = threads.getNumThreads();

    // Zero derivatives for all electrode particles for this thread.  Use the
    // fourth component of the threadForce array for the derivatives.
    float* derivatives = &(*threadForce)[threadIndex][3];
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        derivatives[4 * ii] = 0.0f;
    }

    // Direct space.
    while (true) {
        int nextBlock = atomicBlockCounter++;
        if (nextBlock >= neighborList->getNumBlocks()) {
            break;
        }

        getDerivativesBlock(nextBlock, derivatives);
    }

    int start = threadIndex * numElectrodeParticles / numThreads;
    int end = (threadIndex + 1) * numElectrodeParticles / numThreads;

    // Ewald self-interaction and external field contributions, along with
    // Gaussian self-interaction, potential, and Thomas-Fermi contributions.
    for (int ii = start; ii < end; ii++) {
        int i = elecToSys[ii];
        int iePlus1 = elecElec[ii] + 1;
        fvec4 posI((float) posData[i][0], (float) posData[i][1], (float) posData[i][2], 0.0f);
        derivatives[4 * ii] += 2.0f * posq[4 * i + 3] * electrodeSelfScales[iePlus1] - electrodePotentials[iePlus1] - dot3(posI, externalField) + plasmaTerm;
    }
}

void CpuConstantPotentialForce::getExceptionDeltaR(const fvec4& posI, const fvec4& posJ, fvec4& deltaR, float& r2) const {
    deltaR = posJ - posI;
    if (exceptionsArePeriodic) {
        if (triclinic) {
            deltaR -= boxVectorsVec4[2] * floorf(deltaR[2] * recipBoxSize[2] + 0.5f);
            deltaR -= boxVectorsVec4[1] * floorf(deltaR[1] * recipBoxSize[1] + 0.5f);
            deltaR -= boxVectorsVec4[0] * floorf(deltaR[0] * recipBoxSize[0] + 0.5f);
        }
        else {
            fvec4 base = round(deltaR * recipBoxSize) * boxSize;
            deltaR = deltaR - base;
        }
    }
    r2 = dot3(deltaR, deltaR);
}

void CpuConstantPotentialForce::initializeLookupTables(int ie, int je) {
    float * energyTableIJ = &energyLookupTable[(((ie + 1) * (numElectrodes + 1)) + je + 1) * (NUM_TABLE_POINTS + 4)];
    float * forceTableIJ = &forceLookupTable[(((ie + 1) * (numElectrodes + 1)) + je + 1) * (NUM_TABLE_POINTS + 4)];

    bool ieIsElectrode = ie != -1;
    bool jeIsElectrode = je != -1;
    bool involvesElectrode = ieIsElectrode || jeIsElectrode;

    double iWidth = ieIsElectrode ? electrodeParams[ie][GaussianWidthIndex] : 0.0;
    double jWidth = jeIsElectrode ? electrodeParams[je][GaussianWidthIndex] : 0.0;
    double eta = involvesElectrode ? 1.0 / sqrt(iWidth * iWidth + jWidth * jWidth) : 0.0;

    for (int iTable = 0; iTable < NUM_TABLE_POINTS + 4; iTable++) {
        double r = (nonbondedCutoff * iTable) / NUM_TABLE_POINTS;
        double alphaR = ewaldAlpha * r;
        double etaR = eta * r;

        energyTableIJ[iTable] = erfc(alphaR) - (involvesElectrode ? erfc(etaR) : 0.0);
        forceTableIJ[iTable] = energyTableIJ[iTable] + TWO_OVER_SQRT_PI * (alphaR * exp(-alphaR * alphaR) - (involvesElectrode ? etaR * exp(-etaR * etaR) : 0.0));
    }

    if (ie != je) {
        float * energyTableJI = &energyLookupTable[(((je + 1) * (numElectrodes + 1)) + ie + 1) * (NUM_TABLE_POINTS + 4)];
        float * forceTableJI = &forceLookupTable[(((je + 1) * (numElectrodes + 1)) + ie + 1) * (NUM_TABLE_POINTS + 4)];

        for (int iTable = 0; iTable < NUM_TABLE_POINTS + 4; iTable++) {
            energyTableJI[iTable] = energyTableIJ[iTable];
            forceTableJI[iTable] = forceTableIJ[iTable];
        }
    }
}

CpuConstantPotentialPmeIO::CpuConstantPotentialPmeIO(float* posq, float* force, float* chargeDerivatives, int numParticles, int numElectrodeParticles) :
    posq(posq), force(force), chargeDerivatives(chargeDerivatives), numParticles(numParticles), numElectrodeParticles(numElectrodeParticles) {
}

float* CpuConstantPotentialPmeIO::getPosq() {
    return posq;
}

void CpuConstantPotentialPmeIO::setForce(float* f) {
    for (int i = 0; i < numParticles; i++) {
        force[4 * i] += f[4 * i];
        force[4 * i + 1] += f[4 * i + 1];
        force[4 * i + 2] += f[4 * i + 2];
    }
}

void CpuConstantPotentialPmeIO::setChargeDerivatives(float* derivatives) {
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        chargeDerivatives[ii] += derivatives[ii];
    }
}
