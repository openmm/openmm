/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2008-2025 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Evan Pretti                                        *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "openmm/Context.h"
#include "openmm/internal/ConstantPotentialForceImpl.h"
#include "openmm/common/BondedUtilities.h"
#include "openmm/common/CommonCalcConstantPotentialForce.h"
#include "openmm/common/CommonKernelUtilities.h"
#include "openmm/common/ContextSelector.h"
#include "openmm/common/NonbondedUtilities.h"
#include "CommonKernelSources.h"
#include "SimTKOpenMMRealType.h"
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <iterator>
#include <set>
#include "tnt_array2d.h"
#include "jama_cholesky.h"

using namespace OpenMM;
using namespace std;

class CommonCalcConstantPotentialForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const ConstantPotentialForce& force, const vector<int>& sysElec, const vector<mm_double4>& electrodeParams) : force(force), sysElec(sysElec), electrodeParams(electrodeParams) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        int elec1 = sysElec[particle1];
        int elec2 = sysElec[particle2];
        if (elec1 == -1 && elec2 == -1) {
            // Both particles are non-electrode; check their fixed charges.
            double charge1, charge2;
            force.getParticleParameters(particle1, charge1);
            force.getParticleParameters(particle2, charge2);
            return (charge1 == charge2);
        }
        if (elec1 == -1 || elec2 == -1) {
            // One particle is non-electrode (but not both, handled above).
            return false;
        }
        // Both particles are electrode particles.
        mm_double4 electrodeParams1 = electrodeParams[elec1];
        mm_double4 electrodeParams2 = electrodeParams[elec2];
        return electrodeParams1.x == electrodeParams2.x && electrodeParams1.y == electrodeParams2.y && electrodeParams1.z == electrodeParams2.z;
    }
    int getNumParticleGroups() {
        return force.getNumExceptions();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int particle1, particle2;
        double chargeProd;
        force.getExceptionParameters(index, particle1, particle2, chargeProd);
        particles.resize(2);
        particles[0] = particle1;
        particles[1] = particle2;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2;
        double chargeProd1, chargeProd2;
        force.getExceptionParameters(group1, particle1, particle2, chargeProd1);
        force.getExceptionParameters(group2, particle1, particle2, chargeProd2);
        return (chargeProd1 == chargeProd2);
    }
private:
    const ConstantPotentialForce& force;
    const vector<int>& sysElec;
    const vector<mm_double4>& electrodeParams;
};

class CommonCalcConstantPotentialForceKernel::ReorderListener : public ComputeContext::ReorderListener {
public:
    ReorderListener(ComputeContext& cc, int numElectrodeParticles, int paddedProblemSize, const vector<int>& sysToElec, const vector<int>& elecToSys) :
            cc(cc), numElectrodeParticles(numElectrodeParticles), paddedProblemSize(paddedProblemSize), sysToElec(sysToElec), elecToSys(elecToSys) {
        numParticles = cc.getNumAtoms();
        lastOrder.assign(cc.getAtomIndex().begin(), cc.getAtomIndex().end());
    }
    void addChargeArray(ComputeArray& chargeArray) {
        chargeArrays.push_back(&chargeArray);
    }
    void execute() {
        // Reorder guess charges.
        if (chargeArrays.empty()) {
            return;
        }
        const vector<int>& order = cc.getAtomIndex();
        for (int index = 0; index < chargeArrays.size(); index++) {
            ComputeArray* chargeArray = chargeArrays[index];

            vector<double> charges(paddedProblemSize), swap(numElectrodeParticles);
            chargeArray->download(charges, true);
            for (int ii = 0; ii < numElectrodeParticles; ii++) {
                swap[sysToElec[lastOrder[elecToSys[ii]]]] = charges[ii];
            }
            for (int ii = 0; ii < numElectrodeParticles; ii++) {
                charges[ii] = swap[sysToElec[order[elecToSys[ii]]]];
            }
            chargeArray->upload(charges, true);
        }
        lastOrder.assign(order.begin(), order.end());
    }
private:
    ComputeContext& cc;
    int numParticles, numElectrodeParticles, paddedProblemSize;
    const vector<int>& sysToElec;
    const vector<int>& elecToSys;
    vector<int> lastOrder;
    vector<ComputeArray*> chargeArrays;
};

class CommonCalcConstantPotentialForceKernel::InvalidatePostComputation : public ComputeContext::ForcePostComputation {
public:
    InvalidatePostComputation(ComputeContext& cc, CommonConstantPotentialSolver* solver) : cc(cc), solver(solver) {
    }
    double computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        if(!cc.getForcesValid()) {
            solver->discardSavedSolution();
        }
        return 0.0;
    }
private:
    ComputeContext& cc;
    CommonConstantPotentialSolver* solver;
};

CommonConstantPotentialSolver::CommonConstantPotentialSolver(ComputeContext& cc, int numParticles, int numElectrodeParticles, int paddedProblemSize) :
    numParticles(numParticles),
    numElectrodeParticles(numElectrodeParticles),
    paddedProblemSize(paddedProblemSize),
    valid(false),
    hasSavedSolution(false)
{
    savedPositions.initialize(cc, cc.getPosq().getSize(), cc.getPosq().getElementSize(), "savedPositions");
    savedCharges.initialize(cc, paddedProblemSize, cc.getUseDoublePrecision() ? sizeof(double) : sizeof(float), "savedCharges");
    checkSavedPositionsKernelResult.initialize<int>(cc, 1, "checkSavedPositionsKernelResult");
}

CommonConstantPotentialSolver::~CommonConstantPotentialSolver() {
}

void CommonConstantPotentialSolver::compileKernels(CommonCalcConstantPotentialForceKernel& kernel) {
    map<string, string> defines;
    defines["NUM_PARTICLES"] = kernel.cc.intToString(numParticles);
    ComputeProgram program = kernel.cc.compileProgram(CommonKernelSources::constantPotentialSolver, defines);

    checkSavedPositionsKernel = program->createKernel("checkSavedPositions");
    checkSavedPositionsKernel->addArg(kernel.cc.getPosq());
    checkSavedPositionsKernel->addArg(savedPositions);
    checkSavedPositionsKernel->addArg(checkSavedPositionsKernelResult);
}

void CommonConstantPotentialSolver::invalidate() {
    valid = false;
    hasSavedSolution = false;
}

void CommonConstantPotentialSolver::discardSavedSolution() {
    hasSavedSolution = false;
}

void CommonConstantPotentialSolver::solve(CommonCalcConstantPotentialForceKernel& kernel) {
    // If box vectors or positions have not changed, and there is a solution
    // already computed, we can simply reload it instead of solving again.
    if (hasSavedSolution) {
        if (savedBoxVectors[0] != kernel.boxVectors[0] || savedBoxVectors[1] != kernel.boxVectors[1] || savedBoxVectors[2] != kernel.boxVectors[2]) {
            hasSavedSolution = false;
        }
    }
    if (hasSavedSolution) {
        kernel.cc.clearBuffer(checkSavedPositionsKernelResult);
        checkSavedPositionsKernel->execute(numParticles);
        int result;
        checkSavedPositionsKernelResult.download(&result);
        if (result) {
            hasSavedSolution = false;
        }
    }
    if (hasSavedSolution) {
        savedCharges.copyTo(kernel.electrodeCharges);
        kernel.mustUpdateElectrodeCharges = true;
        return;
    }

    solveImpl(kernel);

    hasSavedSolution = true;
    savedBoxVectors[0] = kernel.boxVectors[0];
    savedBoxVectors[1] = kernel.boxVectors[1];
    savedBoxVectors[2] = kernel.boxVectors[2];
    kernel.cc.getPosq().copyTo(savedPositions);
    kernel.electrodeCharges.copyTo(savedCharges);
}

void CommonConstantPotentialSolver::getGuessChargeArrays(vector<ComputeArray*>& arrays) {
    arrays.clear();
}

CommonConstantPotentialMatrixSolver::CommonConstantPotentialMatrixSolver(ComputeContext& cc, int numParticles, int numElectrodeParticles, int paddedProblemSize) : CommonConstantPotentialSolver(cc, numParticles, numElectrodeParticles, paddedProblemSize) {
    int elementSize = cc.getUseDoublePrecision() ? sizeof(double) : sizeof(float);
    electrodePosData.initialize(cc, numElectrodeParticles, cc.getPosq().getElementSize(), "electrodePosData");
    capacitance.initialize(cc, (size_t) paddedProblemSize * paddedProblemSize, elementSize, "capacitance");
    constraintVector.initialize(cc, numElectrodeParticles, elementSize, "constraintVector");
    checkSavedElectrodePositionsKernelResult.initialize<int>(cc, 1, "checkSavedElectrodePositionsKernelResult");
}

void CommonConstantPotentialMatrixSolver::compileKernels(CommonCalcConstantPotentialForceKernel& kernel) {
    CommonConstantPotentialSolver::compileKernels(kernel);

    map<string, string> defines;
    defines["NUM_PARTICLES"] = kernel.cc.intToString(numParticles);
    defines["NUM_ELECTRODE_PARTICLES"] = kernel.cc.intToString(numElectrodeParticles);
    defines["THREAD_BLOCK_SIZE"] = kernel.cc.intToString(kernel.maxThreadBlockSize);
    defines["CHUNK_SIZE"] = kernel.cc.intToString(kernel.chunkSize);
    defines["CHUNK_COUNT"] = kernel.cc.intToString(kernel.chunkCount);
    defines["PADDED_PROBLEM_SIZE"] = kernel.cc.intToString(kernel.paddedProblemSize);
    if (kernel.useChargeConstraint) {
        defines["USE_CHARGE_CONSTRAINT"] = "1";
    }
    ComputeProgram program = kernel.cc.compileProgram(CommonKernelSources::constantPotentialMatrixSolver, defines);

    checkSavedElectrodePositionsKernel = program->createKernel("checkSavedElectrodePositions");
    checkSavedElectrodePositionsKernel->addArg(kernel.cc.getPosq());
    checkSavedElectrodePositionsKernel->addArg(electrodePosData);
    checkSavedElectrodePositionsKernel->addArg(kernel.elecToSys);
    checkSavedElectrodePositionsKernel->addArg(checkSavedElectrodePositionsKernelResult);

    saveElectrodePositionsKernel = program->createKernel("saveElectrodePositions");
    saveElectrodePositionsKernel->addArg(kernel.cc.getPosq());
    saveElectrodePositionsKernel->addArg(electrodePosData);
    saveElectrodePositionsKernel->addArg(kernel.elecToSys);

    solveKernel = program->createKernel("solve");
    solveKernel->addArg(kernel.electrodeCharges);
    solveKernel->addArg(kernel.chargeDerivatives);
    solveKernel->addArg(capacitance);
    if (kernel.useChargeConstraint) {
        solveKernel->addArg(constraintVector);
        solveKernel->addArg(); // chargeTarget
    }
}

void CommonConstantPotentialMatrixSolver::solveImpl(CommonCalcConstantPotentialForceKernel& kernel) {
    ensureValid(kernel);

    // Zero electrode charges and get derivatives at zero charge.
    kernel.cc.clearBuffer(kernel.electrodeCharges);
    kernel.mustUpdateElectrodeCharges = true;
    kernel.doDerivatives();

    // Solve for electrode charges directly using capacitance matrix and
    // calculated derivatives.
    if (kernel.useChargeConstraint) {
        if (kernel.cc.getUseDoublePrecision()) {
            solveKernel->setArg(4, kernel.chargeTarget);
        }
        else {
            solveKernel->setArg(4, (float) kernel.chargeTarget);
        }
    }
    solveKernel->execute(kernel.maxThreadBlockSize, kernel.maxThreadBlockSize);
    kernel.mustUpdateElectrodeCharges = true;
}

void CommonConstantPotentialMatrixSolver::ensureValid(CommonCalcConstantPotentialForceKernel& kernel) {
    // Initializes or updates the precomputed capacitance matrix if this is its
    // first use or electrode parameters have changed since its initialization.

    // Check for changes to box vectors or electrode positions that might
    // invalidate a matrix that is currently marked valid.
    if (valid) {
        if (boxVectors[0] != kernel.boxVectors[0] || boxVectors[1] != kernel.boxVectors[1] || boxVectors[2] != kernel.boxVectors[2]) {
            valid = false;
        }
    }
    if (valid) {
        kernel.cc.clearBuffer(checkSavedElectrodePositionsKernelResult);
        checkSavedElectrodePositionsKernel->execute(numElectrodeParticles);
        int result;
        checkSavedElectrodePositionsKernelResult.download(&result);
        if (result) {
            valid = false;
        }
    }
    if (valid) {
        return;
    }

    // We must have a valid neighbor list before populating the matrix.
    kernel.ensureValidNeighborList();

    // Store the current box vectors and electrode positions before updating the
    // capacitance matrix.
    valid = true;
    boxVectors[0] = kernel.boxVectors[0];
    boxVectors[1] = kernel.boxVectors[1];
    boxVectors[2] = kernel.boxVectors[2];
    saveElectrodePositionsKernel->execute(numElectrodeParticles);

    TNT::Array2D<double> A(paddedProblemSize, paddedProblemSize);
    vector<double> dUdQ0(numElectrodeParticles);
    vector<double> dUdQ(numElectrodeParticles);

    // Get derivatives when all electrode charges are zeroed.
    kernel.cc.clearBuffer(kernel.electrodeCharges);
    kernel.mustUpdateElectrodeCharges = true;
    kernel.doDerivatives();
    kernel.chargeDerivatives.download(dUdQ0, true);

    vector<double> electrodeCharges(paddedProblemSize);
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        // Get derivatives when one electrode charge is set.
        electrodeCharges[ii] = 1.0;
        kernel.electrodeCharges.upload(electrodeCharges, true);
        kernel.mustUpdateElectrodeCharges = true;
        kernel.doDerivatives();
        kernel.chargeDerivatives.download(dUdQ, true);
        electrodeCharges[ii] = 0.0;

        // Set matrix elements, subtracting zero charge derivatives so that the
        // matrix will end up being the (charge-independent) Hessian.
        for (int jj = 0; jj < ii; jj++) {
            A[ii][jj] = A[jj][ii] = dUdQ[jj] - dUdQ0[jj];
        }
        A[ii][ii] = dUdQ[ii] - dUdQ0[ii];
        for (int jj = numElectrodeParticles; jj < paddedProblemSize; jj++) {
            A[ii][jj] = A[jj][ii] = 0.0;
        }
    }
    for (int ii = numElectrodeParticles; ii < paddedProblemSize; ii++) {
        for (int jj = numElectrodeParticles; jj < paddedProblemSize; jj++) {
            A[ii][jj] = ii == jj ? 1.0 : 0.0;
        }
    }

    // Compute Cholesky decomposition representation of the inverse.
    JAMA::Cholesky<double> choleskyInverse(A);
    if (!choleskyInverse.is_spd()) {
        throw OpenMMException("Electrode matrix not positive definite");
    }

    // Load the results of the Cholesky decomposition into a buffer.
    TNT::Array2D<double> choleskyLower = choleskyInverse.getL();
    vector<double> hostCapacitance(capacitance.getSize());
    size_t index = 0;
    for (int ii = 0; ii < paddedProblemSize; ii++) {
        double scale = 1.0 / choleskyLower[ii][ii];

        // Load the lower triangle.
        for (int jj = 0; jj < ii; jj++) {
            hostCapacitance[index++] = choleskyLower[ii][jj] * scale;
        }

        // Load the reciprocal of the diagonal element.
        hostCapacitance[index++] = scale;

        // Load the transpose of the lower triangle into the upper triangle.
        for (int jj = ii + 1; jj < paddedProblemSize; jj++) {
            hostCapacitance[index++] = choleskyLower[jj][ii] * scale;
        }
    }
    capacitance.upload(hostCapacitance, true);

    // Precompute the appropriate scaling vector to enforce constant total
    // charge if requested.  The vector is parallel to one obtained by solving
    // Aq = b for all q_i = 1 (ensuring the constrained charges will actually be
    // the correct constrained minimum of the quadratic form for the energy),
    // and is scaled so that adding it to a vector of charges increases the
    // total charge by 1 (making it easy to calculate the necessary offset).
    if (kernel.useChargeConstraint) {
        TNT::Array1D<double> solution = choleskyInverse.solve(TNT::Array1D<double>(paddedProblemSize, 1.0));
        vector<double> hostConstraintVector(static_cast<double*>(solution), static_cast<double*>(solution) + numElectrodeParticles);

        double constraintScaleInv = 0.0;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            constraintScaleInv += hostConstraintVector[ii];
        }
        double constraintScale = 1.0 / constraintScaleInv;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            hostConstraintVector[ii] *= constraintScale;
        }

        constraintVector.upload(hostConstraintVector, true);
    }
}

CommonConstantPotentialCGSolver::CommonConstantPotentialCGSolver(ComputeContext& cc, int numParticles, int numElectrodeParticles, int paddedProblemSize, bool precond) :
        CommonConstantPotentialSolver(cc, numParticles, numElectrodeParticles, paddedProblemSize), precondRequested(precond) {
    int elementSize = cc.getUseDoublePrecision() ? sizeof(double) : sizeof(float);
    threadBlockCount = cc.getNumThreadBlocks();
    threadBlockSize = min(cc.getMaxThreadBlockSize(), cc.computeThreadBlockSize(max(5 * elementSize, 4 * elementSize + (int) sizeof(double))));
    q.initialize(cc, paddedProblemSize, elementSize, "q");
    grad.initialize(cc, numElectrodeParticles, elementSize, "grad");
    projGrad.initialize(cc, numElectrodeParticles, elementSize, "projGrad");
    precGrad.initialize(cc, numElectrodeParticles, elementSize, "precGrad");
    qStep.initialize(cc, paddedProblemSize, elementSize, "qStep");
    gradStep.initialize(cc, numElectrodeParticles, elementSize, "gradStep");
    grad0.initialize(cc, numElectrodeParticles, elementSize, "grad0");
    qLast.initialize(cc, paddedProblemSize, elementSize, "qLast");
    blockSums1.initialize(cc, threadBlockCount, 5 * elementSize, "blockSums1");
    blockSums2.initialize(cc, 1, 3 * elementSize, "blockSums2");
    convergedResult.initialize(cc, 1, sizeof(int), "convergedResult");
    if (precondRequested) {
        // If double precision is supported, this will hold double values;
        // otherwise, it will hold mm_float2 values.
        precondVector.initialize(cc, numElectrodeParticles, sizeof(double), "precondVector");
    }

    convergedDownloadStartEvent = cc.createEvent();
    convergedDownloadFinishEvent = cc.createEvent();
    convergedDownloadQueue = cc.createQueue();

    cc.clearBuffer(qLast);
}

void CommonConstantPotentialCGSolver::compileKernels(CommonCalcConstantPotentialForceKernel& kernel) {
    CommonConstantPotentialSolver::compileKernels(kernel);

    map<string, string> defines;
    defines["ERROR_TARGET"] = kernel.cc.doubleToString(kernel.cgErrorTol * kernel.cgErrorTol * numElectrodeParticles);
    defines["NUM_ELECTRODE_PARTICLES"] = kernel.cc.intToString(numElectrodeParticles);
    defines["THREAD_BLOCK_SIZE"] = kernel.cc.intToString(threadBlockSize);
    defines["THREAD_BLOCK_COUNT"] = kernel.cc.intToString(threadBlockCount);
    if (kernel.useChargeConstraint) {
        defines["USE_CHARGE_CONSTRAINT"] = "1";
    }
    if (precondRequested) {
        defines["PRECOND_REQUESTED"] = "1";
    }
    ComputeProgram program = kernel.cc.compileProgram(CommonKernelSources::constantPotentialCGSolver, defines);

    solveInitializeStep1Kernel = program->createKernel("solveInitializeStep1");
    solveInitializeStep1Kernel->addArg(kernel.electrodeCharges);
    solveInitializeStep1Kernel->addArg(qLast);
    if (kernel.useChargeConstraint) {
        solveInitializeStep1Kernel->addArg(); // chargeTarget
    }

    solveInitializeStep2Kernel = program->createKernel("solveInitializeStep2");
    solveInitializeStep2Kernel->addArg(kernel.chargeDerivatives);
    solveInitializeStep2Kernel->addArg(grad);
    solveInitializeStep2Kernel->addArg(projGrad);
    solveInitializeStep2Kernel->addArg(convergedResult);

    solveInitializeStep3Kernel = program->createKernel("solveInitializeStep3");
    solveInitializeStep3Kernel->addArg(kernel.electrodeCharges);
    solveInitializeStep3Kernel->addArg(kernel.chargeDerivatives);
    solveInitializeStep3Kernel->addArg(grad);
    solveInitializeStep3Kernel->addArg(projGrad);
    solveInitializeStep3Kernel->addArg(precGrad);
    solveInitializeStep3Kernel->addArg(qStep);
    solveInitializeStep3Kernel->addArg(grad0);
    if (precondRequested) {
        solveInitializeStep3Kernel->addArg(precondVector);
        solveInitializeStep3Kernel->addArg(); // precondActivated
    }

    solveLoopStep1Kernel = program->createKernel("solveLoopStep1");
    solveLoopStep1Kernel->addArg(kernel.chargeDerivatives);
    solveLoopStep1Kernel->addArg(q);
    solveLoopStep1Kernel->addArg(grad);
    solveLoopStep1Kernel->addArg(qStep);
    solveLoopStep1Kernel->addArg(gradStep);
    solveLoopStep1Kernel->addArg(grad0);
    solveLoopStep1Kernel->addArg(blockSums1);

    solveLoopStep2Kernel = program->createKernel("solveLoopStep2");
    solveLoopStep2Kernel->addArg(blockSums1);
    solveLoopStep2Kernel->addArg(convergedResult);

    solveLoopStep3Kernel = program->createKernel("solveLoopStep3");
    solveLoopStep3Kernel->addArg(q);
    solveLoopStep3Kernel->addArg(grad);
    solveLoopStep3Kernel->addArg(qStep);
    solveLoopStep3Kernel->addArg(gradStep);
    solveLoopStep3Kernel->addArg(blockSums1);
    solveLoopStep3Kernel->addArg(convergedResult);
    if (kernel.useChargeConstraint) {
        solveLoopStep3Kernel->addArg(); // chargeTarget
    }

    solveLoopStep4Kernel = program->createKernel("solveLoopStep4");
    solveLoopStep4Kernel->addArg(grad);
    solveLoopStep4Kernel->addArg(projGrad);
    solveLoopStep4Kernel->addArg(precGrad);
    solveLoopStep4Kernel->addArg(gradStep);
    solveLoopStep4Kernel->addArg(blockSums2);
    solveLoopStep4Kernel->addArg(convergedResult);
    if (precondRequested) {
        solveLoopStep4Kernel->addArg(precondVector);
        solveLoopStep4Kernel->addArg(); // precondActivated
    }

    solveLoopStep5Kernel = program->createKernel("solveLoopStep5");
    solveLoopStep5Kernel->addArg(kernel.electrodeCharges);
    solveLoopStep5Kernel->addArg(precGrad);
    solveLoopStep5Kernel->addArg(qStep);
    solveLoopStep5Kernel->addArg(blockSums1);
    solveLoopStep5Kernel->addArg(blockSums2);
    solveLoopStep5Kernel->addArg(convergedResult);
}

void CommonConstantPotentialCGSolver::solveImpl(CommonCalcConstantPotentialForceKernel& kernel) {
    ensureValid(kernel);

    if (kernel.useChargeConstraint) {
        if (kernel.cc.getUseDoublePrecision()) {
            solveInitializeStep1Kernel->setArg(2, kernel.chargeTarget);
            solveLoopStep3Kernel->setArg(6, kernel.chargeTarget);
        }
        else {
            solveInitializeStep1Kernel->setArg(2, (float) kernel.chargeTarget);
            solveLoopStep3Kernel->setArg(6, (float) kernel.chargeTarget);
        }
    }
    if (precondRequested) {
        solveInitializeStep3Kernel->setArg(8, (int) precondActivated);
        solveLoopStep4Kernel->setArg(7, (int) precondActivated);
    }

    int converged;

    // Evaluate the initial gradient Aq - b.
    solveInitializeStep1Kernel->execute(threadBlockSize, threadBlockSize);
    kernel.mustUpdateElectrodeCharges = true;
    kernel.doDerivatives();

    // Check for convergence at the initial guess charges.
    solveInitializeStep2Kernel->execute(threadBlockSize, threadBlockSize);
    convergedResult.download(&converged);
    if (converged) {
        return;
    }

    // Save the current charges, then evaluate the gradient with zero charges
    // (-b) so that we can later compute Ap as (Ap - b) - (-b).
    kernel.electrodeCharges.copyTo(q);
    kernel.cc.clearBuffer(kernel.electrodeCharges);
    kernel.mustUpdateElectrodeCharges = true;
    kernel.doDerivatives();

    // Prepare for conjugate gradient iterations.
    solveInitializeStep3Kernel->execute(threadBlockSize, threadBlockSize);
    kernel.cc.clearBuffer(blockSums1);

    // Perform conjugate gradient iterations.
    int* convergedPinned = (int*) kernel.cc.getPinnedBuffer();
    for (int iter = 0; iter <= numElectrodeParticles; iter++) {
        // Evaluate the matrix-vector product A qStep.
        kernel.mustUpdateElectrodeCharges = true;
        if (iter > 0) {
            // This is a subsequent iteration; check for convergence.  We can
            // start clearing buffers for the derivatives while we wait for the
            // flag to finish being copied from the device back to the host.
            kernel.initDoDerivatives();
            kernel.initPmeExecute();
            convergedDownloadFinishEvent->wait();
            converged = *convergedPinned;
            if (converged) {
                break;
            }
            kernel.doDerivatives(false);
        } else {
            // This is the first iteration; evaluate derivatives normally.
            kernel.doDerivatives();
        }

        solveLoopStep1Kernel->execute(numElectrodeParticles, threadBlockSize);
        solveLoopStep2Kernel->execute(threadBlockSize, threadBlockSize);
        solveLoopStep3Kernel->execute(numElectrodeParticles, threadBlockSize);

        // Periodically recompute the gradient vector instead of updating it to
        // reduce the accumulation of roundoff error.
        if (iter != 0 && iter % 32 == 0) {
            kernel.mustUpdateElectrodeCharges = true;
            kernel.doDerivatives();
            kernel.chargeDerivatives.copyTo(grad);
        }

        solveLoopStep4Kernel->execute(threadBlockSize, threadBlockSize);

        convergedDownloadStartEvent->enqueue();
        kernel.cc.setCurrentQueue(convergedDownloadQueue);
        convergedDownloadStartEvent->queueWait(convergedDownloadQueue);
        convergedResult.download(convergedPinned, false);
        convergedDownloadFinishEvent->enqueue();
        kernel.cc.restoreDefaultQueue();

        solveLoopStep5Kernel->execute(numElectrodeParticles, threadBlockSize);

        if (iter == numElectrodeParticles) {
            // No more iterations are allowed: download the convergence flag.
            convergedDownloadFinishEvent->wait();
            converged = *convergedPinned;
        }
    }

    if (!converged) {
        throw OpenMMException("Constant potential conjugate gradient iterations not converged");
    }

    // Store the final charges.
    q.copyTo(kernel.electrodeCharges);
    kernel.mustUpdateElectrodeCharges = true;
}

void CommonConstantPotentialCGSolver::getGuessChargeArrays(vector<ComputeArray*>& arrays) {
    CommonConstantPotentialSolver::getGuessChargeArrays(arrays);
    arrays.push_back(&qLast);
}

void CommonConstantPotentialCGSolver::ensureValid(CommonCalcConstantPotentialForceKernel& kernel) {
    // Initializes or updates information for a preconditioner for the conjugate
    // gradient method if this is its first use or electrode parameters have
    // changed since its initialization.

    // No action is required if the box vectors have not changed.
    if (valid && boxVectors[0] == kernel.boxVectors[0] && boxVectors[1] == kernel.boxVectors[1] && boxVectors[2] == kernel.boxVectors[2]) {
        return;
    }

    valid = true;
    boxVectors[0] = kernel.boxVectors[0];
    boxVectors[1] = kernel.boxVectors[1];
    boxVectors[2] = kernel.boxVectors[2];

    precondActivated = false;
    if (precondRequested) {
        // If electrode self-energy contributions differ between electrodes, a
        // preconditioner may help convergence; otherwise, it provides no
        // benefit and may slow convergence due to roundoff error.
        for (int ie = 1; ie < kernel.numElectrodes; ie++) {
            // Note hostElectrodeParams[1].w has the scale for electrode 0, etc.
            if (kernel.hostElectrodeParams[ie + 1].w != kernel.hostElectrodeParams[ie].w) {
                precondActivated = true;
                break;
            }
        }
    }

    float pmeTerm = 0.0f;

    if (precondActivated) {
        // Save all positions and charges.
        vector<mm_double4> posqSave(kernel.cc.getPaddedNumAtoms());
        vector<double> qSave;
        kernel.cc.getPosq().download(posqSave, true);
        if (!kernel.usePosqCharges) {
            qSave.resize(kernel.cc.getPaddedNumAtoms());
            kernel.charges.download(qSave, true);
        }

        vector<mm_double4> posqCopy(posqSave);
        vector<double> qCopy(qSave);

        // Zero all charges.
        if (kernel.usePosqCharges) {
            for (int i = 0; i < numParticles; i++) {
                posqCopy[i].w = 0.0;
            }
        }
        else {
            for (int i = 0; i < numParticles; i++) {
                qCopy[i] = 0.0;
            }
        }

        // Place a unit charge at the origin.
        int i0 = kernel.hostElecToSys[0];
        posqCopy[i0] = mm_double4(0.0, 0.0, 0.0, 1.0);
        if (!kernel.usePosqCharges) {
            qCopy[i0] = 1.0;
        }

        kernel.cc.getPosq().upload(posqCopy, true);
        if (!kernel.usePosqCharges) {
            kernel.charges.upload(qCopy, true);
        }

        // Perform a reference PME calculation with a single charge at the
        // origin to find the constant offset on the preconditioner diagonal due
        // to the PME calculation.  This will actually vary slightly with
        // position but only due to finite accuracy of the PME splines, so it is
        // fine to assume it will be constant for the preconditioner.
        kernel.cc.clearBuffer(kernel.chargeDerivativesFixed);
        kernel.pmeShouldSort = true;
        kernel.pmeExecute(false, false, true);
        kernel.pmeShouldSort = true;
        vector<long> derivatives(numElectrodeParticles);
        kernel.chargeDerivativesFixed.download(derivatives);
        double pmeTerm = derivatives[0] / (double) 0x100000000;

        // Restore all positions and charges.
        kernel.cc.getPosq().upload(posqSave, true);
        if (!kernel.usePosqCharges) {
            kernel.charges.upload(qSave, true);
        }

        // The diagonal has a contribution from reciprocal space, Ewald
        // self-interaction, Ewald neutralizing plasma, Gaussian
        // self-interaction, and Thomas-Fermi contributions.
        double plasmaScale = CommonCalcConstantPotentialForceKernel::PLASMA_SCALE / (boxVectors[0][0] * boxVectors[1][1] * boxVectors[2][2] * kernel.ewaldAlpha * kernel.ewaldAlpha);
        vector<double> hostPrecondVector(numElectrodeParticles);
        double precondScaleInv = 0.0;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            hostPrecondVector[ii] = 1.0 / (2.0 * (kernel.hostElectrodeParams[kernel.hostElecElec[ii] + 1].w - plasmaScale) + pmeTerm);
            precondScaleInv += hostPrecondVector[ii];
        }
        double precondScale = 1.0 / precondScaleInv;
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            hostPrecondVector[ii] *= precondScale;
        }
        if (kernel.cc.getSupportsDoublePrecision()) {
            precondVector.upload(hostPrecondVector);
        }
        else {
            vector<mm_float2> hostPrecondVectorSplit(numElectrodeParticles);
            for (int ii = 0; ii < numElectrodeParticles; ii++) {
                float const x = (float) hostPrecondVector[ii];
                hostPrecondVectorSplit[ii].x = x;
                hostPrecondVectorSplit[ii].y = (float) (hostPrecondVector[ii] - (double) x);
            }
            precondVector.upload(hostPrecondVectorSplit);
        }
    }
}

const double CommonCalcConstantPotentialForceKernel::SELF_ALPHA_SCALE = ONE_4PI_EPS0 / sqrt(PI_M);
const double CommonCalcConstantPotentialForceKernel::SELF_ETA_SCALE = ONE_4PI_EPS0 / sqrt(2.0 * PI_M);
const double CommonCalcConstantPotentialForceKernel::SELF_TF_SCALE = 1.0 / (2.0 * EPSILON0);
const double CommonCalcConstantPotentialForceKernel::PLASMA_SCALE = 1.0 / (8.0 * EPSILON0);

CommonCalcConstantPotentialForceKernel::~CommonCalcConstantPotentialForceKernel() {
    ContextSelector selector(cc);

    if (solver != NULL) {
        delete solver;
    }
}

void CommonCalcConstantPotentialForceKernel::commonInitialize(const System& system, const ConstantPotentialForce& force, bool deviceIsCpu, bool useFixedPointChargeSpreading) {
    ContextSelector selector(cc);
    if (cc.getNumContexts() > 1) {
        throw OpenMMException("ConstantPotentialForce does not support using multiple devices");
    }

    forceGroup = force.getForceGroup();

    maxThreadBlockSize = cc.getMaxThreadBlockSize();
    int elementSize = cc.getUseDoublePrecision() ? sizeof(double) : sizeof(float);

    this->deviceIsCpu = deviceIsCpu;
    this->useFixedPointChargeSpreading = useFixedPointChargeSpreading;
    this->usePosqCharges = cc.requestPosqCharges();

    int forceIndex;
    for (forceIndex = 0; forceIndex < system.getNumForces() && &system.getForce(forceIndex) != &force; ++forceIndex) {
    }
    string prefix = "constantPotential" + cc.intToString(forceIndex) + "_";

    // Get particle parameters.
    numParticles = force.getNumParticles();
    setCharges.resize(numParticles);
    for (int i = 0; i < numParticles; i++) {
        force.getParticleParameters(i, setCharges[i]);
    }

    // Get exceptions and identify "1-4" exceptions (those that don't zero the
    // charge product).
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd;
        force.getExceptionParameters(i, particle1, particle2, chargeProd);
        exclusions.push_back(pair<int, int>(particle1, particle2));
        if (chargeProd != 0.0) {
            exceptionIndex[i] = exceptions.size();
            exceptions.push_back(i);
        }
    }

    // Get a list of all exclusions per particle, including exclusions for
    // self-interaction.
    vector<vector<int> > exclusionList(numParticles);
    for (int i = 0; i < numParticles; i++) {
        exclusionList[i].push_back(i);
    }
    for (auto exclusion : exclusions) {
        exclusionList[exclusion.first].push_back(exclusion.second);
        exclusionList[exclusion.second].push_back(exclusion.first);
    }

    // Get nonbonded parameters.
    cutoff = force.getCutoffDistance();
    ConstantPotentialForceImpl::calcPMEParameters(system, force, ewaldAlpha, gridSizeX, gridSizeY, gridSizeZ);
    gridSizeX = cc.findLegalFFTDimension(gridSizeX);
    gridSizeY = cc.findLegalFFTDimension(gridSizeY);
    gridSizeZ = cc.findLegalFFTDimension(gridSizeZ);

    // Get electrode parameters.  sysToElec will be a map from system particle
    // indices to electrode particle indices (or -1 if the particle is not an
    // electrode particle), while elecToSys will be a map from electrode
    // particle indices to system particle indices.  sysElec will be a map from
    // system particle indices to electrode indices (or -1 if the particle is
    // not an electrode particle), while elecElec will be a map from electrode
    // particle indices to electrode indices.  Precompute and store electrode
    // parameters for electrodes [-1, numElectrodes - 1] as [0, numElectrodes]
    // in hostElectrodeParams for convenient lookup:
    //   x: potential
    //   y: gaussianWidth
    //   z: thomasFermiScale
    //   w: Precomputed self-energy scale
    numElectrodes = force.getNumElectrodes();
    hostSysToElec.resize(cc.getPaddedNumAtoms(), -1);
    hostSysElec.resize(cc.getPaddedNumAtoms(), -1);
    hostElectrodeParams.resize(numElectrodes + 1);
    for (int ie = -1; ie < numElectrodes; ie++) {
        double potential = 0.0;
        double gaussianWidth = 0.0;
        double thomasFermiScale = 0.0;
        double selfScale = -SELF_ALPHA_SCALE * ewaldAlpha;
        if (ie != -1) {
            set<int> electrodeParticles;
            force.getElectrodeParameters(ie, electrodeParticles, potential, gaussianWidth, thomasFermiScale);
            for (int i : electrodeParticles) {
                hostSysToElec[i] = hostElecToSys.size();
                hostSysElec[i] = ie;
                hostElecToSys.push_back(i);
                hostElecElec.push_back(ie);
            }
            selfScale += SELF_ETA_SCALE / gaussianWidth + SELF_TF_SCALE * thomasFermiScale;
        }
        hostElectrodeParams[ie + 1] = mm_double4(potential, gaussianWidth, thomasFermiScale, selfScale);
    }
    numElectrodeParticles = hostElecToSys.size();
    hasElectrodes = (numElectrodeParticles != 0);

    chunkSize = deviceIsCpu ? 1 : 32;
    chunkCount = (numElectrodeParticles + chunkSize - 1) / chunkSize;
    paddedProblemSize = chunkCount * chunkSize;

    // Clear charges on electrode particles.
    for (int ii = 0; ii < numElectrodeParticles; ii++) {
        setCharges[hostElecToSys[ii]] = 0.0;
    }

    // If using a charge constraint, set the charge target to be that on the
    // electrode particles (so that the overall charge is constrained correctly
    // if the non-electrolyte particles are non-neutral).
    useChargeConstraint = force.getUseChargeConstraint();
    if (useChargeConstraint) {
        chargeTarget = force.getChargeConstraintTarget();
        for (int i = 0; i < numParticles; i++) {
            chargeTarget -= setCharges[i];
        }
    }

    // Save external field.
    force.getExternalField(externalField);

    // Set up PME.
    pmeSetup();

    // Create bonded force for exclusion parameters.
    int numExclusions = force.getNumExceptions();
    if (numExclusions > 0) {
        exclusionScales.initialize(cc, numExclusions, elementSize, "exclusionScales");
        vector<vector<int> > hostExclusionAtoms(numExclusions, vector<int>(2));
        vector<double> hostExclusionScales(numExclusions);
        for (int i = 0; i < numExclusions; i++) {
            hostExclusionAtoms[i][0] = exclusions[i].first;
            hostExclusionAtoms[i][1] = exclusions[i].second;
            hostExclusionScales[i] = ONE_4PI_EPS0 * setCharges[exclusions[i].first] * setCharges[exclusions[i].second];
        }
        exclusionScales.upload(hostExclusionScales, true);

        map<string, string> replacements;
        replacements["APPLY_PERIODIC"] = force.getExceptionsUsePeriodicBoundaryConditions() ? "1" : "0";
        replacements["EWALD_ALPHA"] = cc.doubleToString(ewaldAlpha);
        replacements["PARAMS"] = cc.getBondedUtilities().addArgument(exclusionScales, "real");
        replacements["TWO_OVER_SQRT_PI"] = cc.doubleToString(2.0 / sqrt(M_PI));
        cc.getBondedUtilities().addInteraction(hostExclusionAtoms, cc.replaceStrings(CommonKernelSources::constantPotentialExclusions, replacements), forceGroup);
    }

    // Upload 1-4 exception parameters and create bonded force.
    int numExceptions = exceptions.size();
    if (numExceptions > 0) {
        exceptionScales.initialize(cc, numExceptions, elementSize, "exceptionScales");
        vector<vector<int> > hostExceptionAtoms(numExceptions, vector<int>(2));
        vector<double> hostExceptionScales(numExceptions);
        for (int i = 0; i < numExceptions; i++) {
            double chargeProd;
            force.getExceptionParameters(exceptions[i], hostExceptionAtoms[i][0], hostExceptionAtoms[i][1], chargeProd);
            hostExceptionScales[i] = ONE_4PI_EPS0 * chargeProd;
        }
        exceptionScales.upload(hostExceptionScales, true);

        map<string, string> replacements;
        replacements["APPLY_PERIODIC"] = force.getExceptionsUsePeriodicBoundaryConditions() ? "1" : "0";
        replacements["PARAMS"] = cc.getBondedUtilities().addArgument(exceptionScales, "real");
        cc.getBondedUtilities().addInteraction(hostExceptionAtoms, cc.replaceStrings(CommonKernelSources::constantPotentialExceptions, replacements), forceGroup);
    }

    // Upload parameters for electrodes.
    electrodeParams.initialize(cc, numElectrodes + 1, 4 * elementSize, "electrodeParams");
    electrodeParams.upload(hostElectrodeParams, true);

    // Initialize lookup tables for constant potential on the device.
    sysToElec.initialize<int>(cc, cc.getPaddedNumAtoms(), "sysToElec");
    sysElec.initialize<int>(cc, cc.getPaddedNumAtoms(), "sysElec");
    sysToElec.upload(hostSysToElec);
    sysElec.upload(hostSysElec);
    if (hasElectrodes) {
        elecToSys.initialize<int>(cc, numElectrodeParticles, "elecToSys");
        elecElec.initialize<int>(cc, numElectrodeParticles, "elecElec");
        elecToSys.upload(hostElecToSys);
        elecElec.upload(hostElecElec);
    }

    // Initialize nonbonded force.
    map<string, string> nonbondedReplacements;
    nonbondedReplacements["EWALD_ALPHA"] = cc.doubleToString(ewaldAlpha);
    nonbondedReplacements["ONE_4PI_EPS0"] = cc.doubleToString(ONE_4PI_EPS0);
    nonbondedReplacements["TWO_OVER_SQRT_PI"] = cc.doubleToString(2.0 / sqrt(M_PI));
    if (usePosqCharges) {
        nonbondedReplacements["CHARGE1"] = "posq1.w";
        nonbondedReplacements["CHARGE2"] = "posq2.w";
    }
    else {
        nonbondedReplacements["CHARGE1"] = prefix + "charge1";
        nonbondedReplacements["CHARGE2"] = prefix + "charge2";
        cc.getNonbondedUtilities().addParameter(ComputeParameterInfo(charges, prefix + "charge", "real", 1));
    }
    nonbondedReplacements["SYSELEC1"] = prefix + "sysElec1";
    nonbondedReplacements["SYSELEC2"] = prefix + "sysElec2";
    cc.getNonbondedUtilities().addParameter(ComputeParameterInfo(sysElec, prefix + "sysElec", "int", 1));
    nonbondedReplacements["PARAMS"] = prefix + "params";
    cc.getNonbondedUtilities().addArgument(ComputeParameterInfo(electrodeParams, prefix + "params", "real", 4));
    cc.getNonbondedUtilities().addInteraction(true, true, true, force.getCutoffDistance(), exclusionList, cc.replaceStrings(CommonKernelSources::constantPotentialCoulombEnergyForces, nonbondedReplacements), force.getForceGroup(), true, false);

    // Initialize the constant potential solver.
    method = force.getConstantPotentialMethod();
    cgErrorTol = force.getCGErrorTolerance();
    if (method == ConstantPotentialForce::Matrix) {
        if (hasElectrodes) {
            solver = new CommonConstantPotentialMatrixSolver(cc, numParticles, numElectrodeParticles, paddedProblemSize);
        }
    }
    else if (method == ConstantPotentialForce::CG) {
        if (hasElectrodes) {
            solver = new CommonConstantPotentialCGSolver(cc, numParticles, numElectrodeParticles, paddedProblemSize, force.getUsePreconditioner());
        }
    }
    else {
        throw OpenMMException("internal error: invalid constant potential method");
    }

    // Upload fixed charges and initial guesses for electrode charges.
    hostNonElectrodeCharges = setCharges;
    hostNonElectrodeCharges.resize(cc.getPaddedNumAtoms());
    nonElectrodeCharges.initialize(cc, cc.getPaddedNumAtoms(), elementSize, "nonElectrodeCharges");
    nonElectrodeCharges.upload(hostNonElectrodeCharges, true);
    if (hasElectrodes) {
        hostElectrodeCharges.resize(paddedProblemSize);
        electrodeCharges.initialize(cc, paddedProblemSize, elementSize, "electrodeCharges");
        electrodeCharges.upload(hostElectrodeCharges, true);
        chargeDerivatives.initialize(cc, numElectrodeParticles, elementSize, "chargeDerivatives");
        chargeDerivativesFixed.initialize<int64_t>(cc, numElectrodeParticles, "chargeDerivativesFixed");
    }

    // nonElectrodeCharges holds all fixed charges, and all electrode charges in
    // nonElectrodeCharges should be zeroed.  electrodeCharges holds all current
    // electrode charges.  charges holds charges to actually use for evaluation
    // if usePosqCharges is false; otherwise, those charges will be in posq.
    // The mustUpdate flags indicate that something has changed and charges or
    // posq must be updated before energy/force/derivative evaluation.
    charges.initialize(cc, cc.getPaddedNumAtoms(), elementSize, "charges");
    mustUpdateNonElectrodeCharges = true;
    mustUpdateElectrodeCharges = hasElectrodes;

    totalChargeBuffer.initialize(cc, 1, elementSize, "totalChargeBuffer");

    info = new ForceInfo(force, hostSysElec, hostElectrodeParams);
    cc.addForce(info);

    // Initialize cell offsets for finite field computation.
    posCellOffsets.initialize<mm_int4>(cc, cc.getPaddedNumAtoms(), "posCellOffsets");
    hostPosCellOffsets = cc.getPosCellOffsets();
    posCellOffsets.upload(hostPosCellOffsets);

    // Create a reorder listener to swap electrode guess charges.
    ReorderListener* listener = new ReorderListener(cc, numElectrodeParticles, paddedProblemSize, hostSysToElec, hostElecToSys);
    if (hasElectrodes) {
        vector<ComputeArray*> guessChargeArrays;
        solver->getGuessChargeArrays(guessChargeArrays);
        for (ComputeArray* guessChargeArray : guessChargeArrays) {
            listener->addChargeArray(*guessChargeArray);
        }
    }
    cc.addReorderListener(listener);

    // Create a post computation object to avoid caching solutions if forces
    // are invalid (e.g., due to the neighbor list needing to be resized).
    if (hasElectrodes) {
        cc.addPostComputation(new InvalidatePostComputation(cc, solver));
    }
}

void CommonCalcConstantPotentialForceKernel::copyParametersToContext(ContextImpl& context, const ConstantPotentialForce& force, int firstParticle, int lastParticle, int firstException, int lastException, int firstElectrode, int lastElectrode) {
    ContextSelector selector(cc);

    // Get particle parameters.
    if (force.getNumParticles() != numParticles) {
        throw OpenMMException("updateParametersInContext: The number of particles has changed");
    }
    if (firstParticle <= lastParticle) {
        for (int i = firstParticle; i <= lastParticle; i++) {
            if (hostSysElec[i] == -1) {
                force.getParticleParameters(i, setCharges[i]);
                hostNonElectrodeCharges[i] = setCharges[i];
            }
        }
        if (cc.getUseDoublePrecision()) {
            nonElectrodeCharges.uploadSubArray(&hostNonElectrodeCharges[firstParticle], firstParticle, lastParticle - firstParticle + 1);
        }
        else {
            vector<float> hostNonElectrodeChargesFloat(lastParticle - firstParticle + 1);
            for (int i = firstParticle; i <= lastParticle; i++) {
                hostNonElectrodeChargesFloat[i - firstParticle] = (float) hostNonElectrodeCharges[i];
            }
            nonElectrodeCharges.uploadSubArray(&hostNonElectrodeChargesFloat[0], firstParticle, lastParticle - firstParticle + 1);
        }
        mustUpdateNonElectrodeCharges = true;
    }

    // Check exceptions.
    int numExclusions = force.getNumExceptions();
    int numExceptions = exceptions.size();

    vector<int> checkExceptions;
    for (int i = 0; i < numExclusions; i++) {
        int particle1, particle2;
        double chargeProd;
        force.getExceptionParameters(i, particle1, particle2, chargeProd);
        if (exceptionIndex.find(i) == exceptionIndex.end()) {
            if (chargeProd != 0.0) {
                throw OpenMMException("updateParametersInContext: The set of non-excluded exceptions has changed");
            }
        }
        else {
            checkExceptions.push_back(i);
        }
    }
    if (checkExceptions.size() != numExceptions) {
        throw OpenMMException("updateParametersInContext: The set of non-excluded exceptions has changed");
    }

    // Upload exclusion parameters.
    if (numExclusions > 0 && (firstParticle <= lastParticle || firstException <= lastException)) {
        vector<double> hostExclusionScales(numExclusions);
        for (int i = 0; i < numExclusions; i++) {
            hostExclusionScales[i] = ONE_4PI_EPS0 * setCharges[exclusions[i].first] * setCharges[exclusions[i].second];
        }
        exclusionScales.upload(hostExclusionScales, true);
    }

    // Upload exception parameters.
    if (numExceptions > 0 && firstException <= lastException) {
        vector<double> hostExceptionScales(numExceptions);
        for (int i = 0; i < numExceptions; i++) {
            int particle1, particle2;
            double chargeProd;
            force.getExceptionParameters(exceptions[i], particle1, particle2, chargeProd);
            hostExceptionScales[i] = ONE_4PI_EPS0 * chargeProd;
        }
        exceptionScales.upload(hostExceptionScales, true);
    }

    // Get electrode parameters.
    set<int> allElectrodeParticles;
    for (int ie = 0; ie < force.getNumElectrodes(); ie++) {
        set<int> electrodeParticles;
        double potential;
        double gaussianWidth;
        double thomasFermiScale;
        force.getElectrodeParameters(ie, electrodeParticles, potential, gaussianWidth, thomasFermiScale);
        for (int i : electrodeParticles) {
            if (hostSysElec[i] != ie) {
                throw OpenMMException("updateParametersInContext: The electrode assignment of a particle has changed");
            }
            allElectrodeParticles.insert(i);
        }
        if (ie >= firstElectrode && ie <= lastElectrode) {
            double selfScale = SELF_ETA_SCALE / gaussianWidth + SELF_TF_SCALE * thomasFermiScale - SELF_ALPHA_SCALE * ewaldAlpha;
            hostElectrodeParams[ie + 1] = mm_double4(potential, gaussianWidth, thomasFermiScale, selfScale);
        }
    }
    if (allElectrodeParticles.size() != numElectrodeParticles) {
        throw OpenMMException("updateParametersInContext: The electrode state of a particle has changed");
    }
    if (firstElectrode <= lastElectrode) {
        electrodeParams.upload(hostElectrodeParams, true);
    }
    
    // Update charge target.
    if (useChargeConstraint) {
        chargeTarget = force.getChargeConstraintTarget();
        for (int i = 0; i < numParticles; i++) {
            chargeTarget -= setCharges[i];
        }
    }

    // Update external field.
    force.getExternalField(externalField);

    cc.invalidateMolecules(info, firstParticle <= lastParticle || firstElectrode <= lastElectrode, firstException <= lastException);

    if (solver != NULL) {
        solver->discardSavedSolution();
        if (firstElectrode <= lastElectrode) {
            solver->invalidate();
        }
    }
}

void CommonCalcConstantPotentialForceKernel::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    alpha = ewaldAlpha;
    nx = gridSizeX;
    ny = gridSizeY;
    nz = gridSizeZ;
}

double CommonCalcConstantPotentialForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    ContextSelector selector(cc);

    ensureInitialized(context);

    cc.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    setKernelInputs(includeEnergy, includeForces);

    pmeShouldSort = true;

    if (solver != NULL) {
        solver->solve(*this);
    }

    return doEnergyForces(includeForces, includeEnergy);
}

void CommonCalcConstantPotentialForceKernel::getCharges(ContextImpl& context, vector<double>& chargesOut) {
    ContextSelector selector(cc);

    ensureInitialized(context);
    ensureValidNeighborList();

    cc.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    setKernelInputs(false, false);

    pmeShouldSort = true;

    if (solver != NULL) {
        solver->solve(*this);
    }

    // Preserve fixed charges exactly and load solved values of fluctuating
    // charges.
    chargesOut = setCharges;
    if (hasElectrodes) {
        const vector<int>& order = cc.getAtomIndex();
        electrodeCharges.download(hostElectrodeCharges, true);
        for (int ii = 0; ii < numElectrodeParticles; ii++) {
            chargesOut[order[hostElecToSys[ii]]] = hostElectrodeCharges[ii];
        }
    }
}

void CommonCalcConstantPotentialForceKernel::ensureInitialized(ContextImpl& context) {
    if (hasInitializedKernel) {
        return;
    }
    
    NonbondedUtilities& nb = cc.getNonbondedUtilities();

    map<string, string> defines;
    defines["CUTOFF"] = cc.doubleToString(cutoff);
    defines["CUTOFF_SQUARED"] = cc.doubleToString(cutoff * cutoff);
    defines["EWALD_ALPHA"] = cc.doubleToString(ewaldAlpha);
    defines["ONE_4PI_EPS0"] = cc.doubleToString(ONE_4PI_EPS0);
    defines["NUM_PARTICLES"] = cc.intToString(numParticles);
    defines["NUM_ELECTRODE_PARTICLES"] = cc.intToString(numElectrodeParticles);
    defines["NUM_EXCLUSION_TILES"] = cc.intToString(nb.getExclusionTiles().getSize());
    defines["PADDED_NUM_ATOMS"] = cc.intToString(cc.getPaddedNumAtoms());
    defines["PLASMA_SCALE"] = cc.doubleToString(PLASMA_SCALE);
    defines["THREAD_BLOCK_SIZE"] = cc.intToString(maxThreadBlockSize);
    defines["TILE_SIZE"] = cc.intToString(ComputeContext::TileSize);
    defines["WORK_GROUP_SIZE"] = cc.intToString(nb.getForceThreadBlockSize());
    if (usePosqCharges) {
        defines["USE_POSQ_CHARGES"] = "1";
    }
    if (deviceIsCpu) {
        defines["DEVICE_IS_CPU"] = "1";
    }
    ComputeProgram program = cc.compileProgram(CommonKernelSources::constantPotential, defines);

    updateNonElectrodeChargesKernel = program->createKernel("updateNonElectrodeCharges");
    updateNonElectrodeChargesKernel->addArg(cc.getPosq());
    updateNonElectrodeChargesKernel->addArg(charges);
    updateNonElectrodeChargesKernel->addArg(nonElectrodeCharges);
    updateNonElectrodeChargesKernel->addArg(sysElec);

    if (hasElectrodes) {
        updateElectrodeChargesKernel = program->createKernel("updateElectrodeCharges");
        updateElectrodeChargesKernel->addArg(cc.getPosq());
        updateElectrodeChargesKernel->addArg(charges);
        updateElectrodeChargesKernel->addArg(electrodeCharges);
        updateElectrodeChargesKernel->addArg(elecToSys);
    }

    getTotalChargeKernel = program->createKernel("getTotalCharge");
    getTotalChargeKernel->addArg(cc.getPosq());
    getTotalChargeKernel->addArg(charges);
    getTotalChargeKernel->addArg(totalChargeBuffer);

    evaluateSelfEnergyForcesKernel = program->createKernel("evaluateSelfEnergyForces");
    evaluateSelfEnergyForcesKernel->addArg(cc.getPosq());
    evaluateSelfEnergyForcesKernel->addArg(charges);
    evaluateSelfEnergyForcesKernel->addArg(sysElec);
    evaluateSelfEnergyForcesKernel->addArg(electrodeParams);
    evaluateSelfEnergyForcesKernel->addArg(totalChargeBuffer);
    evaluateSelfEnergyForcesKernel->addArg(posCellOffsets);
    evaluateSelfEnergyForcesKernel->addArg(); // periodicBoxVecX
    evaluateSelfEnergyForcesKernel->addArg(); // periodicBoxVecY
    evaluateSelfEnergyForcesKernel->addArg(); // periodicBoxVecZ
    evaluateSelfEnergyForcesKernel->addArg(); // externalField
    evaluateSelfEnergyForcesKernel->addArg(cc.getEnergyBuffer());
    evaluateSelfEnergyForcesKernel->addArg(cc.getLongForceBuffer());

    if (hasElectrodes) {
        evaluateDirectDerivativesKernel = program->createKernel("evaluateDirectDerivatives");
        evaluateDirectDerivativesKernel->addArg(cc.getPosq());
        evaluateDirectDerivativesKernel->addArg(charges);
        evaluateDirectDerivativesKernel->addArg(sysToElec);
        evaluateDirectDerivativesKernel->addArg(sysElec);
        evaluateDirectDerivativesKernel->addArg(electrodeParams);
        evaluateDirectDerivativesKernel->addArg(chargeDerivativesFixed);
        evaluateDirectDerivativesKernel->addArg(); // periodicBoxSize
        evaluateDirectDerivativesKernel->addArg(); // invPeriodicBoxSize
        evaluateDirectDerivativesKernel->addArg(); // periodicBoxVecX
        evaluateDirectDerivativesKernel->addArg(); // periodicBoxVecY
        evaluateDirectDerivativesKernel->addArg(); // periodicBoxVecZ
        evaluateDirectDerivativesKernel->addArg(nb.getExclusionTiles());
        evaluateDirectDerivativesKernel->addArg(nb.getInteractingTiles());
        evaluateDirectDerivativesKernel->addArg(nb.getInteractionCount());
        evaluateDirectDerivativesKernel->addArg(nb.getBlockCenters());
        evaluateDirectDerivativesKernel->addArg(nb.getBlockBoundingBoxes());
        evaluateDirectDerivativesKernel->addArg(nb.getInteractingAtoms());
        evaluateDirectDerivativesKernel->addArg(); // maxTiles

        finishDerivativesKernel = program->createKernel("finishDerivatives");
        finishDerivativesKernel->addArg(cc.getPosq());
        finishDerivativesKernel->addArg(charges);
        finishDerivativesKernel->addArg(elecToSys);
        finishDerivativesKernel->addArg(elecElec);
        finishDerivativesKernel->addArg(electrodeParams);
        finishDerivativesKernel->addArg(totalChargeBuffer);
        finishDerivativesKernel->addArg(posCellOffsets);
        finishDerivativesKernel->addArg(); // periodicBoxVecX
        finishDerivativesKernel->addArg(); // periodicBoxVecY
        finishDerivativesKernel->addArg(); // periodicBoxVecZ
        finishDerivativesKernel->addArg(); // externalField
        finishDerivativesKernel->addArg(chargeDerivatives);
        finishDerivativesKernel->addArg(chargeDerivativesFixed);
    }

    pmeCompileKernels();
    if (solver != NULL) {
        solver->compileKernels(*this);
    }

    hasInitializedKernel = true;
}

double CommonCalcConstantPotentialForceKernel::doEnergyForces(bool includeForces, bool includeEnergy) {
    if (mustUpdateNonElectrodeCharges) {
        updateNonElectrodeChargesKernel->execute(numParticles);
        mustUpdateNonElectrodeCharges = false;
    }
    if (mustUpdateElectrodeCharges) {
        if (hasElectrodes) {
            updateElectrodeChargesKernel->execute(numElectrodeParticles);
        }
        mustUpdateElectrodeCharges = false;
    }

    pmeExecute(includeEnergy, includeForces, false);

    // Ewald neutralizing plasma and per-particle energy.
    if (includeEnergy || includeForces) {
        getTotalChargeKernel->execute(maxThreadBlockSize, maxThreadBlockSize);
        if (cc.getUseDoublePrecision()) {
            evaluateSelfEnergyForcesKernel->setArg(9, mm_double4(externalField[0], externalField[1], externalField[2], 0));
        }
        else {
            evaluateSelfEnergyForcesKernel->setArg(9, mm_float4((float) externalField[0], (float) externalField[1], (float) externalField[2], 0));
        }
        evaluateSelfEnergyForcesKernel->execute(numParticles);
    }

    return 0.0;
}

void CommonCalcConstantPotentialForceKernel::initDoDerivatives() {
    if (mustUpdateNonElectrodeCharges) {
        updateNonElectrodeChargesKernel->execute(numParticles);
        mustUpdateNonElectrodeCharges = false;
    }
    if (mustUpdateElectrodeCharges) {
        if (hasElectrodes) {
            updateElectrodeChargesKernel->execute(numElectrodeParticles);
        }
        mustUpdateElectrodeCharges = false;
    }

    cc.clearBuffer(chargeDerivatives);
    cc.clearBuffer(chargeDerivativesFixed);
}

void CommonCalcConstantPotentialForceKernel::doDerivatives(bool init) {
    if (init) {
        initDoDerivatives();
    }

    pmeExecute(false, false, true, init);

    NonbondedUtilities& nb = cc.getNonbondedUtilities();
    evaluateDirectDerivativesKernel->execute(nb.getNumForceThreadBlocks() * nb.getForceThreadBlockSize(), nb.getForceThreadBlockSize());

    // Ewald neutralizing plasma and per-particle derivatives.
    getTotalChargeKernel->execute(maxThreadBlockSize, maxThreadBlockSize);
    finishDerivativesKernel->execute(numElectrodeParticles);
}

void CommonCalcConstantPotentialForceKernel::pmeSetup() {
    pmeDefines["PME_ORDER"] = cc.intToString(PmeOrder);
    pmeDefines["NUM_ATOMS"] = cc.intToString(numParticles);
    pmeDefines["PADDED_NUM_ATOMS"] = cc.intToString(cc.getPaddedNumAtoms());
    pmeDefines["NUM_INDICES"] = cc.intToString(numElectrodeParticles);
    pmeDefines["RECIP_EXP_FACTOR"] = cc.doubleToString(M_PI * M_PI / (ewaldAlpha * ewaldAlpha));
    pmeDefines["GRID_SIZE_X"] = cc.intToString(gridSizeX);
    pmeDefines["GRID_SIZE_Y"] = cc.intToString(gridSizeY);
    pmeDefines["GRID_SIZE_Z"] = cc.intToString(gridSizeZ);
    pmeDefines["EPSILON_FACTOR"] = cc.doubleToString(sqrt(ONE_4PI_EPS0));
    pmeDefines["M_PI"] = cc.doubleToString(M_PI);
    if (deviceIsCpu) {
        pmeDefines["DEVICE_IS_CPU"] = "1";
    }
    if (useFixedPointChargeSpreading) {
        pmeDefines["USE_FIXED_POINT_CHARGE_SPREADING"] = "1";
    }

    // Create required data structures.

    int elementSize = (cc.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
    int gridElements = gridSizeX*gridSizeY*gridSizeZ;
    pmeGrid1.initialize(cc, gridElements, 2*elementSize, "pmeGrid1");
    pmeGrid2.initialize(cc, gridElements, 2*elementSize, "pmeGrid2");
    pmeBsplineModuliX.initialize(cc, gridSizeX, elementSize, "pmeBsplineModuliX");
    pmeBsplineModuliY.initialize(cc, gridSizeY, elementSize, "pmeBsplineModuliY");
    pmeBsplineModuliZ.initialize(cc, gridSizeZ, elementSize, "pmeBsplineModuliZ");
    pmeAtomGridIndex.initialize<mm_int2>(cc, numParticles, "pmeAtomGridIndex");
    int energyElementSize = (cc.getUseDoublePrecision() || cc.getUseMixedPrecision() ? sizeof(double) : sizeof(float));
    sort = cc.createSort(new SortTrait(), cc.getNumAtoms());
    fft = cc.createFFT(gridSizeX, gridSizeY, gridSizeZ, true);

    // Initialize the b-spline moduli.

    int maxSize = max(max(gridSizeX, gridSizeY), gridSizeZ);
    vector<double> data(PmeOrder);
    vector<double> ddata(PmeOrder);
    vector<double> bsplines_data(maxSize);
    data[PmeOrder-1] = 0.0;
    data[1] = 0.0;
    data[0] = 1.0;
    for (int i = 3; i < PmeOrder; i++) {
        double div = 1.0/(i-1.0);
        data[i-1] = 0.0;
        for (int j = 1; j < (i-1); j++)
            data[i-j-1] = div*(j*data[i-j-2]+(i-j)*data[i-j-1]);
        data[0] = div*data[0];
    }

    // Differentiate.

    ddata[0] = -data[0];
    for (int i = 1; i < PmeOrder; i++)
        ddata[i] = data[i-1]-data[i];
    double div = 1.0/(PmeOrder-1);
    data[PmeOrder-1] = 0.0;
    for (int i = 1; i < (PmeOrder-1); i++)
        data[PmeOrder-i-1] = div*(i*data[PmeOrder-i-2]+(PmeOrder-i)*data[PmeOrder-i-1]);
    data[0] = div*data[0];
    for (int i = 0; i < maxSize; i++)
        bsplines_data[i] = 0.0;
    for (int i = 1; i <= PmeOrder; i++)
        bsplines_data[i] = data[i-1];

    // Evaluate the actual bspline moduli for X/Y/Z.

    for (int dim = 0; dim < 3; dim++) {
        int ndata = (dim == 0 ? gridSizeX : dim == 1 ? gridSizeY : gridSizeZ);
        vector<double> moduli(ndata);
        for (int i = 0; i < ndata; i++) {
            double sc = 0.0;
            double ss = 0.0;
            for (int j = 0; j < ndata; j++) {
                double arg = (2.0*M_PI*i*j)/ndata;
                sc += bsplines_data[j]*cos(arg);
                ss += bsplines_data[j]*sin(arg);
            }
            moduli[i] = sc*sc+ss*ss;
        }
        for (int i = 0; i < ndata; i++)
            if (moduli[i] < 1.0e-7)
                moduli[i] = (moduli[(i-1+ndata)%ndata]+moduli[(i+1)%ndata])*0.5;
        if (dim == 0)
            pmeBsplineModuliX.upload(moduli, true);
        else if (dim == 1)
            pmeBsplineModuliY.upload(moduli, true);
        else
            pmeBsplineModuliZ.upload(moduli, true);
    }
}

void CommonCalcConstantPotentialForceKernel::pmeCompileKernels() {
    map<string, string> replacements;
    replacements["CHARGE"] = (usePosqCharges ? "pos.w" : "charges[atom]");
    ComputeProgram program = cc.compileProgram(cc.replaceStrings(CommonKernelSources::pme, replacements), pmeDefines);
    
    pmeGridIndexKernel = program->createKernel("findAtomGridIndex");
    pmeGridIndexKernel->addArg(cc.getPosq());
    pmeGridIndexKernel->addArg(pmeAtomGridIndex);
    for (int i = 0; i < 8; i++) {
        pmeGridIndexKernel->addArg();
    }

    pmeSpreadChargeKernel = program->createKernel("gridSpreadCharge");
    pmeSpreadChargeKernel->addArg(cc.getPosq());
    if (useFixedPointChargeSpreading) {
        pmeSpreadChargeKernel->addArg(pmeGrid2);
    }
    else {
        pmeSpreadChargeKernel->addArg(pmeGrid1);
    }
    for (int i = 0; i < 8; i++) {
        pmeSpreadChargeKernel->addArg();
    }
    pmeSpreadChargeKernel->addArg(pmeAtomGridIndex);
    pmeSpreadChargeKernel->addArg(charges);

    pmeConvolutionKernel = program->createKernel("reciprocalConvolution");
    pmeConvolutionKernel->addArg(pmeGrid2);
    pmeConvolutionKernel->addArg(pmeBsplineModuliX);
    pmeConvolutionKernel->addArg(pmeBsplineModuliY);
    pmeConvolutionKernel->addArg(pmeBsplineModuliZ);
    for (int i = 0; i < 3; i++) {
        pmeConvolutionKernel->addArg();
    }

    pmeEvalEnergyKernel = program->createKernel("gridEvaluateEnergy");
    pmeEvalEnergyKernel->addArg(pmeGrid2);
    pmeEvalEnergyKernel->addArg(cc.getEnergyBuffer());
    pmeEvalEnergyKernel->addArg(pmeBsplineModuliX);
    pmeEvalEnergyKernel->addArg(pmeBsplineModuliY);
    pmeEvalEnergyKernel->addArg(pmeBsplineModuliZ);
    for (int i = 0; i < 3; i++) {
        pmeEvalEnergyKernel->addArg();
    }

    pmeInterpolateForceKernel = program->createKernel("gridInterpolateForce");
    pmeInterpolateForceKernel->addArg(cc.getPosq());
    pmeInterpolateForceKernel->addArg(cc.getLongForceBuffer());
    pmeInterpolateForceKernel->addArg(pmeGrid1);
    for (int i = 0; i < 8; i++) {
        pmeInterpolateForceKernel->addArg();
    }
    pmeInterpolateForceKernel->addArg(pmeAtomGridIndex);
    pmeInterpolateForceKernel->addArg(charges);

    if (hasElectrodes) {
        pmeInterpolateChargeDerivativesKernel = program->createKernel("gridInterpolateChargeDerivatives");
        pmeInterpolateChargeDerivativesKernel->addArg(cc.getPosq());
        pmeInterpolateChargeDerivativesKernel->addArg(chargeDerivativesFixed);
        pmeInterpolateChargeDerivativesKernel->addArg(pmeGrid1);
        for (int i = 0; i < 8; i++) {
            pmeInterpolateChargeDerivativesKernel->addArg();
        }
        pmeInterpolateChargeDerivativesKernel->addArg(elecToSys);
    }

    if (useFixedPointChargeSpreading) {
        pmeFinishSpreadChargeKernel = program->createKernel("finishSpreadCharge");
        pmeFinishSpreadChargeKernel->addArg(pmeGrid2);
        pmeFinishSpreadChargeKernel->addArg(pmeGrid1);
    }
}

void CommonCalcConstantPotentialForceKernel::initPmeExecute() {
    if (useFixedPointChargeSpreading) {
        cc.clearBuffer(pmeGrid2);
    }
    else {
        cc.clearBuffer(pmeGrid1);
    }
}

void CommonCalcConstantPotentialForceKernel::pmeExecute(bool includeEnergy, bool includeForces, bool includeChargeDerivatives, bool init) {
    if (init) {
        initPmeExecute();
    }
    if (pmeShouldSort) {
        pmeGridIndexKernel->execute(cc.getNumAtoms());
        sort->sort(pmeAtomGridIndex);
        pmeShouldSort = false;
    }
    pmeSpreadChargeKernel->execute(cc.getNumAtoms());
    if (useFixedPointChargeSpreading) {
        pmeFinishSpreadChargeKernel->execute(gridSizeX*gridSizeY*gridSizeZ);
    }
    fft->execFFT(pmeGrid1, pmeGrid2, true);

    if (includeEnergy) {
        pmeEvalEnergyKernel->execute(gridSizeX*gridSizeY*gridSizeZ);
    }

    if (includeForces || includeChargeDerivatives) {
        pmeConvolutionKernel->execute(gridSizeX*gridSizeY*gridSizeZ);
        fft->execFFT(pmeGrid2, pmeGrid1, false);
    
        if (includeForces) {
            if (deviceIsCpu) {
                pmeInterpolateForceKernel->execute(cc.getNumThreadBlocks(), 1);
            }
            else {
                pmeInterpolateForceKernel->execute(cc.getNumAtoms());
            }
        }
    
        if (includeChargeDerivatives && hasElectrodes) {
            if (deviceIsCpu) {
                pmeInterpolateChargeDerivativesKernel->execute(cc.getNumThreadBlocks(), 1);
            }
            else {
                pmeInterpolateChargeDerivativesKernel->execute(numElectrodeParticles);
            }
        }
    }
}

void CommonCalcConstantPotentialForceKernel::setKernelInputs(bool includeEnergy, bool includeForces) {
    double determinant = boxVectors[0][0] * boxVectors[1][1] * boxVectors[2][2];
    double scale = 1.0 / determinant;
    mm_double4 recipBoxVectors[3];
    recipBoxVectors[0] = mm_double4(boxVectors[1][1] * boxVectors[2][2] * scale, 0, 0, 0);
    recipBoxVectors[1] = mm_double4(-boxVectors[1][0] * boxVectors[2][2] * scale, boxVectors[0][0] * boxVectors[2][2] * scale, 0, 0);
    recipBoxVectors[2] = mm_double4((boxVectors[1][0] * boxVectors[2][1] - boxVectors[1][1] * boxVectors[2][0]) * scale, -boxVectors[0][0] * boxVectors[2][1] * scale, boxVectors[0][0] * boxVectors[1][1] * scale, 0);
    mm_float4 recipBoxVectorsFloat[3];
    for (int i = 0; i < 3; i++) {
        recipBoxVectorsFloat[i] = mm_float4((float) recipBoxVectors[i].x, (float) recipBoxVectors[i].y, (float) recipBoxVectors[i].z, 0);
    }

    setPeriodicBoxArgs(cc, pmeGridIndexKernel, 2);
    setPeriodicBoxArgs(cc, pmeSpreadChargeKernel, 2);
    if (includeEnergy || includeForces) {
        if (cc.getUseDoublePrecision()) {
            evaluateSelfEnergyForcesKernel->setArg(6, mm_double4(boxVectors[0][0], boxVectors[0][1], boxVectors[0][2], 0));
            evaluateSelfEnergyForcesKernel->setArg(7, mm_double4(boxVectors[1][0], boxVectors[1][1], boxVectors[1][2], 0));
            evaluateSelfEnergyForcesKernel->setArg(8, mm_double4(boxVectors[2][0], boxVectors[2][1], boxVectors[2][2], 0));
        }
        else {
            evaluateSelfEnergyForcesKernel->setArg(6, mm_float4((float) boxVectors[0][0], (float) boxVectors[0][1], (float) boxVectors[0][2], 0));
            evaluateSelfEnergyForcesKernel->setArg(7, mm_float4((float) boxVectors[1][0], (float) boxVectors[1][1], (float) boxVectors[1][2], 0));
            evaluateSelfEnergyForcesKernel->setArg(8, mm_float4((float) boxVectors[2][0], (float) boxVectors[2][1], (float) boxVectors[2][2], 0));
        }
    }
    if (includeForces) {
        setPeriodicBoxArgs(cc, pmeInterpolateForceKernel, 3);
    }
    if (hasElectrodes) {
        setPeriodicBoxArgs(cc, pmeInterpolateChargeDerivativesKernel, 3);
        setPeriodicBoxArgs(cc, evaluateDirectDerivativesKernel, 6);

        evaluateDirectDerivativesKernel->setArg(17, (unsigned int) cc.getNonbondedUtilities().getInteractingTiles().getSize());

        if (cc.getUseDoublePrecision()) {
            finishDerivativesKernel->setArg(7, mm_double4(boxVectors[0][0], boxVectors[0][1], boxVectors[0][2], 0));
            finishDerivativesKernel->setArg(8, mm_double4(boxVectors[1][0], boxVectors[1][1], boxVectors[1][2], 0));
            finishDerivativesKernel->setArg(9, mm_double4(boxVectors[2][0], boxVectors[2][1], boxVectors[2][2], 0));

            finishDerivativesKernel->setArg(10, mm_double4(externalField[0], externalField[1], externalField[2], 0));
        }
        else {
            finishDerivativesKernel->setArg(7, mm_float4((float) boxVectors[0][0], (float) boxVectors[0][1], (float) boxVectors[0][2], 0));
            finishDerivativesKernel->setArg(8, mm_float4((float) boxVectors[1][0], (float) boxVectors[1][1], (float) boxVectors[1][2], 0));
            finishDerivativesKernel->setArg(9, mm_float4((float) boxVectors[2][0], (float) boxVectors[2][1], (float) boxVectors[2][2], 0));

            finishDerivativesKernel->setArg(10, mm_float4((float) externalField[0], (float) externalField[1], (float) externalField[2], 0));
        }
    }

    if (cc.getUseDoublePrecision()) {
        pmeGridIndexKernel->setArg(7, recipBoxVectors[0]);
        pmeGridIndexKernel->setArg(8, recipBoxVectors[1]);
        pmeGridIndexKernel->setArg(9, recipBoxVectors[2]);

        pmeSpreadChargeKernel->setArg(7, recipBoxVectors[0]);
        pmeSpreadChargeKernel->setArg(8, recipBoxVectors[1]);
        pmeSpreadChargeKernel->setArg(9, recipBoxVectors[2]);

        if (includeEnergy) {
            pmeEvalEnergyKernel->setArg<mm_double4>(5, recipBoxVectors[0]);
            pmeEvalEnergyKernel->setArg<mm_double4>(6, recipBoxVectors[1]);
            pmeEvalEnergyKernel->setArg<mm_double4>(7, recipBoxVectors[2]);
        }

        if (includeForces || hasElectrodes) {
            pmeConvolutionKernel->setArg<mm_double4>(4, recipBoxVectors[0]);
            pmeConvolutionKernel->setArg<mm_double4>(5, recipBoxVectors[1]);
            pmeConvolutionKernel->setArg<mm_double4>(6, recipBoxVectors[2]);

            if (includeForces) {
                pmeInterpolateForceKernel->setArg(8, recipBoxVectors[0]);
                pmeInterpolateForceKernel->setArg(9, recipBoxVectors[1]);
                pmeInterpolateForceKernel->setArg(10, recipBoxVectors[2]);
            }

            if (hasElectrodes) {
                pmeInterpolateChargeDerivativesKernel->setArg(8, recipBoxVectors[0]);
                pmeInterpolateChargeDerivativesKernel->setArg(9, recipBoxVectors[1]);
                pmeInterpolateChargeDerivativesKernel->setArg(10, recipBoxVectors[2]);
            }
        }
    }
    else {
        pmeGridIndexKernel->setArg(7, recipBoxVectorsFloat[0]);
        pmeGridIndexKernel->setArg(8, recipBoxVectorsFloat[1]);
        pmeGridIndexKernel->setArg(9, recipBoxVectorsFloat[2]);

        pmeSpreadChargeKernel->setArg(7, recipBoxVectorsFloat[0]);
        pmeSpreadChargeKernel->setArg(8, recipBoxVectorsFloat[1]);
        pmeSpreadChargeKernel->setArg(9, recipBoxVectorsFloat[2]);

        if (includeEnergy) {
            pmeEvalEnergyKernel->setArg<mm_float4>(5, recipBoxVectorsFloat[0]);
            pmeEvalEnergyKernel->setArg<mm_float4>(6, recipBoxVectorsFloat[1]);
            pmeEvalEnergyKernel->setArg<mm_float4>(7, recipBoxVectorsFloat[2]);
        }

        if (includeForces || hasElectrodes) {
            pmeConvolutionKernel->setArg<mm_float4>(4, recipBoxVectorsFloat[0]);
            pmeConvolutionKernel->setArg<mm_float4>(5, recipBoxVectorsFloat[1]);
            pmeConvolutionKernel->setArg<mm_float4>(6, recipBoxVectorsFloat[2]);

            if (includeForces) {
                pmeInterpolateForceKernel->setArg(8, recipBoxVectorsFloat[0]);
                pmeInterpolateForceKernel->setArg(9, recipBoxVectorsFloat[1]);
                pmeInterpolateForceKernel->setArg(10, recipBoxVectorsFloat[2]);
            }

            if (hasElectrodes) {
                pmeInterpolateChargeDerivativesKernel->setArg(8, recipBoxVectorsFloat[0]);
                pmeInterpolateChargeDerivativesKernel->setArg(9, recipBoxVectorsFloat[1]);
                pmeInterpolateChargeDerivativesKernel->setArg(10, recipBoxVectorsFloat[2]);
            }
        }
    }

    const vector<mm_int4>& newPosCellOffsets = cc.getPosCellOffsets();
    bool mustUpdatePosCellOffsets = false;
    for (int i = 0; i < numParticles; i++) {
        mm_int4 oldOffset = hostPosCellOffsets[i];
        mm_int4 newOffset = newPosCellOffsets[i];
        if (newOffset.x != oldOffset.x || newOffset.y != oldOffset.y || newOffset.z != oldOffset.z) {
            mustUpdatePosCellOffsets = true;
            break;
        }
    }
    if (mustUpdatePosCellOffsets) {
        hostPosCellOffsets.assign(newPosCellOffsets.begin(), newPosCellOffsets.end());
        posCellOffsets.upload(hostPosCellOffsets);
    }
}

void CommonCalcConstantPotentialForceKernel::ensureValidNeighborList() {
    // Save the forcesValid flag since we use it to monitor the neighbor list build.
    bool oldForcesValid = cc.getForcesValid();

    do {
        // If we need to try to build the neighbor list again (i.e., it needs to be made bigger),
        // getForcesValid() will return false after computeInteractions() completes.
        cc.setForcesValid(true);
        cc.getNonbondedUtilities().prepareInteractions(1 << forceGroup);
        cc.getNonbondedUtilities().computeInteractions(1 << forceGroup, false, false);
    } while(!cc.getForcesValid());

    if (hasElectrodes) {
        evaluateDirectDerivativesKernel->setArg(17, (unsigned int) cc.getNonbondedUtilities().getInteractingTiles().getSize());
    }

    // Restore the old value of the flag.
    cc.setForcesValid(oldForcesValid);
}
