/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2026 Stanford University and the Authors.           *
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

#include "openmm/common/CommonMinimizeKernel.h"
#include "openmm/common/ContextSelector.h"
#include "CommonKernelSources.h"
#include <map>

using namespace OpenMM;
using namespace std;

// Constants for constraint handling.

const double CommonMinimizeKernel::minConstraintTol = 1e-4;
const double CommonMinimizeKernel::kRestraintScale = 100;
const double CommonMinimizeKernel::prevMaxErrorInit = 1e10;
const double CommonMinimizeKernel::kRestraintScaleUp = 10;
const double CommonMinimizeKernel::constraintTolScale = 100;

// Constants for L-BFGS.

const double CommonMinimizeKernel::fTol = 1e-4;
const double CommonMinimizeKernel::wolfeParam = 0.9;
const double CommonMinimizeKernel::stepScaleDown = 0.5;
const double CommonMinimizeKernel::stepScaleUp = 2.1;
const double CommonMinimizeKernel::minStep = 1e-20;
const double CommonMinimizeKernel::maxStep = 1e20;
const int CommonMinimizeKernel::numVectors = 6;
const int CommonMinimizeKernel::maxLineSearchIterations = 40;

CommonMinimizeKernel::~CommonMinimizeKernel() {
    if (cpuContext != NULL) {
        delete cpuContext;
    }
}

void CommonMinimizeKernel::initialize(const System& system) {
    numParticles = system.getNumParticles();
    numVariables = numParticles * 3;
    numConstraints = system.getNumConstraints();

    hostPositions.resize(numParticles);
    hostX.resize(numVariables);
    hostGrad.resize(numVariables);
    hostConstraintIndices.resize(numConstraints);
    hostConstraintDistances.resize(numConstraints);

    for (int i = 0; i < numConstraints; i++) {
        system.getConstraintParameters(i, hostConstraintIndices[i].x, hostConstraintIndices[i].y, hostConstraintDistances[i]);
    }
}

void CommonMinimizeKernel::execute(ContextImpl& context, double tolerance, int maxIterations, MinimizationReporter* reporter) {
    ContextSelector selector(cc);

    if (!isSetup) {
        // Load system and integrator information.

        mixedIsDouble = cc.getUseMixedPrecision() || cc.getUseDoublePrecision();
        elementSize = mixedIsDouble ? sizeof(double) : sizeof(float);
        threadBlockSize = cc.getMaxThreadBlockSize();
        numVariableBlocks = (numVariables + threadBlockSize - 1) / threadBlockSize;
        numConstraintBlocks = (numConstraints + threadBlockSize - 1) / threadBlockSize;

        // Initialize all device-side arrays and compile kernels.

        setup(context);
        isSetup = true;
    }

    const Integrator& integrator = context.getIntegrator();
    forceGroups = integrator.getIntegrationForceGroups();
    constraintTol = integrator.getConstraintTolerance();

    this->tolerance = tolerance * sqrt((double) numParticles);
    this->maxIterations = maxIterations;
    this->reporter = reporter;

    double workingConstraintTol = max(minConstraintTol, constraintTol);
    kRestraint = kRestraintScale / workingConstraintTol;
    context.applyConstraints(workingConstraintTol);

    recordInitialPosKernel->execute(numParticles);

    double prevMaxError = prevMaxErrorInit;
    while (true) {
        lbfgs(context);

        if (!numConstraints) {
            // There are no constraints, so we are finished.
            break;
        }

        getConstraintErrorKernel->execute(threadBlockSize, threadBlockSize);
        double maxError = downloadReturnValue();
        if (maxError <= workingConstraintTol) {
            // All constraints are satisfied.
            break;
        }
        restorePosKernel->setArg(2, xInit);
        restorePosKernel->execute(numParticles);
        if (maxError >= prevMaxError) {
            // Further tightening the springs doesn't seem to be helping, so just give up.
            break;
        }
        prevMaxError = maxError;
        kRestraint *= kRestraintScaleUp;
        if (maxError > constraintTolScale * workingConstraintTol) {
            // We've gotten far enough from a valid state that we might have
            // trouble getting back, so reset to the original positions.
            xInit.copyTo(x);
        }
    }

    if (constraintTol < workingConstraintTol) {
        context.applyConstraints(constraintTol);
    }
}

void CommonMinimizeKernel::setup(ContextImpl& context) {
    // Initialize arrays.

    if (numConstraints) {
        constraintIndices.initialize<mm_int2>(cc, numConstraints, "constraintIndices");
        constraintDistances.initialize(cc, numConstraints, elementSize, "constraintDistances");
    }
    xInit.initialize(cc, numVariables, elementSize, "xInit");
    x.initialize(cc, numVariables, elementSize, "x");
    xPrev.initialize(cc, numVariables, elementSize, "xPrev");
    grad.initialize(cc, numVariables, elementSize, "grad");
    gradPrev.initialize(cc, numVariables, elementSize, "gradPrev");
    dir.initialize(cc, numVariables, elementSize, "dir");
    alpha.initialize(cc, numVectors, elementSize, "alpha");
    scale.initialize(cc, numVectors, elementSize, "scale");
    xDiff.initialize(cc, numVectors * numVariables, elementSize, "xDiff");
    gradDiff.initialize(cc, numVectors * numVariables, elementSize, "gradDiff");
    reduceBuffer.initialize(cc, max(2 * numVariableBlocks, numConstraintBlocks), elementSize, "reduceBuffer");
    returnFlag.initialize<int>(cc, 1, "returnFlag");
    returnValue.initialize(cc, 1, elementSize, "returnValue");

    // Compile kernels and set arguments.

    map<string, string> defines;
    defines["THREAD_BLOCK_SIZE"] = cc.intToString(threadBlockSize);
    defines["NUM_CONSTRAINT_BLOCKS"] = cc.intToString(numConstraintBlocks);
    defines["NUM_CONSTRAINTS"] = cc.intToString(numConstraints);
    defines["NUM_PADDED"] = cc.intToString(cc.getPaddedNumAtoms());
    defines["NUM_PARTICLES"] = cc.intToString(numParticles);
    defines["NUM_VARIABLE_BLOCKS"] = cc.intToString(numVariableBlocks);
    defines["NUM_VARIABLES"] = cc.intToString(numVariables);
    defines["NUM_VECTORS"] = cc.intToString(numVectors);

    ComputeProgram program = cc.compileProgram(CommonKernelSources::minimize, defines);

    recordInitialPosKernel = program->createKernel("recordInitialPos");
    recordInitialPosKernel->addArg(cc.getPosq());
    recordInitialPosKernel->addArg(cc.getAtomIndexArray());
    recordInitialPosKernel->addArg(xInit);
    recordInitialPosKernel->addArg(x);
    if (cc.getUseMixedPrecision()) {
        recordInitialPosKernel->addArg(cc.getPosqCorrection());
    }

    restorePosKernel = program->createKernel("restorePos");
    restorePosKernel->addArg(cc.getPosq());
    restorePosKernel->addArg(cc.getAtomIndexArray());
    restorePosKernel->addArg(); // x (could also be launched with xInit)
    if (cc.getUseMixedPrecision()) {
        restorePosKernel->addArg(cc.getPosqCorrection());
    }

    convertForcesKernel = program->createKernel("convertForces");
    convertForcesKernel->addArg(cc.getVelm());
    convertForcesKernel->addArg(cc.getLongForceBuffer());
    convertForcesKernel->addArg(cc.getAtomIndexArray());
    convertForcesKernel->addArg(grad);
    convertForcesKernel->addArg(returnFlag);

    if (numConstraints) {
        getConstraintEnergyForcesKernel = program->createKernel("getConstraintEnergyForces");
        getConstraintEnergyForcesKernel->addArg(cc.getLongForceBuffer());
        getConstraintEnergyForcesKernel->addArg(cc.getAtomIndexArray());
        getConstraintEnergyForcesKernel->addArg(constraintIndices);
        getConstraintEnergyForcesKernel->addArg(constraintDistances);
        getConstraintEnergyForcesKernel->addArg(x);
        getConstraintEnergyForcesKernel->addArg(reduceBuffer);
        getConstraintEnergyForcesKernel->addArg(); // kRestraint

        reduceConstraintEnergyKernel = program->createKernel("reduceConstraintEnergy");
        reduceConstraintEnergyKernel->addArg(reduceBuffer);
        reduceConstraintEnergyKernel->addArg(returnValue);

        getConstraintErrorKernel = program->createKernel("getConstraintError");
        getConstraintErrorKernel->addArg(constraintIndices);
        getConstraintErrorKernel->addArg(constraintDistances);
        getConstraintErrorKernel->addArg(x);
        getConstraintErrorKernel->addArg(returnValue);
    }

    initializeDirKernel = program->createKernel("initializeDir");
    initializeDirKernel->addArg(grad);
    initializeDirKernel->addArg(dir);
    initializeDirKernel->addArg(returnValue);

    gradNormPart1Kernel = program->createKernel("gradNormPart1");
    gradNormPart1Kernel->addArg(grad);
    gradNormPart1Kernel->addArg(reduceBuffer);

    gradNormPart2Kernel = program->createKernel("gradNormPart2");
    gradNormPart2Kernel->addArg(grad);
    gradNormPart2Kernel->addArg(reduceBuffer);
    gradNormPart2Kernel->addArg(returnValue);

    getDiffKernel = program->createKernel("getDiff");
    getDiffKernel->addArg(x);
    getDiffKernel->addArg(xPrev);
    getDiffKernel->addArg(grad);
    getDiffKernel->addArg(gradPrev);
    getDiffKernel->addArg(xDiff);
    getDiffKernel->addArg(gradDiff);
    getDiffKernel->addArg(); // end

    getScalePart1Kernel = program->createKernel("getScalePart1");
    getScalePart1Kernel->addArg(xDiff);
    getScalePart1Kernel->addArg(gradDiff);
    getScalePart1Kernel->addArg(reduceBuffer);
    getScalePart1Kernel->addArg(); // end

    getScalePart2Kernel = program->createKernel("getScalePart2");
    getScalePart2Kernel->addArg(xDiff);
    getScalePart2Kernel->addArg(gradDiff);
    getScalePart2Kernel->addArg(scale);
    getScalePart2Kernel->addArg(reduceBuffer);
    getScalePart2Kernel->addArg(returnValue);
    getScalePart2Kernel->addArg(); // end

    reinitializeDirKernel = program->createKernel("reinitializeDir");
    reinitializeDirKernel->addArg(grad);
    reinitializeDirKernel->addArg(dir);
    reinitializeDirKernel->addArg(xDiff);
    reinitializeDirKernel->addArg(reduceBuffer);
    reinitializeDirKernel->addArg(); // vectorIndex

    reduceAlphaKernel = program->createKernel("reduceAlpha");
    reduceAlphaKernel->addArg(alpha);
    reduceAlphaKernel->addArg(scale);
    reduceAlphaKernel->addArg(reduceBuffer);
    reduceAlphaKernel->addArg(); // vectorIndex

    updateDirAlphaKernel = program->createKernel("updateDirAlpha");
    updateDirAlphaKernel->addArg(dir);
    updateDirAlphaKernel->addArg(alpha);
    updateDirAlphaKernel->addArg(xDiff);
    updateDirAlphaKernel->addArg(gradDiff);
    updateDirAlphaKernel->addArg(reduceBuffer);
    updateDirAlphaKernel->addArg(); // vectorIndex

    scaleDirKernel = program->createKernel("scaleDir");
    scaleDirKernel->addArg(dir);
    scaleDirKernel->addArg(alpha);
    scaleDirKernel->addArg(gradDiff);
    scaleDirKernel->addArg(reduceBuffer);
    scaleDirKernel->addArg(returnValue);
    scaleDirKernel->addArg(); // vectorIndex

    reduceBetaKernel = program->createKernel("reduceBeta");
    reduceBetaKernel->addArg(scale);
    reduceBetaKernel->addArg(reduceBuffer);
    reduceBetaKernel->addArg(returnValue);
    reduceBetaKernel->addArg(); // vectorIndex

    updateDirBetaKernel = program->createKernel("updateDirBeta");
    updateDirBetaKernel->addArg(dir);
    updateDirBetaKernel->addArg(alpha);
    updateDirBetaKernel->addArg(xDiff);
    updateDirBetaKernel->addArg(gradDiff);
    updateDirBetaKernel->addArg(reduceBuffer);
    updateDirBetaKernel->addArg(returnValue);
    updateDirBetaKernel->addArg(); // vectorIndex

    updateDirFinalKernel = program->createKernel("updateDirFinal");
    updateDirFinalKernel->addArg(dir);
    updateDirFinalKernel->addArg(alpha);
    updateDirFinalKernel->addArg(xDiff);
    updateDirFinalKernel->addArg(returnValue);
    updateDirFinalKernel->addArg(); // vectorIndex

    lineSearchStepKernel = program->createKernel("lineSearchStep");
    lineSearchStepKernel->addArg(x);
    lineSearchStepKernel->addArg(xPrev);
    lineSearchStepKernel->addArg(dir);
    lineSearchStepKernel->addArg(); // step

    lineSearchDotPart1Kernel = program->createKernel("lineSearchDotPart1");
    lineSearchDotPart1Kernel->addArg(grad);
    lineSearchDotPart1Kernel->addArg(dir);
    lineSearchDotPart1Kernel->addArg(reduceBuffer);

    lineSearchDotPart2Kernel = program->createKernel("lineSearchDotPart2");
    lineSearchDotPart2Kernel->addArg(reduceBuffer);
    lineSearchDotPart2Kernel->addArg(returnValue);

    // Upload constraint data.

    if (numConstraints) {
        constraintIndices.upload(hostConstraintIndices);
        constraintDistances.upload(hostConstraintDistances, true);
    }
}

void CommonMinimizeKernel::lbfgs(ContextImpl& context) {
    // Evaluate the gradient at the starting point and get a starting step size.

    bool overflow;
    energy = evaluate(context, overflow);
    if (overflow) {
        throw OpenMMException("Energy or force at minimization starting point is infinite or NaN.");
    }
    gradNormPart1Kernel->execute(numVariables, threadBlockSize);
    gradNormPart2Kernel->execute(threadBlockSize, threadBlockSize);
    if (downloadReturnValue() <= tolerance) {
        return;
    }
    initializeDirKernel->execute(numVariables);
    double step = 1.0;

    for (int iteration = 1, end = 0;;) {
        // Try a line search (if it fails, revert to the position and gradient
        // at the start of the search and abort the optimization).

        x.copyTo(xPrev);
        grad.copyTo(gradPrev);
        if (!lineSearch(context, step)) {
            xPrev.copyTo(x);
            gradPrev.copyTo(grad);
            return;
        }

        // Check for convergence or exceeding the maximum number of steps.

        if (reporter != NULL && report(context, iteration)) {
            return;
        }
        gradNormPart1Kernel->execute(numVariables, threadBlockSize);
        gradNormPart2Kernel->execute(threadBlockSize, threadBlockSize);
        if (downloadReturnValue() <= tolerance || (maxIterations && maxIterations < iteration + 1)) {
            return;
        }

        // Do L-BFGS update of search direction.

        getDiffKernel->setArg(6, end);
        getDiffKernel->execute(numVariables);
        getScalePart1Kernel->setArg(3, end);
        getScalePart1Kernel->execute(numVariables, threadBlockSize);
        getScalePart2Kernel->setArg(5, end);
        getScalePart2Kernel->execute(threadBlockSize, threadBlockSize);

        int limit = min(numVectors, iteration++);
        if (++end >= numVectors) {
            end -= numVectors;
        }

        int vectorIndex = (end ? end : numVectors) - 1;
        reinitializeDirKernel->setArg(4, vectorIndex);
        reinitializeDirKernel->execute(numVariables, threadBlockSize);

        for (int vector = 0; vector < limit; vector++) {
            if (vector && --vectorIndex < 0) {
                vectorIndex += numVectors;
            }

            reduceAlphaKernel->setArg(3, vectorIndex);
            reduceAlphaKernel->execute(threadBlockSize, threadBlockSize);

            if (vector < limit - 1) {
                updateDirAlphaKernel->setArg(5, vectorIndex);
                updateDirAlphaKernel->execute(numVariables, threadBlockSize);
            }
        }

        scaleDirKernel->setArg(5, vectorIndex);
        scaleDirKernel->execute(numVariables);

        for (int vector = 0; vector < limit; vector++) {
            reduceBetaKernel->setArg(3, vectorIndex);
            reduceBetaKernel->execute(threadBlockSize, threadBlockSize);

            if (vector < limit - 1) {
                updateDirBetaKernel->setArg(6, vectorIndex);
                updateDirBetaKernel->execute(numVariables, threadBlockSize);

                if (++vectorIndex >= numVectors) {
                    vectorIndex -= numVectors;
                }
            }
        }

        updateDirFinalKernel->setArg(4, vectorIndex);
        updateDirFinalKernel->execute(numVariables);

        step = 1.0;
    }
}

bool CommonMinimizeKernel::lineSearch(ContextImpl& context, double& step) {
    // Check state at starting point for line search.

    lineSearchDotPart1Kernel->execute(numVariables, threadBlockSize);
    lineSearchDotPart2Kernel->execute(threadBlockSize, threadBlockSize);
    double dotStart = downloadReturnValue();
    if (dotStart > 0) {
        return false;
    }
    double energyStart = energy;
    double stepScale;

    // Take line search steps.

    for (int count = 0;;) {
        if (mixedIsDouble) {
            lineSearchStepKernel->setArg(3, step);
        }
        else {
            lineSearchStepKernel->setArg(3, (float) step);
        }
        lineSearchStepKernel->execute(numVariables);
        bool overflow;
        energy = evaluate(context, overflow);
        count++;

        if (overflow || energy > energyStart + step * fTol * dotStart) {
            stepScale = stepScaleDown;
        }
        else {
            lineSearchDotPart1Kernel->execute(numVariables, threadBlockSize);
            lineSearchDotPart2Kernel->execute(threadBlockSize, threadBlockSize);
            double dot = downloadReturnValue();
            if (dot < wolfeParam * dotStart) {
                stepScale = stepScaleUp;
            }
            else if (dot > -wolfeParam * dotStart) {
                stepScale = stepScaleDown;
            }
            else {
                // Strong Wolfe condition satisfied.

                return true;
            }
        }

        if (step < minStep || step > maxStep || count >= maxLineSearchIterations) {
            return false;
        }

        step *= stepScale;
    }
}

double CommonMinimizeKernel::evaluate(ContextImpl& context, bool& overflow) {
    overflow = false;
    restorePosKernel->setArg(2, x);
    restorePosKernel->execute(numParticles);
    context.computeVirtualSites();
    double hostEnergy = context.calcForcesAndEnergy(true, true, forceGroups);

    if (numConstraints) {
        if (mixedIsDouble) {
            getConstraintEnergyForcesKernel->setArg(6, kRestraint);
        }
        else {
            getConstraintEnergyForcesKernel->setArg(6, (float) kRestraint);
        }
        getConstraintEnergyForcesKernel->execute(numConstraints, threadBlockSize);
        reduceConstraintEnergyKernel->execute(threadBlockSize, threadBlockSize);
        hostEnergy += downloadReturnValue();
    }

    // Copy forces into the gradient buffer; if they are too large, we will
    // fall back to downloading positions and computing forces on the CPU.

    cc.clearBuffer(returnFlag);
    convertForcesKernel->execute(numParticles);
    if (downloadReturnFlag()) {
        hostEnergy = evaluateCpu(context, overflow);
    }
    if (!isfinite(hostEnergy)) {
        overflow = true;
    }
    return hostEnergy;
}

double CommonMinimizeKernel::evaluateCpu(ContextImpl& context, bool& overflow) {
    // Create a CPU context if one has not already been created.

    const System& system = context.getSystem();

    if (cpuContext == NULL) {
        Platform* cpuPlatform;
        try {
            cpuPlatform = &Platform::getPlatformByName("CPU");
        }
        catch (...) {
            cpuPlatform = &Platform::getPlatformByName("Reference");
        }
        cpuContext = new Context(system, cpuIntegrator, *cpuPlatform);
    }

    // Download positions and evaluate forces on the CPU.

    x.download(hostX, true);
    for (int i = 0; i < numParticles; i++) {
        hostPositions[i] = Vec3(hostX[3 * i], hostX[3 * i + 1], hostX[3 * i + 2]);
    }
    cpuContext->setState(context.getOwner().getState(State::Parameters));
    cpuContext->setPositions(hostPositions);
    cpuContext->computeVirtualSites();
    State state = cpuContext->getState(State::Energy | State::Forces, false, forceGroups);
    double hostEnergy = state.getPotentialEnergy();
    const vector<Vec3>& hostForces = state.getForces();

    // Prepare the gradient to send back to the optimizer.

    for (int i = 0; i < numParticles; i++) {
        if (system.getParticleMass(i) != 0) {
            hostGrad[3 * i] = -hostForces[i][0];
            hostGrad[3 * i + 1] = -hostForces[i][1];
            hostGrad[3 * i + 2] = -hostForces[i][2];
        }
    }

    // Apply harmonic forces for constraints.

    for (int i = 0; i < numConstraints; i++) {
        mm_int2 indices = hostConstraintIndices[i];
        double distance = hostConstraintDistances[i];
        Vec3 delta = hostPositions[indices.y] - hostPositions[indices.x];
        double r2 = delta.dot(delta);
        double r = sqrt(r2);
        delta *= 1 / r;
        double dr = r - distance;
        double kdr = kRestraint * dr;
        hostEnergy += 0.5 * kdr * dr;
        if (system.getParticleMass(indices.x) != 0) {
            hostGrad[3 * indices.y] -= kdr * delta[0];
            hostGrad[3 * indices.y + 1] -= kdr * delta[1];
            hostGrad[3 * indices.y + 2] -= kdr * delta[2];
        }
        if (system.getParticleMass(indices.y) != 0) {
            hostGrad[3 * indices.x] += kdr * delta[0];
            hostGrad[3 * indices.x + 1] += kdr * delta[1];
            hostGrad[3 * indices.x + 2] += kdr * delta[2];
        }
    }

    // Check for overflow of the forces.  We need to check for this here and
    // either back up the line search or abort, since the GPU platforms won't
    // check for NaN positions, and the optimizer could get stuck otherwise.

    if (mixedIsDouble) {
        for (int i = 0; i < numVariables; i++) {
            if (!isfinite(hostGrad[i])) {
                overflow = true;
                return hostEnergy;
            }
        }
    }
    else {
        for (int i = 0; i < numVariables; i++) {
            if (!isfinite((float) hostGrad[i])) {
                overflow = true;
                return hostEnergy;
            }
        }
    }

    // Upload forces and return the final energy.

    grad.upload(hostGrad, true);
    return hostEnergy;
}

bool CommonMinimizeKernel::report(ContextImpl& context, int iteration) {
    x.download(hostX, true);
    grad.download(hostGrad, true);

    double restraintEnergy = 0.0, maxError = 0.0;
    for (int i = 0; i < numConstraints; i++) {
        mm_int2 indices = hostConstraintIndices[i];
        double distance = hostConstraintDistances[i];
        Vec3 delta = Vec3(hostX[3 * indices.y] - hostX[3 * indices.x], hostX[3 * indices.y + 1] - hostX[3 * indices.x + 1], hostX[3 * indices.y + 2] - hostX[3 * indices.x + 2]);
        double r2 = delta.dot(delta);
        double r = sqrt(r2);
        double dr = r - distance;
        restraintEnergy += 0.5 * kRestraint * dr * dr;
        maxError = max(maxError, fabs(dr) / distance);
    }

    map<string, double> args;
    args["restraint energy"] = restraintEnergy;
    args["system energy"] = energy - restraintEnergy;
    args["restraint strength"] = kRestraint;
    args["max constraint error"] = maxError;
    return reporter->report(iteration - 1, hostX, hostGrad, args);
}

int CommonMinimizeKernel::downloadReturnFlag() {
    int hostReturnFlag;
    returnFlag.download(&hostReturnFlag);
    return hostReturnFlag;
}

double CommonMinimizeKernel::downloadReturnValue() {
    if (mixedIsDouble) {
        double hostReturnValue;
        returnValue.download(&hostReturnValue);
        return hostReturnValue;
    }
    else {
        float hostReturnValue;
        returnValue.download(&hostReturnValue);
        return (double) hostReturnValue;
    }
}
