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
        pinnedMemory = cc.getPinnedBuffer();

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

    if (mixedIsDouble) {
        getDiffKernel->setArg(11, this->tolerance);
    }
    else {
        getDiffKernel->setArg(11, (float) this->tolerance);
    }

    double workingConstraintTol = max(minConstraintTol, constraintTol);
    kRestraint = kRestraintScale / workingConstraintTol;
    context.applyConstraints(workingConstraintTol);

    recordInitialPosKernel->execute(numParticles);

    double prevMaxError1 = prevMaxErrorInit, prevMaxError2 = prevMaxErrorInit;
    while (true) {
        lbfgs(context);

        if (!numConstraints) {
            // There are no constraints, so we are finished.
            break;
        }

        getConstraintErrorKernel->execute(threadBlockSize, threadBlockSize);
        double maxError = downloadReturnValueSync();
        if (maxError <= workingConstraintTol) {
            // All constraints are satisfied.
            break;
        }
        restorePosKernel->setArg(2, xInit);
        restorePosKernel->execute(numParticles);
        if (maxError >= prevMaxError2) {
            // Further tightening the springs doesn't seem to be helping, so just give up.
            break;
        }
        prevMaxError2 = prevMaxError1;
        prevMaxError1 = maxError;
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
    alpha.initialize(cc, numVectors + 1, elementSize, "alpha");
    scale.initialize(cc, numVectors, elementSize, "scale");
    xDiff.initialize(cc, numVectors * numVariables, elementSize, "xDiff");
    gradDiff.initialize(cc, numVectors * numVariables, elementSize, "gradDiff");
    reduceBuffer.initialize(cc, max(2 * numVariableBlocks, numConstraintBlocks), elementSize, "reduceBuffer");
    returnFlag.initialize<int>(cc, 1, "returnFlag");
    returnValue.initialize(cc, 1, elementSize, "returnValue");
    gradNorm.initialize(cc, 1, elementSize, "gradNorm");
    lineSearchData.initialize(cc, 4, elementSize, "lineSearchData");
    lineSearchDataBackup.initialize(cc, 4, elementSize, "lineSearchDataBackup");

    // Compile kernels and set arguments.

    map<string, string> defines;
    defines["THREAD_BLOCK_SIZE"] = cc.intToString(threadBlockSize);
    defines["NUM_VECTORS"] = cc.intToString(numVectors);
    defines["LBFGS_FTOL"] = cc.doubleToString(fTol);
    defines["LBFGS_WOLFE"] = cc.doubleToString(wolfeParam);
    defines["LBFGS_SCALE_DOWN"] = cc.doubleToString(stepScaleDown);
    defines["LBFGS_SCALE_UP"] = cc.doubleToString(stepScaleUp);
    defines["LBFGS_MIN_STEP"] = cc.doubleToString(minStep);
    defines["LBFGS_MAX_STEP"] = cc.doubleToString(maxStep);

    ComputeProgram program = cc.compileProgram(CommonKernelSources::minimize, defines);

    recordInitialPosKernel = program->createKernel("recordInitialPos");
    recordInitialPosKernel->addArg(cc.getPosq());
    recordInitialPosKernel->addArg(cc.getAtomIndexArray());
    recordInitialPosKernel->addArg(xInit);
    recordInitialPosKernel->addArg(x);
    recordInitialPosKernel->addArg(numParticles);
    if (cc.getUseMixedPrecision()) {
        recordInitialPosKernel->addArg(cc.getPosqCorrection());
    }

    restorePosKernel = program->createKernel("restorePos");
    restorePosKernel->addArg(cc.getPosq());
    restorePosKernel->addArg(cc.getAtomIndexArray());
    restorePosKernel->addArg(); // x (could also be launched with xInit)
    restorePosKernel->addArg(returnValue);
    restorePosKernel->addArg(numParticles);
    if (cc.getUseMixedPrecision()) {
        restorePosKernel->addArg(cc.getPosqCorrection());
    }

    convertForcesKernel = program->createKernel("convertForces");
    convertForcesKernel->addArg(cc.getVelm());
    convertForcesKernel->addArg(cc.getLongForceBuffer());
    convertForcesKernel->addArg(cc.getAtomIndexArray());
    convertForcesKernel->addArg(grad);
    convertForcesKernel->addArg(returnValue);
    convertForcesKernel->addArg(numParticles);
    convertForcesKernel->addArg(cc.getPaddedNumAtoms());

    if (numConstraints) {
        getConstraintEnergyForcesKernel = program->createKernel("getConstraintEnergyForces");
        getConstraintEnergyForcesKernel->addArg(cc.getLongForceBuffer());
        getConstraintEnergyForcesKernel->addArg(cc.getAtomIndexArray());
        getConstraintEnergyForcesKernel->addArg(constraintIndices);
        getConstraintEnergyForcesKernel->addArg(constraintDistances);
        getConstraintEnergyForcesKernel->addArg(x);
        getConstraintEnergyForcesKernel->addArg(returnValue);
        getConstraintEnergyForcesKernel->addArg(cc.getPaddedNumAtoms());
        getConstraintEnergyForcesKernel->addArg(numConstraints);
        getConstraintEnergyForcesKernel->addArg(); // kRestraint

        getConstraintErrorKernel = program->createKernel("getConstraintError");
        getConstraintErrorKernel->addArg(constraintIndices);
        getConstraintErrorKernel->addArg(constraintDistances);
        getConstraintErrorKernel->addArg(x);
        getConstraintErrorKernel->addArg(returnValue);
        getConstraintErrorKernel->addArg(numConstraints);
    }

    initializeDirKernel = program->createKernel("initializeDir");
    initializeDirKernel->addArg(grad);
    initializeDirKernel->addArg(dir);
    initializeDirKernel->addArg(gradNorm);
    initializeDirKernel->addArg(lineSearchData);
    initializeDirKernel->addArg(numVariables);

    gradNormPart1Kernel = program->createKernel("gradNormPart1");
    gradNormPart1Kernel->addArg(grad);
    gradNormPart1Kernel->addArg(reduceBuffer);
    gradNormPart1Kernel->addArg(numVariables);

    gradNormPart2Kernel = program->createKernel("gradNormPart2");
    gradNormPart2Kernel->addArg(grad);
    gradNormPart2Kernel->addArg(reduceBuffer);
    gradNormPart2Kernel->addArg(gradNorm);
    gradNormPart2Kernel->addArg(numVariables);
    gradNormPart2Kernel->addArg(numVariableBlocks);

    getDiffKernel = program->createKernel("getDiff");
    getDiffKernel->addArg(x);
    getDiffKernel->addArg(xPrev);
    getDiffKernel->addArg(grad);
    getDiffKernel->addArg(gradPrev);
    getDiffKernel->addArg(xDiff);
    getDiffKernel->addArg(gradDiff);
    getDiffKernel->addArg(reduceBuffer);
    getDiffKernel->addArg(returnFlag);
    getDiffKernel->addArg(gradNorm);
    getDiffKernel->addArg(numVariables);
    getDiffKernel->addArg(numVariableBlocks);
    getDiffKernel->addArg(); // tolerance
    getDiffKernel->addArg(); // end

    getScaleKernel = program->createKernel("getScale");
    getScaleKernel->addArg(alpha);
    getScaleKernel->addArg(scale);
    getScaleKernel->addArg(xDiff);
    getScaleKernel->addArg(gradDiff);
    getScaleKernel->addArg(reduceBuffer);
    getScaleKernel->addArg(returnFlag);
    getScaleKernel->addArg(returnValue);
    getScaleKernel->addArg(numVariables);
    getScaleKernel->addArg(numVariableBlocks);
    getScaleKernel->addArg(); // end

    reinitializeDirKernel = program->createKernel("reinitializeDir");
    reinitializeDirKernel->addArg(grad);
    reinitializeDirKernel->addArg(dir);
    reinitializeDirKernel->addArg(alpha);
    reinitializeDirKernel->addArg(scale);
    reinitializeDirKernel->addArg(xDiff);
    reinitializeDirKernel->addArg(returnFlag);
    reinitializeDirKernel->addArg(numVariables);
    reinitializeDirKernel->addArg(); // vectorIndex

    updateDirAlphaKernel = program->createKernel("updateDirAlpha");
    updateDirAlphaKernel->addArg(dir);
    updateDirAlphaKernel->addArg(alpha);
    updateDirAlphaKernel->addArg(scale);
    updateDirAlphaKernel->addArg(xDiff);
    updateDirAlphaKernel->addArg(gradDiff);
    updateDirAlphaKernel->addArg(returnFlag);
    updateDirAlphaKernel->addArg(numVariables);
    updateDirAlphaKernel->addArg(); // vectorIndex

    scaleDirKernel = program->createKernel("scaleDir");
    scaleDirKernel->addArg(dir);
    scaleDirKernel->addArg(alpha);
    scaleDirKernel->addArg(scale);
    scaleDirKernel->addArg(gradDiff);
    scaleDirKernel->addArg(returnFlag);
    scaleDirKernel->addArg(returnValue);
    scaleDirKernel->addArg(numVariables);
    scaleDirKernel->addArg(); // vectorIndex

    updateDirBetaKernel = program->createKernel("updateDirBeta");
    updateDirBetaKernel->addArg(dir);
    updateDirBetaKernel->addArg(alpha);
    updateDirBetaKernel->addArg(scale);
    updateDirBetaKernel->addArg(xDiff);
    updateDirBetaKernel->addArg(gradDiff);
    updateDirBetaKernel->addArg(returnFlag);
    updateDirBetaKernel->addArg(numVariables);
    updateDirBetaKernel->addArg(); // vectorIndex
    updateDirBetaKernel->addArg(); // vectorIndexAlpha

    updateDirFinalKernel = program->createKernel("updateDirFinal");
    updateDirFinalKernel->addArg(dir);
    updateDirFinalKernel->addArg(alpha);
    updateDirFinalKernel->addArg(xDiff);
    updateDirFinalKernel->addArg(returnFlag);
    updateDirFinalKernel->addArg(lineSearchData);
    updateDirFinalKernel->addArg(numVariables);
    updateDirFinalKernel->addArg(); // vectorIndex
    updateDirFinalKernel->addArg(); // vectorIndexAlpha

    lineSearchSetupKernel = program->createKernel("lineSearchSetup");
    lineSearchSetupKernel->addArg(x);
    lineSearchSetupKernel->addArg(xPrev);
    lineSearchSetupKernel->addArg(grad);
    lineSearchSetupKernel->addArg(gradPrev);
    lineSearchSetupKernel->addArg(dir);
    lineSearchSetupKernel->addArg(returnFlag);
    lineSearchSetupKernel->addArg(lineSearchData);
    lineSearchSetupKernel->addArg(numVariables);
    lineSearchSetupKernel->addArg(); // energyStart

    lineSearchStepKernel = program->createKernel("lineSearchStep");
    lineSearchStepKernel->addArg(x);
    lineSearchStepKernel->addArg(xPrev);
    lineSearchStepKernel->addArg(grad);
    lineSearchStepKernel->addArg(gradPrev);
    lineSearchStepKernel->addArg(dir);
    lineSearchStepKernel->addArg(reduceBuffer);
    lineSearchStepKernel->addArg(returnFlag);
    lineSearchStepKernel->addArg(lineSearchData);
    lineSearchStepKernel->addArg(lineSearchDataBackup);
    lineSearchStepKernel->addArg(numVariables);

    lineSearchDotKernel = program->createKernel("lineSearchDot");
    lineSearchDotKernel->addArg(grad);
    lineSearchDotKernel->addArg(dir);
    lineSearchDotKernel->addArg(lineSearchData);
    lineSearchDotKernel->addArg(returnFlag);
    lineSearchDotKernel->addArg(returnValue);
    lineSearchDotKernel->addArg(numVariables);
    lineSearchDotKernel->addArg(); // energy

    lineSearchContinueKernel = program->createKernel("lineSearchContinue");
    lineSearchContinueKernel->addArg(returnFlag);
    lineSearchContinueKernel->addArg(lineSearchData);

    downloadStartEvent = cc.createEvent();
    downloadFinishEvent = cc.createEvent();
    downloadQueue = cc.createQueue();

    // Upload constraint data.

    if (numConstraints) {
        constraintIndices.upload(hostConstraintIndices);
        constraintDistances.upload(hostConstraintDistances, true);
    }
}

void CommonMinimizeKernel::lbfgs(ContextImpl& context) {
    // Evaluate the energy and gradient at the starting point.

    evaluateGpu(context);
    energy += downloadReturnValueSync();
    if (!isfinite(energy)) {
        energy = evaluateCpu(context);
    }
    if (!isfinite(energy)) {
        throw OpenMMException("Energy or force at minimization starting point is infinite or NaN.");
    }

    // Check to see if the starting point is already a minimum.

    gradNormPart1Kernel->execute(numVariables, threadBlockSize);
    gradNormPart2Kernel->execute(threadBlockSize, threadBlockSize);
    if (downloadGradNormSync() <= tolerance) {
        return;
    }

    initializeDirKernel->execute(numVariables);
    for (int iteration = 1, end = 0;;) {
        // Prepare for a line search.

        if (mixedIsDouble) {
            lineSearchSetupKernel->setArg(8, energy);
        }
        else {
            lineSearchSetupKernel->setArg(8, (float) energy);
        }
        lineSearchSetupKernel->execute(numVariables);

        // Take line search steps.

        for (int count = 0;; count++) {
            lineSearchStepKernel->execute(numVariables, threadBlockSize);

            if (count) {
                int hostReturnFlag = downloadReturnFlagFinish();
                if (hostReturnFlag == 1) {
                    break; // Line search success
                }
                else if (hostReturnFlag == 0 || count >= maxLineSearchIterations) {
                    xPrev.copyTo(x);
                    gradPrev.copyTo(grad);
                    return; // Line search failure
                }
            }

            // Evaluate the energy and gradient at the new search point, then
            // decide if and how to continue the line search.

            evaluateGpu(context);
            downloadReturnValueStart();
            runLineSearchKernels();

            energy += downloadReturnValueFinish();
            if (!isfinite(energy)) {
                // Overflow on the GPU: try the CPU.

                energy = evaluateCpu(context);

                // We ran the line search kernels with an invalid gradient, so
                // we need to reset the state to before they ran and retry.

                lineSearchDataBackup.copyTo(lineSearchData);

                int hostReturnFlag = 2; // Continue the line search
                returnFlag.upload(&hostReturnFlag);

                // lineSearchDot will try to read any restraint energy from
                // returnValue, but it will be NaN from the failed GPU run, so
                // reset it to 0 (since the restraint energy from the CPU run
                // is included in the return value that will get uploaded).

                if (mixedIsDouble) {
                    double hostEnergy = 0;
                    returnValue.upload(&hostEnergy);
                }
                else {
                    float hostEnergy = 0;
                    returnValue.upload(&hostEnergy);
                }

                runLineSearchKernels();
            }


            downloadReturnFlagStart();
        }

        // Check for convergence or exceeding the maximum number of steps.

        if ((reporter != NULL && report(context, iteration)) || (maxIterations && maxIterations < iteration + 1)) {
            return;
        }

        // Do L-BFGS update of search direction.  Note that the equivalent of
        // gradNormPart1Kernel has already been executed in lineSearchStepKernel
        // if it was detected that the line search succeeded.

        gradNormPart2Kernel->execute(threadBlockSize, threadBlockSize);
        getDiffKernel->setArg(12, end);
        getDiffKernel->execute(numVariables, threadBlockSize);

        downloadReturnFlagStart();

        getScaleKernel->setArg(9, end);
        getScaleKernel->execute(threadBlockSize, threadBlockSize);

        int limit = min(numVectors, iteration++);
        if (++end >= numVectors) {
            end -= numVectors;
        }

        int vectorIndex = (end ? end : numVectors) - 1;
        reinitializeDirKernel->setArg(7, vectorIndex);
        reinitializeDirKernel->execute(numVariables);

        for (int vector = 0; vector < limit; vector++) {
            if (vector && --vectorIndex < 0) {
                vectorIndex += numVectors;
            }

            if (vector < limit - 1) {
                updateDirAlphaKernel->setArg(7, vectorIndex);
                updateDirAlphaKernel->execute(numVariables);
            }
        }

        scaleDirKernel->setArg(7, vectorIndex);
        scaleDirKernel->execute(numVariables);

        for (int vector = 0; vector < limit - 1; vector++) {
            // scaleDirKernel puts its first result in alpha[numVectors], so for the
            // first vector, load the result from here instead of alpha[vectorIndex]

            updateDirBetaKernel->setArg(7, vectorIndex);
            updateDirBetaKernel->setArg(8, vector ? vectorIndex : numVectors);
            updateDirBetaKernel->execute(numVariables);

            if (++vectorIndex >= numVectors) {
                vectorIndex -= numVectors;
            }
        }

        // If this is the first iteration, limit is 1, we did not go through the loop above,
        // and we need to read the output of scaleDirKernel directly from alpha[numVectors]

        updateDirFinalKernel->setArg(6, vectorIndex);
        updateDirFinalKernel->setArg(7, limit > 1 ? vectorIndex : numVectors);
        updateDirFinalKernel->execute(numVariables);

        if (downloadReturnFlagFinish()) {
            return;
        }
    }
}

void CommonMinimizeKernel::evaluateGpu(ContextImpl& context) {
    // Put the current positions in posq and compute virtual site positions.

    restorePosKernel->setArg(2, x);
    restorePosKernel->execute(numParticles);
    context.computeVirtualSites();

    // Evaluate the forces and energy for the desired interactions as well as
    // harmonic restraints to emulate the constraints.

    energy = context.calcForcesAndEnergy(true, true, forceGroups);
    if (numConstraints) {
        if (mixedIsDouble) {
            getConstraintEnergyForcesKernel->setArg(8, kRestraint);
        }
        else {
            getConstraintEnergyForcesKernel->setArg(8, (float) kRestraint);
        }
        getConstraintEnergyForcesKernel->execute(numConstraints);
    }

    // Convert the forces from fixed to floating point format.  If they are too
    // large, the energy in returnValue will be set to NaN to signal that the
    // results are invalid and we must fall back to computing forces on the CPU.

    convertForcesKernel->execute(numParticles);
}

double CommonMinimizeKernel::evaluateCpu(ContextImpl& context) {
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
                return NAN;
            }
        }
    }
    else {
        for (int i = 0; i < numVariables; i++) {
            if (!isfinite((float) hostGrad[i])) {
                return NAN;
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

void CommonMinimizeKernel::downloadReturnFlagStart() {
    downloadStartEvent->enqueue();
    cc.setCurrentQueue(downloadQueue);
    downloadStartEvent->queueWait(downloadQueue);
    returnFlag.download(pinnedMemory, false);
    downloadFinishEvent->enqueue();
    cc.restoreDefaultQueue();
}

void CommonMinimizeKernel::downloadReturnValueStart() {
    downloadStartEvent->enqueue();
    cc.setCurrentQueue(downloadQueue);
    downloadStartEvent->queueWait(downloadQueue);
    returnValue.download(pinnedMemory, false);
    downloadFinishEvent->enqueue();
    cc.restoreDefaultQueue();
}

int CommonMinimizeKernel::downloadReturnFlagFinish() {
    downloadFinishEvent->wait();
    return *(int*) pinnedMemory;
}

double CommonMinimizeKernel::downloadReturnValueFinish() {
    downloadFinishEvent->wait();
    if (mixedIsDouble) {
        return *(double*) pinnedMemory;
    }
    else {
        return (double) *(float*) pinnedMemory;
    }
}

double CommonMinimizeKernel::downloadReturnValueSync() {
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

double CommonMinimizeKernel::downloadGradNormSync() {
    if (mixedIsDouble) {
        double hostGradNorm;
        gradNorm.download(&hostGradNorm);
        return hostGradNorm;
    }
    else {
        float hostGradNorm;
        gradNorm.download(&hostGradNorm);
        return (double) hostGradNorm;
    }
}

void CommonMinimizeKernel::runLineSearchKernels() {
    if (mixedIsDouble) {
        lineSearchDotKernel->setArg(6, energy);
    }
    else {
        lineSearchDotKernel->setArg(6, (float) energy);
    }
    lineSearchDotKernel->execute(numVariables);
    lineSearchContinueKernel->execute(1);
}
