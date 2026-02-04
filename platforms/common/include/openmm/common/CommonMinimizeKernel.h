#ifndef OPENMM_COMMONMINIMIZEKERNEL_H_
#define OPENMM_COMMONMINIMIZEKERNEL_H_

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

#include "openmm/LocalEnergyMinimizer.h"
#include "openmm/kernels.h"
#include "openmm/common/ComputeContext.h"

namespace OpenMM {

/**
 * This kernel performs local energy minimization.
 */
class CommonMinimizeKernel : public MinimizeKernel {
public:
    CommonMinimizeKernel(std::string name, const Platform& platform, ComputeContext& cc) : MinimizeKernel(name, platform), cc(cc), isSetup(false), cpuContext(NULL), cpuIntegrator(1) {
    }
    ~CommonMinimizeKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     */
    void initialize(const System& system);
    /**
     * Perform local energy minimization.
     *
     * @param context        the context with which to perform the minimization
     * @param tolerance      limiting root-mean-square value of all force components in kJ/mol/nm for convergence
     * @param maxIterations  the maximum number of iterations to perform, or 0 to continue until convergence
     * @param reporter       an optional reporter to invoke after each iteration of minimization
     */
    void execute(ContextImpl& context, double tolerance, int maxIterations, MinimizationReporter* reporter);
private:
    static const double minConstraintTol, kRestraintScale, prevMaxErrorInit, kRestraintScaleUp, constraintTolScale;
    static const double fTol, wolfeParam, stepScaleDown, stepScaleUp, minStep, maxStep;
    static const int numVectors, maxLineSearchIterations;

    void setup(ContextImpl& context);
    void lbfgs(ContextImpl& context);
    bool lineSearch(ContextImpl& context, double& step);
    double evaluate(ContextImpl& context, bool& overflow);
    double evaluateCpu(ContextImpl& context, bool& overflow);
    bool report(ContextImpl& context, int iteration);
    int downloadReturnFlag();
    double downloadReturnValue();

    ComputeContext& cc;

    int numParticles, numVariables, numConstraints;

    std::vector<Vec3> hostPositions;
    std::vector<double> hostX;
    std::vector<double> hostGrad;
    std::vector<mm_int2> hostConstraintIndices;
    std::vector<double> hostConstraintDistances;

    bool isSetup, mixedIsDouble;
    int elementSize, threadBlockSize, numVariableBlocks, numConstraintBlocks;

    int forceGroups;
    double constraintTol;

    double tolerance;
    int maxIterations;
    MinimizationReporter* reporter;

    double kRestraint, energy;

    ComputeArray constraintIndices, constraintDistances;
    ComputeArray xInit, x, xPrev, grad, gradPrev, dir;
    ComputeArray alpha, scale, xDiff, gradDiff;
    ComputeArray reduceBuffer, returnFlag, returnValue;

    ComputeKernel recordInitialPosKernel;
    ComputeKernel restorePosKernel;
    ComputeKernel convertForcesKernel;
    ComputeKernel getConstraintEnergyForcesKernel, reduceConstraintEnergyKernel;
    ComputeKernel getConstraintErrorKernel;
    ComputeKernel initializeDirKernel;
    ComputeKernel gradNormPart1Kernel, gradNormPart2Kernel;
    ComputeKernel getDiffKernel;
    ComputeKernel getScaleKernel;
    ComputeKernel reinitializeDirKernel;
    ComputeKernel updateDirAlphaKernel;
    ComputeKernel scaleDirKernel;
    ComputeKernel updateDirBetaKernel;
    ComputeKernel updateDirFinalKernel;
    ComputeKernel lineSearchStepKernel;
    ComputeKernel lineSearchDotPart1Kernel, lineSearchDotPart2Kernel;

    Context* cpuContext;
    VerletIntegrator cpuIntegrator;
};

} // namespace OpenMM

#endif // OPENMM_COMMONMINIMIZEKERNEL_H_
