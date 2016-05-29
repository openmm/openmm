/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2015 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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

#include "openmm/internal/CustomIntegratorUtilities.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ForceImpl.h"
#include "lepton/Operation.h"
#include "lepton/Parser.h"
#include <set>
#include <sstream>

using namespace OpenMM;
using namespace std;

void CustomIntegratorUtilities::parseCondition(const string& expression, string& lhs, string& rhs, Comparison& comparison) {
    string operators[] = {"=", "<", ">", "!=", "<=", ">="};
    for (int i = 5; i >= 0; i--) {
        int index = expression.find(operators[i]);
        if (index != string::npos) {
            lhs = expression.substr(0, index);
            rhs = expression.substr(index+operators[i].size());
            comparison = Comparison(i);
            return;
        }
    }
    throw OpenMMException("No comparison operator found in condition: "+expression);
}

bool CustomIntegratorUtilities::usesVariable(const Lepton::ExpressionTreeNode& node, const string& variable) {
    const Lepton::Operation& op = node.getOperation();
    if (op.getId() == Lepton::Operation::VARIABLE && op.getName() == variable)
        return true;
    for (int i = 0; i < (int) node.getChildren().size(); i++)
        if (usesVariable(node.getChildren()[i], variable))
            return true;
    return false;
}

bool CustomIntegratorUtilities::usesVariable(const Lepton::ParsedExpression& expression, const string& variable) {
    return usesVariable(expression.getRootNode(), variable);
}

void CustomIntegratorUtilities::analyzeComputations(const ContextImpl& context, const CustomIntegrator& integrator, vector<vector<Lepton::ParsedExpression> >& expressions,
            vector<Comparison>& comparisons, vector<int>& blockEnd, vector<bool>& invalidatesForces, vector<bool>& needsForces, vector<bool>& needsEnergy,
            vector<bool>& computeBoth, vector<int>& forceGroup) {
    int numSteps = integrator.getNumComputations();
    expressions.resize(numSteps);
    comparisons.resize(numSteps);
    invalidatesForces.resize(numSteps, false);
    needsForces.resize(numSteps, false);
    needsEnergy.resize(numSteps, false);
    computeBoth.resize(numSteps, false);
    forceGroup.resize(numSteps, -2);
    vector<CustomIntegrator::ComputationType> stepType(numSteps);
    vector<string> stepVariable(numSteps);

    // Parse the expressions.

    for (int step = 0; step < numSteps; step++) {
        string expression;
        integrator.getComputationStep(step, stepType[step], stepVariable[step], expression);
        if (stepType[step] == CustomIntegrator::IfBlockStart || stepType[step] == CustomIntegrator::WhileBlockStart) {
            // This step involves a condition.

            string lhs, rhs;
            parseCondition(expression, lhs, rhs, comparisons[step]);
            expressions[step].push_back(Lepton::Parser::parse(lhs).optimize());
            expressions[step].push_back(Lepton::Parser::parse(rhs).optimize());
        }
        else if (expression.size() > 0)
            expressions[step].push_back(Lepton::Parser::parse(expression).optimize());
    }

    // Identify which steps invalidate the forces.

    set<string> affectsForce;
    affectsForce.insert("x");
    for (vector<ForceImpl*>::const_iterator iter = context.getForceImpls().begin(); iter != context.getForceImpls().end(); ++iter) {
        const map<string, double> params = (*iter)->getDefaultParameters();
        for (map<string, double>::const_iterator param = params.begin(); param != params.end(); ++param)
            affectsForce.insert(param->first);
    }
    for (int i = 0; i < numSteps; i++)
        invalidatesForces[i] = (stepType[i] == CustomIntegrator::ConstrainPositions || affectsForce.find(stepVariable[i]) != affectsForce.end());

    // Make a list of which steps require valid forces or energy to be known.

    vector<string> forceGroupName;
    vector<string> energyGroupName;
    for (int i = 0; i < 32; i++) {
        stringstream fname;
        fname << "f" << i;
        forceGroupName.push_back(fname.str());
        stringstream ename;
        ename << "energy" << i;
        energyGroupName.push_back(ename.str());
    }
    for (int step = 0; step < numSteps; step++) {
        for (int expr = 0; expr < expressions[step].size(); expr++) {
            if (usesVariable(expressions[step][expr], "f")) {
                needsForces[step] = true;
                forceGroup[step] = -1;
            }
            if (usesVariable(expressions[step][expr], "energy")) {
                needsEnergy[step] = true;
                forceGroup[step] = -1;
            }
            for (int i = 0; i < 32; i++) {
                if (usesVariable(expressions[step][expr], forceGroupName[i])) {
                    if (forceGroup[step] != -2)
                        throw OpenMMException("A single computation step cannot depend on multiple force groups");
                    needsForces[step] = true;
                    forceGroup[step] = i;
                }
                if (usesVariable(expressions[step][expr], energyGroupName[i])) {
                    if (forceGroup[step] != -2)
                        throw OpenMMException("A single computation step cannot depend on multiple force groups");
                    needsEnergy[step] = true;
                    forceGroup[step] = i;
                }
            }
        }
    }
    for (int step = numSteps-2; step >= 0; step--)
        if (forceGroup[step] == -2)
            forceGroup[step] = forceGroup[step+1];

    // Find the end point of each block.

    vector<int> blockStart;
    blockEnd.resize(numSteps, -1);
    for (int step = 0; step < numSteps; step++) {
        if (stepType[step] == CustomIntegrator::IfBlockStart || stepType[step] == CustomIntegrator::WhileBlockStart)
            blockStart.push_back(step);
        else if (stepType[step] == CustomIntegrator::BlockEnd) {
            if (blockStart.size() == 0) {
                stringstream error("CustomIntegrator: Unexpected end of block at computation ");
                error << step;
                throw OpenMMException(error.str());
            }
            blockEnd[blockStart.back()] = step;
            blockStart.pop_back();
        }
    }
    if (blockStart.size() > 0)
        throw OpenMMException("CustomIntegrator: Missing EndBlock");

    // If a step requires either forces or energy, and a later step will require the other one, it's most efficient
    // to compute both at the same time.  Figure out whether we should do that.  In principle it's easy: step through
    // the sequence of computations and see if the other one is used before the next time they get invalidated.
    // Unfortunately, flow control makes this much more complicated, because there are many possible paths to
    // consider.
    //
    // The cost of computing both when we really only needed one is much less than the cost of computing only one,
    // then later finding we need to compute the other separately.  So we always err on the side of computing both.
    // If there is any possible path that would lead to us needing it, go ahead and compute it.
    //
    // So we need to enumerate all possible paths.  For each "if" block, there are two possibilities: execute it
    // or don't.  For each "while" block there are three possibilities: don't execute it; execute it and then
    // continue on; or execute it and then jump back to the beginning.  I'm assuming the number of blocks will
    // always remain small.  Otherwise, this could become very expensive!

    vector<int> jumps(numSteps, -1);
    vector<int> stepsInPath;
    enumeratePaths(0, stepsInPath, jumps, blockEnd, stepType, needsForces, needsEnergy, invalidatesForces, forceGroup, computeBoth);
}

void CustomIntegratorUtilities::enumeratePaths(int firstStep, vector<int> steps, vector<int> jumps, const vector<int>& blockEnd,
            const vector<CustomIntegrator::ComputationType>& stepType, const vector<bool>& needsForces, const vector<bool>& needsEnergy,
            const vector<bool>& invalidatesForces, const vector<int>& forceGroup, vector<bool>& computeBoth) {
    int step = firstStep;
    int numSteps = stepType.size();
    while (step < numSteps) {
        steps.push_back(step);
        if (jumps[step] > 0) {
            // Follow the jump and remove it from the list.

            int nextStep = jumps[step];
            jumps[step] = -1;
            step = nextStep;
        }
        else if (stepType[step] == CustomIntegrator::IfBlockStart) {
            // Consider skipping the block.

            enumeratePaths(blockEnd[step]+1, steps, jumps, blockEnd, stepType, needsForces, needsEnergy, invalidatesForces, forceGroup, computeBoth);

            // Continue on to execute the block.

            step++;
        }
        else if (stepType[step] == CustomIntegrator::WhileBlockStart && jumps[step] != -2) {
            // Consider skipping the block.

            enumeratePaths(blockEnd[step]+1, steps, jumps, blockEnd, stepType, needsForces, needsEnergy, invalidatesForces, forceGroup, computeBoth);

            // Consider executing the block once.

            enumeratePaths(step+1, steps, jumps, blockEnd, stepType, needsForces, needsEnergy, invalidatesForces, forceGroup, computeBoth);

            // Continue on to execute the block twice.

            jumps[step] = -2; // Mark this "while" block as already processed.
            jumps[blockEnd[step]] = step;
            step++;
        }
        else
            step++;
    }
    analyzeForceComputationsForPath(steps, needsForces, needsEnergy, invalidatesForces, forceGroup, computeBoth);
}

void CustomIntegratorUtilities::analyzeForceComputationsForPath(vector<int>& steps, const vector<bool>& needsForces, const vector<bool>& needsEnergy,
            const vector<bool>& invalidatesForces, const vector<int>& forceGroup, vector<bool>& computeBoth) {
    vector<int> candidatePoints;
    int currentGroup = -1;
    for (int i = 0; i < (int) steps.size(); i++) {
        int step = steps[i];
        if (invalidatesForces[step] || ((needsForces[step] || needsEnergy[step]) && forceGroup[step] != currentGroup)) {
            // Forces and energies are invalidated at this step, or it changes to a different force group,
            // so anything from this point on won't affect what we do at earlier steps.

            candidatePoints.clear();
        }
        if (needsForces[step] || needsEnergy[step]) {
            // See if this step affects what we do at earlier points.

            for (int j = 0; j < (int) candidatePoints.size(); j++) {
                int candidate = candidatePoints[j];
                if ((needsForces[candidate] && needsEnergy[step]) || (needsEnergy[candidate] && needsForces[step]))
                    computeBoth[candidate] = true;
            }

            // Add this to the list of candidates that might be affected by later steps.

            candidatePoints.push_back(step);
            currentGroup = forceGroup[step];
        }
    }
}