
/* Portions copyright (c) 2011-2016 Stanford University and Simbios.
 * Contributors: Peter Eastman
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "SimTKOpenMMUtilities.h"
#include "ReferenceVirtualSites.h"
#include "ReferenceCustomDynamics.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/ForceImpl.h"
#include "lepton/Operation.h"
#include "lepton/ParsedExpression.h"
#include "lepton/Parser.h"
#include <set>
#include <sstream>

using namespace std;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceCustomDynamics constructor

   @param numberOfAtoms  number of atoms
   @param integrator     the integrator definition to use

   --------------------------------------------------------------------------------------- */

ReferenceCustomDynamics::ReferenceCustomDynamics(int numberOfAtoms, const CustomIntegrator& integrator) : 
           ReferenceDynamics(numberOfAtoms, integrator.getStepSize(), 0.0), integrator(integrator) {
    sumBuffer.resize(numberOfAtoms);
    oldPos.resize(numberOfAtoms);
    stepType.resize(integrator.getNumComputations());
    stepVariable.resize(integrator.getNumComputations());
    for (int i = 0; i < integrator.getNumComputations(); i++) {
        string expression;
        integrator.getComputationStep(i, stepType[i], stepVariable[i], expression);
    }
    kineticEnergyExpression = Lepton::Parser::parse(integrator.getKineticEnergyExpression()).optimize().createCompiledExpression();
    expressionSet.registerExpression(kineticEnergyExpression);
    kineticEnergyNeedsForce = false;
    if (kineticEnergyExpression.getVariables().find("f") != kineticEnergyExpression.getVariables().end())
        kineticEnergyNeedsForce = true;
}

/**---------------------------------------------------------------------------------------

   ReferenceCustomDynamics destructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomDynamics::~ReferenceCustomDynamics() {
}

void ReferenceCustomDynamics::initialize(ContextImpl& context, vector<RealOpenMM>& masses, map<string, RealOpenMM>& globals) {
    // Some initialization can't be done in the constructor, since we need a ContextImpl from which to get the list of
    // Context parameters.  Instead, we do it the first time update() or computeKineticEnergy() is called.

    int numSteps = stepType.size();
    vector<int> forceGroup;
    vector<vector<Lepton::ParsedExpression> > expressions;
    CustomIntegratorUtilities::analyzeComputations(context, integrator, expressions, comparisons, blockEnd, invalidatesForces, needsForces, needsEnergy, computeBothForceAndEnergy, forceGroup);
    stepExpressions.resize(expressions.size());
    for (int i = 0; i < numSteps; i++) {
        stepExpressions[i].resize(expressions[i].size());
        for (int j = 0; j < (int) expressions[i].size(); j++) {
            stepExpressions[i][j] = expressions[i][j].createCompiledExpression();
            expressionSet.registerExpression(stepExpressions[i][j]);
        }
        if (stepType[i] == CustomIntegrator::WhileBlockStart)
            blockEnd[blockEnd[i]] = i; // Record where to branch back to.
    }

    // Record the variable names and flags for the force and energy in each step.

    forceGroupFlags.resize(numSteps, -1);
    fIndex = expressionSet.getVariableIndex("f");
    energyIndex = expressionSet.getVariableIndex("energy");
    forceVariableIndex.resize(numSteps, fIndex);
    energyVariableIndex.resize(numSteps, energyIndex);
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
    for (int i = 0; i < numSteps; i++) {
        if (needsForces[i] && forceGroup[i] > -1)
            forceVariableIndex[i] = expressionSet.getVariableIndex(forceGroupName[forceGroup[i]]);
        if (needsEnergy[i] && forceGroup[i] > -1)
            energyVariableIndex[i] = expressionSet.getVariableIndex(energyGroupName[forceGroup[i]]);
        if (forceGroup[i] > -1)
            forceGroupFlags[i] = 1<<forceGroup[i];
    }

    // Build the list of inverse masses.

    int numberOfAtoms = masses.size();
    inverseMasses.resize(numberOfAtoms);
    for (int i = 0; i < numberOfAtoms; i++) {
        if (masses[i] == 0.0)
            inverseMasses[i] = 0.0;
        else
            inverseMasses[i] = 1.0/masses[i];
    }

    // Record indices of other variables.

    xIndex = expressionSet.getVariableIndex("x");
    vIndex = expressionSet.getVariableIndex("v");
    mIndex = expressionSet.getVariableIndex("m");
    gaussianIndex = expressionSet.getVariableIndex("gaussian");
    uniformIndex = expressionSet.getVariableIndex("uniform");
    for (int i = 0; i < integrator.getNumPerDofVariables(); i++)
        perDofVariableIndex.push_back(expressionSet.getVariableIndex(integrator.getPerDofVariableName(i)));
    for (int i = 0; i < stepVariable.size(); i++)
        stepVariableIndex.push_back(expressionSet.getVariableIndex(stepVariable[i]));
}

/**---------------------------------------------------------------------------------------

   Update -- driver routine for performing Custom dynamics update of coordinates
   and velocities

   @param context             the context this integrator is updating
   @param numberOfAtoms       number of atoms
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param masses              atom masses
   @param globals             a map containing values of global variables
   @param forcesAreValid      whether the current forces are valid or need to be recomputed

   --------------------------------------------------------------------------------------- */

void ReferenceCustomDynamics::update(ContextImpl& context, int numberOfAtoms, vector<RealVec>& atomCoordinates,
                                     vector<RealVec>& velocities, vector<RealVec>& forces, vector<RealOpenMM>& masses,
                                     map<string, RealOpenMM>& globals, vector<vector<RealVec> >& perDof, bool& forcesAreValid, RealOpenMM tolerance) {
    if (invalidatesForces.size() == 0)
        initialize(context, masses, globals);
    int numSteps = stepType.size();
    globals.insert(context.getParameters().begin(), context.getParameters().end());
    for (map<string, RealOpenMM>::const_iterator iter = globals.begin(); iter != globals.end(); ++iter)
        expressionSet.setVariable(expressionSet.getVariableIndex(iter->first), iter->second);
    oldPos = atomCoordinates;
    
    // Loop over steps and execute them.
    
    for (int step = 0; step < numSteps; ) {
        if ((needsForces[step] || needsEnergy[step]) && (!forcesAreValid || context.getLastForceGroups() != forceGroupFlags[step])) {
            // Recompute forces and/or energy.
            
            bool computeForce = needsForces[step] || computeBothForceAndEnergy[step];
            bool computeEnergy = needsEnergy[step] || computeBothForceAndEnergy[step];
            recordChangedParameters(context, globals);
            RealOpenMM e = context.calcForcesAndEnergy(computeForce, computeEnergy, forceGroupFlags[step]);
            if (computeEnergy)
                energy = e;
            forcesAreValid = true;
        }
        expressionSet.setVariable(energyVariableIndex[step], energy);
        
        // Execute the step.

        int nextStep = step+1;
        switch (stepType[step]) {
            case CustomIntegrator::ComputeGlobal: {
                expressionSet.setVariable(uniformIndex, SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber());
                expressionSet.setVariable(gaussianIndex, SimTKOpenMMUtilities::getNormallyDistributedRandomNumber());
                RealOpenMM result = stepExpressions[step][0].evaluate();
                globals[stepVariable[step]] = result;
                expressionSet.setVariable(stepVariableIndex[step], result);
                break;
            }
            case CustomIntegrator::ComputePerDof: {
                vector<RealVec>* results = NULL;
                if (stepVariableIndex[step] == xIndex)
                    results = &atomCoordinates;
                else if (stepVariableIndex[step] == vIndex)
                    results = &velocities;
                else {
                    for (int j = 0; j < integrator.getNumPerDofVariables(); j++)
                        if (stepVariableIndex[step] == perDofVariableIndex[j])
                            results = &perDof[j];
                }
                if (results == NULL)
                    throw OpenMMException("Illegal per-DOF output variable: "+stepVariable[step]);
                computePerDof(numberOfAtoms, *results, atomCoordinates, velocities, forces, masses, perDof, stepExpressions[step][0], forceVariableIndex[step]);
                break;
            }
            case CustomIntegrator::ComputeSum: {
                computePerDof(numberOfAtoms, sumBuffer, atomCoordinates, velocities, forces, masses, perDof, stepExpressions[step][0], forceVariableIndex[step]);
                RealOpenMM sum = 0.0;
                for (int j = 0; j < numberOfAtoms; j++)
                    if (masses[j] != 0.0)
                        sum += sumBuffer[j][0]+sumBuffer[j][1]+sumBuffer[j][2];
                globals[stepVariable[step]] = sum;
                expressionSet.setVariable(stepVariableIndex[step], sum);
                break;
            }
            case CustomIntegrator::ConstrainPositions: {
                getReferenceConstraintAlgorithm()->apply(oldPos, atomCoordinates, inverseMasses, tolerance);
                oldPos = atomCoordinates;
                break;
            }
            case CustomIntegrator::ConstrainVelocities: {
                getReferenceConstraintAlgorithm()->applyToVelocities(oldPos, velocities, inverseMasses, tolerance);
                break;
            }
            case CustomIntegrator::UpdateContextState: {
                recordChangedParameters(context, globals);
                context.updateContextState();
                globals.insert(context.getParameters().begin(), context.getParameters().end());
                for (map<string, RealOpenMM>::const_iterator iter = globals.begin(); iter != globals.end(); ++iter)
                    expressionSet.setVariable(expressionSet.getVariableIndex(iter->first), iter->second);
                break;
            }
            case CustomIntegrator::IfBlockStart: {
                if (!evaluateCondition(step))
                    nextStep = blockEnd[step]+1;
                break;
            }
            case CustomIntegrator::WhileBlockStart: {
                if (!evaluateCondition(step))
                    nextStep = blockEnd[step]+1;
                break;
            }
            case CustomIntegrator::BlockEnd: {
                if (blockEnd[step] != -1)
                    nextStep = blockEnd[step]; // Return to the start of a while block.
                break;
            }
        }
        if (invalidatesForces[step])
            forcesAreValid = false;
        step = nextStep;
    }
    ReferenceVirtualSites::computePositions(context.getSystem(), atomCoordinates);
    incrementTimeStep();
    recordChangedParameters(context, globals);
}

void ReferenceCustomDynamics::computePerDof(int numberOfAtoms, vector<RealVec>& results, const vector<RealVec>& atomCoordinates,
              const vector<RealVec>& velocities, const vector<RealVec>& forces, const vector<RealOpenMM>& masses,
              const vector<vector<RealVec> >& perDof, const Lepton::CompiledExpression& expression, int forceIndex) {
    // Loop over all degrees of freedom.

    for (int i = 0; i < numberOfAtoms; i++) {
        if (masses[i] != 0.0) {
            expressionSet.setVariable(mIndex, masses[i]);
            for (int j = 0; j < 3; j++) {
                // Compute the expression.

                expressionSet.setVariable(xIndex, atomCoordinates[i][j]);
                expressionSet.setVariable(vIndex, velocities[i][j]);
                expressionSet.setVariable(forceIndex, forces[i][j]);
                expressionSet.setVariable(uniformIndex, SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber());
                expressionSet.setVariable(gaussianIndex, SimTKOpenMMUtilities::getNormallyDistributedRandomNumber());
                for (int k = 0; k < (int) perDof.size(); k++)
                    expressionSet.setVariable(perDofVariableIndex[k], perDof[k][i][j]);
                results[i][j] = expression.evaluate();
            }
        }
    }
}

bool ReferenceCustomDynamics::evaluateCondition(int step) {
    expressionSet.setVariable(uniformIndex, SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber());
    expressionSet.setVariable(gaussianIndex, SimTKOpenMMUtilities::getNormallyDistributedRandomNumber());
    double lhs = stepExpressions[step][0].evaluate();
    double rhs = stepExpressions[step][1].evaluate();
    switch (comparisons[step]) {
        case CustomIntegratorUtilities::EQUAL:
            return (lhs == rhs);
        case CustomIntegratorUtilities::LESS_THAN:
            return (lhs < rhs);
        case CustomIntegratorUtilities::GREATER_THAN:
            return (lhs > rhs);
        case CustomIntegratorUtilities::NOT_EQUAL:
            return (lhs != rhs);
        case CustomIntegratorUtilities::LESS_THAN_OR_EQUAL:
            return (lhs <= rhs);
        case CustomIntegratorUtilities::GREATER_THAN_OR_EQUAL:
            return (lhs >= rhs);
    }
    throw OpenMMException("ReferenceCustomDynamics: Invalid comparison operator");
}

/**
 * Check which context parameters have changed and register them with the context.
 */
void ReferenceCustomDynamics::recordChangedParameters(OpenMM::ContextImpl& context, std::map<std::string, RealOpenMM>& globals) {
    for (map<string, double>::const_iterator iter = context.getParameters().begin(); iter != context.getParameters().end(); ++iter) {
        string name = iter->first;
        double value = globals[name];
        if (value != iter->second)
            context.setParameter(name, globals[name]);
    }
}

/**---------------------------------------------------------------------------------------

   Compute the kinetic energy of the system.

   @param context             the context this integrator is updating
   @param numberOfAtoms       number of atoms
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param masses              atom masses
   @param globals             a map containing values of global variables
   @param perDof              the values of per-DOF variables
   @param forcesAreValid      whether the current forces are valid or need to be recomputed

   --------------------------------------------------------------------------------------- */

double ReferenceCustomDynamics::computeKineticEnergy(OpenMM::ContextImpl& context, int numberOfAtoms, std::vector<OpenMM::RealVec>& atomCoordinates,
        std::vector<OpenMM::RealVec>& velocities, std::vector<OpenMM::RealVec>& forces, std::vector<RealOpenMM>& masses,
        std::map<std::string, RealOpenMM>& globals, std::vector<std::vector<OpenMM::RealVec> >& perDof, bool& forcesAreValid) {
    if (invalidatesForces.size() == 0)
        initialize(context, masses, globals);
    globals.insert(context.getParameters().begin(), context.getParameters().end());
    for (map<string, RealOpenMM>::const_iterator iter = globals.begin(); iter != globals.end(); ++iter)
        expressionSet.setVariable(expressionSet.getVariableIndex(iter->first), iter->second);
    if (kineticEnergyNeedsForce) {
        energy = context.calcForcesAndEnergy(true, true, -1);
        forcesAreValid = true;
    }
    computePerDof(numberOfAtoms, sumBuffer, atomCoordinates, velocities, forces, masses, perDof, kineticEnergyExpression, fIndex);
    RealOpenMM sum = 0.0;
    for (int j = 0; j < numberOfAtoms; j++)
        if (masses[j] != 0.0)
            sum += sumBuffer[j][0]+sumBuffer[j][1]+sumBuffer[j][2];
    return sum;
}
