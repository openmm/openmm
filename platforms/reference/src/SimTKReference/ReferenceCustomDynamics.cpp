
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
using namespace Lepton;

class ReferenceCustomDynamics::DerivFunction : public CustomFunction {
public:
    DerivFunction(map<string, double>& energyParamDerivs, const string& param) : energyParamDerivs(energyParamDerivs), param(param) {
    }
    int getNumArguments() const {
        return 0;
    }
    double evaluate(const double* arguments) const {
        return energyParamDerivs[param];
    }
    double evaluateDerivative(const double* arguments, const int* derivOrder) const {
        return 0;
    }
    CustomFunction* clone() const {
        return new DerivFunction(energyParamDerivs, param);
    }
private:
    map<string, double>& energyParamDerivs;
    string param;
};

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
}

/**---------------------------------------------------------------------------------------

   ReferenceCustomDynamics destructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomDynamics::~ReferenceCustomDynamics() {
}

void ReferenceCustomDynamics::initialize(ContextImpl& context, vector<RealOpenMM>& masses, map<string, RealOpenMM>& globals) {
    // Some initialization can't be done in the constructor, since we need a ContextImpl from which to get the list of
    // Context parameters.  Instead, we do it the first time update() or computeKineticEnergy() is called.

    std::map<std::string, double*> variableLocations;
    variableLocations["x"] = &x;
    variableLocations["v"] = &v;
    variableLocations["m"] = &m;
    variableLocations["f"] = &f;
    variableLocations["energy"] = &energy;
    variableLocations["gaussian"] = &gaussian;
    variableLocations["uniform"] = &uniform;
    perDofVariable.resize(integrator.getNumPerDofVariables());
    for (int i = 0; i < integrator.getNumPerDofVariables(); i++)
        variableLocations[integrator.getPerDofVariableName(i)] = &perDofVariable[i];
    for (int i = 0; i < 32; i++) {
        stringstream fname;
        fname << "f" << i;
        variableLocations[fname.str()] = &f;
        stringstream ename;
        ename << "energy" << i;
        variableLocations[ename.str()] = &energy;
    }
    
    // Parse the expressions.
    
    int numSteps = stepType.size();
    vector<int> forceGroup;
    vector<vector<ParsedExpression> > expressions;
    CustomIntegratorUtilities::analyzeComputations(context, integrator, expressions, comparisons, blockEnd, invalidatesForces, needsForces, needsEnergy, computeBothForceAndEnergy, forceGroup);
    stepExpressions.resize(expressions.size());
    for (int i = 0; i < numSteps; i++) {
        stepExpressions[i].resize(expressions[i].size());
        for (int j = 0; j < (int) expressions[i].size(); j++) {
            stepExpressions[i][j] = ParsedExpression(replaceDerivFunctions(expressions[i][j].getRootNode(), context)).createCompiledExpression();
            stepExpressions[i][j].setVariableLocations(variableLocations);
            expressionSet.registerExpression(stepExpressions[i][j]);
        }
        if (stepType[i] == CustomIntegrator::WhileBlockStart)
            blockEnd[blockEnd[i]] = i; // Record where to branch back to.
    }
    kineticEnergyExpression = Parser::parse(integrator.getKineticEnergyExpression()).optimize().createCompiledExpression();
    kineticEnergyExpression.setVariableLocations(variableLocations);
    expressionSet.registerExpression(kineticEnergyExpression);
    kineticEnergyNeedsForce = false;
    if (kineticEnergyExpression.getVariables().find("f") != kineticEnergyExpression.getVariables().end())
        kineticEnergyNeedsForce = true;

    // Record the force group flags for each step.

    forceGroupFlags.resize(numSteps, -1);
    for (int i = 0; i < numSteps; i++)
        if (forceGroup[i] > -1)
            forceGroupFlags[i] = 1<<forceGroup[i];

    // Build the list of inverse masses.

    int numberOfAtoms = masses.size();
    inverseMasses.resize(numberOfAtoms);
    for (int i = 0; i < numberOfAtoms; i++) {
        if (masses[i] == 0.0)
            inverseMasses[i] = 0.0;
        else
            inverseMasses[i] = 1.0/masses[i];
    }

    // Record indices of variables.

    xIndex = expressionSet.getVariableIndex("x");
    vIndex = expressionSet.getVariableIndex("v");
    for (int i = 0; i < integrator.getNumPerDofVariables(); i++)
        perDofVariableIndex.push_back(expressionSet.getVariableIndex(integrator.getPerDofVariableName(i)));
    for (int i = 0; i < stepVariable.size(); i++)
        stepVariableIndex.push_back(expressionSet.getVariableIndex(stepVariable[i]));
}

ExpressionTreeNode ReferenceCustomDynamics::replaceDerivFunctions(const ExpressionTreeNode& node, ContextImpl& context) {
    const Operation& op = node.getOperation();
    if (op.getId() == Operation::CUSTOM && op.getName() == "deriv") {
        string param = node.getChildren()[1].getOperation().getName();
        if (context.getParameters().find(param) == context.getParameters().end())
            throw OpenMMException("The second argument to deriv() must be a context parameter");
        return ExpressionTreeNode(new Operation::Custom("deriv", new DerivFunction(energyParamDerivs, param)));
    }
    else {
        vector<ExpressionTreeNode> children;
        for (int i = 0; i < (int) node.getChildren().size(); i++)
            children.push_back(replaceDerivFunctions(node.getChildren()[i], context));
        return ExpressionTreeNode(op.clone(), children);
    }
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
            if (computeEnergy) {
                energy = e;
                context.getEnergyParameterDerivatives(energyParamDerivs);
            }
            forcesAreValid = true;
        }
        
        // Execute the step.

        int nextStep = step+1;
        switch (stepType[step]) {
            case CustomIntegrator::ComputeGlobal: {
                uniform = SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber();
                gaussian = SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
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
                computePerDof(numberOfAtoms, *results, atomCoordinates, velocities, forces, masses, perDof, stepExpressions[step][0]);
                break;
            }
            case CustomIntegrator::ComputeSum: {
                computePerDof(numberOfAtoms, sumBuffer, atomCoordinates, velocities, forces, masses, perDof, stepExpressions[step][0]);
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
              const vector<vector<RealVec> >& perDof, const CompiledExpression& expression) {
    // Loop over all degrees of freedom.

    for (int i = 0; i < numberOfAtoms; i++) {
        if (masses[i] != 0.0) {
            m = masses[i];
            for (int j = 0; j < 3; j++) {
                // Compute the expression.

                x = atomCoordinates[i][j];
                v = velocities[i][j];
                f = forces[i][j];
                uniform = SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber();
                gaussian = SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
                for (int k = 0; k < (int) perDof.size(); k++)
                    perDofVariable[k] = perDof[k][i][j];
                results[i][j] = expression.evaluate();
            }
        }
    }
}

bool ReferenceCustomDynamics::evaluateCondition(int step) {
    uniform = SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber();
    gaussian = SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
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
    computePerDof(numberOfAtoms, sumBuffer, atomCoordinates, velocities, forces, masses, perDof, kineticEnergyExpression);
    RealOpenMM sum = 0.0;
    for (int j = 0; j < numberOfAtoms; j++)
        if (masses[j] != 0.0)
            sum += sumBuffer[j][0]+sumBuffer[j][1]+sumBuffer[j][2];
    return sum;
}
