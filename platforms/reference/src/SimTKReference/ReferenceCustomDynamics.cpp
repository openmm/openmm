
/* Portions copyright (c) 2011-2017 Stanford University and Simbios.
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
#include "ReferenceTabulatedFunction.h"
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

/**
 * Determine whether a parsed expression involves any vector functions.
 */
static bool isVectorExpression(const ExpressionTreeNode& node) {
    const Lepton::Operation& op = node.getOperation();
    if (op.getId() == Lepton::Operation::CUSTOM)
        if (op.getName() == "dot" || op.getName() == "cross" || op.getName() == "vector" || op.getName() == "_x" || op.getName() == "_y" || op.getName() == "_z")
            return true;
    for (auto& child : node.getChildren())
        if (isVectorExpression(child))
            return true;
    return false;
}

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

void ReferenceCustomDynamics::initialize(ContextImpl& context, vector<double>& masses, map<string, double>& globals) {
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
    
    // Create custom functions for the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < integrator.getNumTabulatedFunctions(); i++)
        functions[integrator.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(integrator.getTabulatedFunction(i));
    
    // Parse the expressions.
    
    int numSteps = stepType.size();
    vector<int> forceGroup;
    vector<vector<ParsedExpression> > expressions;
    CustomIntegratorUtilities::analyzeComputations(context, integrator, expressions, comparisons, blockEnd, invalidatesForces, needsForces, needsEnergy, computeBothForceAndEnergy, forceGroup, functions);
    stepExpressions.resize(expressions.size());
    stepVectorExpressions.resize(expressions.size());
    for (int i = 0; i < numSteps; i++) {
        stepExpressions[i].resize(expressions[i].size());
        for (int j = 0; j < (int) expressions[i].size(); j++) {
            ParsedExpression parsed(replaceDerivFunctions(expressions[i][j].getRootNode(), context));
            if (isVectorExpression(parsed.getRootNode()))
                stepVectorExpressions[i].push_back(VectorExpression(parsed));
            else {
                stepExpressions[i][j] = parsed.createCompiledExpression();
                stepExpressions[i][j].setVariableLocations(variableLocations);
                expressionSet.registerExpression(stepExpressions[i][j]);
            }
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

    // Delete the custom functions.

    for (auto& function : functions)
        delete function.second;

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
        for (auto& child : node.getChildren())
            children.push_back(replaceDerivFunctions(child, context));
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

void ReferenceCustomDynamics::update(ContextImpl& context, int numberOfAtoms, vector<Vec3>& atomCoordinates,
                                     vector<Vec3>& velocities, vector<Vec3>& forces, vector<double>& masses,
                                     map<string, double>& globals, vector<vector<Vec3> >& perDof, bool& forcesAreValid, double tolerance) {
    if (invalidatesForces.size() == 0)
        initialize(context, masses, globals);
    int numSteps = stepType.size();
    globals.insert(context.getParameters().begin(), context.getParameters().end());
    for (auto& global : globals)
        expressionSet.setVariable(expressionSet.getVariableIndex(global.first), global.second);
    oldPos = atomCoordinates;
    map<int, double> groupEnergy;
    map<int, vector<Vec3> > groupForces;
    if (forcesAreValid)
        groupForces[context.getLastForceGroups()] = forces;
    
    // Loop over steps and execute them.
    
    for (int step = 0; step < numSteps; ) {
        int flags = forceGroupFlags[step];
        if ((needsForces[step] && groupForces.find(flags) == groupForces.end()) || (needsEnergy[step] && groupEnergy.find(flags) == groupEnergy.end())) {
            // Recompute forces and/or energy.
            
            bool computeForce = needsForces[step] || computeBothForceAndEnergy[step];
            bool computeEnergy = needsEnergy[step] || computeBothForceAndEnergy[step];
            recordChangedParameters(context, globals);
            double e = context.calcForcesAndEnergy(computeForce, computeEnergy, forceGroupFlags[step]);
            if (computeForce)
                groupForces[flags] = forces;
            if (computeEnergy) {
                groupEnergy[flags] = e;
                context.getEnergyParameterDerivatives(energyParamDerivs);
            }
            forcesAreValid = true;
        }
        
        // Execute the step.

        energy = (needsEnergy[step] ? groupEnergy[flags] : 0);
        vector<Vec3>& stepForces = (needsForces[step] ? groupForces[flags] : forces);
        int nextStep = step+1;
        bool stepInvalidatesForces = invalidatesForces[step];
        switch (stepType[step]) {
            case CustomIntegrator::ComputeGlobal: {
                uniform = SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber();
                gaussian = SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
                double result = stepExpressions[step][0].evaluate();
                globals[stepVariable[step]] = result;
                expressionSet.setVariable(stepVariableIndex[step], result);
                break;
            }
            case CustomIntegrator::ComputePerDof: {
                vector<Vec3>* results = NULL;
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
                if (stepVectorExpressions[step].size() > 0)
                    computePerParticle(numberOfAtoms, *results, atomCoordinates, velocities, stepForces, masses, perDof, globals, stepVectorExpressions[step][0]);
                else
                    computePerDof(numberOfAtoms, *results, atomCoordinates, velocities, stepForces, masses, perDof, stepExpressions[step][0]);
                break;
            }
            case CustomIntegrator::ComputeSum: {
                if (stepVectorExpressions[step].size() > 0)
                    computePerParticle(numberOfAtoms, sumBuffer, atomCoordinates, velocities, stepForces, masses, perDof, globals, stepVectorExpressions[step][0]);
                else
                    computePerDof(numberOfAtoms, sumBuffer, atomCoordinates, velocities, stepForces, masses, perDof, stepExpressions[step][0]);
                double sum = 0.0;
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
                stepInvalidatesForces = context.updateContextState();
                globals.insert(context.getParameters().begin(), context.getParameters().end());
                for (auto& global : globals)
                    expressionSet.setVariable(expressionSet.getVariableIndex(global.first), global.second);
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
        if (stepInvalidatesForces) {
            forcesAreValid = false;
            groupForces.clear();
            groupEnergy.clear();
        }
        step = nextStep;
    }
    ReferenceVirtualSites::computePositions(context.getSystem(), atomCoordinates);
    incrementTimeStep();
    recordChangedParameters(context, globals);
}

void ReferenceCustomDynamics::computePerDof(int numberOfAtoms, vector<Vec3>& results, const vector<Vec3>& atomCoordinates,
              const vector<Vec3>& velocities, const vector<Vec3>& forces, const vector<double>& masses,
              const vector<vector<Vec3> >& perDof, const CompiledExpression& expression) {
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

void ReferenceCustomDynamics::computePerParticle(int numberOfAtoms, vector<Vec3>& results, const vector<Vec3>& atomCoordinates,
              const vector<Vec3>& velocities, const vector<Vec3>& forces, const vector<double>& masses,
              const vector<vector<Vec3> >& perDof, const map<string, double>& globals, const VectorExpression& expression) {
    // Loop over all degrees of freedom.

    map<string, Vec3> variables;
    for (auto& entry : globals)
        variables[entry.first] = Vec3(entry.second, entry.second, entry.second);
    for (int i = 0; i < numberOfAtoms; i++) {
        if (masses[i] != 0.0) {
            variables["m"] = Vec3(masses[i], masses[i], masses[i]);
            variables["x"] = atomCoordinates[i];
            variables["v"] = velocities[i];
            variables["f"] = forces[i];
            variables["uniform"] = Vec3(SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber(),
                    SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber(), SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber());
            variables["gaussian"] = Vec3(SimTKOpenMMUtilities::getNormallyDistributedRandomNumber(),
                    SimTKOpenMMUtilities::getNormallyDistributedRandomNumber(), SimTKOpenMMUtilities::getNormallyDistributedRandomNumber());
            for (int j = 0; j < perDof.size(); j++)
                variables[integrator.getPerDofVariableName(j)] = perDof[j][i];
            results[i] = expression.evaluate(variables);
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
void ReferenceCustomDynamics::recordChangedParameters(OpenMM::ContextImpl& context, std::map<std::string, double>& globals) {
    for (auto& param : context.getParameters()) {
        string name = param.first;
        double value = globals[name];
        if (value != param.second)
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

double ReferenceCustomDynamics::computeKineticEnergy(OpenMM::ContextImpl& context, int numberOfAtoms, std::vector<OpenMM::Vec3>& atomCoordinates,
        std::vector<OpenMM::Vec3>& velocities, std::vector<OpenMM::Vec3>& forces, std::vector<double>& masses,
        std::map<std::string, double>& globals, std::vector<std::vector<OpenMM::Vec3> >& perDof, bool& forcesAreValid) {
    if (invalidatesForces.size() == 0)
        initialize(context, masses, globals);
    globals.insert(context.getParameters().begin(), context.getParameters().end());
    for (auto& global : globals)
        expressionSet.setVariable(expressionSet.getVariableIndex(global.first), global.second);
    if (kineticEnergyNeedsForce) {
        energy = context.calcForcesAndEnergy(true, true, -1);
        forcesAreValid = true;
    }
    computePerDof(numberOfAtoms, sumBuffer, atomCoordinates, velocities, forces, masses, perDof, kineticEnergyExpression);
    double sum = 0.0;
    for (int j = 0; j < numberOfAtoms; j++)
        if (masses[j] != 0.0)
            sum += sumBuffer[j][0]+sumBuffer[j][1]+sumBuffer[j][2];
    return sum;
}
