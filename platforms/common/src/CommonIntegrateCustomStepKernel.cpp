/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2025 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
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

#include "openmm/common/CommonIntegrateCustomStepKernel.h"
#include "openmm/common/CommonKernelUtilities.h"
#include "openmm/common/ContextSelector.h"
#include "openmm/common/ExpressionUtilities.h"
#include "openmm/Context.h"
#include "openmm/internal/ContextImpl.h"
#include "ReferenceTabulatedFunction.h"
#include "SimTKOpenMMUtilities.h"
#include "CommonKernelSources.h"
#include "lepton/CustomFunction.h"
#include "lepton/ExpressionTreeNode.h"
#include "lepton/Operation.h"
#include "lepton/Parser.h"
#include "lepton/ParsedExpression.h"

using namespace OpenMM;
using namespace std;
using namespace Lepton;

class CommonIntegrateCustomStepKernel::ReorderListener : public ComputeContext::ReorderListener {
public:
    ReorderListener(ComputeContext& cc, vector<ComputeArray>& perDofValues, vector<vector<mm_float4> >& localPerDofValuesFloat, vector<vector<mm_double4> >& localPerDofValuesDouble, vector<bool>& deviceValuesAreCurrent) :
            cc(cc), perDofValues(perDofValues), localPerDofValuesFloat(localPerDofValuesFloat), localPerDofValuesDouble(localPerDofValuesDouble), deviceValuesAreCurrent(deviceValuesAreCurrent) {
        int numAtoms = cc.getNumAtoms();
        lastAtomOrder.resize(numAtoms);
        for (int i = 0; i < numAtoms; i++)
            lastAtomOrder[i] = cc.getAtomIndex()[i];
    }
    void execute() {
        // Reorder the per-DOF variables to reflect the new atom order.

        if (perDofValues.size() == 0)
            return;
        int numAtoms = cc.getNumAtoms();
        const vector<int>& order = cc.getAtomIndex();
        for (int index = 0; index < perDofValues.size(); index++) {
            if (cc.getUseDoublePrecision() || cc.getUseMixedPrecision()) {
                if (deviceValuesAreCurrent[index])
                    perDofValues[index].download(localPerDofValuesDouble[index]);
                vector<mm_double4> swap(numAtoms);
                for (int i = 0; i < numAtoms; i++)
                    swap[lastAtomOrder[i]] = localPerDofValuesDouble[index][i];
                for (int i = 0; i < numAtoms; i++)
                    localPerDofValuesDouble[index][i] = swap[order[i]];
                perDofValues[index].upload(localPerDofValuesDouble[index]);
            }
            else {
                if (deviceValuesAreCurrent[index])
                    perDofValues[index].download(localPerDofValuesFloat[index]);
                vector<mm_float4> swap(numAtoms);
                for (int i = 0; i < numAtoms; i++)
                    swap[lastAtomOrder[i]] = localPerDofValuesFloat[index][i];
                for (int i = 0; i < numAtoms; i++)
                    localPerDofValuesFloat[index][i] = swap[order[i]];
                perDofValues[index].upload(localPerDofValuesFloat[index]);
            }
            deviceValuesAreCurrent[index] = true;
        }
        for (int i = 0; i < numAtoms; i++)
            lastAtomOrder[i] = order[i];
    }
private:
    ComputeContext& cc;
    vector<ComputeArray>& perDofValues;
    vector<vector<mm_float4> >& localPerDofValuesFloat;
    vector<vector<mm_double4> >& localPerDofValuesDouble;
    vector<bool>& deviceValuesAreCurrent;
    vector<int> lastAtomOrder;
};

class CommonIntegrateCustomStepKernel::DerivFunction : public CustomFunction {
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

void CommonIntegrateCustomStepKernel::initialize(const System& system, const CustomIntegrator& integrator) {
    cc.initializeContexts();
    ContextSelector selector(cc);
    cc.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
    numGlobalVariables = integrator.getNumGlobalVariables();
    int elementSize = (cc.getUseDoublePrecision() || cc.getUseMixedPrecision() ? sizeof(double) : sizeof(float));
    sumBuffer.initialize(cc, system.getNumParticles(), elementSize, "sumBuffer");
    summedValue.initialize(cc, 1, elementSize, "summedValue");
    perDofValues.resize(integrator.getNumPerDofVariables());
    localPerDofValuesFloat.resize(perDofValues.size());
    localPerDofValuesDouble.resize(perDofValues.size());
    for (int i = 0; i < perDofValues.size(); i++)
        perDofValues[i].initialize(cc, system.getNumParticles(), 4*elementSize, "perDofVariables");
    localValuesAreCurrent.resize(integrator.getNumPerDofVariables(), false);
    deviceValuesAreCurrent.resize(integrator.getNumPerDofVariables(), false);
    cc.addReorderListener(new ReorderListener(cc, perDofValues, localPerDofValuesFloat, localPerDofValuesDouble, deviceValuesAreCurrent));
    SimTKOpenMMUtilities::setRandomNumberSeed(integrator.getRandomNumberSeed());
}

string CommonIntegrateCustomStepKernel::createPerDofComputation(const string& variable, const Lepton::ParsedExpression& expr, CustomIntegrator& integrator,
        const string& forceName, const string& energyName, vector<const TabulatedFunction*>& functions, vector<pair<string, string> >& functionNames) {
    string tempType = (cc.getSupportsDoublePrecision() ? "double3" : "float3");
    map<string, Lepton::ParsedExpression> expressions;
    expressions[tempType+" tempResult = "] = expr;
    map<string, string> variables;
    variables["x"] = "make_"+tempType+"(position.x, position.y, position.z)";
    variables["v"] = "make_"+tempType+"(velocity.x, velocity.y, velocity.z)";
    variables[forceName] = "make_"+tempType+"(f.x, f.y, f.z)";
    variables["gaussian"] = "make_"+tempType+"(gaussian.x, gaussian.y, gaussian.z)";
    variables["uniform"] = "make_"+tempType+"(uniform.x, uniform.y, uniform.z)";
    variables["m"] = "mass";
    variables["dt"] = "stepSize";
    if (energyName != "")
        variables[energyName] = "make_"+tempType+"(energy)";
    for (int i = 0; i < integrator.getNumGlobalVariables(); i++)
        variables[integrator.getGlobalVariableName(i)] = "make_"+tempType+"(globals["+cc.intToString(globalVariableIndex[i])+"])";
    for (int i = 0; i < integrator.getNumPerDofVariables(); i++)
        variables[integrator.getPerDofVariableName(i)] = "convertToTempType3(perDof"+cc.intToString(i)+")";
    for (int i = 0; i < (int) parameterNames.size(); i++)
        variables[parameterNames[i]] = "make_"+tempType+"(globals["+cc.intToString(parameterVariableIndex[i])+"])";
    vector<pair<ExpressionTreeNode, string> > variableNodes;
    findExpressionsForDerivs(expr.getRootNode(), variableNodes);
    for (auto& var : variables)
        variableNodes.push_back(make_pair(ExpressionTreeNode(new Operation::Variable(var.first)), var.second));
    string result = cc.getExpressionUtilities().createExpressions(expressions, variableNodes, functions, functionNames, "temp", tempType);
    if (variable == "x")
        result += "position.x = tempResult.x; position.y = tempResult.y; position.z = tempResult.z;\n";
    else if (variable == "v")
        result += "velocity.x = tempResult.x; velocity.y = tempResult.y; velocity.z = tempResult.z;\n";
    else if (variable == "")
        result += "sum[index] = tempResult.x+tempResult.y+tempResult.z;\n";
    else {
        for (int i = 0; i < integrator.getNumPerDofVariables(); i++)
            if (variable == integrator.getPerDofVariableName(i)) {
                string varName = "perDof"+cc.intToString(i);
                result += varName+".x = tempResult.x; "+varName+".y = tempResult.y; "+varName+".z = tempResult.z;\n";
            }
    }
    return result;
}

void CommonIntegrateCustomStepKernel::prepareForComputation(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid) {
    ContextSelector selector(cc);
    IntegrationUtilities& integration = cc.getIntegrationUtilities();
    int numAtoms = cc.getNumAtoms();
    int numSteps = integrator.getNumComputations();
    bool useDouble = cc.getUseDoublePrecision() || cc.getUseMixedPrecision();
    string tempType = (cc.getSupportsDoublePrecision() ? "double3" : "float3");
    string perDofType = (useDouble ? "double4" : "float4");
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;

        // Initialize various data structures.

        const map<string, double>& params = context.getParameters();
        for (auto& param : params)
            parameterNames.push_back(param.first);
        kernels.resize(integrator.getNumComputations());
        requiredGaussian.resize(integrator.getNumComputations(), 0);
        requiredUniform.resize(integrator.getNumComputations(), 0);
        needsGlobals.resize(numSteps, false);
        globalExpressions.resize(numSteps);
        stepType.resize(numSteps);
        stepTarget.resize(numSteps);
        merged.resize(numSteps, false);
        modifiesParameters = false;
        sumWorkGroupSize = cc.getMaxThreadBlockSize();
        if (sumWorkGroupSize > 512)
            sumWorkGroupSize = 512;
        map<string, string> defines;
        defines["NUM_ATOMS"] = cc.intToString(cc.getNumAtoms());
        defines["PADDED_NUM_ATOMS"] = cc.intToString(cc.getPaddedNumAtoms());
        defines["WORK_GROUP_SIZE"] = cc.intToString(sumWorkGroupSize);

        // Record the tabulated functions.

        map<string, Lepton::CustomFunction*> functions;
        vector<pair<string, string> > functionNames;
        vector<const TabulatedFunction*> functionList;
        vector<string> tableTypes;
        tabulatedFunctions.resize(integrator.getNumTabulatedFunctions());
        for (int i = 0; i < integrator.getNumTabulatedFunctions(); i++) {
            functionList.push_back(&integrator.getTabulatedFunction(i));
            string name = integrator.getTabulatedFunctionName(i);
            string arrayName = "table"+cc.intToString(i);
            functionNames.push_back(make_pair(name, arrayName));
            functions[name] = createReferenceTabulatedFunction(integrator.getTabulatedFunction(i));
            int width;
            vector<float> f = cc.getExpressionUtilities().computeFunctionCoefficients(integrator.getTabulatedFunction(i), width);
            tabulatedFunctions[i].initialize<float>(cc, f.size(), "TabulatedFunction");
            tabulatedFunctions[i].upload(f);
            if (width == 1)
                tableTypes.push_back("float");
            else
                tableTypes.push_back("float"+cc.intToString(width));
        }

        // Record information about all the computation steps.

        vector<string> variable(numSteps);
        vector<int> forceGroup;
        vector<vector<Lepton::ParsedExpression> > expression;
        CustomIntegratorUtilities::analyzeComputations(context, integrator, expression, comparisons, blockEnd, invalidatesForces, needsForces, needsEnergy, computeBothForceAndEnergy, forceGroup, functions);
        for (int step = 0; step < numSteps; step++) {
            string expr;
            integrator.getComputationStep(step, stepType[step], variable[step], expr);
            if (stepType[step] == CustomIntegrator::WhileBlockStart)
                blockEnd[blockEnd[step]] = step; // Record where to branch back to.
            if (stepType[step] == CustomIntegrator::ComputeGlobal || stepType[step] == CustomIntegrator::IfBlockStart || stepType[step] == CustomIntegrator::WhileBlockStart)
                for (auto& expr : expression[step])
                    globalExpressions[step].push_back(ParsedExpression(replaceDerivFunctions(expr.getRootNode(), context)).createCompiledExpression());
        }
        for (int step = 0; step < numSteps; step++) {
            for (auto& expr : globalExpressions[step])
                expressionSet.registerExpression(expr);
        }

        // Record the indices for variables in the CompiledExpressionSet.

        gaussianVariableIndex = expressionSet.getVariableIndex("gaussian");
        uniformVariableIndex = expressionSet.getVariableIndex("uniform");
        dtVariableIndex = expressionSet.getVariableIndex("dt");
        for (int i = 0; i < integrator.getNumGlobalVariables(); i++)
            globalVariableIndex.push_back(expressionSet.getVariableIndex(integrator.getGlobalVariableName(i)));
        for (auto& name : parameterNames)
            parameterVariableIndex.push_back(expressionSet.getVariableIndex(name));

        // Record the variable names and flags for the force and energy in each step.

        forceGroupFlags.resize(numSteps, integrator.getIntegrationForceGroups());
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
        vector<string> forceName(numSteps, "f");
        vector<string> energyName(numSteps, "energy");
        stepEnergyVariableIndex.resize(numSteps, expressionSet.getVariableIndex("energy"));
        for (int step = 0; step < numSteps; step++) {
            if (needsForces[step] && forceGroup[step] > -1)
                forceName[step] = forceGroupName[forceGroup[step]];
            if (needsEnergy[step] && forceGroup[step] > -1) {
                energyName[step] = energyGroupName[forceGroup[step]];
                stepEnergyVariableIndex[step] = expressionSet.getVariableIndex(energyName[step]);
            }
            if (forceGroup[step] > -1)
                forceGroupFlags[step] = 1<<forceGroup[step];
            if (forceGroupFlags[step] == -2 && step > 0)
                forceGroupFlags[step] = forceGroupFlags[step-1];
            if (forceGroupFlags[step] != -2 && savedForces.find(forceGroupFlags[step]) == savedForces.end()) {
                savedForces[forceGroupFlags[step]] = ComputeArray();
                savedForces[forceGroupFlags[step]].initialize(cc, cc.getLongForceBuffer().getSize(), cc.getLongForceBuffer().getElementSize(), "savedForces");
            }
        }

        // Allocate space for storing global values, both on the host and the device.

        localGlobalValues.resize(expressionSet.getNumVariables());
        int elementSize = (cc.getUseDoublePrecision() || cc.getUseMixedPrecision() ? sizeof(double) : sizeof(float));
        globalValues.initialize(cc, expressionSet.getNumVariables(), elementSize, "globalValues");
        for (int i = 0; i < integrator.getNumGlobalVariables(); i++) {
            localGlobalValues[globalVariableIndex[i]] = initialGlobalVariables[i];
            expressionSet.setVariable(globalVariableIndex[i], initialGlobalVariables[i]);
        }
        for (int i = 0; i < (int) parameterVariableIndex.size(); i++) {
            double value = context.getParameter(parameterNames[i]);
            localGlobalValues[parameterVariableIndex[i]] = value;
            expressionSet.setVariable(parameterVariableIndex[i], value);
        }
        int numContextParams = context.getParameters().size();
        localPerDofEnergyParamDerivs.resize(numContextParams);
        perDofEnergyParamDerivs.initialize(cc, max(1, numContextParams), elementSize, "perDofEnergyParamDerivs");

        // Record information about the targets of steps that will be stored in global variables.

        for (int step = 0; step < numSteps; step++) {
            if (stepType[step] == CustomIntegrator::ComputeGlobal || stepType[step] == CustomIntegrator::ComputeSum) {
                if (variable[step] == "dt")
                    stepTarget[step].type = DT;
                for (int i = 0; i < integrator.getNumGlobalVariables(); i++)
                    if (variable[step] == integrator.getGlobalVariableName(i))
                        stepTarget[step].type = VARIABLE;
                for (auto& name : parameterNames)
                    if (variable[step] == name) {
                        stepTarget[step].type = PARAMETER;
                        modifiesParameters = true;
                    }
                stepTarget[step].variableIndex = expressionSet.getVariableIndex(variable[step]);
            }
        }

        // Identify which per-DOF steps are going to require global variables or context parameters.

        for (int step = 0; step < numSteps; step++) {
            if (stepType[step] == CustomIntegrator::ComputePerDof || stepType[step] == CustomIntegrator::ComputeSum) {
                for (int i = 0; i < integrator.getNumGlobalVariables(); i++)
                    if (usesVariable(expression[step][0], integrator.getGlobalVariableName(i)))
                        needsGlobals[step] = true;
                for (auto& name : parameterNames)
                    if (usesVariable(expression[step][0], name))
                        needsGlobals[step] = true;
            }
        }

        // Determine how each step will represent the position (as just a value, or a value plus a delta).

        hasAnyConstraints = (context.getSystem().getNumConstraints() > 0);
        vector<bool> storePosAsDelta(numSteps, false);
        vector<bool> loadPosAsDelta(numSteps, false);
        if (hasAnyConstraints) {
            bool beforeConstrain = false;
            for (int step = numSteps-1; step >= 0; step--) {
                if (stepType[step] == CustomIntegrator::ConstrainPositions)
                    beforeConstrain = true;
                else if (stepType[step] == CustomIntegrator::ComputePerDof && variable[step] == "x" && beforeConstrain)
                    storePosAsDelta[step] = true;
            }
            bool storedAsDelta = false;
            for (int step = 0; step < numSteps; step++) {
                loadPosAsDelta[step] = storedAsDelta;
                if (storePosAsDelta[step] == true)
                    storedAsDelta = true;
                if (stepType[step] == CustomIntegrator::ConstrainPositions)
                    storedAsDelta = false;
            }
        }

        // Identify steps that can be merged into a single kernel.

        for (int step = 1; step < numSteps; step++) {
            if (invalidatesForces[step-1] || forceGroupFlags[step] != forceGroupFlags[step-1])
                continue;
            if (stepType[step-1] == CustomIntegrator::ComputePerDof && stepType[step] == CustomIntegrator::ComputePerDof)
                merged[step] = true;
        }
        for (int step = numSteps-1; step > 0; step--)
            if (merged[step]) {
                needsForces[step-1] = (needsForces[step] || needsForces[step-1]);
                needsEnergy[step-1] = (needsEnergy[step] || needsEnergy[step-1]);
                needsGlobals[step-1] = (needsGlobals[step] || needsGlobals[step-1]);
                computeBothForceAndEnergy[step-1] = (computeBothForceAndEnergy[step] || computeBothForceAndEnergy[step-1]);
            }

        // Loop over all steps and create the kernels for them.

        for (int step = 0; step < numSteps; step++) {
            if ((stepType[step] == CustomIntegrator::ComputePerDof || stepType[step] == CustomIntegrator::ComputeSum) && !merged[step]) {
                // Compute a per-DOF value.

                stringstream compute;
                for (int i = 0; i < perDofValues.size(); i++)
                    compute << tempType<<" perDof"<<cc.intToString(i)<<" = convertToTempType3(perDofValues"<<cc.intToString(i)<<"[index]);\n";
                int numGaussian = 0, numUniform = 0;
                for (int j = step; j < numSteps && (j == step || merged[j]); j++) {
                    numGaussian += numAtoms*usesVariable(expression[j][0], "gaussian");
                    numUniform += numAtoms*usesVariable(expression[j][0], "uniform");
                    compute << "{\n";
                    if (numGaussian > 0)
                        compute << "float4 gaussian = gaussianValues[gaussianIndex+index];\n";
                    if (numUniform > 0)
                        compute << "float4 uniform = uniformValues[uniformIndex+index];\n";
                    compute << createPerDofComputation(stepType[j] == CustomIntegrator::ComputePerDof ? variable[j] : "", expression[j][0], integrator, forceName[j], energyName[j], functionList, functionNames);
                    if (variable[j] == "x") {
                        if (storePosAsDelta[j]) {
                            if (cc.getSupportsDoublePrecision())
                                compute << "posDelta[index] = convertFromDouble4(position-loadPos(posq, posqCorrection, index));\n";
                            else
                                compute << "posDelta[index] = position-posq[index];\n";
                        }
                        else
                            compute << "storePos(posq, posqCorrection, index, position);\n";
                    }
                    else if (variable[j] == "v") {
                        if (cc.getSupportsDoublePrecision())
                            compute << "velm[index] = convertFromDouble4(velocity);\n";
                        else
                            compute << "velm[index] = velocity;\n";
                    }
                    else {
                        for (int i = 0; i < perDofValues.size(); i++)
                            compute << "perDofValues"<<cc.intToString(i)<<"[index] = make_"<<perDofType<<"(perDof"<<cc.intToString(i)<<".x, perDof"<<cc.intToString(i)<<".y, perDof"<<cc.intToString(i)<<".z, 0);\n";
                    }
                    if (numGaussian > 0)
                        compute << "gaussianIndex += NUM_ATOMS;\n";
                    if (numUniform > 0)
                        compute << "uniformIndex += NUM_ATOMS;\n";
                    compute << "}\n";
                }
                map<string, string> replacements;
                replacements["COMPUTE_STEP"] = compute.str();
                stringstream args;
                for (int i = 0; i < perDofValues.size(); i++) {
                    string valueName = "perDofValues"+cc.intToString(i);
                    args << ", GLOBAL " << perDofType << "* RESTRICT " << valueName;
                }
                for (int i = 0; i < (int) tableTypes.size(); i++)
                    args << ", GLOBAL const " << tableTypes[i]<< "* RESTRICT table" << i;
                replacements["PARAMETER_ARGUMENTS"] = args.str();
                if (loadPosAsDelta[step])
                    defines["LOAD_POS_AS_DELTA"] = "1";
                else if (defines.find("LOAD_POS_AS_DELTA") != defines.end())
                    defines.erase("LOAD_POS_AS_DELTA");
                ComputeProgram program = cc.compileProgram(cc.replaceStrings(CommonKernelSources::customIntegratorPerDof, replacements), defines);
                ComputeKernel kernel = program->createKernel("computePerDof");
                kernels[step].push_back(kernel);
                requiredGaussian[step] = numGaussian;
                requiredUniform[step] = numUniform;
                kernel->addArg(cc.getPosq());
                if (cc.getUseMixedPrecision())
                    kernel->addArg(cc.getPosqCorrection());
                else
                    kernel->addArg(nullptr);
                kernel->addArg(integration.getPosDelta());
                kernel->addArg(cc.getVelm());
                kernel->addArg(cc.getLongForceBuffer());
                kernel->addArg(integration.getStepSize());
                kernel->addArg(globalValues);
                kernel->addArg(sumBuffer);
                for (int i = 0; i < 4; i++)
                    kernel->addArg();
                kernel->addArg(perDofEnergyParamDerivs);
                for (auto& array : perDofValues)
                    kernel->addArg(array);
                for (auto& array : tabulatedFunctions)
                    kernel->addArg(array);
                if (stepType[step] == CustomIntegrator::ComputeSum) {
                    // Create a second kernel for this step that sums the values.

                    program = cc.compileProgram(CommonKernelSources::customIntegrator, defines);
                    kernel = program->createKernel(useDouble ? "computeDoubleSum" : "computeFloatSum");
                    kernels[step].push_back(kernel);
                    kernel->addArg(sumBuffer);
                    kernel->addArg(summedValue);
                    kernel->addArg(numAtoms);
                }
            }
            else if (stepType[step] == CustomIntegrator::ConstrainPositions) {
                // Apply position constraints.

                ComputeProgram program = cc.compileProgram(CommonKernelSources::customIntegrator, defines);
                ComputeKernel kernel = program->createKernel("applyPositionDeltas");
                kernels[step].push_back(kernel);
                kernel->addArg(cc.getPosq());
                if (cc.getUseMixedPrecision())
                    kernel->addArg(cc.getPosqCorrection());
                else
                    kernel->addArg(nullptr);
                kernel->addArg(integration.getPosDelta());
            }
        }

        // Initialize the random number generator.

        int maxUniformRandoms = 1;
        for (int required : requiredUniform)
            maxUniformRandoms = max(maxUniformRandoms, required);
        uniformRandoms.initialize<mm_float4>(cc, maxUniformRandoms, "uniformRandoms");
        randomSeed.initialize<mm_int4>(cc, cc.getNumThreadBlocks()*64, "randomSeed");
        vector<mm_int4> seed(randomSeed.getSize());
        int rseed = integrator.getRandomNumberSeed();
        // A random seed of 0 means use a unique one
        if (rseed == 0)
            rseed = osrngseed();
        unsigned int r = (unsigned int) (rseed+1);
        for (auto& s : seed) {
            s.x = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
            s.y = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
            s.z = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
            s.w = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
        }
        randomSeed.upload(seed);
        ComputeProgram randomProgram = cc.compileProgram(CommonKernelSources::customIntegrator, defines);
        randomKernel = randomProgram->createKernel("generateRandomNumbers");
        randomKernel->addArg(maxUniformRandoms);
        randomKernel->addArg(uniformRandoms);
        randomKernel->addArg(randomSeed);

        // Create the kernel for computing kinetic energy.

        stringstream computeKE;
        for (int i = 0; i < perDofValues.size(); i++)
            computeKE << tempType<<" perDof"<<cc.intToString(i)<<" = convertToTempType3(perDofValues"<<cc.intToString(i)<<"[index]);\n";
        Lepton::ParsedExpression keExpression = Lepton::Parser::parse(integrator.getKineticEnergyExpression()).optimize();
        computeKE << createPerDofComputation("", keExpression, integrator, "f", "", functionList, functionNames);
        map<string, string> replacements;
        replacements["COMPUTE_STEP"] = computeKE.str();
        stringstream args;
        for (int i = 0; i < perDofValues.size(); i++) {
            string valueName = "perDofValues"+cc.intToString(i);
            args << ", GLOBAL " << perDofType << "* RESTRICT " << valueName;
        }
        for (int i = 0; i < (int) tableTypes.size(); i++)
            args << ", GLOBAL const " << tableTypes[i]<< "* RESTRICT table" << i;
        replacements["PARAMETER_ARGUMENTS"] = args.str();
        if (defines.find("LOAD_POS_AS_DELTA") != defines.end())
            defines.erase("LOAD_POS_AS_DELTA");
        ComputeProgram program = cc.compileProgram(cc.replaceStrings(CommonKernelSources::customIntegratorPerDof, replacements), defines);
        kineticEnergyKernel = program->createKernel("computePerDof");
        kineticEnergyKernel->addArg(cc.getPosq());
        if (cc.getUseMixedPrecision())
            kineticEnergyKernel->addArg(cc.getPosqCorrection());
        else
            kineticEnergyKernel->addArg(nullptr);
        kineticEnergyKernel->addArg(integration.getPosDelta());
        kineticEnergyKernel->addArg(cc.getVelm());
        kineticEnergyKernel->addArg(cc.getLongForceBuffer());
        kineticEnergyKernel->addArg(integration.getStepSize());
        kineticEnergyKernel->addArg(globalValues);
        kineticEnergyKernel->addArg(sumBuffer);
        kineticEnergyKernel->addArg();
        kineticEnergyKernel->addArg();
        kineticEnergyKernel->addArg(uniformRandoms);
        if (cc.getUseDoublePrecision() || cc.getUseMixedPrecision())
            kineticEnergyKernel->addArg(0.0);
        else
            kineticEnergyKernel->addArg(0.0f);
        kineticEnergyKernel->addArg(perDofEnergyParamDerivs);
        for (auto& array : perDofValues)
            kineticEnergyKernel->addArg(array);
        for (auto& array : tabulatedFunctions)
            kineticEnergyKernel->addArg(array);

        // Create a second kernel to sum the values.

        program = cc.compileProgram(CommonKernelSources::customIntegrator, defines);
        sumKineticEnergyKernel = program->createKernel(useDouble ? "computeDoubleSum" : "computeFloatSum");
        sumKineticEnergyKernel->addArg(sumBuffer);
        sumKineticEnergyKernel->addArg(summedValue);
        sumKineticEnergyKernel->addArg(numAtoms);

        // Delete the custom functions.

        for (auto& function : functions)
            delete function.second;
    }

    // Make sure all values (variables, parameters, etc.) are up to date.

    for (int i = 0; i < perDofValues.size(); i++) {
        if (!deviceValuesAreCurrent[i]) {
            if (useDouble)
                perDofValues[i].upload(localPerDofValuesDouble[i]);
            else
                perDofValues[i].upload(localPerDofValuesFloat[i]);
            deviceValuesAreCurrent[i] = true;
        }
        localValuesAreCurrent[i] = false;
    }
    double stepSize = integrator.getStepSize();
    recordGlobalValue(stepSize, GlobalTarget(DT, dtVariableIndex), integrator);
    for (int i = 0; i < (int) parameterNames.size(); i++) {
        double value = context.getParameter(parameterNames[i]);
        if (value != localGlobalValues[parameterVariableIndex[i]]) {
            expressionSet.setVariable(parameterVariableIndex[i], value);
            localGlobalValues[parameterVariableIndex[i]] = value;
            deviceGlobalsAreCurrent = false;
        }
    }
}

ExpressionTreeNode CommonIntegrateCustomStepKernel::replaceDerivFunctions(const ExpressionTreeNode& node, ContextImpl& context) {
    // This is called recursively to identify calls to the deriv() function inside global expressions,
    // and replace them with a custom function that returns the correct value.

    const Operation& op = node.getOperation();
    if (op.getId() == Operation::CUSTOM && op.getName() == "deriv") {
        string param = node.getChildren()[1].getOperation().getName();
        if (context.getParameters().find(param) == context.getParameters().end())
            throw OpenMMException("The second argument to deriv() must be a context parameter");
        needsEnergyParamDerivs = true;
        return ExpressionTreeNode(new Operation::Custom("deriv", new DerivFunction(energyParamDerivs, param)));
    }
    else {
        vector<ExpressionTreeNode> children;
        for (auto& child : node.getChildren())
            children.push_back(replaceDerivFunctions(child, context));
        return ExpressionTreeNode(op.clone(), children);
    }
}

void CommonIntegrateCustomStepKernel::findExpressionsForDerivs(const ExpressionTreeNode& node, vector<pair<ExpressionTreeNode, string> >& variableNodes) {
    // This is called recursively to identify calls to the deriv() function inside per-DOF expressions,
    // and record the code to replace them with.

    const Operation& op = node.getOperation();
    if (op.getId() == Operation::CUSTOM && op.getName() == "deriv") {
        string param = node.getChildren()[1].getOperation().getName();
        int index;
        for (index = 0; index < perDofEnergyParamDerivNames.size() && param != perDofEnergyParamDerivNames[index]; index++)
            ;
        if (index == perDofEnergyParamDerivNames.size())
            perDofEnergyParamDerivNames.push_back(param);
        string tempType = (cc.getSupportsDoublePrecision() ? "double3" : "float3");
        variableNodes.push_back(make_pair(node, "make_"+tempType+"(energyParamDerivs["+cc.intToString(index)+"])"));
        needsEnergyParamDerivs = true;
    }
    else {
        for (auto& child : node.getChildren())
            findExpressionsForDerivs(child, variableNodes);
    }
}

void CommonIntegrateCustomStepKernel::execute(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid) {
    ContextSelector selector(cc);
    prepareForComputation(context, integrator, forcesAreValid);
    IntegrationUtilities& integration = cc.getIntegrationUtilities();
    int numAtoms = cc.getNumAtoms();
    int numSteps = integrator.getNumComputations();
    if (!forcesAreValid)
        savedEnergy.clear();

    // Loop over computation steps in the integrator and execute them.

    for (int step = 0; step < numSteps; ) {
        int nextStep = step+1;
        int forceGroups = forceGroupFlags[step];
        int lastForceGroups = context.getLastForceGroups();
        bool haveForces = (!needsForces[step] || (forcesAreValid && lastForceGroups == forceGroups));
        bool haveEnergy = (!needsEnergy[step] || savedEnergy.find(forceGroups) != savedEnergy.end());
        if (!haveForces || !haveEnergy) {
            if (forcesAreValid) {
                if (savedForces.find(lastForceGroups) != savedForces.end() && validSavedForces.find(lastForceGroups) == validSavedForces.end()) {
                    // The forces are still valid.  We just need a different force group right now.  Save the old
                    // forces in case we need them again.

                    cc.getLongForceBuffer().copyTo(savedForces[lastForceGroups]);
                    validSavedForces.insert(lastForceGroups);
                }
            }
            else
                validSavedForces.clear();

            // Recompute forces and/or energy.  Figure out what is actually needed
            // between now and the next time they get invalidated again.

            bool computeForce = (needsForces[step] || computeBothForceAndEnergy[step]);
            bool computeEnergy = (needsEnergy[step] || computeBothForceAndEnergy[step]);
            if (!computeEnergy && validSavedForces.find(forceGroups) != validSavedForces.end()) {
                // We can just restore the forces we saved earlier.

                savedForces[forceGroups].copyTo(cc.getLongForceBuffer());
                context.getLastForceGroups() = forceGroups;
            }
            else {
                recordChangedParameters(context);
                energy = context.calcForcesAndEnergy(computeForce, computeEnergy, forceGroups);
                savedEnergy[forceGroups] = energy;
                if (needsEnergyParamDerivs) {
                    context.getEnergyParameterDerivatives(energyParamDerivs);
                    if (perDofEnergyParamDerivNames.size() > 0) {
                        for (int i = 0; i < perDofEnergyParamDerivNames.size(); i++)
                            localPerDofEnergyParamDerivs[i] = energyParamDerivs[perDofEnergyParamDerivNames[i]];
                        perDofEnergyParamDerivs.upload(localPerDofEnergyParamDerivs, true);
                    }
                }
            }
            forcesAreValid = true;
        }
        if (needsEnergy[step])
            energy = savedEnergy[forceGroups];
        if (needsGlobals[step] && !deviceGlobalsAreCurrent) {
            // Upload the global values to the device.

            globalValues.upload(localGlobalValues, true);
            deviceGlobalsAreCurrent = true;
        }
        bool stepInvalidatesForces = invalidatesForces[step];
        if (stepType[step] == CustomIntegrator::ComputePerDof && !merged[step]) {
            kernels[step][0]->setArg(9, integration.prepareRandomNumbers(requiredGaussian[step]));
            kernels[step][0]->setArg(8, integration.getRandom());
            kernels[step][0]->setArg(10, uniformRandoms);
            if (cc.getUseDoublePrecision() || cc.getUseMixedPrecision())
                kernels[step][0]->setArg(11, energy);
            else
                kernels[step][0]->setArg(11, (float) energy);
            if (requiredUniform[step] > 0)
                randomKernel->execute(numAtoms, 64);
            kernels[step][0]->execute(numAtoms, 128);
        }
        else if (stepType[step] == CustomIntegrator::ComputeGlobal) {
            expressionSet.setVariable(uniformVariableIndex, SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber());
            expressionSet.setVariable(gaussianVariableIndex, SimTKOpenMMUtilities::getNormallyDistributedRandomNumber());
            expressionSet.setVariable(stepEnergyVariableIndex[step], energy);
            recordGlobalValue(globalExpressions[step][0].evaluate(), stepTarget[step], integrator);
        }
        else if (stepType[step] == CustomIntegrator::ComputeSum) {
            kernels[step][0]->setArg(9, integration.prepareRandomNumbers(requiredGaussian[step]));
            kernels[step][0]->setArg(8, integration.getRandom());
            kernels[step][0]->setArg(10, uniformRandoms);
            if (cc.getUseDoublePrecision() || cc.getUseMixedPrecision())
                kernels[step][0]->setArg(11, energy);
            else
                kernels[step][0]->setArg(11, (float) energy);
            if (requiredUniform[step] > 0)
                randomKernel->execute(numAtoms, 64);
            cc.clearBuffer(sumBuffer);
            kernels[step][0]->execute(numAtoms, 128);
            kernels[step][1]->execute(sumWorkGroupSize, sumWorkGroupSize);
            if (cc.getUseDoublePrecision() || cc.getUseMixedPrecision()) {
                double value;
                summedValue.download(&value);
                recordGlobalValue(value, stepTarget[step], integrator);
            }
            else {
                float value;
                summedValue.download(&value);
                recordGlobalValue(value, stepTarget[step], integrator);
            }
        }
        else if (stepType[step] == CustomIntegrator::UpdateContextState) {
            recordChangedParameters(context);
            stepInvalidatesForces = context.updateContextState();
        }
        else if (stepType[step] == CustomIntegrator::ConstrainPositions) {
            if (hasAnyConstraints) {
                cc.getIntegrationUtilities().applyConstraints(integrator.getConstraintTolerance());
                kernels[step][0]->execute(numAtoms);
            }
            cc.getIntegrationUtilities().computeVirtualSites();
        }
        else if (stepType[step] == CustomIntegrator::ConstrainVelocities) {
            cc.getIntegrationUtilities().applyVelocityConstraints(integrator.getConstraintTolerance());
        }
        else if (stepType[step] == CustomIntegrator::IfBlockStart) {
            if (!evaluateCondition(step))
                nextStep = blockEnd[step]+1;
        }
        else if (stepType[step] == CustomIntegrator::WhileBlockStart) {
            if (!evaluateCondition(step))
                nextStep = blockEnd[step]+1;
        }
        else if (stepType[step] == CustomIntegrator::BlockEnd) {
            if (blockEnd[step] != -1)
                nextStep = blockEnd[step]; // Return to the start of a while block.
        }
        if (stepInvalidatesForces) {
            forcesAreValid = false;
            savedEnergy.clear();
        }
        step = nextStep;
    }
    recordChangedParameters(context);

    // Update the time and step count.

    cc.setTime(cc.getTime()+integrator.getStepSize());
    cc.setStepCount(cc.getStepCount()+1);
    cc.reorderAtoms();
    if (cc.getAtomsWereReordered()) {
        forcesAreValid = false;
        validSavedForces.clear();
    }

    // Reduce UI lag.

    flushPeriodically(cc);
}

bool CommonIntegrateCustomStepKernel::evaluateCondition(int step) {
    expressionSet.setVariable(uniformVariableIndex, SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber());
    expressionSet.setVariable(gaussianVariableIndex, SimTKOpenMMUtilities::getNormallyDistributedRandomNumber());
    expressionSet.setVariable(stepEnergyVariableIndex[step], energy);
    double lhs = globalExpressions[step][0].evaluate();
    double rhs = globalExpressions[step][1].evaluate();
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
    throw OpenMMException("Invalid comparison operator");
}

double CommonIntegrateCustomStepKernel::computeKineticEnergy(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid) {
    ContextSelector selector(cc);
    prepareForComputation(context, integrator, forcesAreValid);
    cc.clearBuffer(sumBuffer);
    kineticEnergyKernel->setArg(8, cc.getIntegrationUtilities().getRandom());
    kineticEnergyKernel->setArg(9, 0);
    kineticEnergyKernel->execute(cc.getNumAtoms());
    sumKineticEnergyKernel->execute(sumWorkGroupSize, sumWorkGroupSize);
    if (cc.getUseDoublePrecision() || cc.getUseMixedPrecision()) {
        double ke;
        summedValue.download(&ke);
        return ke;
    }
    else {
        float ke;
        summedValue.download(&ke);
        return ke;
    }
}

void CommonIntegrateCustomStepKernel::recordGlobalValue(double value, GlobalTarget target, CustomIntegrator& integrator) {
    switch (target.type) {
        case DT:
            if (value != localGlobalValues[dtVariableIndex])
                deviceGlobalsAreCurrent = false;
            expressionSet.setVariable(dtVariableIndex, value);
            localGlobalValues[dtVariableIndex] = value;
            cc.getIntegrationUtilities().setNextStepSize(value);
            integrator.setStepSize(value);
            break;
        case VARIABLE:
        case PARAMETER:
            expressionSet.setVariable(target.variableIndex, value);
            localGlobalValues[target.variableIndex] = value;
            deviceGlobalsAreCurrent = false;
            break;
    }
}

void CommonIntegrateCustomStepKernel::recordChangedParameters(ContextImpl& context) {
    if (!modifiesParameters)
        return;
    for (int i = 0; i < (int) parameterNames.size(); i++) {
        double value = context.getParameter(parameterNames[i]);
        if (value != localGlobalValues[parameterVariableIndex[i]])
            context.setParameter(parameterNames[i], localGlobalValues[parameterVariableIndex[i]]);
    }
}

void CommonIntegrateCustomStepKernel::getGlobalVariables(ContextImpl& context, vector<double>& values) const {
    if (!globalValues.isInitialized()) {
        // The data structures haven't been created yet, so just return the list of values that was given earlier.

        values = initialGlobalVariables;
        return;
    }
    values.resize(numGlobalVariables);
    for (int i = 0; i < numGlobalVariables; i++)
        values[i] = localGlobalValues[globalVariableIndex[i]];
}

void CommonIntegrateCustomStepKernel::setGlobalVariables(ContextImpl& context, const vector<double>& values) {
    if (numGlobalVariables == 0)
        return;
    if (!globalValues.isInitialized()) {
        // The data structures haven't been created yet, so just store the list of values.

        initialGlobalVariables = values;
        return;
    }
    for (int i = 0; i < numGlobalVariables; i++) {
        localGlobalValues[globalVariableIndex[i]] = values[i];
        expressionSet.setVariable(globalVariableIndex[i], values[i]);
    }
    deviceGlobalsAreCurrent = false;
}

void CommonIntegrateCustomStepKernel::getPerDofVariable(ContextImpl& context, int variable, vector<Vec3>& values) const {
    ContextSelector selector(cc);
    values.resize(perDofValues[variable].getSize());
    const vector<int>& order = cc.getAtomIndex();
    if (cc.getUseDoublePrecision() || cc.getUseMixedPrecision()) {
        if (!localValuesAreCurrent[variable]) {
            perDofValues[variable].download(localPerDofValuesDouble[variable]);
            localValuesAreCurrent[variable] = true;
        }
        for (int i = 0; i < (int) values.size(); i++) {
            values[order[i]][0] = localPerDofValuesDouble[variable][i].x;
            values[order[i]][1] = localPerDofValuesDouble[variable][i].y;
            values[order[i]][2] = localPerDofValuesDouble[variable][i].z;
        }
    }
    else {
        if (!localValuesAreCurrent[variable]) {
            perDofValues[variable].download(localPerDofValuesFloat[variable]);
            localValuesAreCurrent[variable] = true;
        }
        for (int i = 0; i < (int) values.size(); i++) {
            values[order[i]][0] = localPerDofValuesFloat[variable][i].x;
            values[order[i]][1] = localPerDofValuesFloat[variable][i].y;
            values[order[i]][2] = localPerDofValuesFloat[variable][i].z;
        }
    }
}

void CommonIntegrateCustomStepKernel::setPerDofVariable(ContextImpl& context, int variable, const vector<Vec3>& values) {
    const vector<int>& order = cc.getAtomIndex();
    localValuesAreCurrent[variable] = true;
    deviceValuesAreCurrent[variable] = false;
    if (cc.getUseDoublePrecision() || cc.getUseMixedPrecision()) {
        localPerDofValuesDouble[variable].resize(values.size());
        for (int i = 0; i < (int) values.size(); i++)
            localPerDofValuesDouble[variable][i] = mm_double4(values[order[i]][0], values[order[i]][1], values[order[i]][2], 0);
    }
    else {
        localPerDofValuesFloat[variable].resize(values.size());
        for (int i = 0; i < (int) values.size(); i++)
            localPerDofValuesFloat[variable][i] = mm_float4(values[order[i]][0], values[order[i]][1], values[order[i]][2], 0);
    }
}
