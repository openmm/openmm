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

#include "openmm/common/CommonCalcCustomGBForceKernel.h"
#include "openmm/common/CommonKernelUtilities.h"
#include "openmm/common/ContextSelector.h"
#include "openmm/common/ExpressionUtilities.h"
#include "openmm/Context.h"
#include "openmm/internal/ContextImpl.h"
#include "CommonKernelSources.h"
#include "lepton/CustomFunction.h"
#include "lepton/ExpressionTreeNode.h"
#include "lepton/Operation.h"
#include "lepton/Parser.h"
#include "lepton/ParsedExpression.h"

using namespace OpenMM;
using namespace std;
using namespace Lepton;

class CommonCalcCustomGBForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const CustomGBForce& force) : force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        thread_local static vector<double> params1, params2;
        force.getParticleParameters(particle1, params1);
        force.getParticleParameters(particle2, params2);
        for (int i = 0; i < (int) params1.size(); i++)
            if (params1[i] != params2[i])
                return false;
        return true;
    }
    int getNumParticleGroups() {
        return force.getNumExclusions();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int particle1, particle2;
        force.getExclusionParticles(index, particle1, particle2);
        particles.resize(2);
        particles[0] = particle1;
        particles[1] = particle2;
    }
    bool areGroupsIdentical(int group1, int group2) {
        return true;
    }
private:
    const CustomGBForce& force;
};

CommonCalcCustomGBForceKernel::~CommonCalcCustomGBForceKernel() {
    ContextSelector selector(cc);
    if (params != NULL)
        delete params;
    if (computedValues != NULL)
        delete computedValues;
    if (energyDerivs != NULL)
        delete energyDerivs;
    if (energyDerivChain != NULL)
        delete energyDerivChain;
    for (auto d : dValuedParam)
        delete d;
}

void CommonCalcCustomGBForceKernel::initialize(const System& system, const CustomGBForce& force) {
    ContextSelector selector(cc);
    if (cc.getNumContexts() > 1)
        throw OpenMMException("CustomGBForce does not support using multiple devices");
    NonbondedUtilities& nb = cc.getNonbondedUtilities();
    cutoff = force.getCutoffDistance();
    bool useExclusionsForValue = false;
    numComputedValues = force.getNumComputedValues();
    vector<string> computedValueNames(numComputedValues);
    vector<string> computedValueExpressions(numComputedValues);
    if (numComputedValues > 0) {
        CustomGBForce::ComputationType type;
        force.getComputedValueParameters(0, computedValueNames[0], computedValueExpressions[0], type);
        if (type == CustomGBForce::SingleParticle)
            throw OpenMMException("The first computed value for a CustomGBForce must be of type ParticlePair or ParticlePairNoExclusions.");
        useExclusionsForValue = (type == CustomGBForce::ParticlePair);
        for (int i = 1; i < numComputedValues; i++) {
            force.getComputedValueParameters(i, computedValueNames[i], computedValueExpressions[i], type);
            if (type != CustomGBForce::SingleParticle)
                throw OpenMMException("A CustomGBForce may only have one computed value of type ParticlePair or ParticlePairNoExclusions.");
        }
    }
    int forceIndex;
    for (forceIndex = 0; forceIndex < system.getNumForces() && &system.getForce(forceIndex) != &force; ++forceIndex)
        ;
    string prefix = "custom"+cc.intToString(forceIndex)+"_";

    // Record parameters and exclusions.

    int numParticles = force.getNumParticles();
    int paddedNumParticles = cc.getPaddedNumAtoms();
    int numParams = force.getNumPerParticleParameters();
    params = new ComputeParameterSet(cc, force.getNumPerParticleParameters(), paddedNumParticles, "customGBParameters", true);
    computedValues = new ComputeParameterSet(cc, numComputedValues, paddedNumParticles, "customGBComputedValues", true, cc.getUseDoublePrecision());
    if (force.getNumGlobalParameters() > 0)
        globals.initialize<float>(cc, force.getNumGlobalParameters(), "customGBGlobals");
    vector<vector<float> > paramVector(paddedNumParticles, vector<float>(numParams, 0));
    vector<vector<int> > exclusionList(numParticles);
    for (int i = 0; i < numParticles; i++) {
        vector<double> parameters;
        force.getParticleParameters(i, parameters);
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
        exclusionList[i].push_back(i);
    }
    for (int i = 0; i < force.getNumExclusions(); i++) {
        int particle1, particle2;
        force.getExclusionParticles(i, particle1, particle2);
        exclusionList[particle1].push_back(particle2);
        exclusionList[particle2].push_back(particle1);
    }
    params->setParameterValues(paramVector);

    // Record the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    vector<pair<string, string> > functionDefinitions;
    vector<const TabulatedFunction*> functionList;
    stringstream tableArgs;
    tabulatedFunctionArrays.resize(force.getNumTabulatedFunctions());
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        functionList.push_back(&force.getTabulatedFunction(i));
        string name = force.getTabulatedFunctionName(i);
        tabulatedFunctionUpdateCount[name] = force.getTabulatedFunction(i).getUpdateCount();
        string arrayName = prefix+"table"+cc.intToString(i);
        functionDefinitions.push_back(make_pair(name, arrayName));
        functions[name] = cc.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cc.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
        tabulatedFunctionArrays[i].initialize<float>(cc, f.size(), "TabulatedFunction");
        tabulatedFunctionArrays[i].upload(f);
        nb.addArgument(ComputeParameterInfo(tabulatedFunctionArrays[i], arrayName, "float", width));
        tableArgs << ", GLOBAL const float";
        if (width > 1)
            tableArgs << width;
        tableArgs << "* RESTRICT " << arrayName;
    }

    // Record the global parameters.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
    }
    if (globals.isInitialized())
        globals.upload(globalParamValues);

    // Record derivatives of expressions needed for the chain rule terms.

    vector<vector<Lepton::ParsedExpression> > valueGradientExpressions(numComputedValues);
    vector<vector<Lepton::ParsedExpression> > valueDerivExpressions(numComputedValues);
    vector<vector<Lepton::ParsedExpression> > valueParamDerivExpressions(numComputedValues);
    needParameterGradient = false;
    for (int i = 0; i < numComputedValues; i++) {
        Lepton::ParsedExpression ex = Lepton::Parser::parse(computedValueExpressions[i], functions).optimize();
        if (i > 0) {
            valueGradientExpressions[i].push_back(ex.differentiate("x").optimize());
            valueGradientExpressions[i].push_back(ex.differentiate("y").optimize());
            valueGradientExpressions[i].push_back(ex.differentiate("z").optimize());
            if (!isZeroExpression(valueGradientExpressions[i][0]) || !isZeroExpression(valueGradientExpressions[i][1]) || !isZeroExpression(valueGradientExpressions[i][2]))
                needParameterGradient = true;
            for (int j = 0; j < i; j++)
                valueDerivExpressions[i].push_back(ex.differentiate(computedValueNames[j]).optimize());
        }
        for (int j = 0; j < force.getNumEnergyParameterDerivatives(); j++)
            valueParamDerivExpressions[i].push_back(ex.differentiate(force.getEnergyParameterDerivativeName(j)).optimize());
    }
    vector<vector<Lepton::ParsedExpression> > energyDerivExpressions(force.getNumEnergyTerms());
    vector<vector<Lepton::ParsedExpression> > energyParamDerivExpressions(force.getNumEnergyTerms());
    vector<bool> needChainForValue(numComputedValues, false);
    for (int i = 0; i < force.getNumEnergyTerms(); i++) {
        string expression;
        CustomGBForce::ComputationType type;
        force.getEnergyTermParameters(i, expression, type);
        Lepton::ParsedExpression ex = Lepton::Parser::parse(expression, functions).optimize();
        for (int j = 0; j < numComputedValues; j++) {
            if (type == CustomGBForce::SingleParticle) {
                energyDerivExpressions[i].push_back(ex.differentiate(computedValueNames[j]).optimize());
                if (!isZeroExpression(energyDerivExpressions[i].back()))
                    needChainForValue[j] = true;
            }
            else {
                energyDerivExpressions[i].push_back(ex.differentiate(computedValueNames[j]+"1").optimize());
                if (!isZeroExpression(energyDerivExpressions[i].back()))
                    needChainForValue[j] = true;
                energyDerivExpressions[i].push_back(ex.differentiate(computedValueNames[j]+"2").optimize());
                if (!isZeroExpression(energyDerivExpressions[i].back()))
                    needChainForValue[j] = true;
            }
        }
        for (int j = 0; j < force.getNumEnergyParameterDerivatives(); j++)
            energyParamDerivExpressions[i].push_back(ex.differentiate(force.getEnergyParameterDerivativeName(j)).optimize());
    }
    bool deviceIsCpu = cc.getIsCPU();
    int elementSize = (cc.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
    valueBuffers.initialize<long long>(cc, cc.getPaddedNumAtoms(), "customGBValueBuffers");
    longEnergyDerivs.initialize<long long>(cc, numComputedValues*cc.getPaddedNumAtoms(), "customGBLongEnergyDerivatives");
    energyDerivs = new ComputeParameterSet(cc, numComputedValues, cc.getPaddedNumAtoms(), "customGBEnergyDerivatives", true);
    cc.addAutoclearBuffer(valueBuffers);
    energyDerivChain = new ComputeParameterSet(cc, numComputedValues, cc.getPaddedNumAtoms(), "customGBEnergyDerivativeChain", true);
    needEnergyParamDerivs = (force.getNumEnergyParameterDerivatives() > 0);
    dValue0dParam.resize(force.getNumEnergyParameterDerivatives());
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        dValuedParam.push_back(new ComputeParameterSet(cc, numComputedValues, cc.getPaddedNumAtoms(), "dValuedParam", true, cc.getUseDoublePrecision()));
        dValue0dParam[i].initialize<long long>(cc, cc.getPaddedNumAtoms(), "dValue0dParam");
        cc.addAutoclearBuffer(dValue0dParam[i]);
        string name = force.getEnergyParameterDerivativeName(i);
        cc.addEnergyParameterDerivative(name);
    }

    // Create the kernels.

    bool useCutoff = (force.getNonbondedMethod() != CustomGBForce::NoCutoff);
    bool usePeriodic = (force.getNonbondedMethod() != CustomGBForce::NoCutoff && force.getNonbondedMethod() != CustomGBForce::CutoffNonPeriodic);
    int numAtomBlocks = cc.getPaddedNumAtoms()/32;
    {
        // Create the N2 value kernel.

        vector<pair<ExpressionTreeNode, string> > variables;
        map<string, string> rename;
        ExpressionTreeNode rnode(new Operation::Variable("r"));
        variables.push_back(make_pair(rnode, "r"));
        variables.push_back(make_pair(ExpressionTreeNode(new Operation::Square(), rnode), "r2"));
        variables.push_back(make_pair(ExpressionTreeNode(new Operation::Reciprocal(), rnode), "invR"));
        for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
            const string& name = force.getPerParticleParameterName(i);
            variables.push_back(makeVariable(name+"1", "((real) params"+params->getParameterSuffix(i, "1)")));
            variables.push_back(makeVariable(name+"2", "((real) params"+params->getParameterSuffix(i, "2)")));
            rename[name+"1"] = name+"2";
            rename[name+"2"] = name+"1";
        }
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = "globals["+cc.intToString(i)+"]";
            variables.push_back(makeVariable(name, value));
        }
        map<string, Lepton::ParsedExpression> n2ValueExpressions;
        stringstream n2ValueSource;
        Lepton::ParsedExpression ex = Lepton::Parser::parse(computedValueExpressions[0], functions).optimize();
        n2ValueExpressions["tempValue1 = "] = ex;
        n2ValueExpressions["tempValue2 = "] = ex.renameVariables(rename);
        for (int i = 0; i < valueParamDerivExpressions[0].size(); i++) {
            string variableBase = "temp_dValue0dParam"+cc.intToString(i+1);
            if (!isZeroExpression(valueParamDerivExpressions[0][i])) {
                n2ValueExpressions[variableBase+"_1 = "] = valueParamDerivExpressions[0][i];
                n2ValueExpressions[variableBase+"_2 = "] = valueParamDerivExpressions[0][i].renameVariables(rename);
            }
        }
        n2ValueSource << cc.getExpressionUtilities().createExpressions(n2ValueExpressions, variables, functionList, functionDefinitions, "temp");
        map<string, string> replacements;
        string n2ValueStr = n2ValueSource.str();
        replacements["COMPUTE_VALUE"] = n2ValueStr;
        stringstream extraArgs, atomParams, loadLocal1, loadLocal2, load1, load2, tempDerivs1, tempDerivs2, storeDeriv1, storeDeriv2;
        if (force.getNumGlobalParameters() > 0)
            extraArgs << ", GLOBAL const float* globals";
        pairValueUsesParam.resize(params->getParameterInfos().size(), false);
        for (int i = 0; i < (int) params->getParameterInfos().size(); i++) {
            ComputeParameterInfo& buffer = params->getParameterInfos()[i];
            string paramName = "params"+cc.intToString(i+1);
            if (n2ValueStr.find(paramName+"1") != n2ValueStr.npos || n2ValueStr.find(paramName+"2") != n2ValueStr.npos) {
                extraArgs << ", GLOBAL const " << buffer.getType() << "* RESTRICT global_" << paramName;
                atomParams << "LOCAL " << buffer.getType() << " local_" << paramName << "[LOCAL_BUFFER_SIZE];\n";
                loadLocal1 << "local_" << paramName << "[localAtomIndex] = " << paramName << "1;\n";
                loadLocal2 << "local_" << paramName << "[localAtomIndex] = global_" << paramName << "[j];\n";
                load1 << buffer.getType() << " " << paramName << "1 = global_" << paramName << "[atom1];\n";
                load2 << buffer.getType() << " " << paramName << "2 = local_" << paramName << "[atom2];\n";
                pairValueUsesParam[i] = true;
            }
        }
        for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
            string derivName = "dValue0dParam"+cc.intToString(i+1);
            extraArgs << ", GLOBAL mm_ulong* RESTRICT global_" << derivName;
            atomParams << "LOCAL real local_" << derivName << "[LOCAL_BUFFER_SIZE];\n";
            loadLocal2 << "local_" << derivName << "[localAtomIndex] = 0;\n";
            load1 << "real " << derivName << " = 0;\n";
            if (!isZeroExpression(valueParamDerivExpressions[0][i])) {
                load2 << "real temp_" << derivName << "_1 = 0;\n";
                load2 << "real temp_" << derivName << "_2 = 0;\n";
                tempDerivs1 << derivName << " += temp_" << derivName << "_1;\n";
                if (deviceIsCpu)
                    tempDerivs2 << "local_" << derivName << "[j] += temp_" << derivName << "_2;\n";
                else
                    tempDerivs2 << "local_" << derivName << "[tbx+tj] += temp_" << derivName << "_2;\n";
                storeDeriv1 << "ATOMIC_ADD(&global_" << derivName << "[offset1], (mm_ulong) realToFixedPoint(" << derivName << "));\n";
                if (deviceIsCpu)
                    storeDeriv2 << "ATOMIC_ADD(&global_" << derivName << "[offset2], (mm_ulong) realToFixedPoint(local_" << derivName << "[tgx]));\n";
                else
                    storeDeriv2 << "ATOMIC_ADD(&global_" << derivName << "[offset2], (mm_ulong) realToFixedPoint(local_" << derivName << "[LOCAL_ID]));\n";
            }
        }
        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
        replacements["ATOM_PARAMETER_DATA"] = atomParams.str();
        replacements["LOAD_LOCAL_PARAMETERS_FROM_1"] = loadLocal1.str();
        replacements["LOAD_LOCAL_PARAMETERS_FROM_GLOBAL"] = loadLocal2.str();
        replacements["LOAD_ATOM1_PARAMETERS"] = load1.str();
        replacements["LOAD_ATOM2_PARAMETERS"] = load2.str();
        replacements["ADD_TEMP_DERIVS1"] = tempDerivs1.str();
        replacements["ADD_TEMP_DERIVS2"] = tempDerivs2.str();
        replacements["STORE_PARAM_DERIVS1"] = storeDeriv1.str();
        replacements["STORE_PARAM_DERIVS2"] = storeDeriv2.str();
        if (useCutoff)
            pairValueDefines["USE_CUTOFF"] = "1";
        if (usePeriodic)
            pairValueDefines["USE_PERIODIC"] = "1";
        if (useExclusionsForValue)
            pairValueDefines["USE_EXCLUSIONS"] = "1";
        pairValueDefines["LOCAL_BUFFER_SIZE"] = cc.intToString(deviceIsCpu ? 32 : nb.getForceThreadBlockSize());
        pairValueDefines["CUTOFF_SQUARED"] = cc.doubleToString(cutoff*cutoff);
        pairValueDefines["NUM_ATOMS"] = cc.intToString(cc.getNumAtoms());
        pairValueDefines["PADDED_NUM_ATOMS"] = cc.intToString(cc.getPaddedNumAtoms());
        pairValueDefines["NUM_BLOCKS"] = cc.intToString(numAtomBlocks);
        pairValueDefines["TILE_SIZE"] = "32";
        string file;
        if (deviceIsCpu)
            file = CommonKernelSources::customGBValueN2_cpu;
        else
            file = CommonKernelSources::customGBValueN2;
        pairValueSrc = cc.replaceStrings(file, replacements);
    }
    {
        // Create the kernel to reduce the N2 value and calculate other values.

        stringstream reductionSource, extraArgs, deriv0;
        if (force.getNumGlobalParameters() > 0)
            extraArgs << ", GLOBAL const float* globals";
        for (int i = 0; i < (int) params->getParameterInfos().size(); i++) {
            ComputeParameterInfo& buffer = params->getParameterInfos()[i];
            string paramName = "params"+cc.intToString(i+1);
            extraArgs << ", GLOBAL const " << buffer.getType() << "* RESTRICT " << paramName;
        }
        for (int i = 0; i < (int) computedValues->getParameterInfos().size(); i++) {
            ComputeParameterInfo& buffer = computedValues->getParameterInfos()[i];
            string valueName = "values"+cc.intToString(i+1);
            extraArgs << ", GLOBAL " << buffer.getType() << "* RESTRICT global_" << valueName;
            reductionSource << buffer.getType() << " local_" << valueName << ";\n";
        }
        for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
            string variableName = "dValuedParam_0_"+cc.intToString(i);
            extraArgs << ", GLOBAL const mm_long* RESTRICT dValue0dParam" << i;
            deriv0 << "real " << variableName << " = RECIP((real) 0x100000000)*dValue0dParam" << i << "[index];\n";
            for (int j = 0; j < dValuedParam[i]->getParameterInfos().size(); j++)
                extraArgs << ", GLOBAL real* RESTRICT global_dValuedParam_" << j << "_" << i;
            deriv0 << "global_dValuedParam_0_" << i << "[index] = dValuedParam_0_" << i << ";\n";
        }
        reductionSource << "local_values" << computedValues->getParameterSuffix(0) << " = sum;\n";
        map<string, string> variables;
        variables["x"] = "pos.x";
        variables["y"] = "pos.y";
        variables["z"] = "pos.z";
        for (int i = 0; i < force.getNumPerParticleParameters(); i++)
            variables[force.getPerParticleParameterName(i)] = "params"+params->getParameterSuffix(i, "[index]");
        for (int i = 0; i < force.getNumGlobalParameters(); i++)
            variables[force.getGlobalParameterName(i)] = "globals["+cc.intToString(i)+"]";
        for (int i = 1; i < numComputedValues; i++) {
            variables[computedValueNames[i-1]] = "local_values"+computedValues->getParameterSuffix(i-1);
            map<string, Lepton::ParsedExpression> valueExpressions;
            valueExpressions["local_values"+computedValues->getParameterSuffix(i)+" = "] = Lepton::Parser::parse(computedValueExpressions[i], functions).optimize();
            reductionSource << cc.getExpressionUtilities().createExpressions(valueExpressions, variables, functionList, functionDefinitions, "value"+cc.intToString(i)+"_temp");
        }
        for (int i = 0; i < (int) computedValues->getParameterInfos().size(); i++) {
            string valueName = "values"+cc.intToString(i+1);
            reductionSource << "global_" << valueName << "[index] = local_" << valueName << ";\n";
        }
        if (needEnergyParamDerivs) {
            map<string, Lepton::ParsedExpression> derivExpressions;
            for (int i = 1; i < numComputedValues; i++) {
                for (int j = 0; j < valueParamDerivExpressions[i].size(); j++)
                    derivExpressions["real dValuedParam_"+cc.intToString(i)+"_"+cc.intToString(j)+" = "] = valueParamDerivExpressions[i][j];
                for (int j = 0; j < i; j++)
                    derivExpressions["real dVdV_"+cc.intToString(i)+"_"+cc.intToString(j)+" = "] = valueDerivExpressions[i][j];
            }
            reductionSource << cc.getExpressionUtilities().createExpressions(derivExpressions, variables, functionList, functionDefinitions, "derivChain_temp");
            for (int i = 1; i < numComputedValues; i++) {
                for (int j = 0; j < i; j++)
                    for (int k = 0; k < valueParamDerivExpressions[i].size(); k++)
                        reductionSource << "dValuedParam_" << i << "_" << k << " += dVdV_" << i << "_" << j << "*dValuedParam_" << j <<"_" << k << ";\n";
                for (int j = 0; j < valueParamDerivExpressions[i].size(); j++)
                    reductionSource << "global_dValuedParam_" << i << "_" << j << "[index] = dValuedParam_" << i << "_" << j << ";\n";
            }
        }
        map<string, string> replacements;
        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
        replacements["REDUCE_PARAM0_DERIV"] = deriv0.str();
        replacements["COMPUTE_VALUES"] = reductionSource.str();
        map<string, string> defines;
        defines["NUM_ATOMS"] = cc.intToString(cc.getNumAtoms());
        ComputeProgram program = cc.compileProgram(cc.replaceStrings(CommonKernelSources::customGBValuePerParticle, replacements), defines);
        perParticleValueKernel = program->createKernel("computePerParticleValues");
    }
    {
        // Create the N2 energy kernel.

        vector<pair<ExpressionTreeNode, string> > variables;
        ExpressionTreeNode rnode(new Operation::Variable("r"));
        variables.push_back(make_pair(rnode, "r"));
        variables.push_back(make_pair(ExpressionTreeNode(new Operation::Square(), rnode), "r2"));
        variables.push_back(make_pair(ExpressionTreeNode(new Operation::Reciprocal(), rnode), "invR"));
        for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
            const string& name = force.getPerParticleParameterName(i);
            variables.push_back(makeVariable(name+"1", "((real) params"+params->getParameterSuffix(i, "1)")));
            variables.push_back(makeVariable(name+"2", "((real) params"+params->getParameterSuffix(i, "2)")));
        }
        for (int i = 0; i < numComputedValues; i++) {
            variables.push_back(makeVariable(computedValueNames[i]+"1", "values"+computedValues->getParameterSuffix(i, "1")));
            variables.push_back(makeVariable(computedValueNames[i]+"2", "values"+computedValues->getParameterSuffix(i, "2")));
        }
        for (int i = 0; i < force.getNumGlobalParameters(); i++)
            variables.push_back(makeVariable(force.getGlobalParameterName(i), "globals["+cc.intToString(i)+"]"));
        stringstream n2EnergySource;
        bool anyExclusions = (force.getNumExclusions() > 0);
        for (int i = 0; i < force.getNumEnergyTerms(); i++) {
            string expression;
            CustomGBForce::ComputationType type;
            force.getEnergyTermParameters(i, expression, type);
            if (type == CustomGBForce::SingleParticle)
                continue;
            bool exclude = (anyExclusions && type == CustomGBForce::ParticlePair);
            map<string, Lepton::ParsedExpression> n2EnergyExpressions;
            n2EnergyExpressions["tempEnergy += "] = Lepton::Parser::parse(expression, functions).optimize();
            n2EnergyExpressions["dEdR += "] = Lepton::Parser::parse(expression, functions).differentiate("r").optimize();
            for (int j = 0; j < numComputedValues; j++) {
                if (needChainForValue[j]) {
                    string index = cc.intToString(j+1);
                    n2EnergyExpressions["/*"+cc.intToString(i+1)+"*/ deriv"+index+"_1 += "] = energyDerivExpressions[i][2*j];
                    n2EnergyExpressions["/*"+cc.intToString(i+1)+"*/ deriv"+index+"_2 += "] = energyDerivExpressions[i][2*j+1];
                }
            }
            for (int j = 0; j < force.getNumEnergyParameterDerivatives(); j++)
                n2EnergyExpressions["energyParamDeriv"+cc.intToString(j)+" += interactionScale*"] = energyParamDerivExpressions[i][j];
            if (exclude)
                n2EnergySource << "if (!isExcluded) {\n";
            n2EnergySource << cc.getExpressionUtilities().createExpressions(n2EnergyExpressions, variables, functionList, functionDefinitions, "temp");
            if (exclude)
                n2EnergySource << "}\n";
        }
        map<string, string> replacements;
        string n2EnergyStr = n2EnergySource.str();
        replacements["COMPUTE_INTERACTION"] = n2EnergyStr;
        stringstream extraArgs, atomParams, loadLocal1, loadLocal2, clearLocal, load1, load2, declare1, recordDeriv, storeDerivs1, storeDerivs2, initParamDerivs, saveParamDerivs;
        if (force.getNumGlobalParameters() > 0)
            extraArgs << ", GLOBAL const float* globals";
        pairEnergyUsesParam.resize(params->getParameterInfos().size(), false);
        for (int i = 0; i < (int) params->getParameterInfos().size(); i++) {
            ComputeParameterInfo& buffer = params->getParameterInfos()[i];
            string paramName = "params"+cc.intToString(i+1);
            if (n2EnergyStr.find(paramName+"1") != n2EnergyStr.npos || n2EnergyStr.find(paramName+"2") != n2EnergyStr.npos) {
                extraArgs << ", GLOBAL const " << buffer.getType() << "* RESTRICT global_" << paramName;
                atomParams << "LOCAL " << buffer.getType() << " local_" << paramName << "[LOCAL_BUFFER_SIZE];\n";
                loadLocal1 << "local_" << paramName << "[localAtomIndex] = " << paramName << "1;\n";
                loadLocal2 << "local_" << paramName << "[localAtomIndex] = global_" << paramName << "[j];\n";
                load1 << buffer.getType() << " " << paramName << "1 = global_" << paramName << "[atom1];\n";
                load2 << buffer.getType() << " " << paramName << "2 = local_" << paramName << "[atom2];\n";
                pairEnergyUsesParam[i] = true;
            }
        }
        pairEnergyUsesValue.resize(computedValues->getParameterInfos().size(), false);
        for (int i = 0; i < (int) computedValues->getParameterInfos().size(); i++) {
            ComputeParameterInfo& buffer = computedValues->getParameterInfos()[i];
            string valueName = "values"+cc.intToString(i+1);
            if (n2EnergyStr.find(valueName+"1") != n2EnergyStr.npos || n2EnergyStr.find(valueName+"2") != n2EnergyStr.npos) {
                extraArgs << ", GLOBAL const " << buffer.getType() << "* RESTRICT global_" << valueName;
                atomParams << "LOCAL " << buffer.getType() << " local_" << valueName << "[LOCAL_BUFFER_SIZE];\n";
                loadLocal1 << "local_" << valueName << "[localAtomIndex] = " << valueName << "1;\n";
                loadLocal2 << "local_" << valueName << "[localAtomIndex] = global_" << valueName << "[j];\n";
                load1 << buffer.getType() << " " << valueName << "1 = global_" << valueName << "[atom1];\n";
                load2 << buffer.getType() << " " << valueName << "2 = local_" << valueName << "[atom2];\n";
                pairEnergyUsesValue[i] = true;
            }
        }
        extraArgs << ", GLOBAL mm_ulong* RESTRICT derivBuffers";
        for (int i = 0; i < numComputedValues; i++) {
            string index = cc.intToString(i+1);
            atomParams << "LOCAL real local_deriv" << index << "[LOCAL_BUFFER_SIZE];\n";
            clearLocal << "local_deriv" << index << "[localAtomIndex] = 0.0f;\n";
            declare1 << "real deriv" << index << "_1 = 0;\n";
            load2 << "real deriv" << index << "_2 = 0;\n";
            recordDeriv << "local_deriv" << index << "[atom2] += deriv" << index << "_2;\n";
            storeDerivs1 << "STORE_DERIVATIVE_1(" << index << ")\n";
            storeDerivs2 << "STORE_DERIVATIVE_2(" << index << ")\n";
        }
        if (needEnergyParamDerivs) {
            extraArgs << ", GLOBAL mixed* RESTRICT energyParamDerivs";
            const vector<string>& allParamDerivNames = cc.getEnergyParamDerivNames();
            int numDerivs = allParamDerivNames.size();
            for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
                initParamDerivs << "mixed energyParamDeriv" << i << " = 0;\n";
                for (int index = 0; index < numDerivs; index++)
                    if (allParamDerivNames[index] == force.getEnergyParameterDerivativeName(i))
                        saveParamDerivs << "energyParamDerivs[GLOBAL_ID*" << numDerivs << "+" << index << "] += energyParamDeriv" << i << ";\n";
            }
        }
        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
        replacements["ATOM_PARAMETER_DATA"] = atomParams.str();
        replacements["LOAD_LOCAL_PARAMETERS_FROM_1"] = loadLocal1.str();
        replacements["LOAD_LOCAL_PARAMETERS_FROM_GLOBAL"] = loadLocal2.str();
        replacements["CLEAR_LOCAL_DERIVATIVES"] = clearLocal.str();
        replacements["LOAD_ATOM1_PARAMETERS"] = load1.str();
        replacements["LOAD_ATOM2_PARAMETERS"] = load2.str();
        replacements["DECLARE_ATOM1_DERIVATIVES"] = declare1.str();
        replacements["RECORD_DERIVATIVE_2"] = recordDeriv.str();
        replacements["STORE_DERIVATIVES_1"] = storeDerivs1.str();
        replacements["STORE_DERIVATIVES_2"] = storeDerivs2.str();
        replacements["INIT_PARAM_DERIVS"] = initParamDerivs.str();
        replacements["SAVE_PARAM_DERIVS"] = saveParamDerivs.str();
        if (useCutoff)
            pairEnergyDefines["USE_CUTOFF"] = "1";
        if (usePeriodic)
            pairEnergyDefines["USE_PERIODIC"] = "1";
        if (anyExclusions)
            pairEnergyDefines["USE_EXCLUSIONS"] = "1";
        pairEnergyDefines["LOCAL_BUFFER_SIZE"] = cc.intToString(deviceIsCpu ? 32 : nb.getForceThreadBlockSize());
        pairEnergyDefines["CUTOFF_SQUARED"] = cc.doubleToString(cutoff*cutoff);
        pairEnergyDefines["NUM_ATOMS"] = cc.intToString(cc.getNumAtoms());
        pairEnergyDefines["PADDED_NUM_ATOMS"] = cc.intToString(cc.getPaddedNumAtoms());
        pairEnergyDefines["NUM_BLOCKS"] = cc.intToString(numAtomBlocks);
        pairEnergyDefines["TILE_SIZE"] = "32";
        string file;
        if (deviceIsCpu)
            file = CommonKernelSources::customGBEnergyN2_cpu;
        else
            file = CommonKernelSources::customGBEnergyN2;
        pairEnergySrc = cc.replaceStrings(file, replacements);
    }
    {
        // Create the kernel to reduce the derivatives and calculate per-particle energy terms.

        stringstream compute, extraArgs, reduce, initParamDerivs, saveParamDerivs;
        if (force.getNumGlobalParameters() > 0)
            extraArgs << ", GLOBAL const float* globals";
        for (int i = 0; i < (int) params->getParameterInfos().size(); i++) {
            ComputeParameterInfo& buffer = params->getParameterInfos()[i];
            string paramName = "params"+cc.intToString(i+1);
            extraArgs << ", GLOBAL const " << buffer.getType() << "* RESTRICT " << paramName;
        }
        for (int i = 0; i < (int) computedValues->getParameterInfos().size(); i++) {
            ComputeParameterInfo& buffer = computedValues->getParameterInfos()[i];
            string valueName = "values"+cc.intToString(i+1);
            extraArgs << ", GLOBAL const " << buffer.getType() << "* RESTRICT " << valueName;
        }
        for (int i = 0; i < (int) energyDerivs->getParameterInfos().size(); i++) {
            ComputeParameterInfo& buffer = energyDerivs->getParameterInfos()[i];
            string index = cc.intToString(i+1);
            extraArgs << ", GLOBAL " << buffer.getType() << "* RESTRICT derivBuffers" << index;
            compute << buffer.getType() << " deriv" << index << " = derivBuffers" << index << "[index];\n";
        }
        for (int i = 0; i < (int) energyDerivChain->getParameterInfos().size(); i++) {
            ComputeParameterInfo& buffer = energyDerivChain->getParameterInfos()[i];
            string index = cc.intToString(i+1);
            extraArgs << ", GLOBAL " << buffer.getType() << "* RESTRICT derivChain" << index;
        }
        extraArgs << ", GLOBAL const mm_long* RESTRICT derivBuffersIn";
        for (int i = 0; i < energyDerivs->getNumParameters(); ++i)
            reduce << "derivBuffers" << energyDerivs->getParameterSuffix(i, "[index]") <<
                    " = RECIP((real) 0x100000000)*derivBuffersIn[index+PADDED_NUM_ATOMS*" << cc.intToString(i) << "];\n";
        if (needEnergyParamDerivs) {
            extraArgs << ", GLOBAL mixed* RESTRICT energyParamDerivs";
            const vector<string>& allParamDerivNames = cc.getEnergyParamDerivNames();
            int numDerivs = allParamDerivNames.size();
            for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
                initParamDerivs << "mixed energyParamDeriv" << i << " = 0;\n";
                for (int index = 0; index < numDerivs; index++)
                    if (allParamDerivNames[index] == force.getEnergyParameterDerivativeName(i))
                        saveParamDerivs << "energyParamDerivs[GLOBAL_ID*" << numDerivs << "+" << index << "] += energyParamDeriv" << i << ";\n";
            }
        }

        // Compute the various expressions.

        map<string, string> variables;
        variables["x"] = "pos.x";
        variables["y"] = "pos.y";
        variables["z"] = "pos.z";
        for (int i = 0; i < force.getNumPerParticleParameters(); i++)
            variables[force.getPerParticleParameterName(i)] = "params"+params->getParameterSuffix(i, "[index]");
        for (int i = 0; i < force.getNumGlobalParameters(); i++)
            variables[force.getGlobalParameterName(i)] = "globals["+cc.intToString(i)+"]";
        for (int i = 0; i < numComputedValues; i++)
            variables[computedValueNames[i]] = "values"+computedValues->getParameterSuffix(i, "[index]");
        map<string, Lepton::ParsedExpression> expressions;
        for (int i = 0; i < force.getNumEnergyTerms(); i++) {
            string expression;
            CustomGBForce::ComputationType type;
            force.getEnergyTermParameters(i, expression, type);
            if (type != CustomGBForce::SingleParticle)
                continue;
            Lepton::ParsedExpression parsed = Lepton::Parser::parse(expression, functions).optimize();
            expressions["/*"+cc.intToString(i+1)+"*/ energy += "] = parsed;
            for (int j = 0; j < numComputedValues; j++)
                expressions["/*"+cc.intToString(i+1)+"*/ deriv"+energyDerivs->getParameterSuffix(j)+" += "] = energyDerivExpressions[i][j];
            Lepton::ParsedExpression gradx = parsed.differentiate("x").optimize();
            Lepton::ParsedExpression grady = parsed.differentiate("y").optimize();
            Lepton::ParsedExpression gradz = parsed.differentiate("z").optimize();
            if (!isZeroExpression(gradx))
                expressions["/*"+cc.intToString(i+1)+"*/ force.x -= "] = gradx;
            if (!isZeroExpression(grady))
                expressions["/*"+cc.intToString(i+1)+"*/ force.y -= "] = grady;
            if (!isZeroExpression(gradz))
                expressions["/*"+cc.intToString(i+1)+"*/ force.z -= "] = gradz;
            for (int j = 0; j < force.getNumEnergyParameterDerivatives(); j++)
                expressions["/*"+cc.intToString(i+1)+"*/ energyParamDeriv"+cc.intToString(j)+" += "] = energyParamDerivExpressions[i][j];
        }
        for (int i = 1; i < numComputedValues; i++)
            for (int j = 0; j < i; j++)
                expressions["real dV"+cc.intToString(i)+"dV"+cc.intToString(j)+" = "] = valueDerivExpressions[i][j];
        compute << cc.getExpressionUtilities().createExpressions(expressions, variables, functionList, functionDefinitions, "temp");

        // Record values.

        for (int i = 0; i < (int) energyDerivs->getParameterInfos().size(); i++) {
            string index = cc.intToString(i+1);
            compute << "derivBuffers" << index << "[index] = deriv" << index << ";\n";
        }
        compute << "forceBuffers[index] += realToFixedPoint(force.x);\n";
        compute << "forceBuffers[index+PADDED_NUM_ATOMS] += realToFixedPoint(force.y);\n";
        compute << "forceBuffers[index+PADDED_NUM_ATOMS*2] += realToFixedPoint(force.z);\n";
        for (int i = 1; i < numComputedValues; i++) {
            compute << "real totalDeriv"<<i<<" = dV"<<i<<"dV0";
            for (int j = 1; j < i; j++)
                compute << " + totalDeriv"<<j<<"*dV"<<i<<"dV"<<j;
            compute << ";\n";
            compute << "deriv"<<(i+1)<<" *= totalDeriv"<<i<<";\n";
        }
        for (int i = 0; i < (int) energyDerivs->getParameterInfos().size(); i++) {
            string index = cc.intToString(i+1);
            compute << "derivChain" << index << "[index] = deriv" << index << ";\n";
        }
        map<string, string> replacements;
        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
        replacements["REDUCE_DERIVATIVES"] = reduce.str();
        replacements["COMPUTE_ENERGY"] = compute.str();
        replacements["INIT_PARAM_DERIVS"] = initParamDerivs.str();
        replacements["SAVE_PARAM_DERIVS"] = saveParamDerivs.str();
        map<string, string> defines;
        defines["NUM_ATOMS"] = cc.intToString(cc.getNumAtoms());
        defines["PADDED_NUM_ATOMS"] = cc.intToString(cc.getPaddedNumAtoms());
        ComputeProgram program = cc.compileProgram(cc.replaceStrings(CommonKernelSources::customGBEnergyPerParticle, replacements), defines);
        perParticleEnergyKernel = program->createKernel("computePerParticleEnergy");
    }
    if (needParameterGradient || needEnergyParamDerivs) {
        // Create the kernel to compute chain rule terms for computed values that depend explicitly on particle coordinates, and for
        // derivatives with respect to global parameters.

        stringstream compute, extraArgs, initParamDerivs, saveParamDerivs;
        if (force.getNumGlobalParameters() > 0)
            extraArgs << ", GLOBAL const float* globals";
        for (int i = 0; i < (int) params->getParameterInfos().size(); i++) {
            ComputeParameterInfo& buffer = params->getParameterInfos()[i];
            string paramName = "params"+cc.intToString(i+1);
            extraArgs << ", GLOBAL const " << buffer.getType() << "* RESTRICT " << paramName;
        }
        for (int i = 0; i < (int) computedValues->getParameterInfos().size(); i++) {
            ComputeParameterInfo& buffer = computedValues->getParameterInfos()[i];
            string valueName = "values"+cc.intToString(i+1);
            extraArgs << ", GLOBAL const " << buffer.getType() << "* RESTRICT " << valueName;
        }
        for (int i = 0; i < (int) energyDerivs->getParameterInfos().size(); i++) {
            ComputeParameterInfo& buffer = energyDerivs->getParameterInfos()[i];
            string index = cc.intToString(i+1);
            extraArgs << ", GLOBAL " << buffer.getType() << "* RESTRICT derivBuffers" << index;
            compute << buffer.getType() << " deriv" << index << " = derivBuffers" << index << "[index];\n";
        }
        if (needEnergyParamDerivs) {
            extraArgs << ", GLOBAL mixed* RESTRICT energyParamDerivs";
            const vector<string>& allParamDerivNames = cc.getEnergyParamDerivNames();
            int numDerivs = allParamDerivNames.size();
            for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
                for (int j = 0; j < dValuedParam[i]->getParameterInfos().size(); j++)
                    extraArgs << ", GLOBAL real* RESTRICT dValuedParam_" << j << "_" << i;
                initParamDerivs << "mixed energyParamDeriv" << i << " = 0;\n";
                for (int index = 0; index < numDerivs; index++)
                    if (allParamDerivNames[index] == force.getEnergyParameterDerivativeName(i))
                        saveParamDerivs << "energyParamDerivs[GLOBAL_ID*" << numDerivs << "+" << index << "] += energyParamDeriv" << i << ";\n";
            }
        }
        map<string, string> variables;
        variables["x"] = "pos.x";
        variables["y"] = "pos.y";
        variables["z"] = "pos.z";
        for (int i = 0; i < force.getNumPerParticleParameters(); i++)
            variables[force.getPerParticleParameterName(i)] = "params"+params->getParameterSuffix(i, "[index]");
        for (int i = 0; i < force.getNumGlobalParameters(); i++)
            variables[force.getGlobalParameterName(i)] = "globals["+cc.intToString(i)+"]";
        for (int i = 0; i < numComputedValues; i++)
            variables[computedValueNames[i]] = "values"+computedValues->getParameterSuffix(i, "[index]");
        if (needParameterGradient) {
            for (int i = 1; i < numComputedValues; i++) {
                string is = cc.intToString(i);
                compute << "real3 dV"<<is<<"dR = make_real3(0);\n";
                for (int j = 1; j < i; j++) {
                    if (!isZeroExpression(valueDerivExpressions[i][j])) {
                        map<string, Lepton::ParsedExpression> derivExpressions;
                        string js = cc.intToString(j);
                        derivExpressions["real dV"+is+"dV"+js+" = "] = valueDerivExpressions[i][j];
                        compute << cc.getExpressionUtilities().createExpressions(derivExpressions, variables, functionList, functionDefinitions, "temp_"+is+"_"+js);
                        compute << "dV"<<is<<"dR += dV"<<is<<"dV"<<js<<"*dV"<<js<<"dR;\n";
                    }
                }
                map<string, Lepton::ParsedExpression> gradientExpressions;
                if (!isZeroExpression(valueGradientExpressions[i][0]))
                    gradientExpressions["dV"+is+"dR.x += "] = valueGradientExpressions[i][0];
                if (!isZeroExpression(valueGradientExpressions[i][1]))
                    gradientExpressions["dV"+is+"dR.y += "] = valueGradientExpressions[i][1];
                if (!isZeroExpression(valueGradientExpressions[i][2]))
                    gradientExpressions["dV"+is+"dR.z += "] = valueGradientExpressions[i][2];
                compute << cc.getExpressionUtilities().createExpressions(gradientExpressions, variables, functionList, functionDefinitions, "gradtemp_"+is);
            }
            for (int i = 1; i < numComputedValues; i++)
                compute << "force -= deriv"<<energyDerivs->getParameterSuffix(i)<<"*dV"<<i<<"dR;\n";
        }
        if (needEnergyParamDerivs)
            for (int i = 0; i < numComputedValues; i++)
                for (int j = 0; j < dValuedParam.size(); j++)
                    compute << "energyParamDeriv"<<j<<" += deriv"<<energyDerivs->getParameterSuffix(i)<<"*dValuedParam_"<<i<<"_"<<j<<"[index];\n";
        map<string, string> replacements;
        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
        replacements["COMPUTE_FORCES"] = compute.str();
        replacements["INIT_PARAM_DERIVS"] = initParamDerivs.str();
        replacements["SAVE_PARAM_DERIVS"] = saveParamDerivs.str();
        map<string, string> defines;
        defines["NUM_ATOMS"] = cc.intToString(cc.getNumAtoms());
        defines["PADDED_NUM_ATOMS"] = cc.intToString(cc.getPaddedNumAtoms());
        ComputeProgram program = cc.compileProgram(cc.replaceStrings(CommonKernelSources::customGBGradientChainRule, replacements), defines);
        gradientChainRuleKernel = program->createKernel("computeGradientChainRuleTerms");
    }
    {
        // Create the code to calculate chain rule terms as part of the default nonbonded kernel.

        vector<pair<ExpressionTreeNode, string> > globalVariables;
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = "globals["+cc.intToString(i)+"]";
            globalVariables.push_back(makeVariable(name, prefix+value));
        }
        vector<pair<ExpressionTreeNode, string> > variables = globalVariables;
        map<string, string> rename;
        ExpressionTreeNode rnode(new Operation::Variable("r"));
        variables.push_back(make_pair(rnode, "r"));
        variables.push_back(make_pair(ExpressionTreeNode(new Operation::Square(), rnode), "r2"));
        variables.push_back(make_pair(ExpressionTreeNode(new Operation::Reciprocal(), rnode), "invR"));
        for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
            const string& name = force.getPerParticleParameterName(i);
            variables.push_back(makeVariable(name+"1", "((real) "+prefix+"params"+params->getParameterSuffix(i, "1)")));
            variables.push_back(makeVariable(name+"2", "((real) "+prefix+"params"+params->getParameterSuffix(i, "2)")));
            rename[name+"1"] = name+"2";
            rename[name+"2"] = name+"1";
        }
        map<string, Lepton::ParsedExpression> derivExpressions;
        stringstream chainSource;
        Lepton::ParsedExpression dVdR = Lepton::Parser::parse(computedValueExpressions[0], functions).differentiate("r").optimize();
        derivExpressions["real dV0dR1 = "] = dVdR;
        derivExpressions["real dV0dR2 = "] = dVdR.renameVariables(rename);
        chainSource << cc.getExpressionUtilities().createExpressions(derivExpressions, variables, functionList, functionDefinitions, prefix+"temp0_");
        if (needChainForValue[0]) {
            if (useExclusionsForValue)
                chainSource << "if (!isExcluded) {\n";
            chainSource << "tempForce -= dV0dR1*" << prefix << "dEdV" << energyDerivs->getParameterSuffix(0, "1") << ";\n";
            chainSource << "tempForce -= dV0dR2*" << prefix << "dEdV" << energyDerivs->getParameterSuffix(0, "2") << ";\n";
            if (useExclusionsForValue)
                chainSource << "}\n";
        }
        for (int i = 1; i < numComputedValues; i++) {
            if (needChainForValue[i]) {
                chainSource << "tempForce -= dV0dR1*" << prefix << "dEdV" << energyDerivs->getParameterSuffix(i, "1") << ";\n";
                chainSource << "tempForce -= dV0dR2*" << prefix << "dEdV" << energyDerivs->getParameterSuffix(i, "2") << ";\n";
            }
        }
        map<string, string> replacements;
        string chainStr = chainSource.str();
        replacements["COMPUTE_FORCE"] = chainStr;
        string source = cc.replaceStrings(CommonKernelSources::customGBChainRule, replacements);
        vector<ComputeParameterInfo> parameters;
        vector<ComputeParameterInfo> arguments;
        for (int i = 0; i < (int) params->getParameterInfos().size(); i++) {
            ComputeParameterInfo& buffer = params->getParameterInfos()[i];
            string paramName = prefix+"params"+cc.intToString(i+1);
            if (chainStr.find(paramName+"1") != chainStr.npos || chainStr.find(paramName+"2") != chainStr.npos)
                parameters.push_back(ComputeParameterInfo(buffer.getArray(), paramName, buffer.getComponentType(), buffer.getNumComponents()));
        }
        for (int i = 0; i < (int) computedValues->getParameterInfos().size(); i++) {
            ComputeParameterInfo& buffer = computedValues->getParameterInfos()[i];
            string paramName = prefix+"values"+cc.intToString(i+1);
            if (chainStr.find(paramName+"1") != chainStr.npos || chainStr.find(paramName+"2") != chainStr.npos)
                parameters.push_back(ComputeParameterInfo(buffer.getArray(), paramName, buffer.getComponentType(), buffer.getNumComponents()));
        }
        for (int i = 0; i < (int) energyDerivChain->getParameterInfos().size(); i++) {
            if (needChainForValue[i]) {
                ComputeParameterInfo& buffer = energyDerivChain->getParameterInfos()[i];
                string paramName = prefix+"dEdV"+cc.intToString(i+1);
                parameters.push_back(ComputeParameterInfo(buffer.getArray(), paramName, buffer.getComponentType(), buffer.getNumComponents()));
            }
        }
        if (globals.isInitialized()) {
            globals.upload(globalParamValues);
            arguments.push_back(ComputeParameterInfo(globals, prefix+"globals", "float", 1));
        }
        nb.addInteraction(useCutoff, usePeriodic, force.getNumExclusions() > 0, cutoff, exclusionList, source, force.getForceGroup());
        for (auto param : parameters)
            nb.addParameter(param);
        for (auto arg : arguments)
            nb.addArgument(arg);
    }
    info = new ForceInfo(force);
    cc.addForce(info);
    cc.addAutoclearBuffer(longEnergyDerivs);
}

double CommonCalcCustomGBForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    ContextSelector selector(cc);
    bool deviceIsCpu = cc.getIsCPU();
    NonbondedUtilities& nb = cc.getNonbondedUtilities();
    int elementSize = (cc.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;

        // These two kernels can't be compiled in initialize(), because the nonbonded utilities object
        // has not yet been initialized then.

        {
            int numExclusionTiles = nb.getExclusionTiles().getSize();
            pairValueDefines["NUM_TILES_WITH_EXCLUSIONS"] = cc.intToString(numExclusionTiles);
            int numContexts = cc.getNumContexts();
            int startExclusionIndex = cc.getContextIndex()*numExclusionTiles/numContexts;
            int endExclusionIndex = (cc.getContextIndex()+1)*numExclusionTiles/numContexts;
            pairValueDefines["FIRST_EXCLUSION_TILE"] = cc.intToString(startExclusionIndex);
            pairValueDefines["LAST_EXCLUSION_TILE"] = cc.intToString(endExclusionIndex);
            pairValueDefines["CUTOFF"] = cc.doubleToString(cutoff);
            ComputeProgram program = cc.compileProgram(pairValueSrc, pairValueDefines);
            pairValueKernel = program->createKernel("computeN2Value");
            pairValueSrc = "";
            pairValueDefines.clear();
        }
        {
            int numExclusionTiles = nb.getExclusionTiles().getSize();
            pairEnergyDefines["NUM_TILES_WITH_EXCLUSIONS"] = cc.intToString(numExclusionTiles);
            int numContexts = cc.getNumContexts();
            int startExclusionIndex = cc.getContextIndex()*numExclusionTiles/numContexts;
            int endExclusionIndex = (cc.getContextIndex()+1)*numExclusionTiles/numContexts;
            pairEnergyDefines["FIRST_EXCLUSION_TILE"] = cc.intToString(startExclusionIndex);
            pairEnergyDefines["LAST_EXCLUSION_TILE"] = cc.intToString(endExclusionIndex);
            pairEnergyDefines["CUTOFF"] = cc.doubleToString(cutoff);
            ComputeProgram program = cc.compileProgram(pairEnergySrc, pairEnergyDefines);
            pairEnergyKernel = program->createKernel("computeN2Energy");
            pairEnergySrc = "";
            pairEnergyDefines.clear();
        }

        // Set arguments for kernels.

        maxTiles = (nb.getUseCutoff() ? nb.getInteractingTiles().getSize() : 0);
        int numAtomBlocks = cc.getPaddedNumAtoms()/32;
        pairValueKernel->addArg(cc.getPosq());
        pairValueKernel->addArg(cc.getNonbondedUtilities().getExclusions());
        pairValueKernel->addArg(cc.getNonbondedUtilities().getExclusionTiles());
        pairValueKernel->addArg(valueBuffers);
        if (nb.getUseCutoff()) {
            pairValueKernel->addArg(nb.getInteractingTiles());
            pairValueKernel->addArg(nb.getInteractionCount());
            for (int i = 0; i < 5; i++)
                pairValueKernel->addArg(); // Periodic box size arguments are set when the kernel is executed.
            pairValueKernel->addArg(maxTiles);
            pairValueKernel->addArg(nb.getBlockCenters());
            pairValueKernel->addArg(nb.getBlockBoundingBoxes());
            pairValueKernel->addArg(nb.getInteractingAtoms());
        }
        else
            pairValueKernel->addArg(numAtomBlocks*(numAtomBlocks+1)/2);
        if (globals.isInitialized())
            pairValueKernel->addArg(globals);
        for (int i = 0; i < (int) params->getParameterInfos().size(); i++) {
            if (pairValueUsesParam[i]) {
                ComputeParameterInfo& buffer = params->getParameterInfos()[i];
                pairValueKernel->addArg(buffer.getArray());
            }
        }
        for (auto& d : dValue0dParam)
            pairValueKernel->addArg(d);
        for (auto& function : tabulatedFunctionArrays)
            pairValueKernel->addArg(function);
        perParticleValueKernel->addArg(cc.getPosq());
        perParticleValueKernel->addArg(valueBuffers);
        if (globals.isInitialized())
            perParticleValueKernel->addArg(globals);
        for (auto& buffer : params->getParameterInfos())
            perParticleValueKernel->addArg(buffer.getArray());
        for (auto& buffer : computedValues->getParameterInfos())
            perParticleValueKernel->addArg(buffer.getArray());
        for (int i = 0; i < dValuedParam.size(); i++) {
            perParticleValueKernel->addArg(dValue0dParam[i]);
            for (int j = 0; j < dValuedParam[i]->getParameterInfos().size(); j++)
                perParticleValueKernel->addArg(dValuedParam[i]->getParameterInfos()[j].getArray());
        }
        for (auto& function : tabulatedFunctionArrays)
            perParticleValueKernel->addArg(function);
        pairEnergyKernel->addArg(cc.getLongForceBuffer());
        pairEnergyKernel->addArg(cc.getEnergyBuffer());
        pairEnergyKernel->addArg(cc.getPosq());
        pairEnergyKernel->addArg(cc.getNonbondedUtilities().getExclusions());
        pairEnergyKernel->addArg(cc.getNonbondedUtilities().getExclusionTiles());
        pairEnergyKernel->addArg(); // Whether to include energy.
        if (nb.getUseCutoff()) {
            pairEnergyKernel->addArg(nb.getInteractingTiles());
            pairEnergyKernel->addArg(nb.getInteractionCount());
            for (int i = 0; i < 5; i++)
                pairEnergyKernel->addArg(); // Periodic box size arguments are set when the kernel is executed.
            pairEnergyKernel->addArg(maxTiles);
            pairEnergyKernel->addArg(nb.getBlockCenters());
            pairEnergyKernel->addArg(nb.getBlockBoundingBoxes());
            pairEnergyKernel->addArg(nb.getInteractingAtoms());
        }
        else
            pairEnergyKernel->addArg(numAtomBlocks*(numAtomBlocks+1)/2);
        if (globals.isInitialized())
            pairEnergyKernel->addArg(globals);
        for (int i = 0; i < (int) params->getParameterInfos().size(); i++) {
            if (pairEnergyUsesParam[i]) {
                ComputeParameterInfo& buffer = params->getParameterInfos()[i];
                pairEnergyKernel->addArg(buffer.getArray());
            }
        }
        for (int i = 0; i < (int) computedValues->getParameterInfos().size(); i++) {
            if (pairEnergyUsesValue[i]) {
                ComputeParameterInfo& buffer = computedValues->getParameterInfos()[i];
                pairEnergyKernel->addArg(buffer.getArray());
            }
        }
        pairEnergyKernel->addArg(longEnergyDerivs);
        if (needEnergyParamDerivs)
            pairEnergyKernel->addArg(cc.getEnergyParamDerivBuffer());
        for (auto& function : tabulatedFunctionArrays)
            pairEnergyKernel->addArg(function);
        perParticleEnergyKernel->addArg(cc.getEnergyBuffer());
        perParticleEnergyKernel->addArg(cc.getPosq());
        perParticleEnergyKernel->addArg(cc.getLongForceBuffer());
        if (globals.isInitialized())
            perParticleEnergyKernel->addArg(globals);
        for (auto& buffer : params->getParameterInfos())
            perParticleEnergyKernel->addArg(buffer.getArray());
        for (auto& buffer : computedValues->getParameterInfos())
            perParticleEnergyKernel->addArg(buffer.getArray());
        for (auto& buffer : energyDerivs->getParameterInfos())
            perParticleEnergyKernel->addArg(buffer.getArray());
        for (auto& buffer : energyDerivChain->getParameterInfos())
            perParticleEnergyKernel->addArg(buffer.getArray());
        perParticleEnergyKernel->addArg(longEnergyDerivs);
        if (needEnergyParamDerivs)
            perParticleEnergyKernel->addArg(cc.getEnergyParamDerivBuffer());
        for (auto& function : tabulatedFunctionArrays)
            perParticleEnergyKernel->addArg(function);
        if (needParameterGradient || needEnergyParamDerivs) {
            gradientChainRuleKernel->addArg(cc.getPosq());
            gradientChainRuleKernel->addArg(cc.getLongForceBuffer());
            if (globals.isInitialized())
                gradientChainRuleKernel->addArg(globals);
            for (auto& buffer : params->getParameterInfos())
                gradientChainRuleKernel->addArg(buffer.getArray());
            for (auto& buffer : computedValues->getParameterInfos())
                gradientChainRuleKernel->addArg(buffer.getArray());
            for (auto& buffer : energyDerivs->getParameterInfos())
                gradientChainRuleKernel->addArg(buffer.getArray());
            if (needEnergyParamDerivs) {
                gradientChainRuleKernel->addArg(cc.getEnergyParamDerivBuffer());
                for (auto d : dValuedParam)
                    for (auto& buffer : d->getParameterInfos())
                        gradientChainRuleKernel->addArg(buffer.getArray());
            }
            for (auto& function : tabulatedFunctionArrays)
                gradientChainRuleKernel->addArg(function);
        }
    }
    if (globals.isInitialized()) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            float value = (float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed)
            globals.upload(globalParamValues);
    }
    pairEnergyKernel->setArg(5, (int) includeEnergy);
    if (nb.getUseCutoff()) {
        setPeriodicBoxArgs(cc, pairValueKernel, 6);
        setPeriodicBoxArgs(cc, pairEnergyKernel, 8);
        if (maxTiles < nb.getInteractingTiles().getSize()) {
            maxTiles = nb.getInteractingTiles().getSize();
            pairValueKernel->setArg(11, maxTiles);
            pairEnergyKernel->setArg(13, maxTiles);
        }
    }
    pairValueKernel->execute(nb.getNumForceThreadBlocks()*nb.getForceThreadBlockSize(), nb.getForceThreadBlockSize());
    perParticleValueKernel->execute(cc.getPaddedNumAtoms());
    pairEnergyKernel->execute(nb.getNumForceThreadBlocks()*nb.getForceThreadBlockSize(), nb.getForceThreadBlockSize());
    perParticleEnergyKernel->execute(cc.getPaddedNumAtoms());
    if (needParameterGradient || needEnergyParamDerivs)
        gradientChainRuleKernel->execute(cc.getPaddedNumAtoms());
    return 0.0;
}

void CommonCalcCustomGBForceKernel::copyParametersToContext(ContextImpl& context, const CustomGBForce& force) {
    ContextSelector selector(cc);
    int numParticles = force.getNumParticles();
    if (numParticles != cc.getNumAtoms())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");

    // Record the per-particle parameters.

    vector<vector<float> > paramVector(cc.getPaddedNumAtoms(), vector<float>(force.getNumPerParticleParameters(), 0));
    vector<double> parameters;
    for (int i = 0; i < numParticles; i++) {
        force.getParticleParameters(i, parameters);
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
    }
    params->setParameterValues(paramVector);

    // See if any tabulated functions have changed.

    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        string name = force.getTabulatedFunctionName(i);
        if (force.getTabulatedFunction(i).getUpdateCount() != tabulatedFunctionUpdateCount[name]) {
            tabulatedFunctionUpdateCount[name] = force.getTabulatedFunction(i).getUpdateCount();
            int width;
            vector<float> f = cc.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
            tabulatedFunctionArrays[i].upload(f);
        }
    }

    // Mark that the current reordering may be invalid.

    cc.invalidateMolecules(info, true, false);
}
