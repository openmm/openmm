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

#include "openmm/common/CommonCalcCustomNonbondedForceKernel.h"
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

class CommonCalcCustomNonbondedForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const CustomNonbondedForce& force) : force(force) {
        if (force.getNumInteractionGroups() > 0) {
            groupsForParticle.resize(force.getNumParticles());
            for (int i = 0; i < force.getNumInteractionGroups(); i++) {
                set<int> set1, set2;
                force.getInteractionGroupParameters(i, set1, set2);
                for (int p : set1)
                    groupsForParticle[p].insert(2*i);
                for (int p : set2)
                    groupsForParticle[p].insert(2*i+1);
            }
        }
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        thread_local static vector<double> params1, params2;
        force.getParticleParameters(particle1, params1);
        force.getParticleParameters(particle2, params2);
        for (int i = 0; i < (int) params1.size(); i++)
            if (params1[i] != params2[i])
                return false;
        if (groupsForParticle.size() > 0 && groupsForParticle[particle1] != groupsForParticle[particle2])
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
    const CustomNonbondedForce& force;
    vector<set<int> > groupsForParticle;
};

class CommonCalcCustomNonbondedForceKernel::LongRangePostComputation : public ComputeContext::ForcePostComputation {
public:
    LongRangePostComputation(ComputeContext& cc, double& longRangeCoefficient, vector<double>& longRangeCoefficientDerivs, CustomNonbondedForce* force) :
            cc(cc), longRangeCoefficient(longRangeCoefficient), longRangeCoefficientDerivs(longRangeCoefficientDerivs), force(force) {
    }
    double computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        if ((groups&(1<<force->getForceGroup())) == 0)
            return 0;
        if (!cc.getWorkThread().isCurrentThread())
            cc.getWorkThread().flush();
        Vec3 a, b, c;
        cc.getPeriodicBoxVectors(a, b, c);
        double volume = a[0]*b[1]*c[2];
        map<string, double>& derivs = cc.getEnergyParamDerivWorkspace();
        for (int i = 0; i < longRangeCoefficientDerivs.size(); i++)
            derivs[force->getEnergyParameterDerivativeName(i)] += longRangeCoefficientDerivs[i]/volume;
        return longRangeCoefficient/volume;
    }
private:
    ComputeContext& cc;
    double& longRangeCoefficient;
    vector<double>& longRangeCoefficientDerivs;
    CustomNonbondedForce* force;
};

class CommonCalcCustomNonbondedForceKernel::LongRangeTask : public ComputeContext::WorkTask {
public:
    LongRangeTask(ComputeContext& cc, Context& context, CustomNonbondedForceImpl::LongRangeCorrectionData& data,
                  double& longRangeCoefficient, vector<double>& longRangeCoefficientDerivs, CustomNonbondedForce* force) :
                        cc(cc), context(context), data(data), longRangeCoefficient(longRangeCoefficient),
                        longRangeCoefficientDerivs(longRangeCoefficientDerivs), force(force) {
    }
    void execute() {
        CustomNonbondedForceImpl::calcLongRangeCorrection(*force, data, context, longRangeCoefficient, longRangeCoefficientDerivs, cc.getThreadPool());
    }
private:
    ComputeContext& cc;
    Context& context;
    CustomNonbondedForceImpl::LongRangeCorrectionData& data;
    double& longRangeCoefficient;
    vector<double>& longRangeCoefficientDerivs;
    CustomNonbondedForce* force;
};

CommonCalcCustomNonbondedForceKernel::~CommonCalcCustomNonbondedForceKernel() {
    ContextSelector selector(cc);
    if (params != NULL)
        delete params;
    if (computedValues != NULL)
        delete computedValues;
    if (forceCopy != NULL)
        delete forceCopy;
}

void CommonCalcCustomNonbondedForceKernel::initialize(const System& system, const CustomNonbondedForce& force) {
    ContextSelector selector(cc);
    int forceIndex;
    for (forceIndex = 0; forceIndex < system.getNumForces() && &system.getForce(forceIndex) != &force; ++forceIndex)
        ;
    string prefix = (force.getNumInteractionGroups() == 0 ? "custom"+cc.intToString(forceIndex)+"_" : "");

    // Record parameters and exclusions.

    int numParticles = force.getNumParticles();
    int paddedNumParticles = cc.getPaddedNumAtoms();
    int numParams = force.getNumPerParticleParameters();
    params = new ComputeParameterSet(cc, numParams, paddedNumParticles, "customNonbondedParameters", true);
    if (force.getNumGlobalParameters() > 0)
        globals.initialize<float>(cc, force.getNumGlobalParameters(), "customNonbondedGlobals");
    vector<vector<float> > paramVector(paddedNumParticles, vector<float>(numParams, 0));
    vector<vector<int> > exclusionList(numParticles);
    for (int i = 0; i < numParticles; i++) {
        vector<double> parameters;
        force.getParticleParameters(i, parameters);
        paramVector[i].resize(parameters.size());
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
    vector<string> tableTypes;
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
        if (force.getNumInteractionGroups() == 0)
            cc.getNonbondedUtilities().addArgument(ComputeParameterInfo(tabulatedFunctionArrays[i], arrayName, "float", width));
        if (width == 1)
            tableTypes.push_back("float");
        else
            tableTypes.push_back("float"+cc.intToString(width));
        tableArgs << ", GLOBAL const float";
        if (width > 1)
            tableArgs << width;
        tableArgs << "* RESTRICT " << arrayName;
    }

    // Record information for the expressions.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
    }
    if (globals.isInitialized())
        globals.upload(globalParamValues);
    bool useCutoff = (force.getNonbondedMethod() != CustomNonbondedForce::NoCutoff);
    bool usePeriodic = (force.getNonbondedMethod() != CustomNonbondedForce::NoCutoff && force.getNonbondedMethod() != CustomNonbondedForce::CutoffNonPeriodic);
    Lepton::ParsedExpression energyExpression = Lepton::Parser::parse(force.getEnergyFunction(), functions).optimize();
    Lepton::ParsedExpression forceExpression = energyExpression.differentiate("r").optimize();
    map<string, Lepton::ParsedExpression> forceExpressions;
    forceExpressions["real customEnergy = "] = energyExpression;
    forceExpressions["tempForce -= "] = forceExpression;

    // Record which per-particle parameters and computed values appear in the energy expression.

    if (force.getNumComputedValues() > 0)
        computedValues = new ComputeParameterSet(cc, force.getNumComputedValues(), paddedNumParticles, "customNonbondedComputedValues", true);
    for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
        string name = force.getPerParticleParameterName(i);
        if (usesVariable(energyExpression, name+"1") || usesVariable(energyExpression, name+"2")) {
            paramNames.push_back(name);
            paramBuffers.push_back(params->getParameterInfos()[i]);
        }
    }
    for (int i = 0; i < force.getNumComputedValues(); i++) {
        string name, expression;
        force.getComputedValueParameters(i, name, expression);
        if (usesVariable(energyExpression, name+"1") || usesVariable(energyExpression, name+"2")) {
            computedValueNames.push_back(name);
            computedValueBuffers.push_back(computedValues->getParameterInfos()[i]);
        }
    }

    // Create the kernels.

    vector<pair<ExpressionTreeNode, string> > variables;
    ExpressionTreeNode rnode(new Operation::Variable("r"));
    variables.push_back(make_pair(rnode, "r"));
    variables.push_back(make_pair(ExpressionTreeNode(new Operation::Square(), rnode), "r2"));
    variables.push_back(make_pair(ExpressionTreeNode(new Operation::Reciprocal(), rnode), "invR"));
    for (int i = 0; i < paramNames.size(); i++) {
        variables.push_back(makeVariable(paramNames[i]+"1", "((real) "+prefix+"params"+cc.intToString(i+1)+"1)"));
        variables.push_back(makeVariable(paramNames[i]+"2", "((real) "+prefix+"params"+cc.intToString(i+1)+"2)"));
    }
    for (int i = 0; i < computedValueNames.size(); i++) {
        variables.push_back(makeVariable(computedValueNames[i]+"1", prefix+"values"+cc.intToString(i+1)+"1"));
        variables.push_back(makeVariable(computedValueNames[i]+"2", prefix+"values"+cc.intToString(i+1)+"2"));
    }
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        const string& name = force.getGlobalParameterName(i);
        string value = "globals["+cc.intToString(i)+"]";
        variables.push_back(makeVariable(name, prefix+value));
    }
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string paramName = force.getEnergyParameterDerivativeName(i);
        string derivVariable = cc.getNonbondedUtilities().addEnergyParameterDerivative(paramName);
        Lepton::ParsedExpression derivExpression = energyExpression.differentiate(paramName).optimize();
        forceExpressions[derivVariable+" += interactionScale*switchValue*"] = derivExpression;
    }
    stringstream compute;
    compute << cc.getExpressionUtilities().createExpressions(forceExpressions, variables, functionList, functionDefinitions, prefix+"temp");
    map<string, string> replacements;
    replacements["COMPUTE_FORCE"] = compute.str();
    replacements["USE_SWITCH"] = (useCutoff && force.getUseSwitchingFunction() ? "1" : "0");
    if (force.getUseSwitchingFunction()) {
        // Compute the switching coefficients.

        replacements["SWITCH_CUTOFF"] = cc.doubleToString(force.getSwitchingDistance());
        replacements["SWITCH_C3"] = cc.doubleToString(10/pow(force.getSwitchingDistance()-force.getCutoffDistance(), 3.0));
        replacements["SWITCH_C4"] = cc.doubleToString(15/pow(force.getSwitchingDistance()-force.getCutoffDistance(), 4.0));
        replacements["SWITCH_C5"] = cc.doubleToString(6/pow(force.getSwitchingDistance()-force.getCutoffDistance(), 5.0));
    }
    string source = cc.replaceStrings(CommonKernelSources::customNonbonded, replacements);
    if (force.getNumInteractionGroups() > 0)
        initInteractionGroups(force, source, tableTypes);
    else {
        cc.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, true, force.getCutoffDistance(), exclusionList, source, force.getForceGroup(), numParticles > 2000);
        for (int i = 0; i < paramBuffers.size(); i++)
            cc.getNonbondedUtilities().addParameter(ComputeParameterInfo(paramBuffers[i].getArray(), prefix+"params"+cc.intToString(i+1),
                    paramBuffers[i].getComponentType(), paramBuffers[i].getNumComponents()));
        for (int i = 0; i < computedValueBuffers.size(); i++)
            cc.getNonbondedUtilities().addParameter(ComputeParameterInfo(computedValueBuffers[i].getArray(), prefix+"values"+cc.intToString(i+1),
                    computedValueBuffers[i].getComponentType(), computedValueBuffers[i].getNumComponents()));
        if (globals.isInitialized()) {
            globals.upload(globalParamValues);
            cc.getNonbondedUtilities().addArgument(ComputeParameterInfo(globals, prefix+"globals", "float", 1));
        }
    }
    if (force.getNumComputedValues() > 0) {
        // Create the kernel to calculate computed values.

        stringstream valuesSource, args;
        for (int i = 0; i < computedValues->getParameterInfos().size(); i++) {
            ComputeParameterInfo& buffer = computedValues->getParameterInfos()[i];
            string valueName = "values"+cc.intToString(i+1);
            if (i > 0)
                args << ", ";
            args << "GLOBAL " << buffer.getType() << "* RESTRICT global_" << valueName;
            valuesSource << buffer.getType() << " local_" << valueName << ";\n";
        }
        if (force.getNumGlobalParameters() > 0)
            args << ", GLOBAL const float* globals";
        for (int i = 0; i < params->getParameterInfos().size(); i++) {
            ComputeParameterInfo& buffer = params->getParameterInfos()[i];
            string paramName = "params"+cc.intToString(i+1);
            args << ", GLOBAL const " << buffer.getType() << "* RESTRICT " << paramName;
        }
        map<string, string> variables;
        for (int i = 0; i < force.getNumPerParticleParameters(); i++)
            variables[force.getPerParticleParameterName(i)] = "params"+params->getParameterSuffix(i, "[index]");
        for (int i = 0; i < force.getNumGlobalParameters(); i++)
            variables[force.getGlobalParameterName(i)] = "globals["+cc.intToString(i)+"]";
        for (int i = 0; i < force.getNumComputedValues(); i++) {
            string name, expression;
            force.getComputedValueParameters(i, name, expression);
            variables[name] = "local_values"+computedValues->getParameterSuffix(i);
            map<string, Lepton::ParsedExpression> valueExpressions;
            valueExpressions["local_values"+computedValues->getParameterSuffix(i)+" = "] = Lepton::Parser::parse(expression, functions).optimize();
            valuesSource << cc.getExpressionUtilities().createExpressions(valueExpressions, variables, functionList, functionDefinitions, "value"+cc.intToString(i)+"_temp");
        }
        for (int i = 0; i < (int) computedValues->getParameterInfos().size(); i++) {
            string valueName = "values"+cc.intToString(i+1);
            valuesSource << "global_" << valueName << "[index] = local_" << valueName << ";\n";
        }
        map<string, string> replacements;
        replacements["PARAMETER_ARGUMENTS"] = args.str()+tableArgs.str();
        replacements["COMPUTE_VALUES"] = valuesSource.str();
        map<string, string> defines;
        defines["NUM_ATOMS"] = cc.intToString(cc.getNumAtoms());
        ComputeProgram program = cc.compileProgram(cc.replaceStrings(CommonKernelSources::customNonbondedComputedValues, replacements), defines);
        computedValuesKernel = program->createKernel("computePerParticleValues");
        for (auto& value : computedValues->getParameterInfos())
            computedValuesKernel->addArg(value.getArray());
        if (globals.isInitialized())
            computedValuesKernel->addArg(globals);
        for (auto& parameter : params->getParameterInfos())
            computedValuesKernel->addArg(parameter.getArray());
        for (auto& function : tabulatedFunctionArrays)
            computedValuesKernel->addArg(function);
    }
    info = new ForceInfo(force);
    cc.addForce(info);

    // Record information for the long range correction.

    if (force.getNonbondedMethod() == CustomNonbondedForce::CutoffPeriodic && force.getUseLongRangeCorrection() && cc.getContextIndex() == 0) {
        forceCopy = new CustomNonbondedForce(force);
        longRangeCorrectionData = CustomNonbondedForceImpl::prepareLongRangeCorrection(force, cc.getThreadPool().getNumThreads());
        cc.addPostComputation(new LongRangePostComputation(cc, longRangeCoefficient, longRangeCoefficientDerivs, forceCopy));
        hasInitializedLongRangeCorrection = false;
    }
    else {
        longRangeCoefficient = 0.0;
        hasInitializedLongRangeCorrection = true;
    }
}

void CommonCalcCustomNonbondedForceKernel::initInteractionGroups(const CustomNonbondedForce& force, const string& interactionSource, const vector<string>& tableTypes) {
    // Process groups to form tiles.

    vector<vector<int> > atomLists;
    vector<pair<int, int> > tiles;
    vector<int> tileGroup;
    vector<vector<int> > duplicateAtomsForGroup;
    for (int group = 0; group < force.getNumInteractionGroups(); group++) {
        // Get the list of atoms in this group and sort them.

        set<int> set1, set2;
        force.getInteractionGroupParameters(group, set1, set2);
        vector<int> atoms1, atoms2;
        atoms1.insert(atoms1.begin(), set1.begin(), set1.end());
        atoms2.insert(atoms2.begin(), set2.begin(), set2.end());
        sort(atoms1.begin(), atoms1.end());
        sort(atoms2.begin(), atoms2.end());
        duplicateAtomsForGroup.push_back(vector<int>());
        set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(),
                inserter(duplicateAtomsForGroup[group], duplicateAtomsForGroup[group].begin()));
        sort(duplicateAtomsForGroup[group].begin(), duplicateAtomsForGroup[group].end());

        // Find how many tiles we will create for this group.

        int tileWidth = min(min(32, (int) atoms1.size()), (int) atoms2.size());
        if (tileWidth == 0)
            continue;
        int numBlocks1 = (atoms1.size()+tileWidth-1)/tileWidth;
        int numBlocks2 = (atoms2.size()+tileWidth-1)/tileWidth;

        // Add the tiles.

        int firstTile = tiles.size();
        for (int i = 0; i < numBlocks1; i++)
            for (int j = 0; j < numBlocks2; j++) {
                tiles.push_back(make_pair(atomLists.size()+i, atomLists.size()+numBlocks1+j));
                tileGroup.push_back(group);
            }

        // Add the atom lists.

        for (int i = 0; i < numBlocks1; i++) {
            vector<int> atoms;
            int first = i*tileWidth;
            int last = min((i+1)*tileWidth, (int) atoms1.size());
            for (int j = first; j < last; j++)
                atoms.push_back(atoms1[j]);
            atomLists.push_back(atoms);
        }
        for (int i = 0; i < numBlocks2; i++) {
            vector<int> atoms;
            int first = i*tileWidth;
            int last = min((i+1)*tileWidth, (int) atoms2.size());
            for (int j = first; j < last; j++)
                atoms.push_back(atoms2[j]);
            atomLists.push_back(atoms);
        }
    }

    // Build a lookup table for quickly identifying excluded interactions.

    vector<set<int> > exclusions(force.getNumParticles());
    for (int i = 0; i < force.getNumExclusions(); i++) {
        int p1, p2;
        force.getExclusionParticles(i, p1, p2);
        exclusions[p1].insert(p2);
        exclusions[p2].insert(p1);
    }

    // Build the exclusion flags for each tile.  While we're at it, filter out tiles
    // where all interactions are excluded, and sort the tiles by size.

    vector<vector<int> > exclusionFlags(tiles.size());
    vector<pair<int, int> > tileOrder;
    for (int tile = 0; tile < tiles.size(); tile++) {
        bool swapped = false;
        if (atomLists[tiles[tile].first].size() < atomLists[tiles[tile].second].size()) {
            // For efficiency, we want the first axis to be the larger one.

            int swap = tiles[tile].first;
            tiles[tile].first = tiles[tile].second;
            tiles[tile].second = swap;
            swapped = true;
        }
        vector<int>& atoms1 = atomLists[tiles[tile].first];
        vector<int>& atoms2 = atomLists[tiles[tile].second];
        vector<int>& duplicateAtoms = duplicateAtomsForGroup[tileGroup[tile]];
        vector<int>& flags = exclusionFlags[tile];
        flags.resize(atoms1.size(), (int) (1LL<<atoms2.size())-1);
        int numExcluded = 0;
        for (int i = 0; i < (int) atoms1.size(); i++) {
            int a1 = atoms1[i];
            bool a1IsDuplicate = binary_search(duplicateAtoms.begin(), duplicateAtoms.end(), a1);
            for (int j = 0; j < (int) atoms2.size(); j++) {
                int a2 = atoms2[j];
                bool isExcluded = false;
                if (a1 == a2 || exclusions[a1].find(a2) != exclusions[a1].end())
                    isExcluded = true; // This is an excluded interaction.
                else if ((a1 > a2) == swapped && a1IsDuplicate && binary_search(duplicateAtoms.begin(), duplicateAtoms.end(), a2))
                    isExcluded = true; // Both atoms are in both sets, so skip duplicate interactions.
                if (isExcluded) {
                    flags[i] &= -1-(1<<j);
                    numExcluded++;
                }
            }
        }
        if (numExcluded == atoms1.size()*atoms2.size())
            continue; // All interactions are excluded.
        tileOrder.push_back(make_pair(-((int)atoms2.size()), tile));
    }
    sort(tileOrder.begin(), tileOrder.end());

    // Merge tiles to get as close as possible to 32 along the first axis of each one.

    vector<int> tileSetStart;
    tileSetStart.push_back(0);
    int tileSetSize = 0;
    for (int i = 0; i < tileOrder.size(); i++) {
        int tile = tileOrder[i].second;
        int size = atomLists[tiles[tile].first].size();
        if (tileSetSize+size > 32) {
            tileSetStart.push_back(i);
            tileSetSize = 0;
        }
        tileSetSize += size;
    }
    tileSetStart.push_back(tileOrder.size());

    // Build the data structures.

    int numTileSets = tileSetStart.size()-1;
    vector<mm_int4> groupData;
    for (int tileSet = 0; tileSet < numTileSets; tileSet++) {
        int indexInTileSet = 0;
        int minSize = 0;
        if (cc.getSIMDWidth() < 32) {
            // We need to include a barrier inside the inner loop, so ensure that all
            // threads will loop the same number of times.

            for (int i = tileSetStart[tileSet]; i < tileSetStart[tileSet+1]; i++)
                minSize = max(minSize, (int) atomLists[tiles[tileOrder[i].second].first].size());
        }
        for (int i = tileSetStart[tileSet]; i < tileSetStart[tileSet+1]; i++) {
            int tile = tileOrder[i].second;
            vector<int>& atoms1 = atomLists[tiles[tile].first];
            vector<int>& atoms2 = atomLists[tiles[tile].second];
            int range = indexInTileSet + ((indexInTileSet+max(minSize, (int) atoms1.size()))<<16);
            int allFlags = (1<<atoms2.size())-1;
            for (int j = 0; j < (int) atoms1.size(); j++) {
                int a1 = atoms1[j];
                int a2 = (j < atoms2.size() ? atoms2[j] : 0);
                int flags = (exclusionFlags[tile].size() > 0 ? exclusionFlags[tile][j] : allFlags);
                groupData.push_back(mm_int4(a1, a2, range, flags<<indexInTileSet));
            }
            indexInTileSet += atoms1.size();
        }
        for (; indexInTileSet < 32; indexInTileSet++)
            groupData.push_back(mm_int4(0, 0, minSize<<16, 0));
    }
    interactionGroupData.initialize<mm_int4>(cc, groupData.size(), "interactionGroupData");
    interactionGroupData.upload(groupData);
    numGroupTiles.initialize<int>(cc, 1, "numGroupTiles");

    // Allocate space for a neighbor list, if necessary.

    if (force.getNonbondedMethod() != CustomNonbondedForce::NoCutoff && groupData.size() > cc.getNumThreadBlocks()) {
        filteredGroupData.initialize<mm_int4>(cc, groupData.size(), "filteredGroupData");
        interactionGroupData.copyTo(filteredGroupData);
        int numTiles = groupData.size()/32;
        numGroupTiles.upload(&numTiles);
    }

    // Create the kernel.

    hasParamDerivs = (force.getNumEnergyParameterDerivatives() > 0);
    map<string, string> replacements;
    replacements["COMPUTE_INTERACTION"] = interactionSource;
    const string suffixes[] = {"x", "y", "z", "w"};
    stringstream localData;
    int localDataSize = 0;
    for (int i = 0; i < paramBuffers.size(); i++) {
        localData<<paramBuffers[i].getComponentType()<<" params"<<(i+1)<<";\n";
        localDataSize += paramBuffers[i].getSize();
    }
    for (int i = 0; i < computedValueBuffers.size(); i++) {
        localData<<computedValueBuffers[i].getComponentType()<<" values"<<(i+1)<<";\n";
        localDataSize += computedValueBuffers[i].getSize();
    }
    replacements["ATOM_PARAMETER_DATA"] = localData.str();
    stringstream args;
    for (int i = 0; i < paramBuffers.size(); i++)
        args<<", GLOBAL const "<<paramBuffers[i].getType()<<"* RESTRICT global_params"<<(i+1);
    for (int i = 0; i < computedValueBuffers.size(); i++)
        args<<", GLOBAL const "<<computedValueBuffers[i].getType()<<"* RESTRICT global_values"<<(i+1);
    for (int i = 0; i < tabulatedFunctionArrays.size(); i++)
        args << ", GLOBAL const " << tableTypes[i]<< "* RESTRICT table" << i;
    if (globals.isInitialized())
        args<<", GLOBAL const float* RESTRICT globals";
    if (hasParamDerivs)
        args << ", GLOBAL mixed* RESTRICT energyParamDerivs";
    replacements["PARAMETER_ARGUMENTS"] = args.str();
    stringstream load1;
    for (int i = 0; i < paramBuffers.size(); i++)
        load1<<paramBuffers[i].getType()<<" params"<<(i+1)<<"1 = global_params"<<(i+1)<<"[atom1];\n";
    for (int i = 0; i < computedValueBuffers.size(); i++)
        load1<<computedValueBuffers[i].getType()<<" values"<<(i+1)<<"1 = global_values"<<(i+1)<<"[atom1];\n";
    replacements["LOAD_ATOM1_PARAMETERS"] = load1.str();
    stringstream loadLocal2;
    for (int i = 0; i < paramBuffers.size(); i++)
        loadLocal2<<"localData[LOCAL_ID].params"<<(i+1)<<" = global_params"<<(i+1)<<"[atom2];\n";
    for (int i = 0; i < computedValueBuffers.size(); i++)
        loadLocal2<<"localData[LOCAL_ID].values"<<(i+1)<<" = global_values"<<(i+1)<<"[atom2];\n";
    replacements["LOAD_LOCAL_PARAMETERS"] = loadLocal2.str();
    stringstream load2;
    for (int i = 0; i < paramBuffers.size(); i++)
        load2<<paramBuffers[i].getType()<<" params"<<(i+1)<<"2 = localData[localIndex].params"<<(i+1)<<";\n";
    for (int i = 0; i < computedValueBuffers.size(); i++)
        load2<<computedValueBuffers[i].getType()<<" values"<<(i+1)<<"2 = localData[localIndex].values"<<(i+1)<<";\n";
    replacements["LOAD_ATOM2_PARAMETERS"] = load2.str();
    stringstream initDerivs, saveDerivs;
    const vector<string>& allParamDerivNames = cc.getEnergyParamDerivNames();
    int numDerivs = allParamDerivNames.size();
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string paramName = force.getEnergyParameterDerivativeName(i);
        string derivVariable = cc.getNonbondedUtilities().addEnergyParameterDerivative(paramName);
        initDerivs<<"mixed "<<derivVariable<<" = 0;\n";
        for (int index = 0; index < numDerivs; index++)
            if (allParamDerivNames[index] == paramName)
                saveDerivs<<"energyParamDerivs[GLOBAL_ID*numDerivatives+"<<index<<"] += "<<derivVariable<<";\n";
    }
    replacements["INIT_DERIVATIVES"] = initDerivs.str();
    replacements["SAVE_DERIVATIVES"] = saveDerivs.str();
    map<string, string> defines;
    if (force.getNonbondedMethod() != CustomNonbondedForce::NoCutoff)
        defines["USE_CUTOFF"] = "1";
    if (force.getNonbondedMethod() == CustomNonbondedForce::CutoffPeriodic)
        defines["USE_PERIODIC"] = "1";
    int localMemorySize = max(32, cc.getNonbondedUtilities().getForceThreadBlockSize());
    defines["LOCAL_MEMORY_SIZE"] = cc.intToString(localMemorySize);
    defines["WARPS_IN_BLOCK"] = cc.intToString(localMemorySize/32);
    double cutoff = force.getCutoffDistance();
    defines["CUTOFF_SQUARED"] = cc.doubleToString(cutoff*cutoff);
    double paddedCutoff = cc.getNonbondedUtilities().padCutoff(cutoff);
    defines["PADDED_CUTOFF_SQUARED"] = cc.doubleToString(paddedCutoff*paddedCutoff);
    defines["PADDED_NUM_ATOMS"] = cc.intToString(cc.getPaddedNumAtoms());
    defines["TILE_SIZE"] = "32";
    defines["NUM_TILES"] = cc.intToString(numTileSets);
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*numTileSets/numContexts;
    int endIndex = (cc.getContextIndex()+1)*numTileSets/numContexts;
    defines["FIRST_TILE"] = cc.intToString(startIndex);
    defines["LAST_TILE"] = cc.intToString(endIndex);
    if ((localDataSize/4)%2 == 0 && !cc.getUseDoublePrecision())
        defines["PARAMETER_SIZE_IS_EVEN"] = "1";
    ComputeProgram program = cc.compileProgram(cc.replaceStrings(CommonKernelSources::customNonbondedGroups, replacements), defines);
    interactionGroupKernel = program->createKernel("computeInteractionGroups");
    prepareNeighborListKernel = program->createKernel("prepareToBuildNeighborList");
    buildNeighborListKernel = program->createKernel("buildNeighborList");
    numGroupThreadBlocks = cc.getNonbondedUtilities().getNumForceThreadBlocks();
}

double CommonCalcCustomNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    useNeighborList = (filteredGroupData.isInitialized() && cc.getNonbondedUtilities().getUseCutoff());
    if (useNeighborList && cc.getContextIndex() > 0) {
        // When using a neighbor list, run the whole calculation on a single device.
        return 0.0;
    }
    ContextSelector selector(cc);
    bool recomputeLongRangeCorrection = !hasInitializedLongRangeCorrection;
    if (globals.isInitialized()) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            float value = (float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed) {
            globals.upload(globalParamValues);
            if (forceCopy != NULL)
                recomputeLongRangeCorrection = true;
        }
    }
    if (recomputeLongRangeCorrection) {
        if (includeEnergy || forceCopy->getNumEnergyParameterDerivatives() > 0) {
            cc.getWorkThread().addTask(new LongRangeTask(cc, context.getOwner(), longRangeCorrectionData, longRangeCoefficient, longRangeCoefficientDerivs, forceCopy));
            hasInitializedLongRangeCorrection = true;
        }
        else
            hasInitializedLongRangeCorrection = false;
    }
    if (computedValues != NULL)
        computedValuesKernel->execute(cc.getNumAtoms());
    if (interactionGroupData.isInitialized()) {
        if (!hasInitializedKernel) {
            hasInitializedKernel = true;
            interactionGroupKernel->addArg(cc.getLongForceBuffer());
            interactionGroupKernel->addArg(cc.getEnergyBuffer());
            interactionGroupKernel->addArg(cc.getPosq());
            interactionGroupKernel->addArg((useNeighborList ? filteredGroupData : interactionGroupData));
            interactionGroupKernel->addArg(numGroupTiles);
            interactionGroupKernel->addArg((int) useNeighborList);
            for (int i = 0; i < 5; i++)
                interactionGroupKernel->addArg(); // Periodic box information will be set just before it is executed.
            interactionGroupKernel->addArg((int) cc.getEnergyParamDerivNames().size());
            for (auto& buffer : paramBuffers)
                interactionGroupKernel->addArg(buffer.getArray());
            for (auto& buffer : computedValueBuffers)
                interactionGroupKernel->addArg(buffer.getArray());
            for (auto& function : tabulatedFunctionArrays)
                interactionGroupKernel->addArg(function);
            if (globals.isInitialized())
                interactionGroupKernel->addArg(globals);
            if (hasParamDerivs)
                interactionGroupKernel->addArg(cc.getEnergyParamDerivBuffer());
            if (useNeighborList) {
                // Initialize kernels for building the interaction group neighbor list.

                prepareNeighborListKernel->addArg(cc.getNonbondedUtilities().getRebuildNeighborList());
                prepareNeighborListKernel->addArg(numGroupTiles);
                buildNeighborListKernel->addArg(cc.getNonbondedUtilities().getRebuildNeighborList());
                buildNeighborListKernel->addArg(numGroupTiles);
                buildNeighborListKernel->addArg(cc.getPosq());
                buildNeighborListKernel->addArg(interactionGroupData);
                buildNeighborListKernel->addArg(filteredGroupData);
                for (int i = 0; i < 5; i++)
                    buildNeighborListKernel->addArg(); // Periodic box information will be set just before it is executed.
            }
        }
        int forceThreadBlockSize = max(32, cc.getNonbondedUtilities().getForceThreadBlockSize());
        if (useNeighborList) {
            // Rebuild the neighbor list, if necessary.

            setPeriodicBoxArgs(cc, buildNeighborListKernel, 5);
            prepareNeighborListKernel->execute(1, 1);
            buildNeighborListKernel->execute(numGroupThreadBlocks*forceThreadBlockSize, forceThreadBlockSize);
        }
        setPeriodicBoxArgs(cc, interactionGroupKernel, 6);
        interactionGroupKernel->execute(numGroupThreadBlocks*forceThreadBlockSize, forceThreadBlockSize);
    }
    return 0;
}

void CommonCalcCustomNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const CustomNonbondedForce& force, int firstParticle, int lastParticle) {
    ContextSelector selector(cc);
    int numParticles = force.getNumParticles();
    if (numParticles != cc.getNumAtoms())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");

    // Record the per-particle parameters.

    if (firstParticle <= lastParticle) {
        int numToSet = lastParticle-firstParticle+1;
        int numParams = force.getNumPerParticleParameters();
        vector<vector<float> > paramVector(numToSet, vector<float>(numParams, 0));
        vector<double> parameters;
        for (int i = 0; i < numToSet; i++) {
            force.getParticleParameters(firstParticle+i, parameters);
            paramVector[i].resize(parameters.size());
            for (int j = 0; j < (int) parameters.size(); j++)
                paramVector[i][j] = (float) parameters[j];
        }
        params->setParameterValuesSubset(firstParticle, paramVector);
    }

    // If necessary, recompute the long range correction.

    if (forceCopy != NULL) {
        longRangeCorrectionData = CustomNonbondedForceImpl::prepareLongRangeCorrection(force, cc.getThreadPool().getNumThreads());
        CustomNonbondedForceImpl::calcLongRangeCorrection(force, longRangeCorrectionData, context.getOwner(), longRangeCoefficient, longRangeCoefficientDerivs, cc.getThreadPool());
        hasInitializedLongRangeCorrection = false;
        *forceCopy = force;
    }

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

    cc.invalidateMolecules(info, firstParticle <= lastParticle, false);
}
