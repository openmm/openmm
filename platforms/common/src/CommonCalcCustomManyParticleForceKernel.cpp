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

#include "openmm/common/CommonCalcCustomManyParticleForceKernel.h"
#include "openmm/common/CommonKernelUtilities.h"
#include "openmm/common/ContextSelector.h"
#include "openmm/common/ExpressionUtilities.h"
#include "openmm/Context.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/CustomManyParticleForceImpl.h"
#include "CommonKernelSources.h"
#include "SimTKOpenMMRealType.h"
#include "lepton/CustomFunction.h"
#include "lepton/ExpressionTreeNode.h"
#include "lepton/Operation.h"
#include "lepton/Parser.h"
#include "lepton/ParsedExpression.h"

using namespace OpenMM;
using namespace std;
using namespace Lepton;

class CommonCalcCustomManyParticleForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const CustomManyParticleForce& force) : force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        thread_local static vector<double> params1, params2;
        int type1, type2;
        force.getParticleParameters(particle1, params1, type1);
        force.getParticleParameters(particle2, params2, type2);
        if (type1 != type2)
            return false;
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
    const CustomManyParticleForce& force;
};

CommonCalcCustomManyParticleForceKernel::~CommonCalcCustomManyParticleForceKernel() {
    ContextSelector selector(cc);
    if (params != NULL)
        delete params;
}

void CommonCalcCustomManyParticleForceKernel::initialize(const System& system, const CustomManyParticleForce& force) {
    ContextSelector selector(cc);
    int numParticles = force.getNumParticles();
    int particlesPerSet = force.getNumParticlesPerSet();
    bool centralParticleMode = (force.getPermutationMode() == CustomManyParticleForce::UniqueCentralParticle);
    nonbondedMethod = CalcCustomManyParticleForceKernel::NonbondedMethod(force.getNonbondedMethod());
    forceWorkgroupSize = 128;
    findNeighborsWorkgroupSize = (cc.getSIMDWidth() >= 32 ? 128 : 32);

    // Record parameter values.

    params = new ComputeParameterSet(cc, force.getNumPerParticleParameters(), numParticles, "customManyParticleParameters");
    vector<vector<float> > paramVector(numParticles);
    for (int i = 0; i < numParticles; i++) {
        vector<double> parameters;
        int type;
        force.getParticleParameters(i, parameters, type);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
    }
    params->setParameterValues(paramVector);
    info = new ForceInfo(force);
    cc.addForce(info);

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
        string arrayName = "table"+cc.intToString(i);
        functionDefinitions.push_back(make_pair(name, arrayName));
        functions[name] = cc.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cc.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
        tabulatedFunctionArrays[i].initialize<float>(cc, f.size(), "TabulatedFunction");
        tabulatedFunctionArrays[i].upload(f);
        tableArgs << ", GLOBAL const float";
        if (width > 1)
            tableArgs << width;
        tableArgs << "* RESTRICT " << arrayName;
    }

    // Record information about parameters.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
    }
    vector<pair<ExpressionTreeNode, string> > variables;
    for (int i = 0; i < particlesPerSet; i++) {
        string index = cc.intToString(i+1);
        variables.push_back(makeVariable("x"+index, "pos"+index+".x"));
        variables.push_back(makeVariable("y"+index, "pos"+index+".y"));
        variables.push_back(makeVariable("z"+index, "pos"+index+".z"));
    }
    for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
        const string& name = force.getPerParticleParameterName(i);
        for (int j = 0; j < particlesPerSet; j++) {
            string index = cc.intToString(j+1);
            variables.push_back(makeVariable(name+index, "((real) params"+params->getParameterSuffix(i, index)+")"));
        }
    }
    if (force.getNumGlobalParameters() > 0) {
        globals.initialize<float>(cc, force.getNumGlobalParameters(), "customManyParticleGlobals");
        globals.upload(globalParamValues);
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = "globals["+cc.intToString(i)+"]";
            variables.push_back(makeVariable(name, value));
        }
    }

    // Build data structures for type filters.

    vector<int> particleTypesVec;
    vector<int> orderIndexVec;
    vector<std::vector<int> > particleOrderVec;
    int numTypes;
    CustomManyParticleForceImpl::buildFilterArrays(force, numTypes, particleTypesVec, orderIndexVec, particleOrderVec);
    bool hasTypeFilters = (particleOrderVec.size() > 1);
    if (hasTypeFilters) {
        particleTypes.initialize<int>(cc, particleTypesVec.size(), "customManyParticleTypes");
        orderIndex.initialize<int>(cc, orderIndexVec.size(), "customManyParticleOrderIndex");
        particleOrder.initialize<int>(cc, particleOrderVec.size()*particlesPerSet, "customManyParticleOrder");
        particleTypes.upload(particleTypesVec);
        orderIndex.upload(orderIndexVec);
        vector<int> flattenedOrder(particleOrder.getSize());
        for (int i = 0; i < (int) particleOrderVec.size(); i++)
            for (int j = 0; j < particlesPerSet; j++)
                flattenedOrder[i*particlesPerSet+j] = particleOrderVec[i][j];
        particleOrder.upload(flattenedOrder);
    }

    // Build data structures for exclusions.

    if (force.getNumExclusions() > 0) {
        vector<vector<int> > particleExclusions(numParticles);
        for (int i = 0; i < force.getNumExclusions(); i++) {
            int p1, p2;
            force.getExclusionParticles(i, p1, p2);
            particleExclusions[p1].push_back(p2);
            particleExclusions[p2].push_back(p1);
        }
        vector<int> exclusionsVec;
        vector<int> exclusionStartIndexVec(numParticles+1);
        exclusionStartIndexVec[0] = 0;
        for (int i = 0; i < numParticles; i++) {
            sort(particleExclusions[i].begin(), particleExclusions[i].end());
            exclusionsVec.insert(exclusionsVec.end(), particleExclusions[i].begin(), particleExclusions[i].end());
            exclusionStartIndexVec[i+1] = exclusionsVec.size();
        }
        exclusions.initialize<int>(cc, exclusionsVec.size(), "customManyParticleExclusions");
        exclusionStartIndex.initialize<int>(cc, exclusionStartIndexVec.size(), "customManyParticleExclusionStart");
        exclusions.upload(exclusionsVec);
        exclusionStartIndex.upload(exclusionStartIndexVec);
    }

    // Build data structures for the neighbor list.

    int numAtomBlocks = cc.getPaddedNumAtoms()/32;
    if (nonbondedMethod != NoCutoff) {
        int elementSize = (cc.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
        blockCenter.initialize(cc, numAtomBlocks, 4*elementSize, "blockCenter");
        blockBoundingBox.initialize(cc, numAtomBlocks, 4*elementSize, "blockBoundingBox");
        numNeighborPairs.initialize<int>(cc, 1, "customManyParticleNumNeighborPairs");
        neighborStartIndex.initialize<int>(cc, numParticles+1, "customManyParticleNeighborStartIndex");
        numNeighborsForAtom.initialize<int>(cc, numParticles, "customManyParticleNumNeighborsForAtom");

        // Select a size for the array that holds the neighbor list.  We have to make a fairly
        // arbitrary guess, but if this turns out to be too small we'll increase it later.

        maxNeighborPairs = 150*numParticles;
        neighborPairs.initialize<mm_int2>(cc, maxNeighborPairs, "customManyParticleNeighborPairs");
        neighbors.initialize<int>(cc, maxNeighborPairs, "customManyParticleNeighbors");
    }

    // Generate the kernel.

    Lepton::ParsedExpression energyExpression = CustomManyParticleForceImpl::prepareExpression(force, functions);
    map<string, Lepton::ParsedExpression> forceExpressions;
    stringstream compute;
    for (int i = 0; i < (int) params->getParameterInfos().size(); i++) {
        ComputeParameterInfo& parameter = params->getParameterInfos()[i];
        compute<<parameter.getType()<<" params"<<(i+1)<<" = global_params"<<(i+1)<<"[index];\n";
    }
    forceExpressions["energy += "] = energyExpression;
    vector<string> forceNames;
    for (int i = 0; i < particlesPerSet; i++) {
        string istr = cc.intToString(i+1);
        string forceName = "force"+istr;
        forceNames.push_back(forceName);
        compute<<"real3 "<<forceName<<" = make_real3(0);\n";
        Lepton::ParsedExpression forceExpressionX = energyExpression.differentiate("x"+istr).optimize();
        Lepton::ParsedExpression forceExpressionY = energyExpression.differentiate("y"+istr).optimize();
        Lepton::ParsedExpression forceExpressionZ = energyExpression.differentiate("z"+istr).optimize();
        if (!isZeroExpression(forceExpressionX))
            forceExpressions[forceName+".x -= "] = forceExpressionX;
        if (!isZeroExpression(forceExpressionY))
            forceExpressions[forceName+".y -= "] = forceExpressionY;
        if (!isZeroExpression(forceExpressionZ))
            forceExpressions[forceName+".z -= "] = forceExpressionZ;
    }
    compute << cc.getExpressionUtilities().createExpressions(forceExpressions, variables, functionList, functionDefinitions, "temp", "real", force.usesPeriodicBoundaryConditions());

    // Store forces to global memory.

    for (int i = 0; i < particlesPerSet; i++)
        compute<<"storeForce(atom"<<(i+1)<<", "<<forceNames[i]<<", forceBuffers);\n";

    // Create other replacements that depend on the number of particles per set.

    stringstream numCombinations, atomsForCombination, isValidCombination, permute, loadData, verifyCutoff, verifyExclusions;
    if (hasTypeFilters) {
        permute<<"int particleSet[] = {";
        for (int i = 0; i < particlesPerSet; i++) {
            permute<<"p"<<(i+1);
            if (i < particlesPerSet-1)
                permute<<", ";
        }
        permute<<"};\n";
    }
    for (int i = 0; i < particlesPerSet; i++) {
        if (hasTypeFilters)
            permute<<"int atom"<<(i+1)<<" = particleSet[particleOrder["<<particlesPerSet<<"*order+"<<i<<"]];\n";
        else
            permute<<"int atom"<<(i+1)<<" = p"<<(i+1)<<";\n";
        loadData<<"real3 pos"<<(i+1)<<" = trimTo3(posq[atom"<<(i+1)<<"]);\n";
        for (int j = 0; j < (int) params->getParameterInfos().size(); j++)
            loadData<<params->getParameterInfos()[j].getType()<<" params"<<(j+1)<<(i+1)<<" = global_params"<<(j+1)<<"[atom"<<(i+1)<<"];\n";
    }
    if (centralParticleMode) {
        for (int i = 1; i < particlesPerSet; i++) {
            if (i > 1)
                isValidCombination<<" && p"<<(i+1)<<">p"<<i<<" && ";
            isValidCombination<<"p"<<(i+1)<<"!=p1";
        }
    }
    else {
        for (int i = 2; i < particlesPerSet; i++) {
            if (i > 2)
                isValidCombination<<" && ";
            isValidCombination<<"a"<<(i+1)<<">a"<<i;
        }
    }
    atomsForCombination<<"int tempIndex = index;\n";
    for (int i = 1; i < particlesPerSet; i++) {
        if (i > 1)
            numCombinations<<"*";
        numCombinations<<"numNeighbors";
        if (centralParticleMode)
            atomsForCombination<<"int a"<<(i+1)<<" = tempIndex%numNeighbors;\n";
        else
            atomsForCombination<<"int a"<<(i+1)<<" = 1+tempIndex%numNeighbors;\n";
        if (i < particlesPerSet-1)
            atomsForCombination<<"tempIndex /= numNeighbors;\n";
    }
    if (particlesPerSet > 2) {
        if (centralParticleMode)
            atomsForCombination<<"a2 = (a3%2 == 0 ? a2 : numNeighbors-a2-1);\n";
        else
            atomsForCombination<<"a2 = (a3%2 == 0 ? a2 : numNeighbors-a2+1);\n";
    }
    for (int i = 1; i < particlesPerSet; i++) {
        if (nonbondedMethod == NoCutoff) {
            if (centralParticleMode)
                atomsForCombination<<"int p"<<(i+1)<<" = a"<<(i+1)<<";\n";
            else
                atomsForCombination<<"int p"<<(i+1)<<" = p1+a"<<(i+1)<<";\n";
        }
        else {
            if (centralParticleMode)
                atomsForCombination<<"int p"<<(i+1)<<" = neighbors[firstNeighbor+a"<<(i+1)<<"];\n";
            else
                atomsForCombination<<"int p"<<(i+1)<<" = neighbors[firstNeighbor-1+a"<<(i+1)<<"];\n";
        }
    }
    if (nonbondedMethod != NoCutoff) {
        for (int i = 1; i < particlesPerSet; i++)
            verifyCutoff<<"real3 pos"<<(i+1)<<" = trimTo3(posq[p"<<(i+1)<<"]);\n";
        if (!centralParticleMode) {
            for (int i = 1; i < particlesPerSet; i++) {
                for (int j = i+1; j < particlesPerSet; j++)
                    verifyCutoff<<"includeInteraction &= (delta(pos"<<(i+1)<<", pos"<<(j+1)<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ).w < CUTOFF_SQUARED);\n";
            }
        }
    }
    if (force.getNumExclusions() > 0) {
        int startCheckFrom = (nonbondedMethod == NoCutoff ? 0 : 1);
        for (int i = startCheckFrom; i < particlesPerSet; i++)
            for (int j = i+1; j < particlesPerSet; j++)
                verifyExclusions<<"includeInteraction &= !isInteractionExcluded(p"<<(i+1)<<", p"<<(j+1)<<", exclusions, exclusionStartIndex);\n";
    }
    string computeTypeIndex = "particleTypes[p"+cc.intToString(particlesPerSet)+"]";
    for (int i = particlesPerSet-2; i >= 0; i--)
        computeTypeIndex = "particleTypes[p"+cc.intToString(i+1)+"]+"+cc.intToString(numTypes)+"*("+computeTypeIndex+")";

    // Create replacements for extra arguments.

    stringstream extraArgs;
    if (force.getNumGlobalParameters() > 0)
        extraArgs << ", GLOBAL const float* globals";
    for (int i = 0; i < (int) params->getParameterInfos().size(); i++) {
        ComputeParameterInfo& parameter = params->getParameterInfos()[i];
        extraArgs<<", GLOBAL const "<<parameter.getType()<<"* RESTRICT global_params"<<(i+1);
    }

    // Create the kernels.

    map<string, string> replacements;
    replacements["COMPUTE_INTERACTION"] = compute.str();
    replacements["NUM_CANDIDATE_COMBINATIONS"] = numCombinations.str();
    replacements["FIND_ATOMS_FOR_COMBINATION_INDEX"] = atomsForCombination.str();
    replacements["IS_VALID_COMBINATION"] = isValidCombination.str();
    replacements["VERIFY_CUTOFF"] = verifyCutoff.str();
    replacements["VERIFY_EXCLUSIONS"] = verifyExclusions.str();
    replacements["PERMUTE_ATOMS"] = permute.str();
    replacements["LOAD_PARTICLE_DATA"] = loadData.str();
    replacements["COMPUTE_TYPE_INDEX"] = computeTypeIndex;
    replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
    map<string, string> defines;
    if (nonbondedMethod != NoCutoff)
        defines["USE_CUTOFF"] = "1";
    if (nonbondedMethod == CutoffPeriodic)
        defines["USE_PERIODIC"] = "1";
    if (centralParticleMode)
        defines["USE_CENTRAL_PARTICLE"] = "1";
    if (hasTypeFilters)
        defines["USE_FILTERS"] = "1";
    if (force.getNumExclusions() > 0)
        defines["USE_EXCLUSIONS"] = "1";
    defines["NUM_ATOMS"] = cc.intToString(cc.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cc.intToString(cc.getPaddedNumAtoms());
    defines["M_PI"] = cc.doubleToString(M_PI);
    defines["CUTOFF_SQUARED"] = cc.doubleToString(force.getCutoffDistance()*force.getCutoffDistance());
    defines["TILE_SIZE"] = cc.intToString(32);
    defines["NUM_BLOCKS"] = cc.intToString(numAtomBlocks);
    defines["FIND_NEIGHBORS_WORKGROUP_SIZE"] = cc.intToString(findNeighborsWorkgroupSize);
    ComputeProgram program = cc.compileProgram(cc.replaceStrings(CommonKernelSources::pointFunctions+CommonKernelSources::customManyParticle, replacements), defines);
    forceKernel = program->createKernel("computeInteraction");
    blockBoundsKernel = program->createKernel("findBlockBounds");
    neighborsKernel = program->createKernel("findNeighbors");
    startIndicesKernel = program->createKernel("computeNeighborStartIndices");
    copyPairsKernel = program->createKernel("copyPairsToNeighborList");
    event = cc.createEvent();
}

double CommonCalcCustomManyParticleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    ContextSelector selector(cc);
    if (!hasInitializedKernel) {
        hasInitializedKernel = true;

        // Set arguments for the force kernel.

        forceKernel->addArg(cc.getLongForceBuffer());
        forceKernel->addArg(cc.getEnergyBuffer());
        forceKernel->addArg(cc.getPosq());
        for (int i = 0; i < 5; i++)
            forceKernel->addArg();
        setPeriodicBoxArgs(cc, forceKernel, 3);
        if (nonbondedMethod != NoCutoff) {
            forceKernel->addArg(neighbors);
            forceKernel->addArg(neighborStartIndex);
        }
        if (particleTypes.isInitialized()) {
            forceKernel->addArg(particleTypes);
            forceKernel->addArg(orderIndex);
            forceKernel->addArg(particleOrder);
        }
        if (exclusions.isInitialized()) {
            forceKernel->addArg(exclusions);
            forceKernel->addArg(exclusionStartIndex);
        }
        if (globals.isInitialized())
            forceKernel->addArg(globals);
        for (auto& parameter : params->getParameterInfos())
            forceKernel->addArg(parameter.getArray());
        for (auto& function : tabulatedFunctionArrays)
            forceKernel->addArg(function);

        if (nonbondedMethod != NoCutoff) {
            // Set arguments for the block bounds kernel.

            for (int i = 0; i < 5; i++)
                blockBoundsKernel->addArg(); // Periodic box information will be set just before it is executed.
            blockBoundsKernel->addArg(cc.getPosq());
            blockBoundsKernel->addArg(blockCenter);
            blockBoundsKernel->addArg(blockBoundingBox);
            blockBoundsKernel->addArg(numNeighborPairs);

            // Set arguments for the neighbor list kernel.

            for (int i = 0; i < 5; i++)
                neighborsKernel->addArg(); // Periodic box information will be set just before it is executed.
            neighborsKernel->addArg(cc.getPosq());
            neighborsKernel->addArg(blockCenter);
            neighborsKernel->addArg(blockBoundingBox);
            neighborsKernel->addArg(neighborPairs);
            neighborsKernel->addArg(numNeighborPairs);
            neighborsKernel->addArg(numNeighborsForAtom);
            neighborsKernel->addArg(maxNeighborPairs);
            if (exclusions.isInitialized()) {
                neighborsKernel->addArg(exclusions);
                neighborsKernel->addArg(exclusionStartIndex);
            }

            // Set arguments for the kernel to find neighbor list start indices.

            startIndicesKernel->addArg(numNeighborsForAtom);
            startIndicesKernel->addArg(neighborStartIndex);
            startIndicesKernel->addArg(numNeighborPairs);
            startIndicesKernel->addArg(maxNeighborPairs);

            // Set arguments for the kernel to assemble the final neighbor list.

            copyPairsKernel->addArg(neighborPairs);
            copyPairsKernel->addArg(neighbors);
            copyPairsKernel->addArg(numNeighborPairs);
            copyPairsKernel->addArg(maxNeighborPairs);
            copyPairsKernel->addArg(numNeighborsForAtom);
            copyPairsKernel->addArg(neighborStartIndex);
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
    while (true) {
        int* numPairs = (int*) cc.getPinnedBuffer();
        if (nonbondedMethod != NoCutoff) {
            setPeriodicBoxArgs(cc, forceKernel, 3);
            setPeriodicBoxArgs(cc, blockBoundsKernel, 0);
            setPeriodicBoxArgs(cc, neighborsKernel, 0);
            blockBoundsKernel->execute(cc.getPaddedNumAtoms()/32);
            neighborsKernel->execute(cc.getNumAtoms(), findNeighborsWorkgroupSize);

            // We need to make sure there was enough memory for the neighbor list.  Download the
            // information asynchronously so kernels can be running at the same time.

            numNeighborPairs.download(numPairs, false);
            event->enqueue();
            startIndicesKernel->execute(256, 256);
            copyPairsKernel->execute(maxNeighborPairs);
        }
        int maxThreads = min(cc.getNumAtoms()*forceWorkgroupSize, (int) cc.getEnergyBuffer().getSize());
        forceKernel->execute(maxThreads, forceWorkgroupSize);
        if (nonbondedMethod != NoCutoff) {
            // Make sure there was enough memory for the neighbor list.

            event->wait();
            if (*numPairs > maxNeighborPairs) {
                // Resize the arrays and run the calculation again.

                maxNeighborPairs = (int) (1.1*(*numPairs));
                neighborPairs.resize(maxNeighborPairs);
                neighbors.resize(maxNeighborPairs);
                neighborsKernel->setArg(11, maxNeighborPairs);
                startIndicesKernel->setArg(3, maxNeighborPairs);
                copyPairsKernel->setArg(3, maxNeighborPairs);
                continue;
            }
        }
        break;
    }
    return 0.0;
}

void CommonCalcCustomManyParticleForceKernel::copyParametersToContext(ContextImpl& context, const CustomManyParticleForce& force) {
    ContextSelector selector(cc);
    int numParticles = force.getNumParticles();
    if (numParticles != cc.getNumAtoms())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");

    // Record the per-particle parameters.

    vector<vector<float> > paramVector(numParticles);
    vector<double> parameters;
    int type;
    for (int i = 0; i < numParticles; i++) {
        force.getParticleParameters(i, parameters, type);
        paramVector[i].resize(parameters.size());
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

    cc.invalidateMolecules(info);
}
