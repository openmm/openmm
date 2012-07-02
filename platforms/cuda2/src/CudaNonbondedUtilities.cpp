/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2012 Stanford University and the Authors.      *
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

#include "openmm/OpenMMException.h"
#include "CudaNonbondedUtilities.h"
#include "CudaArray.h"
#include "CudaKernelSources.h"
#include "CudaExpressionUtilities.h"
#include <map>
#include <set>
#include <utility>

using namespace OpenMM;
using namespace std;

#define CHECK_RESULT(result) \
    if (result != CUDA_SUCCESS) { \
        std::stringstream m; \
        m<<errorMessage<<": "<<context.getErrorString(result)<<" ("<<result<<")"<<" at "<<__FILE__<<":"<<__LINE__; \
        throw OpenMMException(m.str());\
    }

CudaNonbondedUtilities::CudaNonbondedUtilities(CudaContext& context) : context(context), cutoff(-1.0), useCutoff(false), anyExclusions(false),
        exclusionIndices(NULL), exclusionRowIndices(NULL), exclusions(NULL), interactingTiles(NULL), interactionFlags(NULL),
        interactionCount(NULL), blockCenter(NULL), blockBoundingBox(NULL), pinnedInteractionCount(NULL), nonbondedForceGroup(0) {
    // Decide how many thread blocks to use.

    string errorMessage = "Error initializing nonbonded utilities";
    int multiprocessors;
    CHECK_RESULT(cuDeviceGetAttribute(&multiprocessors, CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT, context.getDevice()));
    numForceThreadBlocks = 3*multiprocessors;
    forceThreadBlockSize = (context.getComputeCapability() < 2.0 ? 128 : 256);
}

CudaNonbondedUtilities::~CudaNonbondedUtilities() {
    if (exclusionIndices != NULL)
        delete exclusionIndices;
    if (exclusionRowIndices != NULL)
        delete exclusionRowIndices;
    if (exclusions != NULL)
        delete exclusions;
    if (interactingTiles != NULL)
        delete interactingTiles;
    if (interactionFlags != NULL)
        delete interactionFlags;
    if (interactionCount != NULL)
        delete interactionCount;
    if (blockCenter != NULL)
        delete blockCenter;
    if (blockBoundingBox != NULL)
        delete blockBoundingBox;
    if (pinnedInteractionCount != NULL)
        cuMemFreeHost(pinnedInteractionCount);
}

void CudaNonbondedUtilities::addInteraction(bool usesCutoff, bool usesPeriodic, bool usesExclusions, double cutoffDistance, const vector<vector<int> >& exclusionList, const string& kernel, int forceGroup) {
    if (cutoff != -1.0) {
        if (usesCutoff != useCutoff)
            throw OpenMMException("All Forces must agree on whether to use a cutoff");
        if (usesPeriodic != usePeriodic)
            throw OpenMMException("All Forces must agree on whether to use periodic boundary conditions");
        if (cutoffDistance != cutoff)
            throw OpenMMException("All Forces must use the same cutoff distance");
        if (forceGroup != nonbondedForceGroup)
            throw OpenMMException("All nonbonded forces must be in the same force group");
    }
    if (usesExclusions)
        requestExclusions(exclusionList);
    useCutoff = usesCutoff;
    usePeriodic = usesPeriodic;
    cutoff = cutoffDistance;
    kernelSource += kernel+"\n";
    nonbondedForceGroup = forceGroup;
}

void CudaNonbondedUtilities::addParameter(const ParameterInfo& parameter) {
    parameters.push_back(parameter);
}

void CudaNonbondedUtilities::addArgument(const ParameterInfo& parameter) {
    arguments.push_back(parameter);
}

void CudaNonbondedUtilities::requestExclusions(const vector<vector<int> >& exclusionList) {
    if (anyExclusions) {
        bool sameExclusions = (exclusionList.size() == atomExclusions.size());
        for (int i = 0; i < (int) exclusionList.size() && sameExclusions; i++) {
            if (exclusionList[i].size() != atomExclusions[i].size())
                sameExclusions = false;
            for (int j = 0; j < (int) exclusionList[i].size(); j++)
                if (exclusionList[i][j] != atomExclusions[i][j])
                    sameExclusions = false;
        }
        if (!sameExclusions)
            throw OpenMMException("All Forces must have identical exceptions");
    }
    else {
        atomExclusions = exclusionList;
        anyExclusions = true;
    }
}

void CudaNonbondedUtilities::initialize(const System& system) {
    if (cutoff == -1.0)
        return; // There are no nonbonded interactions in the System.
    string errorMessage = "Error initializing nonbonded utilities";    
    if (atomExclusions.size() == 0) {
        // No exclusions were specifically requested, so just mark every atom as not interacting with itself.
        
        atomExclusions.resize(context.getNumAtoms());
        for (int i = 0; i < (int) atomExclusions.size(); i++)
            atomExclusions[i].push_back(i);
    }

    // Create the list of tiles.

    numAtoms = context.getNumAtoms();
    int numAtomBlocks = context.getNumAtomBlocks();
    int totalTiles = numAtomBlocks*(numAtomBlocks+1)/2;
    int numContexts = context.getPlatformData().contexts.size();
    startTileIndex = context.getContextIndex()*totalTiles/numContexts;
    int endTileIndex = (context.getContextIndex()+1)*totalTiles/numContexts;
    numTiles = endTileIndex-startTileIndex;

    // Build a list of indices for the tiles with exclusions.

    set<pair<int, int> > tilesWithExclusions;
    for (int atom1 = 0; atom1 < (int) atomExclusions.size(); ++atom1) {
        int x = atom1/CudaContext::TileSize;
        for (int j = 0; j < (int) atomExclusions[atom1].size(); ++j) {
            int atom2 = atomExclusions[atom1][j];
            int y = atom2/CudaContext::TileSize;
            tilesWithExclusions.insert(make_pair(max(x, y), min(x, y)));
        }
    }
    if (context.getPaddedNumAtoms() > context.getNumAtoms()) {
        for (int i = 0; i < numAtomBlocks; ++i)
            tilesWithExclusions.insert(make_pair(numAtomBlocks-1, i));
    }
    vector<unsigned int> exclusionRowIndicesVec(numAtomBlocks+1, 0);
    vector<unsigned int> exclusionIndicesVec;
    int currentRow = 0;
    for (set<pair<int, int> >::const_iterator iter = tilesWithExclusions.begin(); iter != tilesWithExclusions.end(); ++iter) {
        while (iter->first != currentRow)
            exclusionRowIndicesVec[++currentRow] = exclusionIndicesVec.size();
        exclusionIndicesVec.push_back(iter->second);
    }
    exclusionRowIndicesVec[++currentRow] = exclusionIndicesVec.size();
    exclusionIndices = CudaArray::create<unsigned int>(context, exclusionIndicesVec.size(), "exclusionIndices");
    exclusionRowIndices = CudaArray::create<unsigned int>(context, exclusionRowIndicesVec.size(), "exclusionRowIndices");
    exclusionIndices->upload(exclusionIndicesVec);
    exclusionRowIndices->upload(exclusionRowIndicesVec);

    // Record the exclusion data.

    exclusions = CudaArray::create<unsigned int>(context, tilesWithExclusions.size()*CudaContext::TileSize, "exclusions");
    vector<unsigned int> exclusionVec(exclusions->getSize());
    for (int i = 0; i < exclusions->getSize(); ++i)
        exclusionVec[i] = 0xFFFFFFFF;
    for (int atom1 = 0; atom1 < (int) atomExclusions.size(); ++atom1) {
        int x = atom1/CudaContext::TileSize;
        int offset1 = atom1-x*CudaContext::TileSize;
        for (int j = 0; j < (int) atomExclusions[atom1].size(); ++j) {
            int atom2 = atomExclusions[atom1][j];
            int y = atom2/CudaContext::TileSize;
            int offset2 = atom2-y*CudaContext::TileSize;
            if (x > y) {
                int index = findExclusionIndex(x, y, exclusionIndicesVec, exclusionRowIndicesVec);
                exclusionVec[index+offset1] &= 0xFFFFFFFF-(1<<offset2);
            }
            else {
                int index = findExclusionIndex(y, x, exclusionIndicesVec, exclusionRowIndicesVec);
                exclusionVec[index+offset2] &= 0xFFFFFFFF-(1<<offset1);
            }
        }
    }

    // Mark all interactions that involve a padding atom as being excluded.

    for (int atom1 = context.getNumAtoms(); atom1 < context.getPaddedNumAtoms(); ++atom1) {
        int x = atom1/CudaContext::TileSize;
        int offset1 = atom1-x*CudaContext::TileSize;
        for (int atom2 = 0; atom2 < context.getPaddedNumAtoms(); ++atom2) {
            int y = atom2/CudaContext::TileSize;
            int offset2 = atom2-y*CudaContext::TileSize;
            if (x >= y) {
                int index = findExclusionIndex(x, y, exclusionIndicesVec, exclusionRowIndicesVec);
                exclusionVec[index+offset1] &= 0xFFFFFFFF-(1<<offset2);
            }
            if (y >= x) {
                int index = findExclusionIndex(y, x, exclusionIndicesVec, exclusionRowIndicesVec);
                exclusionVec[index+offset2] &= 0xFFFFFFFF-(1<<offset1);
            }
        }
    }
    atomExclusions.clear(); // We won't use this again, so free the memory it used
    exclusions->upload(exclusionVec);

    // Create data structures for the neighbor list.

    if (useCutoff) {
        // Select a size for the arrays that hold the neighbor list.  This estimate is intentionally very
        // high, because if it ever is too small, we have to fall back to the N^2 algorithm.

        double4 boxSize = context.getPeriodicBoxSize();
        maxTiles = (int) (numTiles*(cutoff/boxSize.x+cutoff/boxSize.y+cutoff/boxSize.z));
        if (maxTiles > numTiles)
            maxTiles = numTiles;
        if (maxTiles < 1)
            maxTiles = 1;
        interactingTiles = CudaArray::create<ushort2>(context, maxTiles, "interactingTiles");
        interactionFlags = CudaArray::create<unsigned int>(context, maxTiles, "interactionFlags");
        interactionCount = CudaArray::create<unsigned int>(context, 1, "interactionCount");
        blockCenter = CudaArray::create<float4>(context, numAtomBlocks, "blockCenter");
        blockBoundingBox = CudaArray::create<float4>(context, numAtomBlocks, "blockBoundingBox");
        CHECK_RESULT(cuMemHostAlloc((void**) &pinnedInteractionCount, sizeof(unsigned int), 0));
        pinnedInteractionCount[0] = 0;
        interactionCount->upload(pinnedInteractionCount);
    }

    // Create kernels.

    forceKernel = createInteractionKernel(kernelSource, parameters, arguments, true, true);
    if (useCutoff) {
        map<string, string> defines;
        defines["NUM_BLOCKS"] = context.intToString(context.getNumAtomBlocks());
        defines["CUTOFF_SQUARED"] = context.doubleToString(cutoff*cutoff);
        if (usePeriodic)
            defines["USE_PERIODIC"] = "1";
        CUmodule interactingBlocksProgram = context.createModule(CudaKernelSources::vectorOps+CudaKernelSources::findInteractingBlocks, defines);
        findBlockBoundsKernel = context.getKernel(interactingBlocksProgram, "findBlockBounds");
        findBlockBoundsArgs.push_back(&numAtoms);
        findBlockBoundsArgs.push_back(context.getPeriodicBoxSizePointer());
        findBlockBoundsArgs.push_back(context.getInvPeriodicBoxSizePointer());
        findBlockBoundsArgs.push_back(&context.getPosq().getDevicePointer());
        findBlockBoundsArgs.push_back(&blockCenter->getDevicePointer());
        findBlockBoundsArgs.push_back(&blockBoundingBox->getDevicePointer());
        findBlockBoundsArgs.push_back(&interactionCount->getDevicePointer());
        findInteractingBlocksKernel = context.getKernel(interactingBlocksProgram, "findBlocksWithInteractions");
        findInteractingBlocksArgs.push_back(context.getPeriodicBoxSizePointer());
        findInteractingBlocksArgs.push_back(context.getInvPeriodicBoxSizePointer());
        findInteractingBlocksArgs.push_back(&blockCenter->getDevicePointer());
        findInteractingBlocksArgs.push_back(&blockBoundingBox->getDevicePointer());
        findInteractingBlocksArgs.push_back(&interactionCount->getDevicePointer());
        findInteractingBlocksArgs.push_back(&interactingTiles->getDevicePointer());
        findInteractingBlocksArgs.push_back(&interactionFlags->getDevicePointer());
        findInteractingBlocksArgs.push_back(&context.getPosq().getDevicePointer());
        findInteractingBlocksArgs.push_back(&maxTiles);
        findInteractingBlocksArgs.push_back(&startTileIndex);
        findInteractingBlocksArgs.push_back(&numTiles);
        findInteractionsWithinBlocksKernel = context.getKernel(interactingBlocksProgram, "findInteractionsWithinBlocks");
        findInteractionsWithinBlocksArgs.push_back(context.getPeriodicBoxSizePointer());
        findInteractionsWithinBlocksArgs.push_back(context.getInvPeriodicBoxSizePointer());
        findInteractionsWithinBlocksArgs.push_back(&context.getPosq().getDevicePointer());
        findInteractionsWithinBlocksArgs.push_back(&interactingTiles->getDevicePointer());
        findInteractionsWithinBlocksArgs.push_back(&blockCenter->getDevicePointer());
        findInteractionsWithinBlocksArgs.push_back(&blockBoundingBox->getDevicePointer());
        findInteractionsWithinBlocksArgs.push_back(&interactionFlags->getDevicePointer());
        findInteractionsWithinBlocksArgs.push_back(&interactionCount->getDevicePointer());
        findInteractionsWithinBlocksArgs.push_back(&maxTiles);
    }
}

int CudaNonbondedUtilities::findExclusionIndex(int x, int y, const vector<unsigned int>& exclusionIndices, const vector<unsigned int>& exclusionRowIndices) {
    int start = exclusionRowIndices[x];
    int end = exclusionRowIndices[x+1];
    for (int i = start; i < end; i++)
        if (exclusionIndices[i] == y)
            return i*CudaContext::TileSize;
    throw OpenMMException("Internal error: exclusion in unexpected tile");
}

void CudaNonbondedUtilities::prepareInteractions() {
    if (!useCutoff)
        return;
    if (usePeriodic) {
        double4 box = context.getPeriodicBoxSize();
        double minAllowedSize = 1.999999*cutoff;
        if (box.x < minAllowedSize || box.y < minAllowedSize || box.z < minAllowedSize)
            throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
    }

    // Compute the neighbor list.

    context.executeKernel(findBlockBoundsKernel, &findBlockBoundsArgs[0], context.getNumAtoms());
    context.executeKernel(findInteractingBlocksKernel, &findInteractingBlocksArgs[0], context.getNumAtoms());
    context.executeKernel(findInteractionsWithinBlocksKernel, &findInteractionsWithinBlocksArgs[0], context.getNumAtoms(), 128);
}

void CudaNonbondedUtilities::computeInteractions() {
    if (cutoff != -1.0)
        context.executeKernel(forceKernel, &forceArgs[0], numForceThreadBlocks*forceThreadBlockSize, forceThreadBlockSize);
}

void CudaNonbondedUtilities::updateNeighborListSize() {
    if (!useCutoff)
        return;
    interactionCount->download(pinnedInteractionCount);
    if (pinnedInteractionCount[0] <= (unsigned int) maxTiles)
        return;

    // The most recent timestep had too many interactions to fit in the arrays.  Make the arrays bigger to prevent
    // this from happening in the future.

    maxTiles = (int) (1.2*pinnedInteractionCount[0]);
    int numTiles = context.getNumAtomBlocks()*(context.getNumAtomBlocks()+1)/2;
    if (maxTiles > numTiles)
        maxTiles = numTiles;
    delete interactingTiles;
    interactingTiles = CudaArray::create<ushort2>(context, maxTiles, "interactingTiles");
    forceArgs[8] = &interactingTiles->getDevicePointer();
    findInteractingBlocksArgs[5] = &interactingTiles->getDevicePointer();
    delete interactionFlags;
    interactionFlags = CudaArray::create<unsigned int>(context, maxTiles, "interactionFlags");
    forceArgs[13] = &interactionFlags->getDevicePointer();
    findInteractingBlocksArgs[6] = &interactionFlags->getDevicePointer();
    findInteractionsWithinBlocksArgs[3] = &interactingTiles->getDevicePointer();
    findInteractionsWithinBlocksArgs[6] = &interactionFlags->getDevicePointer();
}

void CudaNonbondedUtilities::setTileRange(int startTileIndex, int numTiles) {
    this->startTileIndex = startTileIndex;
    this->numTiles = numTiles;
}

CUfunction CudaNonbondedUtilities::createInteractionKernel(const string& source, vector<ParameterInfo>& params, vector<ParameterInfo>& arguments, bool useExclusions, bool isSymmetric) {
    map<string, string> replacements;
    replacements["COMPUTE_INTERACTION"] = source;
    const string suffixes[] = {"x", "y", "z", "w"};
    stringstream localData;
    int localDataSize = 0;
    for (int i = 0; i < (int) params.size(); i++) {
        if (params[i].getNumComponents() == 1)
            localData<<params[i].getType()<<" "<<params[i].getName()<<";\n";
        else {
            for (int j = 0; j < params[i].getNumComponents(); ++j)
                localData<<params[i].getComponentType()<<" "<<params[i].getName()<<"_"<<suffixes[j]<<";\n";
        }
        localDataSize += params[i].getSize();
    }
    replacements["ATOM_PARAMETER_DATA"] = localData.str();
    stringstream args;
    for (int i = 0; i < (int) params.size(); i++) {
        args << ", const ";
        args << params[i].getType();
        args << "* __restrict__ global_";
        args << params[i].getName();
    }
    for (int i = 0; i < (int) arguments.size(); i++) {
        args << ", const ";
        args << arguments[i].getType();
        args << "* __restrict__ ";
        args << arguments[i].getName();
    }
    replacements["PARAMETER_ARGUMENTS"] = args.str();
    stringstream loadLocal1;
    for (int i = 0; i < (int) params.size(); i++) {
        if (params[i].getNumComponents() == 1) {
            loadLocal1<<"localData[localAtomIndex]."<<params[i].getName()<<" = "<<params[i].getName()<<"1;\n";
        }
        else {
            for (int j = 0; j < params[i].getNumComponents(); ++j)
                loadLocal1<<"localData[localAtomIndex]."<<params[i].getName()<<"_"<<suffixes[j]<<" = "<<params[i].getName()<<"1."<<suffixes[j]<<";\n";
        }
    }
    replacements["LOAD_LOCAL_PARAMETERS_FROM_1"] = loadLocal1.str();
    stringstream loadLocal2;
    for (int i = 0; i < (int) params.size(); i++) {
        if (params[i].getNumComponents() == 1) {
            loadLocal2<<"localData[localAtomIndex]."<<params[i].getName()<<" = global_"<<params[i].getName()<<"[j];\n";
        }
        else {
            loadLocal2<<params[i].getType()<<" temp_"<<params[i].getName()<<" = global_"<<params[i].getName()<<"[j];\n";
            for (int j = 0; j < params[i].getNumComponents(); ++j)
                loadLocal2<<"localData[localAtomIndex]."<<params[i].getName()<<"_"<<suffixes[j]<<" = temp_"<<params[i].getName()<<"."<<suffixes[j]<<";\n";
        }
    }
    replacements["LOAD_LOCAL_PARAMETERS_FROM_GLOBAL"] = loadLocal2.str();
    stringstream load1;
    for (int i = 0; i < (int) params.size(); i++) {
        load1 << params[i].getType();
        load1 << " ";
        load1 << params[i].getName();
        load1 << "1 = global_";
        load1 << params[i].getName();
        load1 << "[atom1];\n";
    }
    replacements["LOAD_ATOM1_PARAMETERS"] = load1.str();
    stringstream load2j;
    for (int i = 0; i < (int) params.size(); i++) {
        if (params[i].getNumComponents() == 1) {
            load2j<<params[i].getType()<<" "<<params[i].getName()<<"2 = localData[atom2]."<<params[i].getName()<<";\n";
        }
        else {
            load2j<<params[i].getType()<<" "<<params[i].getName()<<"2 = make_"<<params[i].getType()<<"(";
            for (int j = 0; j < params[i].getNumComponents(); ++j) {
                if (j > 0)
                    load2j<<", ";
                load2j<<"localData[atom2]."<<params[i].getName()<<"_"<<suffixes[j];
            }
            load2j<<");\n";
        }
    }
    replacements["LOAD_ATOM2_PARAMETERS"] = load2j.str();
    map<string, string> defines;
    if (useCutoff)
        defines["USE_CUTOFF"] = "1";
    if (usePeriodic)
        defines["USE_PERIODIC"] = "1";
    if (useExclusions)
        defines["USE_EXCLUSIONS"] = "1";
    if (isSymmetric)
        defines["USE_SYMMETRIC"] = "1";
    defines["THREAD_BLOCK_SIZE"] = context.intToString(forceThreadBlockSize);
    defines["CUTOFF_SQUARED"] = context.doubleToString(cutoff*cutoff);
    defines["NUM_ATOMS"] = context.intToString(context.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = context.intToString(context.getPaddedNumAtoms());
    defines["NUM_BLOCKS"] = context.intToString(context.getNumAtomBlocks());
    if ((localDataSize/4)%2 == 0 && !context.getUseDoublePrecision())
        defines["PARAMETER_SIZE_IS_EVEN"] = "1";
    CUmodule program = context.createModule(CudaKernelSources::vectorOps+context.replaceStrings(CudaKernelSources::nonbonded, replacements), defines);
    CUfunction kernel = context.getKernel(program, "computeNonbonded");

    // Set arguments to the Kernel.

    int index = 0;
    forceArgs.push_back(&context.getForce().getDevicePointer());
    forceArgs.push_back(&context.getEnergyBuffer().getDevicePointer());
    forceArgs.push_back(&context.getPosq().getDevicePointer());
    forceArgs.push_back(&exclusions->getDevicePointer());
    forceArgs.push_back(&exclusionIndices->getDevicePointer());
    forceArgs.push_back(&exclusionRowIndices->getDevicePointer());
    forceArgs.push_back(&startTileIndex);
    forceArgs.push_back(&numTiles);
    if (useCutoff) {
        forceArgs.push_back(&interactingTiles->getDevicePointer());
        forceArgs.push_back(&interactionCount->getDevicePointer());
        forceArgs.push_back(context.getPeriodicBoxSizePointer());
        forceArgs.push_back(context.getInvPeriodicBoxSizePointer());
        forceArgs.push_back(&maxTiles);
        forceArgs.push_back(&interactionFlags->getDevicePointer());
    }
    for (int i = 0; i < (int) params.size(); i++)
        forceArgs.push_back(&params[i].getMemory());
    for (int i = 0; i < (int) arguments.size(); i++)
        forceArgs.push_back(&arguments[i].getMemory());
    return kernel;
}
