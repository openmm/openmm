/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2023 Stanford University and the Authors.      *
 * Portions copyright (c) 2020-2023 Advanced Micro Devices, Inc.              *
 * Authors: Peter Eastman, Nicholas Curtis                                    *
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
#include "HipNonbondedUtilities.h"
#include "HipArray.h"
#include "HipContext.h"
#include "HipKernelSources.h"
#include "HipExpressionUtilities.h"
#include "HipSort.h"
#include <algorithm>
#include <map>
#include <set>
#include <utility>

using namespace OpenMM;
using namespace std;

#define CHECK_RESULT(result) \
    if (result != hipSuccess) { \
        std::stringstream m; \
        m<<errorMessage<<": "<<context.getErrorString(result)<<" ("<<result<<")"<<" at "<<__FILE__<<":"<<__LINE__; \
        throw OpenMMException(m.str());\
    }


class HipNonbondedUtilities::BlockSortTrait : public HipSort::SortTrait {
public:
    BlockSortTrait() {}
    int getDataSize() const {return sizeof(int);}
    int getKeySize() const {return sizeof(int);}
    const char* getDataType() const {return "unsigned int";}
    const char* getKeyType() const {return "unsigned int";}
    const char* getMinKey() const {return "0";}
    const char* getMaxKey() const {return "0xFFFFFFFFu";}
    const char* getMaxValue() const {return "0xFFFFFFFFu";}
    const char* getSortKey() const {return "value";}
};

HipNonbondedUtilities::HipNonbondedUtilities(HipContext& context) : context(context), useCutoff(false), usePeriodic(false), useNeighborList(false), anyExclusions(false), usePadding(true),
        blockSorter(NULL), pinnedCountBuffer(NULL), forceRebuildNeighborList(true), lastCutoff(0.0), groupFlags(0), canUsePairList(true), tilesAfterReorder(0) {
    // Decide how many thread blocks to use.

    string errorMessage = "Error initializing nonbonded utilities";
    CHECK_RESULT(hipEventCreateWithFlags(&downloadCountEvent, context.getEventFlags()));
    CHECK_RESULT(hipHostMalloc((void**) &pinnedCountBuffer, 2*sizeof(unsigned int), context.getHostMallocFlags()));
    numForceThreadBlocks = 5*4*context.getMultiprocessors();
    forceThreadBlockSize = 64;
    findInteractingBlocksThreadBlockSize = context.getSIMDWidth();

    // When building the neighbor list, we can optionally use large blocks (32 * warpSize atoms) to
    // accelerate the process.  This makes building the neighbor list faster, but it prevents
    // us from sorting atom blocks by size, which leads to a slightly less efficient neighbor
    // list.  We guess based on system size which will be faster.

    useLargeBlocks = (context.getNumAtoms() > 90000);
    setKernelSource(HipKernelSources::nonbonded);
}

HipNonbondedUtilities::~HipNonbondedUtilities() {
    if (blockSorter != NULL)
        delete blockSorter;
    if (pinnedCountBuffer != NULL)
        hipHostFree(pinnedCountBuffer);
    hipEventDestroy(downloadCountEvent);
}

void HipNonbondedUtilities::addInteraction(bool usesCutoff, bool usesPeriodic, bool usesExclusions, double cutoffDistance, const vector<vector<int> >& exclusionList, const string& kernel, int forceGroup, bool usesNeighborList) {
    addInteraction(usesCutoff, usesPeriodic, usesExclusions, cutoffDistance, exclusionList, kernel, forceGroup, usesNeighborList, false);
}

void HipNonbondedUtilities::addInteraction(bool usesCutoff, bool usesPeriodic, bool usesExclusions, double cutoffDistance, const vector<vector<int> >& exclusionList, const string& kernel, int forceGroup, bool usesNeighborList, bool supportsPairList) {
    if (groupCutoff.size() > 0) {
        if (usesCutoff != useCutoff)
            throw OpenMMException("All Forces must agree on whether to use a cutoff");
        if (usesPeriodic != usePeriodic)
            throw OpenMMException("All Forces must agree on whether to use periodic boundary conditions");
        if (usesCutoff && groupCutoff.find(forceGroup) != groupCutoff.end() && groupCutoff[forceGroup] != cutoffDistance)
            throw OpenMMException("All Forces in a single force group must use the same cutoff distance");
    }
    if (usesExclusions)
        requestExclusions(exclusionList);
    useCutoff = usesCutoff;
    usePeriodic = usesPeriodic;
    useNeighborList |= (usesNeighborList && useCutoff);
    groupCutoff[forceGroup] = cutoffDistance;
    groupFlags |= 1<<forceGroup;
    canUsePairList &= supportsPairList;
    if (kernel.size() > 0) {
        if (groupKernelSource.find(forceGroup) == groupKernelSource.end())
            groupKernelSource[forceGroup] = "";
        map<string, string> replacements;
        replacements["CUTOFF"] = "CUTOFF_"+context.intToString(forceGroup);
        replacements["CUTOFF_SQUARED"] = "CUTOFF_"+context.intToString(forceGroup)+"_SQUARED";
        groupKernelSource[forceGroup] += context.replaceStrings(kernel, replacements)+"\n";
    }
}

void HipNonbondedUtilities::addParameter(ComputeParameterInfo parameter) {
    parameters.push_back(ParameterInfo(parameter.getName(), parameter.getComponentType(), parameter.getNumComponents(),
            parameter.getSize(), context.unwrap(parameter.getArray()).getDevicePointer(), parameter.isConstant()));
}

void HipNonbondedUtilities::addParameter(const ParameterInfo& parameter) {
    parameters.push_back(parameter);
}

void HipNonbondedUtilities::addArgument(ComputeParameterInfo parameter) {
    arguments.push_back(ParameterInfo(parameter.getName(), parameter.getComponentType(), parameter.getNumComponents(),
            parameter.getSize(), context.unwrap(parameter.getArray()).getDevicePointer(), parameter.isConstant()));
}

void HipNonbondedUtilities::addArgument(const ParameterInfo& parameter) {
    arguments.push_back(parameter);
}

string HipNonbondedUtilities::addEnergyParameterDerivative(const string& param) {
    // See if the parameter has already been added.

    int index;
    for (index = 0; index < energyParameterDerivatives.size(); index++)
        if (param == energyParameterDerivatives[index])
            break;
    if (index == energyParameterDerivatives.size())
        energyParameterDerivatives.push_back(param);
    context.addEnergyParameterDerivative(param);
    return string("energyParamDeriv")+context.intToString(index);
}

void HipNonbondedUtilities::requestExclusions(const vector<vector<int> >& exclusionList) {
    if (anyExclusions) {
        bool sameExclusions = (exclusionList.size() == atomExclusions.size());
        for (int i = 0; i < (int) exclusionList.size() && sameExclusions; i++) {
             if (exclusionList[i].size() != atomExclusions[i].size())
                 sameExclusions = false;
            set<int> expectedExclusions;
            expectedExclusions.insert(atomExclusions[i].begin(), atomExclusions[i].end());
            for (int j = 0; j < (int) exclusionList[i].size(); j++)
                if (expectedExclusions.find(exclusionList[i][j]) == expectedExclusions.end())
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

static bool compareInt2(int2 a, int2 b) {
    return ((a.y < b.y) || (a.y == b.y && a.x < b.x));
}

static bool compareInt2LargeSIMD(int2 a, int2 b) {
    // This version is used on devices with SIMD width greater than tile size.  It puts diagonal tiles before off-diagonal
    // ones to reduce thread divergence.

    if (a.x == a.y) {
        if (b.x == b.y)
            return (a.x < b.x);
        return true;
    }
    if (b.x == b.y)
        return false;
    return ((a.y < b.y) || (a.y == b.y && a.x < b.x));
}

void HipNonbondedUtilities::initialize(const System& system) {
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
    int numContexts = context.getPlatformData().contexts.size();
    setAtomBlockRange(context.getContextIndex()/(double) numContexts, (context.getContextIndex()+1)/(double) numContexts);

    // Build a list of tiles that contain exclusions.

    set<pair<int, int> > tilesWithExclusions;
    for (int atom1 = 0; atom1 < (int) atomExclusions.size(); ++atom1) {
        int x = atom1/HipContext::TileSize;
        for (int j = 0; j < (int) atomExclusions[atom1].size(); ++j) {
            int atom2 = atomExclusions[atom1][j];
            int y = atom2/HipContext::TileSize;
            tilesWithExclusions.insert(make_pair(max(x, y), min(x, y)));
        }
    }
    vector<int2> exclusionTilesVec;
    for (set<pair<int, int> >::const_iterator iter = tilesWithExclusions.begin(); iter != tilesWithExclusions.end(); ++iter)
        exclusionTilesVec.push_back(make_int2(iter->first, iter->second));
    sort(exclusionTilesVec.begin(), exclusionTilesVec.end(), context.getSIMDWidth() <= 32 || !useNeighborList ? compareInt2 : compareInt2LargeSIMD);
    exclusionTiles.initialize<int2>(context, exclusionTilesVec.size(), "exclusionTiles");
    exclusionTiles.upload(exclusionTilesVec);
    map<pair<int, int>, int> exclusionTileMap;
    for (int i = 0; i < (int) exclusionTilesVec.size(); i++) {
        int2 tile = exclusionTilesVec[i];
        exclusionTileMap[make_pair(tile.x, tile.y)] = i;
    }
    vector<vector<int> > exclusionBlocksForBlock(numAtomBlocks);
    for (set<pair<int, int> >::const_iterator iter = tilesWithExclusions.begin(); iter != tilesWithExclusions.end(); ++iter) {
        exclusionBlocksForBlock[iter->first].push_back(iter->second);
        if (iter->first != iter->second)
            exclusionBlocksForBlock[iter->second].push_back(iter->first);
    }
    vector<unsigned int> exclusionRowIndicesVec(numAtomBlocks+1, 0);
    vector<unsigned int> exclusionIndicesVec;
    for (int i = 0; i < numAtomBlocks; i++) {
        exclusionIndicesVec.insert(exclusionIndicesVec.end(), exclusionBlocksForBlock[i].begin(), exclusionBlocksForBlock[i].end());
        exclusionRowIndicesVec[i+1] = exclusionIndicesVec.size();
    }
    maxExclusions = 0;
    for (int i = 0; i < (int) exclusionBlocksForBlock.size(); i++)
        maxExclusions = (maxExclusions > exclusionBlocksForBlock[i].size() ? maxExclusions : exclusionBlocksForBlock[i].size());
    exclusionIndices.initialize<unsigned int>(context, exclusionIndicesVec.size(), "exclusionIndices");
    exclusionRowIndices.initialize<unsigned int>(context, exclusionRowIndicesVec.size(), "exclusionRowIndices");
    exclusionIndices.upload(exclusionIndicesVec);
    exclusionRowIndices.upload(exclusionRowIndicesVec);

    // Record the exclusion data.

    exclusions.initialize<tileflags>(context, tilesWithExclusions.size()*HipContext::TileSize, "exclusions");
    tileflags allFlags = (tileflags) -1;
    vector<tileflags> exclusionVec(exclusions.getSize(), allFlags);
    for (int atom1 = 0; atom1 < (int) atomExclusions.size(); ++atom1) {
        int x = atom1/HipContext::TileSize;
        int offset1 = atom1-x*HipContext::TileSize;
        for (int j = 0; j < (int) atomExclusions[atom1].size(); ++j) {
            int atom2 = atomExclusions[atom1][j];
            int y = atom2/HipContext::TileSize;
            int offset2 = atom2-y*HipContext::TileSize;
            if (x > y) {
                int index = exclusionTileMap[make_pair(x, y)]*HipContext::TileSize;
                exclusionVec[index+offset1] &= allFlags-(1<<offset2);
            }
            else {
                int index = exclusionTileMap[make_pair(y, x)]*HipContext::TileSize;
                exclusionVec[index+offset2] &= allFlags-(1<<offset1);
            }
        }
    }
    atomExclusions.clear(); // We won't use this again, so free the memory it used
    exclusions.upload(exclusionVec);

    // Create data structures for the neighbor list.

    if (useCutoff) {
        // Select a size for the arrays that hold the neighbor list.  We have to make a fairly
        // arbitrary guess, but if this turns out to be too small we'll increase it later.

        maxTiles = 20*numAtomBlocks;
        if (maxTiles > numTiles)
            maxTiles = numTiles;
        if (maxTiles < 1)
            maxTiles = 1;
        maxSinglePairs = 5*numAtoms;
        // HIP-TODO: This may require tuning
        numTilesInBatch = numAtomBlocks < 2000 ? 4 : 1;
        interactingTiles.initialize<int>(context, maxTiles, "interactingTiles");
        interactingAtoms.initialize<int>(context, HipContext::TileSize*maxTiles, "interactingAtoms");
        interactionCount.initialize<unsigned int>(context, 2, "interactionCount");
        singlePairs.initialize<int2>(context, maxSinglePairs, "singlePairs");
        int elementSize = (context.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
        blockCenter.initialize(context, numAtomBlocks, 4*elementSize, "blockCenter");
        blockBoundingBox.initialize(context, numAtomBlocks, 4*elementSize, "blockBoundingBox");
        sortedBlocks.initialize<unsigned int>(context, numAtomBlocks, "sortedBlocks");
        sortedBlockCenter.initialize(context, numAtomBlocks+1, 4*elementSize, "sortedBlockCenter");
        sortedBlockBoundingBox.initialize(context, numAtomBlocks+1, 4*elementSize, "sortedBlockBoundingBox");
        blockSizeRange.initialize(context, 2, elementSize, "blockSizeRange");
        largeBlockCenter.initialize(context, numAtomBlocks, 4*elementSize, "largeBlockCenter");
        largeBlockBoundingBox.initialize(context, numAtomBlocks*4, elementSize, "largeBlockBoundingBox");
        oldPositions.initialize(context, numAtoms, 4*elementSize, "oldPositions");
        rebuildNeighborList.initialize<int>(context, 1, "rebuildNeighborList");
        blockSorter = new HipSort(context, new BlockSortTrait(), numAtomBlocks, false);
        vector<unsigned int> count(2, 0);
        interactionCount.upload(count);
        rebuildNeighborList.upload(&count[0]);
        if (context.getUseDoublePrecision()) {
            blockSizeRange.upload(vector<double>{1e38, 0});
        } else {
            blockSizeRange.upload(vector<float>{1e38, 0});
        }
    }

    // Record arguments for kernels.

    forceArgs.push_back(&context.getForce().getDevicePointer());
    forceArgs.push_back(&context.getEnergyBuffer().getDevicePointer());
    forceArgs.push_back(&context.getPosq().getDevicePointer());
    forceArgs.push_back(&exclusions.getDevicePointer());
    forceArgs.push_back(&exclusionTiles.getDevicePointer());
    forceArgs.push_back(&startTileIndex);
    forceArgs.push_back(&numTiles);
    if (useCutoff) {
        forceArgs.push_back(&interactingTiles.getDevicePointer());
        forceArgs.push_back(&interactionCount.getDevicePointer());
        forceArgs.push_back(context.getPeriodicBoxSizePointer());
        forceArgs.push_back(context.getInvPeriodicBoxSizePointer());
        forceArgs.push_back(context.getPeriodicBoxVecXPointer());
        forceArgs.push_back(context.getPeriodicBoxVecYPointer());
        forceArgs.push_back(context.getPeriodicBoxVecZPointer());
        forceArgs.push_back(&maxTiles);
        forceArgs.push_back(&blockCenter.getDevicePointer());
        forceArgs.push_back(&blockBoundingBox.getDevicePointer());
        forceArgs.push_back(&interactingAtoms.getDevicePointer());
        forceArgs.push_back(&maxSinglePairs);
        forceArgs.push_back(&singlePairs.getDevicePointer());
    }
    for (int i = 0; i < (int) parameters.size(); i++)
        forceArgs.push_back(&parameters[i].getMemory());
    for (ParameterInfo& arg : arguments)
        forceArgs.push_back(&arg.getMemory());
    if (energyParameterDerivatives.size() > 0)
        forceArgs.push_back(&context.getEnergyParamDerivBuffer().getDevicePointer());
    if (useCutoff) {
        findBlockBoundsArgs.push_back(&numAtoms);
        findBlockBoundsArgs.push_back(context.getPeriodicBoxSizePointer());
        findBlockBoundsArgs.push_back(context.getInvPeriodicBoxSizePointer());
        findBlockBoundsArgs.push_back(context.getPeriodicBoxVecXPointer());
        findBlockBoundsArgs.push_back(context.getPeriodicBoxVecYPointer());
        findBlockBoundsArgs.push_back(context.getPeriodicBoxVecZPointer());
        findBlockBoundsArgs.push_back(&context.getPosq().getDevicePointer());
        findBlockBoundsArgs.push_back(&blockCenter.getDevicePointer());
        findBlockBoundsArgs.push_back(&blockBoundingBox.getDevicePointer());
        findBlockBoundsArgs.push_back(&rebuildNeighborList.getDevicePointer());
        findBlockBoundsArgs.push_back(&blockSizeRange.getDevicePointer());
        computeSortKeysArgs.push_back(&blockBoundingBox.getDevicePointer());
        computeSortKeysArgs.push_back(&sortedBlocks.getDevicePointer());
        computeSortKeysArgs.push_back(&blockSizeRange.getDevicePointer());
        sortBoxDataArgs.push_back(&sortedBlocks.getDevicePointer());
        sortBoxDataArgs.push_back(&blockCenter.getDevicePointer());
        sortBoxDataArgs.push_back(&blockBoundingBox.getDevicePointer());
        sortBoxDataArgs.push_back(&sortedBlockCenter.getDevicePointer());
        sortBoxDataArgs.push_back(&sortedBlockBoundingBox.getDevicePointer());
        if (useLargeBlocks) {
            sortBoxDataArgs.push_back(&largeBlockCenter.getDevicePointer());
            sortBoxDataArgs.push_back(&largeBlockBoundingBox.getDevicePointer());
            sortBoxDataArgs.push_back(context.getPeriodicBoxSizePointer());
            sortBoxDataArgs.push_back(context.getInvPeriodicBoxSizePointer());
            sortBoxDataArgs.push_back(context.getPeriodicBoxVecXPointer());
            sortBoxDataArgs.push_back(context.getPeriodicBoxVecYPointer());
            sortBoxDataArgs.push_back(context.getPeriodicBoxVecZPointer());
        }
        sortBoxDataArgs.push_back(&context.getPosq().getDevicePointer());
        sortBoxDataArgs.push_back(&oldPositions.getDevicePointer());
        sortBoxDataArgs.push_back(&interactionCount.getDevicePointer());
        sortBoxDataArgs.push_back(&rebuildNeighborList.getDevicePointer());
        sortBoxDataArgs.push_back(&forceRebuildNeighborList);
        sortBoxDataArgs.push_back(&blockSizeRange.getDevicePointer());
        findInteractingBlocksArgs.push_back(context.getPeriodicBoxSizePointer());
        findInteractingBlocksArgs.push_back(context.getInvPeriodicBoxSizePointer());
        findInteractingBlocksArgs.push_back(context.getPeriodicBoxVecXPointer());
        findInteractingBlocksArgs.push_back(context.getPeriodicBoxVecYPointer());
        findInteractingBlocksArgs.push_back(context.getPeriodicBoxVecZPointer());
        findInteractingBlocksArgs.push_back(&interactionCount.getDevicePointer());
        findInteractingBlocksArgs.push_back(&interactingTiles.getDevicePointer());
        findInteractingBlocksArgs.push_back(&interactingAtoms.getDevicePointer());
        findInteractingBlocksArgs.push_back(&singlePairs.getDevicePointer());
        findInteractingBlocksArgs.push_back(&context.getPosq().getDevicePointer());
        findInteractingBlocksArgs.push_back(&maxTiles);
        findInteractingBlocksArgs.push_back(&maxSinglePairs);
        findInteractingBlocksArgs.push_back(&startBlockIndex);
        findInteractingBlocksArgs.push_back(&numBlocks);
        findInteractingBlocksArgs.push_back(&sortedBlocks.getDevicePointer());
        findInteractingBlocksArgs.push_back(&sortedBlockCenter.getDevicePointer());
        findInteractingBlocksArgs.push_back(&sortedBlockBoundingBox.getDevicePointer());
        if (useLargeBlocks) {
            findInteractingBlocksArgs.push_back(&largeBlockCenter.getDevicePointer());
            findInteractingBlocksArgs.push_back(&largeBlockBoundingBox.getDevicePointer());
        }
        findInteractingBlocksArgs.push_back(&exclusionIndices.getDevicePointer());
        findInteractingBlocksArgs.push_back(&exclusionRowIndices.getDevicePointer());
        findInteractingBlocksArgs.push_back(&oldPositions.getDevicePointer());
        findInteractingBlocksArgs.push_back(&rebuildNeighborList.getDevicePointer());
        copyInteractionCountsArgs.push_back(&interactionCount.getDevicePointer());
        copyInteractionCountsArgs.push_back(&pinnedCountBuffer);
    }
}

double HipNonbondedUtilities::getMaxCutoffDistance() {
    double cutoff = 0.0;
    for (map<int, double>::const_iterator iter = groupCutoff.begin(); iter != groupCutoff.end(); ++iter)
        cutoff = max(cutoff, iter->second);
    return cutoff;
}

double HipNonbondedUtilities::padCutoff(double cutoff) {
    double padding = (usePadding ? 0.08*cutoff : 0.0);
    return cutoff+padding;
}

void HipNonbondedUtilities::prepareInteractions(int forceGroups) {
    if ((forceGroups&groupFlags) == 0)
        return;
    if (groupKernels.find(forceGroups) == groupKernels.end())
        createKernelsForGroups(forceGroups);
    KernelSet& kernels = groupKernels[forceGroups];
    if (useCutoff && usePeriodic) {
        double4 box = context.getPeriodicBoxSize();
        double minAllowedSize = 1.999999*kernels.cutoffDistance;
        if (box.x < minAllowedSize || box.y < minAllowedSize || box.z < minAllowedSize)
            throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
    }
    if (!useNeighborList)
        return;
    if (numTiles == 0)
        return;

    // Compute the neighbor list.

    if (lastCutoff != kernels.cutoffDistance)
        forceRebuildNeighborList = true;
    context.executeKernelFlat(kernels.findBlockBoundsKernel, &findBlockBoundsArgs[0], context.getPaddedNumAtoms(), context.getSIMDWidth());
    context.executeKernelFlat(kernels.computeSortKeysKernel, &computeSortKeysArgs[0], context.getNumAtomBlocks());
    blockSorter->sort(sortedBlocks);
    context.executeKernelFlat(kernels.sortBoxDataKernel, &sortBoxDataArgs[0], context.getNumAtoms(), 64);
    context.executeKernelFlat(kernels.findInteractingBlocksKernel, &findInteractingBlocksArgs[0], context.getNumAtomBlocks() * context.getSIMDWidth() * numTilesInBatch, findInteractingBlocksThreadBlockSize);
    forceRebuildNeighborList = false;
    lastCutoff = kernels.cutoffDistance;
    context.executeKernelFlat(kernels.copyInteractionCountsKernel, &copyInteractionCountsArgs[0], 1, 1);
    hipEventRecord(downloadCountEvent, context.getCurrentStream());
}

void HipNonbondedUtilities::computeInteractions(int forceGroups, bool includeForces, bool includeEnergy) {
    if ((forceGroups&groupFlags) == 0)
        return;
    KernelSet& kernels = groupKernels[forceGroups];
    if (kernels.hasForces) {
        hipFunction_t& kernel = (includeForces ? (includeEnergy ? kernels.forceEnergyKernel : kernels.forceKernel) : kernels.energyKernel);
        if (kernel == NULL)
            kernel = createInteractionKernel(kernels.source, parameters, arguments, true, true, forceGroups, includeForces, includeEnergy);
        context.executeKernelFlat(kernel, &forceArgs[0], numForceThreadBlocks*forceThreadBlockSize, forceThreadBlockSize);
    }
    if (useNeighborList && numTiles > 0) {
        hipEventSynchronize(downloadCountEvent);
        updateNeighborListSize();
    }
}

bool HipNonbondedUtilities::updateNeighborListSize() {
    if (!useCutoff)
        return false;
    if (context.getStepsSinceReorder() == 0 || tilesAfterReorder == 0)
        tilesAfterReorder = pinnedCountBuffer[0];
    else if (context.getStepsSinceReorder() > 25 && pinnedCountBuffer[0] > 1.1*tilesAfterReorder)
        context.forceReorder();
    if (pinnedCountBuffer[0] <= maxTiles && pinnedCountBuffer[1] <= maxSinglePairs)
        return false;

    // The most recent timestep had too many interactions to fit in the arrays.  Make the arrays bigger to prevent
    // this from happening in the future.

    if (pinnedCountBuffer[0] > maxTiles) {
        maxTiles = (unsigned int) (1.2*pinnedCountBuffer[0]);
        unsigned int numBlocks = context.getNumAtomBlocks();
        int totalTiles = numBlocks*(numBlocks+1)/2;
        if (maxTiles > totalTiles)
            maxTiles = totalTiles;
        interactingTiles.resize(maxTiles);
        interactingAtoms.resize(HipContext::TileSize*(size_t) maxTiles);
        if (forceArgs.size() > 0)
            forceArgs[7] = &interactingTiles.getDevicePointer();
        findInteractingBlocksArgs[6] = &interactingTiles.getDevicePointer();
        if (forceArgs.size() > 0)
            forceArgs[17] = &interactingAtoms.getDevicePointer();
        findInteractingBlocksArgs[7] = &interactingAtoms.getDevicePointer();
    }
    if (pinnedCountBuffer[1] > maxSinglePairs) {
        maxSinglePairs = (unsigned int) (1.2*pinnedCountBuffer[1]);
        singlePairs.resize(maxSinglePairs);
        if (forceArgs.size() > 0)
            forceArgs[19] = &singlePairs.getDevicePointer();
        findInteractingBlocksArgs[8] = &singlePairs.getDevicePointer();
    }
    forceRebuildNeighborList = true;
    context.setForcesValid(false);
    return true;
}

void HipNonbondedUtilities::setUsePadding(bool padding) {
    usePadding = padding;
}

void HipNonbondedUtilities::setAtomBlockRange(double startFraction, double endFraction) {
    int numAtomBlocks = context.getNumAtomBlocks();
    startBlockIndex = (int) (startFraction*numAtomBlocks);
    numBlocks = (int) (endFraction*numAtomBlocks)-startBlockIndex;
    long long totalTiles = context.getNumAtomBlocks()*((long long)context.getNumAtomBlocks()+1)/2;
    startTileIndex = (int) (startFraction*totalTiles);
    numTiles = (long long) (endFraction*totalTiles)-startTileIndex;
    forceRebuildNeighborList = true;
}

void HipNonbondedUtilities::createKernelsForGroups(int groups) {
    KernelSet kernels;
    double cutoff = 0.0;
    string source;
    for (int i = 0; i < 32; i++) {
        if ((groups&(1<<i)) != 0) {
            cutoff = max(cutoff, groupCutoff[i]);
            source += groupKernelSource[i];
        }
    }
    kernels.hasForces = (source.size() > 0);
    kernels.cutoffDistance = cutoff;
    kernels.source = source;
    kernels.forceKernel = kernels.energyKernel = kernels.forceEnergyKernel = NULL;
    if (useCutoff) {
        double paddedCutoff = padCutoff(cutoff);
        map<string, string> defines;
        defines["TILE_SIZE"] = context.intToString(HipContext::TileSize);
        defines["NUM_BLOCKS"] = context.intToString(context.getNumAtomBlocks());
        defines["NUM_ATOMS"] = context.intToString(context.getNumAtoms());
        defines["PADDED_NUM_ATOMS"] = context.intToString(context.getPaddedNumAtoms());
        defines["PADDING"] = context.doubleToString(paddedCutoff-cutoff);
        defines["PADDED_CUTOFF"] = context.doubleToString(paddedCutoff);
        defines["PADDED_CUTOFF_SQUARED"] = context.doubleToString(paddedCutoff*paddedCutoff);
        defines["NUM_TILES_WITH_EXCLUSIONS"] = context.intToString(exclusionTiles.getSize());
        if (usePeriodic)
            defines["USE_PERIODIC"] = "1";
        if (context.getBoxIsTriclinic())
            defines["TRICLINIC"] = "1";
        if (useLargeBlocks)
            defines["USE_LARGE_BLOCKS"] = "1";
        defines["MAX_EXCLUSIONS"] = context.intToString(maxExclusions);
        int maxBits = 0;
        if (canUsePairList) {
            if (context.getUseDoublePrecision()) {
                maxBits = 4;
            }
            else {
                if (context.getSIMDWidth() > 32) {
                    // CDNA
                    if (context.getNumAtoms() < 100000)
                        maxBits = 4;
                    else // Large systems
                        maxBits = 0;
                }
                else {
                    // RDNA
                    if (context.getNumAtoms() < 100000)
                        maxBits = 4;
                    else if (context.getNumAtoms() < 500000)
                        maxBits = 2;
                    else // Very large systems
                        maxBits = 0;
                }
            }
        }
        defines["MAX_BITS_FOR_PAIRS"] = context.intToString(maxBits);
        defines["NUM_TILES_IN_BATCH"] = context.intToString(numTilesInBatch);
        defines["GROUP_SIZE"] = context.intToString(findInteractingBlocksThreadBlockSize);
        int binShift = 1;
        while (1<<binShift <= context.getNumAtomBlocks())
            binShift++;
        defines["BIN_SHIFT"] = context.intToString(binShift);
        defines["BLOCK_INDEX_MASK"] = context.intToString((1<<binShift)-1);
        hipModule_t interactingBlocksProgram = context.createModule(HipKernelSources::vectorOps+HipKernelSources::findInteractingBlocks, defines);
        kernels.findBlockBoundsKernel = context.getKernel(interactingBlocksProgram, "findBlockBounds");
        kernels.computeSortKeysKernel = context.getKernel(interactingBlocksProgram, "computeSortKeys");
        kernels.sortBoxDataKernel = context.getKernel(interactingBlocksProgram, "sortBoxData");
        kernels.findInteractingBlocksKernel = context.getKernel(interactingBlocksProgram, "findBlocksWithInteractions");
        kernels.copyInteractionCountsKernel = context.getKernel(interactingBlocksProgram, "copyInteractionCounts");
    }
    groupKernels[groups] = kernels;
}

hipFunction_t HipNonbondedUtilities::createInteractionKernel(const string& source, vector<ParameterInfo>& params, vector<ParameterInfo>& arguments, bool useExclusions, bool isSymmetric, int groups, bool includeForces, bool includeEnergy) {
    map<string, string> replacements;
    replacements["COMPUTE_INTERACTION"] = source;
    const string suffixes[] = {"x", "y", "z", "w"};
    stringstream args;
    for (const ParameterInfo& param : params) {
        args << ", ";
        if (param.isConstant())
            args << "const ";
        args << param.getType();
        args << "* __restrict__ global_";
        args << param.getName();
    }
    for (const ParameterInfo& arg : arguments) {
        args << ", ";
        if (arg.isConstant())
            args << "const ";
        args << arg.getType();
        args << "* __restrict__ ";
        args << arg.getName();
    }
    if (energyParameterDerivatives.size() > 0)
        args << ", mixed* __restrict__ energyParamDerivs";
    replacements["PARAMETER_ARGUMENTS"] = args.str();

    stringstream load1;
    for (const ParameterInfo& param : params) {
        load1 << param.getType();
        load1 << " ";
        load1 << param.getName();
        load1 << "1 = global_";
        load1 << param.getName();
        load1 << "[atom1];\n";
    }
    replacements["LOAD_ATOM1_PARAMETERS"] = load1.str();

    // Part 1. Defines for on diagonal exclusion tiles

    stringstream broadcastWarpData;
    broadcastWarpData << "posq2.x = SHFL(shflPosq.x, j);\n";
    broadcastWarpData << "posq2.y = SHFL(shflPosq.y, j);\n";
    broadcastWarpData << "posq2.z = SHFL(shflPosq.z, j);\n";
    broadcastWarpData << "posq2.w = SHFL(shflPosq.w, j);\n";
    for (const ParameterInfo& param : params) {
        broadcastWarpData << param.getType() << " shfl" << param.getName() << ";\n";
        for (int j = 0; j < param.getNumComponents(); j++) {
            if (param.getNumComponents() == 1)
                broadcastWarpData << "shfl" << param.getName() << "=SHFL(" << param.getName() <<"1,j);\n";
            else
                broadcastWarpData << "shfl" << param.getName()+"."+suffixes[j] << "=SHFL(" << param.getName()+"1."+suffixes[j] <<",j);\n";
        }
    }
    replacements["BROADCAST_WARP_DATA"] = broadcastWarpData.str();

    // Part 2. Defines for off-diagonal exclusions, and neighborlist tiles.
    stringstream declareLocal2;
    for (const ParameterInfo& param : params)
        declareLocal2<<param.getType()<<" shfl"<<param.getName()<<";\n";
    replacements["DECLARE_LOCAL_PARAMETERS"] = declareLocal2.str();

    stringstream loadLocal2;
    for (const ParameterInfo& param : params)
        loadLocal2<<"shfl"<<param.getName()<<" = global_"<<param.getName()<<"[j];\n";
    replacements["LOAD_LOCAL_PARAMETERS_FROM_GLOBAL"] = loadLocal2.str();

    stringstream load2j;
    for (const ParameterInfo& param : params)
        load2j<<param.getType()<<" "<<param.getName()<<"2 = shfl"<<param.getName()<<";\n";
    replacements["LOAD_ATOM2_PARAMETERS"] = load2j.str();

    stringstream clearLocal;
    for (const ParameterInfo& param : params) {
        clearLocal<<"shfl";
        clearLocal<<param.getName()<<" = ";
        if (param.getNumComponents() == 1)
            clearLocal<<"0;\n";
        else
            clearLocal<<"make_"<<param.getType()<<"(0);\n";
    }
    replacements["CLEAR_LOCAL_PARAMETERS"] = clearLocal.str();

    stringstream initDerivs;
    for (int i = 0; i < energyParameterDerivatives.size(); i++)
        initDerivs<<"mixed energyParamDeriv"<<i<<" = 0;\n";
    replacements["INIT_DERIVATIVES"] = initDerivs.str();
    stringstream saveDerivs;
    const vector<string>& allParamDerivNames = context.getEnergyParamDerivNames();
    int numDerivs = allParamDerivNames.size();
    for (int i = 0; i < energyParameterDerivatives.size(); i++)
        for (int index = 0; index < numDerivs; index++)
            if (allParamDerivNames[index] == energyParameterDerivatives[i])
                saveDerivs<<"energyParamDerivs[GLOBAL_ID*"<<numDerivs<<"+"<<index<<"] += energyParamDeriv"<<i<<";\n";
    replacements["SAVE_DERIVATIVES"] = saveDerivs.str();

    stringstream shuffleWarpData;
    shuffleWarpData << "shflPosq = warpRotateLeft<TILE_SIZE>(shflPosq);\n";
    shuffleWarpData << "shflForce = warpRotateLeft<TILE_SIZE>(shflForce);\n";
    for (const ParameterInfo& param : params) {
        shuffleWarpData<<"shfl"<<param.getName()<<"=warpRotateLeft<TILE_SIZE>(shfl"<<param.getName()<<");\n";
    }
    replacements["SHUFFLE_WARP_DATA"] = shuffleWarpData.str();

    map<string, string> defines;
    if (useCutoff)
        defines["USE_CUTOFF"] = "1";
    if (usePeriodic)
        defines["USE_PERIODIC"] = "1";
    if (useExclusions)
        defines["USE_EXCLUSIONS"] = "1";
    if (isSymmetric)
        defines["USE_SYMMETRIC"] = "1";
    if (useNeighborList)
        defines["USE_NEIGHBOR_LIST"] = "1";
    defines["ENABLE_SHUFFLE"] = "1"; // Used only in hippoNonbonded.cc
    if (includeForces)
        defines["INCLUDE_FORCES"] = "1";
    if (includeEnergy)
        defines["INCLUDE_ENERGY"] = "1";
    defines["THREAD_BLOCK_SIZE"] = context.intToString(forceThreadBlockSize);
    double maxCutoff = 0.0;
    for (int i = 0; i < 32; i++) {
        if ((groups&(1<<i)) != 0) {
            double cutoff = groupCutoff[i];
            maxCutoff = max(maxCutoff, cutoff);
            defines["CUTOFF_"+context.intToString(i)+"_SQUARED"] = context.doubleToString(cutoff*cutoff);
            defines["CUTOFF_"+context.intToString(i)] = context.doubleToString(cutoff);
        }
    }
    defines["MAX_CUTOFF"] = context.doubleToString(maxCutoff);
    defines["MAX_CUTOFF_SQUARED"] = context.doubleToString(maxCutoff*maxCutoff);
    defines["NUM_ATOMS"] = context.intToString(context.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = context.intToString(context.getPaddedNumAtoms());
    defines["NUM_BLOCKS"] = context.intToString(context.getNumAtomBlocks());
    defines["TILE_SIZE"] = context.intToString(HipContext::TileSize);
    int numExclusionTiles = exclusionTiles.getSize();
    defines["NUM_TILES_WITH_EXCLUSIONS"] = context.intToString(numExclusionTiles);
    int numContexts = context.getPlatformData().contexts.size();
    int startExclusionIndex = context.getContextIndex()*numExclusionTiles/numContexts;
    int endExclusionIndex = (context.getContextIndex()+1)*numExclusionTiles/numContexts;
    defines["FIRST_EXCLUSION_TILE"] = context.intToString(startExclusionIndex);
    defines["LAST_EXCLUSION_TILE"] = context.intToString(endExclusionIndex);
    hipModule_t program = context.createModule(HipKernelSources::vectorOps+context.replaceStrings(kernelSource, replacements), defines);
    hipFunction_t kernel = context.getKernel(program, "computeNonbonded");
    return kernel;
}

void HipNonbondedUtilities::setKernelSource(const string& source) {
    kernelSource = source;
}
