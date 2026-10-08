/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2009-2026 Stanford University and the Authors.      *
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

#include "MetalNonbondedUtilities.h"
#include "MetalArray.h"
#include "MetalContext.h"
#include "MetalKernelSources.h"
#include "openmm/OpenMMException.h"
#include "openmm/common/CommonKernelUtilities.h"
#include <algorithm>
#include <map>
#include <set>
#include <utility>

using namespace OpenMM;
using namespace std;

class MetalNonbondedUtilities::BlockSortTrait : public ComputeSortImpl::SortTrait {
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

MetalNonbondedUtilities::MetalNonbondedUtilities(MetalContext& context) : context(context), useCutoff(false), usePeriodic(false), useNeighborList(false), anyExclusions(false), usePadding(true),
        forceRebuildNeighborList(true), groupFlags(0), canUsePairList(true), tilesAfterReorder(0) {
    // Decide how many thread blocks to use.

    string errorMessage = "Error initializing nonbonded utilities";
    numForceThreadBlocks = 6*context.getNumGPUCores();
    forceThreadBlockSize = 256;

    // When building the neighbor list, we can optionally use large blocks (1024 atoms) to
    // accelerate the process.  This makes building the neighbor list faster, but it prevents
    // us from sorting atom blocks by size, which leads to a slightly less efficient neighbor
    // list.  We guess based on system size which will be faster.

    useLargeBlocks = (context.getNumAtoms() > 90000);
    setKernelSource(MetalKernelSources::nonbonded);
}

MetalNonbondedUtilities::~MetalNonbondedUtilities() {
}

void MetalNonbondedUtilities::addInteraction(bool usesCutoff, bool usesPeriodic, bool usesExclusions, double cutoffDistance,
        const vector<vector<int> >& exclusionList, const string& kernel, int forceGroup, bool useNeighborList, bool supportsPairList) {
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
    this->useNeighborList |= (useNeighborList && useCutoff);
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

void MetalNonbondedUtilities::addParameter(ComputeParameterInfo parameter) {
    parameters.push_back(parameter);
}

void MetalNonbondedUtilities::addArgument(ComputeParameterInfo parameter) {
    arguments.push_back(parameter);
}

string MetalNonbondedUtilities::addEnergyParameterDerivative(const string& param) {
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

void MetalNonbondedUtilities::requestExclusions(const vector<vector<int> >& exclusionList) {
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

static bool compareInt2(mm_int2 a, mm_int2 b) {
    return ((a.y < b.y) || (a.y == b.y && a.x < b.x));
}

void MetalNonbondedUtilities::initialize(const System& system) {
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
        int x = atom1/MetalContext::TileSize;
        for (int j = 0; j < (int) atomExclusions[atom1].size(); ++j) {
            int atom2 = atomExclusions[atom1][j];
            int y = atom2/MetalContext::TileSize;
            tilesWithExclusions.insert(make_pair(max(x, y), min(x, y)));
        }
    }
    vector<mm_int2> exclusionTilesVec;
    for (set<pair<int, int> >::const_iterator iter = tilesWithExclusions.begin(); iter != tilesWithExclusions.end(); ++iter)
        exclusionTilesVec.push_back(mm_int2(iter->first, iter->second));
    sort(exclusionTilesVec.begin(), exclusionTilesVec.end(), compareInt2);
    exclusionTiles.initialize<mm_int2>(context, exclusionTilesVec.size(), "exclusionTiles");
    exclusionTiles.upload(exclusionTilesVec);
    map<pair<int, int>, int> exclusionTileMap;
    for (int i = 0; i < (int) exclusionTilesVec.size(); i++) {
        mm_int2 tile = exclusionTilesVec[i];
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

    exclusions.initialize<tileflags>(context, tilesWithExclusions.size()*MetalContext::TileSize, "exclusions");
    tileflags allFlags = (tileflags) -1;
    vector<tileflags> exclusionVec(exclusions.getSize(), allFlags);
    for (int atom1 = 0; atom1 < (int) atomExclusions.size(); ++atom1) {
        int x = atom1/MetalContext::TileSize;
        int offset1 = atom1-x*MetalContext::TileSize;
        for (int j = 0; j < (int) atomExclusions[atom1].size(); ++j) {
            int atom2 = atomExclusions[atom1][j];
            int y = atom2/MetalContext::TileSize;
            int offset2 = atom2-y*MetalContext::TileSize;
            if (x > y) {
                int index = exclusionTileMap[make_pair(x, y)]*MetalContext::TileSize;
                exclusionVec[index+offset1] &= allFlags-(1<<offset2);
            }
            else {
                int index = exclusionTileMap[make_pair(y, x)]*MetalContext::TileSize;
                exclusionVec[index+offset2] &= allFlags-(1<<offset1);
            }
        }
    }
    atomExclusions.clear(); // We won't use this again, so free the memory it used
    exclusions.upload(exclusionVec);

    // Create data structures for the neighbor list.

    maxCutoff = getMaxCutoffDistance();
    if (useCutoff) {
        // Select a size for the arrays that hold the neighbor list.  We have to make a fairly
        // arbitrary guess, but if this turns out to be too small we'll increase it later.

        maxTiles = 20*numAtomBlocks;
        if (maxTiles > numTiles)
            maxTiles = numTiles;
        if (maxTiles < 1)
            maxTiles = 1;
        maxSinglePairs = 5*numAtoms;
        interactingTiles.initialize<int>(context, maxTiles, "interactingTiles");
        interactingAtoms.initialize<int>(context, MetalContext::TileSize*maxTiles, "interactingAtoms");
        interactionCount.initialize<unsigned int>(context, 2, "interactionCount");
        singlePairs.initialize<mm_int2>(context, maxSinglePairs, "singlePairs");
        blockCenter.initialize<mm_float4>(context, numAtomBlocks, "blockCenter");
        blockBoundingBox.initialize<mm_float4>(context, numAtomBlocks, "blockBoundingBox");
        sortedBlocks.initialize<unsigned int>(context, numAtomBlocks, "sortedBlocks");
        sortedBlockCenter.initialize<mm_float4>(context, numAtomBlocks+1, "sortedBlockCenter");
        sortedBlockBoundingBox.initialize<mm_float4>(context, numAtomBlocks+1, "sortedBlockBoundingBox");
        numBlockSizes = min((context.getNumAtomBlocks()+63)/64, context.getNumThreadBlocks());
        blockSizeRange.initialize<mm_float2>(context, numBlockSizes, "blockSizeRange");
        largeBlockCenter.initialize<mm_float4>(context, numAtomBlocks, "largeBlockCenter");
        largeBlockBoundingBox.initialize<mm_float4>(context, numAtomBlocks, "largeBlockBoundingBox");
        oldPositions.initialize<mm_float4>(context, numAtoms, "oldPositions");
        rebuildNeighborList.initialize<int>(context, 1, "rebuildNeighborList");
        blockSorter = context.createSort(new BlockSortTrait(), numAtomBlocks, false);
        vector<unsigned int> count(2, 0);
        interactionCount.upload(count);
        rebuildNeighborList.upload(&count[0]);
        event = context.createEvent();
    }
}

double MetalNonbondedUtilities::getMaxCutoffDistance() {
    double cutoff = 0.0;
    for (map<int, double>::const_iterator iter = groupCutoff.begin(); iter != groupCutoff.end(); ++iter)
        cutoff = max(cutoff, iter->second);
    return cutoff;
}

double MetalNonbondedUtilities::padCutoff(double cutoff) {
    double padding = (usePadding ? 0.08*cutoff : 0.0);
    return cutoff+padding;
}

void MetalNonbondedUtilities::prepareInteractions(int forceGroups) {
    if ((forceGroups&groupFlags) == 0)
        return;
    if (groupKernels.find(forceGroups) == groupKernels.end())
        createKernelsForGroups(forceGroups);
    KernelSet& kernels = groupKernels[forceGroups];
    if (useCutoff && usePeriodic) {
        mm_double4 box = context.getPeriodicBoxSize();
        double minAllowedSize = 1.999999*maxCutoff;
        if (box.x < minAllowedSize || box.y < minAllowedSize || box.z < minAllowedSize)
            throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
    }
    if (!useNeighborList)
        return;
    if (numTiles == 0)
        return;

    // Compute the neighbor list.

    setPeriodicBoxArgs(context, kernels.findBlockBoundsKernel, 1);
    kernels.findBlockBoundsKernel->execute(context.getNumAtomBlocks());
    kernels.computeSortKeysKernel->execute(context.getNumAtomBlocks());
    if (useLargeBlocks)
        setPeriodicBoxArgs(context, kernels.sortBoxDataKernel, 7);
    blockSorter->sort(sortedBlocks);
    kernels.sortBoxDataKernel->setArg(useLargeBlocks ? 16 : 9, forceRebuildNeighborList);
    kernels.sortBoxDataKernel->execute(context.getNumAtoms());
    setPeriodicBoxArgs(context, kernels.findInteractingBlocksKernel, 0);
    kernels.findInteractingBlocksKernel->execute(context.getNumAtoms(), 256);
    forceRebuildNeighborList = false;
    if (useNeighborList && numTiles > 0)
        event->enqueue();
}

void MetalNonbondedUtilities::computeInteractions(int forceGroups, bool includeForces, bool includeEnergy) {
    if ((forceGroups&groupFlags) == 0)
        return;
    KernelSet& kernels = groupKernels[forceGroups];
    if (kernels.hasForces && (includeForces || includeEnergy)) {
        ComputeKernel& kernel = (includeForces ? (includeEnergy ? kernels.forceEnergyKernel : kernels.forceKernel) : kernels.energyKernel);
        if (kernel.use_count() == 0)
            kernel = createInteractionKernel(kernels.source, parameters, arguments, true, true, forceGroups, includeForces, includeEnergy);
        if (useCutoff)
            setPeriodicBoxArgs(context, kernel, 9);
        kernel->execute(numForceThreadBlocks*forceThreadBlockSize, forceThreadBlockSize);
    }
    if (useNeighborList && numTiles > 0) {
        event->wait();
        updateNeighborListSize();
    }
}

bool MetalNonbondedUtilities::updateNeighborListSize() {
    if (!useCutoff)
        return false;
    unsigned int* countBuffer = (unsigned int*) interactionCount.getBuffer()->contents();
    if (context.getStepsSinceReorder() == 0 || tilesAfterReorder == 0)
        tilesAfterReorder = countBuffer[0];
    else if (context.getStepsSinceReorder() > 25 && countBuffer[0] > 1.1*tilesAfterReorder)
        context.forceReorder();
    if (countBuffer[0] <= maxTiles && countBuffer[1] <= maxSinglePairs)
        return false;

    // The most recent timestep had too many interactions to fit in the arrays.  Make the arrays bigger to prevent
    // this from happening in the future.

    if (countBuffer[0] > maxTiles) {
        maxTiles = (unsigned int) (1.2*countBuffer[0]);
        unsigned int numBlocks = context.getNumAtomBlocks();
        int totalTiles = numBlocks*(numBlocks+1)/2;
        if (maxTiles > totalTiles)
            maxTiles = totalTiles;
        interactingTiles.resize(maxTiles);
        interactingAtoms.resize(MetalContext::TileSize*(size_t) maxTiles);
    }
    if (countBuffer[1] > maxSinglePairs) {
        maxSinglePairs = (unsigned int) (1.2*countBuffer[1]);
        singlePairs.resize(maxSinglePairs);
    }
    for (auto& entry : groupKernels) {
        KernelSet& kernels = entry.second;
        kernels.findInteractingBlocksKernel->setArg(10, maxTiles);
        kernels.findInteractingBlocksKernel->setArg(11, maxSinglePairs);
        for (const ComputeKernel& kernel : {kernels.forceKernel, kernels.energyKernel, kernels.forceEnergyKernel})
            if (kernel != nullptr) {
                kernel->setArg(14, maxTiles);
                kernel->setArg(18, maxSinglePairs);
            }
    }
    forceRebuildNeighborList = true;
    context.setForcesValid(false);
    return true;
}

void MetalNonbondedUtilities::setUsePadding(bool padding) {
    usePadding = padding;
}

void MetalNonbondedUtilities::setAtomBlockRange(double startFraction, double endFraction) {
    int numAtomBlocks = context.getNumAtomBlocks();
    startBlockIndex = (int) (startFraction*numAtomBlocks);
    numBlocks = (int) (endFraction*numAtomBlocks)-startBlockIndex;
    long long totalTiles = context.getNumAtomBlocks()*((long long)context.getNumAtomBlocks()+1)/2;
    startTileIndex = (int) (startFraction*totalTiles);
    numTiles = (long long) (endFraction*totalTiles)-startTileIndex;
    forceRebuildNeighborList = true;
}

void MetalNonbondedUtilities::createKernelsForGroups(int groups) {
    KernelSet kernels;
    string source;
    for (int i = 0; i < 32; i++) {
        if ((groups&(1<<i)) != 0) {
            source += groupKernelSource[i];
        }
    }
    kernels.hasForces = (source.size() > 0);
    kernels.source = source;
    kernels.forceKernel = kernels.energyKernel = kernels.forceEnergyKernel = NULL;
    if (useCutoff) {
        double paddedCutoff = padCutoff(maxCutoff);
        map<string, string> defines;
        defines["TILE_SIZE"] = context.intToString(MetalContext::TileSize);
        defines["NUM_BLOCKS"] = context.intToString(context.getNumAtomBlocks());
        defines["NUM_ATOMS"] = context.intToString(context.getNumAtoms());
        defines["PADDING"] = context.doubleToString(paddedCutoff-maxCutoff);
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
        defines["MAX_BITS_FOR_PAIRS"] = (canUsePairList ? "3" : "0");
        int binShift = 1;
        while (1<<binShift <= context.getNumAtomBlocks())
            binShift++;
        defines["BIN_SHIFT"] = context.intToString(binShift);
        defines["BLOCK_INDEX_MASK"] = context.intToString((1<<binShift)-1);

        // Create the kernels.

        ComputeProgram interactingBlocksProgram = context.compileProgram(MetalKernelSources::findInteractingBlocks, defines);
        kernels.findBlockBoundsKernel = interactingBlocksProgram->createKernel("findBlockBounds");
        kernels.computeSortKeysKernel = interactingBlocksProgram->createKernel("computeSortKeys");
        kernels.sortBoxDataKernel = interactingBlocksProgram->createKernel("sortBoxData");
        kernels.findInteractingBlocksKernel = interactingBlocksProgram->createKernel("findBlocksWithInteractions");

        // Set arguments for kernels.

        kernels.findBlockBoundsKernel->addArg(numAtoms);
        for (int i = 0; i < 5; i++)
            kernels.findBlockBoundsKernel->addArg();
        kernels.findBlockBoundsKernel->addArg(context.getPosq());
        kernels.findBlockBoundsKernel->addArg(blockCenter);
        kernels.findBlockBoundsKernel->addArg(blockBoundingBox);
        kernels.findBlockBoundsKernel->addArg(rebuildNeighborList);
        kernels.findBlockBoundsKernel->addArg(blockSizeRange);
        kernels.computeSortKeysKernel->addArg(blockBoundingBox);
        kernels.computeSortKeysKernel->addArg(sortedBlocks);
        kernels.computeSortKeysKernel->addArg(blockSizeRange);
        kernels.computeSortKeysKernel->addArg(numBlockSizes);
        kernels.sortBoxDataKernel->addArg(sortedBlocks);
        kernels.sortBoxDataKernel->addArg(blockCenter);
        kernels.sortBoxDataKernel->addArg(blockBoundingBox);
        kernels.sortBoxDataKernel->addArg(sortedBlockCenter);
        kernels.sortBoxDataKernel->addArg(sortedBlockBoundingBox);
        if (useLargeBlocks) {
            kernels.sortBoxDataKernel->addArg(largeBlockCenter);
            kernels.sortBoxDataKernel->addArg(largeBlockBoundingBox);
            for (int i = 0; i < 5; i++)
                kernels.sortBoxDataKernel->addArg();
        }
        kernels.sortBoxDataKernel->addArg(context.getPosq());
        kernels.sortBoxDataKernel->addArg(oldPositions);
        kernels.sortBoxDataKernel->addArg(interactionCount);
        kernels.sortBoxDataKernel->addArg(rebuildNeighborList);
        kernels.sortBoxDataKernel->addArg();
        for (int i = 0; i < 5; i++)
            kernels.findInteractingBlocksKernel->addArg();
        kernels.findInteractingBlocksKernel->addArg(interactionCount);
        kernels.findInteractingBlocksKernel->addArg(interactingTiles);
        kernels.findInteractingBlocksKernel->addArg(interactingAtoms);
        kernels.findInteractingBlocksKernel->addArg(singlePairs);
        kernels.findInteractingBlocksKernel->addArg(context.getPosq());
        kernels.findInteractingBlocksKernel->addArg(maxTiles);
        kernels.findInteractingBlocksKernel->addArg(maxSinglePairs);
        kernels.findInteractingBlocksKernel->addArg(startBlockIndex);
        kernels.findInteractingBlocksKernel->addArg(numBlocks);
        kernels.findInteractingBlocksKernel->addArg(sortedBlocks);
        kernels.findInteractingBlocksKernel->addArg(sortedBlockCenter);
        kernels.findInteractingBlocksKernel->addArg(sortedBlockBoundingBox);
        if (useLargeBlocks) {
            kernels.findInteractingBlocksKernel->addArg(largeBlockCenter);
            kernels.findInteractingBlocksKernel->addArg(largeBlockBoundingBox);
        }
        kernels.findInteractingBlocksKernel->addArg(exclusionIndices);
        kernels.findInteractingBlocksKernel->addArg(exclusionRowIndices);
        kernels.findInteractingBlocksKernel->addArg(oldPositions);
        kernels.findInteractingBlocksKernel->addArg(rebuildNeighborList);
    }
    groupKernels[groups] = kernels;
}

ComputeKernel MetalNonbondedUtilities::createInteractionKernel(const string& source, vector<ComputeParameterInfo>& params, vector<ComputeParameterInfo>& arguments, bool useExclusions, bool isSymmetric, int groups, bool includeForces, bool includeEnergy) {
    map<string, string> replacements;
    replacements["COMPUTE_INTERACTION"] = source;
    const string suffixes[] = {"x", "y", "z", "w"};
    stringstream localData;
    int localDataSize = 0;
    for (const ComputeParameterInfo& param : params) {
        if (param.getNumComponents() == 1)
            localData<<param.getType()<<" "<<param.getName()<<";\n";
        else {
            for (int j = 0; j < param.getNumComponents(); ++j)
                localData<<param.getComponentType()<<" "<<param.getName()<<"_"<<suffixes[j]<<";\n";
        }
        localDataSize += param.getSize();
    }
    replacements["ATOM_PARAMETER_DATA"] = localData.str();
    stringstream args;
    for (const ComputeParameterInfo& param : params) {
        args << ", ";
        if (param.isConstant())
            args << "GLOBAL const ";
        args << param.getType();
        args << "* RESTRICT global_";
        args << param.getName();
    }
    for (const ComputeParameterInfo& arg : arguments) {
        args << ", ";
        if (arg.isConstant())
            args << "GLOBAL const ";
        args << arg.getType();
        args << "* RESTRICT ";
        args << arg.getName();
    }
    if (energyParameterDerivatives.size() > 0)
        args << ", GLOBAL mixed* RESTRICT energyParamDerivs";
    replacements["PARAMETER_ARGUMENTS"] = args.str();

    stringstream load1;
    for (const ComputeParameterInfo& param : params) {
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
    replacements["BROADCAST_WARP_DATA"] = broadcastWarpData.str();

    // Part 2. Defines for off-diagonal exclusions, and neighborlist tiles.

    stringstream load2g;
    for (const ComputeParameterInfo& param : params)
        load2g<<param.getType()<<" "<<param.getName()<<"2 = global_"<<param.getName()<<"[atom2];\n";
    replacements["LOAD_ATOM2_PARAMETERS_FROM_GLOBAL"] = load2g.str();

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
    shuffleWarpData << "shflPosq.x = SHFL(shflPosq.x, tgx+1);\n";
    shuffleWarpData << "shflPosq.y = SHFL(shflPosq.y, tgx+1);\n";
    shuffleWarpData << "shflPosq.z = SHFL(shflPosq.z, tgx+1);\n";
    shuffleWarpData << "shflPosq.w = SHFL(shflPosq.w, tgx+1);\n";
    shuffleWarpData << "shflForce.x = SHFL(shflForce.x, tgx+1);\n";
    shuffleWarpData << "shflForce.y = SHFL(shflForce.y, tgx+1);\n";
    shuffleWarpData << "shflForce.z = SHFL(shflForce.z, tgx+1);\n";
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
    defines["ENABLE_SHUFFLE"] = "1";
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
    defines["NUM_ATOMS"] = context.intToString(context.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = context.intToString(context.getPaddedNumAtoms());
    defines["NUM_BLOCKS"] = context.intToString(context.getNumAtomBlocks());
    defines["TILE_SIZE"] = context.intToString(MetalContext::TileSize);
    int numExclusionTiles = exclusionTiles.getSize();
    defines["NUM_TILES_WITH_EXCLUSIONS"] = context.intToString(numExclusionTiles);
    int numContexts = context.getPlatformData().contexts.size();
    int startExclusionIndex = context.getContextIndex()*numExclusionTiles/numContexts;
    int endExclusionIndex = (context.getContextIndex()+1)*numExclusionTiles/numContexts;
    defines["FIRST_EXCLUSION_TILE"] = context.intToString(startExclusionIndex);
    defines["LAST_EXCLUSION_TILE"] = context.intToString(endExclusionIndex);
    if ((localDataSize/4)%2 == 0 && !context.getUseDoublePrecision())
        defines["PARAMETER_SIZE_IS_EVEN"] = "1";
    ComputeProgram program = context.compileProgram(context.replaceStrings(kernelSource, replacements), defines);
    ComputeKernel kernel = program->createKernel("computeNonbonded");

    // Set arguments to the Kernel.

    kernel->addArg(context.getLongForceBuffer());
    kernel->addArg(context.getEnergyBuffer());
    kernel->addArg(context.getPosq());
    kernel->addArg(exclusions);
    kernel->addArg(exclusionTiles);
    kernel->addArg(startTileIndex);
    kernel->addArg(numTiles);
    if (useCutoff) {
        kernel->addArg(interactingTiles);
        kernel->addArg(interactionCount);
        for (int i = 0; i < 5; i++)
            kernel->addArg();
        kernel->addArg(maxTiles);
        kernel->addArg(blockCenter);
        kernel->addArg(blockBoundingBox);
        kernel->addArg(interactingAtoms);
        kernel->addArg(maxSinglePairs);
        kernel->addArg(singlePairs);
    }
    for (ComputeParameterInfo& param : parameters)
        kernel->addArg(param.getArray());
    for (ComputeParameterInfo& arg : arguments)
        kernel->addArg(arg.getArray());
    if (energyParameterDerivatives.size() > 0)
        kernel->addArg(context.getEnergyParamDerivBuffer());
    return kernel;
}

void MetalNonbondedUtilities::setKernelSource(const string& source) {
    kernelSource = source;
}
