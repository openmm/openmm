/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM OpenCL Platform.
 * Source: platforms/opencl/src/OpenCLNonbondedUtilities.cpp
 *
 * Original OpenCL Platform code:
 * Portions copyright (c) 2009-2025 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 *
 * Metal Platform code:
 * Portions copyright (c) 2026 Chun-Chi Hung.
 * Authors: Chun-Chi Hung
 *
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
#include "MetalNonbondedUtilities.h"
#include "openmm/common/ComputeArray.h"
#include "MetalContext.h"
#include "MetalQueue.h"
#include "MetalOpenCLKernelSources.h"
#include "MetalNonbondedSources.h"
#include <cstdint>
#include <limits>
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

MetalNonbondedUtilities::MetalNonbondedUtilities(MetalContext& context) : context(context), downloadedCount(0), countReadbackPending(false),
        useCutoff(false), usePeriodic(false), anyExclusions(false), usePadding(true), useNeighborList(false),
        forceRebuildNeighborList(true), canUsePairList(OPENMM_METAL_FAST_SPARSE_PAIRS), groupFlags(0), tilesAfterReorder(0) {
    numForceThreadBlocks = context.getNumThreadBlocks();
    forceThreadBlockSize = min(256, context.getMaxThreadBlockSize());
    forceThreadBlockSize -= forceThreadBlockSize%MetalContext::TileSize;
    if (forceThreadBlockSize < MetalContext::TileSize)
        throw OpenMMException("The Metal device cannot execute a nonbonded tile");
    useLargeBlocks = (context.getNumAtoms() > 100000);
    setKernelSource(MetalOpenCLKernelSources::nonbonded);
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
    canUsePairList &= supportsPairList;
    groupCutoff[forceGroup] = cutoffDistance;
    groupFlags |= 1u<<forceGroup;
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
    // This version is used on devices with SIMD width of 32 or less.  It sorts tiles to improve cache efficiency.

    return ((a.y < b.y) || (a.y == b.y && a.x < b.x));
}

static bool compareInt2LargeSIMD(mm_int2 a, mm_int2 b) {
    // This version is used on devices with SIMD width greater than 32.  It puts diagonal tiles before off-diagonal
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

void MetalNonbondedUtilities::initialize(const System& system) {
    if (atomExclusions.size() == 0) {
        // No exclusions were specifically requested, so just mark every atom as not interacting with itself.

        atomExclusions.resize(context.getNumAtoms());
        for (int i = 0; i < (int) atomExclusions.size(); i++)
            atomExclusions[i].push_back(i);
    }

    // Create the list of tiles.

    int numAtomBlocks = context.getNumAtomBlocks();
    int numContexts = context.getNumContexts();
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
    sort(exclusionTilesVec.begin(), exclusionTilesVec.end(), context.getSIMDWidth() <= 32 || !useNeighborList ? compareInt2 : compareInt2LargeSIMD);
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

    exclusions.initialize<unsigned int>(context, tilesWithExclusions.size()*MetalContext::TileSize, "exclusions");
    unsigned int allFlags = (unsigned int) -1;
    vector<unsigned int> exclusionVec(exclusions.getSize(), allFlags);
    for (int i = 0; i < exclusions.getSize(); ++i)
        exclusionVec[i] = 0xFFFFFFFF;
    for (int atom1 = 0; atom1 < (int) atomExclusions.size(); ++atom1) {
        int x = atom1/MetalContext::TileSize;
        int offset1 = atom1-x*MetalContext::TileSize;
        for (int j = 0; j < (int) atomExclusions[atom1].size(); ++j) {
            int atom2 = atomExclusions[atom1][j];
            int y = atom2/MetalContext::TileSize;
            int offset2 = atom2-y*MetalContext::TileSize;
            if (x > y) {
                int index = exclusionTileMap[make_pair(x, y)]*MetalContext::TileSize;
                exclusionVec[index+offset1] &= allFlags-(1u<<offset2);
            }
            else {
                int index = exclusionTileMap[make_pair(y, x)]*MetalContext::TileSize;
                exclusionVec[index+offset2] &= allFlags-(1u<<offset1);
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

        int maxTiles = 20*numAtomBlocks;
        if (maxTiles > numTiles)
            maxTiles = numTiles;
        if (maxTiles < 1)
            maxTiles = 1;
        int numAtoms = context.getNumAtoms();
        interactingTiles.initialize<int>(context, maxTiles, "interactingTiles");
        interactingAtoms.initialize<int>(context, MetalContext::TileSize*(size_t) maxTiles, "interactingAtoms");
        canUsePairList &= useNeighborList;
        interactionCount.initialize<unsigned int>(context, canUsePairList ? 2 : 1, "interactionCount");
        if (canUsePairList) {
            uint64_t maxPairs = max(uint64_t(1), uint64_t(5)*numAtoms);
            if (maxPairs > (uint64_t) numeric_limits<int>::max())
                throw OpenMMException("The sparse pair list exceeds the Metal array index limit");
            singlePairs.initialize<mm_int2>(context, (int) maxPairs, "singlePairs");
        }
        int elementSize = sizeof(float);
        int boundsElementSize = OPENMM_METAL_FAST_FP16_BOUNDS ? 2 : elementSize;
        blockCenter.initialize(context, numAtomBlocks, 4*elementSize, "blockCenter");
        blockBoundingBox.initialize(context, numAtomBlocks, 4*elementSize, "blockBoundingBox");
        sortedBlocks.initialize<unsigned int>(context, numAtomBlocks, "sortedBlocks");
        sortedBlockCenter.initialize(context, numAtomBlocks+1, 4*elementSize, "sortedBlockCenter");
        sortedBlockBoundingBox.initialize(context, numAtomBlocks+1, 4*boundsElementSize, "sortedBlockBoundingBox");
        numBlockSizes = min((context.getNumAtomBlocks()+63)/64, context.getNumThreadBlocks());
        blockSizeRange.initialize(context, numBlockSizes, 2*elementSize, "blockSizeRange");
        largeBlockCenter.initialize(context, numAtomBlocks, 4*elementSize, "largeBlockCenter");
        largeBlockBoundingBox.initialize(context, numAtomBlocks, 4*boundsElementSize, "largeBlockBoundingBox");
        oldPositions.initialize(context, numAtoms, 4*elementSize, "oldPositions");
        rebuildNeighborList.initialize<int>(context, 1, "rebuildNeighborList");
        blockSorter = context.createSort(new BlockSortTrait(), numAtomBlocks, false);
        vector<unsigned int> count(interactionCount.getSize(), 0);
        interactionCount.upload(count);
        count.resize(1);
        rebuildNeighborList.upload(count);
    }
}

/** @brief Bind the five single-precision periodic-box arguments used by OpenCL kernels. */
static void setPeriodicBoxArgs(MetalContext& context, ComputeKernel& kernel, int index) {
    Vec3 a, b, c;
    context.getPeriodicBoxVectors(a, b, c);
    kernel->setArg(index++, mm_float4(a[0], b[1], c[2], 0));
    kernel->setArg(index++, mm_float4(1/a[0], 1/b[1], 1/c[2], 0));
    kernel->setArg(index++, mm_float4(a[0], a[1], a[2], 0));
    kernel->setArg(index++, mm_float4(b[0], b[1], b[2], 0));
    kernel->setArg(index, mm_float4(c[0], c[1], c[2], 0));
}

/** @brief Reserve argument slots before binding the reused OpenCL kernel interface. */
static ComputeKernel createKernel(ComputeProgram& program, const string& name, int arguments) {
    ComputeKernel kernel = program->createKernel(name);
    for (int i = 0; i < arguments; i++)
        kernel->addArg();
    return kernel;
}

double MetalNonbondedUtilities::getMaxCutoffDistance() {
    double cutoff = 0.0;
    for (map<int, double>::const_iterator iter = groupCutoff.begin(); iter != groupCutoff.end(); ++iter)
        cutoff = max(cutoff, iter->second);
    return cutoff;
}

double MetalNonbondedUtilities::padCutoff(double cutoff) {
    double padding = (usePadding ? 0.1*cutoff : 0.0);
    return cutoff+padding;
}

void MetalNonbondedUtilities::prepareInteractions(int forceGroups) {
    if ((forceGroups&groupFlags) == 0)
        return;
    if (groupKernels.find(forceGroups) == groupKernels.end())
        createKernelsForGroups(forceGroups);
    KernelSet& kernels = groupKernels[forceGroups];
    if (useCutoff && usePeriodic) {
        Vec3 a, b, c;
        context.getPeriodicBoxVectors(a, b, c);
        mm_float4 box(a[0], b[1], c[2], 0);
        double minAllowedSize = 1.999999*maxCutoff;
        if (box.x < minAllowedSize || box.y < minAllowedSize || box.z < minAllowedSize)
            throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
    }
    if (!useNeighborList)
        return;
    if (numTiles == 0)
        return;

    // Recover a pending count if a force evaluation was interrupted before
    // computeInteractions(). Never overwrite its dedicated readback storage.
    if (countReadbackPending) {
        context.unwrap(interactionCount).finishDownload();
        countReadbackPending = false;
    }

    // Compute the neighbor list.

    setPeriodicBoxArgs(context, kernels.findBlockBoundsKernel, 1);
#if OPENMM_METAL_FAST_BLOCK_BOUNDS
    kernels.findBlockBoundsKernel->execute(numBlockSizes*64, 64);
#else
    kernels.findBlockBoundsKernel->execute(context.getNumAtomBlocks());
#endif
    kernels.computeSortKeysKernel->execute(context.getNumAtomBlocks());
    if (useLargeBlocks)
        setPeriodicBoxArgs(context, kernels.sortBoxDataKernel, 12);
    blockSorter->sort(sortedBlocks);
    kernels.sortBoxDataKernel->setArg(9, (int) (forceRebuildNeighborList));
    kernels.sortBoxDataKernel->execute(context.getNumAtoms());
    setPeriodicBoxArgs(context, kernels.findInteractingBlocksKernel, 0);
    kernels.findInteractingBlocksKernel->execute(context.getNumAtoms(), interactingBlocksThreadBlockSize);
    forceRebuildNeighborList = false;
    // Match OpenCL: read counts before standalone and fused force kernels, and
    // let the GPU start this phase while the host prepares the remaining work.
    context.unwrap(interactionCount).beginDownload();
    countReadbackPending = true;
}

void MetalNonbondedUtilities::computeInteractions(int forceGroups, bool includeForces, bool includeEnergy) {
    if ((forceGroups&groupFlags) == 0)
        return;
    KernelSet& kernels = groupKernels[forceGroups];
    if (kernels.hasForces && (includeForces || includeEnergy)) {
        ComputeKernel& kernel = (includeForces ? (includeEnergy ? kernels.forceEnergyKernel : kernels.forceKernel) : kernels.energyKernel);
        if (!kernel)
            kernel = createInteractionKernel(kernels.source, parameters, arguments, true, true, forceGroups, includeForces, includeEnergy);
        if (useCutoff)
            setPeriodicBoxArgs(context, kernel, 9);
        kernel->execute(numForceThreadBlocks*forceThreadBlockSize, forceThreadBlockSize);
    }
    if (useNeighborList && numTiles > 0) {
        // Keep force work running while waiting only for the earlier count copy.
        context.getCurrentMetalQueue().flush();
        updateNeighborListSize();
    }
}

bool MetalNonbondedUtilities::updateNeighborListSize() {
    if (!useCutoff)
        return false;
    unsigned int counts[2] = {0, 0};
    if (countReadbackPending) {
        const unsigned int* data = static_cast<const unsigned int*>(context.unwrap(interactionCount).finishDownload());
        counts[0] = data[0];
        if (canUsePairList)
            counts[1] = data[1];
        countReadbackPending = false;
    }
    else {
        // Explicit size checks outside the normal force-evaluation sequence.
        interactionCount.download(counts);
    }
    downloadedCount = counts[0];
    if (context.getStepsSinceReorder() == 0 || tilesAfterReorder == 0)
        tilesAfterReorder = downloadedCount;
    else if (context.getStepsSinceReorder() > 25 && downloadedCount > 1.1*tilesAfterReorder)
        context.forceReorder();
    if (downloadedCount <= interactingTiles.getSize() && (!canUsePairList || counts[1] <= singlePairs.getSize()))
        return false;

    // The most recent timestep had too many interactions to fit in the arrays.  Make the arrays bigger to prevent
    // this from happening in the future.

    unsigned int maxTiles = interactingTiles.getSize();
    if (downloadedCount > maxTiles) {
        uint64_t blocks = context.getNumAtomBlocks();
        uint64_t totalTiles = blocks*(blocks+1)/2;
        uint64_t requestedTiles = min((uint64_t) (1.2*downloadedCount), totalTiles);
        if (requestedTiles > (uint64_t) numeric_limits<int>::max()/MetalContext::TileSize)
            throw OpenMMException("The neighbor list exceeds the Metal array index limit");
        maxTiles = (unsigned int) requestedTiles;
        interactingTiles.resize(maxTiles);
        interactingAtoms.resize(MetalContext::TileSize*(size_t) maxTiles);
    }
    if (canUsePairList && counts[1] > singlePairs.getSize()) {
        uint64_t requestedPairs = (uint64_t) (1.2*counts[1]+1);
        if (requestedPairs > (uint64_t) numeric_limits<int>::max())
            throw OpenMMException("The sparse pair list exceeds the Metal array index limit");
        singlePairs.resize((int) requestedPairs);
    }
    for (map<int, KernelSet>::iterator iter = groupKernels.begin(); iter != groupKernels.end(); ++iter) {
        KernelSet& kernels = iter->second;
        if (canUsePairList) {
            for (ComputeKernel* kernel : {&kernels.forceKernel, &kernels.energyKernel, &kernels.forceEnergyKernel}) {
                if (*kernel) {
                    (*kernel)->setArg(18, (unsigned int) singlePairs.getSize());
                    (*kernel)->setArg(19, singlePairs);
                }
            }
            kernels.findInteractingBlocksKernel->setArg(19, (unsigned int) singlePairs.getSize());
            kernels.findInteractingBlocksKernel->setArg(20, singlePairs);
        }
        if (kernels.forceKernel) {
            kernels.forceKernel->setArg(7, interactingTiles);
            kernels.forceKernel->setArg(14, (unsigned int) (maxTiles));
            kernels.forceKernel->setArg(17, interactingAtoms);
        }
        if (kernels.energyKernel) {
            kernels.energyKernel->setArg(7, interactingTiles);
            kernels.energyKernel->setArg(14, (unsigned int) (maxTiles));
            kernels.energyKernel->setArg(17, interactingAtoms);
        }
        if (kernels.forceEnergyKernel) {
            kernels.forceEnergyKernel->setArg(7, interactingTiles);
            kernels.forceEnergyKernel->setArg(14, (unsigned int) (maxTiles));
            kernels.forceEnergyKernel->setArg(17, interactingAtoms);
        }
        kernels.findInteractingBlocksKernel->setArg(6, interactingTiles);
        kernels.findInteractingBlocksKernel->setArg(7, interactingAtoms);
        kernels.findInteractingBlocksKernel->setArg(9, (unsigned int) (maxTiles));
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
    startTileIndex = (int) (startFraction*totalTiles);;
    numTiles = (long long) (endFraction*totalTiles)-startTileIndex;
    if (useCutoff) {
        // We are using a cutoff, and the kernels have already been created.

        for (map<int, KernelSet>::iterator iter = groupKernels.begin(); iter != groupKernels.end(); ++iter) {
            KernelSet& kernels = iter->second;
            if (kernels.forceKernel) {
                kernels.forceKernel->setArg(5, (unsigned int) (startTileIndex));
                kernels.forceKernel->setArg(6, (uint64_t) (numTiles));
            }
            if (kernels.energyKernel) {
                kernels.energyKernel->setArg(5, (unsigned int) (startTileIndex));
                kernels.energyKernel->setArg(6, (uint64_t) (numTiles));
            }
            if (kernels.forceEnergyKernel) {
                kernels.forceEnergyKernel->setArg(5, (unsigned int) (startTileIndex));
                kernels.forceEnergyKernel->setArg(6, (uint64_t) (numTiles));
            }
            kernels.findInteractingBlocksKernel->setArg(10, (unsigned int) (startBlockIndex));
            kernels.findInteractingBlocksKernel->setArg(11, (unsigned int) (numBlocks));
        }
        forceRebuildNeighborList = true;
    }
}

void MetalNonbondedUtilities::createKernelsForGroups(int groups) {
    KernelSet kernels;
    string source;
    for (int i = 0; i < 32; i++) {
        if ((groups&(1u<<i)) != 0) {
            source += groupKernelSource[i];
        }
    }
    kernels.hasForces = (source.size() > 0);
    kernels.source = source;
    if (useCutoff) {
        double paddedCutoff = padCutoff(maxCutoff);
        map<string, string> defines;
#if OPENMM_METAL_FAST_FP16_BOUNDS
        // Half overflow is represented by +Inf, an outward bound. Keep this
        // neighbor program safe even when force kernels select fast math.
        defines["OPENMM_METAL_REQUIRE_SAFE_MATH"] = "1";
#endif
        defines["TILE_SIZE"] = context.intToString(MetalContext::TileSize);
        defines["NUM_ATOMS"] = context.intToString(context.getNumAtoms());
        defines["PADDING"] = context.doubleToString(paddedCutoff-maxCutoff);
        defines["PADDED_CUTOFF"] = context.doubleToString(paddedCutoff);
        defines["PADDED_CUTOFF_SQUARED"] = context.doubleToString(paddedCutoff*paddedCutoff);
        defines["NUM_TILES_WITH_EXCLUSIONS"] = context.intToString(exclusionTiles.getSize());
        defines["NUM_BLOCKS"] = context.intToString(context.getNumAtomBlocks());
        defines["SIMD_WIDTH"] = context.intToString(context.getSIMDWidth());
        if (usePeriodic)
            defines["USE_PERIODIC"] = "1";
        if (context.getBoxIsTriclinic())
            defines["TRICLINIC"] = "1";
        if (useLargeBlocks)
            defines["USE_LARGE_BLOCKS"] = "1";
        defines["MAX_EXCLUSIONS"] = context.intToString(maxExclusions);
        defines["BUFFER_GROUPS"] = "2";
        int binShift = 1;
        while (1<<binShift <= context.getNumAtomBlocks())
            binShift++;
        defines["BIN_SHIFT"] = context.intToString(binShift);
        defines["BLOCK_INDEX_MASK"] = context.intToString((1<<binShift)-1);
        string file = MetalNonbondedSources::neighbors(MetalOpenCLKernelSources::findInteractingBlocks, canUsePairList);
        int groupSize = min(256, context.getMaxThreadBlockSize());
        while (true) {
            defines["GROUP_SIZE"] = context.intToString(groupSize);
            ComputeProgram interactingBlocksProgram = context.compileProgram(file, defines);
            kernels.findBlockBoundsKernel = createKernel(interactingBlocksProgram, "findBlockBounds", 11);
            kernels.findBlockBoundsKernel->setArg(0, (int) (context.getNumAtoms()));
            kernels.findBlockBoundsKernel->setArg(6, context.getPosq());
            kernels.findBlockBoundsKernel->setArg(7, blockCenter);
            kernels.findBlockBoundsKernel->setArg(8, blockBoundingBox);
            kernels.findBlockBoundsKernel->setArg(9, rebuildNeighborList);
            kernels.findBlockBoundsKernel->setArg(10, blockSizeRange);
            kernels.computeSortKeysKernel = createKernel(interactingBlocksProgram, "computeSortKeys", 4);
            kernels.computeSortKeysKernel->setArg(0, blockBoundingBox);
            kernels.computeSortKeysKernel->setArg(1, sortedBlocks);
            kernels.computeSortKeysKernel->setArg(2, blockSizeRange);
            kernels.computeSortKeysKernel->setArg(3, (int) (numBlockSizes));
            kernels.sortBoxDataKernel = createKernel(interactingBlocksProgram, "sortBoxData", useLargeBlocks ? 17 : 10);
            kernels.sortBoxDataKernel->setArg(0, sortedBlocks);
            kernels.sortBoxDataKernel->setArg(1, blockCenter);
            kernels.sortBoxDataKernel->setArg(2, blockBoundingBox);
            kernels.sortBoxDataKernel->setArg(3, sortedBlockCenter);
            kernels.sortBoxDataKernel->setArg(4, sortedBlockBoundingBox);
            kernels.sortBoxDataKernel->setArg(5, context.getPosq());
            kernels.sortBoxDataKernel->setArg(6, oldPositions);
            kernels.sortBoxDataKernel->setArg(7, interactionCount);
            kernels.sortBoxDataKernel->setArg(8, rebuildNeighborList);
            kernels.sortBoxDataKernel->setArg(9, (int) (true));
            if (useLargeBlocks) {
                kernels.sortBoxDataKernel->setArg(10, largeBlockCenter);
                kernels.sortBoxDataKernel->setArg(11, largeBlockBoundingBox);
            }
            kernels.findInteractingBlocksKernel = createKernel(interactingBlocksProgram, "findBlocksWithInteractions", 19+(useLargeBlocks ? 2 : 0)+(canUsePairList ? 2 : 0));
            kernels.findInteractingBlocksKernel->setArg(5, interactionCount);
            kernels.findInteractingBlocksKernel->setArg(6, interactingTiles);
            kernels.findInteractingBlocksKernel->setArg(7, interactingAtoms);
            kernels.findInteractingBlocksKernel->setArg(8, context.getPosq());
            kernels.findInteractingBlocksKernel->setArg(9, (unsigned int) (interactingTiles.getSize()));
            kernels.findInteractingBlocksKernel->setArg(10, (unsigned int) (startBlockIndex));
            kernels.findInteractingBlocksKernel->setArg(11, (unsigned int) (numBlocks));
            kernels.findInteractingBlocksKernel->setArg(12, sortedBlocks);
            kernels.findInteractingBlocksKernel->setArg(13, sortedBlockCenter);
            kernels.findInteractingBlocksKernel->setArg(14, sortedBlockBoundingBox);
            kernels.findInteractingBlocksKernel->setArg(15, exclusionIndices);
            kernels.findInteractingBlocksKernel->setArg(16, exclusionRowIndices);
            kernels.findInteractingBlocksKernel->setArg(17, oldPositions);
            kernels.findInteractingBlocksKernel->setArg(18, rebuildNeighborList);
            if (canUsePairList) {
                kernels.findInteractingBlocksKernel->setArg(19, (unsigned int) singlePairs.getSize());
                kernels.findInteractingBlocksKernel->setArg(20, singlePairs);
            }
            if (useLargeBlocks) {
                kernels.findInteractingBlocksKernel->setArg(canUsePairList ? 21 : 19, largeBlockCenter);
                kernels.findInteractingBlocksKernel->setArg(canUsePairList ? 22 : 20, largeBlockBoundingBox);
            }
            if (kernels.findInteractingBlocksKernel->getMaxBlockSize() < groupSize) {
                // The device can't handle this block size, so reduce it.

                groupSize -= 32;
                if (groupSize < 32)
                    throw OpenMMException("Failed to create findInteractingBlocks kernel");
                continue;
            }
            break;
        }
        interactingBlocksThreadBlockSize = groupSize;
    }
    groupKernels[groups] = kernels;
}

ComputeKernel MetalNonbondedUtilities::createInteractionKernel(const string& source, vector<ComputeParameterInfo>& params, vector<ComputeParameterInfo>& arguments, bool useExclusions, bool isSymmetric, int groups, bool includeForces, bool includeEnergy) {
    // Only the upstream template is adapted. Caller-provided kernels keep their ABI.
    bool shuffle = OPENMM_METAL_FAST_NONBONDED_SHUFFLE && kernelSource == MetalOpenCLKernelSources::nonbonded;
    bool sparsePairs = canUsePairList && useCutoff;
    string sourceTemplate = shuffle ? MetalNonbondedSources::cudaSource() : kernelSource;
    if (sparsePairs && !shuffle)
        MetalNonbondedSources::addOpenCLPairs(sourceTemplate);
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
        args << ", __global ";
        if (param.isConstant())
            args << "const ";
        if (param.getNumComponents() == 3)
            args << param.getComponentType();
        else
            args << param.getType();
        args << "* restrict global_";
        args << param.getName();
    }
    for (ComputeParameterInfo& arg : arguments) {
        args << ", __global ";
        if (arg.isConstant())
            args << "const ";
        args << arg.getType() << "* restrict " << arg.getName();
    }
    if (energyParameterDerivatives.size() > 0)
        args << ", __global mixed* restrict energyParamDerivs";
    replacements["PARAMETER_ARGUMENTS"] = args.str();
    stringstream loadLocal1;
    for (const ComputeParameterInfo& param : params) {
        if (param.getNumComponents() == 1) {
            loadLocal1<<"localData[localAtomIndex]."<<param.getName()<<" = "<<param.getName()<<"1;\n";
        }
        else {
            for (int j = 0; j < param.getNumComponents(); ++j)
                loadLocal1<<"localData[localAtomIndex]."<<param.getName()<<"_"<<suffixes[j]<<" = "<<param.getName()<<"1."<<suffixes[j]<<";\n";
        }
    }
    replacements["LOAD_LOCAL_PARAMETERS_FROM_1"] = loadLocal1.str();
    replacements["DECLARE_LOCAL_PARAMETERS"] = "";
    stringstream loadLocal2;
    for (const ComputeParameterInfo& param : params) {
        if (param.getNumComponents() == 1) {
            loadLocal2<<"localData[localAtomIndex]."<<param.getName()<<" = global_"<<param.getName()<<"[j];\n";
        }
        else {
            if (param.getNumComponents() == 3)
                loadLocal2<<param.getType()<<" temp_"<<param.getName()<<" = make_"<<param.getType()<<"(global_"<<param.getName()<<"[3*j], global_"<<param.getName()<<"[3*j+1], global_"<<param.getName()<<"[3*j+2]);\n";
            else
                loadLocal2<<param.getType()<<" temp_"<<param.getName()<<" = global_"<<param.getName()<<"[j];\n";
            for (int j = 0; j < param.getNumComponents(); ++j)
                loadLocal2<<"localData[localAtomIndex]."<<param.getName()<<"_"<<suffixes[j]<<" = temp_"<<param.getName()<<"."<<suffixes[j]<<";\n";
        }
    }
    replacements["LOAD_LOCAL_PARAMETERS_FROM_GLOBAL"] = loadLocal2.str();
    stringstream load1;
    for (const ComputeParameterInfo& param : params) {
        load1<<param.getType()<<" "<<param.getName()<<"1 = ";
        if (param.getNumComponents() == 3)
            load1<<"make_"<<param.getType()<<"(global_"<<param.getName()<<"[3*atom1], global_"<<param.getName()<<"[3*atom1+1], global_"<<param.getName()<<"[3*atom1+2]);\n";
        else
            load1<<"global_"<<param.getName()<<"[atom1];\n";
    }
    replacements["LOAD_ATOM1_PARAMETERS"] = load1.str();
    stringstream load2j;
    for (const ComputeParameterInfo& param : params) {
        if (param.getNumComponents() == 1) {
            load2j<<param.getType()<<" "<<param.getName()<<"2 = localData[atom2]."<<param.getName()<<";\n";
        }
        else {
            load2j<<param.getType()<<" "<<param.getName()<<"2 = make_"<<param.getType()<<"(";
            for (int j = 0; j < param.getNumComponents(); ++j) {
                if (j > 0)
                    load2j<<", ";
                load2j<<"localData[atom2]."<<param.getName()<<"_"<<suffixes[j];
            }
            load2j<<");\n";
        }
    }
    replacements["LOAD_ATOM2_PARAMETERS"] = load2j.str();
    stringstream clearLocal;
    for (const ComputeParameterInfo& param : params) {
        if (param.getNumComponents() == 1)
            clearLocal<<"localData[localAtomIndex]."<<param.getName()<<" = 0;\n";
        else
            for (int j = 0; j < param.getNumComponents(); ++j)
                clearLocal<<"localData[localAtomIndex]."<<param.getName()<<"_"<<suffixes[j]<<" = 0;\n";
    }
    replacements["CLEAR_LOCAL_PARAMETERS"] = clearLocal.str();
    // Sparse pairs reuse CUDA's direct atom loads even when the tiled algorithm
    // remains OpenCL. Packed three-component parameters retain their 12-byte ABI.
    auto globalParameter = [](const ComputeParameterInfo& param, const string& atom) {
        string base = "global_"+param.getName();
        if (param.getNumComponents() == 3)
            return "make_"+param.getType()+"("+base+"[3*"+atom+"], "+base+"[3*"+atom+"+1], "+base+"[3*"+atom+"+2])";
        return base+"["+atom+"]";
    };
    stringstream load2Global;
    for (const ComputeParameterInfo& param : params)
        load2Global<<param.getType()<<" "<<param.getName()<<"2 = "<<globalParameter(param, "atom2")<<";\n";
    replacements["LOAD_ATOM2_PARAMETERS_FROM_GLOBAL"] = load2Global.str();
    if (shuffle) {
        stringstream broadcast, declare, load, local, clear, rotate;
        for (const string& component : suffixes) {
            broadcast<<"posq2."<<component<<" = real_shfl(shflPosq."<<component<<", j);\n";
            rotate<<"shflPosq."<<component<<" = real_shfl(shflPosq."<<component<<", tgx+1);\n";
            if (component != "w")
                rotate<<"shflForce."<<component<<" = real_shfl(shflForce."<<component<<", tgx+1);\n";
        }
        for (const ComputeParameterInfo& param : params) {
            string name = param.getName();
            broadcast<<param.getType()<<" shfl"<<name<<";\n";
            declare<<param.getType()<<" shfl"<<name<<";\n";
            load<<"shfl"<<name<<" = "<<globalParameter(param, "j")<<";\n";
            local<<param.getType()<<" "<<name<<"2 = shfl"<<name<<";\n";
            clear<<"shfl"<<name<<" = "<<(param.getNumComponents() == 1 ? "0" : "make_"+param.getType()+"(0)")<<";\n";
            for (int j = 0; j < param.getNumComponents(); j++) {
                string component = param.getNumComponents() == 1 ? "" : "."+suffixes[j];
                broadcast<<"shfl"<<name<<component<<" = real_shfl("<<name<<"1"<<component<<", j);\n";
                rotate<<"shfl"<<name<<component<<" = real_shfl(shfl"<<name<<component<<", tgx+1);\n";
            }
        }
        replacements["BROADCAST_WARP_DATA"] = broadcast.str();
        replacements["DECLARE_LOCAL_PARAMETERS"] = declare.str();
        replacements["LOAD_LOCAL_PARAMETERS_FROM_GLOBAL"] = load.str();
        replacements["LOAD_ATOM2_PARAMETERS"] = local.str();
        replacements["CLEAR_LOCAL_PARAMETERS"] = clear.str();
        replacements["SHUFFLE_WARP_DATA"] = rotate.str();
    }
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
    if (sparsePairs)
        defines["USE_SPARSE_PAIRS"] = "1";
    if (useCutoff && context.getSIMDWidth() < 32)
        defines["PRUNE_BY_CUTOFF"] = "1";
    if (includeForces)
        defines["INCLUDE_FORCES"] = "1";
    if (includeEnergy)
        defines["INCLUDE_ENERGY"] = "1";
    defines["THREAD_BLOCK_SIZE"] = context.intToString(forceThreadBlockSize);
    defines["FORCE_WORK_GROUP_SIZE"] = context.intToString(forceThreadBlockSize);
    double maxCutoff = 0.0;
    for (int i = 0; i < 32; i++) {
        if ((groups&(1u<<i)) != 0) {
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
    int numContexts = context.getNumContexts();
    int startExclusionIndex = context.getContextIndex()*numExclusionTiles/numContexts;
    int endExclusionIndex = (context.getContextIndex()+1)*numExclusionTiles/numContexts;
    defines["FIRST_EXCLUSION_TILE"] = context.intToString(startExclusionIndex);
    defines["LAST_EXCLUSION_TILE"] = context.intToString(endExclusionIndex);
    if ((localDataSize/4)%2 == 0)
        defines["PARAMETER_SIZE_IS_EVEN"] = "1";
    ComputeProgram program = context.compileProgram(context.replaceStrings(sourceTemplate, replacements), defines);
    int argumentCount = 7+(useCutoff ? 11 : 0)+(sparsePairs ? 2 : 0)+params.size()+arguments.size()+(energyParameterDerivatives.empty() ? 0 : 1);
    ComputeKernel kernel = createKernel(program, "computeNonbonded", argumentCount);

    // Set arguments to the Kernel.

    int index = 0;
    kernel->setArg(index++, context.getLongForceBuffer());
    kernel->setArg(index++, context.getEnergyBuffer());
    kernel->setArg(index++, context.getPosq());
    kernel->setArg(index++, exclusions);
    kernel->setArg(index++, exclusionTiles);
    kernel->setArg(index++, (unsigned int) startTileIndex);
    kernel->setArg(index++, (uint64_t) (numTiles));
    if (useCutoff) {
        kernel->setArg(index++, interactingTiles);
        kernel->setArg(index++, interactionCount);
        index += 5; // The periodic box size arguments are set when the kernel is executed.
        kernel->setArg(index++, (unsigned int) (interactingTiles.getSize()));
        kernel->setArg(index++, blockCenter);
        kernel->setArg(index++, blockBoundingBox);
        kernel->setArg(index++, interactingAtoms);
        if (sparsePairs) {
            kernel->setArg(index++, (unsigned int) singlePairs.getSize());
            kernel->setArg(index++, singlePairs);
        }
    }
    for (ComputeParameterInfo& param : params)
        kernel->setArg(index++, param.getArray());
    for (ComputeParameterInfo& arg : arguments)
        kernel->setArg(index++, arg.getArray());
    if (energyParameterDerivatives.size() > 0)
        kernel->setArg(index++, context.getEnergyParamDerivBuffer());
    return kernel;
}

void MetalNonbondedUtilities::setKernelSource(const string& source) {
    kernelSource = source;
}
