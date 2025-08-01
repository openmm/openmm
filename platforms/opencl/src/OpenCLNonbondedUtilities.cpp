/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2025 Stanford University and the Authors.      *
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
#include "OpenCLNonbondedUtilities.h"
#include "OpenCLArray.h"
#include "OpenCLContext.h"
#include "OpenCLKernelSources.h"
#include "OpenCLExpressionUtilities.h"
#include <algorithm>
#include <map>
#include <set>
#include <utility>

using namespace OpenMM;
using namespace std;

class OpenCLNonbondedUtilities::BlockSortTrait : public ComputeSortImpl::SortTrait {
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

OpenCLNonbondedUtilities::OpenCLNonbondedUtilities(OpenCLContext& context) : context(context), useCutoff(false), usePeriodic(false), useNeighborList(false), anyExclusions(false), usePadding(true),
        pinnedCountBuffer(NULL), pinnedCountMemory(NULL), forceRebuildNeighborList(true), groupFlags(0), tilesAfterReorder(0) {
    // Decide how many thread blocks and force buffers to use.

    deviceIsCpu = (context.getDevice().getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_CPU);
    if (deviceIsCpu) {
        numForceThreadBlocks = context.getNumThreadBlocks();
        forceThreadBlockSize = 1;
    }
    else if (context.getSIMDWidth() == 32) {
        int blocksPerComputeUnit = 4;
        std::string vendor = context.getDevice().getInfo<CL_DEVICE_VENDOR>();
        if (vendor.size() >= 5 && vendor.substr(0, 5) == "Apple") {
            // 1536 threads per GPU core.
            blocksPerComputeUnit = 6;
        }
        numForceThreadBlocks = blocksPerComputeUnit*context.getDevice().getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
        forceThreadBlockSize = 256;
    }
    else {
        numForceThreadBlocks = context.getNumThreadBlocks();
        forceThreadBlockSize = (context.getSIMDWidth() >= 32 ? OpenCLContext::ThreadBlockSize : 32);
    }
    pinnedCountBuffer = new cl::Buffer(context.getContext(), CL_MEM_ALLOC_HOST_PTR, sizeof(unsigned int));
    pinnedCountMemory = (unsigned int*) context.getQueue().enqueueMapBuffer(*pinnedCountBuffer, CL_TRUE, CL_MAP_READ, 0, sizeof(int));
    
    // When building the neighbor list, we can optionally use large blocks (1024 atoms) to
    // accelerate the process.  This makes building the neighbor list faster, but it prevents
    // us from sorting atom blocks by size, which leads to a slightly less efficient neighbor
    // list.  We guess based on system size which will be faster.

    useLargeBlocks = (!deviceIsCpu && context.getNumAtoms() > 100000);

    std::string vendor = context.getDevice().getInfo<CL_DEVICE_VENDOR>();
    isAMD = !deviceIsCpu && ((vendor.size() >= 3 && vendor.substr(0, 3) == "AMD") || (vendor.size() >= 28 && vendor.substr(0, 28) == "Advanced Micro Devices, Inc."));

    setKernelSource(deviceIsCpu ? OpenCLKernelSources::nonbonded_cpu : OpenCLKernelSources::nonbonded);
}

OpenCLNonbondedUtilities::~OpenCLNonbondedUtilities() {
    if (pinnedCountBuffer != NULL)
        delete pinnedCountBuffer;
}

void OpenCLNonbondedUtilities::addInteraction(bool usesCutoff, bool usesPeriodic, bool usesExclusions, double cutoffDistance,
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
    this->useNeighborList |= ((useNeighborList || deviceIsCpu) && useCutoff);
    groupCutoff[forceGroup] = cutoffDistance;
    groupFlags |= 1<<forceGroup;
    if (kernel.size() > 0) {
        if (groupKernelSource.find(forceGroup) == groupKernelSource.end())
            groupKernelSource[forceGroup] = "";
        map<string, string> replacements;
        replacements["CUTOFF"] = "CUTOFF_"+context.intToString(forceGroup);
        replacements["CUTOFF_SQUARED"] = "CUTOFF_"+context.intToString(forceGroup)+"_SQUARED";
        groupKernelSource[forceGroup] += context.replaceStrings(kernel, replacements)+"\n";
    }
}

void OpenCLNonbondedUtilities::addParameter(ComputeParameterInfo parameter) {
    parameters.push_back(parameter);
}

void OpenCLNonbondedUtilities::addArgument(ComputeParameterInfo parameter) {
    arguments.push_back(parameter);
}

string OpenCLNonbondedUtilities::addEnergyParameterDerivative(const string& param) {
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

void OpenCLNonbondedUtilities::requestExclusions(const vector<vector<int> >& exclusionList) {
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

void OpenCLNonbondedUtilities::initialize(const System& system) {
    if (atomExclusions.size() == 0) {
        // No exclusions were specifically requested, so just mark every atom as not interacting with itself.

        atomExclusions.resize(context.getNumAtoms());
        for (int i = 0; i < (int) atomExclusions.size(); i++)
            atomExclusions[i].push_back(i);
    }

    // Create the list of tiles.

    int numAtomBlocks = context.getNumAtomBlocks();
    int numContexts = context.getPlatformData().contexts.size();
    setAtomBlockRange(context.getContextIndex()/(double) numContexts, (context.getContextIndex()+1)/(double) numContexts);

    // Build a list of tiles that contain exclusions.

    set<pair<int, int> > tilesWithExclusions;
    for (int atom1 = 0; atom1 < (int) atomExclusions.size(); ++atom1) {
        int x = atom1/OpenCLContext::TileSize;
        for (int j = 0; j < (int) atomExclusions[atom1].size(); ++j) {
            int atom2 = atomExclusions[atom1][j];
            int y = atom2/OpenCLContext::TileSize;
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
    vector<cl_uint> exclusionRowIndicesVec(numAtomBlocks+1, 0);
    vector<cl_uint> exclusionIndicesVec;
    for (int i = 0; i < numAtomBlocks; i++) {
        exclusionIndicesVec.insert(exclusionIndicesVec.end(), exclusionBlocksForBlock[i].begin(), exclusionBlocksForBlock[i].end());
        exclusionRowIndicesVec[i+1] = exclusionIndicesVec.size();
    }
    maxExclusions = 0;
    for (int i = 0; i < (int) exclusionBlocksForBlock.size(); i++)
        maxExclusions = (maxExclusions > exclusionBlocksForBlock[i].size() ? maxExclusions : exclusionBlocksForBlock[i].size());
    exclusionIndices.initialize<cl_uint>(context, exclusionIndicesVec.size(), "exclusionIndices");
    exclusionRowIndices.initialize<cl_uint>(context, exclusionRowIndicesVec.size(), "exclusionRowIndices");
    exclusionIndices.upload(exclusionIndicesVec);
    exclusionRowIndices.upload(exclusionRowIndicesVec);

    // Record the exclusion data.

    exclusions.initialize<cl_uint>(context, tilesWithExclusions.size()*OpenCLContext::TileSize, "exclusions");
    cl_uint allFlags = (cl_uint) -1;
    vector<cl_uint> exclusionVec(exclusions.getSize(), allFlags);
    for (int i = 0; i < exclusions.getSize(); ++i)
        exclusionVec[i] = 0xFFFFFFFF;
    for (int atom1 = 0; atom1 < (int) atomExclusions.size(); ++atom1) {
        int x = atom1/OpenCLContext::TileSize;
        int offset1 = atom1-x*OpenCLContext::TileSize;
        for (int j = 0; j < (int) atomExclusions[atom1].size(); ++j) {
            int atom2 = atomExclusions[atom1][j];
            int y = atom2/OpenCLContext::TileSize;
            int offset2 = atom2-y*OpenCLContext::TileSize;
            if (x > y) {
                int index = exclusionTileMap[make_pair(x, y)]*OpenCLContext::TileSize;
                exclusionVec[index+offset1] &= allFlags-(1<<offset2);
            }
            else {
                int index = exclusionTileMap[make_pair(y, x)]*OpenCLContext::TileSize;
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

        int maxTiles = 20*numAtomBlocks;
        if (maxTiles > numTiles)
            maxTiles = numTiles;
        if (maxTiles < 1)
            maxTiles = 1;
        int numAtoms = context.getNumAtoms();
        interactingTiles.initialize<cl_int>(context, maxTiles, "interactingTiles");
        interactingAtoms.initialize<cl_int>(context, OpenCLContext::TileSize*maxTiles, "interactingAtoms");
        interactionCount.initialize<cl_uint>(context, 1, "interactionCount");
        int elementSize = (context.getUseDoublePrecision() ? sizeof(cl_double) : sizeof(cl_float));
        blockCenter.initialize(context, numAtomBlocks, 4*elementSize, "blockCenter");
        blockBoundingBox.initialize(context, numAtomBlocks, 4*elementSize, "blockBoundingBox");
        sortedBlocks.initialize<cl_uint>(context, numAtomBlocks, "sortedBlocks");
        sortedBlockCenter.initialize(context, numAtomBlocks+1, 4*elementSize, "sortedBlockCenter");
        sortedBlockBoundingBox.initialize(context, numAtomBlocks+1, 4*elementSize, "sortedBlockBoundingBox");
        numBlockSizes = min((context.getNumAtomBlocks()+63)/64, context.getNumThreadBlocks());
        blockSizeRange.initialize(context, numBlockSizes, 2*elementSize, "blockSizeRange");
        largeBlockCenter.initialize(context, numAtomBlocks, 4*elementSize, "largeBlockCenter");
        largeBlockBoundingBox.initialize(context, numAtomBlocks, 4*elementSize, "largeBlockBoundingBox");
        oldPositions.initialize(context, numAtoms, 4*elementSize, "oldPositions");
        rebuildNeighborList.initialize<int>(context, 1, "rebuildNeighborList");
        blockSorter = context.createSort(new BlockSortTrait(), numAtomBlocks, false);
        vector<cl_uint> count(1, 0);
        interactionCount.upload(count);
        rebuildNeighborList.upload(count);
    }
}

static void setPeriodicBoxArgs(OpenCLContext& cl, cl::Kernel& kernel, int index) {
    if (cl.getUseDoublePrecision()) {
        kernel.setArg<mm_double4>(index++, cl.getPeriodicBoxSizeDouble());
        kernel.setArg<mm_double4>(index++, cl.getInvPeriodicBoxSizeDouble());
        kernel.setArg<mm_double4>(index++, cl.getPeriodicBoxVecXDouble());
        kernel.setArg<mm_double4>(index++, cl.getPeriodicBoxVecYDouble());
        kernel.setArg<mm_double4>(index, cl.getPeriodicBoxVecZDouble());
    }
    else {
        kernel.setArg<mm_float4>(index++, cl.getPeriodicBoxSize());
        kernel.setArg<mm_float4>(index++, cl.getInvPeriodicBoxSize());
        kernel.setArg<mm_float4>(index++, cl.getPeriodicBoxVecX());
        kernel.setArg<mm_float4>(index++, cl.getPeriodicBoxVecY());
        kernel.setArg<mm_float4>(index, cl.getPeriodicBoxVecZ());
    }
}

double OpenCLNonbondedUtilities::getMaxCutoffDistance() {
    double cutoff = 0.0;
    for (map<int, double>::const_iterator iter = groupCutoff.begin(); iter != groupCutoff.end(); ++iter)
        cutoff = max(cutoff, iter->second);
    return cutoff;
}

double OpenCLNonbondedUtilities::padCutoff(double cutoff) {
    double padding = (usePadding ? 0.1*cutoff : 0.0);
    return cutoff+padding;
}

void OpenCLNonbondedUtilities::prepareInteractions(int forceGroups) {
    if ((forceGroups&groupFlags) == 0)
        return;
    if (groupKernels.find(forceGroups) == groupKernels.end())
        createKernelsForGroups(forceGroups);
    KernelSet& kernels = groupKernels[forceGroups];
    if (useCutoff && usePeriodic) {
        mm_float4 box = context.getPeriodicBoxSize();
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
    context.executeKernel(kernels.findBlockBoundsKernel, context.getNumAtomBlocks());
    context.executeKernel(kernels.computeSortKeysKernel, context.getNumAtomBlocks());
    if (useLargeBlocks)
        setPeriodicBoxArgs(context, kernels.sortBoxDataKernel, 12);
    blockSorter->sort(sortedBlocks);
    kernels.sortBoxDataKernel.setArg<cl_int>(9, forceRebuildNeighborList);
    context.executeKernel(kernels.sortBoxDataKernel, context.getNumAtoms());
    setPeriodicBoxArgs(context, kernels.findInteractingBlocksKernel, 0);
    context.executeKernel(kernels.findInteractingBlocksKernel, context.getNumAtoms(), interactingBlocksThreadBlockSize);
    forceRebuildNeighborList = false;
    context.getQueue().enqueueReadBuffer(interactionCount.getDeviceBuffer(), CL_FALSE, 0, sizeof(int), pinnedCountMemory, NULL, &downloadCountEvent);
    if (isAMD)
        context.getQueue().flush();

    #if __APPLE__ && defined(__aarch64__)
    // Segment the command stream to avoid stalls later.
    if (groupKernels[forceGroups].hasForces)
        context.getQueue().flush();
    #endif
}

void OpenCLNonbondedUtilities::computeInteractions(int forceGroups, bool includeForces, bool includeEnergy) {
    if ((forceGroups&groupFlags) == 0)
        return;
    KernelSet& kernels = groupKernels[forceGroups];
    if (kernels.hasForces) {
        if (isAMD)
            context.getQueue().flush();
        cl::Kernel& kernel = (includeForces ? (includeEnergy ? kernels.forceEnergyKernel : kernels.forceKernel) : kernels.energyKernel);
        if (*reinterpret_cast<cl_kernel*>(&kernel) == NULL)
            kernel = createInteractionKernel(kernels.source, parameters, arguments, true, true, forceGroups, includeForces, includeEnergy);
        if (useCutoff)
            setPeriodicBoxArgs(context, kernel, 9);
        context.executeKernel(kernel, numForceThreadBlocks*forceThreadBlockSize, forceThreadBlockSize);
    }
    if (useNeighborList && numTiles > 0) {
        #if __APPLE__ && defined(__aarch64__)
        // Ensure cached up work executes while you're waiting.
        if (kernels.hasForces)
            context.getQueue().flush();
        #endif
        downloadCountEvent.wait();
        updateNeighborListSize();
    }
}

bool OpenCLNonbondedUtilities::updateNeighborListSize() {
    if (!useCutoff)
        return false;
    if (context.getStepsSinceReorder() == 0 || tilesAfterReorder == 0)
        tilesAfterReorder = pinnedCountMemory[0];
    else if (context.getStepsSinceReorder() > 25 && pinnedCountMemory[0] > 1.1*tilesAfterReorder)
        context.forceReorder();
    if (pinnedCountMemory[0] <= interactingTiles.getSize())
        return false;

    // The most recent timestep had too many interactions to fit in the arrays.  Make the arrays bigger to prevent
    // this from happening in the future.

    unsigned int maxTiles = (unsigned int) (1.2*pinnedCountMemory[0]);
    unsigned int numBlocks = context.getNumAtomBlocks();
    int totalTiles = numBlocks*(numBlocks+1)/2;
    if (maxTiles > totalTiles)
        maxTiles = totalTiles;
    interactingTiles.resize(maxTiles);
    interactingAtoms.resize(OpenCLContext::TileSize*(size_t) maxTiles);
    for (map<int, KernelSet>::iterator iter = groupKernels.begin(); iter != groupKernels.end(); ++iter) {
        KernelSet& kernels = iter->second;
        if (*reinterpret_cast<cl_kernel*>(&kernels.forceKernel) != NULL) {
            kernels.forceKernel.setArg<cl::Buffer>(7, interactingTiles.getDeviceBuffer());
            kernels.forceKernel.setArg<cl_uint>(14, maxTiles);
            kernels.forceKernel.setArg<cl::Buffer>(17, interactingAtoms.getDeviceBuffer());
        }
        if (*reinterpret_cast<cl_kernel*>(&kernels.energyKernel) != NULL) {
            kernels.energyKernel.setArg<cl::Buffer>(7, interactingTiles.getDeviceBuffer());
            kernels.energyKernel.setArg<cl_uint>(14, maxTiles);
            kernels.energyKernel.setArg<cl::Buffer>(17, interactingAtoms.getDeviceBuffer());
        }
        if (*reinterpret_cast<cl_kernel*>(&kernels.forceEnergyKernel) != NULL) {
            kernels.forceEnergyKernel.setArg<cl::Buffer>(7, interactingTiles.getDeviceBuffer());
            kernels.forceEnergyKernel.setArg<cl_uint>(14, maxTiles);
            kernels.forceEnergyKernel.setArg<cl::Buffer>(17, interactingAtoms.getDeviceBuffer());
        }
        kernels.findInteractingBlocksKernel.setArg<cl::Buffer>(6, interactingTiles.getDeviceBuffer());
        kernels.findInteractingBlocksKernel.setArg<cl::Buffer>(7, interactingAtoms.getDeviceBuffer());
        kernels.findInteractingBlocksKernel.setArg<cl_uint>(9, maxTiles);
    }
    forceRebuildNeighborList = true;
    context.setForcesValid(false);
    return true;
}

void OpenCLNonbondedUtilities::setUsePadding(bool padding) {
    usePadding = padding;
}

void OpenCLNonbondedUtilities::setAtomBlockRange(double startFraction, double endFraction) {
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
            if (*reinterpret_cast<cl_kernel*>(&kernels.forceKernel) != NULL) {
                kernels.forceKernel.setArg<cl_uint>(5, startTileIndex);
                kernels.forceKernel.setArg<cl_ulong>(6, numTiles);
            }
            if (*reinterpret_cast<cl_kernel*>(&kernels.energyKernel) != NULL) {
                kernels.energyKernel.setArg<cl_uint>(5, startTileIndex);
                kernels.energyKernel.setArg<cl_ulong>(6, numTiles);
            }
            if (*reinterpret_cast<cl_kernel*>(&kernels.forceEnergyKernel) != NULL) {
                kernels.forceEnergyKernel.setArg<cl_uint>(5, startTileIndex);
                kernels.forceEnergyKernel.setArg<cl_ulong>(6, numTiles);
            }
            kernels.findInteractingBlocksKernel.setArg<cl_uint>(10, startBlockIndex);
            kernels.findInteractingBlocksKernel.setArg<cl_uint>(11, numBlocks);
        }
        forceRebuildNeighborList = true;
    }
}

void OpenCLNonbondedUtilities::createKernelsForGroups(int groups) {
    KernelSet kernels;
    string source;
    for (int i = 0; i < 32; i++) {
        if ((groups&(1<<i)) != 0) {
            source += groupKernelSource[i];
        }
    }
    kernels.hasForces = (source.size() > 0);
    kernels.source = source;
    if (useCutoff) {
        double paddedCutoff = padCutoff(maxCutoff);
        map<string, string> defines;
        defines["TILE_SIZE"] = context.intToString(OpenCLContext::TileSize);
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
        defines["BUFFER_GROUPS"] = (deviceIsCpu ? "4" : "2");
        int binShift = 1;
        while (1<<binShift <= context.getNumAtomBlocks())
            binShift++;
        defines["BIN_SHIFT"] = context.intToString(binShift);
        defines["BLOCK_INDEX_MASK"] = context.intToString((1<<binShift)-1);
        string file = (deviceIsCpu ? OpenCLKernelSources::findInteractingBlocks_cpu : OpenCLKernelSources::findInteractingBlocks);
        int groupSize = (deviceIsCpu || context.getSIMDWidth() < 32 ? 32 : 256);
        while (true) {
            defines["GROUP_SIZE"] = context.intToString(groupSize);
            cl::Program interactingBlocksProgram = context.createProgram(file, defines);
            kernels.findBlockBoundsKernel = cl::Kernel(interactingBlocksProgram, "findBlockBounds");
            kernels.findBlockBoundsKernel.setArg<cl_int>(0, context.getNumAtoms());
            kernels.findBlockBoundsKernel.setArg<cl::Buffer>(6, context.getPosq().getDeviceBuffer());
            kernels.findBlockBoundsKernel.setArg<cl::Buffer>(7, blockCenter.getDeviceBuffer());
            kernels.findBlockBoundsKernel.setArg<cl::Buffer>(8, blockBoundingBox.getDeviceBuffer());
            kernels.findBlockBoundsKernel.setArg<cl::Buffer>(9, rebuildNeighborList.getDeviceBuffer());
            kernels.findBlockBoundsKernel.setArg<cl::Buffer>(10, blockSizeRange.getDeviceBuffer());
            kernels.computeSortKeysKernel = cl::Kernel(interactingBlocksProgram, "computeSortKeys");
            kernels.computeSortKeysKernel.setArg<cl::Buffer>(0, blockBoundingBox.getDeviceBuffer());
            kernels.computeSortKeysKernel.setArg<cl::Buffer>(1, sortedBlocks.getDeviceBuffer());
            kernels.computeSortKeysKernel.setArg<cl::Buffer>(2, blockSizeRange.getDeviceBuffer());
            kernels.computeSortKeysKernel.setArg<cl_int>(3, numBlockSizes);
            kernels.sortBoxDataKernel = cl::Kernel(interactingBlocksProgram, "sortBoxData");
            kernels.sortBoxDataKernel.setArg<cl::Buffer>(0, sortedBlocks.getDeviceBuffer());
            kernels.sortBoxDataKernel.setArg<cl::Buffer>(1, blockCenter.getDeviceBuffer());
            kernels.sortBoxDataKernel.setArg<cl::Buffer>(2, blockBoundingBox.getDeviceBuffer());
            kernels.sortBoxDataKernel.setArg<cl::Buffer>(3, sortedBlockCenter.getDeviceBuffer());
            kernels.sortBoxDataKernel.setArg<cl::Buffer>(4, sortedBlockBoundingBox.getDeviceBuffer());
            kernels.sortBoxDataKernel.setArg<cl::Buffer>(5, context.getPosq().getDeviceBuffer());
            kernels.sortBoxDataKernel.setArg<cl::Buffer>(6, oldPositions.getDeviceBuffer());
            kernels.sortBoxDataKernel.setArg<cl::Buffer>(7, interactionCount.getDeviceBuffer());
            kernels.sortBoxDataKernel.setArg<cl::Buffer>(8, rebuildNeighborList.getDeviceBuffer());
            kernels.sortBoxDataKernel.setArg<cl_int>(9, true);
            if (useLargeBlocks) {
                kernels.sortBoxDataKernel.setArg<cl::Buffer>(10, largeBlockCenter.getDeviceBuffer());
                kernels.sortBoxDataKernel.setArg<cl::Buffer>(11, largeBlockBoundingBox.getDeviceBuffer());
            }
            kernels.findInteractingBlocksKernel = cl::Kernel(interactingBlocksProgram, "findBlocksWithInteractions");
            kernels.findInteractingBlocksKernel.setArg<cl::Buffer>(5, interactionCount.getDeviceBuffer());
            kernels.findInteractingBlocksKernel.setArg<cl::Buffer>(6, interactingTiles.getDeviceBuffer());
            kernels.findInteractingBlocksKernel.setArg<cl::Buffer>(7, interactingAtoms.getDeviceBuffer());
            kernels.findInteractingBlocksKernel.setArg<cl::Buffer>(8, context.getPosq().getDeviceBuffer());
            kernels.findInteractingBlocksKernel.setArg<cl_uint>(9, interactingTiles.getSize());
            kernels.findInteractingBlocksKernel.setArg<cl_uint>(10, startBlockIndex);
            kernels.findInteractingBlocksKernel.setArg<cl_uint>(11, numBlocks);
            kernels.findInteractingBlocksKernel.setArg<cl::Buffer>(12, sortedBlocks.getDeviceBuffer());
            kernels.findInteractingBlocksKernel.setArg<cl::Buffer>(13, sortedBlockCenter.getDeviceBuffer());
            kernels.findInteractingBlocksKernel.setArg<cl::Buffer>(14, sortedBlockBoundingBox.getDeviceBuffer());
            kernels.findInteractingBlocksKernel.setArg<cl::Buffer>(15, exclusionIndices.getDeviceBuffer());
            kernels.findInteractingBlocksKernel.setArg<cl::Buffer>(16, exclusionRowIndices.getDeviceBuffer());
            kernels.findInteractingBlocksKernel.setArg<cl::Buffer>(17, oldPositions.getDeviceBuffer());
            kernels.findInteractingBlocksKernel.setArg<cl::Buffer>(18, rebuildNeighborList.getDeviceBuffer());
            if (useLargeBlocks) {
                kernels.findInteractingBlocksKernel.setArg<cl::Buffer>(19, largeBlockCenter.getDeviceBuffer());
                kernels.findInteractingBlocksKernel.setArg<cl::Buffer>(20, largeBlockBoundingBox.getDeviceBuffer());
            }
            if (kernels.findInteractingBlocksKernel.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(context.getDevice()) < groupSize) {
                // The device can't handle this block size, so reduce it.

                groupSize -= 32;
                if (groupSize < 32)
                    throw OpenMMException("Failed to create findInteractingBlocks kernel");
                continue;
            }
            break;
        }
        interactingBlocksThreadBlockSize = (deviceIsCpu ? 1 : groupSize);
    }
    groupKernels[groups] = kernels;
}

cl::Kernel OpenCLNonbondedUtilities::createInteractionKernel(const string& source, vector<ComputeParameterInfo>& params, vector<ComputeParameterInfo>& arguments, bool useExclusions, bool isSymmetric, int groups, bool includeForces, bool includeEnergy) {
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
        if (context.unwrap(arg.getArray()).getDeviceBuffer().getInfo<CL_MEM_TYPE>() == CL_MEM_OBJECT_IMAGE2D) {
            args << ", __read_only image2d_t ";
            args << arg.getName();
        }
        else {
            if ((context.unwrap(arg.getArray()).getDeviceBuffer().getInfo<CL_MEM_FLAGS>() & CL_MEM_READ_ONLY) == 0) {
                args << ", __global ";
                if (arg.isConstant())
                    args << "const ";
            }
            else
                args << ", __constant ";
            args << arg.getType();
            args << "* restrict ";
            args << arg.getName();
        }
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
            load2j<<param.getType()<<" "<<param.getName()<<"2 = ("<<param.getType()<<") (";
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
    defines["TILE_SIZE"] = context.intToString(OpenCLContext::TileSize);
    int numExclusionTiles = exclusionTiles.getSize();
    defines["NUM_TILES_WITH_EXCLUSIONS"] = context.intToString(numExclusionTiles);
    int numContexts = context.getPlatformData().contexts.size();
    int startExclusionIndex = context.getContextIndex()*numExclusionTiles/numContexts;
    int endExclusionIndex = (context.getContextIndex()+1)*numExclusionTiles/numContexts;
    defines["FIRST_EXCLUSION_TILE"] = context.intToString(startExclusionIndex);
    defines["LAST_EXCLUSION_TILE"] = context.intToString(endExclusionIndex);
    if ((localDataSize/4)%2 == 0)
        defines["PARAMETER_SIZE_IS_EVEN"] = "1";
    cl::Program program = context.createProgram(context.replaceStrings(kernelSource, replacements), defines);
    cl::Kernel kernel(program, "computeNonbonded");

    // Set arguments to the Kernel.

    int index = 0;
    kernel.setArg<cl::Memory>(index++, context.getLongForceBuffer().getDeviceBuffer());
    kernel.setArg<cl::Buffer>(index++, context.getEnergyBuffer().getDeviceBuffer());
    kernel.setArg<cl::Buffer>(index++, context.getPosq().getDeviceBuffer());
    kernel.setArg<cl::Buffer>(index++, exclusions.getDeviceBuffer());
    kernel.setArg<cl::Buffer>(index++, exclusionTiles.getDeviceBuffer());
    kernel.setArg<cl_uint>(index++, startTileIndex);
    kernel.setArg<cl_ulong>(index++, numTiles);
    if (useCutoff) {
        kernel.setArg<cl::Buffer>(index++, interactingTiles.getDeviceBuffer());
        kernel.setArg<cl::Buffer>(index++, interactionCount.getDeviceBuffer());
        index += 5; // The periodic box size arguments are set when the kernel is executed.
        kernel.setArg<cl_uint>(index++, interactingTiles.getSize());
        kernel.setArg<cl::Buffer>(index++, blockCenter.getDeviceBuffer());
        kernel.setArg<cl::Buffer>(index++, blockBoundingBox.getDeviceBuffer());
        kernel.setArg<cl::Buffer>(index++, interactingAtoms.getDeviceBuffer());
    }
    for (ComputeParameterInfo& param : params)
        kernel.setArg<cl::Memory>(index++, context.unwrap(param.getArray()).getDeviceBuffer());
    for (ComputeParameterInfo& arg : arguments)
        kernel.setArg<cl::Memory>(index++, context.unwrap(arg.getArray()).getDeviceBuffer());
    if (energyParameterDerivatives.size() > 0)
        kernel.setArg<cl::Memory>(index++, context.getEnergyParamDerivBuffer().getDeviceBuffer());
    return kernel;
}

void OpenCLNonbondedUtilities::setKernelSource(const string& source) {
    kernelSource = source;
}
