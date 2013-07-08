/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2013 Stanford University and the Authors.      *
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
#include "CudaSort.h"
#include <algorithm>
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


class CudaNonbondedUtilities::BlockSortTrait : public CudaSort::SortTrait {
public:
    BlockSortTrait(bool useDouble) : useDouble(useDouble) {
    }
    int getDataSize() const {return useDouble ? sizeof(double2) : sizeof(float2);}
    int getKeySize() const {return useDouble ? sizeof(double) : sizeof(float);}
    const char* getDataType() const {return "real2";}
    const char* getKeyType() const {return "real";}
    const char* getMinKey() const {return "-3.40282e+38f";}
    const char* getMaxKey() const {return "3.40282e+38f";}
    const char* getMaxValue() const {return "make_real2(3.40282e+38f, 3.40282e+38f)";}
    const char* getSortKey() const {return "value.x";}
private:
    bool useDouble;
};

CudaNonbondedUtilities::CudaNonbondedUtilities(CudaContext& context) : context(context), cutoff(-1.0), useCutoff(false), anyExclusions(false), usePadding(true),
        exclusionIndices(NULL), exclusionRowIndices(NULL), exclusionTiles(NULL), exclusions(NULL), interactingTiles(NULL), interactingAtoms(NULL),
        interactionCount(NULL), blockCenter(NULL), blockBoundingBox(NULL), sortedBlocks(NULL), sortedBlockCenter(NULL), sortedBlockBoundingBox(NULL),
        oldPositions(NULL), rebuildNeighborList(NULL), blockSorter(NULL), nonbondedForceGroup(0) {
    // Decide how many thread blocks to use.

    string errorMessage = "Error initializing nonbonded utilities";
    int multiprocessors;
    CHECK_RESULT(cuDeviceGetAttribute(&multiprocessors, CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT, context.getDevice()));
    numForceThreadBlocks = 4*multiprocessors;
    forceThreadBlockSize = (context.getComputeCapability() < 2.0 ? 128 : 256);
}

CudaNonbondedUtilities::~CudaNonbondedUtilities() {
    if (exclusionIndices != NULL)
        delete exclusionIndices;
    if (exclusionRowIndices != NULL)
        delete exclusionRowIndices;
    if (exclusionTiles != NULL)
        delete exclusionTiles;
    if (exclusions != NULL)
        delete exclusions;
    if (interactingTiles != NULL)
        delete interactingTiles;
    if (interactingAtoms != NULL)
        delete interactingAtoms;
    if (interactionCount != NULL)
        delete interactionCount;
    if (blockCenter != NULL)
        delete blockCenter;
    if (blockBoundingBox != NULL)
        delete blockBoundingBox;
    if (sortedBlocks != NULL)
        delete sortedBlocks;
    if (sortedBlockCenter != NULL)
        delete sortedBlockCenter;
    if (sortedBlockBoundingBox != NULL)
        delete sortedBlockBoundingBox;
    if (oldPositions != NULL)
        delete oldPositions;
    if (rebuildNeighborList != NULL)
        delete rebuildNeighborList;
    if (blockSorter != NULL)
        delete blockSorter;
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
    if (kernel.size() > 0)
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

static bool compareUshort2(ushort2 a, ushort2 b) {
    return ((a.y < b.y) || (a.y == b.y && a.x < b.x));
}

void CudaNonbondedUtilities::initialize(const System& system) {
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
        int x = atom1/CudaContext::TileSize;
        for (int j = 0; j < (int) atomExclusions[atom1].size(); ++j) {
            int atom2 = atomExclusions[atom1][j];
            int y = atom2/CudaContext::TileSize;
            tilesWithExclusions.insert(make_pair(max(x, y), min(x, y)));
        }
    }
    vector<ushort2> exclusionTilesVec;
    for (set<pair<int, int> >::const_iterator iter = tilesWithExclusions.begin(); iter != tilesWithExclusions.end(); ++iter)
        exclusionTilesVec.push_back(make_ushort2((unsigned short) iter->first, (unsigned short) iter->second));
    sort(exclusionTilesVec.begin(), exclusionTilesVec.end(), compareUshort2);
    exclusionTiles = CudaArray::create<ushort2>(context, exclusionTilesVec.size(), "exclusionTiles");
    exclusionTiles->upload(exclusionTilesVec);
    map<pair<int, int>, int> exclusionTileMap;
    for (int i = 0; i < (int) exclusionTilesVec.size(); i++) {
        ushort2 tile = exclusionTilesVec[i];
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
    exclusionIndices = CudaArray::create<unsigned int>(context, exclusionIndicesVec.size(), "exclusionIndices");
    exclusionRowIndices = CudaArray::create<unsigned int>(context, exclusionRowIndicesVec.size(), "exclusionRowIndices");
    exclusionIndices->upload(exclusionIndicesVec);
    exclusionRowIndices->upload(exclusionRowIndicesVec);

    // Record the exclusion data.

    exclusions = CudaArray::create<tileflags>(context, tilesWithExclusions.size()*CudaContext::TileSize, "exclusions");
    tileflags allFlags = (tileflags) -1;
    vector<tileflags> exclusionVec(exclusions->getSize(), allFlags);
    for (int atom1 = 0; atom1 < (int) atomExclusions.size(); ++atom1) {
        int x = atom1/CudaContext::TileSize;
        int offset1 = atom1-x*CudaContext::TileSize;
        for (int j = 0; j < (int) atomExclusions[atom1].size(); ++j) {
            int atom2 = atomExclusions[atom1][j];
            int y = atom2/CudaContext::TileSize;
            int offset2 = atom2-y*CudaContext::TileSize;
            if (x > y) {
                int index = exclusionTileMap[make_pair(x, y)]*CudaContext::TileSize;
                exclusionVec[index+offset1] &= allFlags-(1<<offset2);
            }
            else {
                int index = exclusionTileMap[make_pair(y, x)]*CudaContext::TileSize;
                exclusionVec[index+offset2] &= allFlags-(1<<offset1);
            }
        }
    }
    atomExclusions.clear(); // We won't use this again, so free the memory it used
    exclusions->upload(exclusionVec);

    // Create data structures for the neighbor list.

    if (useCutoff) {
        // Select a size for the arrays that hold the neighbor list.  We have to make a fairly
        // arbitrary guess, but if this turns out to be too small we'll increase it later.

        maxTiles = 20*numAtomBlocks;
        if (maxTiles > numTiles)
            maxTiles = numTiles;
        if (maxTiles < 1)
            maxTiles = 1;
        interactingTiles = CudaArray::create<int>(context, maxTiles, "interactingTiles");
        interactingAtoms = CudaArray::create<int>(context, CudaContext::TileSize*maxTiles, "interactingAtoms");
        interactionCount = CudaArray::create<unsigned int>(context, 1, "interactionCount");
        int elementSize = (context.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
        blockCenter = new CudaArray(context, numAtomBlocks, 4*elementSize, "blockCenter");
        blockBoundingBox = new CudaArray(context, numAtomBlocks, 4*elementSize, "blockBoundingBox");
        sortedBlocks = new CudaArray(context, numAtomBlocks, 2*elementSize, "sortedBlocks");
        sortedBlockCenter = new CudaArray(context, numAtomBlocks+1, 4*elementSize, "sortedBlockCenter");
        sortedBlockBoundingBox = new CudaArray(context, numAtomBlocks+1, 4*elementSize, "sortedBlockBoundingBox");
        oldPositions = new CudaArray(context, numAtoms, 4*elementSize, "oldPositions");
        if (context.getUseDoublePrecision()) {
            vector<double4> oldPositionsVec(numAtoms, make_double4(1e30, 1e30, 1e30, 0));
            oldPositions->upload(oldPositionsVec);
        }
        else {
            vector<float4> oldPositionsVec(numAtoms, make_float4(1e30f, 1e30f, 1e30f, 0));
            oldPositions->upload(oldPositionsVec);
        }
        rebuildNeighborList = CudaArray::create<int>(context, 1, "rebuildNeighborList");
        blockSorter = new CudaSort(context, new BlockSortTrait(context.getUseDoublePrecision()), numAtomBlocks);
        vector<unsigned int> count(1, 0);
        interactionCount->upload(count);
    }

    // Create kernels.

    if (kernelSource.size() > 0)
        forceKernel = createInteractionKernel(kernelSource, parameters, arguments, true, true);
    if (useCutoff) {
        double padding = (usePadding ? 0.1*cutoff : 0.0);
        double paddedCutoff = cutoff+padding;
        map<string, string> defines;
        defines["TILE_SIZE"] = context.intToString(CudaContext::TileSize);
        defines["NUM_BLOCKS"] = context.intToString(context.getNumAtomBlocks());
        defines["NUM_ATOMS"] = context.intToString(context.getNumAtoms());
        defines["PADDING"] = context.doubleToString(padding);
        defines["PADDED_CUTOFF"] = context.doubleToString(paddedCutoff);
        defines["PADDED_CUTOFF_SQUARED"] = context.doubleToString(paddedCutoff*paddedCutoff);
        defines["NUM_TILES_WITH_EXCLUSIONS"] = context.intToString(exclusionTiles->getSize());
        if (usePeriodic)
            defines["USE_PERIODIC"] = "1";
        int maxExclusions = 0;
        for (int i = 0; i < (int) exclusionBlocksForBlock.size(); i++)
            maxExclusions = (maxExclusions > exclusionBlocksForBlock[i].size() ? maxExclusions : exclusionBlocksForBlock[i].size());
        defines["MAX_EXCLUSIONS"] = context.intToString(maxExclusions);
        CUmodule interactingBlocksProgram = context.createModule(CudaKernelSources::vectorOps+CudaKernelSources::findInteractingBlocks, defines);
        findBlockBoundsKernel = context.getKernel(interactingBlocksProgram, "findBlockBounds");
        findBlockBoundsArgs.push_back(&numAtoms);
        findBlockBoundsArgs.push_back(context.getPeriodicBoxSizePointer());
        findBlockBoundsArgs.push_back(context.getInvPeriodicBoxSizePointer());
        findBlockBoundsArgs.push_back(&context.getPosq().getDevicePointer());
        findBlockBoundsArgs.push_back(&blockCenter->getDevicePointer());
        findBlockBoundsArgs.push_back(&blockBoundingBox->getDevicePointer());
        findBlockBoundsArgs.push_back(&rebuildNeighborList->getDevicePointer());
        findBlockBoundsArgs.push_back(&sortedBlocks->getDevicePointer());
        sortBoxDataKernel = context.getKernel(interactingBlocksProgram, "sortBoxData");
        sortBoxDataArgs.push_back(&sortedBlocks->getDevicePointer());
        sortBoxDataArgs.push_back(&blockCenter->getDevicePointer());
        sortBoxDataArgs.push_back(&blockBoundingBox->getDevicePointer());
        sortBoxDataArgs.push_back(&sortedBlockCenter->getDevicePointer());
        sortBoxDataArgs.push_back(&sortedBlockBoundingBox->getDevicePointer());
        sortBoxDataArgs.push_back(&context.getPosq().getDevicePointer());
        sortBoxDataArgs.push_back(&oldPositions->getDevicePointer());
        sortBoxDataArgs.push_back(&interactionCount->getDevicePointer());
        sortBoxDataArgs.push_back(&rebuildNeighborList->getDevicePointer());
        findInteractingBlocksKernel = context.getKernel(interactingBlocksProgram, "findBlocksWithInteractions");
        findInteractingBlocksArgs.push_back(context.getPeriodicBoxSizePointer());
        findInteractingBlocksArgs.push_back(context.getInvPeriodicBoxSizePointer());
        findInteractingBlocksArgs.push_back(&interactionCount->getDevicePointer());
        findInteractingBlocksArgs.push_back(&interactingTiles->getDevicePointer());
        findInteractingBlocksArgs.push_back(&interactingAtoms->getDevicePointer());
        findInteractingBlocksArgs.push_back(&context.getPosq().getDevicePointer());
        findInteractingBlocksArgs.push_back(&maxTiles);
        findInteractingBlocksArgs.push_back(&startBlockIndex);
        findInteractingBlocksArgs.push_back(&numBlocks);
        findInteractingBlocksArgs.push_back(&sortedBlocks->getDevicePointer());
        findInteractingBlocksArgs.push_back(&sortedBlockCenter->getDevicePointer());
        findInteractingBlocksArgs.push_back(&sortedBlockBoundingBox->getDevicePointer());
        findInteractingBlocksArgs.push_back(&exclusionIndices->getDevicePointer());
        findInteractingBlocksArgs.push_back(&exclusionRowIndices->getDevicePointer());
        findInteractingBlocksArgs.push_back(&oldPositions->getDevicePointer());
        findInteractingBlocksArgs.push_back(&rebuildNeighborList->getDevicePointer());
    }
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
    blockSorter->sort(*sortedBlocks);
    context.executeKernel(sortBoxDataKernel, &sortBoxDataArgs[0], context.getNumAtoms());
    context.executeKernel(findInteractingBlocksKernel, &findInteractingBlocksArgs[0], context.getNumAtoms(), 256);
}

void CudaNonbondedUtilities::computeInteractions() {
    if (kernelSource.size() > 0) {
        context.executeKernel(forceKernel, &forceArgs[0], numForceThreadBlocks*forceThreadBlockSize, forceThreadBlockSize);
        if (context.getComputeForceCount() == 1)
            updateNeighborListSize(); // This is the first time step, so check whether our initial guess was large enough.
    }
}

void CudaNonbondedUtilities::updateNeighborListSize() {
    if (!useCutoff)
        return;
    unsigned int* pinnedInteractionCount = (unsigned int*) context.getPinnedBuffer();
    interactionCount->download(pinnedInteractionCount);
    if (pinnedInteractionCount[0] <= (unsigned int) maxTiles)
        return;

    // The most recent timestep had too many interactions to fit in the arrays.  Make the arrays bigger to prevent
    // this from happening in the future.

    maxTiles = (int) (1.2*pinnedInteractionCount[0]);
    int totalTiles = context.getNumAtomBlocks()*(context.getNumAtomBlocks()+1)/2;
    if (maxTiles > totalTiles)
        maxTiles = totalTiles;
    delete interactingTiles;
    delete interactingAtoms;
    interactingTiles = NULL; // Avoid an error in the destructor if the following allocation fails
    interactingAtoms = NULL;
    interactingTiles = CudaArray::create<int>(context, maxTiles, "interactingTiles");
    interactingAtoms = CudaArray::create<int>(context, CudaContext::TileSize*maxTiles, "interactingAtoms");
    if (forceArgs.size() > 0)
        forceArgs[7] = &interactingTiles->getDevicePointer();
    findInteractingBlocksArgs[3] = &interactingTiles->getDevicePointer();
    if (forceArgs.size() > 0)
        forceArgs[14] = &interactingAtoms->getDevicePointer();
    findInteractingBlocksArgs[4] = &interactingAtoms->getDevicePointer();
    if (context.getUseDoublePrecision()) {
        vector<double4> oldPositionsVec(numAtoms, make_double4(1e30, 1e30, 1e30, 0));
        oldPositions->upload(oldPositionsVec);
    }
    else {
        vector<float4> oldPositionsVec(numAtoms, make_float4(1e30f, 1e30f, 1e30f, 0));
        oldPositions->upload(oldPositionsVec);
    }
}

void CudaNonbondedUtilities::setUsePadding(bool padding) {
    usePadding = padding;
}

void CudaNonbondedUtilities::setAtomBlockRange(double startFraction, double endFraction) {
    int numAtomBlocks = context.getNumAtomBlocks();
    startBlockIndex = (int) (startFraction*numAtomBlocks);
    numBlocks = (int) (endFraction*numAtomBlocks)-startBlockIndex;
    int totalTiles = context.getNumAtomBlocks()*(context.getNumAtomBlocks()+1)/2;
    startTileIndex = (int) (startFraction*totalTiles);;
    numTiles = (int) (endFraction*totalTiles)-startTileIndex;
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

    int cudaVersion;
    cuDriverGetVersion(&cudaVersion);
    bool useShuffle = (context.getComputeCapability() >= 3.0 && cudaVersion >= 5050);

    // Part 1. Defines for on diagonal exclusion tiles
    stringstream loadLocal1;
    if(useShuffle) {
        // not needed if using shuffles as we can directly fetch from register
    } else {
        for (int i = 0; i < (int) params.size(); i++) {
            if (params[i].getNumComponents() == 1) {
                loadLocal1<<"localData[threadIdx.x]."<<params[i].getName()<<" = "<<params[i].getName()<<"1;\n";
            }
            else {
                for (int j = 0; j < params[i].getNumComponents(); ++j)
                    loadLocal1<<"localData[threadIdx.x]."<<params[i].getName()<<"_"<<suffixes[j]<<" = "<<params[i].getName()<<"1."<<suffixes[j]<<";\n";
            }
        }
    }
    replacements["LOAD_LOCAL_PARAMETERS_FROM_1"] = loadLocal1.str();

    stringstream broadcastWarpData;
    if(useShuffle) {
        broadcastWarpData << "posq2.x = real_shfl(shflPosq.x, j);\n";
        broadcastWarpData << "posq2.y = real_shfl(shflPosq.y, j);\n";
        broadcastWarpData << "posq2.z = real_shfl(shflPosq.z, j);\n";
        broadcastWarpData << "posq2.w = real_shfl(shflPosq.w, j);\n";
        for(int i=0; i< (int) params.size();i++) {
            broadcastWarpData << params[i].getType() << " shfl" << params[i].getName() << ";\n";
            for(int j=0; j < params[i].getNumComponents(); j++) {
                string name;
                if (params[i].getNumComponents() == 1) {
                    broadcastWarpData << "shfl" << params[i].getName() << "=real_shfl(" << params[i].getName() <<"1,j);\n";

                } else {
                    broadcastWarpData << "shfl" << params[i].getName()+"."+suffixes[j] << "=real_shfl(" << params[i].getName()+"1."+suffixes[j] <<",j);\n";
                }
            }
        }
    } else {
        // not used if not shuffling
    }
    replacements["BROADCAST_WARP_DATA"] = broadcastWarpData.str();
    
    // Part 2. Defines for off-diagonal exclusions, and neighborlist tiles. 
    stringstream declareLocal2;
    if(useShuffle) {
        for(int i=0; i< (int) params.size(); i++) {
            declareLocal2<<params[i].getType()<<" shfl"<<params[i].getName()<<";\n";
        }
    } else {
        // not used if using shared memory
    }
    replacements["DECLARE_LOCAL_PARAMETERS"] = declareLocal2.str();

    stringstream loadLocal2;
    if(useShuffle) {
        for(int i=0; i< (int) params.size(); i++) {
            loadLocal2<<"shfl"<<params[i].getName()<<" = global_"<<params[i].getName()<<"[j];\n";
        }
    } else {
        for (int i = 0; i < (int) params.size(); i++) {
            if (params[i].getNumComponents() == 1) {
                loadLocal2<<"localData[threadIdx.x]."<<params[i].getName()<<" = global_"<<params[i].getName()<<"[j];\n";
            }
            else {
                loadLocal2<<params[i].getType()<<" temp_"<<params[i].getName()<<" = global_"<<params[i].getName()<<"[j];\n";
                for (int j = 0; j < params[i].getNumComponents(); ++j)
                    loadLocal2<<"localData[threadIdx.x]."<<params[i].getName()<<"_"<<suffixes[j]<<" = temp_"<<params[i].getName()<<"."<<suffixes[j]<<";\n";
            }
        }
    }
    replacements["LOAD_LOCAL_PARAMETERS_FROM_GLOBAL"] = loadLocal2.str();
   
    stringstream load2j;
    if(useShuffle) {
        for(int i = 0; i < (int) params.size(); i++)
            load2j<<params[i].getType()<<" "<<params[i].getName()<<"2 = shfl"<<params[i].getName()<<";\n";
    } else {
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
    }
    replacements["LOAD_ATOM2_PARAMETERS"] = load2j.str();

    stringstream shuffleWarpData;
    if(useShuffle) {
        shuffleWarpData << "shflPosq.x = real_shfl(shflPosq.x, tgx+1);\n";
        shuffleWarpData << "shflPosq.y = real_shfl(shflPosq.y, tgx+1);\n";
        shuffleWarpData << "shflPosq.z = real_shfl(shflPosq.z, tgx+1);\n";
        shuffleWarpData << "shflPosq.w = real_shfl(shflPosq.w, tgx+1);\n";
        shuffleWarpData << "shflForce.x = real_shfl(shflForce.x, tgx+1);\n";
        shuffleWarpData << "shflForce.y = real_shfl(shflForce.y, tgx+1);\n";
        shuffleWarpData << "shflForce.z = real_shfl(shflForce.z, tgx+1);\n";
        for(int i=0; i < (int) params.size(); i++) {
            if(params[i].getNumComponents() == 1) {
                shuffleWarpData<<"shfl"<<params[i].getName()<<"=real_shfl(shfl"<<params[i].getName()<<", tgx+1);\n";
            } else {
                for(int j=0;j<params[i].getNumComponents();j++) {
                    // looks something like shflsigmaEpsilon.x = real_shfl(shflsigmaEpsilon.x,tgx+1);
                    shuffleWarpData<<"shfl"<<params[i].getName()
                        <<"."<<suffixes[j]<<"=real_shfl(shfl"
                        <<params[i].getName()<<"."<<suffixes[j]
                        <<", tgx+1);\n";
                }
            }
        }
    } else {
        // not used otherwise
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
    if (useShuffle)
        defines["ENABLE_SHUFFLE"] = "1";
    defines["THREAD_BLOCK_SIZE"] = context.intToString(forceThreadBlockSize);
    defines["CUTOFF_SQUARED"] = context.doubleToString(cutoff*cutoff);
    defines["CUTOFF"] = context.doubleToString(cutoff);
    defines["NUM_ATOMS"] = context.intToString(context.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = context.intToString(context.getPaddedNumAtoms());
    defines["NUM_BLOCKS"] = context.intToString(context.getNumAtomBlocks());
    defines["TILE_SIZE"] = context.intToString(CudaContext::TileSize);
    int numExclusionTiles = exclusionTiles->getSize();
    defines["NUM_TILES_WITH_EXCLUSIONS"] = context.intToString(numExclusionTiles);
    int numContexts = context.getPlatformData().contexts.size();
    int startExclusionIndex = context.getContextIndex()*numExclusionTiles/numContexts;
    int endExclusionIndex = (context.getContextIndex()+1)*numExclusionTiles/numContexts;
    defines["FIRST_EXCLUSION_TILE"] = context.intToString(startExclusionIndex);
    defines["LAST_EXCLUSION_TILE"] = context.intToString(endExclusionIndex);
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
    forceArgs.push_back(&exclusionTiles->getDevicePointer());
    forceArgs.push_back(&startTileIndex);
    forceArgs.push_back(&numTiles);
    if (useCutoff) {
        forceArgs.push_back(&interactingTiles->getDevicePointer());
        forceArgs.push_back(&interactionCount->getDevicePointer());
        forceArgs.push_back(context.getPeriodicBoxSizePointer());
        forceArgs.push_back(context.getInvPeriodicBoxSizePointer());
        forceArgs.push_back(&maxTiles);
        forceArgs.push_back(&blockCenter->getDevicePointer());
        forceArgs.push_back(&blockBoundingBox->getDevicePointer());
        forceArgs.push_back(&interactingAtoms->getDevicePointer());
    }
    for (int i = 0; i < (int) params.size(); i++)
        forceArgs.push_back(&params[i].getMemory());
    for (int i = 0; i < (int) arguments.size(); i++)
        forceArgs.push_back(&arguments[i].getMemory());
    return kernel;
}
