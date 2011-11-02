/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2011 Stanford University and the Authors.      *
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

#include "OpenCLNonbondedUtilities.h"
#include "OpenCLArray.h"
#include "OpenCLKernelSources.h"
#include "OpenCLExpressionUtilities.h"
#include <map>
#include <set>
#include <utility>

using namespace OpenMM;
using namespace std;

OpenCLNonbondedUtilities::OpenCLNonbondedUtilities(OpenCLContext& context) : context(context), cutoff(-1.0), useCutoff(false), anyExclusions(false),
        numForceBuffers(0), exclusionIndices(NULL), exclusionRowIndices(NULL), exclusions(NULL), interactingTiles(NULL), interactionFlags(NULL),
        interactionCount(NULL), blockCenter(NULL), blockBoundingBox(NULL) {
    // Decide how many thread blocks and force buffers to use.

    deviceIsCpu = (context.getDevice().getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_CPU);
    forceBufferPerAtomBlock = false;
    if (deviceIsCpu) {
        numForceThreadBlocks = context.getNumThreadBlocks();
        forceThreadBlockSize = 1;
        numForceBuffers = numForceThreadBlocks;
    }
    else if (context.getSIMDWidth() == 32) {
        if (context.getSupports64BitGlobalAtomics()) {
            numForceThreadBlocks = 2*context.getDevice().getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
            forceThreadBlockSize = 256;
            numForceBuffers = 2;
        }
        else {
            numForceThreadBlocks = 4*context.getDevice().getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
            forceThreadBlockSize = 128;
            numForceBuffers = numForceThreadBlocks;
        }
    }
    else {
        numForceThreadBlocks = context.getNumThreadBlocks();
        forceThreadBlockSize = OpenCLContext::ThreadBlockSize;
        numForceBuffers = numForceThreadBlocks;
        if (numForceBuffers >= context.getNumAtomBlocks()) {
            // For small systems, it is more efficient to have one force buffer per block of 32 atoms instead of one per warp.

            forceBufferPerAtomBlock = true;
            numForceBuffers = context.getNumAtomBlocks();
        }
    }
}

OpenCLNonbondedUtilities::~OpenCLNonbondedUtilities() {
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
}

void OpenCLNonbondedUtilities::addInteraction(bool usesCutoff, bool usesPeriodic, bool usesExclusions, double cutoffDistance, const vector<vector<int> >& exclusionList, const string& kernel) {
    if (cutoff != -1.0) {
        if (usesCutoff != useCutoff)
            throw OpenMMException("All Forces must agree on whether to use a cutoff");
        if (usesPeriodic != usePeriodic)
            throw OpenMMException("All Forces must agree on whether to use periodic boundary conditions");
        if (cutoffDistance != cutoff)
            throw OpenMMException("All Forces must use the same cutoff distance");
    }
    if (usesExclusions)
        requestExclusions(exclusionList);
    useCutoff = usesCutoff;
    usePeriodic = usesPeriodic;
    cutoff = cutoffDistance;
    kernelSource += kernel+"\n";
}

void OpenCLNonbondedUtilities::addParameter(const ParameterInfo& parameter) {
    parameters.push_back(parameter);
}

void OpenCLNonbondedUtilities::addArgument(const ParameterInfo& parameter) {
    arguments.push_back(parameter);
}

void OpenCLNonbondedUtilities::requestExclusions(const vector<vector<int> >& exclusionList) {
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

void OpenCLNonbondedUtilities::initialize(const System& system) {
    if (cutoff == -1.0)
        return; // There are no nonbonded interactions in the System.
    
    if (atomExclusions.size() == 0) {
        // No exclusions were specifically requested, so just mark every atom as not interacting with itself.
        
        atomExclusions.resize(context.getNumAtoms());
        for (int i = 0; i < (int) atomExclusions.size(); i++)
            atomExclusions[i].push_back(i);
    }

    // Create the list of tiles.

    int numAtomBlocks = context.getNumAtomBlocks();
    int totalTiles = numAtomBlocks*(numAtomBlocks+1)/2;
    int numContexts = context.getPlatformData().contexts.size();
    startTileIndex = context.getContextIndex()*totalTiles/numContexts;
    int endTileIndex = (context.getContextIndex()+1)*totalTiles/numContexts;
    numTiles = endTileIndex-startTileIndex;

    // Build a list of indices for the tiles with exclusions.

    set<pair<int, int> > tilesWithExclusions;
    for (int atom1 = 0; atom1 < (int) atomExclusions.size(); ++atom1) {
        int x = atom1/OpenCLContext::TileSize;
        for (int j = 0; j < (int) atomExclusions[atom1].size(); ++j) {
            int atom2 = atomExclusions[atom1][j];
            int y = atom2/OpenCLContext::TileSize;
            tilesWithExclusions.insert(make_pair(max(x, y), min(x, y)));
        }
    }
    if (context.getPaddedNumAtoms() > context.getNumAtoms()) {
        for (int i = 0; i < numAtomBlocks; ++i)
            tilesWithExclusions.insert(make_pair(numAtomBlocks-1, i));
    }
    vector<cl_uint> exclusionRowIndicesVec(numAtomBlocks+1, 0);
    vector<cl_uint> exclusionIndicesVec;
    int currentRow = 0;
    for (set<pair<int, int> >::const_iterator iter = tilesWithExclusions.begin(); iter != tilesWithExclusions.end(); ++iter) {
        while (iter->first != currentRow)
            exclusionRowIndicesVec[++currentRow] = exclusionIndicesVec.size();
        exclusionIndicesVec.push_back(iter->second);
    }
    exclusionRowIndicesVec[++currentRow] = exclusionIndicesVec.size();
    exclusionIndices = new OpenCLArray<cl_uint>(context, exclusionIndicesVec.size(), "exclusionIndices");
    exclusionRowIndices = new OpenCLArray<cl_uint>(context, exclusionRowIndicesVec.size(), "exclusionRowIndices");
    exclusionIndices->upload(exclusionIndicesVec);
    exclusionRowIndices->upload(exclusionRowIndicesVec);

    // Record the exclusion data.

    exclusions = new OpenCLArray<cl_uint>(context, tilesWithExclusions.size()*OpenCLContext::TileSize, "exclusions");
    vector<cl_uint> exclusionVec(exclusions->getSize());
    for (int i = 0; i < exclusions->getSize(); ++i)
        exclusionVec[i] = 0xFFFFFFFF;
    for (int atom1 = 0; atom1 < (int) atomExclusions.size(); ++atom1) {
        int x = atom1/OpenCLContext::TileSize;
        int offset1 = atom1-x*OpenCLContext::TileSize;
        for (int j = 0; j < (int) atomExclusions[atom1].size(); ++j) {
            int atom2 = atomExclusions[atom1][j];
            int y = atom2/OpenCLContext::TileSize;
            int offset2 = atom2-y*OpenCLContext::TileSize;
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
        int x = atom1/OpenCLContext::TileSize;
        int offset1 = atom1-x*OpenCLContext::TileSize;
        for (int atom2 = 0; atom2 < context.getPaddedNumAtoms(); ++atom2) {
            int y = atom2/OpenCLContext::TileSize;
            int offset2 = atom2-y*OpenCLContext::TileSize;
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

        mm_float4 boxSize = context.getPeriodicBoxSize();
        int maxInteractingTiles = (int) (numTiles*(cutoff/boxSize.x+cutoff/boxSize.y+cutoff/boxSize.z));
        if (maxInteractingTiles > numTiles)
            maxInteractingTiles = numTiles;
        if (maxInteractingTiles < 1)
            maxInteractingTiles = 1;
        interactingTiles = new OpenCLArray<mm_ushort2>(context, maxInteractingTiles, "interactingTiles");
        interactionFlags = new OpenCLArray<cl_uint>(context, context.getSIMDWidth() == 32 ? maxInteractingTiles : (deviceIsCpu ? 2*maxInteractingTiles : 1), "interactionFlags");
        interactionCount = new OpenCLArray<cl_uint>(context, 1, "interactionCount", true);
        blockCenter = new OpenCLArray<mm_float4>(context, numAtomBlocks, "blockCenter");
        blockBoundingBox = new OpenCLArray<mm_float4>(context, numAtomBlocks, "blockBoundingBox");
        interactionCount->set(0, 0);
        interactionCount->upload();
    }

    // Create kernels.

    forceKernel = createInteractionKernel(kernelSource, parameters, arguments, true, true);
    if (useCutoff) {
        map<string, string> defines;
        defines["NUM_BLOCKS"] = OpenCLExpressionUtilities::intToString(context.getNumAtomBlocks());
        if (forceBufferPerAtomBlock)
            defines["USE_OUTPUT_BUFFER_PER_BLOCK"] = "1";
        if (usePeriodic)
            defines["USE_PERIODIC"] = "1";
        string file = (deviceIsCpu ? OpenCLKernelSources::findInteractingBlocks_cpu : OpenCLKernelSources::findInteractingBlocks);
        cl::Program interactingBlocksProgram = context.createProgram(file, defines);
        findBlockBoundsKernel = cl::Kernel(interactingBlocksProgram, "findBlockBounds");
        findBlockBoundsKernel.setArg<cl_int>(0, context.getNumAtoms());
        findBlockBoundsKernel.setArg<cl::Buffer>(3, context.getPosq().getDeviceBuffer());
        findBlockBoundsKernel.setArg<cl::Buffer>(4, blockCenter->getDeviceBuffer());
        findBlockBoundsKernel.setArg<cl::Buffer>(5, blockBoundingBox->getDeviceBuffer());
        findBlockBoundsKernel.setArg<cl::Buffer>(6, interactionCount->getDeviceBuffer());
        findInteractingBlocksKernel = cl::Kernel(interactingBlocksProgram, "findBlocksWithInteractions");
        findInteractingBlocksKernel.setArg<cl_float>(0, (cl_float) (cutoff*cutoff));
        findInteractingBlocksKernel.setArg<cl::Buffer>(3, blockCenter->getDeviceBuffer());
        findInteractingBlocksKernel.setArg<cl::Buffer>(4, blockBoundingBox->getDeviceBuffer());
        findInteractingBlocksKernel.setArg<cl::Buffer>(5, interactionCount->getDeviceBuffer());
        findInteractingBlocksKernel.setArg<cl::Buffer>(6, interactingTiles->getDeviceBuffer());
        findInteractingBlocksKernel.setArg<cl::Buffer>(7, interactionFlags->getDeviceBuffer());
        findInteractingBlocksKernel.setArg<cl::Buffer>(8, context.getPosq().getDeviceBuffer());
        findInteractingBlocksKernel.setArg<cl_uint>(9, interactingTiles->getSize());
        findInteractingBlocksKernel.setArg<cl_uint>(10, startTileIndex);
        findInteractingBlocksKernel.setArg<cl_uint>(11, startTileIndex+numTiles);
        if (context.getSIMDWidth() == 32 && !deviceIsCpu) {
            findInteractionsWithinBlocksKernel = cl::Kernel(interactingBlocksProgram, "findInteractionsWithinBlocks");
            findInteractionsWithinBlocksKernel.setArg<cl_float>(0, (cl_float) (cutoff*cutoff));
            findInteractionsWithinBlocksKernel.setArg<cl::Buffer>(3, context.getPosq().getDeviceBuffer());
            findInteractionsWithinBlocksKernel.setArg<cl::Buffer>(4, interactingTiles->getDeviceBuffer());
            findInteractionsWithinBlocksKernel.setArg<cl::Buffer>(5, blockCenter->getDeviceBuffer());
            findInteractionsWithinBlocksKernel.setArg<cl::Buffer>(6, blockBoundingBox->getDeviceBuffer());
            findInteractionsWithinBlocksKernel.setArg<cl::Buffer>(7, interactionFlags->getDeviceBuffer());
            findInteractionsWithinBlocksKernel.setArg<cl::Buffer>(8, interactionCount->getDeviceBuffer());
            findInteractionsWithinBlocksKernel.setArg(9, 128*sizeof(cl_uint), NULL);
            findInteractionsWithinBlocksKernel.setArg<cl_uint>(10, interactingTiles->getSize());
        }
    }
}

int OpenCLNonbondedUtilities::findExclusionIndex(int x, int y, const vector<cl_uint>& exclusionIndices, const vector<cl_uint>& exclusionRowIndices) {
    int start = exclusionRowIndices[x];
    int end = exclusionRowIndices[x+1];
    for (int i = start; i < end; i++)
        if (exclusionIndices[i] == y)
            return i*OpenCLContext::TileSize;
    throw OpenMMException("Internal error: exclusion in unexpected tile");
}

void OpenCLNonbondedUtilities::prepareInteractions() {
    if (!useCutoff)
        return;

    // Compute the neighbor list.

    findBlockBoundsKernel.setArg<mm_float4>(1, context.getPeriodicBoxSize());
    findBlockBoundsKernel.setArg<mm_float4>(2, context.getInvPeriodicBoxSize());
    context.executeKernel(findBlockBoundsKernel, context.getNumAtoms());
    findInteractingBlocksKernel.setArg<mm_float4>(1, context.getPeriodicBoxSize());
    findInteractingBlocksKernel.setArg<mm_float4>(2, context.getInvPeriodicBoxSize());
    context.executeKernel(findInteractingBlocksKernel, context.getNumAtoms(), deviceIsCpu ? 1 : -1);
    if (context.getSIMDWidth() == 32 && !deviceIsCpu) {
        findInteractionsWithinBlocksKernel.setArg<mm_float4>(1, context.getPeriodicBoxSize());
        findInteractionsWithinBlocksKernel.setArg<mm_float4>(2, context.getInvPeriodicBoxSize());
        context.executeKernel(findInteractionsWithinBlocksKernel, context.getNumAtoms(), 128);
    }
}

void OpenCLNonbondedUtilities::computeInteractions() {
    if (cutoff != -1.0) {
        if (useCutoff) {
            forceKernel.setArg<mm_float4>(11, context.getPeriodicBoxSize());
            forceKernel.setArg<mm_float4>(12, context.getInvPeriodicBoxSize());
        }
        context.executeKernel(forceKernel, numForceThreadBlocks*forceThreadBlockSize, forceThreadBlockSize);
    }
}

void OpenCLNonbondedUtilities::updateNeighborListSize() {
    if (!useCutoff)
        return;
    interactionCount->download();
    if (interactionCount->get(0) <= (unsigned int) interactingTiles->getSize())
        return;

    // The most recent timestep had too many interactions to fit in the arrays.  Make the arrays bigger to prevent
    // this from happening in the future.

    int newSize = (int) (1.2*interactionCount->get(0));
    int numTiles = context.getNumAtomBlocks()*(context.getNumAtomBlocks()+1)/2;
    if (newSize > numTiles)
        newSize = numTiles;
    delete interactingTiles;
    interactingTiles = new OpenCLArray<mm_ushort2>(context, newSize, "interactingTiles");
    forceKernel.setArg<cl::Buffer>(9, interactingTiles->getDeviceBuffer());
    forceKernel.setArg<cl_uint>(13, newSize);
    findInteractingBlocksKernel.setArg<cl::Buffer>(6, interactingTiles->getDeviceBuffer());
    findInteractingBlocksKernel.setArg<cl_uint>(9, newSize);
    if (context.getSIMDWidth() == 32 || deviceIsCpu) {
        delete interactionFlags;
        interactionFlags = new OpenCLArray<cl_uint>(context, deviceIsCpu ? 2*newSize : newSize, "interactionFlags");
        forceKernel.setArg<cl::Buffer>(14, interactionFlags->getDeviceBuffer());
        findInteractingBlocksKernel.setArg<cl::Buffer>(7, interactionFlags->getDeviceBuffer());
        findInteractionsWithinBlocksKernel.setArg<cl::Buffer>(4, interactingTiles->getDeviceBuffer());
        findInteractionsWithinBlocksKernel.setArg<cl::Buffer>(7, interactionFlags->getDeviceBuffer());
        findInteractionsWithinBlocksKernel.setArg<cl_uint>(10, newSize);
    }
}

void OpenCLNonbondedUtilities::setTileRange(int startTileIndex, int numTiles) {
    this->startTileIndex = startTileIndex;
    this->numTiles = numTiles;
    if (cutoff == -1.0)
        return; // There are no nonbonded interactions in the System.
    forceKernel.setArg<cl_uint>(7, startTileIndex);
    forceKernel.setArg<cl_uint>(8, startTileIndex+numTiles);
    if (useCutoff) {
        findInteractingBlocksKernel.setArg<cl_uint>(10, startTileIndex);
        findInteractingBlocksKernel.setArg<cl_uint>(11, startTileIndex+numTiles);
    }
}

cl::Kernel OpenCLNonbondedUtilities::createInteractionKernel(const string& source, const vector<ParameterInfo>& params, const vector<ParameterInfo>& arguments, bool useExclusions, bool isSymmetric) const {
    map<string, string> replacements;
    replacements["COMPUTE_INTERACTION"] = source;
    int localDataSize = 7*sizeof(cl_float);
    const string suffixes[] = {"x", "y", "z", "w"};
    stringstream localData;
    for (int i = 0; i < (int) params.size(); i++) {
        if (params[i].getNumComponents() == 1)
            localData<<params[i].getType()<<" "<<params[i].getName()<<";\n";
        else {
            for (int j = 0; j < params[i].getNumComponents(); ++j)
                localData<<params[i].getComponentType()<<" "<<params[i].getName()<<"_"<<suffixes[j]<<";\n";
        }
        localDataSize += params[i].getSize();
    }
    if ((localDataSize/4)%2 == 0) {
        localData << "float padding;\n";
        localDataSize += 4;
    }
    replacements["ATOM_PARAMETER_DATA"] = localData.str();
    stringstream args;
    for (int i = 0; i < (int) params.size(); i++) {
        args << ", __global const ";
        args << params[i].getType();
        args << "* restrict global_";
        args << params[i].getName();
    }
    for (int i = 0; i < (int) arguments.size(); i++) {
        if (arguments[i].getMemory().getInfo<CL_MEM_TYPE>() == CL_MEM_OBJECT_IMAGE2D) {
            args << ", __read_only image2d_t ";
            args << arguments[i].getName();
        }
        else {
            if ((arguments[i].getMemory().getInfo<CL_MEM_FLAGS>() & CL_MEM_READ_ONLY) == 0)
                args << ", __global const ";
            else
                args << ", __constant ";
            args << arguments[i].getType();
            args << "* restrict ";
            args << arguments[i].getName();
        }
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
            load2j<<params[i].getType()<<" "<<params[i].getName()<<"2 = ("<<params[i].getType()<<") (";
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
    if (forceBufferPerAtomBlock)
        defines["USE_OUTPUT_BUFFER_PER_BLOCK"] = "1";
    if (useCutoff)
        defines["USE_CUTOFF"] = "1";
    if (usePeriodic)
        defines["USE_PERIODIC"] = "1";
    if (useExclusions)
        defines["USE_EXCLUSIONS"] = "1";
    if (isSymmetric)
        defines["USE_SYMMETRIC"] = "1";
    defines["FORCE_WORK_GROUP_SIZE"] = OpenCLExpressionUtilities::intToString(forceThreadBlockSize);
    defines["CUTOFF_SQUARED"] = OpenCLExpressionUtilities::doubleToString(cutoff*cutoff);
    defines["NUM_ATOMS"] = OpenCLExpressionUtilities::intToString(context.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = OpenCLExpressionUtilities::intToString(context.getPaddedNumAtoms());
    defines["NUM_BLOCKS"] = OpenCLExpressionUtilities::intToString(context.getNumAtomBlocks());
    string file;
    if (deviceIsCpu)
        file = OpenCLKernelSources::nonbonded_cpu;
    else if (context.getSIMDWidth() == 32)
        file = OpenCLKernelSources::nonbonded_nvidia;
    else
        file = OpenCLKernelSources::nonbonded_default;
    cl::Program program = context.createProgram(context.replaceStrings(file, replacements), defines);
    cl::Kernel kernel(program, "computeNonbonded");

    // Set arguments to the Kernel.

    int index = 0;
    if (context.getSupports64BitGlobalAtomics())
        kernel.setArg<cl::Memory>(index++, context.getLongForceBuffer().getDeviceBuffer());
    else
        kernel.setArg<cl::Buffer>(index++, context.getForceBuffers().getDeviceBuffer());
    kernel.setArg<cl::Buffer>(index++, context.getEnergyBuffer().getDeviceBuffer());
    kernel.setArg<cl::Buffer>(index++, context.getPosq().getDeviceBuffer());
    kernel.setArg<cl::Buffer>(index++, exclusions->getDeviceBuffer());
    kernel.setArg<cl::Buffer>(index++, exclusionIndices->getDeviceBuffer());
    kernel.setArg<cl::Buffer>(index++, exclusionRowIndices->getDeviceBuffer());
    kernel.setArg(index++, (deviceIsCpu ? OpenCLContext::TileSize*localDataSize : forceThreadBlockSize*localDataSize), NULL);
    kernel.setArg<cl_uint>(index++, startTileIndex);
    kernel.setArg<cl_uint>(index++, startTileIndex+numTiles);
    if (useCutoff) {
        kernel.setArg<cl::Buffer>(index++, interactingTiles->getDeviceBuffer());
        kernel.setArg<cl::Buffer>(index++, interactionCount->getDeviceBuffer());
        index += 2; // The periodic box size arguments are set when the kernel is executed.
        kernel.setArg<cl_uint>(index++, interactingTiles->getSize());
        kernel.setArg<cl::Buffer>(index++, interactionFlags->getDeviceBuffer());
    }
    else {
        kernel.setArg<cl_uint>(index++, numTiles);
    }
    for (int i = 0; i < (int) params.size(); i++) {
        kernel.setArg<cl::Memory>(index++, params[i].getMemory());
    }
    for (int i = 0; i < (int) arguments.size(); i++) {
        kernel.setArg<cl::Memory>(index++, arguments[i].getMemory());
    }
    return kernel;
}
