/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
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
#include "OpenCLCompact.h"
#include <map>

using namespace OpenMM;
using namespace std;

OpenCLNonbondedUtilities::OpenCLNonbondedUtilities(OpenCLContext& context) : context(context), cutoff(-1.0), useCutoff(false),
        numForceBuffers(0), tiles(NULL), exclusionIndex(NULL), exclusions(NULL), interactingTiles(NULL), interactionFlags(NULL),
        interactionCount(NULL), blockCenter(NULL), blockBoundingBox(NULL), compact(NULL) {
    // Decide how many force buffers to use.

    forceBufferPerAtomBlock = false;
    numForceBuffers = context.getNumThreadBlocks()*OpenCLContext::ThreadBlockSize/OpenCLContext::TileSize;
    if (numForceBuffers >= context.getNumAtomBlocks()) {
        // For small systems, it is more efficient to have one force buffer per block of 32 atoms instead of one per warp.

        forceBufferPerAtomBlock = true;
        numForceBuffers = context.getNumAtomBlocks();
    }
}

OpenCLNonbondedUtilities::~OpenCLNonbondedUtilities() {
    if (tiles != NULL)
        delete tiles;
    if (exclusionIndex != NULL)
        delete exclusionIndex;
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
    if (compact != NULL)
        delete compact;
}

void OpenCLNonbondedUtilities::addInteraction(bool usesCutoff, bool usesPeriodic, double cutoffDistance, const vector<vector<int> >& exclusionList, const string& kernel) {
    if (cutoff != -1.0) {
        if (usesCutoff != useCutoff)
            throw OpenMMException("All Forces must agree on whether to use a cutoff");
        if (usesPeriodic != usePeriodic)
            throw OpenMMException("All Forces must agree on whether to use periodic boundary conditions");
        if (cutoffDistance != cutoff)
            throw OpenMMException("All Forces must use the same cutoff distance");
        bool sameExclusions = (exclusionList.size() == atomExclusions.size());
        for (int i = 0; i < exclusionList.size() && sameExclusions; i++) {
            if (exclusionList[i].size() != atomExclusions[i].size())
                sameExclusions = false;
            for (int j = 0; j < exclusionList[i].size(); j++)
                if (exclusionList[i][j] != atomExclusions[i][j])
                    sameExclusions = false;
        }
        if (!sameExclusions)
            throw OpenMMException("All Forces must have identical exceptions");
    }
    else {
        useCutoff = usesCutoff;
        usePeriodic = usesPeriodic;
        cutoff = cutoffDistance;
        atomExclusions = exclusionList;
        kernelSource += kernel+"\n";
    }
}

void OpenCLNonbondedUtilities::addParameter(const ParameterInfo& parameter) {
    parameters.push_back(parameter);
}

void OpenCLNonbondedUtilities::initialize(const System& system) {
    if (cutoff == -1.0)
        return; // There are no nonbonded interactions in the System.
    
    // Create the list of tiles.

    int numAtomBlocks = context.getNumAtomBlocks();
    int numTiles = numAtomBlocks*(numAtomBlocks+1)/2;
    tiles = new OpenCLArray<cl_uint>(context, numTiles, "tiles");
    vector<cl_uint> tileVec(tiles->getSize());
    unsigned int count = 0;
    for (unsigned int y = 0; y < numAtomBlocks; y++)
        for (unsigned int x = y; x < numAtomBlocks; x++)
            tileVec[count++] = (x << 17) | (y << 2);

    // Create kernels.

    forceKernel = createInteractionKernel(kernelSource, parameters, true);
    map<string, string> defines;
    if (forceBufferPerAtomBlock)
        defines["USE_OUTPUT_BUFFER_PER_BLOCK"] = "1";
    if (useCutoff)
        defines["USE_CUTOFF"] = "1";
    if (usePeriodic)
        defines["USE_PERIODIC"] = "1";
    cl::Program interactingBlocksProgram = context.createProgram(context.loadSourceFromFile("findInteractingBlocks.cl"), defines);
    findBlockBoundsKernel = cl::Kernel(interactingBlocksProgram, "findBlockBounds");
    findInteractingBlocksKernel = cl::Kernel(interactingBlocksProgram, "findBlocksWithInteractions");
    findInteractionsWithinBlocksKernel = cl::Kernel(interactingBlocksProgram, "findInteractionsWithinBlocks");

    // Mark which tiles have exclusions.

    for (int atom1 = 0; atom1 < (int) atomExclusions.size(); ++atom1) {
        int x = atom1/OpenCLContext::TileSize;
        for (int j = 0; j < (int) atomExclusions[atom1].size(); ++j) {
            int atom2 = atomExclusions[atom1][j];
            int y = atom2/OpenCLContext::TileSize;
            int index = (x > y ? x+y*numAtomBlocks-y*(y+1)/2 : y+x*numAtomBlocks-x*(x+1)/2);
            tileVec[index] |= 1;
        }
    }
    if (context.getPaddedNumAtoms() > context.getNumAtoms()) {
        int lastTile = context.getNumAtoms()/OpenCLContext::TileSize;
        for (int i = 0; i < numTiles; ++i) {
            int x = tileVec[i]>>17;
            int y = (tileVec[i]>>2)&0x7FFF;
            if (x == lastTile || y == lastTile)
                tileVec[i] |= 1;
        }
    }

    // Build a list of indices for the tiles with exclusions.

    exclusionIndex = new OpenCLArray<cl_uint>(context, numTiles, "exclusionIndex");
    vector<cl_uint> exclusionIndexVec(exclusionIndex->getSize());
    int numWithExclusions = 0;
    for (int i = 0; i < numTiles; ++i)
        if ((tileVec[i]&1) == 1)
            exclusionIndexVec[i] = (numWithExclusions++)*OpenCLContext::TileSize;

    // Record the exclusion data.

    exclusions = new OpenCLArray<cl_uint>(context, numWithExclusions*OpenCLContext::TileSize, "exclusions");
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
                int tile = x+y*numAtomBlocks-y*(y+1)/2;
                exclusionVec[exclusionIndexVec[tile]+offset1] &= 0xFFFFFFFF-(1<<offset2);
            }
            else {
                int tile = y+x*numAtomBlocks-x*(x+1)/2;
                exclusionVec[exclusionIndexVec[tile]+offset2] &= 0xFFFFFFFF-(1<<offset1);
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
                int tile = x+y*numAtomBlocks-y*(y+1)/2;
                exclusionVec[exclusionIndexVec[tile]+offset1] &= 0xFFFFFFFF-(1<<offset2);
            }
            if (y >= x) {
                int tile = y+x*numAtomBlocks-x*(x+1)/2;
                exclusionVec[exclusionIndexVec[tile]+offset2] &= 0xFFFFFFFF-(1<<offset1);
            }
        }
    }
    atomExclusions.clear(); // We won't use this again, so free the memory it used
    tiles->upload(tileVec);
    exclusions->upload(exclusionVec);
    exclusionIndex->upload(exclusionIndexVec);

    // Record the periodic box size.

    Vec3 boxVectors[3];
    system.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    periodicBoxSize = (mm_float4) {(float) boxVectors[0][0], (float) boxVectors[1][1], (float) boxVectors[2][2], 0.0f};

    // Create data structures for the neighbor list.

    if (useCutoff) {
        interactingTiles = new OpenCLArray<cl_uint>(context, numTiles, "interactingTiles");
        interactionFlags = new OpenCLArray<cl_uint>(context, numTiles, "interactionFlags");
        interactionCount = new OpenCLArray<cl_uint>(context, 1, "interactionCount");
        blockCenter = new OpenCLArray<mm_float4>(context, numAtomBlocks, "blockCenter");
        blockBoundingBox = new OpenCLArray<mm_float4>(context, numAtomBlocks, "blockBoundingBox");
        compact = new OpenCLCompact(context);
    }
}

void OpenCLNonbondedUtilities::prepareInteractions() {
    if (!useCutoff)
        return;

    // Compute the neighbor list.

    findBlockBoundsKernel.setArg<cl_int>(0, context.getNumAtoms());
    findBlockBoundsKernel.setArg<mm_float4>(1, periodicBoxSize);
    findBlockBoundsKernel.setArg<cl::Buffer>(2, context.getPosq().getDeviceBuffer());
    findBlockBoundsKernel.setArg<cl::Buffer>(3, blockCenter->getDeviceBuffer());
    findBlockBoundsKernel.setArg<cl::Buffer>(4, blockBoundingBox->getDeviceBuffer());
    context.executeKernel(findBlockBoundsKernel, context.getNumAtoms());

    findInteractingBlocksKernel.setArg<cl_int>(0, tiles->getSize());
    findInteractingBlocksKernel.setArg<cl_float>(1, cutoff*cutoff);
    findInteractingBlocksKernel.setArg<mm_float4>(2, periodicBoxSize);
    findInteractingBlocksKernel.setArg<cl::Buffer>(3, tiles->getDeviceBuffer());
    findInteractingBlocksKernel.setArg<cl::Buffer>(4, blockCenter->getDeviceBuffer());
    findInteractingBlocksKernel.setArg<cl::Buffer>(5, blockBoundingBox->getDeviceBuffer());
    findInteractingBlocksKernel.setArg<cl::Buffer>(6, interactionFlags->getDeviceBuffer());
    context.executeKernel(findInteractingBlocksKernel, context.getNumAtoms());

    compact->compactStream(*interactingTiles, *tiles, *interactionFlags, *interactionCount);

    findInteractionsWithinBlocksKernel.setArg<cl_float>(0, cutoff*cutoff);
    findInteractionsWithinBlocksKernel.setArg<mm_float4>(1, periodicBoxSize);
    findInteractionsWithinBlocksKernel.setArg<cl::Buffer>(2, context.getPosq().getDeviceBuffer());
    findInteractionsWithinBlocksKernel.setArg<cl::Buffer>(3, interactingTiles->getDeviceBuffer());
    findInteractionsWithinBlocksKernel.setArg<cl::Buffer>(4, blockCenter->getDeviceBuffer());
    findInteractionsWithinBlocksKernel.setArg<cl::Buffer>(5, blockBoundingBox->getDeviceBuffer());
    findInteractionsWithinBlocksKernel.setArg<cl::Buffer>(6, interactionFlags->getDeviceBuffer());
    findInteractionsWithinBlocksKernel.setArg<cl::Buffer>(7, interactionCount->getDeviceBuffer());
    findInteractionsWithinBlocksKernel.setArg(8, OpenCLContext::ThreadBlockSize*sizeof(cl_uint), NULL);
    context.executeKernel(findInteractionsWithinBlocksKernel, context.getNumAtoms());
}

void OpenCLNonbondedUtilities::computeInteractions() {
    invokeInteractionKernel(forceKernel, parameters);
}

cl::Kernel OpenCLNonbondedUtilities::createInteractionKernel(const string& source, const vector<ParameterInfo> params, bool useExclusions) const {
    map<string, string> replacements;
    replacements["COMPUTE_INTERACTION"] = source;
    stringstream args;
    for (int i = 0; i < params.size(); i++) {
        args << ", __global ";
        args << params[i].getType();
        args << "* global_";
        args << params[i].getName();
        args << ", __local ";
        args << params[i].getType();
        args << "* local_";
        args << params[i].getName();
    }
    replacements["PARAMETER_ARGUMENTS"] = args.str();
//            local_sigmaEpsilon[get_local_id(0)] = sigmaEpsilon1;
    stringstream loadLocal1;
    for (int i = 0; i < params.size(); i++) {
        loadLocal1 << "local_";
        loadLocal1 << params[i].getName();
        loadLocal1 << "[get_local_id(0)] = ";
        loadLocal1 << params[i].getName();
        loadLocal1 << "1;";
    }
    replacements["LOAD_LOCAL_PARAMETERS_FROM_1"] = loadLocal1.str();
//                local_sigmaEpsilon[get_local_id(0)] = global_sigmaEpsilon[j];
    stringstream loadLocal2;
    for (int i = 0; i < params.size(); i++) {
        loadLocal2 << "local_";
        loadLocal2 << params[i].getName();
        loadLocal2 << "[get_local_id(0)] = global_";
        loadLocal2 << params[i].getName();
        loadLocal2 << "[j];";
    }
    replacements["LOAD_LOCAL_PARAMETERS_FROM_GLOBAL"] = loadLocal2.str();
    stringstream load1;
    for (int i = 0; i < params.size(); i++) {
        load1 << params[i].getType();
        load1 << " ";
        load1 << params[i].getName();
        load1 << "1 = global_";
        load1 << params[i].getName();
        load1 << "[i];";
    }
    replacements["LOAD_ATOM1_PARAMETERS"] = load1.str();
    stringstream load2j;
    for (int i = 0; i < params.size(); i++) {
        load2j << params[i].getType();
        load2j << " ";
        load2j << params[i].getName();
        load2j << "2 = local_";
        load2j << params[i].getName();
        load2j << "[tbx+j];";
    }
    replacements["LOAD_ATOM2_PARAMETERS_J"] = load2j.str();
    stringstream load2tj;
    for (int i = 0; i < params.size(); i++) {
        load2tj << params[i].getType();
        load2tj << " ";
        load2tj << params[i].getName();
        load2tj << "2 = local_";
        load2tj << params[i].getName();
        load2tj << "[tbx+tj];";
    }
    replacements["LOAD_ATOM2_PARAMETERS_TJ"] = load2tj.str();
    map<string, string> defines;
    if (forceBufferPerAtomBlock)
        defines["USE_OUTPUT_BUFFER_PER_BLOCK"] = "1";
    if (useCutoff)
        defines["USE_CUTOFF"] = "1";
    if (usePeriodic)
        defines["USE_PERIODIC"] = "1";
    if (useExclusions)
        defines["USE_EXCLUSIONS"] = "1";
    cl::Program program = context.createProgram(context.loadSourceFromFile("nonbonded.cl", replacements), defines);
    return cl::Kernel(program, "computeNonbonded");
}

void OpenCLNonbondedUtilities::invokeInteractionKernel(cl::Kernel kernel, const vector<ParameterInfo>& params) {
    kernel.setArg<cl_int>(0, context.getPaddedNumAtoms());
    kernel.setArg<cl::Buffer>(1, context.getForceBuffers().getDeviceBuffer());
    kernel.setArg<cl::Buffer>(2, context.getEnergyBuffer().getDeviceBuffer());
    kernel.setArg<cl::Buffer>(3, context.getPosq().getDeviceBuffer());
    kernel.setArg<cl::Buffer>(4, exclusions->getDeviceBuffer());
    kernel.setArg<cl::Buffer>(5, exclusionIndex->getDeviceBuffer());
    kernel.setArg(6, OpenCLContext::ThreadBlockSize*sizeof(cl_float4), NULL);
    kernel.setArg(7, OpenCLContext::ThreadBlockSize*sizeof(cl_float4), NULL);
    int paramBase = 10;
    if (useCutoff) {
        paramBase = 14;
        kernel.setArg<cl::Buffer>(8, interactingTiles->getDeviceBuffer());
        kernel.setArg<cl_float>(9, cutoff*cutoff);
        kernel.setArg<mm_float4>(10, periodicBoxSize);
        kernel.setArg<cl::Buffer>(11, interactionFlags->getDeviceBuffer());
        kernel.setArg<cl::Buffer>(12, interactionCount->getDeviceBuffer());
        kernel.setArg(13, OpenCLContext::ThreadBlockSize*sizeof(cl_float4), NULL);
    }
    else {
        kernel.setArg<cl::Buffer>(8, tiles->getDeviceBuffer());
        kernel.setArg<cl_uint>(9, tiles->getSize());
    }
    for (int i = 0; i < (int) params.size(); i++) {
        kernel.setArg<cl::Buffer>(i*2+paramBase, params[i].getBuffer());
        kernel.setArg(i*2+paramBase+1, OpenCLContext::ThreadBlockSize*params[i].getSize(), NULL);
    }
    context.executeKernel(kernel, tiles->getSize()*OpenCLContext::TileSize);
}
