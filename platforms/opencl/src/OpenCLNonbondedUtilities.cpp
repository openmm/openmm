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
#include <map>

using namespace OpenMM;
using namespace std;

OpenCLNonbondedUtilities::OpenCLNonbondedUtilities(OpenCLContext& context) : context(context), cutoff(-1.0), useCutoff(false),
        numForceBuffers(0), tiles(NULL), exclusionIndex(NULL), exclusions(NULL) {
}

OpenCLNonbondedUtilities::~OpenCLNonbondedUtilities() {
    if (tiles != NULL)
        delete tiles;
    if (exclusionIndex != NULL)
        delete exclusionIndex;
    if (exclusions != NULL)
        delete exclusions;
}

void OpenCLNonbondedUtilities::addInteraction(bool usesCutoff, bool usesPeriodic, double cutoffDistance, const std::vector<std::vector<int> >& exclusionList) {
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
    }
}

void OpenCLNonbondedUtilities::addParameter(const string& name, const string& type, int size, cl::Buffer& buffer) {
    parameters.push_back(ParameterInfo(name, type, size, buffer));
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

    // Decide how many force buffers to use.

    bool forceBufferPerAtomBlock = false;
    numForceBuffers = context.getNumThreadBlocks()*OpenCLContext::ThreadBlockSize/OpenCLContext::TileSize;
    if (numForceBuffers >= numAtomBlocks) {
        // For small systems, it is more efficient to have one force buffer per block of 32 atoms instead of one per warp.

        forceBufferPerAtomBlock = true;
        numForceBuffers = numAtomBlocks;
    }

    // Create kernels.

    cl::Program forceProgram = context.createProgram(context.loadSourceFromFile("nonbonded.cl"));
    forceKernel = cl::Kernel(forceProgram, "computeNonbonded");

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
    Vec3 boxVectors[3];
    system.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    periodicBoxSize = (mm_float4) {(float) boxVectors[0][0], (float) boxVectors[1][1], (float) boxVectors[2][2], 0.0f};
}

void OpenCLNonbondedUtilities::prepareInteractions() {
    hasComputedInteractions = false;
    if (!useCutoff)
        return;
    // TODO compute the neighbor list
}

void OpenCLNonbondedUtilities::computeInteractions() {
    if (hasComputedInteractions)
        return;
    hasComputedInteractions = true;
    forceKernel.setArg<cl_int>(0, tiles->getSize());
    forceKernel.setArg<cl_int>(1, context.getPaddedNumAtoms());
    forceKernel.setArg<cl_float>(2, cutoff*cutoff);
    forceKernel.setArg<mm_float4>(3, periodicBoxSize);
    forceKernel.setArg<cl::Buffer>(4, context.getForceBuffers().getDeviceBuffer());
    forceKernel.setArg<cl::Buffer>(5, context.getEnergyBuffer().getDeviceBuffer());
    forceKernel.setArg<cl::Buffer>(6, context.getPosq().getDeviceBuffer());
    forceKernel.setArg<cl::Buffer>(7, tiles->getDeviceBuffer());
    forceKernel.setArg<cl::Buffer>(8, exclusions->getDeviceBuffer());
    forceKernel.setArg<cl::Buffer>(9, exclusionIndex->getDeviceBuffer());
    forceKernel.setArg(10, OpenCLContext::ThreadBlockSize*sizeof(cl_float4), NULL);
    forceKernel.setArg(11, OpenCLContext::ThreadBlockSize*sizeof(cl_float4), NULL);
    for (int i = 0; i < (int) parameters.size(); i++) {
        forceKernel.setArg<cl::Buffer>(i*2+12, *parameters[i].buffer);
        forceKernel.setArg(i*2+13, OpenCLContext::ThreadBlockSize*parameters[i].size, NULL);
    }
    context.executeKernel(forceKernel, tiles->getSize()*OpenCLContext::TileSize);
}
