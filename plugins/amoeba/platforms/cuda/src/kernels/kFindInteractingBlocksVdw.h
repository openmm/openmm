/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Scott Le Grand, Peter Eastman                                     *
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

/**
 * Compare the bounding boxes for each pair of blocks.  If they are sufficiently far apart,
 * mark them as non-interacting.
 */
__global__ void METHOD_NAME(kFindBlocksWithInteractionsVdw, _kernel)()
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.workUnits)
    {
        // Extract cell coordinates from appropriate work unit

        unsigned int x = cSim.pWorkUnit[pos];
        unsigned int y = ((x >> 2) & 0x7fff);
        x = (x >> 17);

        // Find the distance between the bounding boxes of the two cells.

        float4 centera = cSim.pGridCenter[x];
        float4 centerb = cSim.pGridCenter[y];
        float dx = centera.x-centerb.x;
        float dy = centera.y-centerb.y;
        float dz = centera.z-centerb.z;
#ifdef USE_PERIODIC
        dx -= floorf(dx*cSim.invPeriodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
        dy -= floorf(dy*cSim.invPeriodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
        dz -= floorf(dz*cSim.invPeriodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;
#endif
        float4 boxSizea = cSim.pGridBoundingBox[x];
        float4 boxSizeb = cSim.pGridBoundingBox[y];
        dx = max(0.0f, abs(dx)-boxSizea.x-boxSizeb.x);
        dy = max(0.0f, abs(dy)-boxSizea.y-boxSizeb.y);
        dz = max(0.0f, abs(dz)-boxSizea.z-boxSizeb.z);
        cSim.pInteractionFlag[pos] = (dx*dx+dy*dy+dz*dz > cAmoebaSim.vdwCutoff2 ? 0 : 1);
        pos += gridDim.x*blockDim.x;
    }
}

/**
 * Compare each atom in one block to the bounding box of another block, and set
 * flags for which ones are interacting.
 */
__global__ void METHOD_NAME(kFindInteractionsWithinBlocksVdw, _kernel)(unsigned int* workUnit)
{
    extern __shared__ volatile unsigned int flags[];
    unsigned int totalWarps = cSim.nonbond_blocks*cSim.nonbond_threads_per_block/GRID;
    unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    unsigned int numWorkUnits = cSim.pInteractionCount[0];
    unsigned int pos = warp*numWorkUnits/totalWarps;
    unsigned int end = (warp+1)*numWorkUnits/totalWarps;
    unsigned int index = threadIdx.x & (GRID - 1);

    unsigned int lasty = 0xFFFFFFFF;
    float4 apos;
    while (pos < end)
    {
        // Extract cell coordinates from appropriate work unit
        unsigned int x = workUnit[pos];
        unsigned int y = ((x >> 2) & 0x7fff);
        bool bExclusionFlag = (x & 0x1);
        x = (x >> 17);
        if (x == y || bExclusionFlag)
        {
            // Assume this block will be dense.

            if (index == 0)
                cSim.pInteractionFlag[pos] = 0xFFFFFFFF;
        }
        else
        {
            // Load the bounding box for x and the atom positions for y.

            float4 center = cSim.pGridCenter[x];
            float4 boxSize = cSim.pGridBoundingBox[x];
            if (y != lasty)
            {
                apos = cSim.pPosq[(y<<GRIDBITS)+index];
            }

            // Find the distance of the atom from the bounding box.

            float dx = apos.x-center.x;
            float dy = apos.y-center.y;
            float dz = apos.z-center.z;
#ifdef USE_PERIODIC
            dx -= floorf(dx*cSim.invPeriodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
            dy -= floorf(dy*cSim.invPeriodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
            dz -= floorf(dz*cSim.invPeriodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;
#endif
            dx = max(0.0f, abs(dx)-boxSize.x);
            dy = max(0.0f, abs(dy)-boxSize.y);
            dz = max(0.0f, abs(dz)-boxSize.z);
            flags[threadIdx.x] = (dx*dx+dy*dy+dz*dz > cAmoebaSim.vdwCutoff2 ? 0 : 1 << index);

            // Sum the flags.

            if (index % 2 == 0)
                flags[threadIdx.x] += flags[threadIdx.x+1];
            if (index % 4 == 0)
                flags[threadIdx.x] += flags[threadIdx.x+2];
            if (index % 8 == 0)
                flags[threadIdx.x] += flags[threadIdx.x+4];
            if (index % 16 == 0)
                flags[threadIdx.x] += flags[threadIdx.x+8];
            if (index == 0)
            {
                unsigned int allFlags = flags[threadIdx.x] + flags[threadIdx.x+16];

                // Count how many flags are set, and based on that decide whether to compute all interactions
                // or only a fraction of them.

                unsigned int bits = (allFlags&0x55555555) + ((allFlags>>1)&0x55555555);
                bits = (bits&0x33333333) + ((bits>>2)&0x33333333);
                bits = (bits&0x0F0F0F0F) + ((bits>>4)&0x0F0F0F0F);
                bits = (bits&0x00FF00FF) + ((bits>>8)&0x00FF00FF);
                bits = (bits&0x0000FFFF) + ((bits>>16)&0x0000FFFF);
                cSim.pInteractionFlag[pos] = (bits > 12 ? 0xFFFFFFFF : allFlags);
            }
            lasty = y;
        }
        pos++;
    }
}
