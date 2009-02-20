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
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/**
 * This file contains the kernels for identifying interacting blocks.  It is included
 * several times in kCalculateCDLJForces.cu with different #defines to generate
 * different versions of the kernels.
 */

/**
 * Find a bounding box for the atoms in each block.
 */
__global__ void METHOD_NAME(kFindBlockBounds, _kernel)()
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int base = pos << GRIDBITS;
    if (base < cSim.atoms)
    {
        float4 apos = cSim.pPosq[base];
#ifdef USE_PERIODIC
        apos.x -= floor(apos.x/cSim.periodicBoxSizeX)*cSim.periodicBoxSizeX;
        apos.y -= floor(apos.y/cSim.periodicBoxSizeY)*cSim.periodicBoxSizeY;
        apos.z -= floor(apos.z/cSim.periodicBoxSizeZ)*cSim.periodicBoxSizeZ;
        float4 firstPoint = apos;
#endif
        float minx = apos.x;
        float maxx = apos.x;
        float miny = apos.y;
        float maxy = apos.y;
        float minz = apos.z;
        float maxz = apos.z;
        for (unsigned int i = 1; i < GRID; i++)
        {
            apos = cSim.pPosq[base+i];
#ifdef USE_PERIODIC
            apos.x -= floor((apos.x-firstPoint.x)/cSim.periodicBoxSizeX+0.5)*cSim.periodicBoxSizeX;
            apos.y -= floor((apos.y-firstPoint.y)/cSim.periodicBoxSizeY+0.5)*cSim.periodicBoxSizeY;
            apos.z -= floor((apos.z-firstPoint.z)/cSim.periodicBoxSizeZ+0.5)*cSim.periodicBoxSizeZ;
#endif
            minx = min(minx, apos.x);
            maxx = max(maxx, apos.x);
            miny = min(miny, apos.y);
            maxy = max(maxy, apos.y);
            minz = min(minz, apos.z);
            maxz = max(maxz, apos.z);
        }
        cSim.pGridBoundingBox[pos] = make_float4(0.5f*(maxx-minx), 0.5f*(maxy-miny), 0.5f*(maxz-minz), 0);
        cSim.pGridCenter[pos] = make_float4(0.5f*(maxx+minx), 0.5f*(maxy+miny), 0.5f*(maxz+minz), 0);
    }
}

/**
 * Compare the bounding boxes for each pair of blocks.  If they are sufficiently far apart,
 * mark them as non-interacting.
 */
__global__ void METHOD_NAME(kFindBlocksWithInteractions, _kernel)()
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
    if (pos < cSim.workUnits)
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
        dx -= floor(dx/cSim.periodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
        dy -= floor(dy/cSim.periodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
        dz -= floor(dz/cSim.periodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;
#endif
        float4 boxSizea = cSim.pGridBoundingBox[x];
        float4 boxSizeb = cSim.pGridBoundingBox[y];
        dx = max(0.0f, abs(dx)-boxSizea.x-boxSizeb.x);
        dy = max(0.0f, abs(dy)-boxSizea.y-boxSizeb.y);
        dz = max(0.0f, abs(dz)-boxSizea.z-boxSizeb.z);
        cSim.pInteractionFlag[pos] = (dx*dx+dy*dy+dz*dz > cSim.nonbondedCutoffSqr ? 0 : 1);
    }
}
