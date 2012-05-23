__kernel void updateBsplines(__global const float4* restrict posq, __global float4* restrict pmeBsplineTheta, __local float4* restrict bsplinesCache,
        __global int2* restrict pmeAtomGridIndex, float4 periodicBoxSize, float4 invPeriodicBoxSize) {
    const float4 scale = 1.0f/(PME_ORDER-1);
    for (int i = get_global_id(0); i < NUM_ATOMS; i += get_global_size(0)) {
        __local float4* data = &bsplinesCache[get_local_id(0)*PME_ORDER];
        float4 pos = posq[i];
        pos.x -= floor(pos.x*invPeriodicBoxSize.x)*periodicBoxSize.x;
        pos.y -= floor(pos.y*invPeriodicBoxSize.y)*periodicBoxSize.y;
        pos.z -= floor(pos.z*invPeriodicBoxSize.z)*periodicBoxSize.z;
        float4 t = (float4) ((pos.x*invPeriodicBoxSize.x)*GRID_SIZE_X,
                             (pos.y*invPeriodicBoxSize.y)*GRID_SIZE_Y,
                             (pos.z*invPeriodicBoxSize.z)*GRID_SIZE_Z, 0.0f);
        float4 dr = (float4) (t.x-(int) t.x, t.y-(int) t.y, t.z-(int) t.z, 0.0f);
        int4 gridIndex = (int4) (((int) t.x) % GRID_SIZE_X,
                                 ((int) t.y) % GRID_SIZE_Y,
                                 ((int) t.z) % GRID_SIZE_Z, 0);
        pmeAtomGridIndex[i] = (int2) (i, gridIndex.x*GRID_SIZE_Y*GRID_SIZE_Z+gridIndex.y*GRID_SIZE_Z+gridIndex.z);
        data[PME_ORDER-1] = 0.0f;
        data[1] = dr;
        data[0] = 1.0f-dr;
        for (int j = 3; j < PME_ORDER; j++) {
            float div = 1.0f/(j-1.0f);
            data[j-1] = div*dr*data[j-2];
            for (int k = 1; k < (j-1); k++)
                data[j-k-1] = div*((dr+(float4) k) *data[j-k-2] + (-dr+(float4) (j-k))*data[j-k-1]);
            data[0] = div*(- dr+1.0f)*data[0];
        }
        data[PME_ORDER-1] = scale*dr*data[PME_ORDER-2];
        for (int j = 1; j < (PME_ORDER-1); j++)
            data[PME_ORDER-j-1] = scale*((dr+(float4) j)*data[PME_ORDER-j-2] + (-dr+(float4) (PME_ORDER-j))*data[PME_ORDER-j-1]);
        data[0] = scale*(-dr+1.0f)*data[0];
        for (int j = 0; j < PME_ORDER; j++) {
            data[j].w = pos.w; // Storing the charge here improves cache coherency in the charge spreading kernel
            pmeBsplineTheta[i+j*NUM_ATOMS] = data[j];
        }
    }
}

/**
 * For each grid point, find the range of sorted atoms associated with that point.
 */
__kernel void findAtomRangeForGrid(__global int2* restrict pmeAtomGridIndex, __global int* restrict pmeAtomRange, __global const float4* restrict posq, float4 periodicBoxSize, float4 invPeriodicBoxSize) {
    int start = (NUM_ATOMS*get_global_id(0))/get_global_size(0);
    int end = (NUM_ATOMS*(get_global_id(0)+1))/get_global_size(0);
    int last = (start == 0 ? -1 : pmeAtomGridIndex[start-1].y);
    for (int i = start; i < end; ++i) {
        int2 atomData = pmeAtomGridIndex[i];
        int gridIndex = atomData.y;
        if (gridIndex != last) {
            for (int j = last+1; j <= gridIndex; ++j)
                pmeAtomRange[j] = i;
            last = gridIndex;
        }
    }

    // Fill in values beyond the last atom.

    if (get_global_id(0) == get_global_size(0)-1) {
        int gridSize = GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z;
        for (int j = last+1; j <= gridSize; ++j)
            pmeAtomRange[j] = NUM_ATOMS;
    }
}

/**
 * The grid index won't be needed again.  Reuse that component to hold the z index, thus saving
 * some work in the charge spreading kernel.
 */
__kernel void recordZIndex(__global int2* restrict pmeAtomGridIndex, __global const float4* restrict posq, float4 periodicBoxSize, float4 invPeriodicBoxSize) {
    int start = (NUM_ATOMS*get_global_id(0))/get_global_size(0);
    int end = (NUM_ATOMS*(get_global_id(0)+1))/get_global_size(0);
    for (int i = start; i < end; ++i) {
        float posz = posq[pmeAtomGridIndex[i].x].z;
        posz -= floor(posz*invPeriodicBoxSize.z)*periodicBoxSize.z;
        int z = ((int) ((posz*invPeriodicBoxSize.z)*GRID_SIZE_Z)) % GRID_SIZE_Z;
        pmeAtomGridIndex[i].y = z;
    }
}

#ifdef SUPPORTS_64_BIT_ATOMICS
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable

#define BUFFER_SIZE (PME_ORDER*PME_ORDER*PME_ORDER)

__kernel __attribute__((reqd_work_group_size(BUFFER_SIZE, 1, 1)))
void gridSpreadCharge(__global const float4* restrict posq, __global const int2* restrict pmeAtomGridIndex, __global const int* restrict pmeAtomRange,
        __global long* restrict pmeGrid, __global const float4* restrict pmeBsplineTheta, float4 periodicBoxSize, float4 invPeriodicBoxSize) {
    int ix = get_local_id(0)/(PME_ORDER*PME_ORDER);
    int remainder = get_local_id(0)-ix*PME_ORDER*PME_ORDER;
    int iy = remainder/PME_ORDER;
    int iz = remainder-iy*PME_ORDER;
    __local float4 theta[PME_ORDER];
    __local float charge[BUFFER_SIZE];
    __local int basex[BUFFER_SIZE];
    __local int basey[BUFFER_SIZE];
    __local int basez[BUFFER_SIZE];
    if (ix < PME_ORDER) {
        for (int baseIndex = get_group_id(0)*BUFFER_SIZE; baseIndex < NUM_ATOMS; baseIndex += get_num_groups(0)*BUFFER_SIZE) {
            // Load the next block of atoms into the buffers.

            if (get_local_id(0) < BUFFER_SIZE) {
                int atomIndex = baseIndex+get_local_id(0);
                if (atomIndex < NUM_ATOMS) {
                    float4 pos = posq[atomIndex];
                    charge[get_local_id(0)] = pos.w;
                    pos.x -= floor(pos.x*invPeriodicBoxSize.x)*periodicBoxSize.x;
                    pos.y -= floor(pos.y*invPeriodicBoxSize.y)*periodicBoxSize.y;
                    pos.z -= floor(pos.z*invPeriodicBoxSize.z)*periodicBoxSize.z;
                    basex[get_local_id(0)] = (int) ((pos.x*invPeriodicBoxSize.x)*GRID_SIZE_X);
                    basey[get_local_id(0)] = (int) ((pos.y*invPeriodicBoxSize.y)*GRID_SIZE_Y);
                    basez[get_local_id(0)] = (int) ((pos.z*invPeriodicBoxSize.z)*GRID_SIZE_Z);
                }
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            int lastIndex = min(BUFFER_SIZE, NUM_ATOMS-baseIndex);
            for (int index = 0; index < lastIndex; index++) {
                int atomIndex = index+baseIndex;
                if (get_local_id(0) < PME_ORDER)
                    theta[get_local_id(0)] = pmeBsplineTheta[atomIndex+get_local_id(0)*NUM_ATOMS];
                barrier(CLK_LOCAL_MEM_FENCE);
                float add = charge[index]*theta[ix].x*theta[iy].y*theta[iz].z;
                int x = basex[index]+ix;
                int y = basey[index]+iy;
                int z = basez[index]+iz;
                x -= (x >= GRID_SIZE_X ? GRID_SIZE_X : 0);
                y -= (y >= GRID_SIZE_Y ? GRID_SIZE_Y : 0);
                z -= (z >= GRID_SIZE_Z ? GRID_SIZE_Z : 0);
                atom_add(&pmeGrid[x*GRID_SIZE_Y*GRID_SIZE_Z+y*GRID_SIZE_Z+z], (long) (add*0xFFFFFFFF));
            }
        }
    }
}

__kernel void finishSpreadCharge(__global long* restrict pmeGrid) {
    __global float2* floatGrid = (__global float2*) pmeGrid;
    const unsigned int gridSize = GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z;
    float scale = EPSILON_FACTOR/(float) 0xFFFFFFFF;
    for (int index = get_global_id(0); index < gridSize; index += get_global_size(0)) {
        long value = pmeGrid[index];
        float2 floatValue = (float2) ((float) (value*scale), 0.0f);
        floatGrid[index] = floatValue;
    }
}
#else
__kernel void gridSpreadCharge(__global const float4* restrict posq, __global const int2* restrict pmeAtomGridIndex, __global const int* restrict pmeAtomRange,
        __global float2* restrict pmeGrid, __global const float4* restrict pmeBsplineTheta) {
    unsigned int numGridPoints = GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z;
    for (int gridIndex = get_global_id(0); gridIndex < numGridPoints; gridIndex += get_global_size(0)) {
        // Compute the charge on a grid point.

        int4 gridPoint;
        gridPoint.x = gridIndex/(GRID_SIZE_Y*GRID_SIZE_Z);
        int remainder = gridIndex-gridPoint.x*GRID_SIZE_Y*GRID_SIZE_Z;
        gridPoint.y = remainder/GRID_SIZE_Z;
        gridPoint.z = remainder-gridPoint.y*GRID_SIZE_Z;
        float result = 0.0f;

        // Loop over all atoms that affect this grid point.

        for (int ix = 0; ix < PME_ORDER; ++ix) {
            int x = gridPoint.x-ix+(gridPoint.x >= ix ? 0 : GRID_SIZE_X);
            for (int iy = 0; iy < PME_ORDER; ++iy) {
                int y = gridPoint.y-iy+(gridPoint.y >= iy ? 0 : GRID_SIZE_Y);
                int z1 = gridPoint.z-PME_ORDER+1;
                z1 += (z1 >= 0 ? 0 : GRID_SIZE_Z);
                int z2 = (z1 < gridPoint.z ? gridPoint.z : GRID_SIZE_Z-1);
                int gridIndex1 = x*GRID_SIZE_Y*GRID_SIZE_Z+y*GRID_SIZE_Z+z1;
                int gridIndex2 = x*GRID_SIZE_Y*GRID_SIZE_Z+y*GRID_SIZE_Z+z2;
                int firstAtom = pmeAtomRange[gridIndex1];
                int lastAtom = pmeAtomRange[gridIndex2+1];
                for (int i = firstAtom; i < lastAtom; ++i)
                {
                    int2 atomData = pmeAtomGridIndex[i];
                    int atomIndex = atomData.x;
                    int z = atomData.y;
                    int iz = gridPoint.z-z+(gridPoint.z >= z ? 0 : GRID_SIZE_Z);
                    float atomCharge = pmeBsplineTheta[atomIndex+ix*NUM_ATOMS].w;
                    result += atomCharge*pmeBsplineTheta[atomIndex+ix*NUM_ATOMS].x*pmeBsplineTheta[atomIndex+iy*NUM_ATOMS].y*pmeBsplineTheta[atomIndex+iz*NUM_ATOMS].z;
                }
                if (z1 > gridPoint.z)
                {
                    gridIndex1 = x*GRID_SIZE_Y*GRID_SIZE_Z+y*GRID_SIZE_Z;
                    gridIndex2 = x*GRID_SIZE_Y*GRID_SIZE_Z+y*GRID_SIZE_Z+gridPoint.z;
                    firstAtom = pmeAtomRange[gridIndex1];
                    lastAtom = pmeAtomRange[gridIndex2+1];
                    for (int i = firstAtom; i < lastAtom; ++i)
                    {
                        int2 atomData = pmeAtomGridIndex[i];
                        int atomIndex = atomData.x;
                        int z = atomData.y;
                        int iz = gridPoint.z-z+(gridPoint.z >= z ? 0 : GRID_SIZE_Z);
                        float atomCharge = pmeBsplineTheta[atomIndex+ix*NUM_ATOMS].w;
                        result += atomCharge*pmeBsplineTheta[atomIndex+ix*NUM_ATOMS].x*pmeBsplineTheta[atomIndex+iy*NUM_ATOMS].y*pmeBsplineTheta[atomIndex+iz*NUM_ATOMS].z;
                    }
                }
            }
        }
        pmeGrid[gridIndex] = (float2) (result*EPSILON_FACTOR, 0.0f);
    }
}
#endif

__kernel void reciprocalConvolution(__global float2* restrict pmeGrid, __global float* restrict energyBuffer, __global const float* restrict pmeBsplineModuliX,
        __global const float* restrict pmeBsplineModuliY, __global const float* restrict pmeBsplineModuliZ, float4 invPeriodicBoxSize, float recipScaleFactor) {
    const unsigned int gridSize = GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z;
    float energy = 0.0f;
    for (int index = get_global_id(0); index < gridSize; index += get_global_size(0)) {
        int kx = index/(GRID_SIZE_Y*GRID_SIZE_Z);
        int remainder = index-kx*GRID_SIZE_Y*GRID_SIZE_Z;
        int ky = remainder/GRID_SIZE_Z;
        int kz = remainder-ky*GRID_SIZE_Z;
        if (kx == 0 && ky == 0 && kz == 0)
            continue;
        int mx = (kx < (GRID_SIZE_X+1)/2) ? kx : (kx-GRID_SIZE_X);
        int my = (ky < (GRID_SIZE_Y+1)/2) ? ky : (ky-GRID_SIZE_Y);
        int mz = (kz < (GRID_SIZE_Z+1)/2) ? kz : (kz-GRID_SIZE_Z);
        float mhx = mx*invPeriodicBoxSize.x;
        float mhy = my*invPeriodicBoxSize.y;
        float mhz = mz*invPeriodicBoxSize.z;
        float bx = pmeBsplineModuliX[kx];
        float by = pmeBsplineModuliY[ky];
        float bz = pmeBsplineModuliZ[kz];
        float2 grid = pmeGrid[index];
        float m2 = mhx*mhx+mhy*mhy+mhz*mhz;
        float denom = m2*bx*by*bz;
        float eterm = recipScaleFactor*EXP(-RECIP_EXP_FACTOR*m2)/denom;
        pmeGrid[index] = (float2) (grid.x*eterm, grid.y*eterm);
        energy += eterm*(grid.x*grid.x + grid.y*grid.y);
    }
    energyBuffer[get_global_id(0)] += 0.5f*energy;
}

__kernel void gridInterpolateForce(__global const float4* restrict posq, __global float4* restrict forceBuffers, __global const float2* restrict pmeGrid,
        float4 periodicBoxSize, float4 invPeriodicBoxSize, __local float4* restrict bsplinesCache) {
    const float4 scale = 1.0f/(PME_ORDER-1);
    __local float4* data = &bsplinesCache[get_local_id(0)*PME_ORDER];
    __local float4* ddata = &bsplinesCache[get_local_id(0)*PME_ORDER + get_local_size(0)*PME_ORDER];
    for (int atom = get_global_id(0); atom < NUM_ATOMS; atom += get_global_size(0)) {
        float4 force = 0.0f;
        float4 pos = posq[atom];
        pos.x -= floor(pos.x*invPeriodicBoxSize.x)*periodicBoxSize.x;
        pos.y -= floor(pos.y*invPeriodicBoxSize.y)*periodicBoxSize.y;
        pos.z -= floor(pos.z*invPeriodicBoxSize.z)*periodicBoxSize.z;
        float4 t = (float4) ((pos.x*invPeriodicBoxSize.x)*GRID_SIZE_X,
                             (pos.y*invPeriodicBoxSize.y)*GRID_SIZE_Y,
                             (pos.z*invPeriodicBoxSize.z)*GRID_SIZE_Z, 0.0f);
        int4 gridIndex = (int4) (((int) t.x) % GRID_SIZE_X,
                                 ((int) t.y) % GRID_SIZE_Y,
                                 ((int) t.z) % GRID_SIZE_Z, 0);

        // Since we need the full set of thetas, it's faster to compute them here than load them
        // from global memory.

        float4 dr = (float4) (t.x-(int) t.x, t.y-(int) t.y, t.z-(int) t.z, 0.0f);
        data[PME_ORDER-1] = 0.0f;
        data[1] = dr;
        data[0] = 1.0f-dr;
        for (int j = 3; j < PME_ORDER; j++) {
            float div = 1.0f/(j-1.0f);
            data[j-1] = div*dr*data[j-2];
            for (int k = 1; k < (j-1); k++)
                data[j-k-1] = div*((dr+(float4) k) *data[j-k-2] + (-dr+(float4) (j-k))*data[j-k-1]);
            data[0] = div*(- dr+1.0f)*data[0];
        }
        ddata[0] = -data[0];
        for (int j = 1; j < PME_ORDER; j++)
            ddata[j] = data[j-1]-data[j];
        data[PME_ORDER-1] = scale*dr*data[PME_ORDER-2];
        for (int j = 1; j < (PME_ORDER-1); j++)
            data[PME_ORDER-j-1] = scale*((dr+(float4) j)*data[PME_ORDER-j-2] + (-dr+(float4) (PME_ORDER-j))*data[PME_ORDER-j-1]);
        data[0] = scale*(-dr+1.0f)*data[0];

        // Compute the force on this atom.

        for (int ix = 0; ix < PME_ORDER; ix++) {
            int xindex = gridIndex.x+ix;
            xindex -= (xindex >= GRID_SIZE_X ? GRID_SIZE_X : 0);
            for (int iy = 0; iy < PME_ORDER; iy++) {
                int yindex = gridIndex.y+iy;
                yindex -= (yindex >= GRID_SIZE_Y ? GRID_SIZE_Y : 0);
                for (int iz = 0; iz < PME_ORDER; iz++) {
                    int zindex = gridIndex.z+iz;
                    zindex -= (zindex >= GRID_SIZE_Z ? GRID_SIZE_Z : 0);
                    int index = xindex*GRID_SIZE_Y*GRID_SIZE_Z + yindex*GRID_SIZE_Z + zindex;
                    float gridvalue = pmeGrid[index].x;
                    force.x += ddata[ix].x*data[iy].y*data[iz].z*gridvalue;
                    force.y += data[ix].x*ddata[iy].y*data[iz].z*gridvalue;
#ifndef MAC_AMD_WORKAROUND
                    force.z += data[ix].x*data[iy].y*ddata[iz].z*gridvalue;
#endif
                }
            }
        }
#ifdef MAC_AMD_WORKAROUND
        for (int ix = 0; ix < PME_ORDER; ix++) {
            int xindex = gridIndex.x+ix;
            xindex -= (xindex >= GRID_SIZE_X ? GRID_SIZE_X : 0);
            for (int iy = 0; iy < PME_ORDER; iy++) {
                int yindex = gridIndex.y+iy;
                yindex -= (yindex >= GRID_SIZE_Y ? GRID_SIZE_Y : 0);
                for (int iz = 0; iz < PME_ORDER; iz++) {
                    int zindex = gridIndex.z+iz;
                    zindex -= (zindex >= GRID_SIZE_Z ? GRID_SIZE_Z : 0);
                    int index = xindex*GRID_SIZE_Y*GRID_SIZE_Z + yindex*GRID_SIZE_Z + zindex;
                    float gridvalue = pmeGrid[index].x;
                    force.z += data[ix].x*data[iy].y*ddata[iz].z*gridvalue;
                }
            }
        }
#endif
        float4 totalForce = forceBuffers[atom];
        float q = pos.w*EPSILON_FACTOR;
        totalForce.x -= q*force.x*GRID_SIZE_X*invPeriodicBoxSize.x;
        totalForce.y -= q*force.y*GRID_SIZE_Y*invPeriodicBoxSize.y;
        totalForce.z -= q*force.z*GRID_SIZE_Z*invPeriodicBoxSize.z;
        forceBuffers[atom] = totalForce;
    }
}
