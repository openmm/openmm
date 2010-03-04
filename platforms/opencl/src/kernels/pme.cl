__kernel void updateGridIndexAndFraction(__global float4* posq, __global float2* pmeAtomGridIndex) {
    for (int i = get_global_id(0); i < NUM_ATOMS; i += get_global_size(0)) {
        float4 pos = posq[i];
        float4 t = (float4) ((pos.x/PERIODIC_BOX_SIZE_X+1.0f)*GRID_SIZE_X,
                             (pos.y/PERIODIC_BOX_SIZE_Y+1.0f)*GRID_SIZE_Y,
                             (pos.z/PERIODIC_BOX_SIZE_Z+1.0f)*GRID_SIZE_Z, 0.0f);
        int4 gridIndex = (int4) (((int) t.x) % GRID_SIZE_X,
                                 ((int) t.y) % GRID_SIZE_Y,
                                 ((int) t.z) % GRID_SIZE_Z, 0);
        pmeAtomGridIndex[i] = (float2) (i, gridIndex.x*GRID_SIZE_Y*GRID_SIZE_Z+gridIndex.y*GRID_SIZE_Z+gridIndex.z);
    }
}

/**
 * For each grid point, find the range of sorted atoms associated with that point.
 */

__kernel void findAtomRangeForGrid(__global float4* posq, __global float2* pmeAtomGridIndex, __global int* pmeAtomRange) {
    int start = (NUM_ATOMS*get_global_id(0))/get_global_size(0);
    int end = (NUM_ATOMS*(get_global_id(0)+1))/get_global_size(0);
    int last = (start == 0 ? -1 : pmeAtomGridIndex[start-1].y);
    for (int i = start; i < end; ++i) {
        float2 atomData = pmeAtomGridIndex[i];
        int gridIndex = atomData.y;
        if (gridIndex != last) {
            for (int j = last+1; j <= gridIndex; ++j)
                pmeAtomRange[j] = i;
            last = gridIndex;
        }

        // The grid index won't be needed again.  Reuse that component to hold the atom charge, thus saving
        // an extra load operation in the charge spreading kernel.

        pmeAtomGridIndex[i].y = posq[(int) atomData.x].w;
    }

    // Fill in values beyond the last atom.

    if (get_global_id(0) == get_global_size(0)-1) {
        int gridSize = GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z;
        for (int j = last+1; j <= gridSize; ++j)
            pmeAtomRange[j] = NUM_ATOMS;
    }
}

__kernel void updateBsplines(__global float4* posq, __global float4* pmeBsplineTheta, __global float4* pmeBsplineDTheta, __local float4* bsplinesCache) {
    const float4 scale = 1.0f/(PME_ORDER-1);
    for (int i = get_global_id(0); i < NUM_ATOMS; i += get_global_size(0)) {
        __local float4* data = &bsplinesCache[get_local_id(0)*PME_ORDER];
        __local float4* ddata = &bsplinesCache[get_local_id(0)*PME_ORDER + get_local_size(0)*PME_ORDER];
        for (int j = 0; j < PME_ORDER; j++) {
	    data[j] = 0.0f;
            ddata[j] = 0.0f;
        }
        float4 pos = posq[i];
        float4 t = (float4) ((pos.x/PERIODIC_BOX_SIZE_X+1.0f)*GRID_SIZE_X,
                             (pos.y/PERIODIC_BOX_SIZE_Y+1.0f)*GRID_SIZE_Y,
                             (pos.z/PERIODIC_BOX_SIZE_Z+1.0f)*GRID_SIZE_Z, 0.0f);
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
        for (int j = 0; j < PME_ORDER; j++) {
            pmeBsplineTheta[i+j*NUM_ATOMS] = data[j];
            pmeBsplineDTheta[i+j*NUM_ATOMS] = ddata[j];
        }
    }
}

__kernel void gridSpreadCharge(__global float2* pmeAtomGridIndex, __global int* pmeAtomRange, __global float2* pmeGrid, __global float4* pmeBsplineTheta) {
    unsigned int numGridPoints = GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z;
    for (int gridIndex = get_global_id(0); gridIndex < numGridPoints; gridIndex += get_global_size(0)) {
        int4 gridPoint;
        gridPoint.x = gridIndex/(GRID_SIZE_Y*GRID_SIZE_Z);
        int remainder = gridIndex-gridPoint.x*GRID_SIZE_Y*GRID_SIZE_Z;
        gridPoint.y = remainder/GRID_SIZE_Z;
        gridPoint.z = remainder-gridPoint.y*GRID_SIZE_Z;
        gridPoint.x += GRID_SIZE_X;
        gridPoint.y += GRID_SIZE_Y;
        gridPoint.z += GRID_SIZE_Z;
        float result = 0.0f;
        for (int ix = 0; ix < PME_ORDER; ++ix)
            for (int iy = 0; iy < PME_ORDER; ++iy)
                for (int iz = 0; iz < PME_ORDER; ++iz) {
                    int x = (gridPoint.x-ix)%GRID_SIZE_X;
                    int y = (gridPoint.y-iy)%GRID_SIZE_Y;
                    int z = (gridPoint.z-iz)%GRID_SIZE_Z;
                    int gridIndex = x*GRID_SIZE_Y*GRID_SIZE_Z+y*GRID_SIZE_Z+z;
                    int firstAtom = pmeAtomRange[gridIndex];
                    int lastAtom = pmeAtomRange[gridIndex+1];
                    for (int i = firstAtom; i < lastAtom; ++i) {
                        float2 atomData = pmeAtomGridIndex[i];
                        int atomIndex = atomData.x;
                        float atomCharge = atomData.y;
                        result += atomCharge*pmeBsplineTheta[atomIndex+ix*NUM_ATOMS].x*pmeBsplineTheta[atomIndex+iy*NUM_ATOMS].y*pmeBsplineTheta[atomIndex+iz*NUM_ATOMS].z;
                    }
                }
        pmeGrid[gridIndex] = (float2) (result*EPSILON_FACTOR, 0.0f);
    }
}

__kernel void reciprocalConvolution(__global float2* pmeGrid, __global float* energyBuffer, __global float* pmeBsplineModuliX,
        __global float* pmeBsplineModuliY, __global float* pmeBsplineModuliZ) {
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
        float mhx = mx/PERIODIC_BOX_SIZE_X;
        float mhy = my/PERIODIC_BOX_SIZE_Y;
        float mhz = mz/PERIODIC_BOX_SIZE_Z;
        float bx = pmeBsplineModuliX[kx];
        float by = pmeBsplineModuliY[ky];
        float bz = pmeBsplineModuliZ[kz];
        float2 grid = pmeGrid[index];
        float m2 = mhx*mhx+mhy*mhy+mhz*mhz;
        float denom = m2*bx*by*bz;
        float eterm = RECIP_SCALE_FACTOR*exp(-RECIP_EXP_FACTOR*m2)/denom;
        pmeGrid[index] = (float2) (grid.x*eterm, grid.y*eterm);
        energy += eterm*(grid.x*grid.x + grid.y*grid.y);
    }
    energyBuffer[get_global_id(0)] += 0.5f*energy;
}

__kernel void gridInterpolateForce(__global float4* posq, __global float4* forceBuffers, __global float4* pmeBsplineTheta, __global float4* pmeBsplineDTheta, __global float2* pmeGrid) {
    for (int atom = get_global_id(0); atom < NUM_ATOMS; atom += get_global_size(0)) {
        float4 force = 0.0f;
        float4 pos = posq[atom];
        float4 t = (float4) ((pos.x/PERIODIC_BOX_SIZE_X+1.0f)*GRID_SIZE_X,
                             (pos.y/PERIODIC_BOX_SIZE_Y+1.0f)*GRID_SIZE_Y,
                             (pos.z/PERIODIC_BOX_SIZE_Z+1.0f)*GRID_SIZE_Z, 0.0f);
        int4 gridIndex = (int4) (((int) t.x) % GRID_SIZE_X,
                                 ((int) t.y) % GRID_SIZE_Y,
                                 ((int) t.z) % GRID_SIZE_Z, 0);
        for (int ix = 0; ix < PME_ORDER; ix++) {
            int xindex = (gridIndex.x + ix) % GRID_SIZE_X;
            float tx = pmeBsplineTheta[atom+ix*NUM_ATOMS].x;
            float dtx = pmeBsplineDTheta[atom+ix*NUM_ATOMS].x;
            for (int iy = 0; iy < PME_ORDER; iy++) {
                int yindex = (gridIndex.y + iy) % GRID_SIZE_Y;
                float ty = pmeBsplineTheta[atom+iy*NUM_ATOMS].y;
                float dty = pmeBsplineDTheta[atom+iy*NUM_ATOMS].y;
                for (int iz = 0; iz < PME_ORDER; iz++) {
                    int zindex = (gridIndex.z + iz) % GRID_SIZE_Z;
                    float tz = pmeBsplineTheta[atom+iz*NUM_ATOMS].z;
                    float dtz = pmeBsplineDTheta[atom+iz*NUM_ATOMS].z;
                    int index = xindex*GRID_SIZE_Y*GRID_SIZE_Z + yindex*GRID_SIZE_Z + zindex;
                    float gridvalue = pmeGrid[index].x;
                    force.x += dtx*ty*tz*gridvalue;
                    force.y += tx*dty*tz*gridvalue;
                    force.z += tx*ty*dtz*gridvalue;
                }
            }
        }
        float4 totalForce = forceBuffers[atom];
        float q = pos.w*EPSILON_FACTOR;
        totalForce.x -= q*force.x*GRID_SIZE_X/PERIODIC_BOX_SIZE_X;
        totalForce.y -= q*force.y*GRID_SIZE_Y/PERIODIC_BOX_SIZE_Y;
        totalForce.z -= q*force.z*GRID_SIZE_Z/PERIODIC_BOX_SIZE_Z;
        forceBuffers[atom] = totalForce;
    }
}
