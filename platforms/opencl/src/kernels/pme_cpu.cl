__kernel void updateBsplines(__global const real4* restrict posq, __global real4* restrict pmeBsplineTheta, __local real4* restrict bsplinesCache, __global int2* restrict pmeAtomGridIndex, real4 periodicBoxSize, real4 invPeriodicBoxSize, __global real4* restrict pmeBsplineDTheta) {
    const real4 scale = 1.0f/(PME_ORDER-1);
    for (int i = get_global_id(0); i < NUM_ATOMS; i += get_global_size(0)) {
        __local real4* data = &bsplinesCache[get_local_id(0)*PME_ORDER];
        __local real4* ddata = &bsplinesCache[get_local_id(0)*PME_ORDER + get_local_size(0)*PME_ORDER];
        for (int j = 0; j < PME_ORDER; j++) {
	    data[j] = 0;
            ddata[j] = 0;
        }
        real4 pos = posq[i];
        pos.x -= floor(pos.x*invPeriodicBoxSize.x)*periodicBoxSize.x;
        pos.y -= floor(pos.y*invPeriodicBoxSize.y)*periodicBoxSize.y;
        pos.z -= floor(pos.z*invPeriodicBoxSize.z)*periodicBoxSize.z;
        real4 t = (real4) ((pos.x*invPeriodicBoxSize.x)*GRID_SIZE_X,
                             (pos.y*invPeriodicBoxSize.y)*GRID_SIZE_Y,
                             (pos.z*invPeriodicBoxSize.z)*GRID_SIZE_Z, 0);
        real4 dr = (real4) (t.x-(int) t.x, t.y-(int) t.y, t.z-(int) t.z, 0);
        data[PME_ORDER-1] = 0;
        data[1] = dr;
        data[0] = 1.0f-dr;
        for (int j = 3; j < PME_ORDER; j++) {
            real div = 1.0f/(j-1.0f);
            data[j-1] = div*dr*data[j-2];
            for (int k = 1; k < (j-1); k++)
                data[j-k-1] = div*((dr+(real4) k) *data[j-k-2] + (-dr+(real4) (j-k))*data[j-k-1]);
            data[0] = div*(- dr+1.0f)*data[0];
        }
        ddata[0] = -data[0];
        for (int j = 1; j < PME_ORDER; j++)
            ddata[j] = data[j-1]-data[j];
        data[PME_ORDER-1] = scale*dr*data[PME_ORDER-2];
        for (int j = 1; j < (PME_ORDER-1); j++)
            data[PME_ORDER-j-1] = scale*((dr+(real4) j)*data[PME_ORDER-j-2] + (-dr+(real4) (PME_ORDER-j))*data[PME_ORDER-j-1]);
        data[0] = scale*(-dr+1.0f)*data[0];
        for (int j = 0; j < PME_ORDER; j++) {
            pmeBsplineTheta[i+j*NUM_ATOMS] = data[j];
            pmeBsplineDTheta[i+j*NUM_ATOMS] = ddata[j];
        }
    }
}

/**
 * This kernel is not actually used when running on a CPU.
 */
__kernel void findAtomRangeForGrid(__global const int2* restrict pmeAtomGridIndex, __global int* restrict pmeAtomRange, __global const real4* restrict posq, real4 periodicBoxSize, real4 invPeriodicBoxSize) {
}

__kernel void gridSpreadCharge(__global const real4* restrict posq, __global const int2* restrict pmeAtomGridIndex, __global const int* restrict pmeAtomRange, __global real2* restrict pmeGrid, __global const real4* restrict pmeBsplineTheta, real4 periodicBoxSize, real4 invPeriodicBoxSize) {
    const int firstx = get_global_id(0)*GRID_SIZE_X/get_global_size(0);
    const int lastx = (get_global_id(0)+1)*GRID_SIZE_X/get_global_size(0);
    for (int gridIndex = firstx*GRID_SIZE_Y*GRID_SIZE_Z; gridIndex < lastx*GRID_SIZE_Y*GRID_SIZE_Z; gridIndex++)
        pmeGrid[gridIndex] = (real2) 0;
    for (int atom = 0; atom < NUM_ATOMS; atom++) {
        real4 pos = posq[atom];
        pos.x -= floor(pos.x*invPeriodicBoxSize.x)*periodicBoxSize.x;
        pos.y -= floor(pos.y*invPeriodicBoxSize.y)*periodicBoxSize.y;
        pos.z -= floor(pos.z*invPeriodicBoxSize.z)*periodicBoxSize.z;
        real4 t = (real4) ((pos.x*invPeriodicBoxSize.x)*GRID_SIZE_X,
                             (pos.y*invPeriodicBoxSize.y)*GRID_SIZE_Y,
                             (pos.z*invPeriodicBoxSize.z)*GRID_SIZE_Z, 0);
        real4 dr = (real4) (t.x-(int) t.x, t.y-(int) t.y, t.z-(int) t.z, 0);
        int4 gridIndex = (int4) (((int) t.x) % GRID_SIZE_X,
                                 ((int) t.y) % GRID_SIZE_Y,
                                 ((int) t.z) % GRID_SIZE_Z, 0);
        real atomCharge = pos.w*EPSILON_FACTOR;
        for (int ix = 0; ix < PME_ORDER; ix++) {
            int xindex = gridIndex.x+ix;
            xindex -= (xindex >= GRID_SIZE_X ? GRID_SIZE_X : 0);
            if (xindex < firstx || xindex >= lastx)
                continue;
            for (int iy = 0; iy < PME_ORDER; iy++) {
                int yindex = gridIndex.y+iy;
                yindex -= (yindex >= GRID_SIZE_Y ? GRID_SIZE_Y : 0);
                for(int iz = 0; iz < PME_ORDER; iz++) {
                    int zindex = gridIndex.z+iz;
                    zindex -= (zindex >= GRID_SIZE_Z ? GRID_SIZE_Z : 0);
                    int index = xindex*GRID_SIZE_Y*GRID_SIZE_Z + yindex*GRID_SIZE_Z + zindex;
                    pmeGrid[index].x += atomCharge*pmeBsplineTheta[atom+ix*NUM_ATOMS].x*pmeBsplineTheta[atom+iy*NUM_ATOMS].y*pmeBsplineTheta[atom+iz*NUM_ATOMS].z;
                }
            }
        }
    }
}

__kernel void reciprocalConvolution(__global real2* restrict pmeGrid, __global real* restrict energyBuffer, __global const real* restrict pmeBsplineModuliX,
        __global const real* restrict pmeBsplineModuliY, __global const real* restrict pmeBsplineModuliZ, real4 invPeriodicBoxSize, real recipScaleFactor) {
    const unsigned int gridSize = GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z;
    real energy = 0;
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
        real mhx = mx*invPeriodicBoxSize.x;
        real mhy = my*invPeriodicBoxSize.y;
        real mhz = mz*invPeriodicBoxSize.z;
        real bx = pmeBsplineModuliX[kx];
        real by = pmeBsplineModuliY[ky];
        real bz = pmeBsplineModuliZ[kz];
        real2 grid = pmeGrid[index];
        real m2 = mhx*mhx+mhy*mhy+mhz*mhz;
        real denom = m2*bx*by*bz;
        real eterm = recipScaleFactor*EXP(-RECIP_EXP_FACTOR*m2)/denom;
        pmeGrid[index] = (real2) (grid.x*eterm, grid.y*eterm);
        energy += eterm*(grid.x*grid.x + grid.y*grid.y);
    }
    energyBuffer[get_global_id(0)] += 0.5f*energy;
}

__kernel void gridInterpolateForce(__global const real4* restrict posq, __global real4* restrict forceBuffers, __global const real2* restrict pmeGrid, real4 periodicBoxSize, real4 invPeriodicBoxSize, __global const real4* restrict pmeBsplineTheta, __global const real4* restrict pmeBsplineDTheta) {
    for (int atom = get_global_id(0); atom < NUM_ATOMS; atom += get_global_size(0)) {
        real4 force = 0;
        real4 pos = posq[atom];
        pos.x -= floor(pos.x*invPeriodicBoxSize.x)*periodicBoxSize.x;
        pos.y -= floor(pos.y*invPeriodicBoxSize.y)*periodicBoxSize.y;
        pos.z -= floor(pos.z*invPeriodicBoxSize.z)*periodicBoxSize.z;
        real4 t = (real4) ((pos.x*invPeriodicBoxSize.x)*GRID_SIZE_X,
                           (pos.y*invPeriodicBoxSize.y)*GRID_SIZE_Y,
                           (pos.z*invPeriodicBoxSize.z)*GRID_SIZE_Z, 0);
        int4 gridIndex = (int4) (((int) t.x) % GRID_SIZE_X,
                                 ((int) t.y) % GRID_SIZE_Y,
                                 ((int) t.z) % GRID_SIZE_Z, 0);
        for (int ix = 0; ix < PME_ORDER; ix++) {
            int xindex = gridIndex.x+ix;
            xindex -= (xindex >= GRID_SIZE_X ? GRID_SIZE_X : 0);
            real tx = pmeBsplineTheta[atom+ix*NUM_ATOMS].x;
            real dtx = pmeBsplineDTheta[atom+ix*NUM_ATOMS].x;
            for (int iy = 0; iy < PME_ORDER; iy++) {
                int yindex = gridIndex.y+iy;
                yindex -= (yindex >= GRID_SIZE_Y ? GRID_SIZE_Y : 0);
                real ty = pmeBsplineTheta[atom+iy*NUM_ATOMS].y;
                real dty = pmeBsplineDTheta[atom+iy*NUM_ATOMS].y;
                for (int iz = 0; iz < PME_ORDER; iz++) {
                    int zindex = gridIndex.z+iz;
                    zindex -= (zindex >= GRID_SIZE_Z ? GRID_SIZE_Z : 0);
                    real tz = pmeBsplineTheta[atom+iz*NUM_ATOMS].z;
                    real dtz = pmeBsplineDTheta[atom+iz*NUM_ATOMS].z;
                    int index = xindex*GRID_SIZE_Y*GRID_SIZE_Z + yindex*GRID_SIZE_Z + zindex;
                    real gridvalue = pmeGrid[index].x;
                    force.x += dtx*ty*tz*gridvalue;
                    force.y += tx*dty*tz*gridvalue;
                    force.z += tx*ty*dtz*gridvalue;
                }
            }
        }
        real4 totalForce = forceBuffers[atom];
        real q = pos.w*EPSILON_FACTOR;
        totalForce.x -= q*force.x*GRID_SIZE_X*invPeriodicBoxSize.x;
        totalForce.y -= q*force.y*GRID_SIZE_Y*invPeriodicBoxSize.y;
        totalForce.z -= q*force.z*GRID_SIZE_Z*invPeriodicBoxSize.z;
        forceBuffers[atom] = totalForce;
    }
}
