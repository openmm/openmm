__kernel void updateBsplines(__global const real4* restrict posq, __global real4* restrict pmeBsplineTheta, __local real4* restrict bsplinesCache,
        __global int2* restrict pmeAtomGridIndex, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        real4 recipBoxVecX, real4 recipBoxVecY, real4 recipBoxVecZ
#ifdef USE_LJPME
        , __global const float2* restrict sigmaEpsilon
#else
        , __global const real* restrict charges
#endif
    ) {
    const real4 scale = 1/(real) (PME_ORDER-1);
    for (int atom = get_global_id(0); atom < NUM_ATOMS; atom += get_global_size(0)) {
        __local real4* data = &bsplinesCache[get_local_id(0)*PME_ORDER];
        real4 pos = posq[atom];
        APPLY_PERIODIC_TO_POS(pos)
        real3 t = (real3) (pos.x*recipBoxVecX.x+pos.y*recipBoxVecY.x+pos.z*recipBoxVecZ.x,
                           pos.y*recipBoxVecY.y+pos.z*recipBoxVecZ.y,
                           pos.z*recipBoxVecZ.z);
        t.x = (t.x-floor(t.x))*GRID_SIZE_X;
        t.y = (t.y-floor(t.y))*GRID_SIZE_Y;
        t.z = (t.z-floor(t.z))*GRID_SIZE_Z;
        real4 dr = (real4) (t.x-(int) t.x, t.y-(int) t.y, t.z-(int) t.z, 0.0f);
        int4 gridIndex = (int4) (((int) t.x) % GRID_SIZE_X,
                                 ((int) t.y) % GRID_SIZE_Y,
                                 ((int) t.z) % GRID_SIZE_Z, 0);
        pmeAtomGridIndex[atom] = (int2) (atom, gridIndex.x*GRID_SIZE_Y*GRID_SIZE_Z+gridIndex.y*GRID_SIZE_Z+gridIndex.z);
#ifndef SUPPORTS_64_BIT_ATOMICS
        data[PME_ORDER-1] = 0.0f;
        data[1] = dr;
        data[0] = 1.0f-dr;
        for (int j = 3; j < PME_ORDER; j++) {
            real div = RECIP(j-1.0f);
            data[j-1] = div*dr*data[j-2];
            for (int k = 1; k < (j-1); k++)
                data[j-k-1] = div*((dr+(real4) k) *data[j-k-2] + (-dr+(real4) (j-k))*data[j-k-1]);
            data[0] = div*(- dr+1.0f)*data[0];
        }
        data[PME_ORDER-1] = scale*dr*data[PME_ORDER-2];
        for (int j = 1; j < (PME_ORDER-1); j++)
            data[PME_ORDER-j-1] = scale*((dr+(real4) j)*data[PME_ORDER-j-2] + (-dr+(real4) (PME_ORDER-j))*data[PME_ORDER-j-1]);
        data[0] = scale*(-dr+1.0f)*data[0];
        for (int j = 0; j < PME_ORDER; j++) {
#ifdef USE_LJPME
            const float2 sigEps = sigmaEpsilon[atom];
            const real charge = 8*sigEps.x*sigEps.x*sigEps.x*sigEps.y;
#else
            const real charge = CHARGE;
#endif
            data[j].w = charge; // Storing the charge here improves cache coherency in the charge spreading kernel
            pmeBsplineTheta[atom+j*NUM_ATOMS] = data[j];
        }
#endif
    }
}

/**
 * For each grid point, find the range of sorted atoms associated with that point.
 */
__kernel void findAtomRangeForGrid(__global int2* restrict pmeAtomGridIndex, __global int* restrict pmeAtomRange, __global const real4* restrict posq) {
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
__kernel void recordZIndex(__global int2* restrict pmeAtomGridIndex, __global const real4* restrict posq, real4 periodicBoxSize, real4 recipBoxVecZ) {
    int start = (NUM_ATOMS*get_global_id(0))/get_global_size(0);
    int end = (NUM_ATOMS*(get_global_id(0)+1))/get_global_size(0);
    for (int i = start; i < end; ++i) {
        real posz = posq[pmeAtomGridIndex[i].x].z;
        posz -= floor(posz*recipBoxVecZ.z)*periodicBoxSize.z;
        int z = ((int) ((posz*recipBoxVecZ.z)*GRID_SIZE_Z)) % GRID_SIZE_Z;
        pmeAtomGridIndex[i].y = z;
    }
}

#ifdef SUPPORTS_64_BIT_ATOMICS
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable

__kernel void gridSpreadCharge(__global const real4* restrict posq, __global const int2* restrict pmeAtomGridIndex, __global const int* restrict pmeAtomRange,
        __global long* restrict pmeGrid, __global const real4* restrict pmeBsplineTheta, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, real4 recipBoxVecX, real4 recipBoxVecY, real4 recipBoxVecZ
#ifdef USE_LJPME
        , __global const float2* restrict sigmaEpsilon
#else
        , __global const real* restrict charges
#endif
    ) {
    // To improve memory efficiency, we divide indices along the z axis into
    // PME_ORDER blocks, where the data for each block is stored together.  We
    // can ensure that all threads write to the same block at the same time,
    // which leads to better coalescing of writes.
    
    __local int zindexTable[GRID_SIZE_Z+PME_ORDER];
    int blockSize = (int) ceil(GRID_SIZE_Z/(real) PME_ORDER);
    for (int i = get_local_id(0); i < GRID_SIZE_Z+PME_ORDER; i += get_local_size(0)) {
        int zindex = i % GRID_SIZE_Z;
	int block = zindex % PME_ORDER;
        zindexTable[i] = zindex/PME_ORDER + block*GRID_SIZE_X*GRID_SIZE_Y*blockSize;
    }
    barrier(CLK_LOCAL_MEM_FENCE);
        
    // Process the atoms in spatially sorted order.  This improves efficiency when writing
    // the grid values.
    
    const real scale = 1/(real) (PME_ORDER-1);
    real4 data[PME_ORDER];
    for (int i = get_global_id(0); i < NUM_ATOMS; i += get_global_size(0)) {
        int atom = pmeAtomGridIndex[i].x;
        real4 pos = posq[atom];
#ifdef USE_LJPME
        const float2 sigEps = sigmaEpsilon[atom];
        const real charge = 8*sigEps.x*sigEps.x*sigEps.x*sigEps.y;
#else
        const real charge = (CHARGE)*EPSILON_FACTOR;
#endif
        if (charge == 0)
            continue;
        APPLY_PERIODIC_TO_POS(pos)
        real3 t = (real3) (pos.x*recipBoxVecX.x+pos.y*recipBoxVecY.x+pos.z*recipBoxVecZ.x,
                           pos.y*recipBoxVecY.y+pos.z*recipBoxVecZ.y,
                           pos.z*recipBoxVecZ.z);
        t.x = (t.x-floor(t.x))*GRID_SIZE_X;
        t.y = (t.y-floor(t.y))*GRID_SIZE_Y;
        t.z = (t.z-floor(t.z))*GRID_SIZE_Z;
        int4 gridIndex = (int4) (((int) t.x) % GRID_SIZE_X,
                                 ((int) t.y) % GRID_SIZE_Y,
                                 ((int) t.z) % GRID_SIZE_Z, 0);

        // Since we need the full set of thetas, it's faster to compute them here than load them
        // from global memory.

        real4 dr = (real4) (t.x-(int) t.x, t.y-(int) t.y, t.z-(int) t.z, 0.0f);
        data[PME_ORDER-1] = 0.0f;
        data[1] = dr;
        data[0] = 1.0f-dr;
        for (int j = 3; j < PME_ORDER; j++) {
            real div = RECIP(j-1.0f);
            data[j-1] = div*dr*data[j-2];
            for (int k = 1; k < (j-1); k++)
                data[j-k-1] = div*((dr+(real4) k) *data[j-k-2] + (-dr+(real4) (j-k))*data[j-k-1]);
            data[0] = div*(-dr+1.0f)*data[0];
        }
        data[PME_ORDER-1] = scale*dr*data[PME_ORDER-2];
        for (int j = 1; j < (PME_ORDER-1); j++)
            data[PME_ORDER-j-1] = scale*((dr+(real4) j)*data[PME_ORDER-j-2] + (-dr+(real4) (PME_ORDER-j))*data[PME_ORDER-j-1]);
        data[0] = scale*(-dr+1.0f)*data[0];

        // Spread the charge from this atom onto each grid point.

	int izoffset = (PME_ORDER-(gridIndex.z%PME_ORDER)) % PME_ORDER;
        for (int ix = 0; ix < PME_ORDER; ix++) {
            int xbase = gridIndex.x+ix;
            xbase -= (xbase >= GRID_SIZE_X ? GRID_SIZE_X : 0);
            xbase = xbase*GRID_SIZE_Y;
            real dx = charge*data[ix].x;
            for (int iy = 0; iy < PME_ORDER; iy++) {
                int ybase = gridIndex.y+iy;
                ybase -= (ybase >= GRID_SIZE_Y ? GRID_SIZE_Y : 0);
                ybase = (xbase+ybase)*blockSize;
                real dxdy = dx*data[iy].y;
                for (int i = 0; i < PME_ORDER; i++) {
		    int iz = (i+izoffset) % PME_ORDER;
                    int zindex = gridIndex.z+iz;
                    int index = ybase + zindexTable[zindex];
                    real add = dxdy*data[iz].z;
                    atom_add(&pmeGrid[index], (long) (add*0x100000000));
                }
            }
        }
    }
}

__kernel void finishSpreadCharge(__global long* restrict grid1, __global real* restrict grid2) {
    // During charge spreading, we shuffled the order of indices along the z
    // axis to make memory access more efficient.  We now need to unshuffle
    // them and convert fixed point values to floating point.

    __local int zindexTable[GRID_SIZE_Z];
    int blockSize = (int) ceil(GRID_SIZE_Z/(real) PME_ORDER);
    for (int i = get_local_id(0); i < GRID_SIZE_Z; i += get_local_size(0)) {
	int block = i % PME_ORDER;
        zindexTable[i] = i/PME_ORDER + block*GRID_SIZE_X*GRID_SIZE_Y*blockSize;
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    const unsigned int gridSize = GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z;
    real scale = 1/(real) 0x100000000;
    for (int index = get_global_id(0); index < gridSize; index += get_global_size(0)) {
        int zindex = index%GRID_SIZE_Z;
        int loadIndex = zindexTable[zindex] + blockSize*(int) (index/GRID_SIZE_Z);
        grid2[index] = scale*grid1[loadIndex];
    }
}
#elif defined(DEVICE_IS_CPU)
__kernel void gridSpreadCharge(__global const real4* restrict posq, __global const int2* restrict pmeAtomGridIndex, __global const int* restrict pmeAtomRange,
        __global real* restrict pmeGrid, __global const real4* restrict pmeBsplineTheta, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, real4 recipBoxVecX, real4 recipBoxVecY, real4 recipBoxVecZ
#ifdef USE_LJPME
        , __global const float2* restrict sigmaEpsilon
#else
        , __global const real* restrict charges
#endif
    ) {
    const int firstx = get_global_id(0)*GRID_SIZE_X/get_global_size(0);
    const int lastx = (get_global_id(0)+1)*GRID_SIZE_X/get_global_size(0);
    if (firstx == lastx)
        return;
    const real4 scale = 1/(real) (PME_ORDER-1);
    real4 data[PME_ORDER];
    
    // Process the atoms in spatially sorted order.  This improves efficiency when writing
    // the grid values.
    
    for (int i = 0; i < NUM_ATOMS; i++) {
        int atom = i;//pmeAtomGridIndex[i].x;
        real4 pos = posq[atom];
        APPLY_PERIODIC_TO_POS(pos)
        real3 t = (real3) (pos.x*recipBoxVecX.x+pos.y*recipBoxVecY.x+pos.z*recipBoxVecZ.x,
                           pos.y*recipBoxVecY.y+pos.z*recipBoxVecZ.y,
                           pos.z*recipBoxVecZ.z);
        t.x = (t.x-floor(t.x))*GRID_SIZE_X;
        t.y = (t.y-floor(t.y))*GRID_SIZE_Y;
        t.z = (t.z-floor(t.z))*GRID_SIZE_Z;
        int4 gridIndex = (int4) (((int) t.x) % GRID_SIZE_X,
                                 ((int) t.y) % GRID_SIZE_Y,
                                 ((int) t.z) % GRID_SIZE_Z, 0);

        // Spread the charge from this atom onto each grid point.

#ifdef USE_LJPME
        const float2 sigEps = sigmaEpsilon[atom];
        const real charge = 8*sigEps.x*sigEps.x*sigEps.x*sigEps.y;
#else
        const real charge = (CHARGE)*EPSILON_FACTOR;
#endif
        if (charge == 0)
            continue;
        bool hasComputedThetas = false;
        for (int ix = 0; ix < PME_ORDER; ix++) {
            int xindex = gridIndex.x+ix;
            xindex -= (xindex >= GRID_SIZE_X ? GRID_SIZE_X : 0);
            if (xindex < firstx || xindex >= lastx)
                continue;
            if (!hasComputedThetas) {
                hasComputedThetas = true;
                
                // Since we need the full set of thetas, it's faster to compute them here than load them
                // from global memory.

                real4 dr = (real4) (t.x-(int) t.x, t.y-(int) t.y, t.z-(int) t.z, 0.0f);
                data[PME_ORDER-1] = 0.0f;
                data[1] = dr;
                data[0] = 1.0f-dr;
                for (int j = 3; j < PME_ORDER; j++) {
                    real div = RECIP(j-1.0f);
                    data[j-1] = div*dr*data[j-2];
                    for (int k = 1; k < (j-1); k++)
                        data[j-k-1] = div*((dr+(real4) k) *data[j-k-2] + (-dr+(real4) (j-k))*data[j-k-1]);
                    data[0] = div*(- dr+1.0f)*data[0];
                }
                data[PME_ORDER-1] = scale*dr*data[PME_ORDER-2];
                for (int j = 1; j < (PME_ORDER-1); j++)
                    data[PME_ORDER-j-1] = scale*((dr+(real4) j)*data[PME_ORDER-j-2] + (-dr+(real4) (PME_ORDER-j))*data[PME_ORDER-j-1]);
                data[0] = scale*(-dr+1.0f)*data[0];
            }
            for (int iy = 0; iy < PME_ORDER; iy++) {
                int yindex = gridIndex.y+iy;
                yindex -= (yindex >= GRID_SIZE_Y ? GRID_SIZE_Y : 0);
                for (int iz = 0; iz < PME_ORDER; iz++) {
                    int zindex = gridIndex.z+iz;
                    zindex -= (zindex >= GRID_SIZE_Z ? GRID_SIZE_Z : 0);
                    int index = xindex*GRID_SIZE_Y*GRID_SIZE_Z + yindex*GRID_SIZE_Z + zindex;
                    pmeGrid[index] += charge*data[ix].x*data[iy].y*data[iz].z;
                }
            }
        }
    }
}
#else
__kernel void gridSpreadCharge(__global const real4* restrict posq, __global const int2* restrict pmeAtomGridIndex, __global const int* restrict pmeAtomRange,
        __global real* restrict pmeGrid, __global const real4* restrict pmeBsplineTheta
#ifdef USE_LJPME
        , __global const float2* restrict sigmaEpsilon
#else
        , __global const real* restrict charges
#endif
    ) {
    unsigned int numGridPoints = GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z;
    for (int gridIndex = get_global_id(0); gridIndex < numGridPoints; gridIndex += get_global_size(0)) {
        // Compute the charge on a grid point.

        int4 gridPoint;
        gridPoint.x = gridIndex/(GRID_SIZE_Y*GRID_SIZE_Z);
        int remainder = gridIndex-gridPoint.x*GRID_SIZE_Y*GRID_SIZE_Z;
        gridPoint.y = remainder/GRID_SIZE_Z;
        gridPoint.z = remainder-gridPoint.y*GRID_SIZE_Z;
        real result = 0.0f;

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
                    real atomCharge = pmeBsplineTheta[atomIndex+ix*NUM_ATOMS].w;
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
                        real atomCharge = pmeBsplineTheta[atomIndex+ix*NUM_ATOMS].w;
                        result += atomCharge*pmeBsplineTheta[atomIndex+ix*NUM_ATOMS].x*pmeBsplineTheta[atomIndex+iy*NUM_ATOMS].y*pmeBsplineTheta[atomIndex+iz*NUM_ATOMS].z;
                    }
                }
            }
        }
        pmeGrid[gridIndex] = result*EPSILON_FACTOR;
    }
}
#endif

__kernel void reciprocalConvolution(__global real2* restrict pmeGrid, __global const real* restrict pmeBsplineModuliX,
        __global const real* restrict pmeBsplineModuliY, __global const real* restrict pmeBsplineModuliZ, real4 recipBoxVecX, real4 recipBoxVecY, real4 recipBoxVecZ) {
    // R2C stores into a half complex matrix where the last dimension is cut by half
    const unsigned int gridSize = GRID_SIZE_X*GRID_SIZE_Y*(GRID_SIZE_Z/2+1);
#ifdef USE_LJPME
    const real recipScaleFactor = -(2*M_PI/6)*SQRT(M_PI)*recipBoxVecX.x*recipBoxVecY.y*recipBoxVecZ.z;
    real bfac = M_PI / EWALD_ALPHA;
    real fac1 = 2*M_PI*M_PI*M_PI*SQRT(M_PI);
    real fac2 = EWALD_ALPHA*EWALD_ALPHA*EWALD_ALPHA;
    real fac3 = -2*EWALD_ALPHA*M_PI*M_PI;
#else
    const real recipScaleFactor = (1.0f/M_PI)*recipBoxVecX.x*recipBoxVecY.y*recipBoxVecZ.z;
#endif

    for (int index = get_global_id(0); index < gridSize; index += get_global_size(0)) {
        // real indices
        int kx = index/(GRID_SIZE_Y*(GRID_SIZE_Z/2+1));
        int remainder = index-kx*GRID_SIZE_Y*(GRID_SIZE_Z/2+1);
        int ky = remainder/(GRID_SIZE_Z/2+1);
        int kz = remainder-ky*(GRID_SIZE_Z/2+1);
        int mx = (kx < (GRID_SIZE_X+1)/2) ? kx : (kx-GRID_SIZE_X);
        int my = (ky < (GRID_SIZE_Y+1)/2) ? ky : (ky-GRID_SIZE_Y);
        int mz = (kz < (GRID_SIZE_Z+1)/2) ? kz : (kz-GRID_SIZE_Z);
        real mhx = mx*recipBoxVecX.x;
        real mhy = mx*recipBoxVecY.x+my*recipBoxVecY.y;
        real mhz = mx*recipBoxVecZ.x+my*recipBoxVecZ.y+mz*recipBoxVecZ.z;
        real bx = pmeBsplineModuliX[kx];
        real by = pmeBsplineModuliY[ky];
        real bz = pmeBsplineModuliZ[kz];
        real2 grid = pmeGrid[index];
        real m2 = mhx*mhx+mhy*mhy+mhz*mhz;
#ifdef USE_LJPME
        real denom = recipScaleFactor/(bx*by*bz);
        real m = SQRT(m2);
        real m3 = m*m2;
        real b = bfac*m;
        real expfac = -b*b;
        real expterm = EXP(expfac);
        real erfcterm = erfc(b);
        real eterm = (fac1*erfcterm*m3 + expterm*(fac2 + fac3*m2)) * denom;
        pmeGrid[index] = (real2) (grid.x*eterm, grid.y*eterm);
#else
        real denom = m2*bx*by*bz;
        real eterm = recipScaleFactor*EXP(-RECIP_EXP_FACTOR*m2)/denom;
        if (kx != 0 || ky != 0 || kz != 0) {
            pmeGrid[index] = (real2) (grid.x*eterm, grid.y*eterm);
        }
#endif
    }
}

__kernel void gridEvaluateEnergy(__global real2* restrict pmeGrid, __global mixed* restrict energyBuffer,
                      __global const real* restrict pmeBsplineModuliX, __global const real* restrict pmeBsplineModuliY, __global const real* restrict pmeBsplineModuliZ,
                      real4 recipBoxVecX, real4 recipBoxVecY, real4 recipBoxVecZ) {
    // R2C stores into a half complex matrix where the last dimension is cut by half
    const unsigned int gridSize = GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z;
 #ifdef USE_LJPME
    const real recipScaleFactor = -(2*M_PI/6)*SQRT(M_PI)*recipBoxVecX.x*recipBoxVecY.y*recipBoxVecZ.z;
    real bfac = M_PI / EWALD_ALPHA;
    real fac1 = 2*M_PI*M_PI*M_PI*SQRT(M_PI);
    real fac2 = EWALD_ALPHA*EWALD_ALPHA*EWALD_ALPHA;
    real fac3 = -2*EWALD_ALPHA*M_PI*M_PI;
#else
    const real recipScaleFactor = (1.0f/M_PI)*recipBoxVecX.x*recipBoxVecY.y*recipBoxVecZ.z;
#endif
 
    mixed energy = 0;
    for (int index = get_global_id(0); index < gridSize; index += get_global_size(0)) {
        // real indices
        int kx = index/(GRID_SIZE_Y*(GRID_SIZE_Z));
        int remainder = index-kx*GRID_SIZE_Y*(GRID_SIZE_Z);
        int ky = remainder/(GRID_SIZE_Z);
        int kz = remainder-ky*(GRID_SIZE_Z);
        int mx = (kx < (GRID_SIZE_X+1)/2) ? kx : (kx-GRID_SIZE_X);
        int my = (ky < (GRID_SIZE_Y+1)/2) ? ky : (ky-GRID_SIZE_Y);
        int mz = (kz < (GRID_SIZE_Z+1)/2) ? kz : (kz-GRID_SIZE_Z);
        real mhx = mx*recipBoxVecX.x;
        real mhy = mx*recipBoxVecY.x+my*recipBoxVecY.y;
        real mhz = mx*recipBoxVecZ.x+my*recipBoxVecZ.y+mz*recipBoxVecZ.z;
        real m2 = mhx*mhx+mhy*mhy+mhz*mhz;
        real bx = pmeBsplineModuliX[kx];
        real by = pmeBsplineModuliY[ky];
        real bz = pmeBsplineModuliZ[kz];
#ifdef USE_LJPME
        real denom = recipScaleFactor/(bx*by*bz);
        real m = SQRT(m2);
        real m3 = m*m2;
        real b = bfac*m;
        real expfac = -b*b;
        real expterm = EXP(expfac);
        real erfcterm = erfc(b);
        real eterm = (fac1*erfcterm*m3 + expterm*(fac2 + fac3*m2)) * denom;
#else
        real denom = m2*bx*by*bz;
        real eterm = recipScaleFactor*EXP(-RECIP_EXP_FACTOR*m2)/denom;
#endif
        if (kz >= (GRID_SIZE_Z/2+1)) {
            kx = ((kx == 0) ? kx : GRID_SIZE_X-kx);
            ky = ((ky == 0) ? ky : GRID_SIZE_Y-ky);
            kz = GRID_SIZE_Z-kz;
        } 
        int indexInHalfComplexGrid = kz + ky*(GRID_SIZE_Z/2+1)+kx*(GRID_SIZE_Y*(GRID_SIZE_Z/2+1));
        real2 grid = pmeGrid[indexInHalfComplexGrid];
#ifndef USE_LJPME
        if (kx != 0 || ky != 0 || kz != 0)
#endif
            energy += eterm*(grid.x*grid.x + grid.y*grid.y);
    }
#if defined(USE_PME_STREAM) && !defined(USE_LJPME)
    energyBuffer[get_global_id(0)] = 0.5f*energy;
#else
    energyBuffer[get_global_id(0)] += 0.5f*energy;
#endif
}

__kernel void gridInterpolateForce(__global const real4* restrict posq, __global real4* restrict forceBuffers, __global const real* restrict pmeGrid,
        real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, real4 recipBoxVecX,
        real4 recipBoxVecY, real4 recipBoxVecZ, __global int2* restrict pmeAtomGridIndex
#ifdef USE_LJPME
        , __global const float2* restrict sigmaEpsilon
#else
        , __global const real* restrict charges
#endif
    ) {
    const real scale = 1/(real) (PME_ORDER-1);
    real4 data[PME_ORDER];
    real4 ddata[PME_ORDER];
    
    // Process the atoms in spatially sorted order.  This improves cache performance when loading
    // the grid values.
    
    for (int i = get_global_id(0); i < NUM_ATOMS; i += get_global_size(0)) {
        int atom = pmeAtomGridIndex[i].x;
        real4 force = 0.0f;
        real4 pos = posq[atom];
        APPLY_PERIODIC_TO_POS(pos)
        real3 t = (real3) (pos.x*recipBoxVecX.x+pos.y*recipBoxVecY.x+pos.z*recipBoxVecZ.x,
                           pos.y*recipBoxVecY.y+pos.z*recipBoxVecZ.y,
                           pos.z*recipBoxVecZ.z);
        t.x = (t.x-floor(t.x))*GRID_SIZE_X;
        t.y = (t.y-floor(t.y))*GRID_SIZE_Y;
        t.z = (t.z-floor(t.z))*GRID_SIZE_Z;
        int4 gridIndex = (int4) (((int) t.x) % GRID_SIZE_X,
                                 ((int) t.y) % GRID_SIZE_Y,
                                 ((int) t.z) % GRID_SIZE_Z, 0);

        // Since we need the full set of thetas, it's faster to compute them here than load them
        // from global memory.

        real4 dr = (real4) (t.x-(int) t.x, t.y-(int) t.y, t.z-(int) t.z, 0.0f);
        data[PME_ORDER-1] = 0.0f;
        data[1] = dr;
        data[0] = 1.0f-dr;
        for (int j = 3; j < PME_ORDER; j++) {
            real div = RECIP(j-1.0f);
            data[j-1] = div*dr*data[j-2];
            for (int k = 1; k < (j-1); k++)
                data[j-k-1] = div*((dr+(real4) k) *data[j-k-2] + (-dr+(real4) (j-k))*data[j-k-1]);
            data[0] = div*(-dr+1.0f)*data[0];
        }
        ddata[0] = -data[0];
        for (int j = 1; j < PME_ORDER; j++)
            ddata[j] = data[j-1]-data[j];
        data[PME_ORDER-1] = scale*dr*data[PME_ORDER-2];
        for (int j = 1; j < (PME_ORDER-1); j++)
            data[PME_ORDER-j-1] = scale*((dr+(real4) j)*data[PME_ORDER-j-2] + (-dr+(real4) (PME_ORDER-j))*data[PME_ORDER-j-1]);
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
                    real gridvalue = pmeGrid[index];
                    force.x += ddata[ix].x*data[iy].y*data[iz].z*gridvalue;
                    force.y += data[ix].x*ddata[iy].y*data[iz].z*gridvalue;
                    force.z += data[ix].x*data[iy].y*ddata[iz].z*gridvalue;
                }
            }
        }
        real4 totalForce = forceBuffers[atom];
#ifdef USE_LJPME
        const float2 sigEps = sigmaEpsilon[atom];
        real q = 8*sigEps.x*sigEps.x*sigEps.x*sigEps.y;
#else
        real q = CHARGE*EPSILON_FACTOR;
#endif
        totalForce.x -= q*(force.x*GRID_SIZE_X*recipBoxVecX.x);
        totalForce.y -= q*(force.x*GRID_SIZE_X*recipBoxVecY.x+force.y*GRID_SIZE_Y*recipBoxVecY.y);
        totalForce.z -= q*(force.x*GRID_SIZE_X*recipBoxVecZ.x+force.y*GRID_SIZE_Y*recipBoxVecZ.y+force.z*GRID_SIZE_Z*recipBoxVecZ.z);
        forceBuffers[atom] = totalForce;
    }
}

__kernel void addForces(__global const real4* restrict forces, __global real4* restrict forceBuffers) {
    for (int atom = get_global_id(0); atom < NUM_ATOMS; atom += get_global_size(0))
        forceBuffers[atom] += forces[atom];
}

__kernel void addEnergy(__global const mixed* restrict pmeEnergyBuffer, __global mixed* restrict energyBuffer, int bufferSize) {
    for (int i = get_global_id(0); i < bufferSize; i += get_global_size(0))
        energyBuffer[i] += pmeEnergyBuffer[i];
}
