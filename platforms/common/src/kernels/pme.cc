KERNEL void findAtomGridIndex(GLOBAL const real4* RESTRICT posq, GLOBAL int2* RESTRICT pmeAtomGridIndex,
        real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        real4 recipBoxVecX, real4 recipBoxVecY, real4 recipBoxVecZ
#ifndef SUPPORTS_64_BIT_ATOMICS
        , GLOBAL real4* RESTRICT pmeBsplineTheta, LOCAL real4* RESTRICT bsplinesCache,
#ifdef CHARGE_FROM_SIGEPS
        GLOBAL const float2* RESTRICT sigmaEpsilon
#else
        GLOBAL const real* RESTRICT charges
#endif
#endif
    ) {
    // Compute the index of the grid point each atom is associated with.

    for (int atom = GLOBAL_ID; atom < NUM_ATOMS; atom += GLOBAL_SIZE) {
        real4 pos = posq[atom];
        APPLY_PERIODIC_TO_POS(pos)
        real3 t = make_real3(pos.x*recipBoxVecX.x+pos.y*recipBoxVecY.x+pos.z*recipBoxVecZ.x,
                             pos.y*recipBoxVecY.y+pos.z*recipBoxVecZ.y,
                             pos.z*recipBoxVecZ.z);
        t.x = (t.x-floor(t.x))*GRID_SIZE_X;
        t.y = (t.y-floor(t.y))*GRID_SIZE_Y;
        t.z = (t.z-floor(t.z))*GRID_SIZE_Z;
        int3 gridIndex = make_int3(((int) t.x) % GRID_SIZE_X,
                                   ((int) t.y) % GRID_SIZE_Y,
                                   ((int) t.z) % GRID_SIZE_Z);
        pmeAtomGridIndex[atom] = make_int2(atom, gridIndex.x*GRID_SIZE_Y*GRID_SIZE_Z+gridIndex.y*GRID_SIZE_Z+gridIndex.z);
#ifndef SUPPORTS_64_BIT_ATOMICS
        // Compute B-splines here for use in the charge spreading kernel.
        const real4 scale = 1/(real) (PME_ORDER-1);
        LOCAL real4* data = &bsplinesCache[LOCAL_ID*PME_ORDER];
        real4 dr = (real4) (t.x-(int) t.x, t.y-(int) t.y, t.z-(int) t.z, 0.0f);
        data[PME_ORDER-1] = 0.0f;
        data[1] = dr;
        data[0] = 1.0f-dr;
        for (int j = 3; j < PME_ORDER; j++) {
            real div = RECIP(j-1.0f);
            data[j-1] = div*dr*data[j-2];
            for (int k = 1; k < (j-1); k++)
                data[j-k-1] = div*((dr+make_real4(k))*data[j-k-2] + (-dr+make_real4(j-k))*data[j-k-1]);
            data[0] = div*(- dr+1.0f)*data[0];
        }
        data[PME_ORDER-1] = scale*dr*data[PME_ORDER-2];
        for (int j = 1; j < (PME_ORDER-1); j++)
            data[PME_ORDER-j-1] = scale*((dr+make_real4(j))*data[PME_ORDER-j-2] + (-dr+make_real4(PME_ORDER-j))*data[PME_ORDER-j-1]);
        data[0] = scale*(-dr+1.0f)*data[0];
        for (int j = 0; j < PME_ORDER; j++) {
#ifdef CHARGE_FROM_SIGEPS
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

#ifdef SUPPORTS_64_BIT_ATOMICS
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable

KERNEL void gridSpreadCharge(GLOBAL const real4* RESTRICT posq,
#ifdef USE_FIXED_POINT_CHARGE_SPREADING
        GLOBAL mm_ulong* RESTRICT pmeGrid,
#else
        GLOBAL real* RESTRICT pmeGrid,
#endif
        real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        real4 recipBoxVecX, real4 recipBoxVecY, real4 recipBoxVecZ, GLOBAL const int2* RESTRICT pmeAtomGridIndex,
#ifdef CHARGE_FROM_SIGEPS
        GLOBAL const float2* RESTRICT sigmaEpsilon
#else
        GLOBAL const real* RESTRICT charges
#endif
        ) {
    // To improve memory efficiency, we divide indices along the z axis into
    // PME_ORDER blocks, where the data for each block is stored together.  We
    // can ensure that all threads write to the same block at the same time,
    // which leads to better coalescing of writes.
    
    LOCAL int zindexTable[GRID_SIZE_Z+PME_ORDER];
    int blockSize = (int) ceil(GRID_SIZE_Z/(real) PME_ORDER);
    for (int i = LOCAL_ID; i < GRID_SIZE_Z+PME_ORDER; i += LOCAL_SIZE) {
        int zindex = i % GRID_SIZE_Z;
	int block = zindex % PME_ORDER;
        zindexTable[i] = zindex/PME_ORDER + block*GRID_SIZE_X*GRID_SIZE_Y*blockSize;
    }
    SYNC_THREADS;
    
    // Process the atoms in spatially sorted order.  This improves efficiency when writing
    // the grid values.
    
    real3 data[PME_ORDER];
    const real scale = RECIP((real) (PME_ORDER-1));
    for (int i = GLOBAL_ID; i < NUM_ATOMS; i += GLOBAL_SIZE) {
        int atom = pmeAtomGridIndex[i].x;
        real4 pos = posq[atom];
#ifdef CHARGE_FROM_SIGEPS
        const float2 sigEps = sigmaEpsilon[atom];
        const real charge = 8*sigEps.x*sigEps.x*sigEps.x*sigEps.y;
#else
        const real charge = (CHARGE)*EPSILON_FACTOR;
#endif
        APPLY_PERIODIC_TO_POS(pos)
        real3 t = make_real3(pos.x*recipBoxVecX.x+pos.y*recipBoxVecY.x+pos.z*recipBoxVecZ.x,
                             pos.y*recipBoxVecY.y+pos.z*recipBoxVecZ.y,
                             pos.z*recipBoxVecZ.z);
        t.x = (t.x-floor(t.x))*GRID_SIZE_X;
        t.y = (t.y-floor(t.y))*GRID_SIZE_Y;
        t.z = (t.z-floor(t.z))*GRID_SIZE_Z;
        int3 gridIndex = make_int3(((int) t.x) % GRID_SIZE_X,
                                   ((int) t.y) % GRID_SIZE_Y,
                                   ((int) t.z) % GRID_SIZE_Z);
        if (charge == 0)
            continue;

        // Since we need the full set of thetas, it's faster to compute them here than load them
        // from global memory.

        real3 dr = make_real3(t.x-(int) t.x, t.y-(int) t.y, t.z-(int) t.z);
        data[PME_ORDER-1] = make_real3(0);
        data[1] = dr;
        data[0] = make_real3(1)-dr;
        for (int j = 3; j < PME_ORDER; j++) {
            real div = RECIP((real) (j-1));
            data[j-1] = div*dr*data[j-2];
            for (int k = 1; k < (j-1); k++)
                data[j-k-1] = div*((dr+make_real3(k))*data[j-k-2] + (make_real3(j-k)-dr)*data[j-k-1]);
            data[0] = div*(make_real3(1)-dr)*data[0];
        }
        data[PME_ORDER-1] = scale*dr*data[PME_ORDER-2];
        for (int j = 1; j < (PME_ORDER-1); j++)
            data[PME_ORDER-j-1] = scale*((dr+make_real3(j))*data[PME_ORDER-j-2] + (make_real3(PME_ORDER-j)-dr)*data[PME_ORDER-j-1]);
        data[0] = scale*(make_real3(1)-dr)*data[0];

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
#ifdef USE_FIXED_POINT_CHARGE_SPREADING
                    ATOMIC_ADD(&pmeGrid[index], (mm_ulong) ((mm_long) (add*0x100000000)));
#else
                    ATOMIC_ADD(&pmeGrid[index], add);
#endif
                }
            }
        }
    }
}

KERNEL void finishSpreadCharge(
#ifdef USE_FIXED_POINT_CHARGE_SPREADING
        GLOBAL const mm_long* RESTRICT grid1,
#else
        GLOBAL const real* RESTRICT grid1,
#endif
        GLOBAL real* RESTRICT grid2) {
    // During charge spreading, we shuffled the order of indices along the z
    // axis to make memory access more efficient.  We now need to unshuffle
    // them.  If the values were accumulated as fixed point, we also need to
    // convert them to floating point.

    LOCAL int zindexTable[GRID_SIZE_Z];
    int blockSize = (int) ceil(GRID_SIZE_Z/(real) PME_ORDER);
    for (int i = LOCAL_ID; i < GRID_SIZE_Z; i += LOCAL_SIZE) {
	int block = i % PME_ORDER;
        zindexTable[i] = i/PME_ORDER + block*GRID_SIZE_X*GRID_SIZE_Y*blockSize;
    }
    SYNC_THREADS;
    const unsigned int gridSize = GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z;
    real scale = 1/(real) 0x100000000;
    for (int index = GLOBAL_ID; index < gridSize; index += GLOBAL_SIZE) {
        int zindex = index%GRID_SIZE_Z;
        int loadIndex = zindexTable[zindex] + blockSize*(int) (index/GRID_SIZE_Z);
#ifdef USE_FIXED_POINT_CHARGE_SPREADING
        grid2[index] = scale*grid1[loadIndex];
#else
        grid2[index] = grid1[loadIndex];
#endif
    }
}
#elif defined(DEVICE_IS_CPU)
KERNEL void gridSpreadCharge(GLOBAL const real4* RESTRICT posq, GLOBAL real* RESTRICT pmeGrid,
        real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        real4 recipBoxVecX, real4 recipBoxVecY, real4 recipBoxVecZ,
#ifdef CHARGE_FROM_SIGEPS
        GLOBAL const float2* RESTRICT sigmaEpsilon
#else
        GLOBAL const real* RESTRICT charges
#endif
    ) {
    const int firstx = GLOBAL_ID*GRID_SIZE_X/GLOBAL_SIZE;
    const int lastx = (GLOBAL_ID+1)*GRID_SIZE_X/GLOBAL_SIZE;
    if (firstx == lastx)
        return;
    const real4 scale = 1/(real) (PME_ORDER-1);
    real4 data[PME_ORDER];
    
    // Process the atoms in spatially sorted order.  This improves efficiency when writing
    // the grid values.
    
    for (int i = 0; i < NUM_ATOMS; i++) {
        int atom = i;
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

#ifdef CHARGE_FROM_SIGEPS
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
/**
 * For each grid point, find the range of sorted atoms associated with that point.
 */
KERNEL void findAtomRangeForGrid(GLOBAL int2* RESTRICT pmeAtomGridIndex, GLOBAL int* RESTRICT pmeAtomRange, GLOBAL const real4* RESTRICT posq) {
    int start = (NUM_ATOMS*GLOBAL_ID)/GLOBAL_SIZE;
    int end = (NUM_ATOMS*(GLOBAL_ID+1))/GLOBAL_SIZE;
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

    if (GLOBAL_ID == GLOBAL_SIZE-1) {
        int gridSize = GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z;
        for (int j = last+1; j <= gridSize; ++j)
            pmeAtomRange[j] = NUM_ATOMS;
    }
}

/**
 * The grid index won't be needed again.  Reuse that component to hold the z index, thus saving
 * some work in the charge spreading kernel.
 */
KERNEL void recordZIndex(GLOBAL int2* RESTRICT pmeAtomGridIndex, GLOBAL const real4* RESTRICT posq, real4 periodicBoxSize, real4 recipBoxVecZ) {
    int start = (NUM_ATOMS*GLOBAL_ID)/GLOBAL_SIZE;
    int end = (NUM_ATOMS*(GLOBAL_ID+1))/GLOBAL_SIZE;
    for (int i = start; i < end; ++i) {
        real posz = posq[pmeAtomGridIndex[i].x].z;
        posz -= floor(posz*recipBoxVecZ.z)*periodicBoxSize.z;
        int z = ((int) ((posz*recipBoxVecZ.z)*GRID_SIZE_Z)) % GRID_SIZE_Z;
        pmeAtomGridIndex[i].y = z;
    }
}

KERNEL void gridSpreadCharge(GLOBAL const real4* RESTRICT posq, GLOBAL real* RESTRICT pmeGrid,
        GLOBAL const int2* RESTRICT pmeAtomGridIndex, GLOBAL const int* RESTRICT pmeAtomRange,
        GLOBAL const real4* RESTRICT pmeBsplineTheta
#ifdef CHARGE_FROM_SIGEPS
        , GLOBAL const float2* RESTRICT sigmaEpsilon
#else
        , GLOBAL const real* RESTRICT charges
#endif
    ) {
    unsigned int numGridPoints = GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z;
    for (int gridIndex = GLOBAL_ID; gridIndex < numGridPoints; gridIndex += GLOBAL_SIZE) {
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

KERNEL void reciprocalConvolution(GLOBAL real2* RESTRICT pmeGrid, GLOBAL const real* RESTRICT pmeBsplineModuliX,
        GLOBAL const real* RESTRICT pmeBsplineModuliY, GLOBAL const real* RESTRICT pmeBsplineModuliZ,
        real4 recipBoxVecX, real4 recipBoxVecY, real4 recipBoxVecZ) {
    // R2C stores into a half complex matrix where the last dimension is cut by half
    const unsigned int gridSize = GRID_SIZE_X*GRID_SIZE_Y*(GRID_SIZE_Z/2+1);
#ifdef USE_LJPME
    const real recipScaleFactor = -(2*M_PI/6)*SQRT(M_PI)*recipBoxVecX.x*recipBoxVecY.y*recipBoxVecZ.z;
    real bfac = M_PI / EWALD_ALPHA;
    real fac1 = 2*M_PI*M_PI*M_PI*SQRT(M_PI);
    real fac2 = EWALD_ALPHA*EWALD_ALPHA*EWALD_ALPHA;
    real fac3 = -2*EWALD_ALPHA*M_PI*M_PI;
#else
    const real recipScaleFactor = RECIP(M_PI)*recipBoxVecX.x*recipBoxVecY.y*recipBoxVecZ.z;
#endif

    for (int index = GLOBAL_ID; index < gridSize; index += GLOBAL_SIZE) {
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
        real erfcterm = ERFC(b);
        real eterm = (fac1*erfcterm*m3 + expterm*(fac2 + fac3*m2)) * denom;
        pmeGrid[index] = make_real2(grid.x*eterm, grid.y*eterm);
#else
        real denom = m2*bx*by*bz;
        real eterm = recipScaleFactor*EXP(-RECIP_EXP_FACTOR*m2)/denom;
        if (kx != 0 || ky != 0 || kz != 0) {
            pmeGrid[index] = make_real2(grid.x*eterm, grid.y*eterm);
        }
#endif
    }
}

KERNEL void gridEvaluateEnergy(GLOBAL real2* RESTRICT pmeGrid, GLOBAL mixed* RESTRICT energyBuffer,
                      GLOBAL const real* RESTRICT pmeBsplineModuliX, GLOBAL const real* RESTRICT pmeBsplineModuliY, GLOBAL const real* RESTRICT pmeBsplineModuliZ,
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
    const real recipScaleFactor = RECIP(M_PI)*recipBoxVecX.x*recipBoxVecY.y*recipBoxVecZ.z;
#endif

    mixed energy = 0;
    for (int index = GLOBAL_ID; index < gridSize; index += GLOBAL_SIZE) {
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
        real erfcterm = ERFC(b);
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
    energyBuffer[GLOBAL_ID] = 0.5f*energy;
#else
    energyBuffer[GLOBAL_ID] += 0.5f*energy;
#endif
}

KERNEL void gridInterpolateForce(GLOBAL const real4* RESTRICT posq, GLOBAL mm_ulong* RESTRICT forceBuffers, GLOBAL const real* RESTRICT pmeGrid,
        real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        real4 recipBoxVecX, real4 recipBoxVecY, real4 recipBoxVecZ, GLOBAL const int2* RESTRICT pmeAtomGridIndex,
#ifdef CHARGE_FROM_SIGEPS
        GLOBAL const float2* RESTRICT sigmaEpsilon
#else
        GLOBAL const real* RESTRICT charges
#endif
        ) {
    real3 data[PME_ORDER];
    real3 ddata[PME_ORDER];
    const real scale = RECIP((real) (PME_ORDER-1));
    
    // Process the atoms in spatially sorted order.  This improves cache performance when loading
    // the grid values.
    
    for (int i = GLOBAL_ID; i < NUM_ATOMS; i += GLOBAL_SIZE) {
        int atom = pmeAtomGridIndex[i].x;
        real3 force = make_real3(0);
        real4 pos = posq[atom];
        APPLY_PERIODIC_TO_POS(pos)
        real3 t = make_real3(pos.x*recipBoxVecX.x+pos.y*recipBoxVecY.x+pos.z*recipBoxVecZ.x,
                             pos.y*recipBoxVecY.y+pos.z*recipBoxVecZ.y,
                             pos.z*recipBoxVecZ.z);
        t.x = (t.x-floor(t.x))*GRID_SIZE_X;
        t.y = (t.y-floor(t.y))*GRID_SIZE_Y;
        t.z = (t.z-floor(t.z))*GRID_SIZE_Z;
        int3 gridIndex = make_int3(((int) t.x) % GRID_SIZE_X,
                                   ((int) t.y) % GRID_SIZE_Y,
                                   ((int) t.z) % GRID_SIZE_Z);

        // Since we need the full set of thetas, it's faster to compute them here than load them
        // from global memory.

        real3 dr = make_real3(t.x-(int) t.x, t.y-(int) t.y, t.z-(int) t.z);
        data[PME_ORDER-1] = make_real3(0);
        data[1] = dr;
        data[0] = make_real3(1)-dr;
        for (int j = 3; j < PME_ORDER; j++) {
            real div = RECIP((real) (j-1));
            data[j-1] = div*dr*data[j-2];
            for (int k = 1; k < (j-1); k++)
                data[j-k-1] = div*((dr+make_real3(k))*data[j-k-2] + (make_real3(j-k)-dr)*data[j-k-1]);
            data[0] = div*(make_real3(1)-dr)*data[0];
        }
        ddata[0] = -data[0];
        for (int j = 1; j < PME_ORDER; j++)
            ddata[j] = data[j-1]-data[j];
        data[PME_ORDER-1] = scale*dr*data[PME_ORDER-2];
        for (int j = 1; j < (PME_ORDER-1); j++)
            data[PME_ORDER-j-1] = scale*((dr+make_real3(j))*data[PME_ORDER-j-2] + (make_real3(PME_ORDER-j)-dr)*data[PME_ORDER-j-1]);
        data[0] = scale*(make_real3(1)-dr)*data[0];

        // Compute the force on this atom.

        for (int ix = 0; ix < PME_ORDER; ix++) {
            int xbase = gridIndex.x+ix;
            xbase -= (xbase >= GRID_SIZE_X ? GRID_SIZE_X : 0);
            xbase = xbase*GRID_SIZE_Y*GRID_SIZE_Z;
            real dx = data[ix].x;
            real ddx = ddata[ix].x;
            
            for (int iy = 0; iy < PME_ORDER; iy++) {
                int ybase = gridIndex.y+iy;
                ybase -= (ybase >= GRID_SIZE_Y ? GRID_SIZE_Y : 0);
                ybase = xbase + ybase*GRID_SIZE_Z;
                real dy = data[iy].y;
                real ddy = ddata[iy].y;
                
                for (int iz = 0; iz < PME_ORDER; iz++) {
                    int zindex = gridIndex.z+iz;
                    zindex -= (zindex >= GRID_SIZE_Z ? GRID_SIZE_Z : 0);
                    int index = ybase + zindex;
                    real gridvalue = pmeGrid[index];
                    force.x += ddx*dy*data[iz].z*gridvalue;
                    force.y += dx*ddy*data[iz].z*gridvalue;
                    force.z += dx*dy*ddata[iz].z*gridvalue;
                }
            }
        }
#ifdef CHARGE_FROM_SIGEPS
        const float2 sigEps = sigmaEpsilon[atom];
        real q = 8*sigEps.x*sigEps.x*sigEps.x*sigEps.y;
#else
        real q = CHARGE*EPSILON_FACTOR;
#endif
        real forceX = -q*(force.x*GRID_SIZE_X*recipBoxVecX.x);
        real forceY = -q*(force.x*GRID_SIZE_X*recipBoxVecY.x+force.y*GRID_SIZE_Y*recipBoxVecY.y);
        real forceZ = -q*(force.x*GRID_SIZE_X*recipBoxVecZ.x+force.y*GRID_SIZE_Y*recipBoxVecZ.y+force.z*GRID_SIZE_Z*recipBoxVecZ.z);
#ifdef USE_PME_STREAM
        ATOMIC_ADD(&forceBuffers[atom], (mm_ulong) ((mm_long) (forceX*0x100000000)));
        ATOMIC_ADD(&forceBuffers[atom+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (forceY*0x100000000)));
        ATOMIC_ADD(&forceBuffers[atom+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (forceZ*0x100000000)));
#else
        forceBuffers[atom] += (mm_ulong) ((mm_long) (forceX*0x100000000));
        forceBuffers[atom+PADDED_NUM_ATOMS] += (mm_ulong) ((mm_long) (forceY*0x100000000));
        forceBuffers[atom+2*PADDED_NUM_ATOMS] += (mm_ulong) ((mm_long) (forceZ*0x100000000));
#endif
    }
}

KERNEL void addForces(GLOBAL const real4* RESTRICT forces, GLOBAL mm_long* RESTRICT forceBuffers) {
    for (int atom = GLOBAL_ID; atom < NUM_ATOMS; atom += GLOBAL_SIZE) {
        real4 f = forces[atom];
        forceBuffers[atom] += (mm_long) (f.x*0x100000000);
        forceBuffers[atom+PADDED_NUM_ATOMS] += (mm_long) (f.y*0x100000000);
        forceBuffers[atom+2*PADDED_NUM_ATOMS] += (mm_long) (f.z*0x100000000);
    }
}

KERNEL void addEnergy(GLOBAL const mixed* RESTRICT pmeEnergyBuffer, GLOBAL mixed* RESTRICT energyBuffer, int bufferSize) {
    for (int i = GLOBAL_ID; i < bufferSize; i += GLOBAL_SIZE)
        energyBuffer[i] += pmeEnergyBuffer[i];
}
