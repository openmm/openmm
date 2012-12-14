extern "C" __global__ void updateBsplines(const real4* __restrict__ posq, real4* __restrict__ pmeBsplineTheta, int2* __restrict__ pmeAtomGridIndex,
            real4 periodicBoxSize, real4 invPeriodicBoxSize) {
    extern __shared__ real3 bsplinesCache[];
    real3* data = &bsplinesCache[threadIdx.x*PME_ORDER];
    const real3 scale = make_real3(RECIP(PME_ORDER-1));
    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < NUM_ATOMS; i += blockDim.x*gridDim.x) {
        real4 pos = posq[i];
        pos.x -= floor(pos.x*invPeriodicBoxSize.x)*periodicBoxSize.x;
        pos.y -= floor(pos.y*invPeriodicBoxSize.y)*periodicBoxSize.y;
        pos.z -= floor(pos.z*invPeriodicBoxSize.z)*periodicBoxSize.z;
        real3 t = make_real3((pos.x*invPeriodicBoxSize.x)*GRID_SIZE_X,
                             (pos.y*invPeriodicBoxSize.y)*GRID_SIZE_Y,
                             (pos.z*invPeriodicBoxSize.z)*GRID_SIZE_Z);
        real3 dr = make_real3(t.x-(int) t.x, t.y-(int) t.y, t.z-(int) t.z);
        int3 gridIndex = make_int3(((int) t.x) % GRID_SIZE_X,
                                 ((int) t.y) % GRID_SIZE_Y,
                                 ((int) t.z) % GRID_SIZE_Z);
        pmeAtomGridIndex[i] = make_int2(i, gridIndex.x*GRID_SIZE_Y*GRID_SIZE_Z+gridIndex.y*GRID_SIZE_Z+gridIndex.z);
        data[PME_ORDER-1] = make_real3(0);
        data[1] = dr;
        data[0] = make_real3(1)-dr;
        for (int j = 3; j < PME_ORDER; j++) {
            real div = RECIP(j-1);
            data[j-1] = div*dr*data[j-2];
            for (int k = 1; k < (j-1); k++)
                data[j-k-1] = div*((dr+make_real3(k)) *data[j-k-2] + (make_real3(j-k)-dr)*data[j-k-1]);
            data[0] = div*(make_real3(1)-dr)*data[0];
        }
        data[PME_ORDER-1] = scale*dr*data[PME_ORDER-2];
        for (int j = 1; j < (PME_ORDER-1); j++)
            data[PME_ORDER-j-1] = scale*((dr+make_real3(j))*data[PME_ORDER-j-2] + (make_real3(PME_ORDER-j)-dr)*data[PME_ORDER-j-1]);
        data[0] = scale*(make_real3(1)-dr)*data[0];
        for (int j = 0; j < PME_ORDER; j++) {
            real3 d = data[j]; // Copy it as a workaround for a bug in CUDA 5.0
            pmeBsplineTheta[i+j*NUM_ATOMS] = make_real4(d.x, d.y, d.z, pos.w);  // Storing the charge here improves cache coherency in the charge spreading kernel
        }
    }
}

/**
 * For each grid point, find the range of sorted atoms associated with that point.
 */
extern "C" __global__ void findAtomRangeForGrid(int2* __restrict__ pmeAtomGridIndex, int* __restrict__ pmeAtomRange, const real4* __restrict__ posq, real4 periodicBoxSize, real4 invPeriodicBoxSize) {
    int start = (NUM_ATOMS*(blockIdx.x*blockDim.x+threadIdx.x))/(blockDim.x*gridDim.x);
    int end = (NUM_ATOMS*(blockIdx.x*blockDim.x+threadIdx.x+1))/(blockDim.x*gridDim.x);
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
    
    if (blockIdx.x == gridDim.x-1 && threadIdx.x == blockDim.x-1) {
        int gridSize = GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z;
        for (int j = last+1; j <= gridSize; ++j)
            pmeAtomRange[j] = NUM_ATOMS;
    }
}

#define BUFFER_SIZE (PME_ORDER*PME_ORDER*PME_ORDER)
extern "C" __global__ void gridSpreadCharge(const real4* __restrict__ posq, real* __restrict__ originalPmeGrid,
        const real4* __restrict__ pmeBsplineTheta, real4 periodicBoxSize, real4 invPeriodicBoxSize) {
    int ix = threadIdx.x/(PME_ORDER*PME_ORDER);
    int remainder = threadIdx.x-ix*PME_ORDER*PME_ORDER;
    int iy = remainder/PME_ORDER;
    int iz = remainder-iy*PME_ORDER;
    __shared__ real4 theta[PME_ORDER];
    __shared__ real charge[BUFFER_SIZE];
    __shared__ int basex[BUFFER_SIZE];
    __shared__ int basey[BUFFER_SIZE];
    __shared__ int basez[BUFFER_SIZE];
    if (ix < PME_ORDER) {
        for (int baseIndex = blockIdx.x*BUFFER_SIZE; baseIndex < NUM_ATOMS; baseIndex += gridDim.x*BUFFER_SIZE) {
            // Load the next block of atoms into the buffers.

            int atomIndex = baseIndex+threadIdx.x;
            if (atomIndex < NUM_ATOMS) {
                real4 pos = posq[atomIndex];
                charge[threadIdx.x] = pos.w;
                pos.x -= floor(pos.x*invPeriodicBoxSize.x)*periodicBoxSize.x;
                pos.y -= floor(pos.y*invPeriodicBoxSize.y)*periodicBoxSize.y;
                pos.z -= floor(pos.z*invPeriodicBoxSize.z)*periodicBoxSize.z;
                basex[threadIdx.x] = (int) ((pos.x*invPeriodicBoxSize.x)*GRID_SIZE_X);
                basey[threadIdx.x] = (int) ((pos.y*invPeriodicBoxSize.y)*GRID_SIZE_Y);
                basez[threadIdx.x] = (int) ((pos.z*invPeriodicBoxSize.z)*GRID_SIZE_Z);
            }
            __syncthreads();
            int lastIndex = min(BUFFER_SIZE, NUM_ATOMS-baseIndex);
            for (int index = 0; index < lastIndex; index++) {
                int atomIndex = index+baseIndex;
                if (threadIdx.x < PME_ORDER)
                    theta[threadIdx.x] = pmeBsplineTheta[atomIndex+threadIdx.x*NUM_ATOMS];
                __syncthreads();
                real add = charge[index]*theta[ix].x*theta[iy].y*theta[iz].z;
                int x = basex[index]+ix;
                int y = basey[index]+iy;
                int z = basez[index]+iz;
                x -= (x >= GRID_SIZE_X ? GRID_SIZE_X : 0);
                y -= (y >= GRID_SIZE_Y ? GRID_SIZE_Y : 0);
                z -= (z >= GRID_SIZE_Z ? GRID_SIZE_Z : 0);
#ifdef USE_DOUBLE_PRECISION
                unsigned long long * ulonglong_p = (unsigned long long *) originalPmeGrid;
                atomicAdd(&ulonglong_p[x*GRID_SIZE_Y*GRID_SIZE_Z+y*GRID_SIZE_Z+z],  static_cast<unsigned long long>((long long) (add*0x100000000)));
#elif __CUDA_ARCH__ < 200
                unsigned long long * ulonglong_p = (unsigned long long *) originalPmeGrid;
                int gridIndex = x*GRID_SIZE_Y*GRID_SIZE_Z+y*GRID_SIZE_Z+z;
                gridIndex = (gridIndex%2 == 0 ? gridIndex/2 : (gridIndex+GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z)/2);
                atomicAdd(&ulonglong_p[gridIndex],  static_cast<unsigned long long>((long long) (add*0x100000000)));
#else
                atomicAdd(&originalPmeGrid[x*GRID_SIZE_Y*GRID_SIZE_Z+y*GRID_SIZE_Z+z], add*EPSILON_FACTOR);
#endif
            }
        }
    }
}

extern "C" __global__ void finishSpreadCharge(long long* __restrict__ originalPmeGrid) {
    real* floatGrid = (real*) originalPmeGrid;
    const unsigned int gridSize = GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z;
    real scale = EPSILON_FACTOR/(real) 0x100000000;
#ifdef USE_DOUBLE_PRECISION
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < gridSize; index += blockDim.x*gridDim.x)
        floatGrid[index] = scale*originalPmeGrid[index];
#else
    for (int index = 2*(blockIdx.x*blockDim.x+threadIdx.x); index < gridSize; index += 2*blockDim.x*gridDim.x) {
        floatGrid[index] = scale*originalPmeGrid[index/2];
        if (index+1 < gridSize)
            floatGrid[index+1] = scale*originalPmeGrid[(index+gridSize+1)/2];
    }
#endif
}

// convolutes on the halfcomplex_pmeGrid, which is of size NX*NY*(NZ/2+1) as F(Q) is conjugate symmetric
extern "C" __global__ void 
reciprocalConvolution(real2* __restrict__ halfcomplex_pmeGrid, real* __restrict__ energyBuffer, 
                      const real* __restrict__ pmeBsplineModuliX,
                      const real* __restrict__ pmeBsplineModuliY, const real* __restrict__ pmeBsplineModuliZ, 
                      real4 periodicBoxSize, real4 invPeriodicBoxSize) {
    // R2C stores into a half complex matrix where the last dimension is cut by half
    const unsigned int gridSize = GRID_SIZE_X*GRID_SIZE_Y*(GRID_SIZE_Z/2+1);
    const real recipScaleFactor = RECIP(M_PI*periodicBoxSize.x*periodicBoxSize.y*periodicBoxSize.z);

    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < gridSize; index += blockDim.x*gridDim.x) {
        // real indices
        int kx = index/(GRID_SIZE_Y*(GRID_SIZE_Z/2+1));
        int remainder = index-kx*GRID_SIZE_Y*(GRID_SIZE_Z/2+1);
        int ky = remainder/(GRID_SIZE_Z/2+1);
        int kz = remainder-ky*(GRID_SIZE_Z/2+1);

        // reciprocal space indices
        int mx = (kx < (GRID_SIZE_X+1)/2) ? kx : (kx-GRID_SIZE_X);
        int my = (ky < (GRID_SIZE_Y+1)/2) ? ky : (ky-GRID_SIZE_Y);
        int mz = (kz < (GRID_SIZE_Z+1)/2) ? kz : (kz-GRID_SIZE_Z);

        // find the coordinates of the reciprocal space vectors
        real mhx = mx*invPeriodicBoxSize.x;
        real mhy = my*invPeriodicBoxSize.y;
        real mhz = mz*invPeriodicBoxSize.z;

        real bx = pmeBsplineModuliX[kx];
        real by = pmeBsplineModuliY[ky];
        real bz = pmeBsplineModuliZ[kz];

        real2 grid = halfcomplex_pmeGrid[index];

        real m2 = mhx*mhx+mhy*mhy+mhz*mhz;
        real denom = m2*bx*by*bz;

        real eterm = recipScaleFactor*EXP(-RECIP_EXP_FACTOR*m2)/denom;

        if (kx != 0 || ky != 0 || kz != 0) {
            halfcomplex_pmeGrid[index] = make_real2(grid.x*eterm, grid.y*eterm);
        }
    }
}

extern "C" __global__
void gridEvaluateEnergy(const real* __restrict__ originalGrid, const real* __restrict__ convolvedGrid, real* __restrict__ energyBuffer) {
    const unsigned int gridSize = GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z;
    real energy = 0;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < gridSize; index += blockDim.x*gridDim.x)
        energy += originalGrid[index]*convolvedGrid[index];
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += 0.5*energy;
}

extern "C" __global__
void gridInterpolateForce(const real4* __restrict__ posq, unsigned long long* __restrict__ forceBuffers, const real* __restrict__ originalPmeGrid,
        real4 periodicBoxSize, real4 invPeriodicBoxSize) {
    real3 data[PME_ORDER];
    real3 ddata[PME_ORDER];
    const real scale = RECIP(PME_ORDER-1);
     
    for (int atom = blockIdx.x*blockDim.x+threadIdx.x; atom < NUM_ATOMS; atom += blockDim.x*gridDim.x) {
        real3 force = make_real3(0);
        real4 pos = posq[atom];
        pos.x -= floor(pos.x*invPeriodicBoxSize.x)*periodicBoxSize.x;
        pos.y -= floor(pos.y*invPeriodicBoxSize.y)*periodicBoxSize.y;
        pos.z -= floor(pos.z*invPeriodicBoxSize.z)*periodicBoxSize.z;
        real3 t = make_real3((pos.x*invPeriodicBoxSize.x)*GRID_SIZE_X,
                             (pos.y*invPeriodicBoxSize.y)*GRID_SIZE_Y,
                             (pos.z*invPeriodicBoxSize.z)*GRID_SIZE_Z);
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
            real div = RECIP(j-1);
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
                    real gridvalue = originalPmeGrid[index];
                    force.x += ddx*dy*data[iz].z*gridvalue;
                    force.y += dx*ddy*data[iz].z*gridvalue;
                    force.z += dx*dy*ddata[iz].z*gridvalue;
                }
            }
        }
        real q = pos.w*EPSILON_FACTOR;
        forceBuffers[atom] +=  static_cast<unsigned long long>((long long) (-q*force.x*GRID_SIZE_X*invPeriodicBoxSize.x*0x100000000));
        forceBuffers[atom+PADDED_NUM_ATOMS] +=  static_cast<unsigned long long>((long long) (-q*force.y*GRID_SIZE_Y*invPeriodicBoxSize.y*0x100000000));
        forceBuffers[atom+2*PADDED_NUM_ATOMS] +=  static_cast<unsigned long long>((long long) (-q*force.z*GRID_SIZE_Z*invPeriodicBoxSize.z*0x100000000));
    }
}
