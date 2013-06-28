extern "C" __global__ void findAtomGridIndex(const real4* __restrict__ posq, int2* __restrict__ pmeAtomGridIndex,
            real4 periodicBoxSize, real4 invPeriodicBoxSize) {
    // Compute the index of the grid point each atom is associated with.
    
    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < NUM_ATOMS; i += blockDim.x*gridDim.x) {
        real4 pos = posq[i];
        pos.x -= floor(pos.x*invPeriodicBoxSize.x)*periodicBoxSize.x;
        pos.y -= floor(pos.y*invPeriodicBoxSize.y)*periodicBoxSize.y;
        pos.z -= floor(pos.z*invPeriodicBoxSize.z)*periodicBoxSize.z;
        real3 t = make_real3((pos.x*invPeriodicBoxSize.x)*GRID_SIZE_X,
                             (pos.y*invPeriodicBoxSize.y)*GRID_SIZE_Y,
                             (pos.z*invPeriodicBoxSize.z)*GRID_SIZE_Z);
        int3 gridIndex = make_int3(((int) t.x) % GRID_SIZE_X,
                                 ((int) t.y) % GRID_SIZE_Y,
                                 ((int) t.z) % GRID_SIZE_Z);
        pmeAtomGridIndex[i] = make_int2(i, gridIndex.x*GRID_SIZE_Y*GRID_SIZE_Z+gridIndex.y*GRID_SIZE_Z+gridIndex.z);
    }
}

extern "C" __global__ void gridSpreadCharge(const real4* __restrict__ posq, real* __restrict__ originalPmeGrid,
        real4 periodicBoxSize, real4 invPeriodicBoxSize, const int2* __restrict__ pmeAtomGridIndex) {
    real3 data[PME_ORDER];
    const real scale = RECIP(PME_ORDER-1);
    
    // Process the atoms in spatially sorted order.  This improves efficiency when writing
    // the grid values.
    
    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < NUM_ATOMS; i += blockDim.x*gridDim.x) {
        int atom = pmeAtomGridIndex[i].x;
        real charge = posq[atom].w;
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
        data[PME_ORDER-1] = scale*dr*data[PME_ORDER-2];
        for (int j = 1; j < (PME_ORDER-1); j++)
            data[PME_ORDER-j-1] = scale*((dr+make_real3(j))*data[PME_ORDER-j-2] + (make_real3(PME_ORDER-j)-dr)*data[PME_ORDER-j-1]);
        data[0] = scale*(make_real3(1)-dr)*data[0];
        
        // Spread the charge from this atom onto each grid point.
         
        for (int ix = 0; ix < PME_ORDER; ix++) {
            int xbase = gridIndex.x+ix;
            xbase -= (xbase >= GRID_SIZE_X ? GRID_SIZE_X : 0);
            xbase = xbase*GRID_SIZE_Y*GRID_SIZE_Z;
            real dx = data[ix].x;
            
            for (int iy = 0; iy < PME_ORDER; iy++) {
                int ybase = gridIndex.y+iy;
                ybase -= (ybase >= GRID_SIZE_Y ? GRID_SIZE_Y : 0);
                ybase = xbase + ybase*GRID_SIZE_Z;
                real dy = data[iy].y;
                
                for (int iz = 0; iz < PME_ORDER; iz++) {
                    int zindex = gridIndex.z+iz;
                    zindex -= (zindex >= GRID_SIZE_Z ? GRID_SIZE_Z : 0);
                    int index = ybase + zindex;

                    real add = charge*dx*dy*data[iz].z;
#ifdef USE_DOUBLE_PRECISION
                    unsigned long long * ulonglong_p = (unsigned long long *) originalPmeGrid;
                    atomicAdd(&ulonglong_p[index],  static_cast<unsigned long long>((long long) (add*0x100000000)));
#elif __CUDA_ARCH__ < 200
                    unsigned long long * ulonglong_p = (unsigned long long *) originalPmeGrid;
                    int gridIndex = index;
                    gridIndex = (gridIndex%2 == 0 ? gridIndex/2 : (gridIndex+GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z)/2);
                    atomicAdd(&ulonglong_p[gridIndex],  static_cast<unsigned long long>((long long) (add*0x100000000)));
#else
                    atomicAdd(&originalPmeGrid[index], add*EPSILON_FACTOR);
#endif
                }
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
        int mx = (kx < (GRID_SIZE_X+1)/2) ? kx : (kx-GRID_SIZE_X);
        int my = (ky < (GRID_SIZE_Y+1)/2) ? ky : (ky-GRID_SIZE_Y);
        int mz = (kz < (GRID_SIZE_Z+1)/2) ? kz : (kz-GRID_SIZE_Z);
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


extern "C" __global__ void
gridEvaluateEnergy(real2* __restrict__ halfcomplex_pmeGrid, real* __restrict__ energyBuffer,
                      const real* __restrict__ pmeBsplineModuliX,
                      const real* __restrict__ pmeBsplineModuliY, const real* __restrict__ pmeBsplineModuliZ,
                      real4 periodicBoxSize, real4 invPeriodicBoxSize) {
    // R2C stores into a half complex matrix where the last dimension is cut by half
    const unsigned int gridSize = GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z;
    const real recipScaleFactor = RECIP(M_PI*periodicBoxSize.x*periodicBoxSize.y*periodicBoxSize.z);
 
    real energy = 0;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < gridSize; index += blockDim.x*gridDim.x) {
        // real indices
        int kx = index/(GRID_SIZE_Y*(GRID_SIZE_Z));
        int remainder = index-kx*GRID_SIZE_Y*(GRID_SIZE_Z);
        int ky = remainder/(GRID_SIZE_Z);
        int kz = remainder-ky*(GRID_SIZE_Z);
        int mx = (kx < (GRID_SIZE_X+1)/2) ? kx : (kx-GRID_SIZE_X);
        int my = (ky < (GRID_SIZE_Y+1)/2) ? ky : (ky-GRID_SIZE_Y);
        int mz = (kz < (GRID_SIZE_Z+1)/2) ? kz : (kz-GRID_SIZE_Z);
        real mhx = mx*invPeriodicBoxSize.x;
        real mhy = my*invPeriodicBoxSize.y;
        real mhz = mz*invPeriodicBoxSize.z;
        real m2 = mhx*mhx+mhy*mhy+mhz*mhz;
        real bx = pmeBsplineModuliX[kx];
        real by = pmeBsplineModuliY[ky];
        real bz = pmeBsplineModuliZ[kz];
        real denom = m2*bx*by*bz;
        real eterm = recipScaleFactor*EXP(-RECIP_EXP_FACTOR*m2)/denom;

        if (kz >= (GRID_SIZE_Z/2+1)) {
            kx = ((kx == 0) ? kx : GRID_SIZE_X-kx);
            ky = ((ky == 0) ? ky : GRID_SIZE_Y-ky);
            kz = GRID_SIZE_Z-kz;
        } 
        int indexInHalfComplexGrid = kz + ky*(GRID_SIZE_Z/2+1)+kx*(GRID_SIZE_Y*(GRID_SIZE_Z/2+1));
        real2 grid = halfcomplex_pmeGrid[indexInHalfComplexGrid];
        if (kx != 0 || ky != 0 || kz != 0) {
            energy += eterm*(grid.x*grid.x + grid.y*grid.y);
        }
    }
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += 0.5f*energy;
}

extern "C" __global__
void gridInterpolateForce(const real4* __restrict__ posq, unsigned long long* __restrict__ forceBuffers, const real* __restrict__ originalPmeGrid,
        real4 periodicBoxSize, real4 invPeriodicBoxSize, const int2* __restrict__ pmeAtomGridIndex) {
    real3 data[PME_ORDER];
    real3 ddata[PME_ORDER];
    const real scale = RECIP(PME_ORDER-1);
    
    // Process the atoms in spatially sorted order.  This improves cache performance when loading
    // the grid values.
    
    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < NUM_ATOMS; i += blockDim.x*gridDim.x) {
        int atom = pmeAtomGridIndex[i].x;
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
        forceBuffers[atom] += static_cast<unsigned long long>((long long) (-q*force.x*GRID_SIZE_X*invPeriodicBoxSize.x*0x100000000));
        forceBuffers[atom+PADDED_NUM_ATOMS] += static_cast<unsigned long long>((long long) (-q*force.y*GRID_SIZE_Y*invPeriodicBoxSize.y*0x100000000));
        forceBuffers[atom+2*PADDED_NUM_ATOMS] += static_cast<unsigned long long>((long long) (-q*force.z*GRID_SIZE_Z*invPeriodicBoxSize.z*0x100000000));
    }
}

extern "C" __global__
void addForces(const real4* __restrict__ forces, unsigned long long* __restrict__ forceBuffers) {
    for (int atom = blockIdx.x*blockDim.x+threadIdx.x; atom < NUM_ATOMS; atom += blockDim.x*gridDim.x) {
        real4 f = forces[atom];
        forceBuffers[atom] += static_cast<unsigned long long>((long long) (f.x*0x100000000));
        forceBuffers[atom+PADDED_NUM_ATOMS] += static_cast<unsigned long long>((long long) (f.y*0x100000000));
        forceBuffers[atom+2*PADDED_NUM_ATOMS] += static_cast<unsigned long long>((long long) (f.z*0x100000000));
    }
}
