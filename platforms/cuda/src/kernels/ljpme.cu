extern "C" __global__ void findAtomDispersionGridIndex(const real4* __restrict__ posq, int2* __restrict__ pmeAtomGridIndex,
            real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
            real3 recipBoxVecX, real3 recipBoxVecY, real3 recipBoxVecZ) {
    // Compute the index of the grid point each atom is associated with.
    
    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < NUM_ATOMS; i += blockDim.x*gridDim.x) {
        real4 pos = posq[i];
        APPLY_PERIODIC_TO_POS(pos)
        real3 t = make_real3(pos.x*recipBoxVecX.x+pos.y*recipBoxVecY.x+pos.z*recipBoxVecZ.x,
                             pos.y*recipBoxVecY.y+pos.z*recipBoxVecZ.y,
                             pos.z*recipBoxVecZ.z);
        t.x = (t.x-floor(t.x))*DISPERSION_GRID_SIZE_X;
        t.y = (t.y-floor(t.y))*DISPERSION_GRID_SIZE_Y;
        t.z = (t.z-floor(t.z))*DISPERSION_GRID_SIZE_Z;
        int3 gridIndex = make_int3(((int) t.x) % DISPERSION_GRID_SIZE_X,
                                   ((int) t.y) % DISPERSION_GRID_SIZE_Y,
                                   ((int) t.z) % DISPERSION_GRID_SIZE_Z);
        pmeAtomGridIndex[i] = make_int2(i, gridIndex.x*DISPERSION_GRID_SIZE_Y*DISPERSION_GRID_SIZE_Z+gridIndex.y*DISPERSION_GRID_SIZE_Z+gridIndex.z);
    }
}

extern "C" __global__ void gridSpreadC6(const real4* __restrict__ posq, real* __restrict__ originalPmeGrid,
        real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        real3 recipBoxVecX, real3 recipBoxVecY, real3 recipBoxVecZ, const int2* __restrict__ pmeAtomGridIndex,
        const real* __restrict__ C6s) {
    real3 data[PME_ORDER];
    const real scale = RECIP(PME_ORDER-1);
    
    // Process the atoms in spatially sorted order.  This improves efficiency when writing
    // the grid values.
    
    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < NUM_ATOMS; i += blockDim.x*gridDim.x) {
        int atom = pmeAtomGridIndex[i].x;
        real4 pos = posq[atom];
        APPLY_PERIODIC_TO_POS(pos)
        real3 t = make_real3(pos.x*recipBoxVecX.x+pos.y*recipBoxVecY.x+pos.z*recipBoxVecZ.x,
                             pos.y*recipBoxVecY.y+pos.z*recipBoxVecZ.y,
                             pos.z*recipBoxVecZ.z);
        t.x = (t.x-floor(t.x))*DISPERSION_GRID_SIZE_X;
        t.y = (t.y-floor(t.y))*DISPERSION_GRID_SIZE_Y;
        t.z = (t.z-floor(t.z))*DISPERSION_GRID_SIZE_Z;
        int3 gridIndex = make_int3(((int) t.x) % DISPERSION_GRID_SIZE_X,
                                   ((int) t.y) % DISPERSION_GRID_SIZE_Y,
                                   ((int) t.z) % DISPERSION_GRID_SIZE_Z);

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
            xbase -= (xbase >= DISPERSION_GRID_SIZE_X ? DISPERSION_GRID_SIZE_X : 0);
            xbase = xbase*DISPERSION_GRID_SIZE_Y*DISPERSION_GRID_SIZE_Z;
            real dx = data[ix].x;
            
            for (int iy = 0; iy < PME_ORDER; iy++) {
                int ybase = gridIndex.y+iy;
                ybase -= (ybase >= DISPERSION_GRID_SIZE_Y ? DISPERSION_GRID_SIZE_Y : 0);
                ybase = xbase + ybase*DISPERSION_GRID_SIZE_Z;
                real dy = data[iy].y;
                
                for (int iz = 0; iz < PME_ORDER; iz++) {
                    int zindex = gridIndex.z+iz;
                    zindex -= (zindex >= DISPERSION_GRID_SIZE_Z ? DISPERSION_GRID_SIZE_Z : 0);
                    int index = ybase + zindex;

                    // We need to grab the C6 coefficient from the array
                    real add = C6s[atom]*dx*dy*data[iz].z;
#ifdef USE_DOUBLE_PRECISION
                    unsigned long long * ulonglong_p = (unsigned long long *) originalPmeGrid;
                    atomicAdd(&ulonglong_p[index],  static_cast<unsigned long long>((long long) (add*0x100000000)));
#elif __CUDA_ARCH__ < 200 || defined(USE_DETERMINISTIC_FORCES)
                    unsigned long long * ulonglong_p = (unsigned long long *) originalPmeGrid;
                    int gridIndex = index;
                    gridIndex = (gridIndex%2 == 0 ? gridIndex/2 : (gridIndex+DISPERSION_GRID_SIZE_X*DISPERSION_GRID_SIZE_Y*DISPERSION_GRID_SIZE_Z)/2);
                    atomicAdd(&ulonglong_p[gridIndex],  static_cast<unsigned long long>((long long) (add*0x100000000)));
#else
                    atomicAdd(&originalPmeGrid[index], add);
#endif

                }
            }
        }
    }
}


extern "C" __global__ void finishSpreadC6(long long* __restrict__ originalPmeGrid) {
    real* floatGrid = (real*) originalPmeGrid;
    const unsigned int gridSize = DISPERSION_GRID_SIZE_X*DISPERSION_GRID_SIZE_Y*DISPERSION_GRID_SIZE_Z;
    real scale = 1.0f/(real) 0x100000000;
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


// convolutes the dispersion grid on the halfcomplex_pmeGrid, which is of size NX*NY*(NZ/2+1) as F(Q) is conjugate symmetric
extern "C" __global__ void 
reciprocalDispersionConvolution(real2* __restrict__ halfcomplex_pmeGrid, mixed* __restrict__ energyBuffer, 
                      const real* __restrict__ pmeBsplineModuliX, const real* __restrict__ pmeBsplineModuliY, const real* __restrict__ pmeBsplineModuliZ, 
                      real4 periodicBoxSize, real3 recipBoxVecX, real3 recipBoxVecY, real3 recipBoxVecZ) {
    // R2C stores into a half complex matrix where the last dimension is cut by half
    const unsigned int gridSize = GRID_SIZE_X*GRID_SIZE_Y*(GRID_SIZE_Z/2+1);
    const real scaleFactor =  -2*M_PI*SQRT(M_PI)*RECIP(6*periodicBoxSize.x*periodicBoxSize.y*periodicBoxSize.z);

    const real alpha = EWALD_DISPERSION_ALPHA;
    real bfac = M_PI / alpha;
    real fac1 = 2*M_PI*M_PI*M_PI*SQRT(M_PI);
    real fac2 = alpha*alpha*alpha;
    real fac3 = -2*alpha*M_PI*M_PI;

    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < gridSize; index += blockDim.x*gridDim.x) {
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
        real denom = scaleFactor/(bx*by*bz);
        real2 grid = halfcomplex_pmeGrid[index];
        real m2 = mhx*mhx+mhy*mhy+mhz*mhz;
        real m = SQRT(m2);
        real m3 = m*m2;
        real b = bfac*m;
        real expfac = -b*b;
        real expterm = EXP(expfac);
#if FAST_ERFC
        // This approximation for erfc is from Abramowitz and Stegun (1964) p. 299.  They cite the following as
        // the original source: C. Hastings, Jr., Approximations for Digital Computers (1955).  It has a maximum
        // error of 1.5e-7.  Stolen by ACS from the CUDA platform's AMOEBA plugin.
        real t = 1.0f/(1.0f+0.3275911f*b);
        real erfcterm = (0.254829592f+(-0.284496736f+(1.421413741f+(-1.453152027f+1.061405429f*t)*t)*t)*t)*t*expterm;
#else
        real erfcterm = ERFC(b);
#endif
        real eterm = (fac1*erfcterm*m3 + expterm*(fac2 + fac3*m2)) * denom;
        halfcomplex_pmeGrid[index] = make_real2(grid.x*eterm, grid.y*eterm);
    }
}


extern "C" __global__ void
gridEvaluateDispersionEnergy(real2* __restrict__ halfcomplex_pmeGrid, mixed* __restrict__ energyBuffer,
                      const real* __restrict__ pmeBsplineModuliX, const real* __restrict__ pmeBsplineModuliY, const real* __restrict__ pmeBsplineModuliZ,
                      real4 periodicBoxSize, real3 recipBoxVecX, real3 recipBoxVecY, real3 recipBoxVecZ) {
    // R2C stores into a half complex matrix where the last dimension is cut by half
    const unsigned int gridSize = DISPERSION_GRID_SIZE_X*DISPERSION_GRID_SIZE_Y*DISPERSION_GRID_SIZE_Z;
    const real scaleFactor =  -2*M_PI*SQRT(M_PI)*RECIP(6*periodicBoxSize.x*periodicBoxSize.y*periodicBoxSize.z);

    const real alpha = EWALD_DISPERSION_ALPHA;
    real bfac = M_PI / alpha;
    real fac1 = 2*M_PI*M_PI*M_PI*SQRT(M_PI);
    real fac2 = alpha*alpha*alpha;
    real fac3 = -2*alpha*M_PI*M_PI;

    mixed energy = 0;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < gridSize; index += blockDim.x*gridDim.x) {
        // real indices
        int kx = index/(DISPERSION_GRID_SIZE_Y*(DISPERSION_GRID_SIZE_Z));
        int remainder = index-kx*DISPERSION_GRID_SIZE_Y*(DISPERSION_GRID_SIZE_Z);
        int ky = remainder/(DISPERSION_GRID_SIZE_Z);
        int kz = remainder-ky*(DISPERSION_GRID_SIZE_Z);
        int mx = (kx < (DISPERSION_GRID_SIZE_X+1)/2) ? kx : (kx-DISPERSION_GRID_SIZE_X);
        int my = (ky < (DISPERSION_GRID_SIZE_Y+1)/2) ? ky : (ky-DISPERSION_GRID_SIZE_Y);
        int mz = (kz < (DISPERSION_GRID_SIZE_Z+1)/2) ? kz : (kz-DISPERSION_GRID_SIZE_Z);
        real mhx = mx*recipBoxVecX.x;
        real mhy = mx*recipBoxVecY.x+my*recipBoxVecY.y;
        real mhz = mx*recipBoxVecZ.x+my*recipBoxVecZ.y+mz*recipBoxVecZ.z;
        real m2 = mhx*mhx+mhy*mhy+mhz*mhz;
        real bx = pmeBsplineModuliX[kx];
        real by = pmeBsplineModuliY[ky];
        real bz = pmeBsplineModuliZ[kz];
        real denom = scaleFactor/(bx*by*bz);
        real m = SQRT(m2);
        real m3 = m*m2;
        real b = bfac*m;
        real expfac = -b*b;
        real expterm = EXP(expfac);
#if FAST_ERFC
        // This approximation for erfc is from Abramowitz and Stegun (1964) p. 299.  They cite the following as
        // the original source: C. Hastings, Jr., Approximations for Digital Computers (1955).  It has a maximum
        // error of 1.5e-7.  Stolen by ACS from the CUDA platform's AMOEBA plugin.
        real t = 1.0f/(1.0f+0.3275911f*b);
        real erfcterm = (0.254829592f+(-0.284496736f+(1.421413741f+(-1.453152027f+1.061405429f*t)*t)*t)*t)*t*expterm;
#else
        real erfcterm = ERFC(b);
#endif
        real eterm = (fac1*erfcterm*m3 + expterm*(fac2 + fac3*m2)) * denom;

        if (kz >= (DISPERSION_GRID_SIZE_Z/2+1)) {
            kx = ((kx == 0) ? kx : DISPERSION_GRID_SIZE_X-kx);
            ky = ((ky == 0) ? ky : DISPERSION_GRID_SIZE_Y-ky);
            kz = DISPERSION_GRID_SIZE_Z-kz;
        } 
        int indexInHalfComplexGrid = kz + ky*(DISPERSION_GRID_SIZE_Z/2+1)+kx*(DISPERSION_GRID_SIZE_Y*(DISPERSION_GRID_SIZE_Z/2+1));
        real2 grid = halfcomplex_pmeGrid[indexInHalfComplexGrid];
        // N.B. We inlcude the 0,0,0 point for dispersion
        energy += eterm*(grid.x*grid.x + grid.y*grid.y);
    }
#ifdef USE_PME_STREAM
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] = 0.5f*energy;
#else
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += 0.5f*energy;
#endif
}


extern "C" __global__
void gridInterpolateDispersionForce(const real4* __restrict__ posq, unsigned long long* __restrict__ forceBuffers, const real* __restrict__ originalPmeGrid,
        real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        real3 recipBoxVecX, real3 recipBoxVecY, real3 recipBoxVecZ, const int2* __restrict__ pmeAtomGridIndex, const real* __restrict__ C6s) {
    real3 data[PME_ORDER];
    real3 ddata[PME_ORDER];
    const real scale = RECIP(PME_ORDER-1);
    
    // Process the atoms in spatially sorted order.  This improves cache performance when loading
    // the grid values.
    
    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < NUM_ATOMS; i += blockDim.x*gridDim.x) {
        int atom = pmeAtomGridIndex[i].x;
        real3 force = make_real3(0);
        real4 pos = posq[atom];
        APPLY_PERIODIC_TO_POS(pos)
        real3 t = make_real3(pos.x*recipBoxVecX.x+pos.y*recipBoxVecY.x+pos.z*recipBoxVecZ.x,
                             pos.y*recipBoxVecY.y+pos.z*recipBoxVecZ.y,
                             pos.z*recipBoxVecZ.z);
        t.x = (t.x-floor(t.x))*DISPERSION_GRID_SIZE_X;
        t.y = (t.y-floor(t.y))*DISPERSION_GRID_SIZE_Y;
        t.z = (t.z-floor(t.z))*DISPERSION_GRID_SIZE_Z;
        int3 gridIndex = make_int3(((int) t.x) % DISPERSION_GRID_SIZE_X,
                                   ((int) t.y) % DISPERSION_GRID_SIZE_Y,
                                   ((int) t.z) % DISPERSION_GRID_SIZE_Z);
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
            xbase -= (xbase >= DISPERSION_GRID_SIZE_X ? DISPERSION_GRID_SIZE_X : 0);
            xbase = xbase*DISPERSION_GRID_SIZE_Y*DISPERSION_GRID_SIZE_Z;
            real dx = data[ix].x;
            real ddx = ddata[ix].x;
            
            for (int iy = 0; iy < PME_ORDER; iy++) {
                int ybase = gridIndex.y+iy;
                ybase -= (ybase >= DISPERSION_GRID_SIZE_Y ? DISPERSION_GRID_SIZE_Y : 0);
                ybase = xbase + ybase*DISPERSION_GRID_SIZE_Z;
                real dy = data[iy].y;
                real ddy = ddata[iy].y;
                
                for (int iz = 0; iz < PME_ORDER; iz++) {
                    int zindex = gridIndex.z+iz;
                    zindex -= (zindex >= DISPERSION_GRID_SIZE_Z ? DISPERSION_GRID_SIZE_Z : 0);
                    int index = ybase + zindex;
                    real gridvalue = originalPmeGrid[index];
                    force.x += ddx*dy*data[iz].z*gridvalue;
                    force.y += dx*ddy*data[iz].z*gridvalue;
                    force.z += dx*dy*ddata[iz].z*gridvalue;
                }
            }
        }
        real q = C6s[atom];
        real forceX = -q*(force.x*DISPERSION_GRID_SIZE_X*recipBoxVecX.x);
        real forceY = -q*(force.x*DISPERSION_GRID_SIZE_X*recipBoxVecY.x+force.y*DISPERSION_GRID_SIZE_Y*recipBoxVecY.y);
        real forceZ = -q*(force.x*DISPERSION_GRID_SIZE_X*recipBoxVecZ.x+force.y*DISPERSION_GRID_SIZE_Y*recipBoxVecZ.y+force.z*DISPERSION_GRID_SIZE_Z*recipBoxVecZ.z);
        atomicAdd(&forceBuffers[atom], static_cast<unsigned long long>((long long) (forceX*0x100000000)));
        atomicAdd(&forceBuffers[atom+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (forceY*0x100000000)));
        atomicAdd(&forceBuffers[atom+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (forceZ*0x100000000)));
    }
}

