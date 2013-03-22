#define ARRAY(x,y) array[(x)-1+((y)-1)*PME_ORDER]

/**
 * Calculate the spline coefficients for a single atom along a single axis.
 */
__device__ void computeBSplinePoint(real4* thetai, real w, real* array) {
    // initialization to get to 2nd order recursion

    ARRAY(2,2) = w;
    ARRAY(2,1) = 1 - w;

    // perform one pass to get to 3rd order recursion

    ARRAY(3,3) = 0.5f * w * ARRAY(2,2);
    ARRAY(3,2) = 0.5f * ((1+w)*ARRAY(2,1)+(2-w)*ARRAY(2,2));
    ARRAY(3,1) = 0.5f * (1-w) * ARRAY(2,1);

    // compute standard B-spline recursion to desired order

    for (int i = 4; i <= PME_ORDER; i++)
    {
        int k = i - 1;
        real denom = RECIP(k);
        ARRAY(i,i) = denom * w * ARRAY(k,k);
        for (int j = 1; j <= i-2; j++)
            ARRAY(i,i-j) = denom * ((w+j)*ARRAY(k,i-j-1)+(i-j-w)*ARRAY(k,i-j));
        ARRAY(i,1) = denom * (1-w) * ARRAY(k,1);
    }

    // get coefficients for the B-spline first derivative

    int k = PME_ORDER - 1;
    ARRAY(k,PME_ORDER) = ARRAY(k,PME_ORDER-1);
    for (int i = PME_ORDER-1; i >= 2; i--)
        ARRAY(k,i) = ARRAY(k,i-1) - ARRAY(k,i);
    ARRAY(k,1) = -ARRAY(k,1);

    // get coefficients for the B-spline second derivative

    k = PME_ORDER - 2;
    ARRAY(k,PME_ORDER-1) = ARRAY(k,PME_ORDER-2);
    for (int i = PME_ORDER-2; i >= 2; i--)
        ARRAY(k,i) = ARRAY(k,i-1) - ARRAY(k,i);
    ARRAY(k,1) = -ARRAY(k,1);
    ARRAY(k,PME_ORDER) = ARRAY(k,PME_ORDER-1);
    for (int i = PME_ORDER-1; i >= 2; i--)
        ARRAY(k,i) = ARRAY(k,i-1) - ARRAY(k,i);
    ARRAY(k,1) = -ARRAY(k,1);

    // get coefficients for the B-spline third derivative

    k = PME_ORDER - 3;
    ARRAY(k,PME_ORDER-2) = ARRAY(k,PME_ORDER-3);
    for (int i = PME_ORDER-3; i >= 2; i--)
        ARRAY(k,i) = ARRAY(k,i-1) - ARRAY(k,i);
    ARRAY(k,1) = -ARRAY(k,1);
    ARRAY(k,PME_ORDER-1) = ARRAY(k,PME_ORDER-2);
    for (int i = PME_ORDER-2; i >= 2; i--)
        ARRAY(k,i) = ARRAY(k,i-1) - ARRAY(k,i);
    ARRAY(k,1) = -ARRAY(k,1);
    ARRAY(k,PME_ORDER) = ARRAY(k,PME_ORDER-1);
    for (int i = PME_ORDER-1; i >= 2; i--)
        ARRAY(k,i) = ARRAY(k,i-1) - ARRAY(k,i);
    ARRAY(k,1) = -ARRAY(k,1);

    // copy coefficients from temporary to permanent storage

    for (int i = 1; i <= PME_ORDER; i++)
        thetai[i-1] = make_real4(ARRAY(PME_ORDER,i), ARRAY(PME_ORDER-1,i), ARRAY(PME_ORDER-2,i), ARRAY(PME_ORDER-3,i));
}

/**
 * Compute the index of the grid point each atom is associated with.
 */
extern "C" __global__ void findAtomGridIndex(const real4* __restrict__ posq, int2* __restrict__ pmeAtomGridIndex,
        real4 periodicBoxSize, real4 invPeriodicBoxSize) {
    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < NUM_ATOMS; i += blockDim.x*gridDim.x) {
        real4 pos = posq[i];
        pos.x -= floor(pos.x*invPeriodicBoxSize.x)*periodicBoxSize.x;
        pos.y -= floor(pos.y*invPeriodicBoxSize.y)*periodicBoxSize.y;
        pos.z -= floor(pos.z*invPeriodicBoxSize.z)*periodicBoxSize.z;

        // First axis.

        real w = pos.x*invPeriodicBoxSize.x;
        real fr = GRID_SIZE_X*(w-(int)(w+0.5f)+0.5f);
        int ifr = (int) fr;
        int igrid1 = ifr-PME_ORDER+1;

        // Second axis.

        w = pos.y*invPeriodicBoxSize.y;
        fr = GRID_SIZE_Y*(w-(int)(w+0.5f)+0.5f);
        ifr = (int) fr;
        int igrid2 = ifr-PME_ORDER+1;

        // Third axis.

        w = pos.z*invPeriodicBoxSize.z;
        fr = GRID_SIZE_Z*(w-(int)(w+0.5f)+0.5f);
        ifr = (int) fr;
        int igrid3 = ifr-PME_ORDER+1;

        // Record the grid point.

        igrid1 += (igrid1 < 0 ? GRID_SIZE_X : 0);
        igrid2 += (igrid2 < 0 ? GRID_SIZE_Y : 0);
        igrid3 += (igrid3 < 0 ? GRID_SIZE_Z : 0);
        pmeAtomGridIndex[i] = make_int2(i, igrid1*GRID_SIZE_Y*GRID_SIZE_Z+igrid2*GRID_SIZE_Z+igrid3);
    }
}

extern "C" __global__ void gridSpreadFixedMultipoles(const real4* __restrict__ posq, const real* __restrict__ labFrameDipole,
        const real* __restrict__ labFrameQuadrupole, real2* __restrict__ pmeGrid, int2* __restrict__ pmeAtomGridIndex,
        real4 periodicBoxSize, real4 invPeriodicBoxSize) {
    const real xscale = GRID_SIZE_X*invPeriodicBoxSize.x;
    const real yscale = GRID_SIZE_Y*invPeriodicBoxSize.y;
    const real zscale = GRID_SIZE_Z*invPeriodicBoxSize.z;
    real array[PME_ORDER*PME_ORDER];
    real4 theta1[PME_ORDER];
    real4 theta2[PME_ORDER];
    real4 theta3[PME_ORDER];
    
    // Process the atoms in spatially sorted order.  This improves cache performance when loading
    // the grid values.
    
    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < NUM_ATOMS; i += blockDim.x*gridDim.x) {
        int m = pmeAtomGridIndex[i].x;
        real4 pos = posq[m];
        pos.x -= floor(pos.x*invPeriodicBoxSize.x)*periodicBoxSize.x;
        pos.y -= floor(pos.y*invPeriodicBoxSize.y)*periodicBoxSize.y;
        pos.z -= floor(pos.z*invPeriodicBoxSize.z)*periodicBoxSize.z;

        // Since we need the full set of thetas, it's faster to compute them here than load them
        // from global memory.

        real w = pos.x*invPeriodicBoxSize.x;
        real fr = GRID_SIZE_X*(w-(int)(w+0.5f)+0.5f);
        int ifr = (int) fr;
        w = fr - ifr;
        int igrid1 = ifr-PME_ORDER+1;
        computeBSplinePoint(theta1, w, array);
        w = pos.y*invPeriodicBoxSize.y;
        fr = GRID_SIZE_Y*(w-(int)(w+0.5f)+0.5f);
        ifr = (int) fr;
        w = fr - ifr;
        int igrid2 = ifr-PME_ORDER+1;
        computeBSplinePoint(theta2, w, array);
        w = pos.z*invPeriodicBoxSize.z;
        fr = GRID_SIZE_Z*(w-(int)(w+0.5f)+0.5f);
        ifr = (int) fr;
        w = fr - ifr;
        int igrid3 = ifr-PME_ORDER+1;
        computeBSplinePoint(theta3, w, array);
        igrid1 += (igrid1 < 0 ? GRID_SIZE_X : 0);
        igrid2 += (igrid2 < 0 ? GRID_SIZE_Y : 0);
        igrid3 += (igrid3 < 0 ? GRID_SIZE_Z : 0);
        
        // Spread the charge from this atom onto each grid point.
         
        for (int ix = 0; ix < PME_ORDER; ix++) {
            int xbase = igrid1+ix;
            xbase -= (xbase >= GRID_SIZE_X ? GRID_SIZE_X : 0);
            xbase = xbase*GRID_SIZE_Y*GRID_SIZE_Z;
            real4 t = theta1[ix];
            
            for (int iy = 0; iy < PME_ORDER; iy++) {
                int ybase = igrid2+iy;
                ybase -= (ybase >= GRID_SIZE_Y ? GRID_SIZE_Y : 0);
                ybase = xbase + ybase*GRID_SIZE_Z;
                real4 u = theta2[iy];
                
                for (int iz = 0; iz < PME_ORDER; iz++) {
                    int zindex = igrid3+iz;
                    zindex -= (zindex >= GRID_SIZE_Z ? GRID_SIZE_Z : 0);
                    int index = ybase + zindex;
                    real4 v = theta3[iz];

                    real atomCharge = pos.w;
                    real atomDipoleX = xscale*labFrameDipole[m*3];
                    real atomDipoleY = yscale*labFrameDipole[m*3+1];
                    real atomDipoleZ = zscale*labFrameDipole[m*3+2];
                    real atomQuadrupoleXX = xscale*xscale*labFrameQuadrupole[m*5];
                    real atomQuadrupoleXY = 2*xscale*yscale*labFrameQuadrupole[m*5+1];
                    real atomQuadrupoleXZ = 2*xscale*zscale*labFrameQuadrupole[m*5+2];
                    real atomQuadrupoleYY = yscale*yscale*labFrameQuadrupole[m*5+3];
                    real atomQuadrupoleYZ = 2*yscale*zscale*labFrameQuadrupole[m*5+4];
                    real atomQuadrupoleZZ = -zscale*zscale*(labFrameQuadrupole[m*5]+labFrameQuadrupole[m*5+3]);
                    real term0 = atomCharge*u.x*v.x + atomDipoleY*u.y*v.x + atomDipoleZ*u.x*v.y + atomQuadrupoleYY*u.z*v.x + atomQuadrupoleZZ*u.x*v.z + atomQuadrupoleYZ*u.y*v.y;
                    real term1 = atomDipoleX*u.x*v.x + atomQuadrupoleXY*u.y*v.x + atomQuadrupoleXZ*u.x*v.y;
                    real term2 = atomQuadrupoleXX * u.x * v.x;
                    real add = term0*t.x + term1*t.y + term2*t.z;
#ifdef USE_DOUBLE_PRECISION
                    unsigned long long * ulonglong_p = (unsigned long long *) pmeGrid;
                    atomicAdd(&ulonglong_p[2*index],  static_cast<unsigned long long>((long long) (add*0x100000000)));
#else
                    atomicAdd(&pmeGrid[index].x, add);
#endif
                }
            }
        }
    }
}

extern "C" __global__ void gridSpreadInducedDipoles(const real4* __restrict__ posq, const real* __restrict__ inducedDipole,
        const real* __restrict__ inducedDipolePolar, real2* __restrict__ pmeGrid, int2* __restrict__ pmeAtomGridIndex,
        real4 periodicBoxSize, real4 invPeriodicBoxSize) {
    const real xscale = GRID_SIZE_X*invPeriodicBoxSize.x;
    const real yscale = GRID_SIZE_Y*invPeriodicBoxSize.y;
    const real zscale = GRID_SIZE_Z*invPeriodicBoxSize.z;
    real array[PME_ORDER*PME_ORDER];
    real4 theta1[PME_ORDER];
    real4 theta2[PME_ORDER];
    real4 theta3[PME_ORDER];
    
    // Process the atoms in spatially sorted order.  This improves cache performance when loading
    // the grid values.
    
    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < NUM_ATOMS; i += blockDim.x*gridDim.x) {
        int m = pmeAtomGridIndex[i].x;
        real4 pos = posq[m];
        pos.x -= floor(pos.x*invPeriodicBoxSize.x)*periodicBoxSize.x;
        pos.y -= floor(pos.y*invPeriodicBoxSize.y)*periodicBoxSize.y;
        pos.z -= floor(pos.z*invPeriodicBoxSize.z)*periodicBoxSize.z;

        // Since we need the full set of thetas, it's faster to compute them here than load them
        // from global memory.

        real w = pos.x*invPeriodicBoxSize.x;
        real fr = GRID_SIZE_X*(w-(int)(w+0.5f)+0.5f);
        int ifr = (int) fr;
        w = fr - ifr;
        int igrid1 = ifr-PME_ORDER+1;
        computeBSplinePoint(theta1, w, array);
        w = pos.y*invPeriodicBoxSize.y;
        fr = GRID_SIZE_Y*(w-(int)(w+0.5f)+0.5f);
        ifr = (int) fr;
        w = fr - ifr;
        int igrid2 = ifr-PME_ORDER+1;
        computeBSplinePoint(theta2, w, array);
        w = pos.z*invPeriodicBoxSize.z;
        fr = GRID_SIZE_Z*(w-(int)(w+0.5f)+0.5f);
        ifr = (int) fr;
        w = fr - ifr;
        int igrid3 = ifr-PME_ORDER+1;
        computeBSplinePoint(theta3, w, array);
        igrid1 += (igrid1 < 0 ? GRID_SIZE_X : 0);
        igrid2 += (igrid2 < 0 ? GRID_SIZE_Y : 0);
        igrid3 += (igrid3 < 0 ? GRID_SIZE_Z : 0);
        
        // Spread the charge from this atom onto each grid point.
         
        for (int ix = 0; ix < PME_ORDER; ix++) {
            int xbase = igrid1+ix;
            xbase -= (xbase >= GRID_SIZE_X ? GRID_SIZE_X : 0);
            xbase = xbase*GRID_SIZE_Y*GRID_SIZE_Z;
            real4 t = theta1[ix];
            
            for (int iy = 0; iy < PME_ORDER; iy++) {
                int ybase = igrid2+iy;
                ybase -= (ybase >= GRID_SIZE_Y ? GRID_SIZE_Y : 0);
                ybase = xbase + ybase*GRID_SIZE_Z;
                real4 u = theta2[iy];
                
                for (int iz = 0; iz < PME_ORDER; iz++) {
                    int zindex = igrid3+iz;
                    zindex -= (zindex >= GRID_SIZE_Z ? GRID_SIZE_Z : 0);
                    int index = ybase + zindex;
                    real4 v = theta3[iz];

                    real inducedDipoleX = xscale*inducedDipole[m*3];
                    real inducedDipoleY = yscale*inducedDipole[m*3+1];
                    real inducedDipoleZ = zscale*inducedDipole[m*3+2];
                    real inducedDipolePolarX = xscale*inducedDipolePolar[m*3];
                    real inducedDipolePolarY = yscale*inducedDipolePolar[m*3+1];
                    real inducedDipolePolarZ = zscale*inducedDipolePolar[m*3+2];
                    real term01 = inducedDipoleY*u.y*v.x + inducedDipoleZ*u.x*v.y;
                    real term11 = inducedDipoleX*u.x*v.x;
                    real term02 = inducedDipolePolarY*u.y*v.x + inducedDipolePolarZ*u.x*v.y;
                    real term12 = inducedDipolePolarX*u.x*v.x;
                    real add1 = term01*t.x + term11*t.y;
                    real add2 = term02*t.x + term12*t.y;
#ifdef USE_DOUBLE_PRECISION
                    unsigned long long * ulonglong_p = (unsigned long long *) pmeGrid;
                    atomicAdd(&ulonglong_p[2*index],  static_cast<unsigned long long>((long long) (add1*0x100000000)));
                    atomicAdd(&ulonglong_p[2*index+1],  static_cast<unsigned long long>((long long) (add2*0x100000000)));
#else
                    atomicAdd(&pmeGrid[index].x, add1);
                    atomicAdd(&pmeGrid[index].y, add2);
#endif
                }
            }
        }
    }
}

/**
 * In double precision, we have to use fixed point to accumulate the grid values, so convert them to floating point.
 */
extern "C" __global__ void finishSpreadCharge(long long* __restrict__ pmeGrid) {
    real* floatGrid = (real*) pmeGrid;
    const unsigned int gridSize = 2*GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z;
    real scale = 1/(real) 0x100000000;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < gridSize; index += blockDim.x*gridDim.x)
        floatGrid[index] = scale*pmeGrid[index];
}

extern "C" __global__ void reciprocalConvolution(real2* __restrict__ pmeGrid, const real* __restrict__ pmeBsplineModuliX,
        const real* __restrict__ pmeBsplineModuliY, const real* __restrict__ pmeBsplineModuliZ, real4 periodicBoxSize, real4 invPeriodicBoxSize) {
    const unsigned int gridSize = GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z;
    real expFactor = M_PI*M_PI/(EWALD_ALPHA*EWALD_ALPHA);
    real scaleFactor = RECIP(M_PI*periodicBoxSize.x*periodicBoxSize.y*periodicBoxSize.z);
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < gridSize; index += blockDim.x*gridDim.x) {
        int kx = index/(GRID_SIZE_Y*GRID_SIZE_Z);
        int remainder = index-kx*GRID_SIZE_Y*GRID_SIZE_Z;
        int ky = remainder/GRID_SIZE_Z;
        int kz = remainder-ky*GRID_SIZE_Z;
        if (kx == 0 && ky == 0 && kz == 0) {
            pmeGrid[index] = make_real2(0, 0);
            continue;
        }
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
        real eterm = scaleFactor*EXP(-expFactor*m2)/denom;
        pmeGrid[index] = make_real2(grid.x*eterm, grid.y*eterm);
    }
}

extern "C" __global__ void computeFixedPotentialFromGrid(const real2* __restrict__ pmeGrid, real* __restrict__ phi,
        long long* __restrict__ fieldBuffers, long long* __restrict__ fieldPolarBuffers,  const real4* __restrict__ posq,
        const real* __restrict__ labFrameDipole, real4 periodicBoxSize, real4 invPeriodicBoxSize, int2* __restrict__ pmeAtomGridIndex) {
    real array[PME_ORDER*PME_ORDER];
    real4 theta1[PME_ORDER];
    real4 theta2[PME_ORDER];
    real4 theta3[PME_ORDER];
    
    // Process the atoms in spatially sorted order.  This improves cache performance when loading
    // the grid values.
    
    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < NUM_ATOMS; i += blockDim.x*gridDim.x) {
        int m = pmeAtomGridIndex[i].x;
        real4 pos = posq[m];
        pos.x -= floor(pos.x*invPeriodicBoxSize.x)*periodicBoxSize.x;
        pos.y -= floor(pos.y*invPeriodicBoxSize.y)*periodicBoxSize.y;
        pos.z -= floor(pos.z*invPeriodicBoxSize.z)*periodicBoxSize.z;

        // Since we need the full set of thetas, it's faster to compute them here than load them
        // from global memory.

        real w = pos.x*invPeriodicBoxSize.x;
        real fr = GRID_SIZE_X*(w-(int)(w+0.5f)+0.5f);
        int ifr = (int) fr;
        w = fr - ifr;
        int igrid1 = ifr-PME_ORDER+1;
        computeBSplinePoint(theta1, w, array);
        w = pos.y*invPeriodicBoxSize.y;
        fr = GRID_SIZE_Y*(w-(int)(w+0.5f)+0.5f);
        ifr = (int) fr;
        w = fr - ifr;
        int igrid2 = ifr-PME_ORDER+1;
        computeBSplinePoint(theta2, w, array);
        w = pos.z*invPeriodicBoxSize.z;
        fr = GRID_SIZE_Z*(w-(int)(w+0.5f)+0.5f);
        ifr = (int) fr;
        w = fr - ifr;
        int igrid3 = ifr-PME_ORDER+1;
        computeBSplinePoint(theta3, w, array);
        igrid1 += (igrid1 < 0 ? GRID_SIZE_X : 0);
        igrid2 += (igrid2 < 0 ? GRID_SIZE_Y : 0);
        igrid3 += (igrid3 < 0 ? GRID_SIZE_Z : 0);

        // Compute the potential from this grid point.

        real tuv000 = 0;
        real tuv001 = 0;
        real tuv010 = 0;
        real tuv100 = 0;
        real tuv200 = 0;
        real tuv020 = 0;
        real tuv002 = 0;
        real tuv110 = 0;
        real tuv101 = 0;
        real tuv011 = 0;
        real tuv300 = 0;
        real tuv030 = 0;
        real tuv003 = 0;
        real tuv210 = 0;
        real tuv201 = 0;
        real tuv120 = 0;
        real tuv021 = 0;
        real tuv102 = 0;
        real tuv012 = 0;
        real tuv111 = 0;
        for (int iz = 0; iz < PME_ORDER; iz++) {
            int k = igrid3+iz-(igrid3+iz >= GRID_SIZE_Z ? GRID_SIZE_Z : 0);
            real4 v = theta3[iz];
            real tu00 = 0;
            real tu10 = 0;
            real tu01 = 0;
            real tu20 = 0;
            real tu11 = 0;
            real tu02 = 0;
            real tu30 = 0;
            real tu21 = 0;
            real tu12 = 0;
            real tu03 = 0;
            for (int iy = 0; iy < PME_ORDER; iy++) {
                int j = igrid2+iy-(igrid2+iy >= GRID_SIZE_Y ? GRID_SIZE_Y : 0);
                real4 u = theta2[iy];
                real4 t = make_real4(0, 0, 0, 0);
                for (int ix = 0; ix < PME_ORDER; ix++) {
                    int i = igrid1+ix-(igrid1+ix >= GRID_SIZE_X ? GRID_SIZE_X : 0);
                    int gridIndex = i*GRID_SIZE_Y*GRID_SIZE_Z + j*GRID_SIZE_Z + k;
                    real tq = pmeGrid[gridIndex].x;
                    real4 tadd = theta1[ix];
                    t.x += tq*tadd.x;
                    t.y += tq*tadd.y;
                    t.z += tq*tadd.z;
                    t.w += tq*tadd.w;
                }
                tu00 += t.x*u.x;
                tu10 += t.y*u.x;
                tu01 += t.x*u.y;
                tu20 += t.z*u.x;
                tu11 += t.y*u.y;
                tu02 += t.x*u.z;
                tu30 += t.w*u.x;
                tu21 += t.z*u.y;
                tu12 += t.y*u.z;
                tu03 += t.x*u.w;
            }
            tuv000 += tu00*v.x;
            tuv100 += tu10*v.x;
            tuv010 += tu01*v.x;
            tuv001 += tu00*v.y;
            tuv200 += tu20*v.x;
            tuv020 += tu02*v.x;
            tuv002 += tu00*v.z;
            tuv110 += tu11*v.x;
            tuv101 += tu10*v.y;
            tuv011 += tu01*v.y;
            tuv300 += tu30*v.x;
            tuv030 += tu03*v.x;
            tuv003 += tu00*v.w;
            tuv210 += tu21*v.x;
            tuv201 += tu20*v.y;
            tuv120 += tu12*v.x;
            tuv021 += tu02*v.y;
            tuv102 += tu10*v.z;
            tuv012 += tu01*v.z;
            tuv111 += tu11*v.y;
        }
        phi[20*m] = tuv000;
        phi[20*m+1] = tuv100;
        phi[20*m+2] = tuv010;
        phi[20*m+3] = tuv001;
        phi[20*m+4] = tuv200;
        phi[20*m+5] = tuv020;
        phi[20*m+6] = tuv002;
        phi[20*m+7] = tuv110;
        phi[20*m+8] = tuv101;
        phi[20*m+9] = tuv011;
        phi[20*m+10] = tuv300;
        phi[20*m+11] = tuv030;
        phi[20*m+12] = tuv003;
        phi[20*m+13] = tuv210;
        phi[20*m+14] = tuv201;
        phi[20*m+15] = tuv120;
        phi[20*m+16] = tuv021;
        phi[20*m+17] = tuv102;
        phi[20*m+18] = tuv012;
        phi[20*m+19] = tuv111;
        real dipoleScale = (4/(real) 3)*(EWALD_ALPHA*EWALD_ALPHA*EWALD_ALPHA)/SQRT_PI;
        long long fieldx = (long long) ((dipoleScale*labFrameDipole[m*3]-GRID_SIZE_X*invPeriodicBoxSize.x*tuv100)*0x100000000);
        fieldBuffers[m] = fieldx;
        fieldPolarBuffers[m] = fieldx;
        long long fieldy = (long long) ((dipoleScale*labFrameDipole[m*3+1]-GRID_SIZE_Y*invPeriodicBoxSize.y*tuv010)*0x100000000);
        fieldBuffers[m+PADDED_NUM_ATOMS] = fieldy;
        fieldPolarBuffers[m+PADDED_NUM_ATOMS] = fieldy;
        long long fieldz = (long long) ((dipoleScale*labFrameDipole[m*3+2]-GRID_SIZE_Z*invPeriodicBoxSize.z*tuv001)*0x100000000);
        fieldBuffers[m+2*PADDED_NUM_ATOMS] = fieldz;
        fieldPolarBuffers[m+2*PADDED_NUM_ATOMS] = fieldz;
    }
}

extern "C" __global__ void computeInducedPotentialFromGrid(const real2* __restrict__ pmeGrid, real* __restrict__ phid,
        real* __restrict__ phip, real* __restrict__ phidp, const real4* __restrict__ posq,
        real4 periodicBoxSize, real4 invPeriodicBoxSize, int2* __restrict__ pmeAtomGridIndex) {
    real array[PME_ORDER*PME_ORDER];
    real4 theta1[PME_ORDER];
    real4 theta2[PME_ORDER];
    real4 theta3[PME_ORDER];
    
    // Process the atoms in spatially sorted order.  This improves cache performance when loading
    // the grid values.
    
    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < NUM_ATOMS; i += blockDim.x*gridDim.x) {
        int m = pmeAtomGridIndex[i].x;
        real4 pos = posq[m];
        pos.x -= floor(pos.x*invPeriodicBoxSize.x)*periodicBoxSize.x;
        pos.y -= floor(pos.y*invPeriodicBoxSize.y)*periodicBoxSize.y;
        pos.z -= floor(pos.z*invPeriodicBoxSize.z)*periodicBoxSize.z;

        // Since we need the full set of thetas, it's faster to compute them here than load them
        // from global memory.

        real w = pos.x*invPeriodicBoxSize.x;
        real fr = GRID_SIZE_X*(w-(int)(w+0.5f)+0.5f);
        int ifr = (int) fr;
        w = fr - ifr;
        int igrid1 = ifr-PME_ORDER+1;
        computeBSplinePoint(theta1, w, array);
        w = pos.y*invPeriodicBoxSize.y;
        fr = GRID_SIZE_Y*(w-(int)(w+0.5f)+0.5f);
        ifr = (int) fr;
        w = fr - ifr;
        int igrid2 = ifr-PME_ORDER+1;
        computeBSplinePoint(theta2, w, array);
        w = pos.z*invPeriodicBoxSize.z;
        fr = GRID_SIZE_Z*(w-(int)(w+0.5f)+0.5f);
        ifr = (int) fr;
        w = fr - ifr;
        int igrid3 = ifr-PME_ORDER+1;
        computeBSplinePoint(theta3, w, array);
        igrid1 += (igrid1 < 0 ? GRID_SIZE_X : 0);
        igrid2 += (igrid2 < 0 ? GRID_SIZE_Y : 0);
        igrid3 += (igrid3 < 0 ? GRID_SIZE_Z : 0);

        // Compute the potential from this grid point.

        real tuv100_1 = 0;
        real tuv010_1 = 0;
        real tuv001_1 = 0;
        real tuv200_1 = 0;
        real tuv020_1 = 0;
        real tuv002_1 = 0;
        real tuv110_1 = 0;
        real tuv101_1 = 0;
        real tuv011_1 = 0;
        real tuv100_2 = 0;
        real tuv010_2 = 0;
        real tuv001_2 = 0;
        real tuv200_2 = 0;
        real tuv020_2 = 0;
        real tuv002_2 = 0;
        real tuv110_2 = 0;
        real tuv101_2 = 0;
        real tuv011_2 = 0;
        real tuv000 = 0;
        real tuv001 = 0;
        real tuv010 = 0;
        real tuv100 = 0;
        real tuv200 = 0;
        real tuv020 = 0;
        real tuv002 = 0;
        real tuv110 = 0;
        real tuv101 = 0;
        real tuv011 = 0;
        real tuv300 = 0;
        real tuv030 = 0;
        real tuv003 = 0;
        real tuv210 = 0;
        real tuv201 = 0;
        real tuv120 = 0;
        real tuv021 = 0;
        real tuv102 = 0;
        real tuv012 = 0;
        real tuv111 = 0;
        for (int iz = 0; iz < PME_ORDER; iz++) {
            int k = igrid3+iz-(igrid3+iz >= GRID_SIZE_Z ? GRID_SIZE_Z : 0);
            real4 v = theta3[iz];
            real tu00_1 = 0;
            real tu01_1 = 0;
            real tu10_1 = 0;
            real tu20_1 = 0;
            real tu11_1 = 0;
            real tu02_1 = 0;
            real tu00_2 = 0;
            real tu01_2 = 0;
            real tu10_2 = 0;
            real tu20_2 = 0;
            real tu11_2 = 0;
            real tu02_2 = 0;
            real tu00 = 0;
            real tu10 = 0;
            real tu01 = 0;
            real tu20 = 0;
            real tu11 = 0;
            real tu02 = 0;
            real tu30 = 0;
            real tu21 = 0;
            real tu12 = 0;
            real tu03 = 0;
            for (int iy = 0; iy < PME_ORDER; iy++) {
                int j = igrid2+iy-(igrid2+iy >= GRID_SIZE_Y ? GRID_SIZE_Y : 0);
                real4 u = theta2[iy];
                real t0_1 = 0;
                real t1_1 = 0;
                real t2_1 = 0;
                real t0_2 = 0;
                real t1_2 = 0;
                real t2_2 = 0;
                real t3 = 0;
                for (int ix = 0; ix < PME_ORDER; ix++) {
                    int i = igrid1+ix-(igrid1+ix >= GRID_SIZE_X ? GRID_SIZE_X : 0);
                    int gridIndex = i*GRID_SIZE_Y*GRID_SIZE_Z + j*GRID_SIZE_Z + k;
                    real2 tq = pmeGrid[gridIndex];
                    real4 tadd = theta1[ix];
                    t0_1 += tq.x*tadd.x;
                    t1_1 += tq.x*tadd.y;
                    t2_1 += tq.x*tadd.z;
                    t0_2 += tq.y*tadd.x;
                    t1_2 += tq.y*tadd.y;
                    t2_2 += tq.y*tadd.z;
                    t3 += (tq.x+tq.y)*tadd.w;
                }
                tu00_1 += t0_1*u.x;
                tu10_1 += t1_1*u.x;
                tu01_1 += t0_1*u.y;
                tu20_1 += t2_1*u.x;
                tu11_1 += t1_1*u.y;
                tu02_1 += t0_1*u.z;
                tu00_2 += t0_2*u.x;
                tu10_2 += t1_2*u.x;
                tu01_2 += t0_2*u.y;
                tu20_2 += t2_2*u.x;
                tu11_2 += t1_2*u.y;
                tu02_2 += t0_2*u.z;
                real t0 = t0_1 + t0_2;
                real t1 = t1_1 + t1_2;
                real t2 = t2_1 + t2_2;
                tu00 += t0*u.x;
                tu10 += t1*u.x;
                tu01 += t0*u.y;
                tu20 += t2*u.x;
                tu11 += t1*u.y;
                tu02 += t0*u.z;
                tu30 += t3*u.x;
                tu21 += t2*u.y;
                tu12 += t1*u.z;
                tu03 += t0*u.w;
            }
            tuv100_1 += tu10_1*v.x;
            tuv010_1 += tu01_1*v.x;
            tuv001_1 += tu00_1*v.y;
            tuv200_1 += tu20_1*v.x;
            tuv020_1 += tu02_1*v.x;
            tuv002_1 += tu00_1*v.z;
            tuv110_1 += tu11_1*v.x;
            tuv101_1 += tu10_1*v.y;
            tuv011_1 += tu01_1*v.y;
            tuv100_2 += tu10_2*v.x;
            tuv010_2 += tu01_2*v.x;
            tuv001_2 += tu00_2*v.y;
            tuv200_2 += tu20_2*v.x;
            tuv020_2 += tu02_2*v.x;
            tuv002_2 += tu00_2*v.z;
            tuv110_2 += tu11_2*v.x;
            tuv101_2 += tu10_2*v.y;
            tuv011_2 += tu01_2*v.y;
            tuv000 += tu00*v.x;
            tuv100 += tu10*v.x;
            tuv010 += tu01*v.x;
            tuv001 += tu00*v.y;
            tuv200 += tu20*v.x;
            tuv020 += tu02*v.x;
            tuv002 += tu00*v.z;
            tuv110 += tu11*v.x;
            tuv101 += tu10*v.y;
            tuv011 += tu01*v.y;
            tuv300 += tu30*v.x;
            tuv030 += tu03*v.x;
            tuv003 += tu00*v.w;
            tuv210 += tu21*v.x;
            tuv201 += tu20*v.y;
            tuv120 += tu12*v.x;
            tuv021 += tu02*v.y;
            tuv102 += tu10*v.z;
            tuv012 += tu01*v.z;
            tuv111 += tu11*v.y;
        }
        phid[10*m]   = 0;
        phid[10*m+1] = tuv100_1;
        phid[10*m+2] = tuv010_1;
        phid[10*m+3] = tuv001_1;
        phid[10*m+4] = tuv200_1;
        phid[10*m+5] = tuv020_1;
        phid[10*m+6] = tuv002_1;
        phid[10*m+7] = tuv110_1;
        phid[10*m+8] = tuv101_1;
        phid[10*m+9] = tuv011_1;

        phip[10*m]   = 0;
        phip[10*m+1] = tuv100_2;
        phip[10*m+2] = tuv010_2;
        phip[10*m+3] = tuv001_2;
        phip[10*m+4] = tuv200_2;
        phip[10*m+5] = tuv020_2;
        phip[10*m+6] = tuv002_2;
        phip[10*m+7] = tuv110_2;
        phip[10*m+8] = tuv101_2;
        phip[10*m+9] = tuv011_2;

        phidp[20*m] = tuv000;
        phidp[20*m+1] = tuv100;
        phidp[20*m+2] = tuv010;
        phidp[20*m+3] = tuv001;
        phidp[20*m+4] = tuv200;
        phidp[20*m+5] = tuv020;
        phidp[20*m+6] = tuv002;
        phidp[20*m+7] = tuv110;
        phidp[20*m+8] = tuv101;
        phidp[20*m+9] = tuv011;
        phidp[20*m+10] = tuv300;
        phidp[20*m+11] = tuv030;
        phidp[20*m+12] = tuv003;
        phidp[20*m+13] = tuv210;
        phidp[20*m+14] = tuv201;
        phidp[20*m+15] = tuv120;
        phidp[20*m+16] = tuv021;
        phidp[20*m+17] = tuv102;
        phidp[20*m+18] = tuv012;
        phidp[20*m+19] = tuv111;
    }
}

extern "C" __global__ void computeFixedMultipoleForceAndEnergy(real4* __restrict__ posq, unsigned long long* __restrict__ forceBuffers,
        long long* __restrict__ torqueBuffers, real* __restrict__ energyBuffer, const real* __restrict__ labFrameDipole,
        const real* __restrict__ labFrameQuadrupole, const real* __restrict__ phi_global, real4 invPeriodicBoxSize) {
    real multipole[10];
    const int deriv1[] = {1, 4, 7, 8, 10, 15, 17, 13, 14, 19};
    const int deriv2[] = {2, 7, 5, 9, 13, 11, 18, 15, 19, 16};
    const int deriv3[] = {3, 8, 9, 6, 14, 16, 12, 19, 17, 18};
    const real xscale = GRID_SIZE_X*invPeriodicBoxSize.x;
    const real yscale = GRID_SIZE_Y*invPeriodicBoxSize.y;
    const real zscale = GRID_SIZE_Z*invPeriodicBoxSize.z;
    real energy = 0;
    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < NUM_ATOMS; i += blockDim.x*gridDim.x) {
        // Compute the torque.

        multipole[0] = posq[i].w;
        multipole[1] = labFrameDipole[i*3];
        multipole[2] = labFrameDipole[i*3+1];
        multipole[3] = labFrameDipole[i*3+2];
        multipole[4] = labFrameQuadrupole[i*5];
        multipole[5] = labFrameQuadrupole[i*5+3];
        multipole[6] = -(multipole[4]+multipole[5]);
        multipole[7] = 2*labFrameQuadrupole[i*5+1];
        multipole[8] = 2*labFrameQuadrupole[i*5+2];
        multipole[9] = 2*labFrameQuadrupole[i*5+4];

        const real* phi = &phi_global[20*i];

        torqueBuffers[i] = (long long) (EPSILON_FACTOR*(multipole[3]*yscale*phi[2] - multipole[2]*zscale*phi[3]
                      + 2*(multipole[6]-multipole[5])*yscale*zscale*phi[9]
                      + multipole[8]*xscale*yscale*phi[7] + multipole[9]*yscale*yscale*phi[5]
                      - multipole[7]*xscale*zscale*phi[8] - multipole[9]*zscale*zscale*phi[6])*0x100000000);

        torqueBuffers[i+PADDED_NUM_ATOMS] = (long long) (EPSILON_FACTOR*(multipole[1]*zscale*phi[3] - multipole[3]*xscale*phi[1]
                      + 2*(multipole[4]-multipole[6])*xscale*zscale*phi[8]
                      + multipole[7]*yscale*zscale*phi[9] + multipole[8]*zscale*zscale*phi[6]
                      - multipole[8]*xscale*xscale*phi[4] - multipole[9]*xscale*yscale*phi[7])*0x100000000);

        torqueBuffers[i+PADDED_NUM_ATOMS*2] = (long long) (EPSILON_FACTOR*(multipole[2]*xscale*phi[1] - multipole[1]*yscale*phi[2]
                      + 2*(multipole[5]-multipole[4])*xscale*yscale*phi[7]
                      + multipole[7]*xscale*xscale*phi[4] + multipole[9]*xscale*zscale*phi[8]
                      - multipole[7]*yscale*yscale*phi[5] - multipole[8]*yscale*zscale*phi[9])*0x100000000);

        // Compute the force and energy.

        multipole[1] *= xscale;
        multipole[2] *= yscale;
        multipole[3] *= zscale;
        multipole[4] *= xscale*xscale;
        multipole[5] *= yscale*yscale;
        multipole[6] *= zscale*zscale;
        multipole[7] *= xscale*yscale;
        multipole[8] *= xscale*zscale;
        multipole[9] *= yscale*zscale;

        real4 f = make_real4(0, 0, 0, 0);
        for (int k = 0; k < 10; k++) {
            energy += multipole[k]*phi[k];
            f.x += multipole[k]*phi[deriv1[k]];
            f.y += multipole[k]*phi[deriv2[k]];
            f.z += multipole[k]*phi[deriv3[k]];
        }
        f.x *= EPSILON_FACTOR*xscale;
        f.y *= EPSILON_FACTOR*yscale;
        f.z *= EPSILON_FACTOR*zscale;
        forceBuffers[i] -= static_cast<unsigned long long>((long long) (f.x*0x100000000));
        forceBuffers[i+PADDED_NUM_ATOMS] -= static_cast<unsigned long long>((long long) (f.y*0x100000000));
        forceBuffers[i+PADDED_NUM_ATOMS*2] -= static_cast<unsigned long long>((long long) (f.z*0x100000000));
    }
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += 0.5f*EPSILON_FACTOR*energy;
}

extern "C" __global__ void computeInducedDipoleForceAndEnergy(real4* __restrict__ posq, unsigned long long* __restrict__ forceBuffers,
        long long* __restrict__ torqueBuffers, real* __restrict__ energyBuffer, const real* __restrict__ labFrameDipole,
        const real* __restrict__ labFrameQuadrupole, const real* __restrict__ inducedDipole_global, const real* __restrict__ inducedDipolePolar_global,
        const real* __restrict__ phi_global, const real* __restrict__ phid_global, const real* __restrict__ phip_global,
        const real* __restrict__ phidp_global, real4 invPeriodicBoxSize) {
    real multipole[10];
    real inducedDipole[3];
    real inducedDipolePolar[3];
    real scales[3];
    const int deriv1[] = {1, 4, 7, 8, 10, 15, 17, 13, 14, 19};
    const int deriv2[] = {2, 7, 5, 9, 13, 11, 18, 15, 19, 16};
    const int deriv3[] = {3, 8, 9, 6, 14, 16, 12, 19, 17, 18};
    const real xscale = GRID_SIZE_X*invPeriodicBoxSize.x;
    const real yscale = GRID_SIZE_Y*invPeriodicBoxSize.y;
    const real zscale = GRID_SIZE_Z*invPeriodicBoxSize.z;
    scales[0] = xscale;
    scales[1] = yscale;
    scales[2] = zscale;
    real energy = 0;
    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < NUM_ATOMS; i += blockDim.x*gridDim.x) {
        // Compute the torque.

        multipole[0] = posq[i].w;
        multipole[1] = labFrameDipole[i*3];
        multipole[2] = labFrameDipole[i*3+1];
        multipole[3] = labFrameDipole[i*3+2];
        multipole[4] = labFrameQuadrupole[i*5];
        multipole[5] = labFrameQuadrupole[i*5+3];
        multipole[6] = -(multipole[4]+multipole[5]);
        multipole[7] = 2*labFrameQuadrupole[i*5+1];
        multipole[8] = 2*labFrameQuadrupole[i*5+2];
        multipole[9] = 2*labFrameQuadrupole[i*5+4];
        const real* phidp = &phidp_global[20*i];
 
        torqueBuffers[i] += (long long) (0.5f*EPSILON_FACTOR*(multipole[3]*yscale*phidp[2] - multipole[2]*zscale*phidp[3]
                      + 2*(multipole[6]-multipole[5])*yscale*zscale*phidp[9]
                      + multipole[8]*xscale*yscale*phidp[7] + multipole[9]*yscale*yscale*phidp[5]
                      - multipole[7]*xscale*zscale*phidp[8] - multipole[9]*zscale*zscale*phidp[6])*0x100000000);

        torqueBuffers[i+PADDED_NUM_ATOMS] += (long long) (0.5f*EPSILON_FACTOR*(multipole[1]*zscale*phidp[3] - multipole[3]*xscale*phidp[1]
                      + 2*(multipole[4]-multipole[6])*xscale*zscale*phidp[8]
                      + multipole[7]*yscale*zscale*phidp[9] + multipole[8]*zscale*zscale*phidp[6]
                      - multipole[8]*xscale*xscale*phidp[4] - multipole[9]*xscale*yscale*phidp[7])*0x100000000);

        torqueBuffers[i+PADDED_NUM_ATOMS*2] += (long long) (0.5f*EPSILON_FACTOR*(multipole[2]*xscale*phidp[1] - multipole[1]*yscale*phidp[2]
                      + 2*(multipole[5]-multipole[4])*xscale*yscale*phidp[7]
                      + multipole[7]*xscale*xscale*phidp[4] + multipole[9]*xscale*zscale*phidp[8]
                      - multipole[7]*yscale*yscale*phidp[5] - multipole[8]*yscale*zscale*phidp[9])*0x100000000);

        // Compute the force and energy.

        multipole[1] *= xscale;
        multipole[2] *= yscale;
        multipole[3] *= zscale;
        multipole[4] *= xscale*xscale;
        multipole[5] *= yscale*yscale;
        multipole[6] *= zscale*zscale;
        multipole[7] *= xscale*yscale;
        multipole[8] *= xscale*zscale;
        multipole[9] *= yscale*zscale;

        inducedDipole[0] = inducedDipole_global[i*3];
        inducedDipole[1] = inducedDipole_global[i*3+1];
        inducedDipole[2] = inducedDipole_global[i*3+2];
        inducedDipolePolar[0] = inducedDipolePolar_global[i*3];
        inducedDipolePolar[1] = inducedDipolePolar_global[i*3+1];
        inducedDipolePolar[2] = inducedDipolePolar_global[i*3+2];
        const real* phi = &phi_global[20*i];
        const real* phip = &phip_global[10*i];
        const real* phid = &phid_global[10*i];
        real4 f = make_real4(0, 0, 0, 0);

        energy += GRID_SIZE_X*invPeriodicBoxSize.x*inducedDipole[0]*phi[1];
        energy += GRID_SIZE_Y*invPeriodicBoxSize.y*inducedDipole[1]*phi[2];
        energy += GRID_SIZE_Z*invPeriodicBoxSize.z*inducedDipole[2]*phi[3];

        for (int k = 0; k < 3; k++) {
            int j1 = deriv1[k+1];
            int j2 = deriv2[k+1];
            int j3 = deriv3[k+1];
            f.x += (inducedDipole[k]+inducedDipolePolar[k])*phi[j1]*(scales[k]/xscale);
            f.y += (inducedDipole[k]+inducedDipolePolar[k])*phi[j2]*(scales[k]/yscale);
            f.z += (inducedDipole[k]+inducedDipolePolar[k])*phi[j3]*(scales[k]/zscale);
#ifndef DIRECT_POLARIZATION
            f.x += (inducedDipole[k]*phip[j1] + inducedDipolePolar[k]*phid[j1])*(scales[k]/xscale);
            f.y += (inducedDipole[k]*phip[j2] + inducedDipolePolar[k]*phid[j2])*(scales[k]/yscale);
            f.z += (inducedDipole[k]*phip[j3] + inducedDipolePolar[k]*phid[j3])*(scales[k]/zscale);
#endif
        }

        f.x *= GRID_SIZE_X*invPeriodicBoxSize.x;
        f.y *= GRID_SIZE_Y*invPeriodicBoxSize.y;
        f.z *= GRID_SIZE_Z*invPeriodicBoxSize.z;
        for (int k = 0; k < 10; k++) {
            f.x += multipole[k]*phidp[deriv1[k]];
            f.y += multipole[k]*phidp[deriv2[k]];
            f.z += multipole[k]*phidp[deriv3[k]];
        }
        f.x *= 0.5f*EPSILON_FACTOR*xscale;
        f.y *= 0.5f*EPSILON_FACTOR*yscale;
        f.z *= 0.5f*EPSILON_FACTOR*zscale;
        forceBuffers[i] -= static_cast<unsigned long long>((long long) (f.x*0x100000000));
        forceBuffers[i+PADDED_NUM_ATOMS] -= static_cast<unsigned long long>((long long) (f.y*0x100000000));
        forceBuffers[i+PADDED_NUM_ATOMS*2] -= static_cast<unsigned long long>((long long) (f.z*0x100000000));
    }
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += 0.5f*EPSILON_FACTOR*energy;
}

extern "C" __global__ void recordInducedFieldDipoles(const real* __restrict__ phid, real* const __restrict__ phip,
        long long* __restrict__ inducedField, long long* __restrict__ inducedFieldPolar, real4 invPeriodicBoxSize) {
    real xscale = GRID_SIZE_X*invPeriodicBoxSize.x*0x100000000;
    real yscale = GRID_SIZE_Y*invPeriodicBoxSize.y*0x100000000;
    real zscale = GRID_SIZE_Z*invPeriodicBoxSize.z*0x100000000;
    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < NUM_ATOMS; i += blockDim.x*gridDim.x) {
        inducedField[i] -= (long long) (xscale*phid[10*i+1]);
        inducedField[i+PADDED_NUM_ATOMS] -= (long long) (yscale*phid[10*i+2]);
        inducedField[i+PADDED_NUM_ATOMS*2] -= (long long) (zscale*phid[10*i+3]);
        inducedFieldPolar[i] -= (long long) (xscale*phip[10*i+1]);
        inducedFieldPolar[i+PADDED_NUM_ATOMS] -= (long long) (yscale*phip[10*i+2]);
        inducedFieldPolar[i+PADDED_NUM_ATOMS*2] -= (long long) (zscale*phip[10*i+3]);
    }
}
