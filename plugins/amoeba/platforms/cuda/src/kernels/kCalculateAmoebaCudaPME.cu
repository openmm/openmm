//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------

#include "amoebaGpuTypes.h"
#include "amoebaCudaKernels.h"
#include "bbsort.h"

//#define AMOEBA_DEBUG

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;

/* Cuda compiler on Windows does not recognized "static const float" values */
#define LOCAL_HACK_PI 3.1415926535897932384626433832795f

void SetCalculateAmoebaPMESim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));
    RTERROR(status, "SetCalculateAmoebaPMESim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));
    RTERROR(status, "SetCalculateAmoebaPMESim: cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

#define ARRAY(x,y) array[(x)-1+((y)-1)*AMOEBA_PME_ORDER]

/**
 * This is called from computeBsplines().  It calculates the spline coefficients for a single atom along a single axis.
 */
__device__ void computeBSplinePoint(float4* thetai, float w, float* array)
{
    // initialization to get to 2nd order recursion

    ARRAY(2,2) = w;
    ARRAY(2,1) = 1.0f - w;

    // perform one pass to get to 3rd order recursion

    ARRAY(3,3) = 0.5f * w * ARRAY(2,2);
    ARRAY(3,2) = 0.5f * ((1.0f+w)*ARRAY(2,1)+(2.0f-w)*ARRAY(2,2));
    ARRAY(3,1) = 0.5f * (1.0f-w) * ARRAY(2,1);

    // compute standard B-spline recursion to desired order

    for (int i = 4; i <= AMOEBA_PME_ORDER; i++)
    {
        int k = i - 1;
        float denom = 1.0f / k;
        ARRAY(i,i) = denom * w * ARRAY(k,k);
        for (int j = 1; j <= i-2; j++)
            ARRAY(i,i-j) = denom * ((w+j)*ARRAY(k,i-j-1)+(i-j-w)*ARRAY(k,i-j));
        ARRAY(i,1) = denom * (1.0f-w) * ARRAY(k,1);
    }

    // get coefficients for the B-spline first derivative

    int k = AMOEBA_PME_ORDER - 1;
    ARRAY(k,AMOEBA_PME_ORDER) = ARRAY(k,AMOEBA_PME_ORDER-1);
    for (int i = AMOEBA_PME_ORDER-1; i >= 2; i--)
        ARRAY(k,i) = ARRAY(k,i-1) - ARRAY(k,i);
    ARRAY(k,1) = -ARRAY(k,1);

    // get coefficients for the B-spline second derivative

    k = AMOEBA_PME_ORDER - 2;
    ARRAY(k,AMOEBA_PME_ORDER-1) = ARRAY(k,AMOEBA_PME_ORDER-2);
    for (int i = AMOEBA_PME_ORDER-2; i >= 2; i--)
        ARRAY(k,i) = ARRAY(k,i-1) - ARRAY(k,i);
    ARRAY(k,1) = -ARRAY(k,1);
    ARRAY(k,AMOEBA_PME_ORDER) = ARRAY(k,AMOEBA_PME_ORDER-1);
    for (int i = AMOEBA_PME_ORDER-1; i >= 2; i--)
        ARRAY(k,i) = ARRAY(k,i-1) - ARRAY(k,i);
    ARRAY(k,1) = -ARRAY(k,1);

    // get coefficients for the B-spline third derivative

    k = AMOEBA_PME_ORDER - 3;
    ARRAY(k,AMOEBA_PME_ORDER-2) = ARRAY(k,AMOEBA_PME_ORDER-3);
    for (int i = AMOEBA_PME_ORDER-3; i >= 2; i--)
        ARRAY(k,i) = ARRAY(k,i-1) - ARRAY(k,i);
    ARRAY(k,1) = -ARRAY(k,1);
    ARRAY(k,AMOEBA_PME_ORDER-1) = ARRAY(k,AMOEBA_PME_ORDER-2);
    for (int i = AMOEBA_PME_ORDER-2; i >= 2; i--)
        ARRAY(k,i) = ARRAY(k,i-1) - ARRAY(k,i);
    ARRAY(k,1) = -ARRAY(k,1);
    ARRAY(k,AMOEBA_PME_ORDER) = ARRAY(k,AMOEBA_PME_ORDER-1);
    for (int i = AMOEBA_PME_ORDER-1; i >= 2; i--)
        ARRAY(k,i) = ARRAY(k,i-1) - ARRAY(k,i);
    ARRAY(k,1) = -ARRAY(k,1);

    // copy coefficients from temporary to permanent storage

    for (int i = 1; i <= AMOEBA_PME_ORDER; i++)
        thetai[i-1] = make_float4(ARRAY(AMOEBA_PME_ORDER,i), ARRAY(AMOEBA_PME_ORDER-1,i), ARRAY(AMOEBA_PME_ORDER-2,i), ARRAY(AMOEBA_PME_ORDER-3,i));
}

/**
 * Compute bspline coefficients.
 */
__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(448, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(160, 1)
#else
__launch_bounds__(160, 1)
#endif
void kComputeAmoebaBsplines_kernel()
{
    extern __shared__ float bsplines_cache[]; // size = block_size*pme_order*pme_order
    float* array = &bsplines_cache[threadIdx.x*AMOEBA_PME_ORDER*AMOEBA_PME_ORDER];

    //  get the B-spline coefficients for each multipole site

    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < cSim.atoms; i += blockDim.x*gridDim.x) {
        float4 posq = cSim.pPosq[i];
        posq.x -= floor(posq.x*cSim.invPeriodicBoxSizeX)*cSim.periodicBoxSizeX;
        posq.y -= floor(posq.y*cSim.invPeriodicBoxSizeY)*cSim.periodicBoxSizeY;
        posq.z -= floor(posq.z*cSim.invPeriodicBoxSizeZ)*cSim.periodicBoxSizeZ;

        // First axis.

        float w = posq.x*cSim.invPeriodicBoxSizeX;
        float fr = cSim.pmeGridSize.x*(w-(int)(w+0.5f)+0.5f);
        int ifr = (int) fr;
        w = fr - ifr;
        int igrid1 = ifr-AMOEBA_PME_ORDER+1;
        computeBSplinePoint(&cAmoebaSim.pThetai1[i*AMOEBA_PME_ORDER], w, array);

        // Second axis.

        w = posq.y*cSim.invPeriodicBoxSizeY;
        fr = cSim.pmeGridSize.y*(w-(int)(w+0.5f)+0.5f);
        ifr = (int) fr;
        w = fr - ifr;
        int igrid2 = ifr-AMOEBA_PME_ORDER+1;
        computeBSplinePoint(&cAmoebaSim.pThetai2[i*AMOEBA_PME_ORDER], w, array);

        // Third axis.

        w = posq.z*cSim.invPeriodicBoxSizeZ;
        fr = cSim.pmeGridSize.z*(w-(int)(w+0.5f)+0.5f);
        ifr = (int) fr;
        w = fr - ifr;
        int igrid3 = ifr-AMOEBA_PME_ORDER+1;
        computeBSplinePoint(&cAmoebaSim.pThetai3[i*AMOEBA_PME_ORDER], w, array);

        // Record the grid point.

        igrid1 += (igrid1 < 0 ? cSim.pmeGridSize.x : 0);
        igrid2 += (igrid2 < 0 ? cSim.pmeGridSize.y : 0);
        igrid3 += (igrid3 < 0 ? cSim.pmeGridSize.z : 0);
        cAmoebaSim.pIgrid[i] = make_int4(igrid1, igrid2, igrid3, 0);
        cSim.pPmeAtomGridIndex[i] = make_int2(i, igrid1*cSim.pmeGridSize.y*cSim.pmeGridSize.z+igrid2*cSim.pmeGridSize.z+igrid3);
    }
}

/**
 * For each grid point, find the range of sorted atoms associated with that point.
 */
__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(1024, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(512, 1)
#else
__launch_bounds__(256, 1)
#endif
void kFindAmoebaAtomRangeForGrid_kernel()
{
    int thread = blockIdx.x*blockDim.x+threadIdx.x;
    int start = (cSim.atoms*thread)/(blockDim.x*gridDim.x);
    int end = (cSim.atoms*(thread+1))/(blockDim.x*gridDim.x);
    int last = (start == 0 ? -1 : cSim.pPmeAtomGridIndex[start-1].y);
    for (int i = start; i < end; ++i)
    {
        int2 atomData = cSim.pPmeAtomGridIndex[i];
        int gridIndex = atomData.y;
        if (gridIndex != last)
        {
            for (int j = last+1; j <= gridIndex; ++j)
                cSim.pPmeAtomRange[j] = i;
            last = gridIndex;
        }

        // The grid index won't be needed again.  Reuse that component to hold the z index, thus saving
        // some work in the charge spreading kernel.

        float posz = cSim.pPosq[atomData.x].z;
        posz -= floor(posz*cSim.invPeriodicBoxSizeZ)*cSim.periodicBoxSizeZ;
        float w = posz*cSim.invPeriodicBoxSizeZ;
        float fr = cSim.pmeGridSize.z*(w-(int)(w+0.5f)+0.5f);
        int z = ((int) fr)-AMOEBA_PME_ORDER+1;
        cSim.pPmeAtomGridIndex[i].y = z;
    }

    // Fill in values beyond the last atom.

    if (thread == blockDim.x*gridDim.x-1)
    {
        int gridSize = cSim.pmeGridSize.x*cSim.pmeGridSize.y*cSim.pmeGridSize.z;
        for (int j = last+1; j <= gridSize; ++j)
            cSim.pPmeAtomRange[j] = cSim.atoms;
    }
}
__global__
__launch_bounds__(64, 10)
void kGridSpreadFixedMultipoles_kernel()
{
    const float xscale = cSim.pmeGridSize.x*cSim.invPeriodicBoxSizeX;
    const float yscale = cSim.pmeGridSize.y*cSim.invPeriodicBoxSizeY;
    const float zscale = cSim.pmeGridSize.z*cSim.invPeriodicBoxSizeZ;
    unsigned int numGridPoints = cSim.pmeGridSize.x*cSim.pmeGridSize.y*cSim.pmeGridSize.z;
    unsigned int numThreads = gridDim.x*blockDim.x;
    for (int gridIndex = blockIdx.x*blockDim.x+threadIdx.x; gridIndex < numGridPoints; gridIndex += numThreads)
    {
        int3 gridPoint;
        gridPoint.x = gridIndex/(cSim.pmeGridSize.y*cSim.pmeGridSize.z);
        int remainder = gridIndex-gridPoint.x*cSim.pmeGridSize.y*cSim.pmeGridSize.z;
        gridPoint.y = remainder/cSim.pmeGridSize.z;
        gridPoint.z = remainder-gridPoint.y*cSim.pmeGridSize.z;
        float result = 0.0f;
        for (int ix = 0; ix < AMOEBA_PME_ORDER; ++ix)
        {
            int x = gridPoint.x-ix+(gridPoint.x >= ix ? 0 : cSim.pmeGridSize.x);
            for (int iy = 0; iy < AMOEBA_PME_ORDER; ++iy)
            {
                int y = gridPoint.y-iy+(gridPoint.y >= iy ? 0 : cSim.pmeGridSize.y);
                int z1 = gridPoint.z-AMOEBA_PME_ORDER+1;
                z1 += (z1 >= 0 ? 0 : cSim.pmeGridSize.z);
                int z2 = (z1 < gridPoint.z ? gridPoint.z : cSim.pmeGridSize.z-1);
                int gridIndex1 = x*cSim.pmeGridSize.y*cSim.pmeGridSize.z+y*cSim.pmeGridSize.z+z1;
                int gridIndex2 = x*cSim.pmeGridSize.y*cSim.pmeGridSize.z+y*cSim.pmeGridSize.z+z2;
                int firstAtom = cSim.pPmeAtomRange[gridIndex1];
                int lastAtom = cSim.pPmeAtomRange[gridIndex2+1];
                for (int i = firstAtom; i < lastAtom; ++i)
                {
                    int2 atomData = cSim.pPmeAtomGridIndex[i];
                    int atomIndex = atomData.x;
                    int z = atomData.y;
                    int iz = gridPoint.z-z+(gridPoint.z >= z ? 0 : cSim.pmeGridSize.z);
                    float atomCharge = cSim.pPosq[atomIndex].w;
                    float atomDipoleX = xscale*cAmoebaSim.pLabFrameDipole[atomIndex*3];
                    float atomDipoleY = yscale*cAmoebaSim.pLabFrameDipole[atomIndex*3+1];
                    float atomDipoleZ = zscale*cAmoebaSim.pLabFrameDipole[atomIndex*3+2];
                    float atomQuadrupoleXX = xscale*xscale*cAmoebaSim.pLabFrameQuadrupole[atomIndex*9];
                    float atomQuadrupoleXY = 2*xscale*yscale*cAmoebaSim.pLabFrameQuadrupole[atomIndex*9+1];
                    float atomQuadrupoleXZ = 2*xscale*zscale*cAmoebaSim.pLabFrameQuadrupole[atomIndex*9+2];
                    float atomQuadrupoleYY = yscale*yscale*cAmoebaSim.pLabFrameQuadrupole[atomIndex*9+4];
                    float atomQuadrupoleYZ = 2*yscale*zscale*cAmoebaSim.pLabFrameQuadrupole[atomIndex*9+5];
                    float atomQuadrupoleZZ = zscale*zscale*cAmoebaSim.pLabFrameQuadrupole[atomIndex*9+8];
                    float4 t = cAmoebaSim.pThetai1[atomIndex*AMOEBA_PME_ORDER+ix];
                    float4 u = cAmoebaSim.pThetai2[atomIndex*AMOEBA_PME_ORDER+iy];
                    float4 v = cAmoebaSim.pThetai3[atomIndex*AMOEBA_PME_ORDER+iz];
                    float term0 = atomCharge*u.x*v.x + atomDipoleY*u.y*v.x + atomDipoleZ*u.x*v.y + atomQuadrupoleYY*u.z*v.x + atomQuadrupoleZZ*u.x*v.z + atomQuadrupoleYZ*u.y*v.y;
                    float term1 = atomDipoleX*u.x*v.x + atomQuadrupoleXY*u.y*v.x + atomQuadrupoleXZ*u.x*v.y;
                    float term2 = atomQuadrupoleXX * u.x * v.x;
                    result += term0*t.x + term1*t.y + term2*t.z;
                }
                if (z1 > gridPoint.z)
                {
                    gridIndex1 = x*cSim.pmeGridSize.y*cSim.pmeGridSize.z+y*cSim.pmeGridSize.z;
                    gridIndex2 = x*cSim.pmeGridSize.y*cSim.pmeGridSize.z+y*cSim.pmeGridSize.z+gridPoint.z;
                    firstAtom = cSim.pPmeAtomRange[gridIndex1];
                    lastAtom = cSim.pPmeAtomRange[gridIndex2+1];
                    for (int i = firstAtom; i < lastAtom; ++i)
                    {
                        int2 atomData = cSim.pPmeAtomGridIndex[i];
                        int atomIndex = atomData.x;
                        int z = atomData.y;
                        int iz = gridPoint.z-z+(gridPoint.z >= z ? 0 : cSim.pmeGridSize.z);
                        float atomCharge = cSim.pPosq[atomIndex].w;
                        float atomDipoleX = xscale*cAmoebaSim.pLabFrameDipole[atomIndex*3];
                        float atomDipoleY = yscale*cAmoebaSim.pLabFrameDipole[atomIndex*3+1];
                        float atomDipoleZ = zscale*cAmoebaSim.pLabFrameDipole[atomIndex*3+2];
                        float atomQuadrupoleXX = xscale*xscale*cAmoebaSim.pLabFrameQuadrupole[atomIndex*9];
                        float atomQuadrupoleXY = 2*xscale*yscale*cAmoebaSim.pLabFrameQuadrupole[atomIndex*9+1];
                        float atomQuadrupoleXZ = 2*xscale*zscale*cAmoebaSim.pLabFrameQuadrupole[atomIndex*9+2];
                        float atomQuadrupoleYY = yscale*yscale*cAmoebaSim.pLabFrameQuadrupole[atomIndex*9+4];
                        float atomQuadrupoleYZ = 2*yscale*zscale*cAmoebaSim.pLabFrameQuadrupole[atomIndex*9+5];
                        float atomQuadrupoleZZ = zscale*zscale*cAmoebaSim.pLabFrameQuadrupole[atomIndex*9+8];
                        float4 t = cAmoebaSim.pThetai1[atomIndex*AMOEBA_PME_ORDER+ix];
                        float4 u = cAmoebaSim.pThetai2[atomIndex*AMOEBA_PME_ORDER+iy];
                        float4 v = cAmoebaSim.pThetai3[atomIndex*AMOEBA_PME_ORDER+iz];
                        float term0 = atomCharge*u.x*v.x + atomDipoleY*u.y*v.x + atomDipoleZ*u.x*v.y + atomQuadrupoleYY*u.z*v.x + atomQuadrupoleZZ*u.x*v.z + atomQuadrupoleYZ*u.y*v.y;
                        float term1 = atomDipoleX*u.x*v.x + atomQuadrupoleXY*u.y*v.x + atomQuadrupoleXZ*u.x*v.y;
                        float term2 = atomQuadrupoleXX * u.x * v.x;
                        result += term0*t.x + term1*t.y + term2*t.z;
                    }
                }
            }
        }
        cSim.pPmeGrid[gridIndex] = make_cuComplex(result, 0.0f);
    }
}

__global__
__launch_bounds__(64, 10)
void kGridSpreadInducedDipoles_kernel()
{
    const float xscale = cSim.pmeGridSize.x*cSim.invPeriodicBoxSizeX;
    const float yscale = cSim.pmeGridSize.y*cSim.invPeriodicBoxSizeY;
    const float zscale = cSim.pmeGridSize.z*cSim.invPeriodicBoxSizeZ;
    unsigned int numGridPoints = cSim.pmeGridSize.x*cSim.pmeGridSize.y*cSim.pmeGridSize.z;
    unsigned int numThreads = gridDim.x*blockDim.x;
    for (int gridIndex = blockIdx.x*blockDim.x+threadIdx.x; gridIndex < numGridPoints; gridIndex += numThreads)
    {
        int3 gridPoint;
        gridPoint.x = gridIndex/(cSim.pmeGridSize.y*cSim.pmeGridSize.z);
        int remainder = gridIndex-gridPoint.x*cSim.pmeGridSize.y*cSim.pmeGridSize.z;
        gridPoint.y = remainder/cSim.pmeGridSize.z;
        gridPoint.z = remainder-gridPoint.y*cSim.pmeGridSize.z;
        cufftComplex result = make_cuComplex(0.0f, 0.0f);
        for (int ix = 0; ix < AMOEBA_PME_ORDER; ++ix)
        {
            int x = gridPoint.x-ix+(gridPoint.x >= ix ? 0 : cSim.pmeGridSize.x);
            for (int iy = 0; iy < AMOEBA_PME_ORDER; ++iy)
            {
                int y = gridPoint.y-iy+(gridPoint.y >= iy ? 0 : cSim.pmeGridSize.y);
                int z1 = gridPoint.z-AMOEBA_PME_ORDER+1;
                z1 += (z1 >= 0 ? 0 : cSim.pmeGridSize.z);
                int z2 = (z1 < gridPoint.z ? gridPoint.z : cSim.pmeGridSize.z-1);
                int gridIndex1 = x*cSim.pmeGridSize.y*cSim.pmeGridSize.z+y*cSim.pmeGridSize.z+z1;
                int gridIndex2 = x*cSim.pmeGridSize.y*cSim.pmeGridSize.z+y*cSim.pmeGridSize.z+z2;
                int firstAtom = cSim.pPmeAtomRange[gridIndex1];
                int lastAtom = cSim.pPmeAtomRange[gridIndex2+1];
                for (int i = firstAtom; i < lastAtom; ++i)
                {
                    int2 atomData = cSim.pPmeAtomGridIndex[i];
                    int atomIndex = atomData.x;
                    int z = atomData.y;
                    int iz = gridPoint.z-z+(gridPoint.z >= z ? 0 : cSim.pmeGridSize.z);
                    float inducedDipoleX = xscale*cAmoebaSim.pInducedDipole[atomIndex*3];
                    float inducedDipoleY = yscale*cAmoebaSim.pInducedDipole[atomIndex*3+1];
                    float inducedDipoleZ = zscale*cAmoebaSim.pInducedDipole[atomIndex*3+2];
                    float inducedDipolePolarX = xscale*cAmoebaSim.pInducedDipolePolar[atomIndex*3];
                    float inducedDipolePolarY = yscale*cAmoebaSim.pInducedDipolePolar[atomIndex*3+1];
                    float inducedDipolePolarZ = zscale*cAmoebaSim.pInducedDipolePolar[atomIndex*3+2];
                    float4 t = cAmoebaSim.pThetai1[atomIndex*AMOEBA_PME_ORDER+ix];
                    float4 u = cAmoebaSim.pThetai2[atomIndex*AMOEBA_PME_ORDER+iy];
                    float4 v = cAmoebaSim.pThetai3[atomIndex*AMOEBA_PME_ORDER+iz];
                    float term01 = inducedDipoleY*u.y*v.x + inducedDipoleZ*u.x*v.y;
                    float term11 = inducedDipoleX*u.x*v.x;
                    float term02 = inducedDipolePolarY*u.y*v.x + inducedDipolePolarZ*u.x*v.y;
                    float term12 = inducedDipolePolarX*u.x*v.x;
                    result.x += term01*t.x + term11*t.y;
                    result.y += term02*t.x + term12*t.y;
                }
                if (z1 > gridPoint.z)
                {
                    gridIndex1 = x*cSim.pmeGridSize.y*cSim.pmeGridSize.z+y*cSim.pmeGridSize.z;
                    gridIndex2 = x*cSim.pmeGridSize.y*cSim.pmeGridSize.z+y*cSim.pmeGridSize.z+gridPoint.z;
                    firstAtom = cSim.pPmeAtomRange[gridIndex1];
                    lastAtom = cSim.pPmeAtomRange[gridIndex2+1];
                    for (int i = firstAtom; i < lastAtom; ++i)
                    {
                        int2 atomData = cSim.pPmeAtomGridIndex[i];
                        int atomIndex = atomData.x;
                        int z = atomData.y;
                        int iz = gridPoint.z-z+(gridPoint.z >= z ? 0 : cSim.pmeGridSize.z);
                        float inducedDipoleX = xscale*cAmoebaSim.pInducedDipole[atomIndex*3];
                        float inducedDipoleY = yscale*cAmoebaSim.pInducedDipole[atomIndex*3+1];
                        float inducedDipoleZ = zscale*cAmoebaSim.pInducedDipole[atomIndex*3+2];
                        float inducedDipolePolarX = xscale*cAmoebaSim.pInducedDipolePolar[atomIndex*3];
                        float inducedDipolePolarY = yscale*cAmoebaSim.pInducedDipolePolar[atomIndex*3+1];
                        float inducedDipolePolarZ = zscale*cAmoebaSim.pInducedDipolePolar[atomIndex*3+2];
                        float4 t = cAmoebaSim.pThetai1[atomIndex*AMOEBA_PME_ORDER+ix];
                        float4 u = cAmoebaSim.pThetai2[atomIndex*AMOEBA_PME_ORDER+iy];
                        float4 v = cAmoebaSim.pThetai3[atomIndex*AMOEBA_PME_ORDER+iz];
                        float term01 = inducedDipoleY*u.y*v.x + inducedDipoleZ*u.x*v.y;
                        float term11 = inducedDipoleX*u.x*v.x;
                        float term02 = inducedDipolePolarY*u.y*v.x + inducedDipolePolarZ*u.x*v.y;
                        float term12 = inducedDipolePolarX*u.x*v.x;
                        result.x += term01*t.x + term11*t.y;
                        result.y += term02*t.x + term12*t.y;
                    }
                }
            }
        }
        cSim.pPmeGrid[gridIndex] = result;
    }
}

__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(768, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(384, 1)
#else
__launch_bounds__(192, 1)
#endif
void kAmoebaReciprocalConvolution_kernel()
{
    const unsigned int gridSize = cSim.pmeGridSize.x*cSim.pmeGridSize.y*cSim.pmeGridSize.z;
    float expFactor = LOCAL_HACK_PI*LOCAL_HACK_PI/(cSim.alphaEwald*cSim.alphaEwald);
    float scaleFactor = 1.0/(LOCAL_HACK_PI*cSim.periodicBoxSizeX*cSim.periodicBoxSizeY*cSim.periodicBoxSizeZ);
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < gridSize; index += blockDim.x*gridDim.x)
    {
        int kx = index/(cSim.pmeGridSize.y*cSim.pmeGridSize.z);
        int remainder = index-kx*cSim.pmeGridSize.y*cSim.pmeGridSize.z;
        int ky = remainder/cSim.pmeGridSize.z;
        int kz = remainder-ky*cSim.pmeGridSize.z;
        if (kx == 0 && ky == 0 && kz == 0)
            continue;
        int mx = (kx < (cSim.pmeGridSize.x+1)/2) ? kx : (kx-cSim.pmeGridSize.x);
        int my = (ky < (cSim.pmeGridSize.y+1)/2) ? ky : (ky-cSim.pmeGridSize.y);
        int mz = (kz < (cSim.pmeGridSize.z+1)/2) ? kz : (kz-cSim.pmeGridSize.z);
        float mhx = mx*cSim.invPeriodicBoxSizeX;
        float mhy = my*cSim.invPeriodicBoxSizeY;
        float mhz = mz*cSim.invPeriodicBoxSizeZ;
        float bx = cSim.pPmeBsplineModuli[0][kx];
        float by = cSim.pPmeBsplineModuli[1][ky];
        float bz = cSim.pPmeBsplineModuli[2][kz];
        cuComplex grid = cSim.pPmeGrid[index];
        float m2 = mhx*mhx+mhy*mhy+mhz*mhz;
        float denom = m2*bx*by*bz;
        float eterm = scaleFactor*exp(-expFactor*m2)/denom;
        cSim.pPmeGrid[index] = make_cuComplex(grid.x*eterm, grid.y*eterm);
    }
}

__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(384, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(192, 1)
#else
__launch_bounds__(96, 1)
#endif
void kComputeFixedPotentialFromGrid_kernel()
{
    // extract the permanent multipole field at each site

    for (int m = blockIdx.x*blockDim.x+threadIdx.x; m < cSim.atoms; m += blockDim.x*gridDim.x) {
        int4 gridPoint = cAmoebaSim.pIgrid[m];
        float tuv000 = 0.0f;
        float tuv001 = 0.0f;
        float tuv010 = 0.0f;
        float tuv100 = 0.0f;
        float tuv200 = 0.0f;
        float tuv020 = 0.0f;
        float tuv002 = 0.0f;
        float tuv110 = 0.0f;
        float tuv101 = 0.0f;
        float tuv011 = 0.0f;
        float tuv300 = 0.0f;
        float tuv030 = 0.0f;
        float tuv003 = 0.0f;
        float tuv210 = 0.0f;
        float tuv201 = 0.0f;
        float tuv120 = 0.0f;
        float tuv021 = 0.0f;
        float tuv102 = 0.0f;
        float tuv012 = 0.0f;
        float tuv111 = 0.0f;
        for (int iz = 0; iz < AMOEBA_PME_ORDER; iz++) {
            int k = gridPoint.z+iz-(gridPoint.z+iz >= cSim.pmeGridSize.z ? cSim.pmeGridSize.z : 0);
            float4 v = cAmoebaSim.pThetai3[m*AMOEBA_PME_ORDER+iz];
            float tu00 = 0.0f;
            float tu10 = 0.0f;
            float tu01 = 0.0f;
            float tu20 = 0.0f;
            float tu11 = 0.0f;
            float tu02 = 0.0f;
            float tu30 = 0.0f;
            float tu21 = 0.0f;
            float tu12 = 0.0f;
            float tu03 = 0.0f;
            for (int iy = 0; iy < AMOEBA_PME_ORDER; iy++) {
                int j = gridPoint.y+iy-(gridPoint.y+iy >= cSim.pmeGridSize.y ? cSim.pmeGridSize.y : 0);
                float4 u = cAmoebaSim.pThetai2[m*AMOEBA_PME_ORDER+iy];
                float4 t = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
                for (int ix = 0; ix < AMOEBA_PME_ORDER; ix++) {
                    int i = gridPoint.x+ix-(gridPoint.x+ix >= cSim.pmeGridSize.x ? cSim.pmeGridSize.x : 0);
                    int gridIndex = i*cSim.pmeGridSize.y*cSim.pmeGridSize.z + j*cSim.pmeGridSize.z + k;
                    float tq = cSim.pPmeGrid[gridIndex].x;
                    float4 tadd = cAmoebaSim.pThetai1[m*AMOEBA_PME_ORDER+ix];
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
        cAmoebaSim.pPhi[20*m] = tuv000;
        cAmoebaSim.pPhi[20*m+1] = tuv100;
        cAmoebaSim.pPhi[20*m+2] = tuv010;
        cAmoebaSim.pPhi[20*m+3] = tuv001;
        cAmoebaSim.pPhi[20*m+4] = tuv200;
        cAmoebaSim.pPhi[20*m+5] = tuv020;
        cAmoebaSim.pPhi[20*m+6] = tuv002;
        cAmoebaSim.pPhi[20*m+7] = tuv110;
        cAmoebaSim.pPhi[20*m+8] = tuv101;
        cAmoebaSim.pPhi[20*m+9] = tuv011;
        cAmoebaSim.pPhi[20*m+10] = tuv300;
        cAmoebaSim.pPhi[20*m+11] = tuv030;
        cAmoebaSim.pPhi[20*m+12] = tuv003;
        cAmoebaSim.pPhi[20*m+13] = tuv210;
        cAmoebaSim.pPhi[20*m+14] = tuv201;
        cAmoebaSim.pPhi[20*m+15] = tuv120;
        cAmoebaSim.pPhi[20*m+16] = tuv021;
        cAmoebaSim.pPhi[20*m+17] = tuv102;
        cAmoebaSim.pPhi[20*m+18] = tuv012;
        cAmoebaSim.pPhi[20*m+19] = tuv111;
    }
}

__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(256, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(128, 1)
#else
__launch_bounds__(64, 1)
#endif
void kComputeInducedPotentialFromGrid_kernel()
{
    // extract the induced dipole field at each site

    for (int m = blockIdx.x*blockDim.x+threadIdx.x; m < cSim.atoms; m += blockDim.x*gridDim.x) {
        int4 gridPoint = cAmoebaSim.pIgrid[m];
        float tuv100_1 = 0.0f;
        float tuv010_1 = 0.0f;
        float tuv001_1 = 0.0f;
        float tuv200_1 = 0.0f;
        float tuv020_1 = 0.0f;
        float tuv002_1 = 0.0f;
        float tuv110_1 = 0.0f;
        float tuv101_1 = 0.0f;
        float tuv011_1 = 0.0f;
        float tuv100_2 = 0.0f;
        float tuv010_2 = 0.0f;
        float tuv001_2 = 0.0f;
        float tuv200_2 = 0.0f;
        float tuv020_2 = 0.0f;
        float tuv002_2 = 0.0f;
        float tuv110_2 = 0.0f;
        float tuv101_2 = 0.0f;
        float tuv011_2 = 0.0f;
        float tuv000 = 0.0f;
        float tuv001 = 0.0f;
        float tuv010 = 0.0f;
        float tuv100 = 0.0f;
        float tuv200 = 0.0f;
        float tuv020 = 0.0f;
        float tuv002 = 0.0f;
        float tuv110 = 0.0f;
        float tuv101 = 0.0f;
        float tuv011 = 0.0f;
        float tuv300 = 0.0f;
        float tuv030 = 0.0f;
        float tuv003 = 0.0f;
        float tuv210 = 0.0f;
        float tuv201 = 0.0f;
        float tuv120 = 0.0f;
        float tuv021 = 0.0f;
        float tuv102 = 0.0f;
        float tuv012 = 0.0f;
        float tuv111 = 0.0f;
        for (int iz = 0; iz < AMOEBA_PME_ORDER; iz++) {
            int k = gridPoint.z+iz-(gridPoint.z+iz >= cSim.pmeGridSize.z ? cSim.pmeGridSize.z : 0);
            float4 v = cAmoebaSim.pThetai3[m*AMOEBA_PME_ORDER+iz];
            float tu00_1 = 0.0f;
            float tu01_1 = 0.0f;
            float tu10_1 = 0.0f;
            float tu20_1 = 0.0f;
            float tu11_1 = 0.0f;
            float tu02_1 = 0.0f;
            float tu00_2 = 0.0f;
            float tu01_2 = 0.0f;
            float tu10_2 = 0.0f;
            float tu20_2 = 0.0f;
            float tu11_2 = 0.0f;
            float tu02_2 = 0.0f;
            float tu00 = 0.0f;
            float tu10 = 0.0f;
            float tu01 = 0.0f;
            float tu20 = 0.0f;
            float tu11 = 0.0f;
            float tu02 = 0.0f;
            float tu30 = 0.0f;
            float tu21 = 0.0f;
            float tu12 = 0.0f;
            float tu03 = 0.0f;
            for (int iy = 0; iy < AMOEBA_PME_ORDER; iy++) {
                int j = gridPoint.y+iy-(gridPoint.y+iy >= cSim.pmeGridSize.y ? cSim.pmeGridSize.y : 0);
                float4 u = cAmoebaSim.pThetai2[m*AMOEBA_PME_ORDER+iy];
                float t0_1 = 0.0f;
                float t1_1 = 0.0f;
                float t2_1 = 0.0f;
                float t0_2 = 0.0f;
                float t1_2 = 0.0f;
                float t2_2 = 0.0f;
                float t3 = 0.0f;
                for (int ix = 0; ix < AMOEBA_PME_ORDER; ix++) {
                    int i = gridPoint.x+ix-(gridPoint.x+ix >= cSim.pmeGridSize.x ? cSim.pmeGridSize.x : 0);
                    int gridIndex = i*cSim.pmeGridSize.y*cSim.pmeGridSize.z + j*cSim.pmeGridSize.z + k;
                    cufftComplex tq = cSim.pPmeGrid[gridIndex];
                    float4 tadd = cAmoebaSim.pThetai1[m*AMOEBA_PME_ORDER+ix];
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
                float t0 = t0_1 + t0_2;
                float t1 = t1_1 + t1_2;
                float t2 = t2_1 + t2_2;
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
        cAmoebaSim.pPhid[10*m+1] = tuv100_1;
        cAmoebaSim.pPhid[10*m+2] = tuv010_1;
        cAmoebaSim.pPhid[10*m+3] = tuv001_1;
        cAmoebaSim.pPhid[10*m+4] = tuv100_1;
        cAmoebaSim.pPhid[10*m+5] = tuv010_1;
        cAmoebaSim.pPhid[10*m+6] = tuv002_1;
        cAmoebaSim.pPhid[10*m+7] = tuv110_1;
        cAmoebaSim.pPhid[10*m+8] = tuv101_1;
        cAmoebaSim.pPhid[10*m+9] = tuv011_1;
        cAmoebaSim.pPhip[10*m+1] = tuv100_2;
        cAmoebaSim.pPhip[10*m+2] = tuv010_2;
        cAmoebaSim.pPhip[10*m+3] = tuv001_2;
        cAmoebaSim.pPhip[10*m+4] = tuv100_2;
        cAmoebaSim.pPhip[10*m+5] = tuv010_2;
        cAmoebaSim.pPhip[10*m+6] = tuv002_2;
        cAmoebaSim.pPhip[10*m+7] = tuv110_2;
        cAmoebaSim.pPhip[10*m+8] = tuv101_2;
        cAmoebaSim.pPhip[10*m+9] = tuv011_2;
        cAmoebaSim.pPhidp[20*m] = tuv000;
        cAmoebaSim.pPhidp[20*m+1] = tuv100;
        cAmoebaSim.pPhidp[20*m+2] = tuv010;
        cAmoebaSim.pPhidp[20*m+3] = tuv001;
        cAmoebaSim.pPhidp[20*m+4] = tuv200;
        cAmoebaSim.pPhidp[20*m+5] = tuv020;
        cAmoebaSim.pPhidp[20*m+6] = tuv002;
        cAmoebaSim.pPhidp[20*m+7] = tuv110;
        cAmoebaSim.pPhidp[20*m+8] = tuv101;
        cAmoebaSim.pPhidp[20*m+9] = tuv011;
        cAmoebaSim.pPhidp[20*m+10] = tuv300;
        cAmoebaSim.pPhidp[20*m+11] = tuv030;
        cAmoebaSim.pPhidp[20*m+12] = tuv003;
        cAmoebaSim.pPhidp[20*m+13] = tuv210;
        cAmoebaSim.pPhidp[20*m+14] = tuv201;
        cAmoebaSim.pPhidp[20*m+15] = tuv120;
        cAmoebaSim.pPhidp[20*m+16] = tuv021;
        cAmoebaSim.pPhidp[20*m+17] = tuv102;
        cAmoebaSim.pPhidp[20*m+18] = tuv012;
        cAmoebaSim.pPhidp[20*m+19] = tuv111;
    }
}

__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(768, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(384, 1)
#else
__launch_bounds__(192, 1)
#endif
void kComputeFixedMultipoleForceAndEnergy_kernel()
{
    float multipole[10];
    const int deriv1[] = {1, 4, 7, 8, 10, 15, 17, 13, 14, 19};
    const int deriv2[] = {2, 7, 5, 9, 13, 11, 18, 15, 19, 16};
    const int deriv3[] = {3, 8, 9, 6, 14, 16, 12, 19, 17, 18};
    const float xscale = cSim.pmeGridSize.x*cSim.invPeriodicBoxSizeX;
    const float yscale = cSim.pmeGridSize.y*cSim.invPeriodicBoxSizeY;
    const float zscale = cSim.pmeGridSize.z*cSim.invPeriodicBoxSizeZ;
    float energy = 0.0f;
    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < cSim.atoms; i += blockDim.x*gridDim.x) {
        // Compute the torque.

        multipole[0] = cSim.pPosq[i].w;
        multipole[1] = cAmoebaSim.pLabFrameDipole[i*3];
        multipole[2] = cAmoebaSim.pLabFrameDipole[i*3+1];
        multipole[3] = cAmoebaSim.pLabFrameDipole[i*3+2];
        multipole[4] = cAmoebaSim.pLabFrameQuadrupole[i*9];
        multipole[5] = cAmoebaSim.pLabFrameQuadrupole[i*9+4];
        multipole[6] = cAmoebaSim.pLabFrameQuadrupole[i*9+8];
        multipole[7] = 2*cAmoebaSim.pLabFrameQuadrupole[i*9+1];
        multipole[8] = 2*cAmoebaSim.pLabFrameQuadrupole[i*9+2];
        multipole[9] = 2*cAmoebaSim.pLabFrameQuadrupole[i*9+5];
        float* phi = &cAmoebaSim.pPhi[20*i];
        cAmoebaSim.pTorque[3*i] = -cAmoebaSim.electric*(multipole[3]*yscale*phi[2] - multipole[2]*zscale*phi[3]
                      + 2.0f*(multipole[6]-multipole[5])*zscale*zscale*phi[9]
                      + multipole[8]*yscale*yscale*phi[7] + multipole[9]*xscale*yscale*phi[5]
                      - multipole[7]*yscale*zscale*phi[8] - multipole[9]*xscale*zscale*phi[6]);
        cAmoebaSim.pTorque[3*i+1] = -cAmoebaSim.electric*(multipole[1]*zscale*phi[3] - multipole[3]*xscale*phi[1]
                      + 2.0f*(multipole[4]-multipole[6])*zscale*zscale*phi[8]
                      + multipole[7]*zscale*zscale*phi[9] + multipole[8]*xscale*zscale*phi[6]
                      - multipole[8]*xscale*xscale*phi[4] - multipole[9]*yscale*yscale*phi[7]);
        cAmoebaSim.pTorque[3*i+2] = -cAmoebaSim.electric*(multipole[2]*xscale*phi[1] - multipole[1]*yscale*phi[2]
                      + 2.0f*(multipole[5]-multipole[4])*yscale*yscale*phi[7]
                      + multipole[7]*xscale*xscale*phi[4] + multipole[9]*yscale*zscale*phi[8]
                      - multipole[7]*xscale*yscale*phi[5] - multipole[8]*zscale*zscale*phi[9]);

        // Compute the force and energy.

        multipole[1] *= xscale;
        multipole[2] *= yscale;
        multipole[3] *= zscale;
        multipole[4] *= xscale*xscale;
        multipole[5] *= xscale*yscale;
        multipole[6] *= xscale*zscale;
        multipole[7] *= yscale*yscale;
        multipole[8] *= yscale*zscale;
        multipole[9] *= zscale*zscale;
        float4 f = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
        for (int k = 0; k < 10; k++) {
            energy += multipole[k]*phi[k];
            f.x += multipole[k]*phi[deriv1[k]];
            f.y += multipole[k]*phi[deriv2[k]];
            f.z += multipole[k]*phi[deriv3[k]];
        }
        f.x *= cAmoebaSim.electric*cSim.pmeGridSize.x*cSim.invPeriodicBoxSizeX;
        f.y *= cAmoebaSim.electric*cSim.pmeGridSize.y*cSim.invPeriodicBoxSizeY;
        f.z *= cAmoebaSim.electric*cSim.pmeGridSize.z*cSim.invPeriodicBoxSizeZ;
        float4 force = cSim.pForce4[i];
        force.x += f.x;
        force.y += f.y;
        force.z += f.z;
        cSim.pForce4[i] = force;


    }
    cSim.pEnergy[blockIdx.x*blockDim.x+threadIdx.x] += 0.5f*cAmoebaSim.electric*energy;
}

__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(768, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(384, 1)
#else
__launch_bounds__(192, 1)
#endif
void kComputeInducedDipoleForceAndEnergy_kernel()
{
    float multipole[10];
    float inducedDipole[3];
    float inducedDipolePolar[3];
    const int deriv1[] = {1, 4, 7, 8, 10, 15, 17, 13, 14, 19};
    const int deriv2[] = {2, 7, 5, 9, 13, 11, 18, 15, 19, 16};
    const int deriv3[] = {3, 8, 9, 6, 14, 16, 12, 19, 17, 18};
    const float xscale = cSim.pmeGridSize.x*cSim.invPeriodicBoxSizeX;
    const float yscale = cSim.pmeGridSize.y*cSim.invPeriodicBoxSizeY;
    const float zscale = cSim.pmeGridSize.z*cSim.invPeriodicBoxSizeZ;
    float energy = 0.0f;
    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < cSim.atoms; i += blockDim.x*gridDim.x) {
        // Compute the torque.

        multipole[0] = cSim.pPosq[i].w;
        multipole[1] = cAmoebaSim.pLabFrameDipole[i*3];
        multipole[2] = cAmoebaSim.pLabFrameDipole[i*3+1];
        multipole[3] = cAmoebaSim.pLabFrameDipole[i*3+2];
        multipole[4] = cAmoebaSim.pLabFrameQuadrupole[i*9];
        multipole[5] = cAmoebaSim.pLabFrameQuadrupole[i*9+4];
        multipole[6] = cAmoebaSim.pLabFrameQuadrupole[i*9+8];
        multipole[7] = 2*cAmoebaSim.pLabFrameQuadrupole[i*9+1];
        multipole[8] = 2*cAmoebaSim.pLabFrameQuadrupole[i*9+2];
        multipole[9] = 2*cAmoebaSim.pLabFrameQuadrupole[i*9+5];
        float* phidp = &cAmoebaSim.pPhidp[20*i];
        cAmoebaSim.pTorque[3*i] = -0.5f*cAmoebaSim.electric*(multipole[3]*yscale*phidp[2] - multipole[2]*zscale*phidp[3]
                      + 2.0f*(multipole[6]-multipole[5])*zscale*zscale*phidp[9]
                      + multipole[8]*yscale*yscale*phidp[7] + multipole[9]*xscale*yscale*phidp[5]
                      - multipole[7]*yscale*zscale*phidp[8] - multipole[9]*xscale*zscale*phidp[6]);
        cAmoebaSim.pTorque[3*i+1] = -0.5f*cAmoebaSim.electric*(multipole[1]*zscale*phidp[3] - multipole[3]*xscale*phidp[1]
                      + 2.0f*(multipole[4]-multipole[6])*zscale*zscale*phidp[8]
                      + multipole[7]*zscale*zscale*phidp[9] + multipole[8]*xscale*zscale*phidp[6]
                      - multipole[8]*xscale*xscale*phidp[4] - multipole[9]*yscale*yscale*phidp[7]);
        cAmoebaSim.pTorque[3*i+2] = -0.5f*cAmoebaSim.electric*(multipole[2]*xscale*phidp[1] - multipole[1]*yscale*phidp[2]
                      + 2.0f*(multipole[5]-multipole[4])*yscale*yscale*phidp[7]
                      + multipole[7]*xscale*xscale*phidp[4] + multipole[9]*yscale*zscale*phidp[8]
                      - multipole[7]*xscale*yscale*phidp[5] - multipole[8]*zscale*zscale*phidp[9]);

        // Compute the force and energy.

        multipole[1] *= xscale;
        multipole[2] *= yscale;
        multipole[3] *= zscale;
        multipole[4] *= xscale*xscale;
        multipole[5] *= xscale*yscale;
        multipole[6] *= xscale*zscale;
        multipole[7] *= yscale*yscale;
        multipole[8] *= yscale*zscale;
        multipole[9] *= zscale*zscale;
        inducedDipole[0] = cAmoebaSim.pInducedDipole[i*3];
        inducedDipole[1] = cAmoebaSim.pInducedDipole[i*3+1];
        inducedDipole[2] = cAmoebaSim.pInducedDipole[i*3+2];
        inducedDipolePolar[0] = cAmoebaSim.pInducedDipolePolar[i*3];
        inducedDipolePolar[1] = cAmoebaSim.pInducedDipolePolar[i*3+1];
        inducedDipolePolar[2] = cAmoebaSim.pInducedDipolePolar[i*3+2];
        float* phi = &cAmoebaSim.pPhi[20*i];
        float* phip = &cAmoebaSim.pPhip[10*i];
        float* phid = &cAmoebaSim.pPhid[10*i];
        float4 f = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
        for (int k = 0; k < 3; k++) {
            int j1 = deriv1[k+1];
            int j2 = deriv2[k+1];
            int j3 = deriv3[k+1];
            energy += inducedDipole[k]*phi[k+1];
            f.x += (inducedDipole[k]+inducedDipolePolar[k])*phi[j1] + inducedDipole[k]*phip[j1] + inducedDipolePolar[k]*phid[j1];
            f.y += (inducedDipole[k]+inducedDipolePolar[k])*phi[j2] + inducedDipole[k]*phip[j2] + inducedDipolePolar[k]*phid[j2];
            f.z += (inducedDipole[k]+inducedDipolePolar[k])*phi[j3] + inducedDipole[k]*phip[j3] + inducedDipolePolar[k]*phid[j3];
        }
        for (int k = 0; k < 10; k++) {
            f.x += multipole[k]*phidp[deriv1[k]];
            f.y += multipole[k]*phidp[deriv2[k]];
            f.z += multipole[k]*phidp[deriv3[k]];
        }
        f.x *= 0.5f*cAmoebaSim.electric*cSim.pmeGridSize.x*cSim.invPeriodicBoxSizeX;
        f.y *= 0.5f*cAmoebaSim.electric*cSim.pmeGridSize.y*cSim.invPeriodicBoxSizeY;
        f.z *= 0.5f*cAmoebaSim.electric*cSim.pmeGridSize.z*cSim.invPeriodicBoxSizeZ;
        float4 force = cSim.pForce4[i];
        force.x += f.x;
        force.y += f.y;
        force.z += f.z;
        cSim.pForce4[i] = force;
    }
    cSim.pEnergy[blockIdx.x*blockDim.x+threadIdx.x] += 0.5f*cAmoebaSim.electric*energy;
}

__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(768, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(384, 1)
#else
__launch_bounds__(192, 1)
#endif
void kRecordFixedMultipoleField_kernel(float* output)
{
    const float xscale = cSim.pmeGridSize.x*cSim.invPeriodicBoxSizeX;
    const float yscale = cSim.pmeGridSize.y*cSim.invPeriodicBoxSizeY;
    const float zscale = cSim.pmeGridSize.z*cSim.invPeriodicBoxSizeZ;
    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < cSim.atoms; i += blockDim.x*gridDim.x) {
        output[3*i] = -xscale*cAmoebaSim.pPhi[20*i+1];
        output[3*i+1] = -yscale*cAmoebaSim.pPhi[20*i+2];
        output[3*i+2] = -zscale*cAmoebaSim.pPhi[20*i+3];
    }
}

__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(768, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(384, 1)
#else
__launch_bounds__(192, 1)
#endif
void kRecordInducedDipoleField_kernel(float* output, float* outputPolar)
{
    const float xscale = cSim.pmeGridSize.x*cSim.invPeriodicBoxSizeX;
    const float yscale = cSim.pmeGridSize.y*cSim.invPeriodicBoxSizeY;
    const float zscale = cSim.pmeGridSize.z*cSim.invPeriodicBoxSizeZ;
    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < cSim.atoms; i += blockDim.x*gridDim.x) {
        output[3*i] -= xscale*cAmoebaSim.pPhid[10*i+1];
        output[3*i+1] -= yscale*cAmoebaSim.pPhid[10*i+2];
        output[3*i+2] -= zscale*cAmoebaSim.pPhid[10*i+3];
        outputPolar[3*i] -= xscale*cAmoebaSim.pPhip[10*i+1];
        outputPolar[3*i+1] -= yscale*cAmoebaSim.pPhip[10*i+2];
        outputPolar[3*i+2] -= zscale*cAmoebaSim.pPhip[10*i+3];
    }
}

extern void cudaComputeAmoebaMapTorquesAndAddTotalForce2(amoebaGpuContext gpu, CUDAStream<float>* psTorque, CUDAStream<float4>* psOutputForce);

/**
 * Compute the potential and forces due to the reciprocal space PME calculation for fixed multipoles.
 */
void kCalculateAmoebaPMEFixedMultipoles(amoebaGpuContext amoebaGpu)
{
    // Compute B-spline coefficients and sort the atoms.

    gpuContext gpu = amoebaGpu->gpuContext;
    int bsplineThreads = (gpu->sm_version >= SM_20 ? 448 : (gpu->sm_version >= SM_12 ? 160 : 160));
    kComputeAmoebaBsplines_kernel<<<gpu->sim.blocks, bsplineThreads, bsplineThreads*AMOEBA_PME_ORDER*AMOEBA_PME_ORDER*sizeof(float)>>>();
    LAUNCHERROR("kComputeAmoebaBsplines");
    bbSort(gpu->psPmeAtomGridIndex->_pDevData, gpu->natoms);
    kFindAmoebaAtomRangeForGrid_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block>>>();
    LAUNCHERROR("kFindAmoebaAtomRangeForGrid");

    // Perform PME for the fixed multipoles.

    kGridSpreadFixedMultipoles_kernel<<<10*gpu->sim.blocks, 64>>>();
    LAUNCHERROR("kGridSpreadFixedMultipoles");
    cufftExecC2C(gpu->fftplan, gpu->psPmeGrid->_pDevData, gpu->psPmeGrid->_pDevData, CUFFT_FORWARD);
    kAmoebaReciprocalConvolution_kernel<<<gpu->sim.blocks, gpu->sim.nonbond_threads_per_block>>>();
    LAUNCHERROR("kAmoebaReciprocalConvolution");
    cufftExecC2C(gpu->fftplan, gpu->psPmeGrid->_pDevData, gpu->psPmeGrid->_pDevData, CUFFT_INVERSE);
    int potentialThreads = (gpu->sm_version >= SM_20 ? 384 : (gpu->sm_version >= SM_12 ? 192 : 96));
    kComputeFixedPotentialFromGrid_kernel<<<gpu->sim.blocks, potentialThreads>>>();
    LAUNCHERROR("kComputeFixedPotentialFromGrid");
    kRecordFixedMultipoleField_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block>>>(amoebaGpu->psE_Field->_pDevData);
    LAUNCHERROR("kRecordFixedMultipoleField");
    kComputeFixedMultipoleForceAndEnergy_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block>>>();
    LAUNCHERROR("kComputeFixedMultipoleForceAndEnergy");
    cudaComputeAmoebaMapTorquesAndAddTotalForce2(amoebaGpu, amoebaGpu->psTorque, gpu->psForce4);
}

/**
 * Compute the potential due to the reciprocal space PME calculation for induced dipoles.
 */
void kCalculateAmoebaPMEInducedDipoleField(amoebaGpuContext amoebaGpu)
{
    // Perform PME for the induced dipoles.

    gpuContext gpu = amoebaGpu->gpuContext;
    kGridSpreadInducedDipoles_kernel<<<10*gpu->sim.blocks, 64>>>();
    LAUNCHERROR("kGridSpreadInducedDipoles");
    cufftExecC2C(gpu->fftplan, gpu->psPmeGrid->_pDevData, gpu->psPmeGrid->_pDevData, CUFFT_FORWARD);
    kAmoebaReciprocalConvolution_kernel<<<gpu->sim.blocks, gpu->sim.nonbond_threads_per_block>>>();
    LAUNCHERROR("kAmoebaReciprocalConvolution");
    cufftExecC2C(gpu->fftplan, gpu->psPmeGrid->_pDevData, gpu->psPmeGrid->_pDevData, CUFFT_INVERSE);
    int potentialThreads = (gpu->sm_version >= SM_20 ? 256 : (gpu->sm_version >= SM_12 ? 128 : 64));
    kComputeInducedPotentialFromGrid_kernel<<<gpu->sim.blocks, potentialThreads>>>();
    LAUNCHERROR("kComputeInducedPotentialFromGrid");
    kRecordInducedDipoleField_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block>>>(amoebaGpu->psWorkVector[0]->_pDevData, amoebaGpu->psWorkVector[1]->_pDevData);
    LAUNCHERROR("kRecordInducedDipoleField");
}

/**
 * Compute the forces due to the reciprocal space PME calculation for induced dipoles.
 */
void kCalculateAmoebaPMEInducedDipoleForces(amoebaGpuContext amoebaGpu)
{
    // Perform PME for the induced dipoles.

    gpuContext gpu = amoebaGpu->gpuContext;
    kGridSpreadInducedDipoles_kernel<<<10*gpu->sim.blocks, 64>>>();
    LAUNCHERROR("kGridSpreadInducedDipoles");
    cufftExecC2C(gpu->fftplan, gpu->psPmeGrid->_pDevData, gpu->psPmeGrid->_pDevData, CUFFT_FORWARD);
    kAmoebaReciprocalConvolution_kernel<<<gpu->sim.blocks, gpu->sim.nonbond_threads_per_block>>>();
    LAUNCHERROR("kAmoebaReciprocalConvolution");
    cufftExecC2C(gpu->fftplan, gpu->psPmeGrid->_pDevData, gpu->psPmeGrid->_pDevData, CUFFT_INVERSE);
    int potentialThreads = (gpu->sm_version >= SM_20 ? 256 : (gpu->sm_version >= SM_12 ? 128 : 64));
    kComputeInducedPotentialFromGrid_kernel<<<gpu->sim.blocks, potentialThreads>>>();
    LAUNCHERROR("kComputeInducedPotentialFromGrid");
    kComputeInducedDipoleForceAndEnergy_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block>>>();
    LAUNCHERROR("kComputeInducedDipoleForceAndEnergy");
    cudaComputeAmoebaMapTorquesAndAddTotalForce2(amoebaGpu, amoebaGpu->psTorque, gpu->psForce4);
}
