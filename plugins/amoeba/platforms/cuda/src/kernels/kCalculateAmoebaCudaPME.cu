//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------

#include "amoebaGpuTypes.h"
#include "amoebaCudaKernels.h"

//#define AMOEBA_DEBUG

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;

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
//      real*8 thetai(4,AMOEBA_PME_ORDER)
//      real*8 array(AMOEBA_PME_ORDER,AMOEBA_PME_ORDER)

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
        thetai[i] = make_float4(ARRAY(AMOEBA_PME_ORDER,i), ARRAY(AMOEBA_PME_ORDER-1,i), ARRAY(AMOEBA_PME_ORDER-2,i), ARRAY(AMOEBA_PME_ORDER-3,i));
}

/**
 * Compute bspline coefficients.
 */
__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(256, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(256, 1)
#else
__launch_bounds__(128, 1)
#endif
void computeBsplines_kernel()
{
    extern __shared__ float bsplines_cache[]; // size = block_size*pme_order*pme_order
    float* array = &bsplines_cache[threadIdx.x*AMOEBA_PME_ORDER*AMOEBA_PME_ORDER];

    //  get the B-spline coefficients for each multipole site

    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < cSim.atoms; i += blockDim.x*gridDim.x) {
        float4 posq = cSim.pPosq[i];
        posq.x -= floor(posq.x*cSim.invPeriodicBoxSizeX)*cSim.periodicBoxSizeX;
        posq.y -= floor(posq.y*cSim.invPeriodicBoxSizeY)*cSim.periodicBoxSizeY;
        posq.z -= floor(posq.z*cSim.invPeriodicBoxSizeZ)*cSim.periodicBoxSizeZ;
        float3 t = make_float3((posq.x*cSim.invPeriodicBoxSizeX)*cSim.pmeGridSize.x,
                               (posq.y*cSim.invPeriodicBoxSizeY)*cSim.pmeGridSize.y,
                               (posq.z*cSim.invPeriodicBoxSizeZ)*cSim.pmeGridSize.z);

        // First axis.

        float w = posq.x*cSim.invPeriodicBoxSizeX;
        float fr = cSim.pmeGridSize.x*(w-(int)(w+0.5f)+0.5f);
        int ifr = (int) fr;
        w = fr - ifr;
        int igrid1 = ifr - AMOEBA_PME_ORDER;
        computeBSplinePoint(&cAmoebaSim.pThetai1[i*AMOEBA_PME_ORDER], w, array);

        // Second axis.

        w = posq.y*cSim.invPeriodicBoxSizeY;
        fr = cSim.pmeGridSize.y*(w-(int)(w+0.5f)+0.5f);
        ifr = (int) fr;
        w = fr - ifr;
        int igrid2 = ifr - AMOEBA_PME_ORDER;
        computeBSplinePoint(&cAmoebaSim.pThetai2[i*AMOEBA_PME_ORDER], w, array);

        // Third axis.

        w = posq.z*cSim.invPeriodicBoxSizeZ;
        fr = cSim.pmeGridSize.z*(w-(int)(w+0.5f)+0.5f);
        ifr = (int) fr;
        w = fr - ifr;
        int igrid3 = ifr - AMOEBA_PME_ORDER;
        computeBSplinePoint(&cAmoebaSim.pThetai3[i*AMOEBA_PME_ORDER], w, array);

        cAmoebaSim.pIgrid[i] = make_int4(igrid1, igrid2, igrid3, 0);
        cSim.pPmeAtomGridIndex[i] = make_int2(i, igrid1*cSim.pmeGridSize.y*cSim.pmeGridSize.z+igrid2*cSim.pmeGridSize.z+igrid3);
    }
}


__global__
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(1024, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(512, 1)
#else
__launch_bounds__(256, 1)
#endif
void kGridSpreadMultipoles_kernel()
{
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
        for (int ix = 0; ix < PME_ORDER; ++ix)
        {
            int x = gridPoint.x-ix+(gridPoint.x >= ix ? 0 : cSim.pmeGridSize.x);
            for (int iy = 0; iy < PME_ORDER; ++iy)
            {
                int y = gridPoint.y-iy+(gridPoint.y >= iy ? 0 : cSim.pmeGridSize.y);
                int z1 = gridPoint.z-PME_ORDER+1;
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
                    float atomDipoleX = cAmoebaSim.pLabFrameDipole[atomIndex+3];
                    float atomDipoleY = cAmoebaSim.pLabFrameDipole[atomIndex+3+1];
                    float atomDipoleZ = cAmoebaSim.pLabFrameDipole[atomIndex+3+2];
                    float atomQuadrupoleXX = cAmoebaSim.pLabFrameQuadrupole[atomIndex+9];
                    float atomQuadrupoleXY = 2*cAmoebaSim.pLabFrameQuadrupole[atomIndex+9+1];
                    float atomQuadrupoleXZ = 2*cAmoebaSim.pLabFrameQuadrupole[atomIndex+9+2];
                    float atomQuadrupoleYY = cAmoebaSim.pLabFrameQuadrupole[atomIndex+9+4];
                    float atomQuadrupoleYZ = 2*cAmoebaSim.pLabFrameQuadrupole[atomIndex+9+5];
                    float atomQuadrupoleZZ = cAmoebaSim.pLabFrameQuadrupole[atomIndex+9+8];
                    float4 t = cAmoebaSim.pThetai1[atomIndex+ix*cSim.atoms];
                    float4 u = cAmoebaSim.pThetai2[atomIndex+iy*cSim.atoms];
                    float4 v = cAmoebaSim.pThetai3[atomIndex+iz*cSim.atoms];
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
                        float atomDipoleX = cAmoebaSim.pLabFrameDipole[atomIndex+3];
                        float atomDipoleY = cAmoebaSim.pLabFrameDipole[atomIndex+3+1];
                        float atomDipoleZ = cAmoebaSim.pLabFrameDipole[atomIndex+3+2];
                        float atomQuadrupoleXX = cAmoebaSim.pLabFrameQuadrupole[atomIndex+9];
                        float atomQuadrupoleXY = 2*cAmoebaSim.pLabFrameQuadrupole[atomIndex+9+1];
                        float atomQuadrupoleXZ = 2*cAmoebaSim.pLabFrameQuadrupole[atomIndex+9+2];
                        float atomQuadrupoleYY = cAmoebaSim.pLabFrameQuadrupole[atomIndex+9+4];
                        float atomQuadrupoleYZ = 2*cAmoebaSim.pLabFrameQuadrupole[atomIndex+9+5];
                        float atomQuadrupoleZZ = cAmoebaSim.pLabFrameQuadrupole[atomIndex+9+8];
                        float4 t = cAmoebaSim.pThetai1[atomIndex+ix*cSim.atoms];
                        float4 u = cAmoebaSim.pThetai2[atomIndex+iy*cSim.atoms];
                        float4 v = cAmoebaSim.pThetai3[atomIndex+iz*cSim.atoms];
                        float term0 = atomCharge*u.x*v.x + atomDipoleY*u.y*v.x + atomDipoleZ*u.x*v.y + atomQuadrupoleYY*u.z*v.x + atomQuadrupoleZZ*u.x*v.z + atomQuadrupoleYZ*u.y*v.y;
                        float term1 = atomDipoleX*u.x*v.x + atomQuadrupoleXY*u.y*v.x + atomQuadrupoleXZ*u.x*v.y;
                        float term2 = atomQuadrupoleXX * u.x * v.x;
                        result += term0*t.x + term1*t.y + term2*t.z;
                    }
                }
            }
        }
        cSim.pPmeGrid[gridIndex] = make_cuComplex(result*sqrt(cSim.epsfac), 0.0f);
    }
}
