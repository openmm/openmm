//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------

#include "amoebaGpuTypes.h"
#include "amoebaCudaKernels.h"
#include "kCalculateAmoebaCudaUtilities.h"
#include "kCalculateAmoebaCudaWcaDispersionParticle.h"
#include "amoebaScaleFactors.h"

#include <stdio.h>

using namespace std;

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;

void SetCalculateAmoebaCudaWcaDispersionSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaWcaDispersionSim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaWcaDispersionFieldSim: cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

void GetCalculateAmoebaCudaWcaDispersionSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaWcaDispersionSim: cudaMemcpyFromSymbol: SetSim copy from cSim failed");
    status = cudaMemcpyFromSymbol(&amoebaGpu->amoebaSim, cAmoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaWcaDispersionSim: cudaMemcpyFromSymbol: SetSim copy from cAmoebaSim failed");
}

//#define AMOEBA_DEBUG
#undef AMOEBA_DEBUG

__device__ void zeroWcaDispersionSharedForce( struct WcaDispersionParticle* sA ) 
{
    // zero shared fields

    sA->force[0]              = 0.0f;
    sA->force[1]              = 0.0f;
    sA->force[2]              = 0.0f;

}

__device__ void loadWcaDispersionShared( struct WcaDispersionParticle* sA, unsigned int atomI,
                                         float4* atomCoord, float2* wcaParameters )
{
    // coordinates, radius, epsilon

    sA->x                        = atomCoord[atomI].x;
    sA->y                        = atomCoord[atomI].y;
    sA->z                        = atomCoord[atomI].z;

    sA->radius                   = wcaParameters[atomI].x;
    sA->epsilon                  = wcaParameters[atomI].y;

}

// load struct and arrays w/ shared data in sA

__device__ void loadWcaDispersionData( struct WcaDispersionParticle* sA,
                                       float4* jCoord, float* jRadius, float* jEpsilon )
{

    // load coordinates, radius, epsilon

    jCoord->x               = sA->x;
    jCoord->y               = sA->y;
    jCoord->z               = sA->z;

    *jRadius                = sA->radius;
    *jEpsilon               = sA->epsilon;
}

__device__ void calculateWcaDispersionInit_kernel( float iRadius,   float iEpsilon,
                                                   float* rmixo,    float* rmixh,
                                                   float* emixo,    float* emixh

 )
{

    float sqrtEps                   = sqrtf(iEpsilon);
    float denominator               = sqrtf(cAmoebaSim.epso) + sqrtEps;
         *emixo                     = 4.0f*cAmoebaSim.epso*iEpsilon / (denominator*denominator);

          denominator               = sqrtf(cAmoebaSim.epsh) + sqrtEps;
         *emixh                     = 4.0f*cAmoebaSim.epsh*iEpsilon / (denominator*denominator);

    float iRadius2                  = iRadius*iRadius;
    float rmino2                    = cAmoebaSim.rmino*cAmoebaSim.rmino; 
         *rmixo                     = 2.0f*(rmino2*cAmoebaSim.rmino + iRadius2*iRadius) / (rmino2 + iRadius2);

    float rminh2                    = cAmoebaSim.rminh*cAmoebaSim.rminh;
         *rmixh                     = 2.0f*(rminh2*cAmoebaSim.rminh + iRadius2*iRadius) / (rminh2+iRadius2);

}

__device__ void calculateWcaDispersionPairIxn_kernel( float4 atomCoordinatesI, float4 atomCoordinatesJ,
                                                      float radiusI,  float radiusJ,
                                                      float rmixo,    float rmixh,
                                                      float emixo,    float emixh,
                                                      float force[3], float* energy

#ifdef AMOEBA_DEBUG
               , float4* debugArray
#endif

 )
{

    const float pi         = 3.1415926535897f;
    const float shctd      = cAmoebaSim.shctd;
    const float awater     = cAmoebaSim.awater;

    // ---------------------------------------------------------------------------------------
    
    // get deltaR, and r between 2 atoms
    
    force[0]                                     = atomCoordinatesJ.x - atomCoordinatesI.x;
    force[1]                                     = atomCoordinatesJ.y - atomCoordinatesI.y;
    force[2]                                     = atomCoordinatesJ.z - atomCoordinatesI.z;

    float r2                                     = force[0]*force[0] + force[1]*force[1] + force[2]*force[2];
    if( r2 <= 0.0f ){
        force[0] = force[1] = force[2] = *energy = 0.0f;
        return;
    }
    float rI                                     = rsqrtf( r2 );
    float r                                      = 1.0f/rI;

    float sk                                     = radiusJ*shctd;
    float sk2                                    = sk*sk;
    if( radiusI >= (r+sk) ){
        force[0] = force[1] = force[2] = *energy = 0.0f;
        return;
    }

    float rmax                                   = radiusI > (r - sk) ? radiusI : (r - sk);
    float lik                                    = rmax;
    float lik2                                   = lik*lik;
    float lik3                                   = lik2*lik;
    float lik4                                   = lik2*lik2;
 
    float uik                                    = (r+sk) < rmixo ? (r+sk) : rmixo;
    float uik2                                   = uik*uik;
    float uik3                                   = uik2*uik;
    float uik4                                   = uik2*uik2;

    // 3453
    float term                                   = 4.0f*pi/(48.f*r)*(3.0f*(lik4-uik4) - 8.0f*r*(lik3-uik3) + 6.0f*(r2-sk2)*(lik2-uik2));

    float r3                                     = r2*r;
    float dl1                                    = lik2*(-lik2 + 2.0f*(r2 + sk2) );
    float dl2                                    = lik*(-lik3 + 4.0f*lik2*r - 6.0f*lik*r2 + 2.0f*lik*sk2 + 4.0f*r3 - 4.0f*r*sk2);
    float dl                                     = radiusI > (r-sk)? dl1 : dl2;

    // 3464

    float du1                                    = uik2*(-uik2 + 2.0f*(r2 + sk2) );
    float du2                                    = uik*(-uik3 + 4.0f*uik2*r - 2.0f*uik*(3.0f*r2 - sk2) + 4.0f*r*(r2 - sk2) );
    //float du2                                    = uik*(uik*( -uik2 + 4.0f*uik*r - 2.0f*(3.0f*r2 - sk2)) + 4.0f*r*(r2 - sk2) );
    float du                                     = (r+sk) > rmixo ? du1 : du2;
          du                                    *= -1.0f;

    float mask2                                  = lik < rmixo ? 1.0f : 0.0f;
    float sum                                    = -mask2*(emixo*term);
    float de                                     = -mask2*emixo*pi*(dl+du)/(4.0f*r2);

    // block at 3476

    uik                                          = (r+sk) < rmixh ? (r+sk) : rmixh;
    uik2                                         = uik*uik;
    uik3                                         = uik2*uik;
    uik4                                         = uik2*uik2;

    // 3481

    term                                         = (pi)/ (12.0f*r) * (3.0f*(lik4-uik4) - 8.0f*r*(lik3-uik3) + 6.0f*(r2-sk2)*(lik2-uik2));

    dl1                                          = lik2*(-lik2 + 2.0f*r2 + 2.0f*sk2);
    dl2                                          = lik*(-lik3 + 4.0f*lik2*r - 6.0f*lik*r2 + 2.0f*lik*sk2 + 4.0f*r3 - 4.0f*r*sk2);
    dl                                           = radiusI > (r-sk) ? dl1 : dl2;

    // 3492

    du1                                          = -uik2*(-uik2 + 2.0f*r2 + 2.0f*sk2);
    du2                                          = -uik*(-uik3 + 4.0f*uik2*r - 6.0f*uik*r2 + 2.0f*uik*sk2 + 4.0f*r3 - 4.0f*r*sk2);
    du                                           = (r+sk) > rmixh ? du1 : du2;

    mask2                                        = lik < rmixh ? 1.0f : 0.0f;
    sum                                         -= mask2*(2.0f*emixh*term);
    de                                          -= mask2*(2.0f*emixh*pi*(dl+du)/(4.0f*r2));

    // 3504

    uik                                          = r + sk;
    uik2                                         = uik*uik;
    uik3                                         = uik2*uik;
    uik4                                         = uik2*uik2;
    float uik5                                   = uik4*uik;
    float uik6                                   = uik3*uik3;
    float uik10                                  = uik5*uik5;
    float uik11                                  = uik10*uik;
    float uik12                                  = uik6*uik6;
    float uik13                                  = uik12*uik;

    lik                                          = rmax > rmixo ? rmax : rmixo;
    lik2                                         = lik*lik;
    lik3                                         = lik2*lik;
    lik4                                         = lik2*lik2;
    float lik5                                   = lik4*lik;
    float lik6                                   = lik3*lik3;
    float lik10                                  = lik5*lik5;
    float lik11                                  = lik10*lik;
    float lik12                                  = lik6*lik6;
    float lik13                                  = lik12*lik;

    // 3525

    term                                         = 4.0f*pi/(120.0f*r*lik5*uik5)*(15.0f*uik*lik*r*(uik4-lik4) - 10.0f*uik2*lik2*(uik3-lik3) + 6.0f*(sk2-r2)*(uik5-lik5));
    dl1                                          = (-5.0f*lik2 + 3.0f*r2 + 3.0f*sk2)/lik5;
    dl2                                          = ( 5.0f*lik3 - 33.0f*lik*r2 - 3.0f*lik*sk2 + 15.0f*(lik2*r+r3-r*sk2))/lik6;
    dl                                           = (radiusI > (r-sk)) || (rmax < rmixo) ? -dl1 : dl2;

    du                                           = (-5.0f*uik3 + 33.0f*uik*r2 + 3.0f*uik*sk2 - 15.0f*(uik2*r+r3-r*sk2))/uik6;

    float rmixo7                                 = rmixo*rmixo*rmixo;
          rmixo7                                 = rmixo7*rmixo7*rmixo;
    float ao                                     = emixo*rmixo7;

    // 3540

    float idisp                                  = -2.0f*ao*term;
    mask2                                        = uik > rmixo ? 1.0f : 0.0f;

    // 3541
    de                                          -= mask2*(2.0f*ao*pi*(dl + du)/(15.0f*r2));

    // 3542

    term                                         = 4.0f*pi/(2640.0f*r*lik12*uik12) * (120.0f*uik*lik*r*(uik11-lik11) - 66.0f*uik2*lik2*(uik10-lik10) + 55.0f*(sk2-r2)*(uik12-lik12));

    // 3546

    dl1                                          = (6.0f*lik2 - 5.0f*r2 - 5.0f*sk2)/lik12;
    dl2                                          = (6.0f*lik3 - 125.0f*lik*r2 - 5.0f*lik*sk2 + 60.0f*(lik2*r+r3-r*sk2))/lik13;
    dl                                           = (radiusI > (r-sk)) || (rmax < rmixo) ? dl1 : dl2;

    // 3554

    du                                           = (-6.0f*uik3 + 125.0f*uik*r2 + 5.0f*uik*sk2 - 60.0f*(uik2*r+r3-r*sk2))/uik13;

    de                                          += mask2*(ao*rmixo7*pi*(dl + du)/(60.0f*r2));
    float irep                                   = ao*rmixo7*term;
    sum                                         += mask2*(irep + idisp);

    // 3562

    lik                                          = rmax > rmixh ? rmax : rmixh;
    lik2                                         = lik*lik;
    lik3                                         = lik2*lik;
    lik4                                         = lik2*lik2;
    lik5                                         = lik4*lik;
    lik6                                         = lik3*lik3;
    lik10                                        = lik5*lik5;
    lik11                                        = lik10*lik;
    lik12                                        = lik6*lik6;
    lik13                                        = lik12*lik;

    // 3572

    term                                         = 4.0f * pi / (120.0f*r*lik5*uik5) * (15.0f*uik*lik*r*(uik4-lik4) -
                                                   10.0f*uik2*lik2*(uik3-lik3) + 6.0f*(sk2-r2)*(uik5-lik5));

    dl1                                          = (-5.0f*lik2 + 3.0f*r2 + 3.0f*sk2)/lik5;
    dl2                                          = (5.0f*lik3 - 33.0f*lik*r2 - 3.0f*lik*sk2+ 15.0f*(lik2*r+r3-r*sk2))/lik6;
    dl                                           = (radiusI > (r-sk)) || (rmax < rmixh) ? -dl1 : dl2;

    du                                           = -(5.0f*uik3 - 33.0f*uik*r2 - 3.0f*uik*sk2 + 15.0f*(uik2*r+r3-r*sk2))/uik6;

    float rmixh7                                 = rmixh*rmixh*rmixh;
          rmixh7                                 = rmixh7*rmixh7*rmixh;
    float ah                                     = emixh * rmixh7;

    // 3587
    idisp                                        = -4.0f * ah * term;

    mask2                                        = uik > rmixh ? 1.0f : 0.0f;
    de                                          -= mask2*(4.0f*ah*pi*(dl + du)/(15.0f*r2));

    term                                         = 4.0f * pi / (2640.0f*r*lik12*uik12) * (120.0f*uik*lik*r*(uik11-lik11) -
                                                   66.0f*uik2*lik2*(uik10-lik10) + 55.0f*(sk2-r2)*(uik12-lik12));

    // 3593

    dl1                                          = -(-6.0f*lik2 + 5.0f*r2 + 5.0f*sk2)/lik12;
    dl2                                          =  (6.0f*lik3 - 125.0f*lik*r2 - 5.0f*lik*sk2 + 60.0f*(lik2*r+r3-r*sk2))/lik13;
    dl                                           = ( (radiusI > (r-sk) ) || (rmax < rmixh) ) ? dl1 : dl2;

    // 3603

    du                                           = -(6.0f*uik3 - 125.0f*uik*r2 -5.0f*uik*sk2 + 60.0f*(uik2*r+r3-r*sk2))/uik13;
    irep                                         = 2.0f*ah*rmixh7*term;

    de                                          += mask2*(ah*rmixh7*pi*(dl+du)/(30.0f*r2));
    sum                                         += mask2*(irep + idisp);

    *energy                                      = sum;

    de                                          *= -(awater/r);
    force[0]                                    *= de;
    force[1]                                    *= de;
    force[2]                                    *= de;

#ifdef AMOEBA_DEBUG
    debugArray[0].x                              = sum;
    debugArray[0].y                              = sum;
    debugArray[0].z                              = sum;
    debugArray[0].w                              = sum;
#if 0
    debugArray[0].x                              = r;
    debugArray[0].y                              = -r*de/awater;
    debugArray[0].z                              = emixo;
    debugArray[0].w                              = mask2;

    debugArray[1].x                              = dl;
    debugArray[1].y                              = du;
    debugArray[1].z                              = lik;
    debugArray[1].w                              = uik;

    debugArray[2].x                              = du1;
    debugArray[2].y                              = du2;
    debugArray[2].z                              = term;
    debugArray[2].w                              = sk;
#endif

#endif
}

// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateAmoebaCudaWcaDispersion.h"

#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateAmoebaCudaWcaDispersion.h"

// reduce psWorkArray_3_1 -> outputArray

static void kReduceWcaDispersion(amoebaGpuContext amoebaGpu, CUDAStream<float>* outputArray )
{
    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                               amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                               amoebaGpu->psWorkArray_3_1->_pDevData, outputArray->_pDevData );
    LAUNCHERROR("kReduceWcaDispersion");
}

// reduce psWorkArray_3_1 -> outputArray

static void kReduceWcaDispersionToFloat4(amoebaGpuContext amoebaGpu, CUDAStream<float4>* outputArray )
{
    kReduceFieldsToFloat4_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                                   amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                                   amoebaGpu->psWorkArray_3_1->_pDevData, outputArray->_pDevData );
    LAUNCHERROR("kReduceWcaDispersion");
}

/**---------------------------------------------------------------------------------------

   Compute WCA dispersion

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

void kCalculateAmoebaWcaDispersionForces( amoebaGpuContext amoebaGpu )
{
  
   // ---------------------------------------------------------------------------------------

    static int threadsPerBlock    = 0;

   // ---------------------------------------------------------------------------------------

     gpuContext gpu    = amoebaGpu->gpuContext;

#ifdef AMOEBA_DEBUG
    static const char* methodName = "kCalculateAmoebaWcaDispersionForces";
    int paddedNumberOfAtoms                    = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
    CUDAStream<float4>* debugArray             = new CUDAStream<float4>(paddedNumberOfAtoms*paddedNumberOfAtoms, 1, "DebugArray");
    memset( debugArray->_pSysData,      0, sizeof( float )*4*paddedNumberOfAtoms*paddedNumberOfAtoms);
    debugArray->Upload();
    int targetAtom                             = 3;
#endif

    // set threads/block first time through

    if( threadsPerBlock == 0 ){
        unsigned int maxThreads;
        if (gpu->sm_version >= SM_20)
            maxThreads = 384;
        else if (gpu->sm_version >= SM_12)
            maxThreads = 192;
        else
            maxThreads = 64;
       threadsPerBlock = std::min(getThreadsPerBlock( amoebaGpu, sizeof(WcaDispersionParticle)), maxThreads);
    }

    if (gpu->bOutputBufferPerWarp){
        kCalculateAmoebaWcaDispersionN2ByWarp_kernel<<<amoebaGpu->nonbondBlocks, threadsPerBlock, sizeof(WcaDispersionParticle)*threadsPerBlock>>>(
                                                            amoebaGpu->psWorkUnit->_pDevData,
                                                            gpu->psPosq4->_pDevData,
                                                            amoebaGpu->psWcaDispersionRadiusEpsilon->_pDevData,
#ifdef AMOEBA_DEBUG
                                                            amoebaGpu->psWorkArray_3_1->_pDevData,
                                                            debugArray->_pDevData, targetAtom );
#else
                                                            amoebaGpu->psWorkArray_3_1->_pDevData );
#endif

    } else {

#ifdef AMOEBA_DEBUG
        (void) fprintf( amoebaGpu->log, "N2 no warp\n" );
        (void) fprintf( amoebaGpu->log, "%s numBlocks=%u numThreads=%u bufferPerWarp=%u atm=%u shrd=%u ixnCt=%u workUnits=%u\n",
                        methodName, amoebaGpu->nonbondBlocks, threadsPerBlock, amoebaGpu->bOutputBufferPerWarp,
                        sizeof(WcaDispersionParticle), sizeof(WcaDispersionParticle)*threadsPerBlock,
                        (*gpu->psInteractionCount)[0], gpu->sim.workUnits );
        (void) fflush( amoebaGpu->log );
#endif

        kCalculateAmoebaWcaDispersionN2_kernel<<<amoebaGpu->nonbondBlocks, threadsPerBlock, sizeof(WcaDispersionParticle)*threadsPerBlock>>>(
                                                            amoebaGpu->psWorkUnit->_pDevData,
                                                            gpu->psPosq4->_pDevData,
                                                            amoebaGpu->psWcaDispersionRadiusEpsilon->_pDevData,
#ifdef AMOEBA_DEBUG
                                                            amoebaGpu->psWorkArray_3_1->_pDevData,
                                                            debugArray->_pDevData, targetAtom );
#else
                                                            amoebaGpu->psWorkArray_3_1->_pDevData );
#endif

    }
    LAUNCHERROR("kCalculateAmoebaWcaDispersion");  

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){

        amoebaGpu->psWorkArray_3_1->Download();
        debugArray->Download();

        (void) fprintf( amoebaGpu->log,"\n" );
        int paddedNumberOfAtoms = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
        double sum = 0.0;
        double sums[8] = { 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0 };
        std::map<int,double> buffers;
        std::vector< std::vector<double> > trackD;
        std::vector< std::vector<double> > trackT;
        std::vector<double> maxD;
        std::vector< std::vector<int> > trackI;
        trackD.resize( gpu->natoms );
        trackT.resize( gpu->natoms );
        maxD.resize( gpu->natoms );
        trackI.resize( gpu->natoms );
        for( int jj = 0; jj < gpu->natoms; jj++ ){
            unsigned int debugIndex = jj;
            (void) fprintf( amoebaGpu->log,"%5d %5d DebugWca\n", targetAtom, jj );
            int block = -1;

            for( int kk = 0; kk < -3; kk++ ){
                if( kk == 1 ){
                    block =  static_cast<int>(debugArray->_pSysData[debugIndex].w + 1.0e-04); 
                    if( buffers.find(block) == buffers.end() ){
                        buffers[block] = 0.0;
                    }
                }
                if( kk == 1  && jj != targetAtom ){
                    sums[0] += debugArray->_pSysData[debugIndex].y;
                    sums[1] += debugArray->_pSysData[debugIndex].z;
                    sums[2] += debugArray->_pSysData[debugIndex].w;
                    double x4 = debugArray->_pSysData[debugIndex].x - (debugArray->_pSysData[debugIndex].y + debugArray->_pSysData[debugIndex].z + debugArray->_pSysData[debugIndex].w);
                    sums[3] += x4;
                    //sum     += debugArray->_pSysData[debugIndex].x;
                    sum     += debugArray->_pSysData[debugIndex].z;
                    buffers[block] += debugArray->_pSysData[debugIndex].z; 
                    (void) fprintf( amoebaGpu->log," %16.9e [%16.9e %16.9e %16.9e %16.9e]\n",
                                    debugArray->_pSysData[debugIndex].x, debugArray->_pSysData[debugIndex].y,
                                    debugArray->_pSysData[debugIndex].z, debugArray->_pSysData[debugIndex].w, x4);
         
                } else if( kk == 2 &&  jj != targetAtom){
                    //sum     += debugArray->_pSysData[debugIndex].x;
                    sum     += debugArray->_pSysData[debugIndex].z;
                    sums[4] += debugArray->_pSysData[debugIndex].y;
                    sums[5] += debugArray->_pSysData[debugIndex].z;
                    sums[6] += debugArray->_pSysData[debugIndex].w;
                    double x4 = debugArray->_pSysData[debugIndex].x - (debugArray->_pSysData[debugIndex].y + debugArray->_pSysData[debugIndex].z + debugArray->_pSysData[debugIndex].w);
                    sums[7] += x4;
                    buffers[block] += debugArray->_pSysData[debugIndex].z; 
                    (void) fprintf( amoebaGpu->log," %16.9e [%16.9e %16.9e %16.9e %16.9e]\n",
                                    debugArray->_pSysData[debugIndex].x, debugArray->_pSysData[debugIndex].y,
                                    debugArray->_pSysData[debugIndex].z, debugArray->_pSysData[debugIndex].w, x4);
                } else {
                    (void) fprintf( amoebaGpu->log,"[%16.9e %16.9e %16.9e %16.9e] %7u\n",
                                    debugArray->_pSysData[debugIndex].x, debugArray->_pSysData[debugIndex].y,
                                    debugArray->_pSysData[debugIndex].z, debugArray->_pSysData[debugIndex].w, debugIndex );
                }
                if( kk == 5 )(void) fprintf( amoebaGpu->log,"\n" );
                debugIndex += paddedNumberOfAtoms;
            }

            block =  static_cast<int>(debugArray->_pSysData[debugIndex+paddedNumberOfAtoms].w + 1.0e-04); 
            if( buffers.find(block) == buffers.end() ){
                buffers[block] = 0.0;
                maxD[block]    = 0.0;
            }
            for( int kk = 0; kk < 3; kk++ ){
                if( kk == 0 && jj != targetAtom ){
                    sum            += debugArray->_pSysData[debugIndex].z;
                    buffers[block] += debugArray->_pSysData[debugIndex].z; 
                    trackI[block].push_back( jj );
                    trackD[block].push_back( debugArray->_pSysData[debugIndex].z );
                    trackT[block].push_back( debugArray->_pSysData[debugIndex].w );
                    if( fabs( debugArray->_pSysData[debugIndex].w ) > maxD[block] ){
                        maxD[block] = fabs( debugArray->_pSysData[debugIndex].w );
                    }
                    (void) fprintf( amoebaGpu->log,"[%16.9e %16.9e %16.9e %16.9e]\n",
                                    debugArray->_pSysData[debugIndex].x, debugArray->_pSysData[debugIndex].y,
                                    debugArray->_pSysData[debugIndex].z, debugArray->_pSysData[debugIndex].w);
         
                } else if( kk == 2 &&  jj != targetAtom){
                    sum             += debugArray->_pSysData[debugIndex].z;
                    buffers[block]  += debugArray->_pSysData[debugIndex].z; 
                    trackI[block].push_back( jj );
                    trackD[block].push_back( debugArray->_pSysData[debugIndex].z );
                    trackT[block].push_back( debugArray->_pSysData[debugIndex].w );
                    if( fabs( debugArray->_pSysData[debugIndex].w ) > maxD[block] ){
                        maxD[block] = fabs( debugArray->_pSysData[debugIndex].w );
                    }
                    (void) fprintf( amoebaGpu->log,"[%16.9e %16.9e %16.9e %16.9e]\n",
                                    debugArray->_pSysData[debugIndex].x, debugArray->_pSysData[debugIndex].y,
                                    debugArray->_pSysData[debugIndex].z, debugArray->_pSysData[debugIndex].w);
                } else {
                    (void) fprintf( amoebaGpu->log,"[%16.9e %16.9e %16.9e %16.9e] %7u\n",
                                    debugArray->_pSysData[debugIndex].x, debugArray->_pSysData[debugIndex].y,
                                    debugArray->_pSysData[debugIndex].z, debugArray->_pSysData[debugIndex].w, debugIndex );
                }
                if( kk == 5 )(void) fprintf( amoebaGpu->log,"\n" );
                debugIndex += paddedNumberOfAtoms;
            }
            (void) fprintf( amoebaGpu->log,"\n" );
        }
        (void) fprintf( amoebaGpu->log,"Total sum=%14.7e\n", sum );
        (void) fprintf( amoebaGpu->log,"DeWW\n" );
        sum = 0.0;
        for( int jj = 0; jj < 4; jj++ ){
            sum += sums[jj] + sums[jj+4];
            (void) fprintf( amoebaGpu->log,"[%14.7e %14.7e] %14.7e\n", sums[jj], sums[jj+4], sums[jj] + sums[jj+4] );
        }
        (void) fprintf( amoebaGpu->log,"Total sum8=%14.7e\n", sum );

        (void) fprintf( amoebaGpu->log,"Buffers\n", sum );
        sum = 0.0;
        for( std::map<int,double>::const_iterator jj = buffers.begin(); jj != buffers.end(); jj++ ){

            (void) fprintf( amoebaGpu->log,"%5d %14.7e", jj->first, jj->second );
            sum += jj->second;

            int block       = jj->first;
            double sumBlock = 0.0;

            for( unsigned int kk = 0; kk < trackI[block].size(); kk++ ){
                 sumBlock +=  trackD[block][kk];
            }
            (void) fprintf( amoebaGpu->log," %5u Total=%14.7e MaxD=%14.7e\n", trackI[block].size(), sumBlock, maxD[block]);
            for( unsigned int kk = 0; kk < trackI[block].size(); kk++ ){
                 (void) fprintf( amoebaGpu->log,"[%5d %14.7e  %14.7e] ", trackI[block][kk], trackD[block][kk], trackT[block][kk]);
                 if( ((kk+1)%3) == 0 )(void) fprintf( amoebaGpu->log,"\n" );
            }
            (void) fprintf( amoebaGpu->log,"\n\n\n" );

        }
        (void) fprintf( amoebaGpu->log,"Total buffer sum=%14.7e\n", sum );
/*
        for( int jj = 0; jj < gpu->natoms; jj++ ){
            int debugIndex   = jj;
            int debugIndex4  = debugIndex + 4*paddedNumberOfAtoms;
            int debugIndex5  = debugIndex + 5*paddedNumberOfAtoms;
            int debugIndex9  = debugIndex + 9*paddedNumberOfAtoms;
            int debugIndex10 = debugIndex + 10*paddedNumberOfAtoms;
            (void) fprintf( amoebaGpu->log,"%6d %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e \n", jj,
                                       debugArray->_pSysData[debugIndex4].x*debugArray->_pSysData[debugIndex4].x,  // r2
                                       debugArray->_pSysData[debugIndex4].y, debugArray->_pSysData[debugIndex9].y, // de
                                       debugArray->_pSysData[debugIndex5].x, debugArray->_pSysData[debugIndex10].x, // dll
                                       debugArray->_pSysData[debugIndex5].y, debugArray->_pSysData[debugIndex10].y, // duu
                                       debugArray->_pSysData[debugIndex4].z, debugArray->_pSysData[debugIndex9].z ); // emxio
        }
*/
    }
#endif

    kReduceWcaDispersionToFloat4( amoebaGpu, gpu->psForce4 );

#ifdef AMOEBA_DEBUG
    if( 0 ){
        gpu->psEnergy->Download();
        double sum = 0.0;
        for (int i = 0; i < gpu->sim.energyOutputBuffers; i++){
            sum  += gpu->psEnergy->_pSysData[i];
            if( fabsf( (*gpu->psEnergy)[i]) > 0.0 )
                (void) fprintf( amoebaGpu->log,"SumQQ %6d %14.7e QQ SUM\n", i, (*gpu->psEnergy)[i] );
        }   
        (void) fprintf( amoebaGpu->log,"%14.7e QQ SUM\n", sum );
    }
#endif

    if( 0 ){
        int paddedNumberOfAtoms             = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
        CUDAStream<float>* psTempForce      = new CUDAStream<float>(4*paddedNumberOfAtoms, 1, "psTempForce");
        kClearFloat( amoebaGpu, 4*paddedNumberOfAtoms, psTempForce );
        kReduceWcaDispersion( amoebaGpu, psTempForce );
        std::vector<int> fileId;
        //fileId.push_back( 0 );
        VectorOfDoubleVectors outputVector;
        cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,   outputVector, NULL, 1.0f );
        cudaLoadCudaFloatArray(  gpu->natoms, 3, psTempForce,    outputVector, NULL, 1.0f );
        cudaWriteVectorOfDoubleVectorsToFile( "CudaWca", fileId, outputVector );
        delete psTempForce;
        //exit(0);
     }

#ifdef AMOEBA_DEBUG
    delete debugArray;
    //exit(0);
#endif

   // ---------------------------------------------------------------------------------------
}
