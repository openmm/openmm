
//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------

#include "amoebaCudaKernels.h"
#include "kCalculateAmoebaCudaUtilities.h"

//#define AMOEBA_DEBUG

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;

void SetCalculateAmoebaCudaPmeFixedEFieldSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaPmeFixedEFieldSim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaPmeFixedEFieldSim: cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

void GetCalculateAmoebaCudaPmeFixedEFieldSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaPmeFixedEFieldSim: cudaMemcpyFromSymbol: SetSim copy from cSim failed");
    status = cudaMemcpyFromSymbol(&amoebaGpu->amoebaSim, cAmoebaSim, sizeof(cudaAmoebaGmxSimulation));         
    RTERROR(status, "GetCalculateAmoebaCudaPmeFixedEFieldSim: cudaMemcpyFromSymbol: SetSim copy from cAmoebaSim failed");
}

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
static void kReducePmeEFieldPolar_kernel( unsigned int fieldComponents, unsigned int outputBuffers, float* EFieldReciprocal,  float* fieldIn, float* fieldOut )
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

    // Reduce field

    const float term = (4.0f/3.0f)*(cSim.alphaEwald*cSim.alphaEwald*cSim.alphaEwald)/cAmoebaSim.sqrtPi;
    //const float term = 0.0f;
    while (pos < fieldComponents)
    {   

        // self-term included here

        float totalField = EFieldReciprocal[pos] + term*cAmoebaSim.pLabFrameDipole[pos];

        float* pFt       = fieldIn + pos;
        unsigned int i   = outputBuffers;
        while (i >= 4)
        {   
            totalField += pFt[0] + pFt[fieldComponents] + pFt[2*fieldComponents] + pFt[3*fieldComponents];
            pFt        += fieldComponents*4;
            i          -= 4;
        }   

        if (i >= 2)
        {   
            totalField += pFt[0] + pFt[fieldComponents];
            pFt        += fieldComponents*2;
            i          -= 2;
        }   

        if (i > 0)
        {   
            totalField += pFt[0];
        }   

        fieldOut[pos]   = totalField;
        pos            += gridDim.x * blockDim.x;
    }   
}


__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
static void kReducePmeEField_kernel( unsigned int fieldComponents, unsigned int outputBuffers,  float* fieldIn, float* fieldOut )
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

    // Reduce field

    const float term = (4.0f/3.0f)*(cSim.alphaEwald*cSim.alphaEwald*cSim.alphaEwald)/cAmoebaSim.sqrtPi;
    //const float term = 0.0f;
    while (pos < fieldComponents)
    {   

        // self-term included here

        float totalField = term*cAmoebaSim.pLabFrameDipole[pos];

        float* pFt       = fieldIn + pos;
        unsigned int i   = outputBuffers;
        while (i >= 4)
        {   
            totalField += pFt[0] + pFt[fieldComponents] + pFt[2*fieldComponents] + pFt[3*fieldComponents];
            pFt        += fieldComponents*4;
            i          -= 4;
        }   

        if (i >= 2)
        {   
            totalField += pFt[0] + pFt[fieldComponents];
            pFt        += fieldComponents*2;
            i          -= 2;
        }   

        if (i > 0)
        {   
            totalField += pFt[0];
        }   

        fieldOut[pos]  += totalField;
        pos            += gridDim.x * blockDim.x;
    }   
}

// reduce psWorkArray_3_1 -> EField
// reduce psWorkArray_3_2 -> EFieldPolar

static void kReducePmeDirectE_Fields(amoebaGpuContext amoebaGpu )
{

    // E_FieldPolar = E_Field (reciprocal) + E_FieldPolar (direct) + self

    kReducePmeEFieldPolar_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                                   amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                                   amoebaGpu->psE_Field->_pDevStream[0], amoebaGpu->psWorkArray_3_2->_pDevStream[0], amoebaGpu->psE_FieldPolar->_pDevStream[0] );
    LAUNCHERROR("kReducePmeE_Fields1");

    // E_Field = E_Field (reciprocal) + E_Field (direct) + self

    kReducePmeEField_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                              amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                              amoebaGpu->psWorkArray_3_1->_pDevStream[0], amoebaGpu->psE_Field->_pDevStream[0] );
    LAUNCHERROR("kReducePmeE_Fields2");
}

// file includes FixedFieldParticle struct definition/load/unload struct and body kernel for fixed E-field

#undef GK
#include "kCalculateAmoebaCudaFixedFieldParticle.h"
__device__ void calculateFixedFieldRealSpacePairIxn_kernel( FixedFieldParticle& atomI, FixedFieldParticle& atomJ,
                                                            float dscale, float pscale, float fields[4][3]
#ifdef AMOEBA_DEBUG
                                                            , float4* pullBack
#endif

 ){

    // compute the real space portion of the Ewald summation
  
    float xr          = atomJ.x - atomI.x;
    float yr          = atomJ.y - atomI.y;
    float zr          = atomJ.z - atomI.z;

    // periodic boundary conditions

    xr               -= floor(xr*cSim.invPeriodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
    yr               -= floor(yr*cSim.invPeriodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
    zr               -= floor(zr*cSim.invPeriodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;

    float r2          = xr*xr + yr* yr + zr*zr;
    float r           = sqrtf(r2);

    // calculate the error function damping terms

    float ralpha      = cSim.alphaEwald*r;
    float bn[4];

    bn[0]             = erfc(ralpha)/r;
    float alsq2       = 2.0f*cSim.alphaEwald*cSim.alphaEwald;
    float alsq2n      = 1.0f/(cAmoebaSim.sqrtPi*cSim.alphaEwald);
    float exp2a       = exp(-(ralpha*ralpha));
    alsq2n           *= alsq2;
    bn[1]             = (bn[0]+alsq2n*exp2a)/r2;

    alsq2n           *= alsq2;
    bn[2]             = (3.0f*bn[1]+alsq2n*exp2a)/r2;

    alsq2n           *= alsq2;
    bn[3]             = (5.0f*bn[2]+alsq2n*exp2a)/r2;

    // compute the error function scaled and unscaled terms

    float scale3      = 1.0f;
    float scale5      = 1.0f;
    float scale7      = 1.0f;
    float damp        = atomI.damp*atomJ.damp;
    if( damp != 0.0f ){

        float ratio  = (r/damp);
              ratio  = ratio*ratio*ratio;

        float pgamma = atomI.thole < atomJ.thole ? atomI.thole : atomJ.thole;

              damp   = -pgamma*ratio;

        if( damp > -50.0f) {
            float expdamp = exp(damp);
            scale3        = 1.0f - expdamp;
            scale5        = 1.0f - expdamp*(1.0f-damp);
            scale7        = 1.0f - expdamp*(1.0f-damp+(0.6f*damp*damp));
        }
    }
    float dsc3        = dscale*scale3;
    float dsc5        = dscale*scale5;
    float dsc7        = dscale*scale7;

    float psc3        = pscale*scale3;
    float psc5        = pscale*scale5;
    float psc7        = pscale*scale7;

    float r3          = (r*r2);
    float r5          = (r3*r2);
    float r7          = (r5*r2);
    float drr3        = (1.0f-dsc3)/r3;
    float drr5        = 3.0f * (1.0f-dsc5)/r5;
    float drr7        = 15.0f * (1.0f-dsc7)/r7;

    float prr3        = (1.0f-psc3) / r3;
    float prr5        = 3.0f *(1.0f-psc5)/r5;
    float prr7        = 15.0f*(1.0f-psc7)/r7;

    float dir         = atomI.labFrameDipole_X*xr + atomI.labFrameDipole_Y*yr + atomI.labFrameDipole_Z*zr;

    float qix         = atomI.labFrameQuadrupole_XX*xr + atomI.labFrameQuadrupole_XY*yr + atomI.labFrameQuadrupole_XZ*zr;
    float qiy         = atomI.labFrameQuadrupole_XY*xr + atomI.labFrameQuadrupole_YY*yr + atomI.labFrameQuadrupole_YZ*zr;
    float qiz         = atomI.labFrameQuadrupole_XZ*xr + atomI.labFrameQuadrupole_YZ*yr + atomI.labFrameQuadrupole_ZZ*zr;

    float qir         = qix*xr + qiy*yr + qiz*zr;

    float dkr         = atomJ.labFrameDipole_X*xr + atomJ.labFrameDipole_Y*yr + atomJ.labFrameDipole_Z*zr;
    float qkx         = atomJ.labFrameQuadrupole_XX*xr + atomJ.labFrameQuadrupole_XY*yr + atomJ.labFrameQuadrupole_XZ*zr;
    float qky         = atomJ.labFrameQuadrupole_XY*xr + atomJ.labFrameQuadrupole_YY*yr + atomJ.labFrameQuadrupole_YZ*zr;
    float qkz         = atomJ.labFrameQuadrupole_XZ*xr + atomJ.labFrameQuadrupole_YZ*yr + atomJ.labFrameQuadrupole_ZZ*zr;
    float qkr         = qkx*xr + qky*yr + qkz*zr;

    float fim[3],fkm[3];
    float fid[3],fkd[3];
    float fip[3],fkp[3];
    fim[0]            = -xr*(bn[1]*atomJ.q-bn[2]*dkr+bn[3]*qkr)
                         - bn[1]*atomJ.labFrameDipole_X + 2.0f*bn[2]*qkx;

    fim[1]            = -yr*(bn[1]*atomJ.q-bn[2]*dkr+bn[3]*qkr)
                         - bn[1]*atomJ.labFrameDipole_Y + 2.0f*bn[2]*qky;

    fim[2]            = -zr*(bn[1]*atomJ.q-bn[2]*dkr+bn[3]*qkr)
                         - bn[1]*atomJ.labFrameDipole_Z + 2.0f*bn[2]*qkz;

    fkm[0]            = xr*(bn[1]*atomI.q+bn[2]*dir+bn[3]*qir)
                         - bn[1]*atomI.labFrameDipole_X - 2.0f*bn[2]*qix;

    fkm[1]            = yr*(bn[1]*atomI.q+bn[2]*dir+bn[3]*qir)
                         - bn[1]*atomI.labFrameDipole_Y - 2.0f*bn[2]*qiy;

    fkm[2]            = zr*(bn[1]*atomI.q+bn[2]*dir+bn[3]*qir)
                         - bn[1]*atomI.labFrameDipole_Z - 2.0f*bn[2]*qiz;

    fid[0]            = -xr*(drr3*atomJ.q-drr5*dkr+drr7*qkr)
                         - drr3*atomJ.labFrameDipole_X + 2.0f*drr5*qkx;

    fid[1]            = -yr*(drr3*atomJ.q-drr5*dkr+drr7*qkr)
                         - drr3*atomJ.labFrameDipole_Y + 2.0f*drr5*qky;

    fid[2]            = -zr*(drr3*atomJ.q-drr5*dkr+drr7*qkr)
                         - drr3*atomJ.labFrameDipole_Z + 2.0f*drr5*qkz;

    fkd[0]            = xr*(drr3*atomI.q+drr5*dir+drr7*qir)
                         - drr3*atomI.labFrameDipole_X - 2.0f*drr5*qix;

    fkd[1]            = yr*(drr3*atomI.q+drr5*dir+drr7*qir)
                         - drr3*atomI.labFrameDipole_Y - 2.0f*drr5*qiy;

    fkd[2]            = zr*(drr3*atomI.q+drr5*dir+drr7*qir)
                         - drr3*atomI.labFrameDipole_Z - 2.0f*drr5*qiz;

    fip[0]            = -xr*(prr3*atomJ.q-prr5*dkr+prr7*qkr)
                         - prr3*atomJ.labFrameDipole_X + 2.0f*prr5*qkx;

    fip[1]            = -yr*(prr3*atomJ.q-prr5*dkr+prr7*qkr)
                         - prr3*atomJ.labFrameDipole_Y + 2.0f*prr5*qky;

    fip[2]            = -zr*(prr3*atomJ.q-prr5*dkr+prr7*qkr)
                         - prr3*atomJ.labFrameDipole_Z + 2.0f*prr5*qkz;

    fkp[0]            = xr*(prr3*atomI.q+prr5*dir+prr7*qir)
                         - prr3*atomI.labFrameDipole_X - 2.0f*prr5*qix;

    fkp[1]            = yr*(prr3*atomI.q+prr5*dir+prr7*qir)
                         - prr3*atomI.labFrameDipole_Y - 2.0f*prr5*qiy;

    fkp[2]            = zr*(prr3*atomI.q+prr5*dir+prr7*qir)
                         - prr3*atomI.labFrameDipole_Z - 2.0f*prr5*qiz;
  
    // increment the field at each site due to this interaction

    if( r2 <= cAmoebaSim.cutoffDistance2 ){

        fields[0][0]       = fim[0] - fid[0];
        fields[0][1]       = fim[1] - fid[1];
        fields[0][2]       = fim[2] - fid[2];

        fields[1][0]       = fkm[0] - fkd[0];
        fields[1][1]       = fkm[1] - fkd[1];
        fields[1][2]       = fkm[2] - fkd[2];

        fields[2][0]       = fim[0] - fip[0];
        fields[2][1]       = fim[1] - fip[1];
        fields[2][2]       = fim[2] - fip[2];

        fields[3][0]       = fkm[0] - fkp[0];
        fields[3][1]       = fkm[1] - fkp[1];
        fields[3][2]       = fkm[2] - fkp[2];
 
    } else {

        fields[0][0]       = 0.0f;
        fields[1][0]       = 0.0f;
        fields[2][0]       = 0.0f;
        fields[3][0]       = 0.0f;
    
        fields[0][1]       = 0.0f;
        fields[1][1]       = 0.0f;
        fields[2][1]       = 0.0f;
        fields[3][1]       = 0.0f;
    
        fields[0][2]       = 0.0f;
        fields[1][2]       = 0.0f;
        fields[2][2]       = 0.0f;
        fields[3][2]       = 0.0f;
    }
#ifdef AMOEBA_DEBUG
    pullBack[0].x = xr;
    pullBack[0].y = yr;
    pullBack[0].z = zr;
    pullBack[0].w = r2;

    pullBack[1].x = atomJ.x - atomI.x;
    pullBack[1].y = atomJ.y - atomI.y;
    pullBack[1].z = atomJ.z - atomI.z;
    pullBack[1].w = (atomJ.x - atomI.x)*(atomJ.x - atomI.x) + (atomJ.y - atomI.y)*(atomJ.y - atomI.y)+ (atomJ.z - atomI.z)*(atomJ.z - atomI.z);

    pullBack[2].x = scale3;
    pullBack[2].y = scale5;
    pullBack[2].z = scale7;
    pullBack[2].w = -1.0;

#endif
}

// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateAmoebaCudaPmeFixedEField.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateAmoebaCudaPmeFixedEField.h"

/**---------------------------------------------------------------------------------------

   Report whether a number is a nan or infinity

   @param number               number to test
   @return 1 if number is  nan or infinity; else return 0

   --------------------------------------------------------------------------------------- */

#ifdef AMOEBA_DEBUG
static int isNanOrInfinity( double number ){
    return (number != number || number == std::numeric_limits<double>::infinity() || number == -std::numeric_limits<double>::infinity()) ? 1 : 0; 
}
#endif

/**---------------------------------------------------------------------------------------

   Compute fixed electric field using PME

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

static void cudaComputeAmoebaPmeDirectFixedEField( amoebaGpuContext amoebaGpu )
{
  
    gpuContext gpu    = amoebaGpu->gpuContext;

#ifdef AMOEBA_DEBUG
    static const char* methodName = "computeCudaAmoebaPmeFixedEField";
    if( amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "\n%s\n", methodName ); (void) fflush( amoebaGpu->log );
    }
    int paddedNumberOfAtoms                    = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;

    // N2 debug array

    CUDAStream<float4>* debugArray             = new CUDAStream<float4>(paddedNumberOfAtoms*paddedNumberOfAtoms, 1, "DebugArray");
    memset( debugArray->_pSysStream[0], 0, sizeof( float )*4*paddedNumberOfAtoms*paddedNumberOfAtoms);
    debugArray->Upload();

    // print intermediate results for the targetAtom 

    unsigned int targetAtom  = 0;

    int maxPrint             = 3002;
    amoebaGpu->psE_Field->Download();
    (void) fprintf( amoebaGpu->log, "Recip EFields In\n" );
    for( int ii = 0; ii < gpu->natoms; ii++ ){
        (void) fprintf( amoebaGpu->log, "%5d ", ii); 

        int indexOffset     = ii*3;

        // E_Field

        int isNan  = isNanOrInfinity( amoebaGpu->psE_Field->_pSysStream[0][indexOffset] );
            isNan += isNanOrInfinity( amoebaGpu->psE_Field->_pSysStream[0][indexOffset+1] );
            isNan += isNanOrInfinity( amoebaGpu->psE_Field->_pSysStream[0][indexOffset+2] );

        (void) fprintf( amoebaGpu->log,"E[%16.9e %16.9e %16.9e] %s\n",
                        amoebaGpu->psE_Field->_pSysStream[0][indexOffset],
                        amoebaGpu->psE_Field->_pSysStream[0][indexOffset+1],
                        amoebaGpu->psE_Field->_pSysStream[0][indexOffset+2], (isNan ? "XXX" :"") );
   
        if( ii == maxPrint && (gpu->natoms - maxPrint) > ii ){
            ii = gpu->natoms - maxPrint;
        }
    }
    (void) fflush( amoebaGpu->log );
    (void) fprintf( amoebaGpu->log, "Recip EFields End\n" );
#endif

    kClearFields_3( amoebaGpu, 2 );

    if (gpu->bOutputBufferPerWarp){
        kCalculateAmoebaPmeDirectFixedE_FieldN2ByWarp_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->nonbondThreadsPerBlock, sizeof(FixedFieldParticle)*amoebaGpu->nonbondThreadsPerBlock>>>(
                                                                           amoebaGpu->psWorkUnit->_pDevStream[0],
                                                                           amoebaGpu->psWorkArray_3_1->_pDevStream[0],
#ifdef AMOEBA_DEBUG
                                                                           amoebaGpu->psWorkArray_3_2->_pDevStream[0],
                                                                           debugArray->_pDevStream[0], targetAtom );
#else
                                                                           amoebaGpu->psWorkArray_3_2->_pDevStream[0] );
#endif
    } else {

        kCalculateAmoebaPmeDirectFixedE_FieldN2_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->nonbondThreadsPerBlock, sizeof(FixedFieldParticle)*amoebaGpu->nonbondThreadsPerBlock>>>(
                                                                           amoebaGpu->psWorkUnit->_pDevStream[0],
                                                                           amoebaGpu->psWorkArray_3_1->_pDevStream[0],
#ifdef AMOEBA_DEBUG
                                                                           amoebaGpu->psWorkArray_3_2->_pDevStream[0],
                                                                           debugArray->_pDevStream[0], targetAtom );
#else
                                                                           amoebaGpu->psWorkArray_3_2->_pDevStream[0] );
#endif
    }
    LAUNCHERROR("kCalculateAmoebaPmeDirectFixedE_Field_kernel");

#if 0
        for( unsigned int ii = 0; ii < amoebaGpu->outputBuffers; ii++ ){
            //float index = 1.0f;
            float index = (float) ii;
            for( unsigned int jj = 0; jj < 3*amoebaGpu->paddedNumberOfAtoms; jj += 3 ){
                unsigned int kk = 3*ii*amoebaGpu->paddedNumberOfAtoms + jj;
                amoebaGpu->psWorkArray_3_1->_pSysStream[0][kk]   = index;
                amoebaGpu->psWorkArray_3_1->_pSysStream[0][kk+1] = index;
                amoebaGpu->psWorkArray_3_1->_pSysStream[0][kk+2] = index;
            }
        }
        amoebaGpu->psWorkArray_3_1->Upload();
#endif

    kReducePmeDirectE_Fields( amoebaGpu );

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){
        gpu->psInteractionCount->Download();
        (void) fprintf( amoebaGpu->log, "AmoebaN2Forces_kernel numBlocks=%u numThreads=%u bufferPerWarp=%u atm=%u shrd=%u Ebuf=%u ixnCt=%u workUnits=%u warp=%d\n",
                        amoebaGpu->nonbondBlocks, amoebaGpu->nonbondThreadsPerBlock, amoebaGpu->bOutputBufferPerWarp,
                        sizeof(FixedFieldParticle), sizeof(FixedFieldParticle)*amoebaGpu->nonbondThreadsPerBlock, amoebaGpu->energyOutputBuffers, 
                        (*gpu->psInteractionCount)[0], gpu->sim.workUnits, gpu->bOutputBufferPerWarp );
        (void) fflush( amoebaGpu->log );
/*
        (void) fprintf( amoebaGpu->log, "Out WorkArray_3_[1,2]  paddedNumberOfAtoms=%d\n",  amoebaGpu->paddedNumberOfAtoms, amoebaGpu->outputBuffers );
        amoebaGpu->psWorkArray_3_1->Download();
        amoebaGpu->psWorkArray_3_2->Download();
        for( int ii = 0; ii < amoebaGpu->paddedNumberOfAtoms; ii++ ){
           (void) fprintf( amoebaGpu->log, "%5d ", ii); 

            int indexOffset     = ii*3;

           // buffer 1

           (void) fprintf( amoebaGpu->log,"WArry1[%16.9e %16.9e %16.9e] ",
                           amoebaGpu->psWorkArray_3_1->_pSysStream[0][indexOffset],
                           amoebaGpu->psWorkArray_3_1->_pSysStream[0][indexOffset+1],
                           amoebaGpu->psWorkArray_3_1->_pSysStream[0][indexOffset+2] );
   
           // buffer 2

           (void) fprintf( amoebaGpu->log,"WArry2[%16.9e %16.9e %16.9e] ",
                           amoebaGpu->psWorkArray_3_2->_pSysStream[0][indexOffset],
                           amoebaGpu->psWorkArray_3_2->_pSysStream[0][indexOffset+1],
                           amoebaGpu->psWorkArray_3_2->_pSysStream[0][indexOffset+2] );

           (void) fprintf( amoebaGpu->log,"\n" );
           if( ii == maxPrint && (gpu->natoms - maxPrint) > ii ){
                ii = gpu->natoms - maxPrint;
           }
        }
        (void) fflush( amoebaGpu->log );
*/
        amoebaGpu->psE_Field->Download();
        amoebaGpu->psE_FieldPolar->Download();
        for( int ii = 0; ii < gpu->natoms; ii++ ){
           (void) fprintf( amoebaGpu->log, "%5d ", ii); 

            int indexOffset     = ii*3;

           // E_Field

           (void) fprintf( amoebaGpu->log,"E[%16.9e %16.9e %16.9e] ",
                           amoebaGpu->psE_Field->_pSysStream[0][indexOffset],
                           amoebaGpu->psE_Field->_pSysStream[0][indexOffset+1],
                           amoebaGpu->psE_Field->_pSysStream[0][indexOffset+2] );
   
           // E_Field polar

           (void) fprintf( amoebaGpu->log,"Epol[%16.9e %16.9e %16.9e] ",
                           amoebaGpu->psE_FieldPolar->_pSysStream[0][indexOffset],
                           amoebaGpu->psE_FieldPolar->_pSysStream[0][indexOffset+1],
                           amoebaGpu->psE_FieldPolar->_pSysStream[0][indexOffset+2] );

           (void) fprintf( amoebaGpu->log,"\n" );
           if( ii == maxPrint && (gpu->natoms - maxPrint) > ii ){
                ii = gpu->natoms - maxPrint;
           }
        }
        (void) fflush( amoebaGpu->log );
        (void) fprintf( amoebaGpu->log, "EFields End\n" );

        (void) fprintf( amoebaGpu->log, "DebugQ\n" );
        debugArray->Download();

        int paddedNumberOfAtoms                    = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
        for( int jj = 0; jj < gpu->natoms; jj++ ){
            int debugIndex = jj;
            (void) fprintf( amoebaGpu->log,"%5d PmeFixedEField\n", jj );
            for( int kk = 0; kk < 10; kk++ ){
                (void) fprintf( amoebaGpu->log,"[%16.9e %16.9e %16.9e %16.9e]\n",
                                debugArray->_pSysStream[0][debugIndex].x, debugArray->_pSysStream[0][debugIndex].y,
                                debugArray->_pSysStream[0][debugIndex].z, debugArray->_pSysStream[0][debugIndex].w );
                debugIndex += paddedNumberOfAtoms;
            }
            (void) fprintf( amoebaGpu->log,"\n" );

        }

        // write results to file

        if( 1 ){
            std::vector<int> fileId;
            //fileId.push_back( 0 );
            VectorOfDoubleVectors outputVector;
            cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,              outputVector );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_Field,      outputVector );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_FieldPolar, outputVector);
            cudaWriteVectorOfDoubleVectorsToFile( "CudaEField", fileId, outputVector );

         }
         delete debugArray;
    }
#endif

        if( 0 ){
            std::vector<int> fileId;
            fileId.push_back( 0 );
            VectorOfDoubleVectors outputVector;
            cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,              outputVector );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_Field,      outputVector );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_FieldPolar, outputVector);
            cudaWriteVectorOfDoubleVectorsToFile( "CudaEField", fileId, outputVector );

         }
}

void cudaComputeAmoebaPmeFixedEField( amoebaGpuContext amoebaGpu )
{
    kCalculateAmoebaPMEFixedMultipoleField( amoebaGpu );
    cudaComputeAmoebaPmeDirectFixedEField( amoebaGpu );
}
