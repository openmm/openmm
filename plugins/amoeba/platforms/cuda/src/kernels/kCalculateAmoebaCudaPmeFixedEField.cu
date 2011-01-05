
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
#undef INCLUDE_FIXED_FIELD_BUFFERS
#define INCLUDE_FIXED_FIELD_BUFFERS
#include "kCalculateAmoebaCudaFixedFieldParticle.h"
#undef INCLUDE_FIXED_FIELD_BUFFERS
__device__ void sumTempBuffer( FixedFieldParticle& atomI, FixedFieldParticle& atomJ ){
    atomI.tempBuffer[0]  += atomJ.tempBuffer[0];
    atomI.tempBuffer[1]  += atomJ.tempBuffer[1];
    atomI.tempBuffer[2]  += atomJ.tempBuffer[2];

    atomI.tempBufferP[0] += atomJ.tempBufferP[0];
    atomI.tempBufferP[1] += atomJ.tempBufferP[1];
    atomI.tempBufferP[2] += atomJ.tempBufferP[2];
}

__device__ void calculateFixedFieldRealSpacePairIxn_kernel( FixedFieldParticle& atomI, FixedFieldParticle& atomJ,
                                                            float dscale, float pscale, float4 fields[3]
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

    float r2          = xr*xr + yr*yr + zr*zr;
    if( r2 <= cSim.nonbondedCutoffSqr ){
        float r           = sqrtf(r2);

        // calculate the error function damping terms

        float ralpha      = cSim.alphaEwald*r;

        float bn0         = erfc(ralpha)/r;
        float alsq2       = 2.0f*cSim.alphaEwald*cSim.alphaEwald;
        float alsq2n      = 1.0f/(cAmoebaSim.sqrtPi*cSim.alphaEwald);
        float exp2a       = exp(-(ralpha*ralpha));
        alsq2n           *= alsq2;
        float bn1         = (bn0+alsq2n*exp2a)/r2;

        alsq2n           *= alsq2;
        float bn2         = (3.0f*bn1+alsq2n*exp2a)/r2;

        alsq2n           *= alsq2;
        float bn3         = (5.0f*bn2+alsq2n*exp2a)/r2;

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

        float fim0            = -xr*(bn1*atomJ.q-bn2*dkr+bn3*qkr)
                             - bn1*atomJ.labFrameDipole_X + 2.0f*bn2*qkx;

        float fim1            = -yr*(bn1*atomJ.q-bn2*dkr+bn3*qkr)
                             - bn1*atomJ.labFrameDipole_Y + 2.0f*bn2*qky;

        float fim2            = -zr*(bn1*atomJ.q-bn2*dkr+bn3*qkr)
                             - bn1*atomJ.labFrameDipole_Z + 2.0f*bn2*qkz;

        float fkm0            = xr*(bn1*atomI.q+bn2*dir+bn3*qir)
                             - bn1*atomI.labFrameDipole_X - 2.0f*bn2*qix;

        float fkm1            = yr*(bn1*atomI.q+bn2*dir+bn3*qir)
                             - bn1*atomI.labFrameDipole_Y - 2.0f*bn2*qiy;

        float fkm2            = zr*(bn1*atomI.q+bn2*dir+bn3*qir)
                             - bn1*atomI.labFrameDipole_Z - 2.0f*bn2*qiz;

        float fid0            = -xr*(drr3*atomJ.q-drr5*dkr+drr7*qkr)
                             - drr3*atomJ.labFrameDipole_X + 2.0f*drr5*qkx;

        float fid1            = -yr*(drr3*atomJ.q-drr5*dkr+drr7*qkr)
                             - drr3*atomJ.labFrameDipole_Y + 2.0f*drr5*qky;

        float fid2            = -zr*(drr3*atomJ.q-drr5*dkr+drr7*qkr)
                             - drr3*atomJ.labFrameDipole_Z + 2.0f*drr5*qkz;

        float fkd0            = xr*(drr3*atomI.q+drr5*dir+drr7*qir)
                             - drr3*atomI.labFrameDipole_X - 2.0f*drr5*qix;

        float fkd1            = yr*(drr3*atomI.q+drr5*dir+drr7*qir)
                             - drr3*atomI.labFrameDipole_Y - 2.0f*drr5*qiy;

        float fkd2            = zr*(drr3*atomI.q+drr5*dir+drr7*qir)
                             - drr3*atomI.labFrameDipole_Z - 2.0f*drr5*qiz;

        float fip0            = -xr*(prr3*atomJ.q-prr5*dkr+prr7*qkr)
                             - prr3*atomJ.labFrameDipole_X + 2.0f*prr5*qkx;

        float fip1            = -yr*(prr3*atomJ.q-prr5*dkr+prr7*qkr)
                             - prr3*atomJ.labFrameDipole_Y + 2.0f*prr5*qky;

        float fip2            = -zr*(prr3*atomJ.q-prr5*dkr+prr7*qkr)
                             - prr3*atomJ.labFrameDipole_Z + 2.0f*prr5*qkz;

        float fkp0            = xr*(prr3*atomI.q+prr5*dir+prr7*qir)
                             - prr3*atomI.labFrameDipole_X - 2.0f*prr5*qix;

        float fkp1            = yr*(prr3*atomI.q+prr5*dir+prr7*qir)
                             - prr3*atomI.labFrameDipole_Y - 2.0f*prr5*qiy;

        float fkp2            = zr*(prr3*atomI.q+prr5*dir+prr7*qir)
                             - prr3*atomI.labFrameDipole_Z - 2.0f*prr5*qiz;

        // increment the field at each site due to this interaction

        fields[0].x       = fim0 - fid0;
        fields[1].x       = fim1 - fid1;
        fields[2].x       = fim2 - fid2;

        fields[0].y       = fkm0 - fkd0;
        fields[1].y       = fkm1 - fkd1;
        fields[2].y       = fkm2 - fkd2;

        fields[0].z       = fim0 - fip0;
        fields[1].z       = fim1 - fip1;
        fields[2].z       = fim2 - fip2;

        fields[0].w       = fkm0 - fkp0;
        fields[1].w       = fkm1 - fkp1;
        fields[2].w       = fkm2 - fkp2;
 
    } else {

        fields[0].x       = 0.0f;
        fields[0].y       = 0.0f;
        fields[0].z       = 0.0f;
        fields[0].w       = 0.0f;
    
        fields[1].x       = 0.0f;
        fields[1].y       = 0.0f;
        fields[1].z       = 0.0f;
        fields[1].w       = 0.0f;
    
        fields[2].x       = 0.0f;
        fields[2].y       = 0.0f;
        fields[2].z       = 0.0f;
        fields[2].w       = 0.0f;
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

#define METHOD_NAME(a, b) a##Cutoff##b
#include "kCalculateAmoebaCudaPmeFixedEField.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##CutoffByWarp##b
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
  
    static unsigned int threadsPerBlock  = 0;
    gpuContext gpu                       = amoebaGpu->gpuContext;

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

    unsigned int targetAtom  = 1280;
#endif

    kClearFields_3( amoebaGpu, 2 );

    // on first pass, set threads/block

    if( threadsPerBlock == 0 ){ 
        unsigned int maxThreads;
        if (gpu->sm_version >= SM_20)
            maxThreads = 384; 
        else if (gpu->sm_version >= SM_12)
            maxThreads = 192;
        else
            maxThreads = 64;
        threadsPerBlock = std::min(getThreadsPerBlock(amoebaGpu, sizeof(FixedFieldParticle)), maxThreads);
    }    

    if (gpu->bOutputBufferPerWarp){
        kCalculateAmoebaPmeDirectFixedE_FieldCutoffByWarp_kernel<<<amoebaGpu->nonbondBlocks, threadsPerBlock, sizeof(FixedFieldParticle)*threadsPerBlock>>>(
                                                                           gpu->sim.pInteractingWorkUnit,
                                                                           amoebaGpu->psWorkArray_3_1->_pDevStream[0],
#ifdef AMOEBA_DEBUG
                                                                           amoebaGpu->psWorkArray_3_2->_pDevStream[0],
                                                                           debugArray->_pDevStream[0], targetAtom );
#else
                                                                           amoebaGpu->psWorkArray_3_2->_pDevStream[0] );
#endif
    } else {
        kCalculateAmoebaPmeDirectFixedE_FieldCutoff_kernel<<<amoebaGpu->nonbondBlocks, threadsPerBlock, sizeof(FixedFieldParticle)*threadsPerBlock>>>(
                                                                           gpu->sim.pInteractingWorkUnit,
                                                                           amoebaGpu->psWorkArray_3_1->_pDevStream[0],
#ifdef AMOEBA_DEBUG
                                                                           amoebaGpu->psWorkArray_3_2->_pDevStream[0],
                                                                           debugArray->_pDevStream[0], targetAtom );
#else
                                                                           amoebaGpu->psWorkArray_3_2->_pDevStream[0] );
#endif
    }
    LAUNCHERROR("kCalculateAmoebaPmeDirectFixedE_Field_kernel");

    kReducePmeDirectE_Fields( amoebaGpu );

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){
        gpu->psInteractionCount->Download();
        (void) fprintf( amoebaGpu->log, "cudaComputeAmoebaPmeDirectFixedEField:  threadsPerBlock=%u getThreadsPerBlock=%d sizeof=%u shrd=%u\n", 
                        threadsPerBlock, getThreadsPerBlock(amoebaGpu, sizeof(FixedFieldParticle)+sizeof(float3)),
                        (sizeof(FixedFieldParticle)+sizeof(float3)), (sizeof(FixedFieldParticle)+sizeof(float3))*threadsPerBlock );
        (void) fprintf( amoebaGpu->log, "AmoebaCutoffForces_kernel numBlocks=%u numThreads=%u bufferPerWarp=%u atm=%u shrd=%u Ebuf=%u ixnCt=%u workUnits=%u warp=%d\n",
                        amoebaGpu->nonbondBlocks, threadsPerBlock, amoebaGpu->bOutputBufferPerWarp,
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
        (void) fprintf( amoebaGpu->log,"E-field (includes self term)" );
        int maxPrint             = 3002;
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
        amoebaGpu->gpuContext->psPosq4->Download();
        for( int jj = 0; jj < gpu->natoms; jj++ ){
            int debugIndex = jj;
if( fabs(debugArray->_pSysStream[0][jj+3*paddedNumberOfAtoms].x) > 0.0 ){
            (void) fprintf( amoebaGpu->log,"%5d PmeFixedEField\n", jj );
            for( int kk = 0; kk < 7; kk++ ){
                (void) fprintf( amoebaGpu->log,"[%16.9e %16.9e %16.9e %16.9e]\n",
                                debugArray->_pSysStream[0][debugIndex].x, debugArray->_pSysStream[0][debugIndex].y,
                                debugArray->_pSysStream[0][debugIndex].z, debugArray->_pSysStream[0][debugIndex].w );
                debugIndex += paddedNumberOfAtoms;
            }
            (void) fprintf( amoebaGpu->log,"\n" );
}

        }

        // write results to file

        if( 1 ){
            std::vector<int> fileId;
            //fileId.push_back( 0 );
            VectorOfDoubleVectors outputVector;
            cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,              outputVector, gpu->psAtomIndex->_pSysData );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_Field,      outputVector, gpu->psAtomIndex->_pSysData );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_FieldPolar, outputVector, gpu->psAtomIndex->_pSysData );
            cudaWriteVectorOfDoubleVectorsToFile( "CudaEField", fileId, outputVector );
         }
         delete debugArray;
    }
#endif

        if( 0 ){
            std::vector<int> fileId;
            fileId.push_back( 0 );
            VectorOfDoubleVectors outputVector;
            cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,              outputVector, gpu->psAtomIndex->_pSysData );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_Field,      outputVector, gpu->psAtomIndex->_pSysData );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_FieldPolar, outputVector, gpu->psAtomIndex->_pSysData);
            cudaWriteVectorOfDoubleVectorsToFile( "CudaEField", fileId, outputVector );

         }
}

void cudaComputeAmoebaPmeFixedEField( amoebaGpuContext amoebaGpu )
{
    kCalculateAmoebaPMEFixedMultipoles( amoebaGpu );
    cudaComputeAmoebaPmeDirectFixedEField( amoebaGpu );

    if( 0 ){
        gpuContext gpu                       = amoebaGpu->gpuContext;
        std::vector<int> fileId;
        fileId.push_back( 0 );
        VectorOfDoubleVectors outputVector;
        cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,              outputVector, gpu->psAtomIndex->_pSysData );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_Field,      outputVector, gpu->psAtomIndex->_pSysData );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psE_FieldPolar, outputVector, gpu->psAtomIndex->_pSysData );
        cudaWriteVectorOfDoubleVectorsToFile( "CudaEField", fileId, outputVector );
    }
}
