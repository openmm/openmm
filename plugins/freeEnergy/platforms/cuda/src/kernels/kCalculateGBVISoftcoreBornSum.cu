/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Scott Le Grand, Peter Eastman                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "GpuGBVISoftcore.h"
#include "GpuFreeEnergyCudaKernels.h"
#include <cuda.h>

struct cudaFreeEnergySimulationGBVI {
    float quinticLowerLimitFactor;
    float quinticUpperLimit;
    float* pSwitchDerivative;
};
struct cudaFreeEnergySimulationGBVI gbviSim;

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaFreeEnergySimulationGBVI gbviSimDev;

void SetCalculateGBVISoftcoreBornSumGpuSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetCalculateGBVISoftcoreBornSumGpuSim copy to cSim failed");
    //(void) fprintf( stderr, "SetCalculateGBVISoftcoreBornSumGpuSim\n" );
}

void GetCalculateGBVISoftcoreBornSumSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

void SetCalculateGBVISoftcoreSupplementarySim( GpuGBVISoftcore* gpuGBVISoftcore )
{
    cudaError_t status;
    gbviSim.pSwitchDerivative        = gpuGBVISoftcore->getGpuSwitchDerivative();
    gbviSim.quinticLowerLimitFactor  = gpuGBVISoftcore->getQuinticLowerLimitFactor();
    gbviSim.quinticUpperLimit        = gpuGBVISoftcore->getQuinticUpperLimit();
    status                           = cudaMemcpyToSymbol(gbviSimDev, &gbviSim, sizeof(cudaFreeEnergySimulationGBVI));
    RTERROR(status, "cudaMemcpyToSymbol: SetCalculateGBVISoftcoreSupplementarySim");

    //(void) fprintf( stderr, "SetCalculateGBVISoftcoreSupplementarySim %14.6e %14.6e swDerv=%p\n",
    //                gbviSim.quinticLowerLimitFactor, gbviSim.quinticUpperLimit, gbviSim.pSwitchDerivative );
}

// create, initialize and enter BornRadiusScaleFactors values (used to scale contribution of atoms to Born sum of other atoms)
// return handle to GpuGBVISoftcore object

extern "C"
GpuGBVISoftcore* gpuSetGBVISoftcoreParameters(gpuContext gpu, float innerDielectric, float solventDielectric, const std::vector<int>& atom,
                                              const std::vector<float>& radius, const std::vector<float>& gamma,
                                              const std::vector<float>& scaledRadii, const std::vector<float>& bornRadiusScaleFactors,
                                              const std::vector<float>& quinticSplineParameters)
{
    static const float electricConstant         = -166.02691f;
    unsigned int atoms                          = atom.size();
    double tau                                  = ((1.0f/innerDielectric)-(1.0f/solventDielectric)); 

    // create gpuGBVISoftcore, load parameters, and track minimum softcore value
    // gpuGBVISoftcore is not really being used (it was in the initial implementation) -- 
    // will be removed in future once confirmed not needed

    GpuGBVISoftcore* gpuGBVISoftcore            = new GpuGBVISoftcore();
    unsigned int numberOfParticles              = radius.size();

    // check if quintic scaling to be applied

    if( quinticSplineParameters.size() == 2 ){
       gpuGBVISoftcore->setBornRadiiScalingMethod( 1 );
       gpuGBVISoftcore->setQuinticLowerLimitFactor( quinticSplineParameters[0] );
       gpuGBVISoftcore->setQuinticUpperLimit(       quinticSplineParameters[1] );
       gpuGBVISoftcore->initializeGpuSwitchDerivative(  gpu->sim.paddedNumberOfAtoms );
    }

    for (unsigned int i = 0; i < bornRadiusScaleFactors.size(); i++) 
    {
            (*gpu->psGBVIData)[i].x = radius[i];
            (*gpu->psGBVIData)[i].y = scaledRadii[i];
            (*gpu->psGBVIData)[i].z = tau*gamma[i];
            (*gpu->psGBVIData)[i].w = bornRadiusScaleFactors[i];

(*gpu->psObcData)[i].x  = radius[i];
(*gpu->psObcData)[i].y  = 0.9f*radius[i];

    }

    // Dummy out extra atom data

    for (unsigned int i = atoms; i < gpu->sim.paddedNumberOfAtoms; i++)
    {
        (*gpu->psBornRadii)[i]      = 0.2f;
        (*gpu->psGBVIData)[i].x     = 0.01f;
        (*gpu->psGBVIData)[i].y     = 0.01f;
        (*gpu->psGBVIData)[i].z     = 0.01f;
        (*gpu->psGBVIData)[i].w     = 1.00f;
    }

#undef DUMP_PARAMETERS
#define DUMP_PARAMETERS 0
#if (DUMP_PARAMETERS == 1)
    (void) fprintf( stderr,"GBVI softcore param %u %u sclMeth=%d LwFct=%8.3f UpLmt=[%12.5e (nm) %12.5e]\nR scaledR gamma*tau= bornRadiusScaleFactor \n",
                    bornRadiusScaleFactors.size(), gpu->sim.paddedNumberOfAtoms,
                    gpuGBVISoftcore->getBornRadiiScalingMethod(), gpuGBVISoftcore->getQuinticLowerLimitFactor(),
                    powf( gpuGBVISoftcore->getQuinticUpperLimit(), -0.3333333f ), gpuGBVISoftcore->getQuinticUpperLimit() );
    int maxPrint = 31;
    for (unsigned int ii = 0; ii < gpu->sim.paddedNumberOfAtoms; ii++) 
    {

        (void) fprintf( stderr,"%6u %14.7e %14.7e %14.7e %14.7e\n",
                        ii, (*gpu->psGBVIData)[ii].x, (*gpu->psGBVIData)[ii].y, (*gpu->psGBVIData)[ii].z, (*gpu->psGBVIData)[ii].w ); 
        if( ii == maxPrint ){
            ii = gpu->sim.paddedNumberOfAtoms - maxPrint;
            if( ii < maxPrint )ii = maxPrint;
        }
    }
#endif

    gpu->psBornRadii->Upload();
    gpu->psGBVIData->Upload();
gpu->psObcData->Upload();
    gpu->sim.preFactor              = 2.0f*electricConstant*((1.0f/innerDielectric)-(1.0f/solventDielectric))*gpu->sim.forceConversionFactor;
    gpuGBVISoftcore->upload( gpu );

#if (DUMP_PARAMETERS == 1)
(void) fprintf( stderr, "gpuSetGBVISoftcoreParameters: preFactor=%14.6e elecCnstnt=%.4f frcCnvrsnFctr=%.4f tau=%.4f.\n",
                gpu->sim.preFactor, 2.0f*electricConstant, gpu->sim.forceConversionFactor, ((1.0f/innerDielectric)-(1.0f/solventDielectric)) );
#endif

    return gpuGBVISoftcore;

}

struct Atom {
    float x;
    float y;
    float z;
    float r;
    float sr;
    float sum;
    float gamma;
    float bornRadiusScaleFactor;
};

__global__ void kClearGBVISoftcoreBornSum_kernel()
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.stride * cSim.nonbondOutputBuffers)
    {   
        ((float*)cSim.pBornSum)[pos] = 0.0f;
        pos += gridDim.x * blockDim.x;
    }   
}

void kClearGBVISoftcoreBornSum(gpuContext gpu) {
    kClearGBVISoftcoreBornSum_kernel<<<gpu->sim.blocks, 384>>>();
}

__global__ void kReduceGBVISoftcoreBornForces_kernel()
{
    unsigned int pos = (blockIdx.x * blockDim.x + threadIdx.x);
    float energy = 0.0f;
    while (pos < cSim.atoms)
    {
        float bornRadius  = cSim.pBornRadii[pos];
        float4 gbviData   = cSim.pGBVIData[pos];
        float totalForce  = 0.0f;
        float* pFt        = cSim.pBornForce + pos;

        int i = cSim.nonbondOutputBuffers;
        while (i >= 4)
        {
            float f1    = *pFt;
            pFt        += cSim.stride;
            float f2    = *pFt;
            pFt        += cSim.stride;
            float f3    = *pFt;
            pFt        += cSim.stride;
            float f4    = *pFt;
            pFt        += cSim.stride;
            totalForce += f1 + f2 + f3 + f4;
            i -= 4;
        }
        if (i >= 2)
        {
            float f1    = *pFt;
            pFt        += cSim.stride;
            float f2    = *pFt;
            pFt        += cSim.stride;
            totalForce += f1 + f2;
            i -= 2;
        }
        if (i > 0)
        {
            totalForce += *pFt;
        }

        float ratio         = (gbviData.x/bornRadius);
        float ratio3        = ratio*ratio*ratio;
        energy             -= gbviData.z*ratio3;
        totalForce         += (3.0f*gbviData.z*ratio3)/bornRadius; // 'cavity' term
        float br2           = bornRadius*bornRadius;
        totalForce         *= (1.0f/3.0f)*br2*br2;

        pFt = cSim.pBornForce + pos;
        *pFt = totalForce;
        pos += gridDim.x * blockDim.x;
    }
    cSim.pEnergy[blockIdx.x * blockDim.x + threadIdx.x] += energy;
}

void kReduceGBVISoftcoreBornForces(gpuContext gpu)
{
    kReduceGBVISoftcoreBornForces_kernel<<<gpu->sim.blocks, gpu->sim.bf_reduce_threads_per_block>>>();
    LAUNCHERROR("kReduceGBVISoftcoreBornForces");

}

__global__ void kReduceGBVISoftcoreBornSum_kernel()
{
    unsigned int pos = (blockIdx.x * blockDim.x + threadIdx.x);
    
    while (pos < cSim.atoms)
    {
        float sum = 0.0f;
        float* pSt = cSim.pBornSum + pos;
        float4 atom = cSim.pGBVIData[pos];
        
        // Get summed Born data
        for (int i = 0; i < cSim.nonbondOutputBuffers; i++)
        {
            sum += *pSt;
       //     printf("%4d %4d A: %9.4f\n", pos, i, *pSt);
            pSt += cSim.stride;
        }
        
        // Now calculate Born radius

        float Rinv           = 1.0f/atom.x;
        sum                  = Rinv*Rinv*Rinv - sum; 
        cSim.pBornRadii[pos] = pow( sum, (-1.0f/3.0f) ); 
        pos += gridDim.x * blockDim.x;
    }   
}

void kReduceGBVISoftcoreBornSum(gpuContext gpu)
{
    //printf("kReduceGBVISoftcoreBornSum\n");
#define GBVISoftcore_DEBUG 0
#if ( GBVISoftcore_DEBUG == 1 )
               gpu->psGBVISoftcoreData->Download();
               gpu->psBornSum->Download();
               gpu->psPosq4->Download();
                (void) fprintf( stderr, "\nkReduceGBVISoftcoreBornSum: Post BornSum %s Born radii & params\n", 
                               (gpu->bIncludeGBVISoftcore ? "GBVI" : "Obc") );
                for( int ii = 0; ii < gpu->natoms; ii++ ){
                   (void) fprintf( stderr, "%d bSum=%14.6e param[%14.6e %14.6e %14.6e] x[%14.6f %14.6f %14.6f %14.6f]\n",
                                   ii, 
                                   gpu->psBornSum->_pSysStream[0][ii],
                                   gpu->psGBVISoftcoreData->_pSysStream[0][ii].x,
                                   gpu->psGBVISoftcoreData->_pSysStream[0][ii].y,
                                   gpu->psGBVISoftcoreData->_pSysStream[0][ii].z,
                                   gpu->psPosq4->_pSysStream[0][ii].x, gpu->psPosq4->_pSysStream[0][ii].y,
                                   gpu->psPosq4->_pSysStream[0][ii].z, gpu->psPosq4->_pSysStream[0][ii].w
                                 );  
                }   
#endif
#undef GBVISoftcore_DEBUG


    kReduceGBVISoftcoreBornSum_kernel<<<gpu->sim.blocks, 384>>>();
    gpu->bRecalculateBornRadii = false;
    LAUNCHERROR("kReduceGBVISoftcoreBornSum");
}

// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateGBVISoftcoreBornSum.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateGBVISoftcoreBornSum.h"

// Include versions of the kernels with cutoffs.

#if 0
#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_CUTOFF
#define METHOD_NAME(a, b) a##Cutoff##b
#include "kCalculateGBVISoftcoreBornSum.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##CutoffByWarp##b
#include "kCalculateGBVISoftcoreBornSum.h"

// Include versions of the kernels with periodic boundary conditions.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_PERIODIC
#define METHOD_NAME(a, b) a##Periodic##b
#include "kCalculateGBVISoftcoreBornSum.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##PeriodicByWarp##b
#include "kCalculateGBVISoftcoreBornSum.h"
#endif

#if 0
__global__ void kClearGBVISoftcoreBornSum_kernel()
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.stride * cSim.nonbondOutputBuffers)
    {
        ((float*)cSim.pBornSum)[pos] = 0.0f;
        pos += gridDim.x * blockDim.x;
    }
}
#endif

 __device__ void quinticSpline( float  x, float rl, float ru, float* outValue, float* outDerivative )
{
   float numerator    = x  - rl;
   float denominator  = ru - rl;
   float ratio        = numerator/denominator;
   float ratio2       = ratio*ratio;
   float ratio3       = ratio2*ratio;

   *outValue          = 1.0f + ratio3*(-10.f + 15.0f*ratio - 6.0f*ratio2);
   *outDerivative     = ratio2*(-30.0f + 60.0f*ratio - 30.0f*ratio2)/denominator;
}

__global__ void kReduceGBVIBornSumQuinticScaling_kernel()
{
    unsigned int pos = (blockIdx.x * blockDim.x + threadIdx.x);
    
    while (pos < cSim.atoms)
    {
        float sum = 0.0f;
        float* pSt = cSim.pBornSum + pos;
        float4 atom = cSim.pGBVIData[pos];
        
        // Get summed Born data
        for (int i = 0; i < cSim.nonbondOutputBuffers; i++)
        {
            sum += *pSt;
       //     printf("%4d %4d A: %9.4f\n", pos, i, *pSt);
            pSt += cSim.stride;
        }
        
        // Now calculate Born radius

        float Rinv           = 1.0f/atom.x;
        float r3             = Rinv*Rinv*Rinv;
        float splineL        = gbviSimDev.quinticLowerLimitFactor*r3;
//float bSum           = sum;
        float switchDeriviative;
        if( sum > splineL ){
            if( sum < r3 ){
                float splineValue, splineDerivative;
                quinticSpline( sum, splineL, r3, &splineValue, &splineDerivative ); 
                switchDeriviative  = splineValue - (r3 - sum)*splineDerivative;
                sum                = (r3 - sum)*splineValue + gbviSimDev.quinticUpperLimit;
            } else {
                sum                = gbviSimDev.quinticUpperLimit;
                switchDeriviative  = 0.0f;
            }
        } else {
            sum               = r3 - sum;
            switchDeriviative = 1.0f;
        }

        cSim.pBornRadii[pos]              = pow( sum, (-1.0f/3.0f) ); 
//cSim.pBornSum[pos]              = bSum;
        gbviSimDev.pSwitchDerivative[pos] = switchDeriviative;
        pos += gridDim.x * blockDim.x;
    }   
}

void kReduceGBVIBornSumQuinticScaling(gpuContext gpu, GpuGBVISoftcore* gpuGBVISoftcore)
{
    //printf("kReduceGBVIBornSumQuinticScaling_kernel\n");
    kReduceGBVIBornSumQuinticScaling_kernel<<<gpu->sim.blocks, 384>>>();
    gpu->bRecalculateBornRadii = false;
    LAUNCHERROR("kReduceGBVIBornSumQuinticScaling_kernel");

#define GBVI_DEBUG 0
#if ( GBVI_DEBUG == 1 )
               gpu->psGBVIData->Download();
               gpu->psBornSum->Download();
               gpu->psBornRadii->Download();
               gpu->psPosq4->Download();
               CUDAStream<float>* psSwitchDerivative = gpuGBVISoftcore->getSwitchDerivative();
                
               psSwitchDerivative->Download();
                (void) fprintf( stderr, "\nkReduceGBVIBornSumQuinticScaling: Post BornSum %s Born radii & params\n", 
                               (gpu->bIncludeGBVI ? "GBVI" : "Obc") );
                for( int ii = 0; ii < gpu->natoms; ii++ ){
                   (void) fprintf( stderr, "%6d bSum=%14.6e bR=%14.6e swDerv=%14.6e param[%14.6e %14.6e %14.6e] x[%14.6f %14.6f %14.6f %14.6f] %s\n",
                                   ii, 
                                   gpu->psBornSum->_pSysStream[0][ii],
                                   gpu->psBornRadii->_pSysStream[0][ii],
                                   psSwitchDerivative->_pSysStream[0][ii],
                                   gpu->psGBVIData->_pSysStream[0][ii].x,
                                   gpu->psGBVIData->_pSysStream[0][ii].y,
                                   gpu->psGBVIData->_pSysStream[0][ii].z,
                                   gpu->psPosq4->_pSysStream[0][ii].x, gpu->psPosq4->_pSysStream[0][ii].y,
                                   gpu->psPosq4->_pSysStream[0][ii].z, gpu->psPosq4->_pSysStream[0][ii].w,
                                   (fabs( psSwitchDerivative->_pSysStream[0][ii] - 1.0 ) > 1.0e-05 ? "SWWWWW" : "")
                                 );  
                }   
#endif
#undef GBVI_DEBUG

}

__global__ void kReduceGBVIBornForcesQuinticScaling_kernel()
{
    unsigned int pos = (blockIdx.x * blockDim.x + threadIdx.x);
    float energy = 0.0f;
    while (pos < cSim.atoms)
    {
        float bornRadius    = cSim.pBornRadii[pos];
        float4 gbviData     = cSim.pGBVIData[pos];
        float  switchDeriv  = gbviSimDev.pSwitchDerivative[pos];
        float totalForce    = 0.0f;
        float* pFt          = cSim.pBornForce + pos;

        int i = cSim.nonbondOutputBuffers;
        while (i >= 4)
        {
            float f1    = *pFt;
            pFt        += cSim.stride;
            float f2    = *pFt;
            pFt        += cSim.stride;
            float f3    = *pFt;
            pFt        += cSim.stride;
            float f4    = *pFt;
            pFt        += cSim.stride;
            totalForce += f1 + f2 + f3 + f4;
            i -= 4;
        }
        if (i >= 2)
        {
            float f1    = *pFt;
            pFt        += cSim.stride;
            float f2    = *pFt;
            pFt        += cSim.stride;
            totalForce += f1 + f2;
            i -= 2;
        }
        if (i > 0)
        {
            totalForce += *pFt;
        }

        float ratio         = (gbviData.x/bornRadius);
        float ratio3        = ratio*ratio*ratio;
        energy             -= gbviData.z*ratio3;
        totalForce         += (3.0f*gbviData.z*ratio3)/bornRadius; // 'cavity' term
        float br2           = bornRadius*bornRadius;
        totalForce         *= (1.0f/3.0f)*br2*br2*switchDeriv;

        pFt = cSim.pBornForce + pos;
        *pFt = totalForce;
        pos += gridDim.x * blockDim.x;
    }
    cSim.pEnergy[blockIdx.x * blockDim.x + threadIdx.x] += energy;
}

void kReduceGBVIBornForcesQuinticScaling(gpuContext gpu)
{
    //printf("kReduceObcGbsaBornForces\n");
    kReduceGBVIBornForcesQuinticScaling_kernel<<<gpu->sim.blocks, gpu->sim.bf_reduce_threads_per_block>>>();
    LAUNCHERROR("kReduceGBVIBornForcesQuinticScaling");
}

void kPrintGBVISoftcore(gpuContext gpu, GpuGBVISoftcore* gpuGBVISoftcore, std::string callId, int call)
{
    int maxPrint = 20;
    (void) fprintf( stderr, "kPrintGBVgSoftcore %s %d\n", callId.c_str(), call );
    gpu->psGBVIData->Download();
    gpu->psBornRadii->Download();
    gpu->psBornForce->Download();
    gpu->psPosq4->Download();
    CUDAStream<float>* switchDeriviative = gpuGBVISoftcore-> getSwitchDerivative( );
    switchDeriviative->Download();

    (void) fprintf( stderr, "BornSum Born radii & params\n" );
    for( int ii = 0; ii < gpu->natoms; ii++ ){
        (void) fprintf( stderr, "%6d prm[%14.6e %14.6e %14.6e] bR=%14.6e bF=%14.6e swDrv=%14.6e x[%14.6f %14.6f %14.6f %14.6f]\n",
                        ii,
                        gpu->psGBVIData->_pSysStream[0][ii].x,
                        gpu->psGBVIData->_pSysStream[0][ii].y,
                        gpu->psGBVIData->_pSysStream[0][ii].z,
                        gpu->psBornRadii->_pSysStream[0][ii],
                        gpu->psBornForce->_pSysStream[0][ii],
                        switchDeriviative->_pSysStream[0][ii],
                        gpu->psPosq4->_pSysStream[0][ii].x, gpu->psPosq4->_pSysStream[0][ii].y,
                        gpu->psPosq4->_pSysStream[0][ii].z, gpu->psPosq4->_pSysStream[0][ii].w );
        if( (ii == maxPrint) && ( ii < (gpu->natoms - maxPrint)) ){
            ii = gpu->natoms - maxPrint;
        }
    }
}

void kCalculateGBVISoftcoreBornSum(gpuContext gpu)
{
    //printf("kCalculateGBVIBornSum\n");
    kClearGBVISoftcoreBornSum( gpu );
    LAUNCHERROR("kClearGBVIBornSum from kCalculateGBVISoftcoreBornSum");

    //size_t numWithInteractions;
    switch (gpu->sim.nonbondedMethod)
    {
        case NO_CUTOFF:

            if (gpu->bOutputBufferPerWarp){
                kCalculateGBVISoftcoreN2ByWarpBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit);
            } else {
                kCalculateGBVISoftcoreN2BornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit);
            }
            break;
#if 0
        case CUTOFF:
            if (gpu->bOutputBufferPerWarp)
                kCalculateGBVICutoffByWarpBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateGBVICutoffBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit );
            break;
        case PERIODIC:
            if (gpu->bOutputBufferPerWarp)
                kCalculateGBVIPeriodicByWarpBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit );
            else
                kCalculateGBVIPeriodicBornSum_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit );
            break;
#endif
    }
    LAUNCHERROR("kCalculateGBVISoftcoreBornSum");
}
