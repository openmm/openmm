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
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "GpuNonbondedSoftcore.h"
#include "GpuFreeEnergyCudaKernels.h"
#include "openmm/OpenMMException.h"
#include <algorithm>

// structure containing array of softcore lambdas

struct cudaFreeEnergySimulationNonBonded {
    float* pParticleSoftCoreLJLambda;
};
struct cudaFreeEnergySimulationNonBonded feSim;

// device handles

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaFreeEnergySimulationNonBonded feSimDev;

// write address of structs to devices

void SetCalculateCDLJSoftcoreGpuSim( gpuContext gpu )
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");

    //(void) fprintf( stderr, "SetCalculateCDLJSoftcoreGpuSim gpu=%p cSim=%p sizeof=%u\n", gpu, &gpu->sim, sizeof(cudaGmxSimulation) ); fflush( stderr );
}

void SetCalculateCDLJSoftcoreSupplementarySim( float* gpuParticleSoftCoreLJLambda)
{
    cudaError_t status;
    feSim.pParticleSoftCoreLJLambda = gpuParticleSoftCoreLJLambda;
    status = cudaMemcpyToSymbol(feSimDev, &feSim, sizeof(cudaFreeEnergySimulationNonBonded));
    RTERROR(status, "cudaMemcpyToSymbol: SetCalculateCDLJSoftcoreSupplementarySim");

    //(void) fprintf( stderr, "SetCalculateCDLJSoftcoreSupplementarySim\n" );
}

void GetCalculateCDLJSoftcoreForcesSim(float* gpuParticleSoftCoreLJLambda)
{
//    cudaError_t status;
//    status = cudaMemcpyFromSymbol(gpuParticleSoftCoreLJLambda, particleSoftCoreLJLambdaDev, sizeof(float*));
//    RTERROR(status, "cudaMemcpyFromSymbol: GetCalculateCDLJSoftcoreForcesSim failed");
}

// create, initialize and entrt SoftCoreLJLambda values
// return handle to GpuNonbondedSoftcore object

static void setSoftcoreExclusions(gpuContext gpu, const std::vector<std::vector<int> >& exclusions) {
    if (gpu->exclusions.size() > 0) { 
        bool ok = (exclusions.size() == gpu->exclusions.size());
        for (unsigned int i = 0; i < exclusions.size() && ok; i++) {
            if (exclusions[i].size() != gpu->exclusions[i].size())
                ok = false;
            else {
                for (unsigned int j = 0; j < exclusions[i].size(); j++) 
                    if (find(gpu->exclusions[i].begin(), gpu->exclusions[i].end(), exclusions[i][j]) == gpu->exclusions[i].end())
                        ok = false;
            }
        }
        if (!ok)
            throw OpenMM::OpenMMException("All nonbonded forces must have identical sets of exceptions");
    }    
    gpu->exclusions = exclusions;
}

extern "C"
GpuNonbondedSoftcore* gpuSetNonbondedSoftcoreParameters(gpuContext gpu, float epsfac, const std::vector<int>& atom, const std::vector<float>& c6,
                                                        const std::vector<float>& c12, const std::vector<float>& q,
                                                        const std::vector<float>& softcoreLJLambdaArray, const std::vector<char>& symbol,
                                                        const std::vector<std::vector<int> >& exclusions, CudaNonbondedMethod method)
{
    unsigned int numberOfParticles     = c6.size();
    gpu->sim.epsfac                    = epsfac;
    gpu->sim.nonbondedMethod           = method;
    if (numberOfParticles > 0)
        setSoftcoreExclusions(gpu, exclusions);
    
    // create gpuNonbondedSoftcore

    GpuNonbondedSoftcore* gpuNonbondedSoftcore = new GpuNonbondedSoftcore();
    gpuNonbondedSoftcore->initializeParticleSoftCoreLJLambda( numberOfParticles );
    float minSoftcore                          = 1.0e+10;
    for (unsigned int i = 0; i < numberOfParticles; i++)
    {
            float p0               = q[i];

            // track min softcore value

            float softcoreLJLambda = softcoreLJLambdaArray[i];
            if( minSoftcore > softcoreLJLambda ){
                minSoftcore = softcoreLJLambda;
            }
            gpuNonbondedSoftcore->setParticleSoftCoreLJLambda( i, softcoreLJLambda );

            float p1 = 0.5f, p2 = 0.0f;               
            if ((c6[i] > 0.0f) && (c12[i] > 0.0f))
            {
                p1 = 0.5f * pow(c12[i] / c6[i], 1.0f / 6.0f);
                p2 = c6[i] * sqrt(1.0f / c12[i]);
            }
            if (symbol.size() > 0)
                gpu->pAtomSymbol[i] = symbol[i];

            (*gpu->psPosq4)[i].w          = p0;
            (*gpu->psSigEps2)[i].x        = p1;
            (*gpu->psSigEps2)[i].y        = p2;
    }
    gpuNonbondedSoftcore->setSoftCoreLJLambda( minSoftcore );

    // Dummy out extra atom data
    for (unsigned int i = numberOfParticles; i < gpu->sim.paddedNumberOfAtoms; i++)
    {
        (*gpu->psPosq4)[i].x                = 100000.0f + i * 10.0f;
        (*gpu->psPosq4)[i].y                = 100000.0f + i * 10.0f;
        (*gpu->psPosq4)[i].z                = 100000.0f + i * 10.0f;
        (*gpu->psPosq4)[i].w                = 0.0f;
        (*gpu->psSigEps2)[i].x              = 0.0f;
        (*gpu->psSigEps2)[i].y              = 0.0f;
    }

#undef DUMP_PARAMETERS
#define DUMP_PARAMETERS 0
#if (DUMP_PARAMETERS == 1)
    (void) fprintf( stderr,"gpuSetNonbondedSoftcoreParameters: %5u epsfac=%14.7e method=%d\n", numberOfParticles, gpu->sim.paddedNumberOfAtoms, epsfac, method );
    int maxPrint = 31;
    for (unsigned int ii = 0; ii < gpu->sim.paddedNumberOfAtoms; ii++){
        (void) fprintf( stderr,"%6u x[%14.7e %14.7e %14.7e %14.7e] sig[%14.7e %14.7e]\n",
                        ii, (*gpu->psPosq4)[ii].x, (*gpu->psPosq4)[ii].y, (*gpu->psPosq4)[ii].z, (*gpu->psPosq4)[ii].w,
                        (*gpu->psSigEps2)[ii].x, (*gpu->psSigEps2)[ii].y );
        if( ii == maxPrint && ii < gpu->sim.paddedNumberOfAtoms - maxPrint ){
           ii = gpu->sim.paddedNumberOfAtoms - maxPrint;
        }
    }
#endif
 
    // upload data to board

    gpuNonbondedSoftcore->upload( gpu );

    gpu->psPosq4->Upload();
    gpu->psSigEps2->Upload();

    return gpuNonbondedSoftcore;
}

// delete gpuNonbondedSoftcore

extern "C"
void gpuDeleteNonbondedSoftcoreParameters( void* gpuNonbondedSoftcore)
{
    GpuNonbondedSoftcore* internalGNonbondedSoftcore = static_cast<GpuNonbondedSoftcore*>(gpuNonbondedSoftcore);
    delete internalGNonbondedSoftcore;
}

extern "C"
bool gpuIsAvailableSoftcore()
{
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    return (deviceCount > 0);
}

struct Atom {
    float x;
    float y;
    float z;
    float q;
    float sig;
    float eps;
    float softCoreLJLambda;
    float fx;
    float fy;
    float fz;
};

#if 0
texture<float, 1, cudaReadModeElementType> tabulatedErfcRef;

__device__ float fastErfc(float r)
{
    float normalized = cSim.tabulatedErfcScale*r;
    int index = (int) normalized;
    float fract2 = normalized-index;
    float fract1 = 1.0f-fract2;
    return fract1*tex1Dfetch(tabulatedErfcRef, index) + fract2*tex1Dfetch(tabulatedErfcRef, index+1);
}

// Include versions of the kernels for N^2 calculations.

#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateNonbondedSoftcore.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateNonbondedSoftcore.h"

#endif

// Include versions of the kernels for N^2 calculations with softcore LJ.

#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2SoftcoreLJ##b
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_SOFTCORE_LJ
#include "kCalculateNonbondedSoftcore.h"

#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2SoftcoreLJByWarp##b
#include "kCalculateNonbondedSoftcore.h"
#undef USE_SOFTCORE_LJ

// Include versions of the kernels with cutoffs.

#if 0
#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_CUTOFF
#define METHOD_NAME(a, b) a##Cutoff##b
#include "kCalculateNonbondedSoftcore.h"
#include "kFindInteractingBlocks.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##CutoffByWarp##b
#include "kCalculateNonbondedSoftcore.h"

// Include versions of the kernels with periodic boundary conditions.

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_PERIODIC
#define METHOD_NAME(a, b) a##Periodic##b
#include "kCalculateNonbondedSoftcore.h"
#include "kFindInteractingBlocks.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##PeriodicByWarp##b
#include "kCalculateNonbondedSoftcore.h"

// Include versions of the kernels for Ewald

#undef METHOD_NAME
#undef USE_OUTPUT_BUFFER_PER_WARP
#define USE_PERIODIC
#define USE_EWALD
#define METHOD_NAME(a, b) a##Ewald##b
#include "kCalculateNonbondedSoftcore.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##EwaldByWarp##b
#include "kCalculateNonbondedSoftcore.h"

// Reciprocal Space Ewald summation is in a separate kernel
#include "kCalculateCDLJEwaldFastReciprocal.h"

void kCalculatePME(gpuContext gpu);
#endif

void kCalculateCDLJSoftcoreForces(gpuContext gpu )
{

    //printf("kCalculateCDLJCutoffForces %d\n", gpu->sim.nonbondedMethod); fflush( stdout );
    switch (gpu->sim.nonbondedMethod)
    {
        case NO_CUTOFF:
           if (gpu->bOutputBufferPerWarp)
                kCalculateCDLJSoftcoreN2SoftcoreLJByWarpForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                         sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit);
           else
                   kCalculateCDLJSoftcoreN2SoftcoreLJForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                           sizeof(Atom)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pWorkUnit );
//(gpu->sim.pWorkUnit, gpuNonbondedSoftcore->getGpuParticleSoftCoreLJLambda());
            LAUNCHERROR("kCalculateCDLJSoftcoreN2Forces");

#if 0
int maxPrint = 31; 
gpu->psWorkUnit->Download();
fprintf( stderr, "kCalculateCDLJSoftcoreForces: bOutputBufferPerWarp=%u blks=%u th/blk=%u wu=%u %u shrd=%u\n", gpu->bOutputBufferPerWarp,
                 gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block, gpu->sim.workUnits, gpu->psWorkUnit->_pSysStream[0][0],
        sizeof(Atom)*gpu->sim.nonbond_threads_per_block );

               gpu->psPosq4->Download();

                (void) fprintf( stderr, "\nkCalculateGBVISoftcoreBornSum: pre BornSum %s Born radii & params\n",
                               (gpu->bIncludeGBVI ? "GBVI" : "Obc") );
                for( int ii = 0; ii < gpu->natoms; ii++ ){
                   (void) fprintf( stderr, "%6d bSum=%14.6e param[%14.6e %14.6e %14.6e] x[%14.6f %14.6f %14.6f %14.6f]\n",
                                   ii,
                                   gpu->psBornSum->_pSysStream[0][ii],
                                   gpu->psGBVIData->_pSysStream[0][ii].x,
                                   gpu->psGBVIData->_pSysStream[0][ii].y,
                                   gpu->psGBVIData->_pSysStream[0][ii].z,
                                   gpu->psPosq4->_pSysStream[0][ii].x, gpu->psPosq4->_pSysStream[0][ii].y,
                                   gpu->psPosq4->_pSysStream[0][ii].z, gpu->psPosq4->_pSysStream[0][ii].w
                                 );
                   if( (ii == maxPrint) && ( ii < (gpu->natoms - maxPrint)) ){
                      ii = gpu->natoms - maxPrint;
                   }
                }

#endif
#undef GBVI


            break;
#if 0
        case CUTOFF:
            kFindBlockBoundsCutoff_kernel<<<(gpu->psGridBoundingBox->_length+63)/64, 64>>>();
            LAUNCHERROR("kFindBlockBoundsCutoff");
            kFindBlocksWithInteractionsCutoff_kernel<<<gpu->sim.interaction_blocks, gpu->sim.interaction_threads_per_block>>>();
            LAUNCHERROR("kFindBlocksWithInteractionsCutoff");
            compactStream(gpu->compactPlan, gpu->sim.pInteractingWorkUnit, gpu->sim.pWorkUnit, gpu->sim.pInteractionFlag, gpu->sim.workUnits, gpu->sim.pInteractionCount);
            kFindInteractionsWithinBlocksCutoff_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                    sizeof(unsigned int)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
#if 0
    static int iteration = 0;
    if (iteration >= 0)
    {
        gpu->psInteractingWorkUnit->Download();
        gpu->psInteractionCount->Download();
/*
    unsigned int totalWarps = cSim.nonbond_blocks*cSim.nonbond_threads_per_block/GRID;
    unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    unsigned int numWorkUnits = cSim.pInteractionCount[0];
    unsigned int pos = warp*numWorkUnits/totalWarps;
    unsigned int end = (warp+1)*numWorkUnits/totalWarps;
*/
 
        printf("# Post kCalculateCDLJCutoffForces %d atoms warps=%d cnt=%u bOutputBufferPerWarp=%d zC=%d\n", 
                gpu->natoms, ((gpu->sim.nonbond_blocks*gpu->sim.nonbond_threads_per_block)/GRID),
                gpu->psInteractionCount->_pSysStream[0][0], gpu->bOutputBufferPerWarp,
                (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block);
        fflush( stdout );
        for (int i = 0; i < gpu->psInteractingWorkUnit->_stride; i++)
        {
            printf("%5d %u\n", i, gpu->psInteractingWorkUnit->_pSysStream[0][i] );
            fflush( stdout );
        }
    }
    iteration++;
#endif


            if (gpu->bOutputBufferPerWarp)
                kCalculateCDLJCutoffByWarpForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateCDLJCutoffForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            LAUNCHERROR("kCalculateCDLJCutoffForces");
            break;
        case PERIODIC:
            kFindBlockBoundsPeriodic_kernel<<<(gpu->psGridBoundingBox->_length+63)/64, 64>>>();
            LAUNCHERROR("kFindBlockBoundsPeriodic");
            kFindBlocksWithInteractionsPeriodic_kernel<<<gpu->sim.interaction_blocks, gpu->sim.interaction_threads_per_block>>>();
            LAUNCHERROR("kFindBlocksWithInteractionsPeriodic");
            compactStream(gpu->compactPlan, gpu->sim.pInteractingWorkUnit, gpu->sim.pWorkUnit, gpu->sim.pInteractionFlag, gpu->sim.workUnits, gpu->sim.pInteractionCount);
            kFindInteractionsWithinBlocksPeriodic_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                    sizeof(unsigned int)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            if (gpu->bOutputBufferPerWarp)
                kCalculateCDLJPeriodicByWarpForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateCDLJPeriodicForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            LAUNCHERROR("kCalculateCDLJPeriodicForces");
            break;
        case EWALD:
        case PARTICLE_MESH_EWALD:
            kFindBlockBoundsPeriodic_kernel<<<(gpu->psGridBoundingBox->_length+63)/64, 64>>>();
            LAUNCHERROR("kFindBlockBoundsPeriodic");
            kFindBlocksWithInteractionsPeriodic_kernel<<<gpu->sim.interaction_blocks, gpu->sim.interaction_threads_per_block>>>();
            LAUNCHERROR("kFindBlocksWithInteractionsPeriodic");
            compactStream(gpu->compactPlan, gpu->sim.pInteractingWorkUnit, gpu->sim.pWorkUnit, gpu->sim.pInteractionFlag, gpu->sim.workUnits, gpu->sim.pInteractionCount);
            kFindInteractionsWithinBlocksPeriodic_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                    sizeof(unsigned int)*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            LAUNCHERROR("kFindInteractionsWithinBlocksPeriodic");
            cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
            cudaBindTexture(NULL, &tabulatedErfcRef, gpu->psTabulatedErfc->_pDevData, &channelDesc, gpu->psTabulatedErfc->_length*sizeof(float));
            if (gpu->bOutputBufferPerWarp)
                kCalculateCDLJEwaldByWarpForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            else
                kCalculateCDLJEwaldForces_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block,
                        (sizeof(Atom)+sizeof(float3))*gpu->sim.nonbond_threads_per_block>>>(gpu->sim.pInteractingWorkUnit);
            LAUNCHERROR("kCalculateCDLJEwaldForces");
            if (gpu->sim.nonbondedMethod == EWALD)
            {
                kCalculateEwaldFastCosSinSums_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.nonbond_threads_per_block>>>();
                LAUNCHERROR("kCalculateEwaldFastCosSinSums");
                kCalculateEwaldFastForces_kernel<<<gpu->sim.blocks, gpu->sim.update_threads_per_block>>>();
                LAUNCHERROR("kCalculateEwaldFastForces");
            }
            else
                kCalculatePME(gpu);
#endif
    }
}

void kPrintForces(gpuContext gpu, std::string idString, int call )
{
 //   printf("kReduceForces\n");
#define GBVI_DEBUG 0
#if ( GBVI_DEBUG == 4 )

                gpu->psBornRadii->Download();
                gpu->psObcData->Download();
                gpu->psObcChain->Download();
                gpu->psBornForce->Download();
                gpu->psForce4->Download();
                gpu->psPosq4->Download();
                int maxPrint = 0; 
int   nanHit       = 0;
int   targetIndex  = 852;
float maxForce     = 3.0e+04;
float maxPosition  = 2.0e+02;
                for( int ii = 0; ii < gpu->natoms; ii++ ){

int   hit  = 0;
float dist = sqrtf( gpu->psPosq4->_pSysStream[0][ii].x*gpu->psPosq4->_pSysStream[0][ii].x + 
                    gpu->psPosq4->_pSysStream[0][ii].y*gpu->psPosq4->_pSysStream[0][ii].y +
                    gpu->psPosq4->_pSysStream[0][ii].z*gpu->psPosq4->_pSysStream[0][ii].z );

if( fabs( gpu->psForce4->_pSysStream[0][ii].x ) > maxForce ||
    fabs( gpu->psForce4->_pSysStream[0][ii].y ) > maxForce ||
    fabs( gpu->psForce4->_pSysStream[0][ii].z ) > maxForce ||
//    gpu->psBornRadii->_pSysStream[0][ii] <= 0.0            ||
    dist > maxPosition                                     ||
    isnan( gpu->psForce4->_pSysStream[0][ii].x )           ||
    isnan( gpu->psForce4->_pSysStream[0][ii].y )           ||
    isnan( gpu->psForce4->_pSysStream[0][ii].z )  ){  
   hit = 1;
} else {
   hit = 0;
}
if( ii == targetIndex || ii == (targetIndex+1) || ii == (targetIndex+2) )hit = 1;
if( isnan( gpu->psBornForce->_pSysStream[0][ii] ) ||
    isnan( gpu->psBornRadii->_pSysStream[0][ii] ) ||
    isnan( gpu->psObcChain->_pSysStream[0][ii]  ) ||
    isnan( gpu->psForce4->_pSysStream[0][ii].x  )  ||
    isnan( gpu->psForce4->_pSysStream[0][ii].y  )  ||
    isnan( gpu->psForce4->_pSysStream[0][ii].z  )  ){  
   hit    = 1;
   nanHit = 1;
}

                //if( hit || ii < maxPrint || ii >= (gpu->natoms - maxPrint) )
                if( hit ){
                    static int firstHit = 1;
                    if( firstHit ){
                       firstHit = 0;
                       (void) fprintf( stderr, "\nkPrintForces: %d [r, scl q] b[r/c/f] f[] x[] Born radii/force (%p %p)\n", call,
                                       gpu->psBornForce, gpu->psBornForce->_pDevStream[0] );
                    }
                    (void) fprintf( stderr, "%6d [%8.3f %8.3f %8.3f] b[%13.6e %13.6e %13.6e] f[%13.6e %13.6e %13.6e] x[%13.6e %13.6e %13.6e] %10.3e %s %s %d\n",
                                    ii, 
                                    (gpu->psObcData->_pSysStream[0][ii].x + 0.009f),
                                    (gpu->psObcData->_pSysStream[0][ii].y/gpu->psObcData->_pSysStream[0][ii].x),
                                    gpu->psPosq4->_pSysStream[0][ii].w,
                                    gpu->psBornRadii->_pSysStream[0][ii],
                                    gpu->psObcChain->_pSysStream[0][ii],
                                    gpu->psBornForce->_pSysStream[0][ii],

                                    gpu->psForce4->_pSysStream[0][ii].x,
                                    gpu->psForce4->_pSysStream[0][ii].y,
                                    gpu->psForce4->_pSysStream[0][ii].z,

                                    gpu->psPosq4->_pSysStream[0][ii].x,
                                    gpu->psPosq4->_pSysStream[0][ii].y,
                                    gpu->psPosq4->_pSysStream[0][ii].z, dist,

                                    (hit ? "XXXXXXX" : "" ), idString.c_str(), call );
                }
                }
                (void) fflush( stderr );
//                if( nanHit )exit(0);
#endif

}

