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

//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------

#include "amoebaCudaKernels.h"
#include "kCalculateAmoebaCudaUtilities.h"
#include <time.h>

//#define AMOEBA_DEBUG

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;

void SetCalculateAmoebaCudaSASAForcesSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaSASASim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaCudaSASASim: cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

void GetCalculateAmoebaCudaSASAForcesSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaSASASim: cudaMemcpyFromSymbol: SetSim copy from cSim failed");
    status = cudaMemcpyFromSymbol(&amoebaGpu->amoebaSim, cAmoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaCudaSASASim: cudaMemcpyFromSymbol: SetSim copy from cAmoebaSim failed");
}


__device__ void insertionSort( unsigned int length, float* A, int* P)
{
    for( unsigned int ii = 1; ii < length; ii++ ){
        float value = A[ii];
        int index   = P[ii];
        int jj      = ii;
        while( jj > 0 && A[jj-1] > value ){
            A[jj] = A[jj-1];
            P[jj] = P[jj-1];
            jj--;
        }
        A[jj] = value;
        P[jj] = index;
    }
}

__device__ void insertionSortArc( unsigned int length, float4* arcList )
{
    for( unsigned int ii = 1; ii < length; ii++ ){
        float4 value = arcList[ii];
        int jj       = ii;
        while( jj > 0 && arcList[jj-1].x > value.x ){
            arcList[jj] = arcList[jj-1];
            jj--;
        }
        arcList[jj] = value;
    }
}

__device__ int setOmitListOrig( float4 atomCoordI, float4* atomCoord, float* grArray, unsigned int ioListLength, int* ioList,
                                int* omitList, float delta )
{

    float pi_f  = 3.1415926535897f;
    for( unsigned int jj = 0; jj < ioListLength-1; jj++ ){

        int atomJ        = ioList[jj];

        float txj        = atomCoord[atomJ].x - atomCoordI.x;
        float tyj        = atomCoord[atomJ].y - atomCoordI.y;
        float tzj        = atomCoord[atomJ].z - atomCoordI.z;

        float bj         = sqrtf( txj*txj + tyj*tyj + tzj*tzj );
        float therj      = 0.5f*pi_f - asin(min(1.0f,max(-1.0f,grArray[jj])));

        for( unsigned int kk = jj + 1; kk < ioListLength; kk++ ){

            int atomK        = ioList[kk];

            float txk        = atomCoord[atomK].x - atomCoordI.x;
            float tyk        = atomCoord[atomK].y - atomCoordI.y;
            float tzk        = atomCoord[atomK].z - atomCoordI.z;

            float bk         = sqrt( txk*txk + tyk*tyk + tzk*tzk );
            float therk      = 0.5f*pi_f - asin(min(1.0f,max(-1.0f,grArray[kk])));
   
            // check to see if kk circle is intersecting jj circle
            // get distance between circle centers and sum of radii
   
            float cc         = (txj*txk + tyj*tyk + tzj*tzk)/(bj*bk);
                  cc         = acos(min(1.0f,max(-1.0f,cc)));
            float td         = therj + therk;
   
            // check to see if circles enclose separate regions
   
            int mask1        = omitList[jj] == 0 && omitList[kk] == 0 && cc < td ? 1 : 0;
            int mask2        = (cc + therk) < therj ? 1 : 0;
            int mask3        = cc <= delta ? 1 : 0;
            omitList[kk]     = omitList[kk] || (mask1 && (mask2 || mask3)) ? 1 : 0;

            if( (cc > delta) && ( (2.0f*pi_f - cc) <= td) ){
                return 1; // done atom
            }
        }
    }
    return 0;
}

__device__ int setOmitList( float4 atomCoordI, float4* atomCoord, float* grArray, unsigned int ioListLength, int* ioList,
                            int* omitList, float delta )
{

    float pi_f                      = 3.1415926535897f;
    unsigned int outputIncludeCount = 0;
    ioList[outputIncludeCount++]    = ioList[0];
    for( unsigned int jj = 0; jj < ioListLength-1; jj++ ){

        int atomJ        = ioList[jj];

        float txj        = atomCoord[atomJ].x - atomCoordI.x;
        float tyj        = atomCoord[atomJ].y - atomCoordI.y;
        float tzj        = atomCoord[atomJ].z - atomCoordI.z;

        float bj         = txj*txj + tyj*tyj + tzj*tzj;
        float therj      = 0.5f*pi_f - asin(min(1.0f,max(-1.0f,grArray[jj])));

        for( unsigned int kk = jj + 1; kk < ioListLength; kk++ ){

            int atomK        = ioList[kk];

            float txk        = atomCoord[atomK].x - atomCoordI.x;
            float tyk        = atomCoord[atomK].y - atomCoordI.y;
            float tzk        = atomCoord[atomK].z - atomCoordI.z;

            float bk         = txk*txk + tyk*tyk + tzk*tzk;
            float therk      = 0.5f*pi_f - asin(min(1.0f,max(-1.0f,grArray[kk])));
   
            // check to see if kk circle is intersecting jj circle
            // get distance between circle centers and sum of radii
   
            float cc         = (txj*txk + tyj*tyk + tzj*tzk)/sqrtf(bj*bk);
                  cc         = acos(min(1.0f,max(-1.0f,cc)));
            float td         = therj + therk;
   
            // check to see if circles enclose separate regions
   
            int mask1        = omitList[jj] == 0 && omitList[kk] == 0 && cc < td ? 1 : 0;
            int mask2        = (cc + therk) < therj ? 1 : 0;
            int mask3        = cc <= delta ? 1 : 0;
            omitList[kk]     = omitList[kk] || (mask1 && (mask2 || mask3)) ? 1 : 0;

            if( (cc > delta) && ( (2.0f*pi_f - cc) <= td) ){
                return 0; // done atom
            }
        }
        if( omitList[jj+1] == 0 ){
            ioList[outputIncludeCount++] = atomJ;
        }
    }
    return outputIncludeCount;
}

__device__ int setArcList( float4 atomCoordI, float4* atomCoord, float* grArray, unsigned int ioListLength, int* ioList,
                           int* omitList, float delta )
{

    float pi_f                      = 3.1415926535897f;
    unsigned int outputIncludeCount = 0;
    ioList[outputIncludeCount++]    = ioList[0];
    for( unsigned int jj = 0; jj < ioListLength-1; jj++ ){

        int atomJ        = ioList[jj];

        float txj        = atomCoord[atomJ].x - atomCoordI.x;
        float tyj        = atomCoord[atomJ].y - atomCoordI.y;
        float tzj        = atomCoord[atomJ].z - atomCoordI.z;

        float bj         = sqrtf( txj*txj + tyj*tyj + tzj*tzj );
        float therj      = 0.5f*pi_f - asin(min(1.0f,max(-1.0f,grArray[jj])));

        for( unsigned int kk = jj + 1; kk < ioListLength; kk++ ){

            int atomK        = ioList[kk];

            float txk        = atomCoord[atomK].x - atomCoordI.x;
            float tyk        = atomCoord[atomK].y - atomCoordI.y;
            float tzk        = atomCoord[atomK].z - atomCoordI.z;

            float bk         = sqrt( txk*txk + tyk*tyk + tzk*tzk );
            float therk      = 0.5f*pi_f - asin(min(1.0f,max(-1.0f,grArray[kk])));
   
            // check to see if kk circle is intersecting jj circle
            // get distance between circle centers and sum of radii
   
            float cc         = (txj*txk + tyj*tyk + tzj*tzk)/(bj*bk);
                  cc         = acos(min(1.0f,max(-1.0f,cc)));
            float td         = therj + therk;
   
            // check to see if circles enclose separate regions
   
            int mask1        = omitList[jj] == 0 && omitList[kk] == 0 && cc < td ? 1 : 0;
            int mask2        = (cc + therk) < therj ? 1 : 0;
            int mask3        = cc <= delta ? 1 : 0;
            omitList[kk]     = omitList[kk] || (mask1 && (mask2 || mask3)) ? 1 : 0;

            if( (cc > delta) && ( (2.0f*pi_f - cc) <= td) ){
                return 0; // done atom
            }
        }
        if( omitList[jj+1] == 0 ){
            ioList[outputIncludeCount++] = atomJ;
        }
    }
    return outputIncludeCount;
}

void insertionSortCpu( unsigned int length, float* A, int* P)
{
    for( unsigned int ii = 1; ii < length; ii++ ){
        float value = A[ii];
        int index   = P[ii];
        int jj      = ii;
        while( jj > 0 && A[jj-1] > value ){
            A[jj] = A[jj-1];
            P[jj] = P[jj-1];
            jj--;
        }
        A[jj] = value;
        P[jj] = index;
    }
}

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_NONBOND_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_NONBOND_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_NONBOND_THREADS_PER_BLOCK, 1)
#endif
void kCalculateAmoebaCudaSASA( float4* atomCoord, float* rPlusProbe,
                                          int* doneAtom, int* ioListCount,
                                          int* workI1Array, int* workI2Array, int* workI3Array, int* workI4Array,
                                          float* workFArray
#ifdef AMOEBA_DEBUG
                                          ,float4* debugArray
#endif
)
{


    unsigned int atomI           = blockIdx.x*blockDim.x+threadIdx.x;
    if( atomI >= cAmoebaSim.numberOfAtoms)
    {
        return;
    }

    unsigned int arrayOffset     = atomI*cAmoebaSim.maxarc;

    float* grArray               = workFArray  + arrayOffset;
    int* ioList                  = workI1Array + arrayOffset;
    int* omitList                = workI2Array + arrayOffset;

    unsigned int atomJ           = 0;
    float4 coordI                = atomCoord[atomI];
    float  delta                 = 1.0e-08f;
    float  rPlusProbeI           = rPlusProbe[atomI];
    unsigned int io              = 0;
    while (atomJ < cAmoebaSim.numberOfAtoms)
    {
        float4 coordJ            = atomCoord[atomJ];
        float  rPlusProbeJ       = rPlusProbe[atomJ];
       
        float rplus              = rPlusProbeI + rPlusProbeJ;
        float tx                 = coordJ.x - coordI.x;
        float ty                 = coordJ.y - coordI.y;
        float tz                 = coordJ.z - coordI.z;
        float ccsq               = tx*tx + ty*ty + tz*tz;
        float cc                 = sqrtf( ccsq );
        float rminus             = rPlusProbeI - rPlusProbeJ;

        int mask1                = (atomI != atomJ) && (fabsf( tx ) < rplus) && (fabsf( ty ) < rplus) && (fabsf( tz ) < rplus) ? 1 : 0;
        int mask2                = ((rplus-cc) > delta) ? mask1 : 0;
        int mask3                = (cc-fabsf(rminus) > delta) ? mask2 : 0; 

        if( mask3 ){

            ioList[io]           = atomJ;
            omitList[io]         = 0;
            grArray[io]          = (ccsq+rplus*rminus) / (2.0f*rPlusProbeI*cc);
            io++;

#ifdef AMOEBA_DEBUG
            debugArray[arrayOffset+io].x = tx;
            debugArray[arrayOffset+io].y = ty;
            debugArray[arrayOffset+io].z = rplus;
            debugArray[arrayOffset+io].w = rPlusProbe[atomJ];
#endif
        }
        atomJ++;
    }

//    ioListCount[atomI] = io;
    insertionSort( io, grArray, ioList );

    // set omit/doneAtom values

    io              = setOmitList( coordI, atomCoord, grArray, io, ioList, omitList, delta );
    doneAtom[atomI] = (io == 0) ? 1 : 0;
}


/**---------------------------------------------------------------------------------------

   Compute SASA force/energy

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

void kCalculateAmoebaSASAForces(amoebaGpuContext amoebaGpu )
{

   // ---------------------------------------------------------------------------------------

#ifdef AMOEBA_DEBUG
    static const char* methodName = "kCalculateAmoebaSASAForces";
    static int timestep = 0;
    std::vector<int> fileId;
    timestep++;
    fileId.resize( 2 );
    fileId[0] = timestep;
    fileId[1] = 1;
#endif

   if( amoebaGpu )return;

   // ---------------------------------------------------------------------------------------

     gpuContext gpu    = amoebaGpu->gpuContext;
    int numOfElems     = amoebaGpu->paddedNumberOfAtoms;
    int numThreads     = min( THREADS_PER_BLOCK, numOfElems );
    int numBlocks      = numOfElems/numThreads;

    if( (numOfElems % numThreads) != 0 )numBlocks++;

#ifdef AMOEBA_DEBUG
    CUDAStream<float4>* debugArray = NULL;
    if( amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "%s %d numOfElems=%d numThreads=%d numBlocks=%d maxarc=%d\n",
                        methodName, gpu->natoms, numOfElems, numThreads, numBlocks, amoebaGpu->amoebaSim.maxarc );
        (void) fflush( amoebaGpu->log );
        int paddedNumberOfAtoms                    = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
        debugArray                                 = new CUDAStream<float4>(paddedNumberOfAtoms*paddedNumberOfAtoms, 1, "DebugArray");
        memset( debugArray->_pSysStream[0],      0, sizeof( float )*4*paddedNumberOfAtoms*paddedNumberOfAtoms);
        debugArray->Upload();

    }   
#endif

    // initialize [ InducedDipole & InducedDipolePolar ] to [ E_Field & E_FieldPolar]*Polarizability

    kCalculateAmoebaCudaSASA<<< numBlocks, numThreads >>>(
         gpu->psPosq4->_pDevStream[0],
         amoebaGpu->psSASA_Radii->_pDevStream[0],
         amoebaGpu->psDoneAtom->_pDevStream[0],
         amoebaGpu->psIoListCount->_pDevStream[0],
         amoebaGpu->psIntWorkArray[0]->_pDevStream[0],
         amoebaGpu->psIntWorkArray[1]->_pDevStream[0],
         amoebaGpu->psIntWorkArray[2]->_pDevStream[0],
         amoebaGpu->psIntWorkArray[3]->_pDevStream[0],
         amoebaGpu->psFloatWorkArray->_pDevStream[0]
#ifdef AMOEBA_DEBUG
, debugArray->_pDevStream[0]
#endif
 );
    LAUNCHERROR("kCalculateAmoebaCudaSASA");

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){

        debugArray->Download();
        amoebaGpu->psIntWorkArray[0]->Download();
        amoebaGpu->psIntWorkArray[1]->Download();
        amoebaGpu->psIoListCount->Download();
        amoebaGpu->psDoneAtom->Download();
        amoebaGpu->psFloatWorkArray->Download();
        (void) fprintf( amoebaGpu->log, "Post kCalculateAmoebaCudaSASA \n" );
        for( unsigned int ii = 0; ii < amoebaGpu->paddedNumberOfAtoms; ii++)
        {
            int index            = ii*amoebaGpu->amoebaSim.maxarc;
            unsigned int length  = (unsigned int) amoebaGpu->psIoListCount->_pSysStream[0][ii];
            int omitSum = 0;
            for( unsigned int jj = 0; jj < length; jj++)
            {
                omitSum += amoebaGpu->psIntWorkArray[1]->_pSysStream[0][index+jj];
            }
            (void) fprintf( amoebaGpu->log, "%5u %5d omitSum=%5d done=%d\n", ii, length, omitSum, amoebaGpu->psDoneAtom->_pSysStream[0][ii] );
            if( ii == 0 ){
                for( unsigned int jj = 0; jj < length; jj++)
                {
                    (void) fprintf( amoebaGpu->log, "    %5d omt=%5d %14.6e %14.6e %14.6e %14.6e %14.6e \n",
                                    amoebaGpu->psIntWorkArray[0]->_pSysStream[0][index+jj],
                                    amoebaGpu->psIntWorkArray[1]->_pSysStream[0][index+jj],
                                    amoebaGpu->psFloatWorkArray->_pSysStream[0][index+jj],
                                    debugArray->_pSysStream[0][index+jj].x, debugArray->_pSysStream[0][index+jj].y, debugArray->_pSysStream[0][index+jj].z, debugArray->_pSysStream[0][index+jj].w);
                }
/*
                int*   permutation   = (int*) malloc( sizeof( int )*length );
                float* values        = (float*) malloc( sizeof( float )*length );
                for( unsigned int jj = 1; jj <= length; jj++)
                {
                    permutation[jj-1] = amoebaGpu->psIntWorkArray[0]->_pSysStream[0][index+jj];
                    values[jj-1]      = amoebaGpu->psFloatWorkArray->_pSysStream[0][index+jj];
                }
                insertionSortCpu( length, values, permutation);
                for( unsigned int jj = 0; jj < length; jj++)
                {
                    (void) fprintf( amoebaGpu->log, "Sortd: %5u %5d %14.6e\n", jj, permutation[jj], values[jj] );
                }
                free( permutation );
                free( values );
*/
            }
        }
        (void) fflush( amoebaGpu->log );
    }   
#endif
 
    if( 0 )
    {
        unsigned int iterations = 5000;
        time_t start            = clock();
        for( unsigned int ii = 0; ii < iterations; ii++)
        {
            kCalculateAmoebaCudaSASA<<< numBlocks, numThreads >>>(
                 gpu->psPosq4->_pDevStream[0],
                 amoebaGpu->psSASA_Radii->_pDevStream[0],
                 amoebaGpu->psDoneAtom->_pDevStream[0],
                 amoebaGpu->psIoListCount->_pDevStream[0],
                 amoebaGpu->psIntWorkArray[0]->_pDevStream[0],
                 amoebaGpu->psIntWorkArray[1]->_pDevStream[0],
                 amoebaGpu->psIntWorkArray[2]->_pDevStream[0],
                 amoebaGpu->psIntWorkArray[3]->_pDevStream[0],
                 amoebaGpu->psFloatWorkArray->_pDevStream[0]
#ifdef AMOEBA_DEBUG
                , debugArray->_pDevStream[0]
#endif
                );
                LAUNCHERROR("kCalculateAmoebaCudaSASA");
        }
        float sasaTime = (float)(clock()-start)/CLOCKS_PER_SEC;
        (void) fprintf( amoebaGpu->log, "Total sasa kernel time=%12.5e time/iter=%12.5e %u %u\n", sasaTime, sasaTime/(float) iterations, (unsigned int) start, (unsigned int) clock()); 

    }
#ifdef AMOEBA_DEBUG
    delete debugArray;
#endif

}
