#ifndef __AMOEBA_GPU_TYPES_H__
#define __AMOEBA_GPU_TYPES_H__

/* -------------------------------------------------------------------------- *
 *                             OpenMMAmoeba                                   *
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

#include "amoebaGpuTypes.h"

#include <string>
#include <vector>

typedef std::vector<std::string> StringVector;
typedef std::vector<StringVector> StringVectorVector;

#define SQRT sqrtf
#define EXP  expf
#define DOT3(u,v) ((u[0])*(v[0]) + (u[1])*(v[1]) + (u[2])*(v[2]))

typedef std::vector<std::vector<double> > VectorOfDoubleVectors;

// local (bond) forces

extern void SetCalculateAmoebaLocalForcesSim(amoebaGpuContext gpu);
extern void GetCalculateAmoebaLocalForcesSim(amoebaGpuContext gpu);
extern void kCalculateAmoebaLocalForces(amoebaGpuContext gpu);

// multipole

extern void SetCalculateAmoebaMultipoleForcesSim(amoebaGpuContext gpu);
extern void GetCalculateAmoebaMultipoleForcesSim(amoebaGpuContext gpu);
extern void kCalculateAmoebaMultipoleForces(amoebaGpuContext amoebaGpu, bool performGk );

// vdw

extern void SetCalculateAmoebaCudaVdw14_7Sim(amoebaGpuContext gpu);
extern void GetCalculateAmoebaCudaVdw14_7Sim(amoebaGpuContext gpu);
extern void kCalculateAmoebaVdw14_7Forces(amoebaGpuContext amoebaGpu, int applyCutoff );

// wca dispersion

extern void SetCalculateAmoebaCudaWcaDispersionSim(amoebaGpuContext gpu);
extern void GetCalculateAmoebaCudaWcaDispersionSim(amoebaGpuContext gpu);
extern void kCalculateAmoebaWcaDispersionForces(amoebaGpuContext amoebaGpu );

// fixed electric field -- no cutoff

extern void SetCalculateAmoebaCudaFixedEFieldSim(amoebaGpuContext gpu);
extern void GetCalculateAmoebaCudaFixedEFieldSim(amoebaGpuContext gpu);
extern void cudaComputeAmoebaFixedEField( amoebaGpuContext gpu);

// fixed electric field  -- PME
extern void SetCalculateAmoebaCudaPmeFixedEFieldSim(amoebaGpuContext gpu);
extern void GetCalculateAmoebaCudaPmeFixedEFieldSim(amoebaGpuContext gpu);
extern void cudaComputeAmoebaPmeFixedEField( amoebaGpuContext gpu);

// fixed electric field and Gk

extern void SetCalculateAmoebaCudaFixedEAndGKFieldsSim(amoebaGpuContext gpu);
extern void GetCalculateAmoebaCudaFixedEAndGKFieldsSim(amoebaGpuContext gpu);
extern void cudaComputeAmoebaFixedEAndGkFields( amoebaGpuContext gpu);

// mutual induced 

extern void SetCalculateAmoebaCudaMutualInducedFieldSim(amoebaGpuContext gpu);
extern void GetCalculateAmoebaCudaMutualInducedFieldSim(amoebaGpuContext gpu);
extern void cudaComputeAmoebaMutualInducedField( amoebaGpuContext gpu);

extern void SetCalculateAmoebaCudaPmeMutualInducedFieldSim(amoebaGpuContext gpu);
extern void GetCalculateAmoebaCudaPmeMutualInducedFieldSim(amoebaGpuContext gpu);
extern void cudaComputeAmoebaPmeMutualInducedField( amoebaGpuContext gpu);

// mutual induced and Gk

extern void SetCalculateAmoebaCudaMutualInducedAndGkFieldsSim(amoebaGpuContext amoebaGpu);
extern void GetCalculateAmoebaCudaMutualInducedAndGkFieldsSim(amoebaGpuContext amoebaGpu);
extern void cudaComputeAmoebaMutualInducedAndGkField( amoebaGpuContext gpu);

extern void cudaComputeAmoebaLabFrameMoments( amoebaGpuContext amoebaGpu );
extern void cudaWriteFloat4AndFloat1ArraysToFile( int numberOfAtoms, const std::string& fname, int timestep, int entriesPerAtom1, CUDAStream<float4>* array1, 
                                                  int entriesPerAtom2, CUDAStream<float>* array2 );

extern void SetCalculateAmoebaElectrostaticSim( amoebaGpuContext amoebaGpu );
extern void GetCalculateAmoebaElectrostaticSim( amoebaGpuContext amoebaGpu );
extern void cudaComputeAmoebaElectrostatic( amoebaGpuContext amoebaGpu, int addTorqueToForce );

extern void SetCalculateAmoebaPmeDirectElectrostaticSim( amoebaGpuContext amoebaGpu );
extern void GetCalculateAmoebaPmeDirectElectrostaticSim( amoebaGpuContext amoebaGpu );
extern void cudaComputeAmoebaPmeElectrostatic( amoebaGpuContext amoebaGpu );

extern void SetCalculateAmoebaCudaMapTorquesSim(amoebaGpuContext gpu);
extern void GetCalculateAmoebaCudaMapTorquesSim(amoebaGpuContext gpu);
extern void cudaComputeAmoebaMapTorqueAndAddToForce( amoebaGpuContext gpu, CUDAStream<float>* psTorque );

extern void SetCalculateAmoebaKirkwoodSim( amoebaGpuContext amoebaGpu );
extern void GetCalculateAmoebaKirkwoodSim( amoebaGpuContext amoebaGpu );
//extern void cudaComputeAmoebaKirkwood( amoebaGpuContext amoebaGpu );
extern void kCalculateAmoebaKirkwood( amoebaGpuContext amoebaGpu );

extern void SetCalculateAmoebaKirkwoodEDiffSim( amoebaGpuContext amoebaGpu );
extern void GetCalculateAmoebaKirkwoodEDiffSim( amoebaGpuContext amoebaGpu );
//extern void cudaComputeAmoebaKirkwoodEDiff( amoebaGpuContext amoebaGpu );
extern void kCalculateAmoebaKirkwoodEDiff( amoebaGpuContext amoebaGpu );

extern void SetCalculateAmoebaObcGbsaBornSumSim( gpuContext gpu );
extern void GetCalculateAmoebaObcGbsaBornSumSim( gpuContext gpu );
extern void cudaComputeAmoebaBornRadii( amoebaGpuContext amoebaGpu );

// OBC -- Part 1
//extern void SetCalculateObcGbsaForces1Sim(gpuContext gpu);
//extern void GetCalculateObcGbsaForces1Sim(gpuContext gpu);
//extern void kCalculateObcGbsaForces1(gpuContext gpu);

extern void SetCalculateAmoebaObcGbsaForces2Sim(amoebaGpuContext amoebaGpu);
extern void GetCalculateAmoebaObcGbsaForces2Sim(amoebaGpuContext amoebaGpu);
extern void kCalculateAmoebaObcGbsaForces2(  amoebaGpuContext amoebaGpu );

extern void  cudaReduceN2ToN( float *N2Array, int N, float *NArray, int includeDiagonal, int offset );
extern float cudaGetSum( int numberOfElements, CUDAStream<float>* array );
extern float cudaGetNorm2( int numberOfElements, CUDAStream<float>* array );
extern int   checkForNansAndInfinities( int numberOfElements, CUDAStream<float>* array );
extern void cudaWriteFloat1AndFloat1ArraysToFile( int numberOfAtoms, const std::string& fname, std::vector<int>& fileId, int entriesPerAtom1, CUDAStream<float>* array1, 
                                                  int entriesPerAtom2, CUDAStream<float>* array2 );
extern void readFile( std::string fileName, StringVectorVector& fileContents );
 
extern void cudaLoadCudaFloatArray(  int numberOfParticles, int entriesPerParticle, CUDAStream<float>*  array, VectorOfDoubleVectors& outputVector, int* order, float conversion );
extern void cudaLoadCudaFloat2Array( int numberOfParticles, int entriesPerParticle, CUDAStream<float2>* array, VectorOfDoubleVectors& outputVector, int* order, float conversion );
extern void cudaLoadCudaFloat4Array( int numberOfParticles, int entriesPerParticle, CUDAStream<float4>* array, VectorOfDoubleVectors& outputVector, int* order, float conversion );
extern void cudaWriteVectorOfDoubleVectorsToFile( const std::string& fname, std::vector<int>& fileId, VectorOfDoubleVectors& outputVector );
extern void initializeCudaFloatArray( int numberOfParticles, int entriesPerParticle, CUDAStream<float>* array, float initValue );
extern void checkForNans( int numberOfParticles, int entriesPerParticle,
                          CUDAStream<float>* array, int* order, int iteration, std::string idString, FILE* log );
extern void checkForNansFloat4( int numberOfParticles, CUDAStream<float4>* array, int* order, int iteration, std::string idString, FILE* log );



extern void kClearFloat( amoebaGpuContext amoebaGpu, unsigned int entries, CUDAStream<float>* fieldToClear );
extern void kClearFloat4( amoebaGpuContext amoebaGpu, unsigned int entries, CUDAStream<float4>* fieldToClear );
extern void kClearFields_1( amoebaGpuContext amoebaGpu );
extern void kClearFields_3( amoebaGpuContext amoebaGpu, unsigned int numberToClear );
extern unsigned int getThreadsPerBlock( amoebaGpuContext amoebaGpu, unsigned int sharedMemoryPerThread, unsigned int sharedMemoryPerBlock );

//extern int isNanOrInfinity( double number );
extern void trackMutualInducedIterations( amoebaGpuContext amoebaGpu, int iteration);
extern void zeroCUDAStreamFloat4( CUDAStream<float4>* streamToCopy );
extern void reduceAndCopyCUDAStreamFloat4( CUDAStream<float4>* streamToCopy, CUDAStream<float>*  outputStream, float conversion );

// PME

extern void SetCalculateAmoebaPMESim( amoebaGpuContext amoebaGpu );
extern void kCalculateAmoebaPMEFixedMultipoles(amoebaGpuContext amoebaGpu);
extern void kCalculateAmoebaPMEInducedDipoleField(amoebaGpuContext amoebaGpu);
extern void kCalculateAmoebaPMEInducedDipoleForces(amoebaGpuContext amoebaGpu);


extern void SetCalculateAmoebaCudaUtilitiesSim( amoebaGpuContext amoebaGpu );
#endif //__AMOEBA_GPU_TYPES_H__

