#ifndef __GPUTYPES_H__
#define __GPUTYPES_H__

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

#include "cudatypes.h"
#include "cudpp.h"
#include <vector>

struct gpuAtomType {
    std::string name;
    char symbol;
    float r;
};

struct gpuMoleculeGroup {
    std::vector<int> atoms;
    std::vector<int> instances;
};

enum SM_VERSION
{
    SM_10,
    SM_11,
    SM_12
};


/* Pointer to this structure will be given 
 * to gromacs functions*/
struct _gpuContext {
    
    
    //Cache this here so that it doesn't
    //have to be repeatedly passed around
    int natoms;
    gpuAtomType* gpAtomTable;
    int gAtomTypes;
    cudaGmxSimulation sim;
    unsigned int* pOutputBufferCounter;
    std::vector<std::vector<int> > exclusions;
    unsigned char* pAtomSymbol;
    std::vector<gpuMoleculeGroup> moleculeGroups;
    float iterations;
    float epsfac;
    float solventDielectric;
    float soluteDielectric;
    int grid;
    bool bCalculateCM;
    bool bRemoveCM;
    bool bRecalculateBornRadii;
    bool bOutputBufferPerWarp;
    bool bIncludeGBSA;
    unsigned long seed;
    SM_VERSION sm_version;
    CUDPPHandle cudpp;
    CUDAStream<float4>* psPosq4;
    CUDAStream<float4>* psPosqP4;
    CUDAStream<float4>* psOldPosq4;
    CUDAStream<float4>* psVelm4;
    CUDAStream<float4>* psForce4;
    CUDAStream<float4>* psxVector4;
    CUDAStream<float4>* psvVector4;
    CUDAStream<float2>* psSigEps2; 
    CUDAStream<float2>* psObcData; 
    CUDAStream<float>* psObcChain;
    CUDAStream<float>* psBornForce;
    CUDAStream<float>* psBornRadii;
    CUDAStream<float>* psBornSum;
    CUDAStream<int4>* psBondID;
    CUDAStream<float2>* psBondParameter;
    CUDAStream<int4>* psBondAngleID1;
    CUDAStream<int2>* psBondAngleID2;
    CUDAStream<float2>* psBondAngleParameter;
    CUDAStream<int4>* psDihedralID1;
    CUDAStream<int4>* psDihedralID2;
    CUDAStream<float4>* psDihedralParameter;
    CUDAStream<int4>* psRbDihedralID1;
    CUDAStream<int4>* psRbDihedralID2;
    CUDAStream<float4>* psRbDihedralParameter1;
    CUDAStream<float2>* psRbDihedralParameter2;
    CUDAStream<int4>* psLJ14ID;
    CUDAStream<float4>* psLJ14Parameter;
    CUDAStream<int>* psNonShakeID;
    CUDAStream<int4>* psShakeID;
    CUDAStream<float4>* psShakeParameter;
    CUDAStream<int4>* psSettleID;
    CUDAStream<float2>* psSettleParameter;
    CUDAStream<unsigned int>* psExclusion;
    CUDAStream<unsigned int>* psWorkUnit;
    CUDAStream<unsigned int>* psInteractingWorkUnit;
    CUDAStream<unsigned int>* psInteractionFlag;
    CUDAStream<size_t>* psInteractionCount;
    CUDAStream<float4>* psRandom4;          // Pointer to sets of 4 random numbers for MD integration
    CUDAStream<float2>* psRandom2;          // Pointer to sets of 2 random numbers for MD integration
    CUDAStream<uint4>* psRandomSeed;        // Pointer to each random seed
    CUDAStream<int>* psRandomPosition;      // Pointer to random number positions
    CUDAStream<float4>* psLinearMomentum;   // Pointer to total linear momentum per CTA
    CUDAStream<int>* psAtomIndex;           // The original index of each atom
    CUDAStream<float4>* psGridBoundingBox;  // The size of each grid cell
    CUDAStream<float4>* psGridCenter;       // The center and radius for each grid cell
};

typedef struct _gpuContext *gpuContext;


// Function prototypes
extern "C"
bool gpuIsAvailable();

extern "C"
int gpuReadBondParameters(gpuContext gpu, char* fname);

extern "C"
void gpuSetBondParameters(gpuContext gpu, const std::vector<int>& atom1, const std::vector<int>& atom2, const std::vector<float>& length, const std::vector<float>& k);

extern "C"
int gpuReadBondAngleParameters(gpuContext gpu, char* fname);

extern "C"
void gpuSetBondAngleParameters(gpuContext gpu, const std::vector<int>& atom1, const std::vector<int>& atom2, const std::vector<int>& atom3,
        const std::vector<float>& angle, const std::vector<float>& k);

extern "C"
int gpuReadDihedralParameters(gpuContext gpu, char* fname);

extern "C"
void gpuSetDihedralParameters(gpuContext gpu, const std::vector<int>& atom1, const std::vector<int>& atom2, const std::vector<int>& atom3, const std::vector<int>& atom4,
        const std::vector<float>& k, const std::vector<float>& phase, const std::vector<int>& periodicity);

extern "C"
int gpuReadRbDihedralParameters(gpuContext gpu, char* fname);

extern "C"
void gpuSetRbDihedralParameters(gpuContext gpu, const std::vector<int>& atom1, const std::vector<int>& atom2, const std::vector<int>& atom3, const std::vector<int>& atom4,
        const std::vector<float>& c0, const std::vector<float>& c1, const std::vector<float>& c2, const std::vector<float>& c3, const std::vector<float>& c4, const std::vector<float>& c5);

extern "C"
int gpuReadLJ14Parameters(gpuContext gpu, char* fname);

extern "C"
void gpuSetLJ14Parameters(gpuContext gpu, float epsfac, float fudge, const std::vector<int>& atom1, const std::vector<int>& atom2,
        const std::vector<float>& c6, const std::vector<float>& c12, const std::vector<float>& q1, const std::vector<float>& q2);

extern "C"
float gpuGetAtomicRadius(gpuContext gpu, std::string s);

extern "C"
unsigned char gpuGetAtomicSymbol(gpuContext gpu, std::string s);

extern "C"
int gpuReadAtomicParameters(gpuContext gpu, char* fname);

extern "C"
int gpuReadCoulombParameters(gpuContext gpu, char* fname);

extern "C"
void gpuSetCoulombParameters(gpuContext gpu, float epsfac, const std::vector<int>& atom, const std::vector<float>& c6, const std::vector<float>& c12, const std::vector<float>& q,
        const std::vector<char>& symbol, const std::vector<std::vector<int> >& exclusions, CudaNonbondedMethod method);

extern "C"
void gpuSetNonbondedCutoff(gpuContext gpu, float cutoffDistance, float solventDielectric);

extern "C"
void gpuSetPeriodicBoxSize(gpuContext gpu, float xsize, float ysize, float zsize);

extern "C"
void gpuSetObcParameters(gpuContext gpu, float innerDielectric, float solventDielectric, const std::vector<int>& atom, const std::vector<float>& radius, const std::vector<float>& scale);

extern "C"
int gpuReadShakeParameters(gpuContext gpu, char* fname);

extern "C"
void gpuSetShakeParameters(gpuContext gpu, const std::vector<int>& atom1, const std::vector<int>& atom2, const std::vector<float>& distance,
        const std::vector<float>& invMass1, const std::vector<float>& invMass2, float tolerance);

extern "C"
int gpuAllocateInitialBuffers(gpuContext gpu);

extern "C"
void gpuReadCoordinates(gpuContext gpu, char* fname);

extern "C"
void gpuSetPositions(gpuContext gpu, const std::vector<float>& x, const std::vector<float>& y, const std::vector<float>& z);

extern "C"
void gpuSetVelocities(gpuContext gpu, const std::vector<float>& x, const std::vector<float>& y, const std::vector<float>& z);

extern "C"
void gpuSetMass(gpuContext gpu, const std::vector<float>& mass);

extern "C"
void gpuInitializeRandoms(gpuContext gpu);

extern "C"
void* gpuInitFromFile(char* fname);

extern "C"
void* gpuInit(int numAtoms);

extern "C"
void gpuSetIntegrationParameters(gpuContext gpu, float tau, float deltaT, float temperature);

extern "C"
void gpuSetVerletIntegrationParameters(gpuContext gpu, float deltaT);

extern "C"
void gpuSetBrownianIntegrationParameters(gpuContext gpu, float tau, float deltaT, float temperature);

extern "C"
void gpuSetAndersenThermostatParameters(gpuContext gpu, float temperature, float collisionProbability);

extern "C"
void gpuShutDown(gpuContext gpu);

extern "C"
int gpuBuildOutputBuffers(gpuContext gpu);

extern "C"
int gpuBuildThreadBlockWorkList(gpuContext gpu);

extern "C"
void gpuBuildExclusionList(gpuContext gpu);

extern "C"
int gpuSetConstants(gpuContext gpu);

extern "C"
void gpuDumpCoordinates(gpuContext gpu);

extern "C"
void gpuDumpPrimeCoordinates(gpuContext gpu);

extern "C"
void gpuDumpForces(gpuContext gpu);

extern "C"
void gpuDumpAtomData(gpuContext gpu);

extern "C"
bool gpuCheckData(gpuContext gpu);

extern "C"
void gpuSetup(void* pVoid);

extern "C"
void kCPUCalculate14(gpuContext gpu);

extern "C"
void kCPUCalculateLocalForces(gpuContext gpu);

extern "C"
void WriteArrayToFile1( gpuContext gpu, char* fname, int step, CUDAStream<float>*  psPos, int numPrint );

extern "C"
void WriteArrayToFile2( gpuContext gpu, char* fname, int step, CUDAStream<float2>* psPos, int numPrint );

extern "C"
void WriteArrayToFile3( gpuContext gpu, char* fname, int step, CUDAStream<float3>* psPos, int numPrint );

extern "C"
void WriteArrayToFile4( gpuContext gpu, char* fname, int step, CUDAStream<float4>* psPos, int numPrint );

extern "C"
void gpuDumpObcInfo(gpuContext gpu);

extern "C"
void gpuDumpObcLoop1(gpuContext gpu); 

extern "C"
void gpuReorderAtoms(gpuContext gpu);

#endif //__GPUTYPES_H__
