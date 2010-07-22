#ifndef AMOEBA_CUDATYPES_H
#define AMOEBA_CUDATYPES_H

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

#include <kernels/cudatypes.h>

#include <stdarg.h>
#include <limits>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cufft.h>
#include <builtin_types.h>
#include <vector_functions.h>

#if 0
#define RTERROR(status, s) \
    if (status != cudaSuccess) { \
        printf("%s %s\n", s, cudaGetErrorString(status)); \
        exit(-1); \
    }

#define LAUNCHERROR(s) \
    { \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            exit(-1); \
        } \
    }
#endif

#if 0
// Pure virtual class to define an interface for objects resident both on GPU and CPU
struct SoADeviceObject {
    virtual void Allocate() = 0;
    virtual void Deallocate() = 0;
    virtual void Upload() = 0;
    virtual void Download() = 0;
};

template <typename T>
struct CUDAStream : public SoADeviceObject
{
    unsigned int    _length;
    unsigned int    _subStreams;
    unsigned int    _stride;
    T**             _pSysStream;
    T**             _pDevStream;
    T*              _pSysData;
    T*              _pDevData;
    std::string     _name;
    CUDAStream(int length, int subStreams = 1, std::string name="");
    CUDAStream(unsigned int length, unsigned int subStreams = 1, std::string name="");
    CUDAStream(unsigned int length, int subStreams = 1, std::string name="");
    CUDAStream(int length, unsigned int subStreams = 1, std::string name="");
    virtual ~CUDAStream();
    void Allocate();
    void Deallocate();
    void Upload();
    void Download();
    void Collapse(unsigned int newstreams = 1, unsigned int interleave = 1);
    T& operator[](int index);
};

float CompareStreams(CUDAStream<float>& s1, CUDAStream<float>& s2, float tolerance, unsigned int maxindex = 0);

template <typename T>
CUDAStream<T>::CUDAStream(int length, unsigned int subStreams, std::string name) : _length(length), _subStreams(subStreams), _stride((length + 0xf) & 0xfffffff0), _name(name)
{
    Allocate();   
}

template <typename T>
CUDAStream<T>::CUDAStream(unsigned int length, int subStreams, std::string name) : _length(length), _subStreams(subStreams), _stride((length + 0xf) & 0xfffffff0), _name(name)
{
    Allocate();   
}

template <typename T>
CUDAStream<T>::CUDAStream(unsigned int length, unsigned int subStreams, std::string name) : _length(length), _subStreams(subStreams), _stride((length + 0xf) & 0xfffffff0), _name(name)
{
    Allocate();   
}

template <typename T>
CUDAStream<T>::CUDAStream(int length, int subStreams, std::string name) : _length(length), _subStreams(subStreams), _stride((length + 0xf) & 0xfffffff0), _name(name)
{
    Allocate();   
}

template <typename T>
CUDAStream<T>::~CUDAStream()
{
    Deallocate();
}

template <typename T>
void CUDAStream<T>::Allocate()
{
    cudaError_t status;
    _pSysStream =   new T*[_subStreams];
    _pDevStream =   new T*[_subStreams];
    _pSysData =     new T[_subStreams * _stride];

    status = cudaMalloc((void **) &_pDevData, _stride * _subStreams * sizeof(T));
    RTERROR(status, (_name+": cudaMalloc in CUDAStream::Allocate failed").c_str());

    for (unsigned int i = 0; i < _subStreams; i++)
    {
        _pSysStream[i] = _pSysData + i * _stride;
        _pDevStream[i] = _pDevData + i * _stride;
    }
}

template <typename T>
void CUDAStream<T>::Deallocate()
{
    cudaError_t status;
    delete[] _pSysStream;
    _pSysStream = NULL;
    delete[] _pDevStream;
    _pDevStream = NULL;
    delete[] _pSysData;
    _pSysData = NULL;
    status = cudaFree(_pDevData);
    RTERROR(status, (_name+": cudaFree in CUDAStream::Deallocate failed").c_str());
}

template <typename T>
void CUDAStream<T>::Upload()
{
    cudaError_t status;
    status = cudaMemcpy(_pDevData, _pSysData, _stride * _subStreams * sizeof(T), cudaMemcpyHostToDevice);
    RTERROR(status, (_name+": cudaMemcpy in CUDAStream::Upload failed").c_str());
}

template <typename T>
void CUDAStream<T>::Download()
{
    cudaError_t status;
    status = cudaMemcpy(_pSysData, _pDevData, _stride * _subStreams * sizeof(T), cudaMemcpyDeviceToHost);
    RTERROR(status, (_name+": cudaMemcpy in CUDAStream::Download failed").c_str());
}

template <typename T>
void CUDAStream<T>::Collapse(unsigned int newstreams, unsigned int interleave)
{
    T* pTemp = new T[_subStreams * _stride];
    unsigned int stream = 0;
    unsigned int pos = 0;
    unsigned int newstride = _stride * _subStreams / newstreams;
    unsigned int newlength = _length * _subStreams / newstreams;

    // Copy data into new format
    for (unsigned int i = 0; i < _length; i++)
    {
        for (unsigned int j = 0; j < _subStreams; j++)
        {
            pTemp[stream * newstride + pos] = _pSysStream[j][i];
            stream++;
            if (stream == newstreams)
            {
                stream = 0;
                pos++;
            }
        }
    }

    // Remap stream pointers;
    for (unsigned int i = 0; i < newstreams; i++)
    {
        _pSysStream[i] = _pSysData + i * newstride;
        _pDevStream[i] = _pDevData + i * newstride;
    }

    // Copy data back intro original stream
    for (unsigned int i = 0; i < newlength; i++)
        for (unsigned int j = 0; j < newstreams; j++)
            _pSysStream[j][i] = pTemp[j * newstride + i];
    
    _stride = newstride;
    _length = newlength;
    _subStreams = newstreams;
    delete[] pTemp;
}

template <typename T>
T& CUDAStream<T>::operator[](int index)
{
    return _pSysData[index];
}

static const unsigned int GRID = 32;
static const unsigned int GRIDBITS = 5;
static const int G8X_NONBOND_THREADS_PER_BLOCK          = 256;
static const int GT2XX_NONBOND_THREADS_PER_BLOCK        = 320;
static const int G8X_BORNFORCE2_THREADS_PER_BLOCK       = 256;
static const int GT2XX_BORNFORCE2_THREADS_PER_BLOCK     = 320;
static const int G8X_SHAKE_THREADS_PER_BLOCK            = 128;
static const int GT2XX_SHAKE_THREADS_PER_BLOCK          = 256;
static const int G8X_UPDATE_THREADS_PER_BLOCK           = 192;
static const int GT2XX_UPDATE_THREADS_PER_BLOCK         = 384;
static const int G8X_LOCALFORCES_THREADS_PER_BLOCK      = 192;
static const int GT2XX_LOCALFORCES_THREADS_PER_BLOCK    = 384;
static const int G8X_THREADS_PER_BLOCK                  = 256;
static const int GT2XX_THREADS_PER_BLOCK                = 256;
static const int G8X_RANDOM_THREADS_PER_BLOCK           = 256;
static const int GT2XX_RANDOM_THREADS_PER_BLOCK         = 384;
static const int G8X_NONBOND_WORKUNITS_PER_SM           = 220;
static const int GT2XX_NONBOND_WORKUNITS_PER_SM         = 256;
static const unsigned int MAX_STACK_SIZE = 8;
static const unsigned int MAX_TABULATED_FUNCTIONS = 4;

static const float PI = 3.14159265358979323846f;

static const int PME_ORDER = 4;

enum CudaNonbondedMethod
{
    NO_CUTOFF,
    CUTOFF,
    PERIODIC,
    EWALD,
    PARTICLE_MESH_EWALD
};

enum ExpressionOp {
    CONSTANT = 0, VARIABLE0, VARIABLE1, VARIABLE2, VARIABLE3, VARIABLE4, VARIABLE5, VARIABLE6, VARIABLE7, VARIABLE8, GLOBAL, CUSTOM, CUSTOM_DERIV, ADD, SUBTRACT, MULTIPLY, DIVIDE,
        POWER, NEGATE, SQRT, EXP, LOG, SIN, COS, SEC, CSC, TAN, COT, ASIN, ACOS, ATAN, SQUARE, CUBE, RECIPROCAL, ADD_CONSTANT, MULTIPLY_CONSTANT, POWER_CONSTANT
};

template<int SIZE>
struct Expression {
    int op[SIZE];
    float arg[SIZE];
    int length, stackSize;
};
#endif

struct cudaAmoebaGmxSimulation {
    // Constants

    unsigned int    amoebaBonds;                    // Number of bonds
    int4*           pAmoebaBondID;                  // Bond atom and output buffer IDs
    float4*         pAmoebaBondParameter;           // Bond parameters
    unsigned int    amoebaBond_offset;              // Offset to end of bonds

    unsigned int    amoebaAngles;                   // Number of bond angles
    int4*           pAmoebaAngleID1;                // Bond angle atom and first output buffer IDs
    int2*           pAmoebaAngleID2;                // Bond angle output buffer IDs
    float2*         pAmoebaAngleParameter;          // Bond angle parameters
    unsigned int    amoebaAngle_offset;             // Offset to end of bond angles

    float amoebaAngleCubicK;                        // cubic factor
    float amoebaAngleQuarticK;                      // quartic factor
    float amoebaAnglePenticK;                       // pentic factor
    float amoebaAngleSexticK;                       // sextic factor

    unsigned int    amoebaInPlaneAngles;            // Number of in-plane angles
    int4*           pAmoebaInPlaneAngleID1;         // Bond angle atom and first output buffer IDs
    int4*           pAmoebaInPlaneAngleID2;         // Bond angle output buffer IDs
    float2*         pAmoebaInPlaneAngleParameter;   // Bond angle parameters
    unsigned int    amoebaInPlaneAngle_offset;      // Offset to end of bond angles

    float amoebaInPlaneAngleCubicK;                 // cubic factor
    float amoebaInPlaneAngleQuarticK;               // quartic factor
    float amoebaInPlaneAnglePenticK;                // pentic factor
    float amoebaInPlaneAngleSexticK;                // sextic factor

    unsigned int    amoebaTorsions;                 // Number of torsions
    int4*           pAmoebaTorsionID1;              // Torsion atom and first output buffer IDs
    int4*           pAmoebaTorsionID2;              // Torsion output buffer IDs
    float4*         pAmoebaTorsionParameter1;       // Torsion parameters
    float2*         pAmoebaTorsionParameter2;       // Torsion parameters
    unsigned int    amoebaTorsion_offset;           // Offset to end of torsions

    unsigned int    amoebaPiTorsions;               // Number of torsions
    int4*           pAmoebaPiTorsionID1;            // PiTorsion atom and first output buffer IDs
    int4*           pAmoebaPiTorsionID2;            // PiTorsion output buffer IDs
    int4*           pAmoebaPiTorsionID3;            // PiTorsion output buffer IDs
    float*          pAmoebaPiTorsionParameter;      // PiTorsion parameters
    unsigned int    amoebaPiTorsion_offset;         // Offset to end of torsions

    unsigned int    amoebaStretchBends;             // Number of stretch bends
    int4*           pAmoebaStretchBendID1;          // stretch bend atoms and first output buffer IDs
    int2*           pAmoebaStretchBendID2;          // stretch bend output buffer IDs
    float4*         pAmoebaStretchBendParameter;    // stretch bend parameters
    unsigned int    amoebaStretchBend_offset;       // Offset to end of stretch bends

    unsigned int    amoebaOutOfPlaneBends;          // Number of stretch bends
    int4*           pAmoebaOutOfPlaneBendID1;       // stretch bend atoms and first output buffer IDs
    int4*           pAmoebaOutOfPlaneBendID2;       // stretch bend output buffer IDs
    float*          pAmoebaOutOfPlaneBendParameter; // stretch bend parameters
    unsigned int    amoebaOutOfPlaneBend_offset;    // Offset to end of stretch bends
    float amoebaOutOfPlaneBendCubicK;               // cubic factor
    float amoebaOutOfPlaneBendQuarticK;             // quartic factor
    float amoebaOutOfPlaneBendPenticK;              // pentic factor
    float amoebaOutOfPlaneBendSexticK;              // sextic factor

    unsigned int    amoebaTorsionTorsions;          // Number of torsion torsions
    int4*           pAmoebaTorsionTorsionID1;       // torsion torsion atoms and first output buffer IDs
    int4*           pAmoebaTorsionTorsionID2;       // torsion torsion output buffer IDs
    int4*           pAmoebaTorsionTorsionID3;       // torsion torsion parameters
    unsigned int    amoebaTorsionTorsion_offset;    // Offset to end of torsion torsions

                                                    // grids
    int   amoebaTorTorGridOffset[4];                // grid offset
    int   amoebaTorTorGridNy[4];                    // 25
    float amoebaTorTorGridBegin[4];                 // -180.0
    float amoebaTorTorGridDelta[4];                 // 15.0
    float4*          pAmoebaTorsionTorsionGrids;    // torsion torsion grids

    unsigned int numberOfAtoms;                     // number of atoms
    unsigned int paddedNumberOfAtoms;               // padded number of atoms
    float scalingDistanceCutoff;                    // scaling cutoff
    float2*         pDampingFactorAndThole;         // Thole & damping factors

    unsigned int amoebaVdwNonReductions;
    int* pAmoebaVdwNonReductionID;

    unsigned int amoebaVdwReductions;
    int4* pAmoebaVdwReductionID;
    float* pAmoebaVdwReduction;
    int* pVdwExclusionIndicesIndex;
    int* pVdwExclusionIndices;

    // WCA constants
    float epso;
    float epsh;
    float rmino;
    float rminh;
    float awater;
    float shctd;
    float dispoff;
 
    float totalMaxWcaDispersionEnergy;
                                                    // scaling indices
    int*            pScaleIndicesIndex;
    int*            pD_ScaleIndices;
    int2*           pP_ScaleIndices;
    int2*           pM_ScaleIndices;

    float electric;   // 3.320637090E+02f;
    float gkc;        // 2.455f;

    float dielec;    // 1.0f;
    float dwater;    // 78.3f;
    float fc;        // electric * 1.0f * (1.0f-dwater)/(0.0f+1.0f*dwater);
    float fd;        // electric * 2.0f * (1.0f-dwater)/(1.0f+2.0f*dwater);
    float fq;        // electric * 3.0f * (1.0f-dwater)/(2.0f+3.0f*dwater);

    // SASA probe radius
    float probeRadius; 
    int maxarc; 

//    texture<float4,2,cudaReadModeElementType> texTorTorGrid;

#if 0
    unsigned int    dihedrals;                      // Number of dihedrals
    int4*           pDihedralID1;                   // Dihedral IDs
    int4*           pDihedralID2;                   // Dihedral output buffer IDs
    float4*         pDihedralParameter;             // Dihedral parameters
    unsigned int    rb_dihedrals;                   // Number of Ryckaert Bellemans dihedrals
    int4*           pRbDihedralID1;                 // Ryckaert Bellemans Dihedral IDs
    int4*           pRbDihedralID2;                 // Ryckaert Bellemans Dihedral output buffer IDs
    float4*         pRbDihedralParameter1;          // Ryckaert Bellemans Dihedral parameters
    float2*         pRbDihedralParameter2;          // Ryckaert Bellemans Dihedral parameters
    unsigned int    LJ14s;                          // Number of Lennard Jones 1-4 interactions
    int4*           pLJ14ID;                        // Lennard Jones 1-4 atom and output buffer IDs
    float4*         pLJ14Parameter;                 // Lennard Jones 1-4 parameters
    float           inverseTotalMass;               // Used in linear momentum removal
    unsigned int    ShakeConstraints;               // Total number of Shake constraints
    unsigned int    settleConstraints;              // Total number of Settle constraints
    unsigned int    ccmaConstraints;                // Total number of CCMA constraints.
    unsigned int    rigidClusters;                  // Total number of rigid clusters
    unsigned int    maxRigidClusterSize;            // The size of the largest rigid cluster
    unsigned int    clusterShakeBlockSize;          // The number of threads to process each rigid cluster
    unsigned int    NonShakeConstraints;            // Total number of NonShake atoms
    unsigned int    maxShakeIterations;             // Maximum shake iterations
    unsigned int    degreesOfFreedom;               // Number of degrees of freedom in system
    float           shakeTolerance;                 // Shake tolerance
    float           InvMassJ;                       // Shake inverse mass for hydrogens
    int*            pNonShakeID;                    // Not Shaking atoms
    int4*           pShakeID;                       // Shake atoms and phase
    float4*         pShakeParameter;                // Shake parameters
    int4*           pSettleID;                      // Settle atoms
    float2*         pSettleParameter;               // Settle parameters
    unsigned int*   pExclusion;                     // Nonbond exclusion data
    unsigned int*   pExclusionIndex;                // Index of exclusion data for each work unit
    unsigned int    dihedral_offset;                // Offset to end of dihedrals
    unsigned int    rb_dihedral_offset;             // Offset to end of Ryckaert Bellemans dihedrals
    unsigned int    LJ14_offset;                    // Offset to end of Lennard Jones 1-4 parameters
    int*            pAtomIndex;                     // The original index of each atom
    float4*         pGridBoundingBox;               // The size of each grid cell
    float4*         pGridCenter;                    // The center of each grid cell
    int2*           pCcmaAtoms;                     // The atoms connected by each CCMA constraint
    float4*         pCcmaDistance;                  // The displacement vector (x, y, z) and constraint distance (w) for each CCMA constraint
    float*          pCcmaDelta1;                    // Workspace for CCMA
    float*          pCcmaDelta2;                    // Workspace for CCMA
    int*            pCcmaAtomConstraints;           // The indices of constraints involving each atom
    int*            pCcmaNumAtomConstraints;        // The number of constraints involving each atom
    short*          pSyncCounter;                   // Used for global thread synchronization
    unsigned int*   pRequiredIterations;            // Used by CCMA to communicate whether iteration has converged
    float*          pCcmaReducedMass;               // The reduced mass for each CCMA constraint
    unsigned int*   pConstraintMatrixColumn;        // The column of each element in the constraint matrix.
    float*          pConstraintMatrixValue;         // The value of each element in the constraint matrix.

    // Mutable stuff
    float4*         pPosq;                          // Pointer to atom positions and charges
    float4*         pPosqP;                         // Pointer to mid-integration atom positions
    float4*         pOldPosq;                       // Pointer to old atom positions
    float4*         pVelm4;                         // Pointer to atom velocity and inverse mass
    float4*         pvVector4;                      // Pointer to atom v Vector
    float4*         pxVector4;                      // Pointer to atom x Vector
    float*          pBornForce;                     // Pointer to Born force data
    float*	    pBornSum;                       // Pointer to Born Radii calculation output buffers
    float*	    pBornRadii;                     // Pointer to Born Radii
    float*          pObcChain;                      // Pointer to OBC chain data
    float4*         pLinearMomentum;                // Pointer to linear momentum
    
    // Random numbers
    float4*         pRandom4a;                      // Pointer to first set of 4 random numbers
    float4*         pRandom4b;                      // Pointer to second set of 4 random numbers
    float2*         pRandom2a;                      // Pointer to first set of 2 random numbers
    float2*         pRandom2b;                      // Pointer to second set of 2 random numbers
    uint4*          pRandomSeed;                    // Pointer to random seeds
    int*            pRandomPosition;                // Pointer to random number positions
    unsigned int    randoms;                        // Number of randoms
    unsigned int    totalRandoms;                   // Number of randoms plus overflow.
    unsigned int    totalRandomsTimesTwo;           // Used for generating randoms
    unsigned int    randomIterations;               // Number of iterations before regenerating randoms
    unsigned int    randomFrames;                   // Number of frames of random numbers
#endif

};

#endif
