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

#include "cudatypes.h"
#include "cudaCompact.h"
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

struct gpuTabulatedFunction {
    gpuTabulatedFunction() : coefficients(NULL) {
    }
    std::string name;
    double min, max;
    CUDAStream<float4>* coefficients;
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
    int device;
    bool useBlockingSync;
    gpuAtomType* gpAtomTable;
    int gAtomTypes;
    cudaGmxSimulation sim;
    unsigned int* pOutputBufferCounter;
    std::vector<std::vector<int> > exclusions;
    unsigned char* pAtomSymbol;
    std::vector<gpuMoleculeGroup> moleculeGroups;
    gpuTabulatedFunction tabulatedFunctions[MAX_TABULATED_FUNCTIONS];
    std::vector<int3> posCellOffsets;
    int iterations;
    float epsfac;
    float solventDielectric;
    float soluteDielectric;
    int grid;
    bool bCalculateCM;
    bool bRemoveCM;
    bool bRecalculateBornRadii;
    bool bOutputBufferPerWarp;
    bool bIncludeGBSA;
    bool bIncludeGBVI;
    bool tabulatedFunctionsChanged;
    unsigned long seed;
    SM_VERSION sm_version;
    compactionPlan compactPlan;
    cufftHandle fftplan;
    CUDAStream<float4>* psPosq4;
    CUDAStream<float4>* psPosqP4;
    CUDAStream<float4>* psOldPosq4;
    CUDAStream<float4>* psVelm4;
    CUDAStream<float4>* psForce4;
    CUDAStream<float>*  psEnergy;           // Energy output buffer
    CUDAStream<float4>* psxVector4;
    CUDAStream<float4>* psvVector4;
    CUDAStream<float2>* psSigEps2; 
    CUDAStream<float4>* psCustomParams;     // Atom parameters for custom nonbonded force
    CUDAStream<int4>* psCustomBondID;             // Atom indices for custom bonds
    CUDAStream<float4>* psCustomBondParams;       // Parameters for custom bonds
    CUDAStream<int4>* psCustomAngleID1;           // Atom indices for custom angles
    CUDAStream<int2>* psCustomAngleID2;           // Atom indices for custom angles
    CUDAStream<float4>* psCustomAngleParams;      // Parameters for custom angles
    CUDAStream<int>* psCustomExternalID;          // Atom indices for custom external force
    CUDAStream<float4>* psCustomExternalParams;   // Parameters for custom external force
    CUDAStream<float4>* psTabulatedFunctionParams; // The min, max, and spacing for each tabulated function
    CUDAStream<float2>* psEwaldCosSinSum;
    CUDAStream<float>* psTabulatedErfc;     // Tabulated values for erfc()
    CUDAStream<cufftComplex>* psPmeGrid;    // Grid points for particle mesh Ewald
    CUDAStream<float>* psPmeBsplineModuli[3];
    CUDAStream<float4>* psPmeBsplineTheta;
    CUDAStream<float4>* psPmeBsplineDtheta;
    CUDAStream<int>* psPmeAtomRange;           // The range of sorted atoms at each grid point
    CUDAStream<float2>* psPmeAtomGridIndex;    // The grid point each atom is at
    CUDAStream<float2>* psObcData;
    CUDAStream<float4>* psGBVIData;
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
    CUDAStream<unsigned int>* psExclusionIndex;
    CUDAStream<unsigned int>* psWorkUnit;
    CUDAStream<unsigned int>* psInteractingWorkUnit;
    CUDAStream<unsigned int>* psInteractionFlag;
    CUDAStream<size_t>* psInteractionCount;
    CUDAStream<float2>* psStepSize;         // The size of the previous and current time steps
    CUDAStream<float>* psLangevinParameters;// Parameters used for Langevin integration
    CUDAStream<float4>* psRandom4;          // Pointer to sets of 4 random numbers for MD integration
    CUDAStream<float2>* psRandom2;          // Pointer to sets of 2 random numbers for MD integration
    CUDAStream<uint4>* psRandomSeed;        // Pointer to each random seed
    CUDAStream<int>* psRandomPosition;      // Pointer to random number positions
    CUDAStream<float4>* psLinearMomentum;   // Pointer to total linear momentum per CTA
    CUDAStream<int>* psAtomIndex;           // The original index of each atom
    CUDAStream<float4>* psGridBoundingBox;  // The size of each grid cell
    CUDAStream<float4>* psGridCenter;       // The center and radius for each grid cell
    CUDAStream<int2>* psCcmaAtoms;          // The atoms connected by each CCMA constraint
    CUDAStream<float4>* psCcmaDistance;     // The displacement vector (x, y, z) and constraint distance (w) for each CCMA constraint
    CUDAStream<int>* psCcmaAtomConstraints; // The indices of constraints involving each atom
    CUDAStream<int>* psCcmaNumAtomConstraints; // The number of constraints involving each atom
    CUDAStream<float>* psCcmaDelta1;        // Workspace for CCMA
    CUDAStream<float>* psCcmaDelta2;        // Workspace for CCMA
    CUDAStream<short>* psSyncCounter;       // Used for global thread synchronization
    CUDAStream<unsigned int>* psRequiredIterations; // Used by SHAKE to communicate whether iteration has converged
    CUDAStream<float>* psCcmaReducedMass;   // The reduced mass for each SHAKE constraint
    CUDAStream<float>* psRigidClusterMatrix;// The inverse constraint matrix for each rigid cluster
    CUDAStream<unsigned int>* psRigidClusterConstraintIndex; // The index of each cluster in the stream containing cluster constraints.
    CUDAStream<unsigned int>* psRigidClusterMatrixIndex; // The index of each cluster in the stream containing cluster matrices.
    CUDAStream<unsigned int>* psConstraintMatrixColumn; // The column of each element in the constraint matrix.
    CUDAStream<float>* psConstraintMatrixValue; // The value of each element in the constraint matrix.
};

typedef struct _gpuContext *gpuContext;


// Function prototypes
extern "C"
bool gpuIsAvailable();

extern "C"
void gpuSetBondParameters(gpuContext gpu, const std::vector<int>& atom1, const std::vector<int>& atom2, const std::vector<float>& length, const std::vector<float>& k);

extern "C"
void gpuSetBondAngleParameters(gpuContext gpu, const std::vector<int>& atom1, const std::vector<int>& atom2, const std::vector<int>& atom3,
        const std::vector<float>& angle, const std::vector<float>& k);

extern "C"
void gpuSetDihedralParameters(gpuContext gpu, const std::vector<int>& atom1, const std::vector<int>& atom2, const std::vector<int>& atom3, const std::vector<int>& atom4,
        const std::vector<float>& k, const std::vector<float>& phase, const std::vector<int>& periodicity);

extern "C"
void gpuSetRbDihedralParameters(gpuContext gpu, const std::vector<int>& atom1, const std::vector<int>& atom2, const std::vector<int>& atom3, const std::vector<int>& atom4,
        const std::vector<float>& c0, const std::vector<float>& c1, const std::vector<float>& c2, const std::vector<float>& c3, const std::vector<float>& c4, const std::vector<float>& c5);

extern "C"
void gpuSetLJ14Parameters(gpuContext gpu, float epsfac, float fudge, const std::vector<int>& atom1, const std::vector<int>& atom2,
        const std::vector<float>& c6, const std::vector<float>& c12, const std::vector<float>& q1, const std::vector<float>& q2);

extern "C"
void gpuSetCoulombParameters(gpuContext gpu, float epsfac, const std::vector<int>& atom, const std::vector<float>& c6, const std::vector<float>& c12, const std::vector<float>& q,
        const std::vector<char>& symbol, const std::vector<std::vector<int> >& exclusions, CudaNonbondedMethod method);

extern "C"
void gpuSetNonbondedCutoff(gpuContext gpu, float cutoffDistance, float solventDielectric);

extern "C"
void gpuSetTabulatedFunction(gpuContext gpu, int index, const std::string& name, const std::vector<double>& values, double min, double max, bool interpolating);

extern "C"
void gpuSetCustomBondParameters(gpuContext gpu, const std::vector<int>& bondAtom1, const std::vector<int>& bondAtom2, const std::vector<std::vector<double> >& bondParams,
            const std::string& energyExp, const std::vector<std::string>& paramNames, const std::vector<std::string>& globalParamNames);

extern "C"
void gpuSetCustomAngleParameters(gpuContext gpu, const std::vector<int>& bondAtom1, const std::vector<int>& bondAtom2, const std::vector<int>& bondAtom3, const std::vector<std::vector<double> >& bondParams,
            const std::string& energyExp, const std::vector<std::string>& paramNames, const std::vector<std::string>& globalParamNames);

extern "C"
void gpuSetCustomExternalParameters(gpuContext gpu, const std::vector<int>& atomIndex, const std::vector<std::vector<double> >& atomParams,
            const std::string& energyExp, const std::vector<std::string>& paramNames, const std::vector<std::string>& globalParamNames);

extern "C"
void gpuSetCustomNonbondedParameters(gpuContext gpu, const std::vector<std::vector<double> >& parameters, const std::vector<std::vector<int> >& exclusions,
            CudaNonbondedMethod method, float cutoffDistance, const std::string& energyExp,
            const std::vector<std::string>& paramNames, const std::vector<std::string>& globalParamNames);

extern "C"
void gpuSetEwaldParameters(gpuContext gpu, float alpha, int kmaxx, int kmaxy, int kmaxz);

extern "C"
void gpuSetPMEParameters(gpuContext gpu, float alpha, int gridSizeX, int gridSizeY, int gridSizeZ);

extern "C"
void gpuSetPeriodicBoxSize(gpuContext gpu, float xsize, float ysize, float zsize);

extern "C"
void gpuSetObcParameters(gpuContext gpu, float innerDielectric, float solventDielectric, const std::vector<float>& radius, const std::vector<float>& scale, const std::vector<float>& charge);

extern "C" 
void gpuSetGBVIParameters(gpuContext gpu, float innerDielectric, float solventDielectric, const std::vector<int>& atom, const std::vector<float>& radius,
                          const std::vector<float>& gammas, const std::vector<float>& scaledRadii);

extern "C"
void gpuSetConstraintParameters(gpuContext gpu, const std::vector<int>& atom1, const std::vector<int>& atom2, const std::vector<float>& distance,
        const std::vector<float>& invMass1, const std::vector<float>& invMass2, float constraintTolerance);

extern "C"
int gpuAllocateInitialBuffers(gpuContext gpu);

extern "C"
void gpuSetPositions(gpuContext gpu, const std::vector<float>& x, const std::vector<float>& y, const std::vector<float>& z);

extern "C"
void gpuSetVelocities(gpuContext gpu, const std::vector<float>& x, const std::vector<float>& y, const std::vector<float>& z);

extern "C"
void gpuSetMass(gpuContext gpu, const std::vector<float>& mass);

extern "C"
void gpuInitializeRandoms(gpuContext gpu);

extern "C"
void* gpuInit(int numAtoms, unsigned int device = 0, bool useBlockingSync = false);

extern "C"
void gpuSetLangevinIntegrationParameters(gpuContext gpu, float tau, float deltaT, float temperature, float errorTol);

extern "C"
void gpuSetVerletIntegrationParameters(gpuContext gpu, float deltaT, float errorTol);

extern "C"
void gpuSetBrownianIntegrationParameters(gpuContext gpu, float tau, float deltaT, float temperature);

extern "C"
void gpuSetAndersenThermostatParameters(gpuContext gpu, float temperature, float collisionFrequency);

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
void gpuReorderAtoms(gpuContext gpu);

extern "C"
void setExclusions(gpuContext gpu, const std::vector<std::vector<int> >& exclusions);

#endif //__GPUTYPES_H__
