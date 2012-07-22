#ifndef __AMOEBA_GPUTYPES_H__
#define __AMOEBA_GPUTYPES_H__

/* -------------------------------------------------------------------------- *
 *                          OpenMMAmoeba                                      *
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

#include "kernels/gputypes.h"
#include "OpenMM.h"
#include "openmm/Vec3.h"
#include "amoebaCudaTypes.h"

#include <map>
typedef std::map<int,float> MapIntFloat;
typedef MapIntFloat::const_iterator MapIntFloatCI;

struct _amoebaGpuContext {
    
    _gpuContext* gpuContext;
    cudaAmoebaGmxSimulation amoebaSim;
 
    FILE* log;

    CUDAStream<int4>*   psAmoebaBondID;
    CUDAStream<float2>* psAmoebaBondParameter;

    CUDAStream<int4>*   psAmoebaUreyBradleyID;
    CUDAStream<float2>* psAmoebaUreyBradleyParameter;

    CUDAStream<int4>*   psAmoebaAngleID1;
    CUDAStream<int2>*   psAmoebaAngleID2;
    CUDAStream<float2>* psAmoebaAngleParameter;

    CUDAStream<int4>*   psAmoebaInPlaneAngleID1;
    CUDAStream<int4>*   psAmoebaInPlaneAngleID2;
    CUDAStream<float2>* psAmoebaInPlaneAngleParameter;

    CUDAStream<int4>*   psAmoebaTorsionID1;
    CUDAStream<int4>*   psAmoebaTorsionID2;
    CUDAStream<float4>* psAmoebaTorsionParameter1;
    CUDAStream<float2>* psAmoebaTorsionParameter2;

    CUDAStream<int4>*   psAmoebaPiTorsionID1;
    CUDAStream<int4>*   psAmoebaPiTorsionID2;
    CUDAStream<int4>*   psAmoebaPiTorsionID3;
    CUDAStream<float>*  psAmoebaPiTorsionParameter;

    CUDAStream<int4>*   psAmoebaStretchBendID1;
    CUDAStream<int2>*   psAmoebaStretchBendID2;
    CUDAStream<float4>* psAmoebaStretchBendParameter;

    CUDAStream<int4>*   psAmoebaOutOfPlaneBendID1;
    CUDAStream<int4>*   psAmoebaOutOfPlaneBendID2;
    CUDAStream<float>*  psAmoebaOutOfPlaneBendParameter;

    CUDAStream<int4>*   psAmoebaTorsionTorsionID1;
    CUDAStream<int4>*   psAmoebaTorsionTorsionID2;
    CUDAStream<int4>*   psAmoebaTorsionTorsionID3;
    CUDAStream<float4>* psAmoebaTorsionTorsionGrids;

    unsigned int workUnits; 

    // workspace arrays

    CUDAStream<float>*  psWorkArray_3_1; 
    CUDAStream<float>*  psWorkArray_3_2; 
    CUDAStream<float>*  psWorkArray_3_3; 
    CUDAStream<float>*  psWorkArray_3_4; 

    CUDAStream<float>*  psWorkArray_1_1; 
    CUDAStream<float>*  psWorkArray_1_2; 

    CUDAStream<int>*  psScalingIndicesIndex; 
    CUDAStream<int>*  ps_D_ScaleIndices; 
    CUDAStream<int2>* ps_P_ScaleIndices; 
    CUDAStream<int2>* ps_M_ScaleIndices; 

    int maxCovalentDegreeSz;
    float solventDielectric;

    // multipole parameters

    CUDAStream<int4>* psMultipoleParticlesIdsAndAxisType;

    // buffer indices used for mapping torques onto forces 

    int torqueMapForce4Delete;
    CUDAStream<int4>*    psMultipoleParticlesTorqueBufferIndices;
    CUDAStream<float4>*  psTorqueMapForce4; 

    CUDAStream<float>* psMolecularDipole;
    CUDAStream<float>* psMolecularQuadrupole;

    
    CUDAStream<unsigned int>* psPotentialWorkUnit;
    CUDAStream<float4>* psPotentialGrid;
    CUDAStream<float>*  psPotential;

    // molecular frame multipoles

    CUDAStream<float>* psLabFrameDipole;
    CUDAStream<float>* psLabFrameQuadrupole;

    // scaling-related parameters

    CUDAStream<float2>*  psDampingFactorAndThole;

    // used to setup scaling constants

    std::vector<int>    covalentDegree;
    std::vector<int>    polarizationDegree;

    // fixed-E field

    CUDAStream<float>*  psE_Field;
    CUDAStream<float>*  psE_FieldPolar;

    int multipoleNonbondedMethod;
    double cutoffDistance;

    // mutual induced field

    int mutualInducedIterativeMethod;
    int mutualInducedMaxIterations;
    int mutualInducedConverged;
    int mutualInducedDone;

    int epsilonThreadsPerBlock;
    float mutualInducedTargetEpsilon;
    float mutualInducedCurrentEpsilon;
    CUDAStream<float>*  psInducedDipole; 
    CUDAStream<float>*  psInducedDipolePolar; 
    CUDAStream<float>*  psPolarizability; 
    CUDAStream<float>*  psCurrentEpsilon; 

    // SOR arrays for mutual induced field

    unsigned int numberOfSorWorkVectors;
    CUDAStream<float>*  psWorkVector[4]; 

    // electrostatic

    CUDAStream<float>*  psTorque; 

    // Kirkwood fields

    CUDAStream<float>*  psGk_Field;
    CUDAStream<float>*  psInducedDipoleS; 
    CUDAStream<float>*  psInducedDipolePolarS; 
    CUDAStream<float>*  psBorn; 
    CUDAStream<float>*  psBornPolar; 

    int includeObcCavityTerm;

    // Vdw fields

    CUDAStream<float2>*  psVdwSigmaEpsilon;

    CUDAStream<int>*     psAmoebaVdwNonReductionID; 
    CUDAStream<int4>*    psAmoebaVdwReductionID; 
    CUDAStream<float>*   psAmoebaVdwReduction; 
    CUDAStream<float4>*  psAmoebaVdwCoordinates; 

    CUDAStream<unsigned int>*   psVdwWorkUnit; 
    CUDAStream<int>* psVdwExclusionIndicesIndex;
    CUDAStream<int>* psVdwExclusionIndices;

    int vdwSigmaCombiningRule;
    int vdwEpsilonCombiningRule;
    std::vector< std::vector<int> > vdwExclusions;

    // Wca dispersion fields

    CUDAStream<float2>*  psWcaDispersionRadiusEpsilon;

    // PME fields

    CUDAStream<float4>* psThetai1;
    CUDAStream<float4>* psThetai2;
    CUDAStream<float4>* psThetai3;
    CUDAStream<int4>* psIgrid;
    CUDAStream<float>* psPhi;
    CUDAStream<float>* psPhid;
    CUDAStream<float>* psPhip;
    CUDAStream<float>* psPhidp;
};

typedef struct _amoebaGpuContext *amoebaGpuContext;

// Function prototypes

extern "C"
amoebaGpuContext amoebaGpuInit( _gpuContext* gpu );

extern "C"
void gpuPrintCudaAmoebaGmxSimulation(amoebaGpuContext gpu,  FILE* log );

extern "C"
void amoebaGpuShutDown(amoebaGpuContext gpu);

extern "C"
void amoebaGpuBuildOutputBuffers( amoebaGpuContext gpu, int hasKirkwood );

extern "C"
int amoebaGpuBuildThreadBlockWorkList( amoebaGpuContext gpu );

extern "C"
void amoebaGpuBuildScalingList( amoebaGpuContext gpu );

extern "C"
void gpuSetAmoebaBondParameters(amoebaGpuContext gpu, const std::vector<int>& atom1, const std::vector<int>& atom2, 
                                const std::vector<float>& length, const std::vector<float>& k, float cubic, float quartic);

extern "C"
void gpuSetAmoebaUreyBradleyParameters(amoebaGpuContext gpu, const std::vector<int>& atom1, const std::vector<int>& atom2, 
                                       const std::vector<float>& length, const std::vector<float>& k, float cubic, float quartic);

extern "C"
void gpuSetAmoebaAngleParameters(amoebaGpuContext gpu, const std::vector<int>& atom1, const std::vector<int>& atom2, const std::vector<int>& atom3,
                                 const std::vector<float>& angle, const std::vector<float>& k, float cubicK,
                                 float quarticK, float penticK, float sexticK);

extern "C"
void gpuSetAmoebaInPlaneAngleParameters(amoebaGpuContext gpu, const std::vector<int>& atom1, const std::vector<int>& atom2,
                                                              const std::vector<int>& atom3, const std::vector<int>& atom4,
                                        const std::vector<float>& angle, const std::vector<float>& k, float cubicK,
                                        float quarticK, float penticK, float sexticK);

extern "C"
void gpuSetAmoebaTorsionParameters(amoebaGpuContext gpu, const std::vector<int>& atom1, const std::vector<int>& atom2,
                                                         const std::vector<int>& atom3, const std::vector<int>& atom4,
                                   const std::vector< std::vector<float> >& torsion1,
                                   const std::vector< std::vector<float> >& torsion2,
                                   const std::vector< std::vector<float> >& torsion3 );

extern "C"
void gpuSetAmoebaPiTorsionParameters(amoebaGpuContext gpu, const std::vector<int>& atom1, const std::vector<int>& atom2,
                                                           const std::vector<int>& atom3, const std::vector<int>& atom4,
                                                           const std::vector<int>& atom5, const std::vector<int>& atom6,
                                                           const std::vector<float>& torsion1 );

extern "C"
void gpuSetAmoebaStretchBendParameters(amoebaGpuContext gpu, const std::vector<int>& atom1, const std::vector<int>& atom2, const std::vector<int>& atom3,
                                       const std::vector<float>& lengthAB,
                                       const std::vector<float>& lengthCB,
                                       const std::vector<float>& angle,
                                       const std::vector<float>& k );

extern "C"
void gpuSetAmoebaOutOfPlaneBendParameters(amoebaGpuContext gpu, const std::vector<int>& atom1, const std::vector<int>& atom2, const std::vector<int>& atom3,
                                          const std::vector<int>& atom4, const std::vector<float>& k,
                                          float cubicK, float quarticK, float penticK, float sexticK );

extern "C"
void gpuSetAmoebaTorsionTorsionParameters(amoebaGpuContext gpu, const std::vector<int>& atom1, const std::vector<int>& atom2, const std::vector<int>& atom3,
                                          const std::vector<int>& atom4, const std::vector<int>& atom5, const std::vector<int>& chiralAtomIndex, const std::vector<int>& gridIndex );

extern "C"
void gpuSetAmoebaTorsionTorsionGrids(amoebaGpuContext gpu, const std::vector< std::vector< std::vector< std::vector<float> > > >& floatGrids );

extern "C"  
void gpuSetAmoebaMultipoleParameters(amoebaGpuContext amoebaGpu, const std::vector<float>& charges, const std::vector<float>& dipoles, const std::vector<float>& quadrupoles,
                                     const std::vector<int>& axisType, const std::vector<int>& multipoleAtomZ, const std::vector<int>& multipoleAtomX,  const std::vector<int>& multipoleAtomY,
                                     const std::vector<float>& tholes, float scalingDistanceCutoff,const std::vector<float>& dampingFactors, const std::vector<float>& polarity,
                                     const std::vector< std::vector< std::vector<int> > >& multipoleAtomCovalentInfo, const std::vector<int>& covalentDegree,
                                     const std::vector<int>& minCovalentIndices,  const std::vector<int>& minCovalentPolarizationIndices, int maxCovalentRange,
                                     int mutualInducedIterationMethod, int mutualInducedMaxIterations, float mutualInducedTargetEpsilon,
                                     int nonbondedMethod, int polarizationType, float cutoffDistance,  float alphaEwald );


extern "C"
void gpuSetupElectrostaticPotentialCalculation( amoebaGpuContext amoebaGpu, const std::vector< OpenMM::Vec3 >& inputGrid );
extern "C"
void gpuLoadElectrostaticPotential( amoebaGpuContext amoebaGpu, unsigned int gridSize, std::vector< double >& outputElectrostaticPotential );
extern "C"
void gpuCleanupElectrostaticPotentialCalculation( amoebaGpuContext amoebaGpu );

extern "C"
void gpuSetAmoebaObcParameters( amoebaGpuContext amoebaGpu , float innerDielectric, float solventDielectric,
                                const std::vector<float>& radius, const std::vector<float>& scale, const std::vector<float>& charge,
                                int includeCavityTerm, float probeRadius, float surfaceAreaFactor);

extern "C"
void gpuSetAmoebaGrycukParameters( amoebaGpuContext amoebaGpu , float innerDielectric, float solventDielectric,
                                   const std::vector<float>& radius, const std::vector<float>& scale, const std::vector<float>& charge,
                                   int includeCavityTerm, float probeRadius, float surfaceAreaFactor);

extern "C"
void gpuSetAmoebaVdwParameters( amoebaGpuContext amoebaGpu,
                                const std::vector<int>& indexIVs, 
                                const std::vector<float>& sigmas,
                                const std::vector<float>& epsilons,
                                const std::vector<float>& reductions,
                                const std::string& sigmaCombiningRule,
                                const std::string& epsilonCombiningRule,
                                const std::vector< std::vector<int> >& allExclusions, int usePBC, float cutoff );
extern "C"
void gpuSetAmoebaPMEParameters(amoebaGpuContext amoebaGpu, float alpha, int gridSizeX, int gridSizeY, int gridSizeZ);

extern "C"
void amoebaGpuBuildVdwExclusionList( amoebaGpuContext amoebaGpu );

extern "C"
void gpuSetAmoebaWcaDispersionParameters( amoebaGpuContext amoebaGpu,
                                const std::vector<float>& radii,
                                const std::vector<float>& epsilons,
                                const float totalMaxWcaDisperionEnergy,
                                const float epso, const float epsh, const float rmino, const float rminh,
                                const float awater, const float shctd, const float dispoff );
 
extern "C"
void amoebaGpuSetConstants(amoebaGpuContext gpu, int updateFlag );

extern "C"
void gpuSetAmoebaBondOffsets(amoebaGpuContext gpu);

#endif //__AMOEBA_GPUTYPES_H__
