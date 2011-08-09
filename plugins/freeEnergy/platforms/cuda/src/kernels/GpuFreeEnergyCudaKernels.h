#ifndef __GPU_FREE_ENERGY_KERNELS_H__
#define __GPU_FREE_ENERGY_KERNELS_H__

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

#include "gputypes.h"
#include "cudatypes.h"
#include "GpuNonbondedSoftcore.h"
#include "GpuLJ14Softcore.h"
#include "GpuObcGbsaSoftcore.h"
#include "GpuGBVISoftcore.h"
#include <vector>
#include <cuda.h>

// Function prototypes

// CDLJ softcore

// setup methods called from CudaFreeEnergyKernels
// nonbonded and 1-4 ixns

extern "C"
bool gpuIsAvailableSoftcore();

extern "C" 
GpuNonbondedSoftcore* gpuSetNonbondedSoftcoreParameters(gpuContext gpu, float epsfac, const std::vector<int>& atom, const std::vector<float>& c6,
                                                        const std::vector<float>& c12, const std::vector<float>& q,
                                                        const std::vector<float>& softcoreLJLambdaArray, const std::vector<char>& symbol,
                                                        const std::vector<std::vector<int> >& exclusions, CudaNonbondedMethod method);

extern "C"
GpuLJ14Softcore* gpuSetLJ14SoftcoreParameters(gpuContext gpu, float epsfac, float fudge, const std::vector<int>& atom1,
                                              const std::vector<int>& atom2, const std::vector<float>& c6, const std::vector<float>& c12,
                                              const std::vector<float>& q1, const std::vector<float>& q2, const std::vector<float>& softcoreLJLambdaArray);

// delete supplemtentary objects, ...

extern "C"
void gpuDeleteNonbondedSoftcoreParameters( void* gpuNonbondedSoftcore);

// write address's to device

extern "C"
void SetCalculateCDLJSoftcoreGpuSim( gpuContext gpu );

extern "C"
void SetCalculateCDLJSoftcoreSupplementarySim( float* gpuParticleSoftCoreLJLambda);

extern "C"
void SetCalculateLocalSoftcoreGpuSim( gpuContext gpu );

// kernel calls to device

extern "C"
void kCalculateCDLJSoftcoreForces(gpuContext gpu );

extern void kCalculateLocalSoftcoreForces( gpuContext gpu );

// GB/VI softcore

// setup method called from CudaFreeEnergyKernels

extern "C" 
GpuGBVISoftcore* gpuSetGBVISoftcoreParameters(gpuContext gpu, float innerDielectric, float solventDielectric, const std::vector<int>& atom, const std::vector<float>& radius, 
                                              const std::vector<float>& gamma, const std::vector<float>& scaledRadii,
                                              const std::vector<float>& bornRadiusScaleFactors, const std::vector<float>& quinticSplineParameters);

// write address's to device

extern "C"
void SetCalculateGBVISoftcoreForcesSim( gpuContext gpu, float* softCoreLJLambda);

extern "C"
void SetCalculateGBVISoftcoreBornSumGpuSim( gpuContext gpu);

extern "C"
void SetCalculateGBVISoftcoreSupplementarySim( GpuGBVISoftcore* gpuGBVISoftcore );

extern "C"
void SetCalculateGBVISoftcoreForces2Sim(gpuContext gpu);

extern "C"
void GetCalculateGBVISoftcoreForces2Sim(gpuContext gpu);

// kernel calls to device

extern void kReduceGBVIBornSumQuinticScaling( gpuContext gpu, GpuGBVISoftcore* gpuGBVISoftcore );
extern void kCalculateGBVISoftcoreBornSum( gpuContext gpu );
extern void kReduceGBVIBornForcesQuinticScaling( gpuContext gpu );
extern void kCalculateGBVISoftcoreForces2( gpuContext gpu );
extern void kReduceGBVISoftcoreBornForces(gpuContext gpu);
extern void kReduceGBVISoftcoreBornSum(gpuContext gpu);
extern void kPrintGBVISoftcore(gpuContext gpu, GpuGBVISoftcore* gpuGBVISoftcore, std::string callId, int call);

extern void kClearSoftcoreBornForces(gpuContext gpu);

// Obc softcore

// setup method called from CudaFreeEnergyKernels

/** 
 * Initialize parameters for Cuda Obc softcore
 * 
 * @param gpu                  gpu context
 * @param innerDielectric      solute dielectric
 * @param solventDielectric    solvent dielectric
 * @param radius               intrinsic Born radii
 * @param scale                Obc scaling factors
 * @param charge               atomic charges (possibly overwritten by other methods?)
 * @param nonPolarScalingFactors non-polar scaling factors
 *
 */

extern "C" 
GpuObcGbsaSoftcore* gpuSetObcSoftcoreParameters(gpuContext gpu, float innerDielectric, float solventDielectric, float nonPolarPrefactor,
                                                const std::vector<float>& radius, const std::vector<float>& scale,
                                                const std::vector<float>& charge, const std::vector<float>& nonPolarScalingFactors);

// delete supplemtentary objects, ...

extern "C"
void gpuDeleteObcSoftcoreParameters( void* gpuNonbondedSoftcore);

// write address's to device

extern "C"
void SetCalculateObcGbsaSoftcoreBornSumSim( gpuContext gpu );

extern "C"
void SetCalculateObcGbsaSoftcoreNonPolarScalingFactorsSim( float* nonPolarScalingFactors );

extern "C"
void SetCalculateObcGbsaSoftcoreNonPolarScalingFactorsObc2Sim( float* nonPolarScalingFactors );

// this method and kCalculateObcGbsaSoftcoreForces2() are being
// used until changes in OpenMM version are made 
extern "C"
void SetCalculateObcGbsaSoftcoreForces2Sim( gpuContext gpu );

// kernel calls to device

extern void kClearObcGbsaSoftcoreBornSum( gpuContext gpu );
extern void kReduceObcGbsaSoftcoreBornForces( gpuContext gpu );
extern void kCalculateObcGbsaSoftcoreBornSum( gpuContext gpu );
extern void kReduceObcGbsaSoftcoreBornSum( gpuContext gpu );


// this method is not needed; the OpenMM version can be used
extern void kCalculateObcGbsaSoftcoreForces2( gpuContext gpu );
extern void kPrintForces( gpuContext gpu, std::string idString, int call );


// shared

extern "C"
void SetCalculateCDLJObcGbsaSoftcoreGpu1Sim( gpuContext gpu );

extern "C"
void SetCalculateCDLJObcGbsaSoftcoreSupplementary1Sim( float* gpuParticleSoftCoreLJLambda);

extern void kCalculateCDLJObcGbsaSoftcoreForces1( gpuContext gpu );

#endif //__GPU_FREE_ENERGY_KERNELS_H__
