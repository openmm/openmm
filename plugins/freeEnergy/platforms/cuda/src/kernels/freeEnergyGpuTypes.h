#ifndef __FREE_ENERGY_GPUTYPES_H__
#define __FREE_ENERGY_GPUTYPES_H__

/* -------------------------------------------------------------------------- *
 *                          OpenMMFreeEnergy                                      *
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
#include "freeEnergyCudaTypes.h"

#include <map>
typedef std::map<int,float> MapIntFloat;
typedef MapIntFloat::const_iterator MapIntFloatCI;

struct _freeEnergyGpuContext {
    
    _gpuContext* gpuContext;
    cudaFreeEnergyGmxSimulation freeEnergySim;
    std::vector<std::vector<int> > exclusions;
    
    CUDAStream<float4>* psSigEps4;
    CUDAStream<int4>*   psLJ14ID;
    CUDAStream<float4>* psLJ14Parameter;
    CUDAStream<float>*  psSwitchDerivative;
    CUDAStream<float>*  psNonPolarScalingFactors;
 
    FILE* log;
};

typedef struct _freeEnergyGpuContext *freeEnergyGpuContext;

// Function prototypes

extern "C" freeEnergyGpuContext freeEnergyGpuInit( _gpuContext* gpu );
extern "C" void freeEnergyGpuShutDown(freeEnergyGpuContext gpu);
extern "C" void freeEnergyGpuSetConstants(freeEnergyGpuContext gpu);
extern "C" unsigned int getThreadsPerBlockFEP( freeEnergyGpuContext freeEnergyGpu, unsigned int sharedMemoryPerThread, unsigned int sharedMemoryPerBlock );

#endif // __FREE_ENERGY_GPUTYPES_H__
