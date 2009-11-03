#ifndef OPENMM_FREE_ENERGY_GPU_LJ14_SOFTCORE_
#define OPENMM_FREE_ENERGY_GPU_LJ14_SOFTCORE_

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

struct cudaFreeEnergySimulationNonbonded14 {

    unsigned int    LJ14s;                          // Number of Lennard Jones 1-4 interactions
    unsigned int    LJ14_offset;                    // Offset to end of Lennard Jones 1-4 parameters
    CudaNonbondedMethod nonbondedMethod;            // How to handle nonbonded interactions
    int4*           pLJ14ID;                        // Lennard Jones 1-4 atom and output buffer IDs
    float4*         pLJ14Parameter;                 // Lennard Jones 1-4 parameters
};

// info related to nonbonded 1-4 softcore

class GpuLJ14Softcore {

    public:

        GpuLJ14Softcore();
        ~GpuLJ14Softcore();

        CUDAStream<int4>* psLJ14SoftcoreID;
        CUDAStream<float4>* psLJ14SoftcoreParameter;
        cudaFreeEnergySimulationNonbonded14 feSim;
        int flipStrides(gpuContext gpu);

    private:
    
};

#endif // OPENMM_FREE_ENERGY_GPU_LJ14_SOFTCORE_
