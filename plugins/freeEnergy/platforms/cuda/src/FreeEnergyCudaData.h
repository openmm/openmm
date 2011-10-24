#ifndef FREE_ENERGY_CUDA_DATA_H_
#define FREE_ENERGY_CUDA_DATA_H_

/* -------------------------------------------------------------------------- *
 *                              OpenMMFreeEnergy                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors:                                                                   *
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

#include "CudaPlatform.h"
#include "kernels/freeEnergyGpuTypes.h"
#include "kernels/cudaKernels.h"
#include "openmm/KernelImpl.h"

namespace OpenMM {

/**
 * Free energy Cuda data
 */
class FreeEnergyCudaData {

public:

    FreeEnergyCudaData( CudaPlatform::PlatformData& data );
    ~FreeEnergyCudaData();

    /**
     * Increment kernel count
     * 
     */
    void incrementKernelCount( void );

    /**
     * Decrement kernel count
     * 
     */
    void decrementKernelCount( void );

    /**
     * Return freeEnergyGpuContext context
     * 
     * @return freeEnergyGpuContext
     */
    freeEnergyGpuContext OPENMMCUDA_EXPORT getFreeEnergyGpu( void ) const;

    /**
     * Set log file reference
     * 
     * @param log file reference; if not set, then no logging
     */
    void setLog( FILE* inputLog );

    /**
     * Get log file reference
     * 
     * @return log file reference
     */
    FILE* getLog( void ) const;

    /**
     * if gpuInitialized is false, write data to board
     * 
     * @param log file reference; if not set, then no logging
     */
    void initializeGpu( void ); 

    /**
     * Set contextImpl
     * 
     * @param contextImpl reference
     */
    void setContextImpl( void* contextImpl ); 

    CudaPlatform::PlatformData& cudaPlatformData;

private:

    freeEnergyGpuContext freeEnergyGpu;
    unsigned int kernelCount;
    void* contextImpl;
    FILE* log;
    bool gpuInitialized;
    double boxDimensions[3];
};

} // namespace OpenMM

#endif /*FREE_ENERGY_CUDA_DATA_H_*/
