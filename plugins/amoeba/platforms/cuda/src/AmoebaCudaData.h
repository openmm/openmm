#ifndef AMOEBA_CUDA_DATA_H_
#define AMOEBA_CUDA_DATA_H_

/* -------------------------------------------------------------------------- *
 *                              AmoebaOpenMM                                  *
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
#include "kernels/amoebaGpuTypes.h"
#include "kernels/cudaKernels.h"
#include "openmm/KernelImpl.h"

namespace OpenMM {

/**
 * This kernel is invoked by AmoebaHarmonicBondForce to calculate the forces acting on the system and the energy of the system.
 */
class AmoebaCudaData {
public:

    AmoebaCudaData( CudaPlatform::PlatformData& data );
    ~AmoebaCudaData();

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
     * Set value of hasAmoebaBonds
     * 
     * @param value of hasAmoebaBonds
     */
    void setHasAmoebaBonds( bool inputHasAmoebaBonds );

    /**
     * Set value of hasAmoebaMultipole
     * 
     * @param value of hasAmoebaMultipole
     */
    void setHasAmoebaMultipole( bool inputHasAmoebaMultipole );

    /**
     * Get value of hasAmoebaMultipole
     * 
     * @return value of hasAmoebaMultipole
     */
    bool getHasAmoebaMultipole( void ) const;

    /**
     * Set value of hasAmoebaGeneralizedKirkwood
     * 
     * @param value of hasAmoebaGeneralizedKirkwood
     */
    void setHasAmoebaGeneralizedKirkwood( bool inputHasAmoebaGeneralizedKirkwood );

    /**
     * Get value of hasAmoebaGeneralizedKirkwood
     * 
     * @return value of hasAmoebaGeneralizedKirkwood
     */
    bool getHasAmoebaGeneralizedKirkwood( void ) const;

    /**
     * Return amoebaGpuContext context
     * 
     * @return amoebaGpuContext
     */
    amoebaGpuContext OPENMMCUDA_EXPORT getAmoebaGpu( void ) const;

    /**
     * Set accessor for LocalForcesKernel
     * 
     * @param inputLocalForceKernel 
     */
    void setAmoebaLocalForcesKernel( KernelImpl* inputLocalForceKernel );

    /**
     * Get accessor for LocalForcesKernel
     * 
     * @return localForceKernel 
     */
    KernelImpl* getAmoebaLocalForcesKernel( void ) const;
    
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

    /**
     * Get multipole force count 
     * 
     * @return multipole force count
     */
    int getMultipoleForceCount( void ) const; 

    /**
     * Get multipole force count 
     * 
     * @return multipole force count
     */
    void incrementMultipoleForceCount( void ); 

    /**
     * Get multipole force count 
     * 
     * @return multipole force count
     */
    int getApplyCutoff( ) const; 

    /**
     * Get multipole force count 
     * 
     * @return multipole force count
     */
    void setApplyCutoff( int applyCutoff ); 

    CudaPlatform::PlatformData& cudaPlatformData;

private:

    amoebaGpuContext amoebaGpu;
    bool hasAmoebaBonds, hasAmoebaGeneralizedKirkwood, hasAmoebaMultipole;
    int multipoleForceCount;
    int applyCutoff;
    KernelImpl* localForceKernel;
    unsigned int kernelCount;
    void* contextImpl;
    FILE* log;
    bool gpuInitialized;
};


} // namespace OpenMM

#endif /*AMOEBA_CUDA_DATA_H_*/
